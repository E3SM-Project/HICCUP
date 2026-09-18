# -------------------------------------------------------------------------------------------------
# Pure-Python vertical remap utilities for HICCUP
#
# Replaces NCO's `ncremap --vrt_fl=...` with an xarray + numpy implementation that streams
# columns through dask, so memory does not blow up on large unstructured grids (ne1024, etc.).
# Top-level entry point is `remap_vertical_py()`; the `_py` suffix disambiguates from the
# pre-existing NCO-based `remap_vertical()` method on `hiccup_data` (which we keep around as
# a fallback). The module is intentionally free of `hiccup_data` coupling so it can be unit
# tested in isolation.
# -------------------------------------------------------------------------------------------------
import os
import numpy as np
import xarray as xr
import dask

from hiccup.hiccup_utilities import tcolor

# default reference pressure (Pa); used only when neither input nor target file carries P0
_DEFAULT_P0 = 1.0e5

# variables that describe the vertical grid itself; never remapped, always pulled from
# the target vert_file and copied through verbatim
_HYBRID_COEF_VARS = ('hyam', 'hybm', 'hyai', 'hybi', 'P0')

# lapse rate used when extrapolating temperature below the lowest source level,
# in K/Pa. 6.5e-4 K/Pa is 6.5 K/(100 hPa), the constant moist adiabatic rate NCO's
# ncremap applies to temperature; keeping the same value means the Python path
# reproduces the vertical profiles HICCUP produced via NCO. Change it to ~9.8e-4
# for a dry adiabat.
_LAPSE_RATE = 6.5e-4

# ---------------------------------------------------------------------------
def _is_temperature(da, name):
  """
  identify temperature fields the way NCO does - by units first, then by a short
  list of conventional names - so lapse-rate extrapolation is applied to T alone
  and not to humidity, winds or tracers
  """
  units = str(da.attrs.get('units','')).strip().lower()
  if units in ('k','kelvin','deg_k','degk','degrees_k'): return True
  return str(name).lower() in ('t','t_mid','ta','temp','temperature')

# Memory tuning knobs. These exist because the remap is driven from a single task on a
# single node - for ne1024-class RRM grids (ncol ~ 3e8) the defaults below are the
# difference between a 60 GB job and an OOM kill. See remap_vertical_py() for details.
#
# target uncompressed bytes held per dask chunk (source levels + target levels for one
# horizontal block); the chunk length along ncol is derived from this
chunk_mem_gb = 16.0
# number of concurrent dask threads used while writing the output; peak memory scales
# with this times chunk_mem_gb, so it is deliberately far below the core count
max_workers = 4
# number of columns processed per pass inside the interpolation kernel; bounds the
# transient target-pressure array so it never scales with the dask chunk size
interp_block_size = 65536

# ---------------------------------------------------------------------------
# pressure-on-grid helpers
# ---------------------------------------------------------------------------
def _resolve_surface_pressure(ds, ps_name, lev_name='lev'):
  """
  return surface pressure as an xarray.DataArray named ps_name
  prefers ds[ps_name]; falls back to exp(ds['lnsp']) for ECMWF IFS layout
  raises ValueError if neither is available
  for the lnsp path, drops only the singleton vertical dim (lev_name) when
  present - using a blanket squeeze() would also strip a singleton time dim
  and misalign PS with the rest of the dataset for single-timestep inputs
  """
  if ps_name in ds.variables:
    return ds[ps_name]
  if 'lnsp' in ds.variables:
    lnsp = ds['lnsp']
    if lev_name in lnsp.dims:
      lnsp = lnsp.isel({lev_name: 0}, drop=True)
    return np.exp(lnsp).rename(ps_name)
  raise ValueError(
    f'input_file must contain surface pressure as {ps_name!r} or as lnsp'
  )

# ---------------------------------------------------------------------------
def _input_pressure_spec(ds, lev_name):
  """
  describe the source vertical coordinate so pressure can be rebuilt one block of
  columns at a time inside the interpolation kernel, instead of materializing the
  full (ncol,lev) pressure array up front
  detects three layouts:
    1. EAM / EAMxx hybrid: hyam is a unitless fraction of P0 -> p = hyam*P0 + hybm*ps
       (P0 defaults to _DEFAULT_P0 when absent - EAM/CAM hyam are unitless
       fractions, so a missing P0 must NOT silently drop us to the bare-lev
       fallback, which mixes hPa lev values with a Pa target grid)
    2. ECMWF IFS hybrid:   hyam is already in Pa -> p = hyam + hybm*ps
    3. pure pressure levels: only the lev coord -> p = ds[lev_name]
       (converted to Pa when the coordinate advertises hPa/millibar units)
  The EAM vs IFS choice is made from the magnitude of hyam, NOT from the presence
  of P0 or lnsp: those cues are ambiguous (an IFS file can carry P0 after
  add_reference_pressure, and an IFS file can supply PS instead of lnsp), and
  keying off them picks the wrong hybrid formula. hyam as a unitless fraction is
  O(1); hyam in Pa is O(1e3-1e4), so the two conventions are cleanly separated by
  several orders of magnitude.
  returns ('plev', da) for layout 3 and ('hybrid', (hyam, hybm, p0)) otherwise,
  where the hybrid form always evaluates as a*p0 + b*ps (p0 is 1.0 for IFS, whose
  hyam is already in Pa)
  """
  variables = set(ds.variables.keys())
  hybrid = {'hyam','hybm'}.issubset(variables)

  if hybrid:
    hyam = ds['hyam']
    hybm = ds['hybm']
    # unitless fraction (EAM/CAM) vs pressure in Pa (ECMWF IFS); hyam is 1-D over
    # levels so this max() is cheap even under dask
    if float(np.asarray(hyam.max())) <= 1.0:
      # EAM / EAMxx hybrid: hyam is a unitless fraction of P0. Default P0 when the
      # file hasn't had add_reference_pressure() applied yet - otherwise we'd fall
      # through to the bare-lev fallback below and use the hPa-valued lev
      # coordinate as if it were Pa, clamping everything below ~10 hPa to a constant.
      p0 = ds['P0'] if 'P0' in variables else _DEFAULT_P0
      return 'hybrid', (hyam, hybm, p0)
    # IFS layout: hyam is already in Pa, no P0 scaling
    return 'hybrid', (hyam, hybm, 1.0)

  if lev_name in ds.coords or lev_name in ds.variables:
    lev = ds[lev_name].astype('float64')
    # pressure-level coordinate: normalize hPa/millibar to Pa so it matches the
    # Pa-valued target grid (p_out = hyam*P0 + hybm*ps)
    units = str(lev.attrs.get('units','')).strip().lower()
    if units in ('hpa','mb','millibar','millibars','mbar'):
      lev = lev*100.0
    return 'plev', lev

  raise ValueError(
    f'cannot infer source vertical grid for lev_name={lev_name!r}; '
    f'expected hybrid (hyam,hybm[,P0|lnsp]) or pressure-level coordinate'
  )

# ---------------------------------------------------------------------------
def _compute_input_pressure(ds, lev_name, ps):
  """
  return an xarray.DataArray of pressure on the source vertical grid
  ps must be the resolved surface pressure DataArray (see _resolve_surface_pressure)
  layout detection lives in _input_pressure_spec()
  """
  kind, spec = _input_pressure_spec(ds, lev_name)
  if kind == 'plev':
    return spec
  hyam, hybm, p0 = spec
  return hyam*p0 + hybm*ps

# ---------------------------------------------------------------------------
def _output_pressure_spec(ds_vert, out_lev_name):
  """
  return (hyam, hybm, p0) for the target hybrid grid, with the coefficients renamed
  onto out_lev_name; target pressure always evaluates as hyam*p0 + hybm*ps
  """
  if not {'hyam','hybm'}.issubset(ds_vert.variables.keys()):
    raise ValueError('ds_vert must contain hyam and hybm')
  p0 = ds_vert['P0'] if 'P0' in ds_vert.variables else _DEFAULT_P0
  hyam = ds_vert['hyam']
  hybm = ds_vert['hybm']
  # ensure hybrid coefs use the canonical output lev dim name
  if hyam.dims[0] != out_lev_name:
    hyam = hyam.rename({hyam.dims[0]: out_lev_name})
    hybm = hybm.rename({hybm.dims[0]: out_lev_name})
  return hyam, hybm, p0

# ---------------------------------------------------------------------------
def _compute_output_pressure(ds_vert, ps, out_lev_name):
  """
  return pressure on the target hybrid grid: p = hyam*P0 + hybm*PS
  ps must already be aligned with the input dataset's non-vertical dims
  """
  hyam, hybm, p0 = _output_pressure_spec(ds_vert, out_lev_name)
  return hyam*p0 + hybm*ps

# ---------------------------------------------------------------------------
# interpolation kernel
# ---------------------------------------------------------------------------
def _interp_block(p_target, p_source, f_source, mode='log_pressure',
                  extrap='constant', lapse_rate=None):
  """
  interpolate a block of columns from p_source -> p_target, fully vectorized
    p_target  (..., n_out)  target pressures, one column per leading index
    p_source  (n_in,) shared by every column, or (..., n_in) per column
    f_source  (..., n_in)   source field values
  - mode:    'log_pressure' (default) or 'linear_pressure'
  - extrap:  'constant' (clamp to the end values), 'linear' (slope from the end
             pair), or 'lapse' (clamp at the top, but extrapolate below the
             lowest source level along a constant lapse rate in pressure)
  - lapse_rate: K/Pa used by extrap='lapse'; defaults to _LAPSE_RATE
  Columns are sorted by pressure defensively, so top-down and bottom-up source
  orderings both work. This replaces a per-column np.interp loop driven by
  np.vectorize: on an ne1024-class grid that loop ran ~3e8 times per variable,
  held the GIL throughout, and dominated the runtime.
  """
  if lapse_rate is None: lapse_rate = _LAPSE_RATE

  p_t = np.asarray(p_target, dtype='float64')
  p_s = np.asarray(p_source, dtype='float64')

  if mode == 'log_pressure':
    x, xp = np.log(p_t), np.log(p_s)
  else:
    x, xp = p_t, p_s

  fp = f_source
  shared_source = (xp.ndim == 1)

  # sort each column by source pressure; np.interp required increasing xp and the
  # gather below has the same requirement. the linear pressures are carried along
  # because the lapse-rate extrapolation below is defined in Pa, not log(Pa)
  if shared_source:
    order = np.argsort(xp)
    xp = xp[order]
    p_s = p_s[order]
    fp = fp[..., order]
  else:
    order = np.argsort(xp, axis=-1)
    xp = np.take_along_axis(xp, order, axis=-1)
    p_s = np.take_along_axis(p_s, order, axis=-1)
    fp = np.take_along_axis(fp, order, axis=-1)

  n_in = xp.shape[-1]

  # a single source level defines no slope and no interval - every target level
  # takes that value, which matches np.interp's behavior for a one-point table
  if n_in == 1:
    return np.broadcast_to(fp[..., :1], x.shape).copy()

  # locate each target pressure in the source column: idx = number of source
  # levels strictly below it, i.e. np.searchsorted(..., side='left')
  if shared_source:
    idx = np.searchsorted(xp, x, side='left')
  else:
    # accumulate the count one source level at a time rather than building the
    # full (..., n_out, n_in) comparison array, which would dwarf the data itself
    idx = np.zeros(x.shape, dtype=np.intp)
    for k in range(n_in):
      idx += (x > xp[..., k:k+1])

  # clamping to a valid interior pair makes the gather below extrapolate along the
  # slope of the first/last pair, which is exactly the 'linear' extrap behavior;
  # 'constant' overrides the out-of-range values afterwards
  i1 = np.clip(idx, 1, n_in-1)
  i0 = i1 - 1

  if shared_source:
    x0, x1 = xp[i0], xp[i1]
  else:
    x0 = np.take_along_axis(xp, i0, axis=-1)
    x1 = np.take_along_axis(xp, i1, axis=-1)
  y0 = np.take_along_axis(fp, i0, axis=-1)
  y1 = np.take_along_axis(fp, i1, axis=-1)

  # duplicate source pressures give a zero-width interval; fall back to the lower
  # end value there instead of dividing by zero
  dx = x1 - x0
  w  = np.divide(x - x0, dx, out=np.zeros_like(x), where=(dx != 0))
  y  = y0 + w*(y1 - y0)

  if extrap in ('constant','lapse'):
    if shared_source:
      lo_p, hi_p = xp[0], xp[-1]
    else:
      lo_p, hi_p = xp[..., :1], xp[..., -1:]
    y = np.where(x < lo_p, fp[..., :1], y)
    y = np.where(x > hi_p, fp[..., -1:], y)

  if extrap == 'lapse':
    # below the lowest source level, follow a constant lapse rate in pressure
    # instead of holding the end value fixed. This reproduces what NCO's ncremap
    # reports doing ("temperature extrapolated toward/into surface assuming
    # constant moist adiabatic lapse rate = 6.5 K/(100 hPa)"), and it matters
    # because the bottom hybrid levels sit below the lowest source pressure
    # wherever surface pressure exceeds it - roughly a third of a global grid,
    # where clamping instead leaves the lowest model levels systematically cool.
    p_bot = p_s[-1] if shared_source else p_s[..., -1:]
    y = np.where(p_t > p_bot, fp[..., -1:] + lapse_rate*(p_t - p_bot), y)

  return y

# ---------------------------------------------------------------------------
def _interp_column(p_target, p_source, f_source, mode='log_pressure',
                   extrap='constant', lapse_rate=None):
  """
  interpolate a single column from p_source -> p_target
  thin wrapper around _interp_block() so the single-column and block paths cannot
  drift apart
  """
  return _interp_block(p_target, p_source, f_source, mode=mode, extrap=extrap,
                       lapse_rate=lapse_rate)

# ---------------------------------------------------------------------------
def _remap_field(field, ps, p_in_spec, p_out_spec, in_lev_name, out_lev_name,
                 mode, extrap, block_size=None):
  """
  remap a single xarray.DataArray from in_lev_name -> out_lev_name
  driven by xarray.apply_ufunc so dask parallelizes over blocks of columns
    field       source field, chunked with the full column in a single chunk
    ps          surface pressure, chunked to match field's horizontal dims
    p_in_spec   source vertical coordinate, from _input_pressure_spec()
    p_out_spec  target (hyam, hybm, p0), from _output_pressure_spec()
  Both pressure coordinates are rebuilt from ps inside the kernel rather than
  passed in as arrays. On an ne1024-class grid the target pressure alone is
  (3e8, 128) float64 = 300 GB; as a dask input it materialized ~1 GB per chunk
  per worker, which is what pushed the job over the node's 512 GB.
  """
  if block_size is None: block_size = interp_block_size

  out_dtype = field.dtype if np.issubdtype(field.dtype, np.floating) else np.float64

  hyam_o, hybm_o, p0_o = p_out_spec
  a_out = np.asarray(hyam_o.values, dtype='float64')
  b_out = np.asarray(hybm_o.values, dtype='float64')
  p0_out = float(np.asarray(p0_o))
  nlev_out = a_out.size

  in_kind, in_spec = p_in_spec
  if in_kind == 'plev':
    p_in_shared = np.asarray(in_spec.values, dtype='float64')
    a_in = b_in = None
    p0_in = None
  else:
    hyam_i, hybm_i, p0_i = in_spec
    p_in_shared = None
    a_in = np.asarray(hyam_i.values, dtype='float64')
    b_in = np.asarray(hybm_i.values, dtype='float64')
    p0_in = float(np.asarray(p0_i))

  def _kernel(f_s, ps_v):
    # apply_ufunc hands us the core dim last: f_s is (..., n_in), ps_v is (...)
    lead  = f_s.shape[:-1]
    n_col = int(np.prod(lead)) if lead else 1
    f_flat  = f_s.reshape(n_col, f_s.shape[-1])
    ps_flat = np.asarray(ps_v, dtype='float64').reshape(-1)
    if ps_flat.size != n_col:
      ps_flat = np.broadcast_to(ps_flat, (n_col,))

    out = np.empty((n_col, nlev_out), dtype=out_dtype)
    for i0 in range(0, n_col, block_size):
      i1 = min(i0+block_size, n_col)
      ps_b = ps_flat[i0:i1, None]
      p_o  = a_out*p0_out + b_out*ps_b
      p_i  = p_in_shared if p_in_shared is not None else a_in*p0_in + b_in*ps_b
      out[i0:i1] = _interp_block(p_o, p_i, f_flat[i0:i1],
                                 mode=mode, extrap=extrap).astype(out_dtype)
    return out.reshape(*lead, nlev_out)

  result = xr.apply_ufunc(
    _kernel,
    field, ps,
    input_core_dims=[[in_lev_name], []],
    output_core_dims=[[out_lev_name]],
    dask='parallelized',
    output_dtypes=[out_dtype],
    dask_gufunc_kwargs={'output_sizes': {out_lev_name: nlev_out}},
  )
  result.attrs = dict(field.attrs)

  # apply_ufunc always moves core dims to the last axis, so a (time,lev,ncol) field
  # would come back as (time,ncol,lev). The NCO path this replaced was order-preserving,
  # and EAM expects the vertical dim in its original slot, so restore the source layout.
  # The trailing Ellipsis is defensive - it parks any dim not present in the source
  # field (from broadcasting against ps) at the end instead of raising.
  dim_order = [out_lev_name if d==in_lev_name else d for d in field.dims]
  result = result.transpose(*dim_order, ...)

  return result

# ---------------------------------------------------------------------------
# chunk selection
# ---------------------------------------------------------------------------
def _storage_chunk_len(ds, dim, var_list):
  """
  return the netCDF storage chunk length along dim, taken as the largest over the
  variables being remapped (0 when they are all contiguous/unchunked)
  """
  store_len = 0
  for v in var_list:
    if v not in ds.variables: continue
    chunksizes = ds[v].encoding.get('chunksizes', None)
    if chunksizes is None: continue
    if dim not in ds[v].dims: continue
    store_len = max(store_len, int(chunksizes[ds[v].dims.index(dim)]))
  return store_len

# ---------------------------------------------------------------------------
def _auto_chunk_len(ds, dim, lev_name, nlev_out, var_list, budget_bytes):
  """
  choose a chunk length along dim that stays inside budget_bytes and is a whole
  multiple of the netCDF storage chunk length
  dask's own 'auto' sizing is unusable here for two reasons: it sizes chunks from
  the source array, ignoring that the target grid can have several times as many
  levels, and it ignores the on-disk chunking - a deflated variable stored in
  (1,2,21428512) chunks gets decompressed once per dask chunk that touches it, so
  unaligned chunks turn one read of the variable into hundreds
  """
  dim_len = int(ds.sizes[dim])
  nlev_in = int(ds.sizes.get(lev_name, 1))

  itemsize = 4
  for v in var_list:
    if v in ds.variables: itemsize = max(itemsize, ds[v].dtype.itemsize)

  # a column costs the source levels read in plus the target levels written out
  bytes_per_col = max(1, (nlev_in + nlev_out)*itemsize)
  budget_cols = max(1, int(budget_bytes // bytes_per_col))

  store_len = _storage_chunk_len(ds, dim, var_list)
  if 0 < store_len < dim_len:
    return int(min(dim_len, max(1, budget_cols//store_len)*store_len))

  return int(min(dim_len, budget_cols))

# ---------------------------------------------------------------------------
# top-level entry point
# ---------------------------------------------------------------------------
def remap_vertical_py(input_file, output_file, vert_file,
                      ps_name='PS', var_list=None, lev_name='lev',
                      mode='log_pressure', extrap='constant',
                      chunks=None, nc_output_format='NETCDF4',
                      mem_gb=None, workers=None, ps_file=None,
                      temperature_extrap='lapse', verbose=False):
  """
  vertically remap fields in input_file onto the hybrid grid defined by vert_file
  and write the result to output_file

  supported source layouts (auto-detected):
    - EAM / EAMxx hybrid:  hyam, hybm, P0, and ps_name present
    - ECMWF IFS hybrid:    hyam (in Pa), hybm, and lnsp present (no P0)
    - pure pressure level: only the lev_name coordinate (e.g. ERA5 plev in Pa)
  surface pressure is taken from ps_name when available, else derived as exp(lnsp);
  in either case the resolved value is written to the output as ps_name

  parameters
    input_file   path to source dataset; must contain lev_name plus a recognized layout
                 (see above) and either ps_name or lnsp for surface pressure
    output_file  path to write remapped dataset
    vert_file    target vertical grid file (must contain hyam, hybm; P0 optional)
    ps_name      name of surface pressure variable in input_file (default 'PS');
                 also the name used for surface pressure in the output
    ps_file      optional separate file to read surface pressure from; lets the
                 multi-file workflow skip copying ps into every per-variable file,
                 which costs a full read and rewrite of each file before any
                 remapping happens. The resolved value is still written to the
                 output as ps_name, so the output layout is unchanged
    var_list     list of variables to remap; if None, remap every var with lev_name in dims
    lev_name     name of source vertical dim (default 'lev')
    mode         'log_pressure' (default) or 'linear_pressure'
    extrap       'constant' (default, clamp), 'linear' (slope from end pair), or
                 'lapse' (clamp at top, lapse rate below the lowest source level)
    temperature_extrap
                 extrapolation used for temperature fields specifically, which
                 NCO treats differently from everything else; defaults to 'lapse'
                 so the output matches what ncremap produced. Set to None to use
                 `extrap` for every variable
    chunks       dict passed to xr.open_dataset; lev_name is forced to -1 (full
                 column required for interp). Any dim left as 'auto' - and every
                 non-lev dim when chunks is None - is sized by _auto_chunk_len()
                 from mem_gb and the file's own storage chunking; pass explicit
                 integers to override
    nc_output_format  netCDF format string passed through to xarray's to_netcdf()
                      (default 'NETCDF4'); set to 'NETCDF3_64BIT' etc. when downstream
                      tools require an older format
    mem_gb       target uncompressed bytes per dask chunk, in GB (default chunk_mem_gb)
    workers      concurrent dask threads used for the write (default max_workers);
                 peak memory is roughly workers*mem_gb, so leaving this at the core
                 count is what turns a large remap into an OOM kill
    verbose      print progress
  """
  if mode not in ('log_pressure','linear_pressure'):
    raise ValueError(f'mode must be log_pressure or linear_pressure, got {mode!r}')
  if extrap not in ('constant','linear','lapse'):
    raise ValueError(f'extrap must be constant, linear or lapse, got {extrap!r}')
  if temperature_extrap is not None \
  and temperature_extrap not in ('constant','linear','lapse'):
    raise ValueError(f'temperature_extrap must be constant, linear, lapse or '
                     f'None, got {temperature_extrap!r}')

  if mem_gb  is None: mem_gb  = chunk_mem_gb
  if workers is None: workers = max_workers

  same_in_out = os.path.abspath(input_file) == os.path.abspath(output_file)
  tmp_file = output_file + '.vrt_tmp.nc' if same_in_out else output_file

  # peek at the source to size chunks: lev_name is always forced to -1 (the full
  # column is needed for interp) and every dim left as 'auto' is sized against the
  # memory budget and the file's storage chunking rather than dask's default
  with xr.open_dataset(input_file) as ds_peek, \
       xr.open_dataset(vert_file) as ds_vert_peek:
    nlev_out = int(ds_vert_peek['hyam'].shape[0])
    if chunks is None:
      chunks = {d: 'auto' for d in ds_peek.dims}
    else:
      chunks = dict(chunks)
    chunks[lev_name] = -1
    chunk_var_list = var_list
    if chunk_var_list is None:
      chunk_var_list = [v for v in ds_peek.data_vars
                        if lev_name in ds_peek[v].dims and v not in _HYBRID_COEF_VARS]
    for d, c in list(chunks.items()):
      if c != 'auto': continue
      if d not in ds_peek.dims:
        chunks.pop(d)
        continue
      chunks[d] = _auto_chunk_len(ds_peek, d, lev_name, nlev_out,
                                  chunk_var_list, mem_gb*(1024**3))

  if verbose:
    print(f'{tcolor.GREEN}  vertical remap: {input_file} -> {output_file}{tcolor.ENDC}')
    print(f'    chunks: {chunks}  workers: {workers}')

  with xr.open_dataset(input_file, chunks=chunks) as ds_in, \
       xr.open_dataset(vert_file) as ds_vert:

    # resolve surface pressure once (handles both PS and IFS lnsp layouts); when
    # ps_file is given it is read from there instead, so the caller does not have
    # to copy ps into the input file first
    ds_ps = None
    if ps_file is None:
      ps = _resolve_surface_pressure(ds_in, ps_name, lev_name=lev_name)
    else:
      ds_ps = xr.open_dataset(ps_file)
      ds_ps = ds_ps.chunk({d:c for d,c in chunks.items() if d in ds_ps.dims})
      ps = _resolve_surface_pressure(ds_ps, ps_name, lev_name=lev_name)

    out_lev_name_native = ds_vert['hyam'].dims[0]
    # if source and target use the same dim name, rename target internally so apply_ufunc
    # can tell them apart; we restore the native name at the end via the output dataset
    out_lev_name = out_lev_name_native
    if out_lev_name == lev_name:
      out_lev_name = f'{lev_name}_target'
      ds_vert = ds_vert.rename({out_lev_name_native: out_lev_name})

    p_in_spec  = _input_pressure_spec(ds_in, lev_name)
    p_out_spec = _output_pressure_spec(ds_vert, out_lev_name)

    # guard against a source/target vertical-unit mismatch. The classic failure
    # is a source file whose hybrid coefficients (or P0) went missing, so the
    # source pressure ends up as a bare hPa lev coordinate (max ~1000) while the
    # target grid is in Pa (max ~1e5). np.interp would then clamp every level
    # below ~10 hPa to a constant instead of raising - fail loudly instead.
    # Only relevant to the bare pressure-coordinate fallback: a hybrid source
    # always yields Pa (hyam*P0 + hybm*ps), so skip the guard there and avoid
    # eagerly reducing p_in/p_out over large dask arrays before apply_ufunc.
    p_in_kind, p_in_plev = p_in_spec
    if p_in_kind == 'plev':
      # p_in_plev is the 1-D lev coordinate (already in Pa if units were hPa/mb)
      p_in_max = float(np.asarray(p_in_plev.max()))
      ps_max = float(np.asarray(ps.max()))
      p0 = float(np.asarray(ds_vert['P0'])) if 'P0' in ds_vert.variables else _DEFAULT_P0
      hyam_max = float(np.asarray(ds_vert['hyam'].max()))
      hybm_max = float(np.asarray(ds_vert['hybm'].max()))
      p_out_max = hyam_max*p0 + hybm_max*ps_max
      if p_in_max > 0.0 and p_in_max < 1.2e3 and p_out_max > 1.0e4:
        raise ValueError(
          f'source vertical pressure (max {p_in_max:.3g}) and target grid '
          f'(max {p_out_max:.3g} Pa) appear to be in different units - the source '
          f'looks like hPa while the target is Pa. This usually means the source '
          f'file is missing its hybrid coefficients (hyam/hybm/P0) at the vertical '
          f'remap step, so lev={lev_name!r} was used directly as pressure. '
          f'Fix by adding P0 (add_reference_pressure), setting lev units to hPa/mb '
          f'(so it can be converted to Pa), or converting lev values to Pa and '
          f'setting units="Pa".'
        )

    # decide which fields get remapped; never remap the hybrid coefficients themselves -
    # those describe the vertical grid and are pulled from vert_file
    if var_list is None:
      var_list = [v for v in ds_in.data_vars
                  if lev_name in ds_in[v].dims and v not in _HYBRID_COEF_VARS]
    else:
      var_list = [v for v in var_list
                  if v in ds_in.data_vars and v not in _HYBRID_COEF_VARS]

    # build output dataset: target hybrid coefs from vert_file, plus passthrough vars
    ds_out = xr.Dataset()
    for v in _HYBRID_COEF_VARS:
      if v in ds_vert.variables:
        ds_out[v] = ds_vert[v]
    # always carry P0 in the output so the file is self-describing - if vert_file
    # didn't provide one, record the default that was used in the pressure calc
    if 'P0' not in ds_out.variables:
      ds_out['P0'] = xr.DataArray(np.float64(_DEFAULT_P0),
                                  attrs={'long_name': 'reference pressure', 'units': 'Pa'})
    ds_out[ps_name] = ps

    # passthrough: any var without the source lev dim that isn't already in ds_out
    for v in ds_in.data_vars:
      if v in ds_out.variables: continue
      if v in var_list:         continue
      if lev_name not in ds_in[v].dims:
        ds_out[v] = ds_in[v]

    # remap each field
    for v in var_list:
      if lev_name not in ds_in[v].dims:
        # caller asked for a field that isn't on the vertical grid; copy it through
        ds_out[v] = ds_in[v]
        continue
      # NCO extrapolates temperature into the surface along a lapse rate but
      # clamps everything else, so the choice is made per variable
      v_extrap = extrap
      if temperature_extrap is not None and _is_temperature(ds_in[v], v):
        v_extrap = temperature_extrap
      if verbose: print(f'    remapping {v} (extrap={v_extrap})')
      ds_out[v] = _remap_field(ds_in[v], ps, p_in_spec, p_out_spec,
                               lev_name, out_lev_name, mode, v_extrap)

    # carry over global attrs and coords that aren't tied to the source vertical dim
    ds_out.attrs = dict(ds_in.attrs)

    # restore the native target lev dim name if we renamed it to avoid collision
    if out_lev_name != out_lev_name_native:
      ds_out = ds_out.rename({out_lev_name: out_lev_name_native})

    # cap the thread pool for the write: the netCDF writer is serialized, so an
    # uncapped scheduler (one thread per core - 256 on a Perlmutter CPU node) just
    # computes chunks faster than they can be drained and holds them all in memory
    with dask.config.set(scheduler='threads', num_workers=workers):
      ds_out.to_netcdf(tmp_file, mode='w', format=nc_output_format)

    if ds_ps is not None: ds_ps.close()

  if same_in_out:
    os.replace(tmp_file, output_file)

  return
# -------------------------------------------------------------------------------------------------
