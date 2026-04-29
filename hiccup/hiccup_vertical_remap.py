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

from hiccup.hiccup_utilities import tcolor

# default reference pressure (Pa); used only when neither input nor target file carries P0
_DEFAULT_P0 = 1.0e5

# variables that describe the vertical grid itself; never remapped, always pulled from
# the target vert_file and copied through verbatim
_HYBRID_COEF_VARS = ('hyam', 'hybm', 'hyai', 'hybi', 'P0')

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
def _compute_input_pressure(ds, lev_name, ps):
  """
  return an xarray.DataArray of pressure on the source vertical grid
  detects three layouts:
    1. EAM / EAMxx hybrid: hyam, hybm, P0 present -> p = hyam*P0 + hybm*ps
    2. ECMWF IFS hybrid:   hyam in Pa, lnsp present (no P0) -> p = hyam + hybm*ps
    3. pure pressure levels: only the lev coord -> p = ds[lev_name]
  ps must be the resolved surface pressure DataArray (see _resolve_surface_pressure)
  """
  variables = set(ds.variables.keys())
  hybrid = {'hyam','hybm'}.issubset(variables)

  if hybrid and 'P0' in variables:
    return ds['hyam']*ds['P0'] + ds['hybm']*ps

  if hybrid and 'lnsp' in variables:
    # IFS layout: hyam is already in Pa
    return ds['hyam'] + ds['hybm']*ps

  if lev_name in ds.coords or lev_name in ds.variables:
    return ds[lev_name].astype('float64')

  raise ValueError(
    f'cannot infer source vertical grid for lev_name={lev_name!r}; '
    f'expected hybrid (hyam,hybm,P0|lnsp) or pressure-level coordinate'
  )

# ---------------------------------------------------------------------------
def _compute_output_pressure(ds_vert, ps, out_lev_name):
  """
  return pressure on the target hybrid grid: p = hyam*P0 + hybm*PS
  ps must already be aligned with the input dataset's non-vertical dims
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
  return hyam*p0 + hybm*ps

# ---------------------------------------------------------------------------
# per-column interpolation kernel
# ---------------------------------------------------------------------------
def _interp_column(p_target, p_source, f_source, mode='log_pressure',
                   extrap='constant'):
  """
  interpolate a single column from p_source -> p_target
  - mode:    'log_pressure' (default) or 'linear_pressure'
  - extrap:  'constant' (clamp, np.interp default) or 'linear' (slope from end pair)
  np.interp requires monotonically increasing xp; we sort defensively per column
  """
  if mode == 'log_pressure':
    x  = np.log(p_target)
    xp = np.log(p_source)
  else:
    x, xp = p_target, p_source

  order = np.argsort(xp)
  xp = xp[order]
  fp = f_source[order]

  y = np.interp(x, xp, fp)

  # linear extrapolation needs a slope from the last two source points; fall back
  # to constant (np.interp's clamp) when the column has fewer than 2 levels
  if extrap == 'linear' and xp.size >= 2:
    left = x < xp[0]
    if left.any():
      slope = (fp[1] - fp[0]) / (xp[1] - xp[0])
      y = np.where(left, fp[0] + (x - xp[0])*slope, y)
    right = x > xp[-1]
    if right.any():
      slope = (fp[-1] - fp[-2]) / (xp[-1] - xp[-2])
      y = np.where(right, fp[-1] + (x - xp[-1])*slope, y)

  return y

# ---------------------------------------------------------------------------
def _remap_field(field, p_in, p_out, in_lev_name, out_lev_name, mode, extrap):
  """
  remap a single xarray.DataArray from in_lev_name -> out_lev_name
  driven by xarray.apply_ufunc with vectorize=True so dask can parallelize over columns
  """
  out_dtype = field.dtype if np.issubdtype(field.dtype, np.floating) else np.float64

  def _kernel(p_t, p_s, f_s):
    return _interp_column(p_t, p_s, f_s, mode=mode, extrap=extrap).astype(out_dtype)

  result = xr.apply_ufunc(
    _kernel,
    p_out, p_in, field,
    input_core_dims=[[out_lev_name], [in_lev_name], [in_lev_name]],
    output_core_dims=[[out_lev_name]],
    vectorize=True,
    dask='parallelized',
    output_dtypes=[out_dtype],
  )
  result.attrs = dict(field.attrs)
  return result

# ---------------------------------------------------------------------------
# top-level entry point
# ---------------------------------------------------------------------------
def remap_vertical_py(input_file, output_file, vert_file,
                      ps_name='PS', var_list=None, lev_name='lev',
                      mode='log_pressure', extrap='constant',
                      chunks=None, nc_output_format='NETCDF4',
                      verbose=False):
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
    var_list     list of variables to remap; if None, remap every var with lev_name in dims
    lev_name     name of source vertical dim (default 'lev')
    mode         'log_pressure' (default) or 'linear_pressure'
    extrap       'constant' (default, clamp) or 'linear' (slope from end pair)
    chunks       dict passed to xr.open_dataset; lev_name is forced to -1 (full
                 column required for interp). When None (default), every non-lev
                 dim is set to 'auto' so dask picks reasonable chunk sizes for
                 large unstructured grids; pass an explicit dict to override
    nc_output_format  netCDF format string passed through to xarray's to_netcdf()
                      (default 'NETCDF4'); set to 'NETCDF3_64BIT' etc. when downstream
                      tools require an older format
    verbose      print progress
  """
  if mode not in ('log_pressure','linear_pressure'):
    raise ValueError(f'mode must be log_pressure or linear_pressure, got {mode!r}')
  if extrap not in ('constant','linear'):
    raise ValueError(f'extrap must be constant or linear, got {extrap!r}')

  same_in_out = os.path.abspath(input_file) == os.path.abspath(output_file)
  tmp_file = output_file + '.vrt_tmp.nc' if same_in_out else output_file

  # build chunks: lev_name is always forced to -1 (full column needed for interp);
  # when chunks is None, default every other dim to dask 'auto' so horizontal/time
  # axes get reasonable chunk sizes instead of one giant chunk per dim
  if chunks is None:
    with xr.open_dataset(input_file) as _ds_peek:
      chunks = {d: 'auto' for d in _ds_peek.dims}
    chunks[lev_name] = -1
  else:
    chunks = dict(chunks)
    chunks[lev_name] = -1

  if verbose:
    print(f'{tcolor.GREEN}  vertical remap: {input_file} -> {output_file}{tcolor.ENDC}')

  with xr.open_dataset(input_file, chunks=chunks) as ds_in, \
       xr.open_dataset(vert_file) as ds_vert:

    # resolve surface pressure once (handles both PS and IFS lnsp layouts)
    ps = _resolve_surface_pressure(ds_in, ps_name, lev_name=lev_name)

    out_lev_name_native = ds_vert['hyam'].dims[0]
    # if source and target use the same dim name, rename target internally so apply_ufunc
    # can tell them apart; we restore the native name at the end via the output dataset
    out_lev_name = out_lev_name_native
    if out_lev_name == lev_name:
      out_lev_name = f'{lev_name}_target'
      ds_vert = ds_vert.rename({out_lev_name_native: out_lev_name})

    p_in  = _compute_input_pressure(ds_in, lev_name, ps)
    p_out = _compute_output_pressure(ds_vert, ps, out_lev_name)

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
      if verbose: print(f'    remapping {v}')
      ds_out[v] = _remap_field(ds_in[v], p_in, p_out, lev_name, out_lev_name, mode, extrap)

    # carry over global attrs and coords that aren't tied to the source vertical dim
    ds_out.attrs = dict(ds_in.attrs)

    # restore the native target lev dim name if we renamed it to avoid collision
    if out_lev_name != out_lev_name_native:
      ds_out = ds_out.rename({out_lev_name: out_lev_name_native})

    ds_out.to_netcdf(tmp_file, mode='w', format=nc_output_format)

  if same_in_out:
    os.replace(tmp_file, output_file)

  return
# -------------------------------------------------------------------------------------------------
