#!/usr/bin/env python
# ===================================================================================================
# Unit testing for hiccup_vertical_remap module
# ===================================================================================================
import os
import shutil
import tempfile
import unittest
import numpy as np
import xarray as xr
from time import perf_counter

from hiccup.hiccup_data_class_timer_methods import print_timer
from hiccup.hiccup_vertical_remap import (remap_vertical_py,
                                          _resolve_surface_pressure,
                                          _compute_input_pressure,
                                          _compute_output_pressure,
                                          _interp_column,
                                          _remap_field,
                                          _DEFAULT_P0)

# ---------------------------------------------------------------------------
# fixture helpers
# ---------------------------------------------------------------------------
def _make_hybrid_coords(nlev, p_top=100.0, p_bot=1.0e5, p0=1.0e5):
  """
  build a simple hybrid sigma-pressure grid with smoothly varying eta
  hyam/hybm at midpoints, hyai/hybi at interfaces; pure-pressure at top, pure-sigma at bottom
  returns (hyam, hybm, hyai, hybi) as numpy arrays
  """
  eta_i = np.linspace(p_top/p0, 1.0, nlev+1)
  # linear ramp: hybi grows from 0 at top to 1 at surface
  hybi = np.clip((eta_i - p_top/p0)/(1.0 - p_top/p0), 0.0, 1.0)
  hyai = eta_i - hybi
  hyam = 0.5*(hyai[:-1] + hyai[1:])
  hybm = 0.5*(hybi[:-1] + hybi[1:])
  return hyam, hybm, hyai, hybi

# ---------------------------------------------------------------------------
def _write_source_dataset(path, ncol=4, ntime=2, nlev_src=20,
                          ps_min=9.5e4, ps_max=1.015e5, lev_name='lev'):
  """
  write a small synthetic source dataset on a hybrid grid:
    - PS varies linearly across columns
    - T(p) = 250 + 30*log(p/1e5)   (smooth in log p)
    - Q(p) = exp(-p/1e5)           (smooth, monotonic)
  layout: (time, lev, ncol)
  """
  hyam, hybm, hyai, hybi = _make_hybrid_coords(nlev_src)
  ps = np.linspace(ps_min, ps_max, ncol)
  ps = np.broadcast_to(ps, (ntime, ncol)).copy()

  # pressure at midpoints: (lev, ncol)
  p = hyam[:, None]*1.0e5 + hybm[:, None]*ps[0, None, :]  # use first time slice for shape
  p = np.broadcast_to(p[None, :, :], (ntime, nlev_src, ncol)).copy()
  # recompute per-time using actual ps for accuracy
  p = (hyam[None, :, None]*1.0e5
       + hybm[None, :, None]*ps[:, None, :])

  T = 250.0 + 30.0*np.log(p/1.0e5)
  Q = np.exp(-p/1.0e5)
  TS = np.full((ntime, ncol), 290.0)  # passthrough scalar field

  ds = xr.Dataset(
    data_vars={
      'T':  ((('time', lev_name, 'ncol')), T.astype(np.float64)),
      'Q':  ((('time', lev_name, 'ncol')), Q.astype(np.float64)),
      'PS': ((('time', 'ncol')), ps.astype(np.float64)),
      'TS': ((('time', 'ncol')), TS.astype(np.float64)),
      'hyam': ((lev_name,), hyam),
      'hybm': ((lev_name,), hybm),
      'P0':   ((), np.float64(1.0e5)),
    },
    coords={
      'time': np.arange(ntime),
      lev_name: np.arange(nlev_src),
      'ncol': np.arange(ncol),
    },
  )
  ds.to_netcdf(path)
  return ds

# ---------------------------------------------------------------------------
def _write_vert_file(path, nlev_tgt=12, lev_name='lev'):
  """
  write a target vertical grid file with a coarser hybrid grid
  """
  hyam, hybm, hyai, hybi = _make_hybrid_coords(nlev_tgt)
  ds = xr.Dataset(
    data_vars={
      'hyam': ((lev_name,), hyam),
      'hybm': ((lev_name,), hybm),
      'hyai': (('ilev',), hyai),
      'hybi': (('ilev',), hybi),
      'P0':   ((), np.float64(1.0e5)),
    },
  )
  ds.to_netcdf(path)
  return ds

# ===========================================================================
class interp_column_test_case(unittest.TestCase):
  """
  Tests for the per-column interpolation kernel
  """
  # ------------------------------------------------------------------------
  def test_log_pressure_exact_for_linear_in_log_p(self):
    """
    a field that is exactly linear in log(p) should round-trip to machine precision
    """
    timer_start = perf_counter()
    p_src = np.linspace(1.0e3, 1.0e5, 30)
    p_tgt = np.linspace(2.0e3, 9.0e4, 12)
    f_src = 250.0 + 30.0*np.log(p_src/1.0e5)
    f_exp = 250.0 + 30.0*np.log(p_tgt/1.0e5)
    f_got = _interp_column(p_tgt, p_src, f_src, mode='log_pressure')
    np.testing.assert_allclose(f_got, f_exp, rtol=1e-12, atol=1e-12)
    print_timer(timer_start, caller='test_log_pressure_exact_for_linear_in_log_p')
  # ------------------------------------------------------------------------
  def test_linear_pressure_exact_for_linear_in_p(self):
    """
    a field that is exactly linear in p should round-trip under linear-in-pressure mode
    """
    timer_start = perf_counter()
    p_src = np.linspace(1.0e3, 1.0e5, 30)
    p_tgt = np.linspace(2.0e3, 9.0e4, 12)
    f_src = 100.0 + 1.5e-3*p_src
    f_exp = 100.0 + 1.5e-3*p_tgt
    f_got = _interp_column(p_tgt, p_src, f_src, mode='linear_pressure')
    np.testing.assert_allclose(f_got, f_exp, rtol=1e-12, atol=1e-12)
    print_timer(timer_start, caller='test_linear_pressure_exact_for_linear_in_p')
  # ------------------------------------------------------------------------
  def test_constant_extrapolation_clamps(self):
    """
    extrap='constant' should clamp output beyond the source range to source endpoints
    """
    timer_start = perf_counter()
    p_src = np.array([1.0e4, 5.0e4, 1.0e5])
    f_src = np.array([10.0, 20.0, 30.0])
    p_tgt = np.array([5.0e3, 1.5e5])  # both outside
    f_got = _interp_column(p_tgt, p_src, f_src, mode='linear_pressure', extrap='constant')
    self.assertAlmostEqual(f_got[0], 10.0)
    self.assertAlmostEqual(f_got[1], 30.0)
    print_timer(timer_start, caller='test_constant_extrapolation_clamps')
  # ------------------------------------------------------------------------
  def test_linear_extrapolation_uses_endpoint_slope(self):
    """
    extrap='linear' should extend with the slope of the last two source points
    """
    timer_start = perf_counter()
    p_src = np.array([1.0e4, 5.0e4, 1.0e5])
    f_src = np.array([10.0, 20.0, 30.0])  # slope = (30-20)/(1e5-5e4) = 2e-4
    p_tgt = np.array([1.5e5])
    f_got = _interp_column(p_tgt, p_src, f_src, mode='linear_pressure', extrap='linear')
    expected = 30.0 + (1.5e5 - 1.0e5)*2.0e-4
    self.assertAlmostEqual(f_got[0], expected)
    print_timer(timer_start, caller='test_linear_extrapolation_uses_endpoint_slope')
  # ------------------------------------------------------------------------
  def test_linear_extrapolation_falls_back_when_single_source_point(self):
    """
    a column with a single source level cannot define a slope; extrap='linear'
    should fall back to constant (clamp) instead of raising IndexError
    """
    timer_start = perf_counter()
    p_src = np.array([5.0e4])
    f_src = np.array([42.0])
    p_tgt = np.array([1.0e4, 5.0e4, 1.0e5])
    f_got = _interp_column(p_tgt, p_src, f_src, mode='linear_pressure', extrap='linear')
    np.testing.assert_array_equal(f_got, np.array([42.0, 42.0, 42.0]))
    print_timer(timer_start, caller='test_linear_extrapolation_falls_back_when_single_source_point')
  # ------------------------------------------------------------------------
  def test_handles_descending_source_pressure(self):
    """
    column ordered top-down (decreasing pressure) should give the same result as ascending
    """
    timer_start = perf_counter()
    p_src = np.linspace(1.0e3, 1.0e5, 30)
    f_src = 250.0 + 30.0*np.log(p_src/1.0e5)
    p_tgt = np.linspace(2.0e3, 9.0e4, 12)

    f_asc  = _interp_column(p_tgt, p_src, f_src, mode='log_pressure')
    f_desc = _interp_column(p_tgt, p_src[::-1], f_src[::-1], mode='log_pressure')
    np.testing.assert_allclose(f_asc, f_desc, rtol=1e-12, atol=1e-12)
    print_timer(timer_start, caller='test_handles_descending_source_pressure')
  # ------------------------------------------------------------------------

# ===========================================================================
class compute_pressure_test_case(unittest.TestCase):
  """
  Tests for _compute_input_pressure / _compute_output_pressure
  """
  # ------------------------------------------------------------------------
  def test_input_pressure_hybrid_eam(self):
    """
    EAM hybrid input (P0 present) should produce p = hyam*P0 + hybm*PS
    """
    timer_start = perf_counter()
    ds = xr.Dataset({
      'hyam': (('lev',), np.array([0.1, 0.2])),
      'hybm': (('lev',), np.array([0.0, 0.5])),
      'P0':   ((), np.float64(1.0e5)),
      'PS':   (('ncol',), np.array([1.0e5, 9.5e4])),
    })
    ps = ds['PS']
    p = _compute_input_pressure(ds, lev_name='lev', ps=ps).values
    expected = np.array([
      [0.1*1.0e5 + 0.0*1.0e5, 0.1*1.0e5 + 0.0*9.5e4],
      [0.2*1.0e5 + 0.5*1.0e5, 0.2*1.0e5 + 0.5*9.5e4],
    ])
    np.testing.assert_allclose(p, expected, rtol=1e-12)
    print_timer(timer_start, caller='test_input_pressure_hybrid_eam')
  # ------------------------------------------------------------------------
  def test_input_pressure_hybrid_eam_defaults_p0_when_missing(self):
    """
    EAM hybrid input with hyam/hybm but NO P0 (file not yet through
    add_reference_pressure) must still use the hybrid formula with a default P0.
    Regression: the old code fell through to the bare-lev fallback and used the
    hPa-valued lev coordinate as Pa, clamping everything below ~10 hPa to a
    constant (only visibly wrong in T, whose profile has structure aloft).
    """
    timer_start = perf_counter()
    ds = xr.Dataset({
      'hyam': (('lev',), np.array([0.1, 0.2])),
      'hybm': (('lev',), np.array([0.0, 0.5])),
      'PS':   (('ncol',), np.array([1.0e5, 9.5e4])),
    }, coords={'lev': np.array([100.0, 800.0])})  # hPa-looking lev, must be ignored
    ps = ds['PS']
    p = _compute_input_pressure(ds, lev_name='lev', ps=ps).values
    expected = np.array([
      [0.1*_DEFAULT_P0 + 0.0*1.0e5, 0.1*_DEFAULT_P0 + 0.0*9.5e4],
      [0.2*_DEFAULT_P0 + 0.5*1.0e5, 0.2*_DEFAULT_P0 + 0.5*9.5e4],
    ])
    np.testing.assert_allclose(p, expected, rtol=1e-12)
    print_timer(timer_start, caller='test_input_pressure_hybrid_eam_defaults_p0_when_missing')
  # ------------------------------------------------------------------------
  def test_input_pressure_pure_pressure_hpa_units_converted(self):
    """
    a pressure-level coordinate that advertises hPa/millibar units should be
    converted to Pa so it matches the Pa-valued target grid
    """
    timer_start = perf_counter()
    plev = np.array([10.0, 100.0, 500.0, 1000.0])  # hPa
    lev = xr.DataArray(plev, dims='plev'); lev.attrs['units'] = 'hPa'
    ds = xr.Dataset({'PS': (('ncol',), np.array([1.0e5]))},
                    coords={'plev': lev})
    p = _compute_input_pressure(ds, lev_name='plev', ps=ds['PS']).values
    np.testing.assert_allclose(p, plev*100.0, rtol=1e-12)
    print_timer(timer_start, caller='test_input_pressure_pure_pressure_hpa_units_converted')
  # ------------------------------------------------------------------------
  def test_input_pressure_hybrid_ifs(self):
    """
    IFS hybrid input (lnsp present, no P0) should produce p = hyam + hybm*ps
    with hyam in Pa
    """
    timer_start = perf_counter()
    hyam_pa = np.array([1.0e3, 2.0e3])
    ds = xr.Dataset({
      'hyam': (('lev',), hyam_pa),
      'hybm': (('lev',), np.array([0.0, 0.5])),
      'lnsp': (('ncol',), np.log(np.array([1.0e5, 9.5e4]))),
    })
    ps = _resolve_surface_pressure(ds, 'PS')
    p = _compute_input_pressure(ds, lev_name='lev', ps=ps).values
    expected = np.array([
      [1.0e3 + 0.0*1.0e5, 1.0e3 + 0.0*9.5e4],
      [2.0e3 + 0.5*1.0e5, 2.0e3 + 0.5*9.5e4],
    ])
    np.testing.assert_allclose(p, expected, rtol=1e-12)
    print_timer(timer_start, caller='test_input_pressure_hybrid_ifs')
  # ------------------------------------------------------------------------
  def test_input_pressure_hybrid_ifs_with_p0_present(self):
    """
    an IFS hybrid input (hyam in Pa) that also carries P0 - e.g. after
    add_reference_pressure - must still use the additive IFS formula
    p = hyam + hybm*ps and NOT scale hyam by P0
    """
    timer_start = perf_counter()
    hyam_pa = np.array([1.0e3, 2.0e3])
    ds = xr.Dataset({
      'hyam': (('lev',), hyam_pa),
      'hybm': (('lev',), np.array([0.0, 0.5])),
      'lnsp': (('ncol',), np.log(np.array([1.0e5, 9.5e4]))),
      'P0':   ((), np.float64(1.0e5)),
    })
    ps = _resolve_surface_pressure(ds, 'PS')
    p = _compute_input_pressure(ds, lev_name='lev', ps=ps).values
    expected = np.array([
      [1.0e3 + 0.0*1.0e5, 1.0e3 + 0.0*9.5e4],
      [2.0e3 + 0.5*1.0e5, 2.0e3 + 0.5*9.5e4],
    ])
    np.testing.assert_allclose(p, expected, rtol=1e-12)
    print_timer(timer_start, caller='test_input_pressure_hybrid_ifs_with_p0_present')
  # ------------------------------------------------------------------------
  def test_input_pressure_hybrid_ifs_with_ps_no_lnsp(self):
    """
    an IFS hybrid input (hyam in Pa) that supplies PS directly instead of lnsp
    must be recognized as IFS from the hyam magnitude, not mis-scaled as EAM
    """
    timer_start = perf_counter()
    hyam_pa = np.array([1.0e3, 2.0e3])
    ds = xr.Dataset({
      'hyam': (('lev',), hyam_pa),
      'hybm': (('lev',), np.array([0.0, 0.5])),
      'PS':   (('ncol',), np.array([1.0e5, 9.5e4])),
    })
    ps = _resolve_surface_pressure(ds, 'PS')
    p = _compute_input_pressure(ds, lev_name='lev', ps=ps).values
    expected = np.array([
      [1.0e3 + 0.0*1.0e5, 1.0e3 + 0.0*9.5e4],
      [2.0e3 + 0.5*1.0e5, 2.0e3 + 0.5*9.5e4],
    ])
    np.testing.assert_allclose(p, expected, rtol=1e-12)
    print_timer(timer_start, caller='test_input_pressure_hybrid_ifs_with_ps_no_lnsp')
  # ------------------------------------------------------------------------
  def test_input_pressure_pure_pressure_levels(self):
    """
    pure pressure-level input should return the lev coord directly
    """
    timer_start = perf_counter()
    plev = np.array([1.0e3, 5.0e3, 1.0e4, 5.0e4, 1.0e5])
    ds = xr.Dataset({'PS': (('ncol',), np.array([1.0e5]))},
                    coords={'plev': plev})
    ps = ds['PS']
    p = _compute_input_pressure(ds, lev_name='plev', ps=ps).values
    np.testing.assert_allclose(p, plev, rtol=1e-12)
    print_timer(timer_start, caller='test_input_pressure_pure_pressure_levels')
  # ------------------------------------------------------------------------
  def test_input_pressure_unknown_layout_raises(self):
    """
    a dataset with neither hybrid coefs nor a recognizable lev coord should raise ValueError
    """
    timer_start = perf_counter()
    ds = xr.Dataset({'PS': (('ncol',), np.array([1.0e5]))})
    with self.assertRaises(ValueError):
      _compute_input_pressure(ds, lev_name='lev', ps=ds['PS'])
    print_timer(timer_start, caller='test_input_pressure_unknown_layout_raises')
  # ------------------------------------------------------------------------
  def test_resolve_surface_pressure_prefers_named_ps(self):
    """
    when both ps_name and lnsp are present, the named variable should win
    """
    timer_start = perf_counter()
    ds = xr.Dataset({
      'PS':   (('ncol',), np.array([1.0e5, 9.5e4])),
      'lnsp': (('ncol',), np.log(np.array([8.0e4, 7.0e4]))),
    })
    ps = _resolve_surface_pressure(ds, 'PS')
    np.testing.assert_array_equal(ps.values, np.array([1.0e5, 9.5e4]))
    print_timer(timer_start, caller='test_resolve_surface_pressure_prefers_named_ps')
  # ------------------------------------------------------------------------
  def test_resolve_surface_pressure_falls_back_to_lnsp(self):
    """
    when only lnsp is present, surface pressure should be derived as exp(lnsp)
    and named ps_name
    """
    timer_start = perf_counter()
    ps_true = np.array([1.0e5, 9.5e4])
    ds = xr.Dataset({'lnsp': (('ncol',), np.log(ps_true))})
    ps = _resolve_surface_pressure(ds, 'PS')
    self.assertEqual(ps.name, 'PS')
    np.testing.assert_allclose(ps.values, ps_true, rtol=1e-12)
    print_timer(timer_start, caller='test_resolve_surface_pressure_falls_back_to_lnsp')
  # ------------------------------------------------------------------------
  def test_resolve_surface_pressure_missing_raises(self):
    """
    a dataset with neither ps_name nor lnsp should raise ValueError
    """
    timer_start = perf_counter()
    ds = xr.Dataset({'foo': (('x',), np.array([1.0]))})
    with self.assertRaises(ValueError):
      _resolve_surface_pressure(ds, 'PS')
    print_timer(timer_start, caller='test_resolve_surface_pressure_missing_raises')
  # ------------------------------------------------------------------------
  def test_resolve_surface_pressure_preserves_singleton_time(self):
    """
    a single-timestep lnsp with a singleton vertical dim should have the lev dim
    dropped but the time dim preserved (regression: blanket squeeze() would drop both)
    """
    timer_start = perf_counter()
    ps_true = np.array([[1.0e5, 9.5e4]])  # shape (time=1, ncol=2)
    lnsp = np.log(ps_true)[:, None, :]    # shape (time=1, lev=1, ncol=2)
    ds = xr.Dataset({'lnsp': (('time', 'lev', 'ncol'), lnsp)})
    ps = _resolve_surface_pressure(ds, 'PS', lev_name='lev')
    self.assertEqual(ps.dims, ('time', 'ncol'))
    self.assertEqual(ps.shape, (1, 2))
    np.testing.assert_allclose(ps.values, ps_true, rtol=1e-12)
    print_timer(timer_start, caller='test_resolve_surface_pressure_preserves_singleton_time')
  # ------------------------------------------------------------------------
  def test_resolve_surface_pressure_lnsp_without_lev_dim(self):
    """
    when lnsp has no vertical dim at all, it should be used as-is (no isel)
    """
    timer_start = perf_counter()
    ps_true = np.array([1.0e5, 9.5e4])
    ds = xr.Dataset({'lnsp': (('ncol',), np.log(ps_true))})
    ps = _resolve_surface_pressure(ds, 'PS', lev_name='lev')
    self.assertEqual(ps.dims, ('ncol',))
    np.testing.assert_allclose(ps.values, ps_true, rtol=1e-12)
    print_timer(timer_start, caller='test_resolve_surface_pressure_lnsp_without_lev_dim')
  # ------------------------------------------------------------------------
  def test_output_pressure_uses_target_hybrid(self):
    """
    target pressure should be hyam*P0 + hybm*PS using the target file's coefs
    """
    timer_start = perf_counter()
    ds_vert = xr.Dataset({
      'hyam': (('lev',), np.array([0.05, 0.3])),
      'hybm': (('lev',), np.array([0.0, 0.7])),
      'P0':   ((), np.float64(1.0e5)),
    })
    ps = xr.DataArray(np.array([1.0e5, 9.0e4]), dims=('ncol',))
    p = _compute_output_pressure(ds_vert, ps, out_lev_name='lev').values
    expected = np.array([
      [0.05*1.0e5, 0.05*1.0e5],
      [0.3*1.0e5 + 0.7*1.0e5, 0.3*1.0e5 + 0.7*9.0e4],
    ])
    np.testing.assert_allclose(p, expected, rtol=1e-12)
    print_timer(timer_start, caller='test_output_pressure_uses_target_hybrid')
  # ------------------------------------------------------------------------
  def test_output_pressure_missing_hybrid_raises(self):
    """
    target file without hyam/hybm should raise ValueError
    """
    timer_start = perf_counter()
    ds_vert = xr.Dataset({'foo': (('x',), np.array([1.0]))})
    ps = xr.DataArray(np.array([1.0e5]), dims=('ncol',))
    with self.assertRaises(ValueError):
      _compute_output_pressure(ds_vert, ps, out_lev_name='lev')
    print_timer(timer_start, caller='test_output_pressure_missing_hybrid_raises')
  # ------------------------------------------------------------------------

# ===========================================================================
class remap_vertical_end_to_end_test_case(unittest.TestCase):
  """
  End-to-end tests that write fixtures, call remap_vertical_py(), and inspect the output
  """
  # ------------------------------------------------------------------------
  def setUp(self):
    self.tmpdir = tempfile.mkdtemp(prefix='hiccup_vrt_test_')
    self.src_file  = os.path.join(self.tmpdir, 'src.nc')
    self.vert_file = os.path.join(self.tmpdir, 'vert.nc')
    self.out_file  = os.path.join(self.tmpdir, 'out.nc')
    self.ds_src  = _write_source_dataset(self.src_file)
    self.ds_vert = _write_vert_file(self.vert_file)
  # ------------------------------------------------------------------------
  def tearDown(self):
    shutil.rmtree(self.tmpdir, ignore_errors=True)
  # ------------------------------------------------------------------------
  def test_smooth_field_remap_is_accurate(self):
    """
    a field that is exactly linear in log(p) should remap to the same analytic field
    on the target grid, to within float64 precision (down to round-off in p computation)
    """
    timer_start = perf_counter()
    remap_vertical_py(self.src_file, self.out_file, self.vert_file, ps_name='PS', lev_name='lev')
    with xr.open_dataset(self.out_file) as ds_out:
      hyam = ds_out['hyam'].values
      hybm = ds_out['hybm'].values
      ps   = ds_out['PS'].values  # (time, ncol)
      p_tgt = hyam[None, :, None]*1.0e5 + hybm[None, :, None]*ps[:, None, :]
      T_exp = 250.0 + 30.0*np.log(p_tgt/1.0e5)
      T_got = ds_out['T'].transpose('time', 'lev', 'ncol').values
      np.testing.assert_allclose(T_got, T_exp, rtol=1e-10, atol=1e-10)
    print_timer(timer_start, caller='test_smooth_field_remap_is_accurate')
  # ------------------------------------------------------------------------
  def test_dimension_order_is_preserved(self):
    """
    the remapped fields must keep the source dimension order (time,lev,ncol).
    xarray.apply_ufunc moves core dims to the last axis, which silently produced
    (time,ncol,lev) output - EAM reads initial conditions positionally and does not
    validate the layout, so a transposed IC blows up in the dycore on the first step
    instead of raising an error. Assert the layout explicitly rather than calling
    .transpose() first, which is what let this regression through before.
    """
    timer_start = perf_counter()
    remap_vertical_py(self.src_file, self.out_file, self.vert_file, ps_name='PS', lev_name='lev')
    with xr.open_dataset(self.out_file) as ds_out:
      for var in ['T','Q']:
        self.assertEqual(ds_out[var].dims, ('time','lev','ncol'),
                         msg=f'{var} dims {ds_out[var].dims} != source order (time,lev,ncol)')
      self.assertEqual(ds_out['PS'].dims, ('time','ncol'))
    print_timer(timer_start, caller='test_dimension_order_is_preserved')
  # ------------------------------------------------------------------------
  def test_dimension_order_is_preserved_for_transposed_source(self):
    """
    a source written as (time,ncol,lev) must come back as (time,ncol,lev) - the remap
    should mirror whatever layout it was handed, not impose one of its own
    """
    timer_start = perf_counter()
    src_file = os.path.join(self.tmpdir, 'src_transposed.nc')
    self.ds_src.transpose('time','ncol','lev').to_netcdf(src_file)
    remap_vertical_py(src_file, self.out_file, self.vert_file, ps_name='PS', lev_name='lev')
    with xr.open_dataset(self.out_file) as ds_out:
      for var in ['T','Q']:
        self.assertEqual(ds_out[var].dims, ('time','ncol','lev'),
                         msg=f'{var} dims {ds_out[var].dims} != source order (time,ncol,lev)')
    print_timer(timer_start, caller='test_dimension_order_is_preserved_for_transposed_source')
  # ------------------------------------------------------------------------
  def test_passthrough_variables_unchanged(self):
    """
    variables without the source lev dim (PS, TS) should be copied through bit-identical
    """
    timer_start = perf_counter()
    remap_vertical_py(self.src_file, self.out_file, self.vert_file, ps_name='PS', lev_name='lev')
    with xr.open_dataset(self.out_file) as ds_out:
      np.testing.assert_array_equal(ds_out['PS'].values, self.ds_src['PS'].values)
      np.testing.assert_array_equal(ds_out['TS'].values, self.ds_src['TS'].values)
    print_timer(timer_start, caller='test_passthrough_variables_unchanged')
  # ------------------------------------------------------------------------
  def test_target_hybrid_coefs_carried_over(self):
    """
    hyam/hybm/hyai/hybi/P0 from the target vert file should appear in the output
    """
    timer_start = perf_counter()
    remap_vertical_py(self.src_file, self.out_file, self.vert_file, ps_name='PS', lev_name='lev')
    with xr.open_dataset(self.out_file) as ds_out:
      np.testing.assert_array_equal(ds_out['hyam'].values, self.ds_vert['hyam'].values)
      np.testing.assert_array_equal(ds_out['hybm'].values, self.ds_vert['hybm'].values)
      np.testing.assert_array_equal(ds_out['hyai'].values, self.ds_vert['hyai'].values)
      np.testing.assert_array_equal(ds_out['hybi'].values, self.ds_vert['hybi'].values)
    print_timer(timer_start, caller='test_target_hybrid_coefs_carried_over')
  # ------------------------------------------------------------------------
  def test_explicit_chunks_produce_correct_output(self):
    """
    explicit horizontal chunking should still produce correct output;
    proves the dask-parallelized interp path works end-to-end
    """
    timer_start = perf_counter()
    remap_vertical_py(self.src_file, self.out_file, self.vert_file, ps_name='PS',
                      lev_name='lev', chunks={'ncol': 2, 'time': 1})
    with xr.open_dataset(self.out_file) as ds_out:
      hyam = ds_out['hyam'].values
      hybm = ds_out['hybm'].values
      ps   = ds_out['PS'].values
      p_tgt = hyam[None, :, None]*1.0e5 + hybm[None, :, None]*ps[:, None, :]
      T_exp = 250.0 + 30.0*np.log(p_tgt/1.0e5)
      T_got = ds_out['T'].transpose('time', 'lev', 'ncol').values
      np.testing.assert_allclose(T_got, T_exp, rtol=1e-10, atol=1e-10)
    print_timer(timer_start, caller='test_explicit_chunks_produce_correct_output')
  # ------------------------------------------------------------------------
  def test_default_chunks_chunk_horizontal_dims(self):
    """
    when chunks=None, the input dataset should be opened with dask backing on
    every dim (lev_name forced to -1) - regression against an earlier default
    that left non-lev dims as a single chunk
    """
    timer_start = perf_counter()
    # build a wider source so 'auto' has something to bite into
    wide_src = os.path.join(self.tmpdir, 'wide_src.nc')
    _write_source_dataset(wide_src, ncol=64, ntime=4, nlev_src=20)
    out = os.path.join(self.tmpdir, 'wide_out.nc')
    remap_vertical_py(wide_src, out, self.vert_file, ps_name='PS', lev_name='lev')
    # confirm correctness on the wider grid
    with xr.open_dataset(out) as ds_out:
      hyam = ds_out['hyam'].values
      hybm = ds_out['hybm'].values
      ps   = ds_out['PS'].values
      p_tgt = hyam[None, :, None]*1.0e5 + hybm[None, :, None]*ps[:, None, :]
      T_exp = 250.0 + 30.0*np.log(p_tgt/1.0e5)
      T_got = ds_out['T'].transpose('time', 'lev', 'ncol').values
      np.testing.assert_allclose(T_got, T_exp, rtol=1e-10, atol=1e-10)
    print_timer(timer_start, caller='test_default_chunks_chunk_horizontal_dims')
  # ------------------------------------------------------------------------
  def test_p0_always_written_even_when_vert_file_lacks_it(self):
    """
    when vert_file has no P0, the output should still carry P0 with the default
    used during pressure computation, so the output file is self-describing
    """
    timer_start = perf_counter()
    no_p0_vert = os.path.join(self.tmpdir, 'vert_no_p0.nc')
    self.ds_vert.drop_vars('P0').to_netcdf(no_p0_vert)
    out = os.path.join(self.tmpdir, 'out_no_p0.nc')
    remap_vertical_py(self.src_file, out, no_p0_vert, ps_name='PS', lev_name='lev')
    with xr.open_dataset(out) as ds_out:
      self.assertIn('P0', ds_out.variables)
      self.assertEqual(float(ds_out['P0'].values), 1.0e5)
    print_timer(timer_start, caller='test_p0_always_written_even_when_vert_file_lacks_it')
  # ------------------------------------------------------------------------
  def test_source_missing_p0_remaps_correctly_not_constant(self):
    """
    a hybrid source that still has hyam/hybm but lost P0 must remap to the correct
    analytic profile - NOT collapse to a constant below ~10 hPa. This is the
    end-to-end guard for the reported bug where only T looked wrong because the
    source pressure silently fell back to the hPa lev coordinate.
    """
    timer_start = perf_counter()
    no_p0_src = os.path.join(self.tmpdir, 'src_no_p0.nc')
    # give lev genuine hPa-looking values so the old fallback would misfire
    hyam = self.ds_src['hyam'].values
    hybm = self.ds_src['hybm'].values
    lev_hpa = (hyam*1.0e5 + hybm*1.0e5)/100.0
    src = self.ds_src.drop_vars('P0').assign_coords(lev=lev_hpa)
    src.to_netcdf(no_p0_src)
    out = os.path.join(self.tmpdir, 'out_no_p0.nc')
    remap_vertical_py(no_p0_src, out, self.vert_file, ps_name='PS', lev_name='lev')
    with xr.open_dataset(out) as ds_out:
      hyam_t = ds_out['hyam'].values
      hybm_t = ds_out['hybm'].values
      ps     = ds_out['PS'].values
      p_tgt  = hyam_t[None, :, None]*1.0e5 + hybm_t[None, :, None]*ps[:, None, :]
      T_exp  = 250.0 + 30.0*np.log(p_tgt/1.0e5)
      T_got  = ds_out['T'].transpose('time', 'lev', 'ncol').values
      np.testing.assert_allclose(T_got, T_exp, rtol=1e-8, atol=1e-8)
      # and it must not be constant in the vertical
      self.assertGreater(np.ptp(T_got, axis=1).min(), 1.0)
    print_timer(timer_start, caller='test_source_missing_p0_remaps_correctly_not_constant')
  # ------------------------------------------------------------------------
  def test_hpa_pa_unit_mismatch_raises(self):
    """
    a source with no hybrid coefs and a bare hPa lev coordinate (no units attr)
    remapped onto a Pa target grid should raise a clear ValueError rather than
    silently clamping everything below ~10 hPa to a constant
    """
    timer_start = perf_counter()
    hpa_src = os.path.join(self.tmpdir, 'src_hpa.nc')
    lev_hpa = np.linspace(1.0, 1000.0, 20)  # hPa, ascending; no units attribute
    ntime, ncol = 2, 4
    T = 250.0 + 30.0*np.log((lev_hpa*100.0)/1.0e5)
    T = np.broadcast_to(T[None, :, None], (ntime, 20, ncol))
    ds = xr.Dataset(
      data_vars={
        'T':  (('time', 'lev', 'ncol'), T.astype(np.float64)),
        'PS': (('time', 'ncol'), np.full((ntime, ncol), 1.0e5)),
      },
      coords={'time': np.arange(ntime), 'lev': lev_hpa, 'ncol': np.arange(ncol)},
    )
    ds.to_netcdf(hpa_src)
    with self.assertRaises(ValueError):
      remap_vertical_py(hpa_src, self.out_file, self.vert_file, ps_name='PS', lev_name='lev')
    print_timer(timer_start, caller='test_hpa_pa_unit_mismatch_raises')
  # ------------------------------------------------------------------------
  def test_var_list_filters_remapped_fields(self):
    """
    explicit var_list should remap only the listed fields; other lev-bearing vars are dropped
    """
    timer_start = perf_counter()
    remap_vertical_py(self.src_file, self.out_file, self.vert_file, ps_name='PS',
                   var_list=['T'], lev_name='lev')
    with xr.open_dataset(self.out_file) as ds_out:
      self.assertIn('T', ds_out.data_vars)
      self.assertNotIn('Q', ds_out.data_vars)
    print_timer(timer_start, caller='test_var_list_filters_remapped_fields')
  # ------------------------------------------------------------------------
  def test_invalid_mode_raises(self):
    """
    invalid interpolation mode should raise ValueError before opening any files
    """
    timer_start = perf_counter()
    with self.assertRaises(ValueError):
      remap_vertical_py(self.src_file, self.out_file, self.vert_file, mode='cubic')
    print_timer(timer_start, caller='test_invalid_mode_raises')
  # ------------------------------------------------------------------------
  def test_invalid_extrap_raises(self):
    """
    invalid extrapolation mode should raise ValueError before opening any files
    """
    timer_start = perf_counter()
    with self.assertRaises(ValueError):
      remap_vertical_py(self.src_file, self.out_file, self.vert_file, extrap='spline')
    print_timer(timer_start, caller='test_invalid_extrap_raises')
  # ------------------------------------------------------------------------
  def test_input_equals_output_overwrites_in_place(self):
    """
    when input_file == output_file, the source file should be overwritten with remapped data
    """
    timer_start = perf_counter()
    in_place = os.path.join(self.tmpdir, 'in_place.nc')
    _write_source_dataset(in_place)
    remap_vertical_py(in_place, in_place, self.vert_file, ps_name='PS', lev_name='lev')
    with xr.open_dataset(in_place) as ds:
      self.assertEqual(ds['T'].sizes['lev'], self.ds_vert['hyam'].sizes['lev'])
    print_timer(timer_start, caller='test_input_equals_output_overwrites_in_place')
  # ------------------------------------------------------------------------
  def test_missing_ps_raises(self):
    """
    a source file lacking BOTH the named surface pressure variable and lnsp
    should raise ValueError
    """
    timer_start = perf_counter()
    no_ps = os.path.join(self.tmpdir, 'no_ps.nc')
    # the synthetic source has no lnsp, so dropping PS leaves no surface-pressure source
    self.ds_src.drop_vars('PS').to_netcdf(no_ps)
    with self.assertRaises(ValueError):
      remap_vertical_py(no_ps, self.out_file, self.vert_file, ps_name='PS', lev_name='lev')
    print_timer(timer_start, caller='test_missing_ps_raises')
  # ------------------------------------------------------------------------
  def test_lnsp_source_remaps_and_writes_ps(self):
    """
    an IFS-style source (hyam in Pa, lnsp instead of PS, no P0) should remap
    successfully and the output should carry PS = exp(lnsp)
    """
    timer_start = perf_counter()
    ifs_src = os.path.join(self.tmpdir, 'ifs_src.nc')
    ntime, ncol, nlev = 2, 4, 20
    # build IFS-style hybrid: hyam in Pa, no P0
    eta_i = np.linspace(100.0/1.0e5, 1.0, nlev+1)
    hybi  = np.clip((eta_i - 100.0/1.0e5)/(1.0 - 100.0/1.0e5), 0.0, 1.0)
    hyai  = (eta_i - hybi)*1.0e5  # convert to Pa
    hyam  = 0.5*(hyai[:-1] + hyai[1:])
    hybm  = 0.5*(hybi[:-1] + hybi[1:])
    ps    = np.broadcast_to(np.linspace(9.5e4, 1.015e5, ncol), (ntime, ncol)).copy()
    lnsp  = np.log(ps)
    p     = hyam[None, :, None] + hybm[None, :, None]*ps[:, None, :]
    T     = 250.0 + 30.0*np.log(p/1.0e5)
    ds = xr.Dataset(
      data_vars={
        'T':    (('time', 'lev', 'ncol'), T.astype(np.float64)),
        'lnsp': (('time', 'ncol'), lnsp.astype(np.float64)),
        'hyam': (('lev',), hyam),
        'hybm': (('lev',), hybm),
      },
      coords={'time': np.arange(ntime), 'lev': np.arange(nlev), 'ncol': np.arange(ncol)},
    )
    ds.to_netcdf(ifs_src)
    remap_vertical_py(ifs_src, self.out_file, self.vert_file, ps_name='PS', lev_name='lev')
    with xr.open_dataset(self.out_file) as ds_out:
      self.assertIn('PS', ds_out.data_vars)
      np.testing.assert_allclose(ds_out['PS'].values, ps, rtol=1e-12)
      # smooth field should round-trip to within float64 round-off
      hyam_t = ds_out['hyam'].values
      hybm_t = ds_out['hybm'].values
      ps_t   = ds_out['PS'].values
      p_tgt  = hyam_t[None, :, None]*1.0e5 + hybm_t[None, :, None]*ps_t[:, None, :]
      T_exp  = 250.0 + 30.0*np.log(p_tgt/1.0e5)
      T_got  = ds_out['T'].transpose('time', 'lev', 'ncol').values
      np.testing.assert_allclose(T_got, T_exp, rtol=1e-10, atol=1e-10)
    print_timer(timer_start, caller='test_lnsp_source_remaps_and_writes_ps')
  # ------------------------------------------------------------------------

# ===========================================================================
if __name__ == '__main__':
  unittest.main()
