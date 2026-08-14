#!/usr/bin/env python
#===================================================================================================
# Unit testing for state_adjustment module
#===================================================================================================
import unittest, tempfile, numpy as np, xarray as xr
from time import perf_counter
from hiccup.hiccup_data_class_timer_methods import print_timer
from hiccup import hiccup_state_adjustment as hsa
from hiccup.standard_atmosphere import standard_atmosphere
from hiccup.hiccup_constants import Rdair
from hiccup.hiccup_constants import Rvapor

verbose_default = False # local verbosity default

#===============================================================================
def remove_supersaturation_test( ds, hybrid_lev=False, pressure_var_name='plev',
                            debug=False, verbose=None, verbose_indent='' ):
  """
  Adjust the surface temperature based on new surace height assumed lapse rate 
    ncol            # columns
    qv              specific humidity
    temperature     temperature at layer mid-points [k]
    pressure        pressure at layer mid-points    (convert to hPa for qv_sat calculation)
  """
  if verbose is None : verbose = verbose_default
  if verbose: print(f'\n{verbose_indent}Removing super saturated data points...')
  if debug: print(f'{verbose_indent}remove_supersaturation: DEBUG MODE ENABLED')

  qv_min = 1.0e-9   # minimum specific humidity value allowed

  if hybrid_lev :
    pressure = get_pressure_from_hybrid(ds)/1e2
  else :
    pressure = ds[pressure_var_name]

  if debug:
    print(); print_stat(pressure,name='pressure in remove_supersaturation')
    print(); print_stat(ds['Q'],name='qv in remove_supersaturation')
    print(); print_stat(ds['T'],name='T in remove_supersaturation')

  # Calculate saturation specific humidity
  qv_sat = calculate_qv_sat_liq(ds['T'],pressure)
  
  if debug:
    print(); print_stat(qv_sat,name='qv_sat in remove_supersaturation')

  # The following check is to avoid the generation of negative values
  # that can occur in the upper stratosphere and mesosphere
  # qv_sat.values = xr.where(qv_sat.values>=0.0,qv_sat,1.0)

  qv_sat = xr.where(qv_sat>=0.0,qv_sat,1.0)

  # Calculate relative humidity for limiter
  rh = ds['Q'] / qv_sat

  if debug:
    print(); print_stat(rh,name='rh in remove_supersaturation')

  # save attributes to restore later
  tmp_attrs = ds['Q'].attrs

  # Apply limiter conditions
  ds['Q'] = xr.where(rh>1.,qv_sat,ds['Q'])
  ds['Q'] = xr.where(rh<0.,qv_min,ds['Q'])
  
  # restore attributes
  ds['Q'].attrs = tmp_attrs

  if debug:
    print(); print_stat(ds['Q'],name='qv in remove_supersaturation after adjustment')

  return

#===============================================================================
def calculate_qv_sat_liq( temperature, pressure ):
  """ 
  calculate saturation specific humidity [kg/kg]
  from temperature [K] and pressure [hPa]
  """

  # Calculate saturation vapor pressure [hPa] over liquid 
  # Bolton, D., 1980: The Computation of Equivalent Potential Temperature, MWR, 108, 1046-1053
  # https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2
  es = 6.112 * np.exp( 17.67*(temperature-273.0)/(temperature-273.0+243.5) ) 

  # Convert to mixing ratio
  r_sat = (Rdair/Rvapor) * es / (pressure - es)

  # Convert mixing ratio to saturation specific humidity
  qv_sat = r_sat / ( 1.0 + r_sat )

  return qv_sat

#===============================================================================
class state_adjustment_test_case(unittest.TestCase):
  """ 
  Tests for hiccup_state_adjustment.py 
  """
  # ----------------------------------------------------------------------------
  def test_adjust_surface_pressure(self):
    """ 
    Does surface pressure interpolation give the right value? 
    """
    timer_start = perf_counter()

    phis_int        = np.array([ 3.5e3, 2.5e3, 1.5e3, 0.5e3, 0.2e3, 0.1e3 ])   # interface altitudes to define plev
    phis_min        = np.min(phis_int)
    phis_new        = np.array([ phis_min-0.5e3, phis_min+0.5e3, phis_min ])
    plev            = len(phis_int)-1
    ilev            = plev+1
    ncol            = len(phis_new)
    isa_int         = standard_atmosphere( phis_int )
    temperature_mid = np.empty((ncol,plev))  
    pressure_mid    = np.empty((ncol,plev))  
    temperature_int = np.empty((ncol,plev+1))
    pressure_int    = np.empty((ncol,plev+1))
    for i in range(ncol):
      for k in range(plev+1):    
        temperature_int[i,k] = isa_int.temperature[k]
        pressure_int[i,k]    = isa_int.pressure[k]
        if k<plev: temperature_mid[i,k] = np.mean( temperature_int[i,k:k+1] )
        if k<plev: pressure_mid[i,k]    = np.mean( pressure_int[i,k:k+1] )
    ps_old        = np.array(pressure_int[:,-1])
    ts_old        = np.array(temperature_int[:,-1])
    phis_old      = np.array([phis_min]*ncol)

    pressure_mid = pressure_mid[0,:]

    # Convert into dataset - then add dummy time dimension 
    # use .copy() to avoid creating a read-only dataset
    ds_data = xr.Dataset({'PHIS':xr.DataArray(phis_old,dims=['ncol'])
                         ,'PS'  :xr.DataArray(ps_old,dims=['ncol'])
                         ,'TS'  :xr.DataArray(ts_old,dims=['ncol'])
                         ,'T'   :xr.DataArray(temperature_mid,dims=['ncol','lev'])
                         # ,'plev':xr.DataArray(pressure_mid,dims=['ncol','lev'])
                         ,'plev':xr.DataArray(pressure_mid,dims=['lev'])
                         # ,'Pint':xr.DataArray(pressure_int,dims=['ncol','lev'])
                         },
                         coords={'ncol':np.arange(ncol)
                                 ,'lev':np.arange(plev)
                                 ,'ilev':np.arange(ilev)}
                        ).expand_dims(time=1,axis=0).copy(deep=True)
    ds_data = ds_data.chunk(2)
    ds_topo = xr.Dataset({'PHIS':(['ncol'],phis_new)})

    ds_data = hsa.adjust_surface_pressure( ds_data, ds_topo, debug=False )

    # Get the adjusted surface pressure for value checking
    ps_new = ds_data['PS'].isel(time=0).values

    # for k in range(plev+1): print(f'  {k}  {pressure_int[0,k]:8.2f}  {phis_int[k]:8.2f} ')
    # for i in range(ncol): print(f'phis_old: {phis_old[i]:04.0f}  phis_new: {phis_new[i]:04.0f}  ps_old: {ps_old[i]:6.2f}  ps_new: {ps_new[i]:6.2f}')
    self.assertTrue( ps_new[0] >ps_old[0] )
    self.assertTrue( ps_new[1] <ps_old[1] )
    self.assertTrue( ps_new[2]==ps_old[2] )

    print_timer(timer_start,caller='test_adjust_surface_pressure')
  # ----------------------------------------------------------------------------
  def test_adjust_surface_temperature(self):
    """ 
    Does surface temperature interpolation give the right value? 
    """
    timer_start = perf_counter()

    phis      = 0.5e3
    phis_old  = np.array([phis]*3)
    phis_new  = np.array([ phis-0.5e3, phis+0.5e3, phis ])
    isa       = standard_atmosphere( phis_old )
    ts_old    = isa.temperature

    ds_data = xr.Dataset({'PHIS':(['ncol'],phis_old)
                         ,'TS'  :(['ncol'],ts_old)})
    ds_topo = xr.Dataset({'PHIS':(['ncol'],phis_new)})

    ds_data = hsa.adjust_surface_temperature( ds_data, ds_topo )

    ts_new = ds_data['TS'].values

    # for i in range(ncol): print(f'phis_old: {phis_old[i]:04.0f}  phis_new: {phis_new[i]:04.0f}  ts_old: {ts_old[i]:6.2f}  ts_new: {ts_new[i]:6.2f}')
    self.assertTrue( ts_new[0] <ts_old[0] )
    self.assertTrue( ts_new[1] >ts_old[1] )
    self.assertTrue( ts_new[2]==ts_old[2] )

    print_timer(timer_start,caller='test_adjust_surface_temperature')
  # ----------------------------------------------------------------------------
  def test_remove_supersaturation(self):
    """ 
    Do supersaturated values get limited correctly? 
    """
    timer_start = perf_counter()

    temperature_in = 300
    pressure_in    = 1010
    qv_sat = hsa.calculate_qv_sat_liq(temperature_in,pressure_in)
    qv = xr.DataArray(np.array([ 1.1*qv_sat , 1.0*qv_sat , 0.9*qv_sat ]))
    ncol = len(qv.values)
    temperature = xr.DataArray([temperature_in]*ncol)
    pressure    = xr.DataArray([pressure_in]   *ncol)

    # Convert into dataset
    ds = xr.Dataset({'Q'  :xr.DataArray(qv,dims=['ncol'])
                    ,'T'  :xr.DataArray(temperature,dims=['ncol'])
                    ,'plev':xr.DataArray(pressure,dims=['ncol'])
                    }, coords={'ncol':np.arange(ncol)} )

    # hsa.remove_supersaturation( ds )
    remove_supersaturation_test(ds)

    rh_out = ds['Q'].values/qv_sat
    expected_answer = np.array([1.0, 1.0, 0.9])
    self.assertTrue( np.all( np.abs(rh_out-expected_answer)<1e-10 ) )

    print_timer(timer_start,caller='test_remove_supersaturation')
  # ----------------------------------------------------------------------------
  def test_remove_supersaturation2(self):
    """ 
    Do supersaturated values get limited correctly? 
    """
    timer_start = perf_counter()

    temperature_in = 300
    pressure_in    = 1010
    qv_sat = hsa.calculate_qv_sat_liq(temperature_in,pressure_in)
    # qv = xr.DataArray(np.array([ 1.1*qv_sat , 1.0*qv_sat , 0.9*qv_sat ]))

    qv = xr.DataArray(np.random.uniform(0.9,1.1,int(1e8)))

    ncol = len(qv.values)
    temperature = xr.DataArray([temperature_in]*ncol)
    pressure    = xr.DataArray([pressure_in]   *ncol)

    # Convert into dataset
    ds = xr.Dataset({'Q'  :xr.DataArray(qv,dims=['ncol'])
                    ,'T'  :xr.DataArray(temperature,dims=['ncol'])
                    ,'plev':xr.DataArray(pressure,dims=['ncol'])
                    }, coords={'ncol':np.arange(ncol)} )

    # hsa.remove_supersaturation( ds )
    remove_supersaturation_test(ds)

    rh_out = ds['Q'].values/qv_sat
    # expected_answer = np.array([1.0, 1.0, 0.9])
    # self.assertTrue( np.all( np.abs(rh_out-expected_answer)<1e-10 ) )

    print_timer(timer_start,caller='test_remove_supersaturation2')
  
  # ----------------------------------------------------------------------------
  def test_adjust_cld_wtr(self):
    """
    Does cloud water and ice get limited properly?
    """
    timer_start = perf_counter()

    vals = [0,-1,1]
    expected_answer = [0,0,1]
    ncol = np.arange(len(vals))
    ds = xr.Dataset({'CLDLIQ':xr.DataArray(vals,dims=['ncol'])
                    ,'CLDICE':xr.DataArray(vals,dims=['ncol'])
                    }, coords={'ncol':ncol} )

    ds = hsa.adjust_cld_wtr( ds )

    self.assertTrue( np.all( np.abs(ds['CLDLIQ'].values-expected_answer)<1e-10 ) )
    self.assertTrue( np.all( np.abs(ds['CLDICE'].values-expected_answer)<1e-10 ) )

    print_timer(timer_start,caller='test_adjust_cld_wtr')
  # ----------------------------------------------------------------------------
  def test_adjust_cloud_fraction(self):
    """
    Does cloud fraction get limited properly?
    """
    timer_start = perf_counter()

    input_vals = [0,-1,1,2]
    expected_answer = [0,0,1,1]
    ncol = np.arange(len(input_vals))
    ds = xr.Dataset({'CLOUD_FRAC':xr.DataArray(input_vals,dims=['ncol'])
                    }, coords={'ncol':ncol} )

    ds = hsa.adjust_cloud_fraction( ds, frac_var_name='CLOUD_FRAC' )
    
    self.assertTrue( np.all( np.abs(ds['CLOUD_FRAC'].values-expected_answer)<1e-10 ) )

    print_timer(timer_start,caller='test_adjust_cloud_fraction')
  # ----------------------------------------------------------------------------
  def test_build_gaussian_smoother(self):
    """
    Does the Gaussian smoother produce a row-normalized operator?
    """
    timer_start = perf_counter()

    lon2d,lat2d = np.meshgrid( np.linspace(0,360,48,endpoint=False)
                             , np.linspace(-80,80,24) )
    lat = lat2d.ravel(); lon = lon2d.ravel()

    S = hsa.build_gaussian_smoother( lat, lon, corr_length_km=2000. )

    # rows should sum to 1 (partition of unity) so the filter preserves the mean
    row_sum = np.asarray(S.sum(axis=1)).ravel()
    self.assertTrue( np.all( np.abs(row_sum-1.)<1e-10 ) )
    # a constant field should be unchanged by the smoother
    const = np.full(lat.size,3.7)
    self.assertTrue( np.all( np.abs((S@const)-3.7)<1e-10 ) )

    print_timer(timer_start,caller='test_build_gaussian_smoother')
  # ----------------------------------------------------------------------------
  def test_apply_correlated_perturbations(self):
    """
    Are correlated perturbations the right magnitude, reproducible, and
    spatially coherent?
    """
    timer_start = perf_counter()

    # synthetic unstructured grid with lat/lon coordinates
    lon2d,lat2d = np.meshgrid( np.linspace(0,360,72,endpoint=False)
                             , np.linspace(-80,80,36) )
    lat = lat2d.ravel(); lon = lon2d.ravel()
    ncol = lat.size
    rng = np.random.default_rng(0)
    base = rng.standard_normal(ncol)*5. + 280.
    def make_ds():
      return xr.Dataset({'T':xr.DataArray(base.copy(),dims=['ncol'])
                        ,'lat':xr.DataArray(lat,dims=['ncol'])
                        ,'lon':xr.DataArray(lon,dims=['ncol'])
                        }, coords={'ncol':np.arange(ncol)} )

    corr_length_km = 2000.
    ds1 = hsa.apply_correlated_perturbations( make_ds(), var_list=['T']
                                            , corr_length_km=corr_length_km, seed=123 )
    pert1 = ds1['T'].values - base

    # magnitude should be ~1% of the field std-dev, matching apply_random_perturbations
    self.assertTrue( np.abs( pert1.std()/base.std() - 0.01 ) < 2e-3 )

    # same seed should give identical results
    ds2 = hsa.apply_correlated_perturbations( make_ds(), var_list=['T']
                                            , corr_length_km=corr_length_km, seed=123 )
    self.assertTrue( np.all( ds2['T'].values == ds1['T'].values ) )

    # perturbations should be spatially coherent: nearby columns more
    # correlated than distant ones
    lat_r=np.deg2rad(lat); lon_r=np.deg2rad(lon)
    xyz=np.column_stack([np.cos(lat_r)*np.cos(lon_r)
                        ,np.cos(lat_r)*np.sin(lon_r),np.sin(lat_r)])
    f=(pert1-pert1.mean())/pert1.std()
    i=rng.integers(0,ncol,200000); j=rng.integers(0,ncol,200000)
    chord=np.linalg.norm(xyz[i]-xyz[j],axis=1)
    gc=6371.*2*np.arcsin(np.clip(chord/2,0,1)); prod=f[i]*f[j]
    corr_near = prod[gc<corr_length_km].mean()
    corr_far  = prod[gc>2*corr_length_km].mean()
    self.assertTrue( corr_near > 0.3 )
    self.assertTrue( corr_near > corr_far )

    print_timer(timer_start,caller='test_apply_correlated_perturbations')
  # ----------------------------------------------------------------------------
  def test_create_perturbed_file(self):
    """
    Does create_perturbed_file preserve non-perturbed data, perturb the
    requested variables, and reproduce results for a given seed?
    """
    timer_start = perf_counter()

    lon2d,lat2d = np.meshgrid( np.linspace(0,360,48,endpoint=False)
                             , np.linspace(-80,80,24) )
    lat = lat2d.ravel(); lon = lon2d.ravel()
    ncol = lat.size
    rng = np.random.default_rng(0)
    ds = xr.Dataset({'T':xr.DataArray(rng.standard_normal(ncol)*5+280,dims=['ncol'])
                    ,'Q':xr.DataArray(rng.standard_normal(ncol),dims=['ncol'])
                    ,'lat':xr.DataArray(lat,dims=['ncol'])
                    ,'lon':xr.DataArray(lon,dims=['ncol'])
                    }, coords={'ncol':np.arange(ncol)} )

    with tempfile.TemporaryDirectory() as tmp:
      in_file  = f'{tmp}/in.nc'
      out_file = f'{tmp}/out.nc'
      rep_file = f'{tmp}/rep.nc'
      ds.to_netcdf(in_file)

      hsa.create_perturbed_file( in_file, out_file, var_list=['T']
                               , corr_length_km=2000., seed=7 )
      orig = xr.open_dataset(in_file)
      pert = xr.open_dataset(out_file)

      # non-perturbed variable must be unchanged
      self.assertTrue( np.array_equal( orig['Q'].values, pert['Q'].values ) )
      # perturbed variable must actually change, ~1% of the field std-dev
      d = pert['T'].values - orig['T'].values
      self.assertTrue( d.std()>0. )
      self.assertTrue( np.abs( d.std()/orig['T'].std().values - 0.01 ) < 3e-3 )

      # same seed reproduces the perturbed file exactly
      hsa.create_perturbed_file( in_file, rep_file, var_list=['T']
                               , corr_length_km=2000., seed=7 )
      rep = xr.open_dataset(rep_file)
      self.assertTrue( np.array_equal( rep['T'].values, pert['T'].values ) )

      # refuse to overwrite an existing file without clobber
      with self.assertRaises(OSError):
        hsa.create_perturbed_file( in_file, out_file, var_list=['T'], seed=1 )

      orig.close(); pert.close(); rep.close()

    print_timer(timer_start,caller='test_create_perturbed_file')
  # ----------------------------------------------------------------------------
  # def test_dry_mass_fixer(self):
  #   """ """
  #   hsa.dry_mass_fixer( ncol, plev, hyai, hybi, wgt, qv, mass_ref, ps_in, ps_out )
  #   self.assertTrue( np.all( ps_in[0] == ps_out[0] ) )
#===============================================================================
if __name__ == '__main__':
    unittest.main()
