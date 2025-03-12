#!/usr/bin/env python
#===================================================================================================
# Unit testing for state_adjustment module
#===================================================================================================
import unittest
import numpy as np
import xarray as xr
from hiccup import hiccup_state_adjustment as hsa

class standard_atmosphere:
  def __init__(self,altitude):
    """ 
    Return a 2-uplet (pressure, temperature) depending on provided altitude.
    Units are SI (m, PA, Kelvin)
    """
    self.altitude = altitude
    self.pressure = np.empty(altitude.shape)
    self.temperature = np.empty(altitude.shape)
    for k in range(len(altitude)):
      if altitude[k]<=11000: 
        # troposphere
        self.pressure[k] = 101325 * (1 - 2.25569E-5 * self.altitude[k])**5.25616
        self.temperature[k] = 288.14 - 0.00649 * self.altitude[k]        
      elif altitude[k]<=20000:  
        self.pressure[k] = 0.223356 * 101325 * np.exp(-0.000157688 * (self.altitude[k] - 11000))  # stratosphere
        self.temperature[k] = 216.66
      else:
        raise ValueError('altitude out of range [0-20000m]')

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
                         }
                         ,coords={'ncol':np.arange(ncol)
                                 ,'lev':np.arange(plev)
                                 ,'ilev':np.arange(ilev)}
                        ).expand_dims(time=1,axis=0).copy(deep=True)
    ds_data = ds_data.chunk(2)
    ds_topo = xr.Dataset({'PHIS':(['ncol'],phis_new)})

    hsa.adjust_surface_pressure( ds_data, ds_topo, debug=False )

    # Get the adjusted surface pressure for value checking
    ps_new = ds_data['PS'].isel(time=0).values

    # for k in range(plev+1): print(f'  {k}  {pressure_int[0,k]:8.2f}  {phis_int[k]:8.2f} ')
    # for i in range(ncol): print(f'phis_old: {phis_old[i]:04.0f}  phis_new: {phis_new[i]:04.0f}  ps_old: {ps_old[i]:6.2f}  ps_new: {ps_new[i]:6.2f}')
    self.assertTrue( ps_new[0] >ps_old[0] )
    self.assertTrue( ps_new[1] <ps_old[1] )
    self.assertTrue( ps_new[2]==ps_old[2] )
  # ----------------------------------------------------------------------------
  def test_adjust_surface_temperature(self):
    """ 
    Does surface temperature interpolation give the right value? 
    """
    phis      = 0.5e3
    phis_old  = np.array([phis]*3)
    phis_new  = np.array([ phis-0.5e3, phis+0.5e3, phis ])
    isa       = standard_atmosphere( phis_old )
    ts_old    = isa.temperature

    ds_data = xr.Dataset({'PHIS':(['ncol'],phis_old)
                         ,'TS'  :(['ncol'],ts_old)})
    ds_topo = xr.Dataset({'PHIS':(['ncol'],phis_new)})

    hsa.adjust_surface_temperature( ds_data, ds_topo )

    ts_new = ds_data['TS'].values

    # for i in range(ncol): print(f'phis_old: {phis_old[i]:04.0f}  phis_new: {phis_new[i]:04.0f}  ts_old: {ts_old[i]:6.2f}  ts_new: {ts_new[i]:6.2f}')
    self.assertTrue( ts_new[0] <ts_old[0] )
    self.assertTrue( ts_new[1] >ts_old[1] )
    self.assertTrue( ts_new[2]==ts_old[2] )
  # ----------------------------------------------------------------------------
  def test_remove_supersaturation(self):
    """ 
    Do supersaturated values get limited correctly? 
    """
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

    hsa.remove_supersaturation( ds )

    rh_out = ds['Q'].values/qv_sat
    expected_answer = np.array([1.0, 1.0, 0.9])
    self.assertTrue( np.all( np.abs(rh_out-expected_answer)<1e-10 ) )
  # ----------------------------------------------------------------------------
  def test_adjust_cld_wtr(self):
    """
    Does cloud water and ice get limited properly?
    """
    vals = [0,-1,1]
    expected_answer = [0,0,1]
    ncol = np.arange(len(vals))
    ds = xr.Dataset({'CLDLIQ':xr.DataArray(vals,dims=['ncol'])
                    ,'CLDICE':xr.DataArray(vals,dims=['ncol'])
                    }, coords={'ncol':ncol} )

    hsa.adjust_cld_wtr( ds )

    self.assertTrue( np.all( np.abs(ds['CLDLIQ'].values-expected_answer)<1e-10 ) )
    self.assertTrue( np.all( np.abs(ds['CLDICE'].values-expected_answer)<1e-10 ) )
  # ----------------------------------------------------------------------------
  def test_adjust_cloud_fraction(self):
    """
    Does cloud fraction get limited properly?
    """
    input_vals = [0,-1,1,2]
    expected_answer = [0,0,1,1]
    ncol = np.arange(len(input_vals))
    ds = xr.Dataset({'CLOUD_FRAC':xr.DataArray(input_vals,dims=['ncol'])
                    }, coords={'ncol':ncol} )

    hsa.adjust_cloud_fraction( ds, frac_var_name='CLOUD_FRAC' )
    
    self.assertTrue( np.all( np.abs(ds['CLOUD_FRAC'].values-expected_answer)<1e-10 ) )
  # ----------------------------------------------------------------------------
  # def test_dry_mass_fixer(self):
  #   """ """
  #   hsa.dry_mass_fixer( ncol, plev, hyai, hybi, wgt, qv, mass_ref, ps_in, ps_out )
  #   self.assertTrue( np.all( ps_in[0] == ps_out[0] ) )
#===============================================================================
if __name__ == '__main__':
    unittest.main()
