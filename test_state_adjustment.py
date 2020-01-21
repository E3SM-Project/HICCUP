#!/usr/bin/env python
#===================================================================================================
# Unit testing for state_adjustment module
#===================================================================================================
import unittest
import numpy as np
import xarray as xr
import hiccup_state_adjustment

class standard_atmosphere:
  def __init__(self,altitude):
    """ return a 2-uplet (pressure, temperature) depending on provided altitude.
    Units are SI (m, PA, Kelvin)"""
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
  """ Tests for state_adjustment.py """
  #-----------------------------------------------------------------------------
  def test_adjust_surface_pressure(self):
    """ Does surface pressure interpolation give the right value? """
    phis_int        = np.array([ 2.5e3, 1.5e3, 0.5e3 ])   # interface altitudes to define plev
    phis_min        = np.min(phis_int)
    phis_new        = np.array([ phis_min-0.5e3, phis_min+0.5e3, phis_min ])
    plev            = len(phis_int)-1
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
    ps_new        = np.empty(ncol)
    ts_new        = np.empty(ncol)

    state_adjustment.adjust_surface_pressure( plev, ncol, temperature_mid,  \
                                              pressure_mid, pressure_int,   \
                                              phis_old, ps_old, phis_new, ps_new )

    # for k in range(plev+1): print(f'  {k}  {pressure_int[0,k]:8.2f}  {phis_int[k]:8.2f} ')
    # for i in range(ncol): print(f'phis_old: {phis_old[i]:04.0f}  phis_new: {phis_new[i]:04.0f}  ps_old: {ps_old[i]:6.2f}  ps_new: {ps_new[i]:6.2f}')
    self.assertTrue( ps_new[0] >ps_old[0] )
    self.assertTrue( ps_new[1] <ps_old[1] )
    self.assertTrue( ps_new[2]==ps_old[2] )
  #-----------------------------------------------------------------------------
  def test_adjust_surface_temperature(self):
    """ Does surface temperature interpolation give the right value? """
    phis      = 0.5e3
    phis_old  = np.array([phis]*3)
    phis_new  = np.array([ phis-0.5e3, phis+0.5e3, phis ])
    isa       = standard_atmosphere( phis_old )
    ts_old    = isa.temperature
    ncol      = len(phis_old)
    ts_new    = np.empty(ncol)

    state_adjustment.adjust_surface_temperature( ncol, phis_old, ts_old, phis_new, ts_new )

    # for i in range(ncol): print(f'phis_old: {phis_old[i]:04.0f}  phis_new: {phis_new[i]:04.0f}  ts_old: {ts_old[i]:6.2f}  ts_new: {ts_new[i]:6.2f}')
    self.assertTrue( ts_new[0] <ts_old[0] )
    self.assertTrue( ts_new[1] >ts_old[1] )
    self.assertTrue( ts_new[2]==ts_old[2] )
  #-----------------------------------------------------------------------------
  def test_remove_supersaturation(self):
    """ do negative values get limited correctly? """
    temperature_in = 300
    pressure_in    = 1010e2
    qv_sat = state_adjustment.calculate_qv_sat(temperature_in,pressure_in/1e2)
    qv = xr.DataArray(np.array([ 1.1*qv_sat , 1.0*qv_sat , 0.9*qv_sat ]))
    ncol = len(qv.values)
    temperature = xr.DataArray([temperature_in]*ncol)
    pressure    = xr.DataArray([pressure_in]   *ncol)
    state_adjustment.remove_supersaturation( qv, temperature, pressure )
    rh_out = qv/qv_sat
    self.assertTrue( np.all( rh_out.values==np.array([1.0, 1.0, 0.9]) ) )
  #-----------------------------------------------------------------------------
  # def test_dry_mass_fixer(self):
  #   """ """
  #   state_adjustment.dry_mass_fixer( ncol, plev, hyai, hybi, wgt, qv, mass_ref, ps_in, ps_out )
  #   self.assertTrue( np.all( ps_in[0] == ps_out[0] ) )
#===============================================================================
if __name__ == '__main__':
    unittest.main()
