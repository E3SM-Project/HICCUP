#!/usr/bin/env python
#===================================================================================================
# Unit testing for state_adjustment module
#===================================================================================================
import unittest
import numpy as np
import state_adjustment

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
  """Tests for state_adjustment.py """
  
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
                                              phis_old, ps_old, ts_old,     \
                                              phis_new, ps_new, ts_new )

    # for k in range(plev+1): print(f'  {k}  {pressure_int[0,k]:8.2f}  {phis_int[k]:8.2f} ')
    # for i in range(ncol): print(f'phis_old: {phis_old[i]:04.0f}  phis_new: {phis_new[i]:04.0f}  ps_old: {ps_old[i]:6.2f}  ps_new: {ps_new[i]:6.2f}')
    self.assertTrue( ps_new[0] >ps_old[0] )
    self.assertTrue( ps_new[1] <ps_old[1] )
    self.assertTrue( ps_new[2]==ps_old[2] )

#===============================================================================
if __name__ == '__main__':
    unittest.main()
