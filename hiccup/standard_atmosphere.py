import numpy as np
#===============================================================================
'''
The "U.S. Standard Atmosphere 1976" is a model of how the pressure, temperature,
density, and viscosity of the Earth's atmosphere changes with altitude.
It is defined as having a temperature of 288.15 K at the sea level (1013.25 hPa)

Geometrical altitude is the scale of elevation we would measure with a tape
measure. Geopotential altitude is based on a scale that relates altitude to
gravitational equipotentials, or surfaces of constant gravitational potential
energy per unit mass.

The atmosphere are divided as follows:
  - Troposphere  - 0  to 11 km (36.000 ft) altitude
  - Stratosphere - 11 to 51 km (167.000 ft) altitude
  - Mesosphere   - 51 to 71 km (232.000 ft) altitude
  - Ionosphere   - above 71 km (232.000 ft) altitude
'''
#===============================================================================
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
