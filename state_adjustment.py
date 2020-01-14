import xarray as xr
import numpy as np
#-------------------------------------------------------------------------------
lapse   = 0.0065        # std. atmosphere lapse rate              ~ -6.5 K/km
gravit  = 9.80616       # acceleration of gravity                 ~ m/s^2
boltz   = 1.38065e-23   # boltzmann's constant                    ~ J/k/molecule
avogad  = 6.02214e26    # avogadro's number                       ~ molecules/kmole
mwdair  = 28.966        # molecular weight of dry air             ~ kg/kmole
rgas    = avogad*boltz  # universal gas constant                  ~ J/k/kmole
rdair   = rgas/mwdair   # gas constant for dry air                ~ J/k/kg

T_ref1    = 290.5       # reference temperature for sfc adjustments
T_ref2    = 255.0       # reference temperature for sfc adjustments
#-------------------------------------------------------------------------------
# Adjust surface pressure and temperature
# Algorithm based on sea-level pressure calculation
# from section 3.1.b of NCAR NT-396 
# "Vertical Interpolation and Truncation of Model-Coordinate Data"
# https://opensky.ucar.edu/islandora/object/technotes%3A168
# similar to components/cam/src/physics/cam/cpslec.F90
#-------------------------------------------------------------------------------
def adjust_surface_pressure( plev, ncol, temperature, pressure_mid, pressure_int,  \
                             phis_old, ps_old, ts_old, phis_new, ps_new, ts_new ):
  """ 
  Adjust the surface pressure based on new surace height assumed lapse rate 
    plev            # levels
    ncol            # columns
    temperature     temperature on level centers    [K]
    pressure_mid    pressure on level centers       [Pa]
    pressure_int    pressure on level interfaces    [Pa]
    phis_old        input old surface height        [m]
    ps_old          input old surface pressure      [Pa]
    ts_old          input old surface temperature   [K]
    phis_new        input new surface height        [m]
    ps_new          output surface pressure         [Pa]
    ts_new          output surface temperature      [K]
  """
  
  phis_threshold = 0.001  # threshold for determining whether to calculate new pressure
  dz_min = 150.           # min distance of first level above sfc for Tbot and Pbot [m]

  for i in range(ncol) :  

    del_phis = phis_old[i] - phis_new[i]

    # interpolate new surface temperature value
    ts_new[i] = ts_old[i] - lapse*(del_phis/gravit)

    # If difference between analysis and model phis is negligible,
    # then set model Ps = analysis
    if(np.abs(del_phis) <= phis_threshold) :

      ps_new[i] = ps_old[i]

    else:

      # move up from surface to define level for Tbot and Pbot
      z = 0.
      for k in range(plev-1,0,-1) :
        hkk    = 0.5*( pressure_int[i,k+1] - pressure_int[i,k] ) / pressure_mid[i,k]
        z_incr = (rdair/gravit)*temperature[i,k]*hkk
        z = z + z_incr
        if ( z > dz_min ) : break
        z = z + z_incr

      # if (k==plev) : exit(f'Error:  could not find model level above {z_min} ')

      # Define Tbot & Pbot
      tbot  = temperature[i,k]
      pbot  = pressure_mid[i,k]

      alpha = lapse*rdair/gravit                           # pg 8 eq 6

      # provisional extrapolated surface temperature
      Tstar = tbot + alpha*tbot*( ps_old[i]/pbot - 1.)      # pg 8 eq 5
      # T0    = Tstar + lapse*phis_new[i]/gravit              # pg 9 eq 13
      T0    = Tstar + lapse*phis_old[i]/gravit              # pg 9 eq 13

      # adjustments for very high (T_ref1) or low (T_ref2) temperatures 
      if (Tstar <= T_ref1) and (T0 > T_ref1) :
        # inhibit low pressure under elevated hot terrain
        alpha = rdair/phis_new[i]*( T_ref1 - Tstar )        # pg 9 eq 14.1
      elif (Tstar > T_ref1) and (T0 > T_ref1) :
        # inhibit low pressure under elevated hot terrain
        Tstar = ( T_ref1 + Tstar )*0.5                      # pg 9 eq 14.2
      if (Tstar < T_ref2) :
        # inhibit unduly high pressure below elevated cold terrain
        Tstar = ( T_ref2 + Tstar )*0.5                      # pg 9 eq 14.3

      # Calculate new surface pressure
      # beta = phis_new[i]/(rdair*Tstar)
      beta = del_phis/(rdair*Tstar)
      temp = beta*(1. - 0.5*alpha*beta + (1./3.)*(alpha*beta)**2. )
      ps_new[i] = ps_old[i] * np.exp( temp )                # pg 9 eq 12

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------