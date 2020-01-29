import xarray as xr
import numpy as np
#-------------------------------------------------------------------------------
lapse   = 0.0065        # std. atmosphere lapse rate              ~ -6.5 K/km
gravit  = 9.80616       # acceleration of gravity                 ~ m/s^2
boltz   = 1.38065e-23   # boltzmann's constant                    ~ J/k/molecule
avogad  = 6.02214e26    # avogadro's number                       ~ molecules/kmole
mwdair  = 28.966        # molecular weight of dry air             ~ kg/kmole
mwvapor = 18.0          # molecular weight of dry air             ~ kg/kmole
Rgas    = avogad*boltz  # universal gas constant                  ~ J/k/kmole
Rdair   = Rgas/mwdair   # gas constant for dry air                ~ J/k/kg
Rvapor  = Rgas/mwvapor
P0      = 1e5           # reference pressure

T_ref1    = 290.5       # reference temperature for sfc adjustments
T_ref2    = 255.0       # reference temperature for sfc adjustments

phis_threshold = 1e-3   # threshold for determining if 2 phis values are different
dz_min = 150.           # min distance [m] from sfc to minimize effects radiation
#-------------------------------------------------------------------------------
# Adjust surface pressure
# Algorithm based on sea-level pressure calculation
# from section 3.1.b of NCAR NT-396 
# "Vertical Interpolation and Truncation of Model-Coordinate Data"
# https://opensky.ucar.edu/islandora/object/technotes%3A168
# similar to components/cam/src/physics/cam/cpslec.F90
# Also see: IFS Documentation Cycle CY23r4, 
# Part VI: Technical and Computational Procedures, 
# Chapter 2 FULL-POS post-processing and interpolation
#-------------------------------------------------------------------------------
def adjust_surface_pressure( plev, ncol, temperature, pressure_mid, pressure_int,  \
                             phis_old, ps_old, phis_new, ps_new ):
  """ 
  Adjust the surface pressure based on new surace height assumed lapse rate 
    plev            # levels
    ncol            # columns
    temperature     temperature on level centers    [K]
    pressure_mid    pressure on level centers       [Pa]
    pressure_int    pressure on level interfaces    [Pa]
    phis_old        input old surface height        [m]
    ps_old          input old surface pressure      [Pa]
    phis_new        input new surface height        [m]
    ps_new          output surface pressure         [Pa]
  """

  for i in range(ncol) :  

    del_phis = phis_old[i] - phis_new[i]

    # If difference between analysis and model phis is negligible,
    # then set model Ps = analysis
    if(np.abs(del_phis) <= phis_threshold) :

      ps_new[i] = ps_old[i]

    else:

      # move up from surface to define k level for Tbot and Pbot
      # zis calculated from the hydrostatic equation
      z = 0.
      for k in range(plev-1,0,-1) :
        hkk    = 0.5*( pressure_int[i,k+1] - pressure_int[i,k] ) / pressure_mid[i,k]
        z_incr = (Rdair/gravit)*temperature[i,k]*hkk
        z = z + z_incr
        if ( z > dz_min ) : break
        z = z + z_incr

      if (k==plev) : exit(f'ERROR: could not find model level {dz_min} m above the surface ')

      # Define Tbot & Pbot
      tbot  = temperature[i,k]
      pbot  = pressure_mid[i,k]

      alpha = lapse*Rdair/gravit                           # pg 8 eq 6

      # provisional extrapolated surface temperature
      Tstar = tbot + alpha*tbot*( ps_old[i]/pbot - 1.)      # pg 8 eq 5
      # T0    = Tstar + lapse*phis_new[i]/gravit              # pg 9 eq 13
      T0    = Tstar + lapse*phis_old[i]/gravit              # pg 9 eq 13

      # adjustments for very high (T_ref1) or low (T_ref2) temperatures 
      if (Tstar <= T_ref1) and (T0 > T_ref1) :
        # inhibit low pressure under elevated hot terrain
        alpha = Rdair/phis_new[i]*( T_ref1 - Tstar )        # pg 9 eq 14.1
      elif (Tstar > T_ref1) and (T0 > T_ref1) :
        # inhibit low pressure under elevated hot terrain
        Tstar = ( T_ref1 + Tstar )*0.5                      # pg 9 eq 14.2
      if (Tstar < T_ref2) :
        # inhibit unduly high pressure below elevated cold terrain
        Tstar = ( T_ref2 + Tstar )*0.5                      # pg 9 eq 14.3

      # Calculate new surface pressure
      # beta = phis_new[i]/(Rdair*Tstar)
      beta = del_phis/(gravit*Tstar)
      temp = del_phis/(Rdair*Tstar)*(1. - 0.5*alpha*beta + (1./3.)*(alpha*beta)**2. )
      ps_new[i] = ps_old[i] * np.exp( temp )                # pg 9 eq 12

#-------------------------------------------------------------------------------
# Adjust surface temperature
# Algorithm based on sea-level pressure calculation
# from section 3.1.b of NCAR NT-396 
# "Vertical Interpolation and Truncation of Model-Coordinate Data"
# https://opensky.ucar.edu/islandora/object/technotes%3A168
# similar to components/cam/src/physics/cam/cpslec.F90
# Also see: IFS Documentation Cycle CY23r4, 
# Part VI: Technical and Computational Procedures, 
# Chapter 2 FULL-POS post-processing and interpolation
#-------------------------------------------------------------------------------
def adjust_surface_temperature( ncol, phis_old, ts_old, phis_new, ts_new ):
  """ 
  Adjust the surface temperature based on new surace height assumed lapse rate 
    ncol            # columns
    ts_old          input old surface temperature   [K]
    ts_new          output surface temperature      [K]
  """
  for i in range(ncol) :  
    del_phis = phis_old[i] - phis_new[i]
    ts_new[i] = ts_old[i] - lapse*(del_phis/gravit)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def remove_supersaturation( qv, temperature, pressure ):
  """
  Adjust the surface temperature based on new surace height assumed lapse rate 
    ncol            # columns
    qv              specific humidity
    temperature     temperature at layer mid-points [k]
    pressure        pressure at layer mid-points    [Pa]
  """
  qv_min = 1.0e-9   # minimum specific humidity value allowed

  # Calculate saturation specific humidity
  qv_sat = calculate_qv_sat_liq(temperature,pressure/1e2)

  # The following check is to avoid the generation of negative values
  # that can occur in the upper stratosphere and mesosphere
  qv_sat = qv_sat.where(qv_sat>=0.0,other=1.0)

  # Calculate relative humidity for limiter
  rh = qv / qv_sat

  # Apply limiter conditions
  qv.values = xr.where(rh>1.,qv_sat,qv)
  qv.values = xr.where(rh<0.,qv_min,qv)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
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
  qv_sat = r_sat / ( 1 + r_sat )

  return qv_sat
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def calculate_qv_sat_ice( temperature, pressure ):
  """ 
  calculate saturation specific humidity [kg/kg]
  from temperature [K] and pressure [hPa]
  """

  # Calculate saturation vapor pressure over ice 
  # Chapter 4 of WMO GUIDE TO METEOROLOGICAL INSTRUMENTS AND METHODS OF OBSERVATION
  # https://www.wmo.int/pages/prog/www/IMOP/CIMO-Guide.html
  ei = 6.112 * np.exp( 22.46*(temperature-273.0)/(temperature-273.0+272.62) )

  # Convert to mixing ratio
  r_sat = (Rdair/Rvapor) * ei / (pressure - ei)

  # Convert mixing ratio to saturation specific humidity
  qv_sat = r_sat / ( 1 + r_sat )

  return qv_sat
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# def dry_mass_fixer( ncol, plev, hyai, hybi, wgt, qv, mass_ref, ps_in, ps_out ):
#   """ 
#   NOT TESTED - THIS APPEARS TO ONLY BE FOR THE SPECTRAL DYCOR (EUL)?
#   Adjust atmospheric mass based upon qv 
#     plev            # levels
#     ncol            # columns
#     hyai            hybrid coefficient for level interfaces
#     hybi            hybrid coefficient for level interfaces
#     wgt             integration weights
#     qv              specific humidity
#     mass_ref        Dry mass of Ref. atmosphere
#     ps_in           input surface pressure
#     ps_out          output adjusted surface pressure
#   """

#   # Compute separate pdel's from "A" and "B" portions of 
#   # hybrid vertical grid for later use in global integrals
#   pdela = np.empty([ncol,plev])
#   pdelb = np.empty([ncol,plev])
#   for i in range(ncol):
#     for k in range(plev):
#       pdela[i,k] = ( hyai[k+1] - hyai[k] )*P0
#       pdelb[i,k] = ( hybi[k+1] - hybi[k] )*ps_in[i]

#   # Compute integrals of mass, moisture, and geopotential height
#   ps_sum  = 0.
#   for i in range(ncol): ps_sum  = ps_sum  + wgt[i]*ps_in[i]
#   mass_init = ps_sum/ncol
#   mass_qv1 = 0.
#   mass_qv2 = 0.

#   # Calculate global integrals needed for water vapor adjustment
#   for k in range(plev):
#     dotproda = 0.
#     dotprodb = 0.
#     for i in range(ncol):
#       dotproda = dotproda + wgt[i]*qv[i,k]*pdela[i,k]
#       dotprodb = dotprodb + wgt[i]*qv[i,k]*pdelb[i,k]
#     mass_qv1 = mass_qv1 + dotproda/ncol
#     mass_qv2 = mass_qv2 + dotprodb/ncol

#   # Normalize average mass, height
#   mass_init = mass_init*0.5/gravit
#   mass_qv1 = mass_qv1*.5/gravit
#   mass_qv2 = mass_qv2*.5/gravit

#   # Compute and apply an initial mass fix factor 
#   # which preserves horizontal gradients of ln(ps)
#   mass_fix = (mass_ref + mass_qv1)/(mass_init - mass_qv2)

#   for i in range(ncol): ps_out[i] = ps_in[i]*mass_fix

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------