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
# Simple routine for chcecking variable values - useful for debugging
#-------------------------------------------------------------------------------
def print_stat(x,name='(no name)',unit='',fmt='f',stat='naxh',indent='  '):
  """ 
  Simple routine for printing various statistics or properites of a variable.
  The characters of the "stat" string variable are used to specify the order 
  and type of quantities to calculate
    n   minimum value
    a   average across all dimensions
    x   maximum value
    s   standard deviation
    h   shape
  """
  if fmt=='f' : fmt = '%f'
  if fmt=='e' : fmt = '%e'
  if unit!='' : unit = '['+str(unit)+']'
  print(indent+name+" "+unit)
  for c in list(stat):
      if c=='n' : print(indent+'min: '+fmt%x.min() )
      if c=='a' : print(indent+'avg: '+fmt%x.mean())
      if c=='x' : print(indent+'max: '+fmt%x.max() )
      if c=='s' : print(indent+'std: '+fmt%x.std() )
      if c=='h' : print(indent+'shp: '+str(x.shape) )
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
def adjust_surface_pressure( ds_data, ds_topo, pressure_var_name='plev',
                             lev_coord_name='lev', debug=False ):
  """ 
  Adjust the surface pressure based on surace height difference 
  and assumed standard atmosphere lapse rate. Input datasets must
  include the following variables:
    time                  time coordinate
    ncol                  column index coordinate
    <lev_coord_name>      vertical level coordinate
    PS                    input old surface pressure      [Pa]
    PHIS                  input surface geopotential      [m]
    T                     temperature on level centers    [K]
    <pressure_var_name>   pressure on level centers       [Pa]
  the target surface geopotential (PHIS) must also be included in ds_topo
  """

  # Check for required variables in input datasets
  for var in ['time','ncol',lev_coord_name] :
    if var not in ds_data.dims : raise KeyError(f'{var} is missing from ds_data')
  for var in ['PHIS','PS','T',pressure_var_name] :
    if var not in ds_data.variables : raise KeyError(f'{var} is missing from ds_data')
  if 'PHIS' not in ds_topo.variables : raise KeyError('PHIS is missing from ds_data')

  if debug :
    # Debugging print statements
    print('Before Adjustment:')
    print_stat(ds_data['PS'],name='PS (old)')

  nlev = len(ds_data[lev_coord_name])

  # Loop over all time and columns to make adjustment
  for t in range(len(ds_data['time'])) :
    for i in range(len(ds_data['ncol'])) :

      del_phis = ds_data['PHIS'].isel(ncol=i) - ds_topo['PHIS'].isel(ncol=i)

      # Only update PHIS if phis difference is not negligible
      if(np.abs(del_phis) > phis_threshold) :

        # grab pressure and temperature of current column to make code more readable
        pressure = ds_data[pressure_var_name]
        if 'time' in pressure.dims : pressure = pressure.isel(time=t)
        if 'ncol' in pressure.dims : pressure = pressure.isel(ncol=i)
        temperature = ds_data['T'].isel(time=t,ncol=i)

        # move up from surface to define k level for Tbot and Pbot
        # z is calculated from the hydrostatic equation
        z = 0.
        first = True
        for k in range(nlev-1,0,-1) :
          if first :
            # For first level use surface pressure to get dp
            dp = ds_data['PS'].isel(time=t,ncol=i) - pressure[k]
            first = False
          else:
            dp = pressure[k+1] - pressure[k]
          rho = pressure[k] / (Rdair*temperature[k])
          dz = dp/(rho/gravit)
          z = z + dz
          if ( z > dz_min ) : break

        if (k==nlev) : exit(f'ERROR: could not find model level {dz_min} m above the surface ')

        # Define Tbot & Pbot
        tbot  = temperature[k]
        pbot  = pressure[k]

        alpha = lapse*Rdair/gravit                                              # pg 8 eq 6

        # provisional extrapolated surface temperature
        Tstar = tbot + alpha*tbot*( ds_data['PS'][:,i]/pbot - 1.)               # pg 8 eq 5
        T0    = Tstar + lapse*ds_data['PHIS'][:,i]/gravit                       # pg 9 eq 13

        # adjustments for very high (T_ref1) or low (T_ref2) temperatures 
        if (Tstar <= T_ref1) and (T0 > T_ref1) :
          # inhibit low pressure under elevated hot terrain
          alpha = Rdair/ds_topo['PHIS'][i]*( T_ref1 - Tstar )                   # pg 9 eq 14.1
        elif (Tstar > T_ref1) and (T0 > T_ref1) :
          # inhibit low pressure under elevated hot terrain
          Tstar = ( T_ref1 + Tstar )*0.5                                        # pg 9 eq 14.2
        if (Tstar < T_ref2) :
          # inhibit unduly high pressure below elevated cold terrain
          Tstar = ( T_ref2 + Tstar )*0.5                                        # pg 9 eq 14.3

        # Calculate new surface pressure
        beta = del_phis/(gravit*Tstar)
        temp = del_phis/(Rdair*Tstar)*(1. - 0.5*alpha*beta + (1./3.)*(alpha*beta)**2. )
        ds_data['PS'][:,i] = ds_data['PS'][:,i].values * np.exp( temp.values )  # pg 9 eq 12


  if debug :
    # Debugging print statements
    print('After Adjustment:')
    print_stat(ds_data['PS'],name='PS (new)')

  return

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
def adjust_surface_temperature( ds_data, ds_topo, debug=False ):
  """ 
  Adjust the surface temperature based on surace height difference 
  and assumed standard atmosphere lapse rate 
    ds        xarray dataset containing surface temperature and 
              surface geopotential on Model-Coordinateel grid 
    ds_topo   xarray dataset containing smoothed model topography 
              (i.e. target topo)
  """
  
  # Check for required variables in input datasets
  if 'TS'   not in ds_data.variables : 
    raise KeyError('sfc temperature (TS) variable is missing from ds_data')
  if 'PHIS' not in ds_data.variables : 
    raise KeyError('sfc geopotential (PHIS) variable is missing from ds_data')
  if 'PHIS' not in ds_topo.variables : 
    raise KeyError('sfc geopotential (PHIS) variable is missing from ds_topo')
  if len(ds_data['ncol'].values) != len(ds_topo['ncol'].values) : 
    raise IndexError('ncol dimensions of input datasets do not match')

  if debug :
    # Debugging print statements
    print('Before Adjustment:')
    print_stat(ds_data['PHIS'],name='PHIS (old)')
    print_stat(ds_topo['PHIS'],name='PHIS (new)')
    print_stat(ds_data['TS'],name='TS (old)')

  ds_data['TS'] = ds_data['TS'] - ( ds_data['PHIS'] - ds_topo['PHIS'] )*lapse/gravit

  if debug :
    # Debugging print statements
    print('After Adjustment:')
    print_stat(ds_data['TS'],name='TS (new)')

  return 

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
  qv_sat = r_sat / ( 1.0 + r_sat )

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