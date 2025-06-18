import xarray as xr, numpy as np, os, datetime, copy, numba
#-------------------------------------------------------------------------------
from hiccup.hiccup_constants import std_lapse, dry_lapse
from hiccup.hiccup_constants import gravit
from hiccup.hiccup_constants import boltz
from hiccup.hiccup_constants import avogad
from hiccup.hiccup_constants import MW_dryair
from hiccup.hiccup_constants import MW_ozone
from hiccup.hiccup_constants import MW_vapor
from hiccup.hiccup_constants import Rgas
from hiccup.hiccup_constants import Rdair
from hiccup.hiccup_constants import Rvapor
from hiccup.hiccup_constants import P0
from hiccup.hiccup_utilities import print_stat
from hiccup.hiccup_utilities import chk_finite

T_ref1    = 290.5       # reference temperature for sfc adjustments
T_ref2    = 255.0       # reference temperature for sfc adjustments

phis_threshold = 1e-3   # threshold for determining if 2 phis values are different
z_min = 150.            # min distance [m] from sfc to minimize effects radiation

lapse_rate = std_lapse
# lapse_rate = dry_lapse

verbose_default = False # local verbosity default

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
                             lev_coord_name='lev', debug=False,
                             verbose=None, verbose_indent='' ):
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
  if verbose is None : verbose = verbose_default
  if verbose: print(f'\n{verbose_indent}Adjusting surface pressure...')
  if debug: print(f'{verbose_indent}adjust_surface_pressure: DEBUG MODE ENABLED')

  # # define minimum threshold to use when dividing by topo height
  # topo_min_value = 10.

  # Make sure to use PHIS_d if file contains both
  if 'PHIS_d' in ds_topo.variables : 
    if 'ncol' in ds_topo.variables: ds_topo = ds_topo.drop(['ncol'])
    if 'PHIS' in ds_topo.data_vars: ds_topo = ds_topo.drop(['PHIS'])
    ds_topo = ds_topo.rename({'PHIS_d':'PHIS','ncol_d':'ncol'})

  rename_ncol = False
  if 'ncol_d' in ds_data.dims: ds_data = ds_data.rename({'ncol_d':'ncol'}) ; rename_ncol = True

  # Check for required variables in input datasets
  for var in ['time','ncol',lev_coord_name] :
    if var not in ds_data.dims : raise KeyError(f'{var} is missing from ds_data')
  for var in ['PS','T',pressure_var_name] :
    if var not in ds_data.variables : raise KeyError(f'{var} is missing from ds_data')
  if 'PHIS' not in ds_data.variables :
    if 'PHIS_d' in ds_data.variables :
      ds_data = ds_data.rename({'PHIS_d':'PHIS'})
    else:
      raise KeyError(f'PHIS is missing from ds_data')

  # Check for required variables in target topography
  if 'PHIS' not in ds_topo.variables : raise KeyError(f'PHIS is missing from ds_data')

  # Check to make sure that [pressure] levels are ordered top to bottom
  if ds_data[lev_coord_name][0] > ds_data[lev_coord_name][-1]:
    raise ValueError(f'The level coordinate ({lev_coord_name}) must be ordered top/low to bottom/high')

  if debug:
    # Debugging print statements
    print(f'{verbose_indent}Before Adjustment:')
    print_stat(ds_data['PS'],name='PS (old)')

  nlev = len(ds_data[lev_coord_name])

  #-----------------------------------------------------------------------------
  # NOTE - Below is the original method that utilizes the minimum altitude (z_min)
  # to use for interpolation, but this is overly costly for high-res grids.
  # This method should be put into a seperate function available as an option.

  # # Make 3D pressure variable with surface pressure field added at the bottom
  # ps_lev_coord = ds_data[pressure_var_name][lev_coord_name]
  # ps_lev_coord = ps_lev_coord.max().values + ps_lev_coord.diff(lev_coord_name).max().values
  # ps_tmp = ds_data['PS'].expand_dims({lev_coord_name:[ps_lev_coord]},axis=-1)
  # pressure = ds_data[pressure_var_name]
  # if 'time' not in pressure.dims : pressure = pressure.expand_dims(time=len(ps_tmp['time']),axis=0)
  # if 'ncol' not in pressure.dims : pressure = pressure.expand_dims(ncol=len(ps_tmp['ncol']),axis=2)
  # # If ps_tmp has extra lat/lon coords they will cause an error, so just drop them
  # if 'lat' in  ps_tmp.coords : ps_tmp = ps_tmp.drop('lat')
  # if 'lon' in  ps_tmp.coords : ps_tmp = ps_tmp.drop('lon')
  # pressure_with_ps = xr.concat( [ pressure, ps_tmp ], dim=lev_coord_name )

  # # calculate pressure thickness
  # dp = pressure_with_ps.isel({lev_coord_name:slice(None,None,-1)}).diff(dim=lev_coord_name)
  # dp = dp.isel({lev_coord_name:slice(None,None,-1)})
  # dp = dp*-1

  # # calculate dz from hydrostatic formula
  # dz = dp / ( gravit * pressure / (Rdair * ds_data['T']) )

  # # integrate dz to get z
  # z = dz.cumsum(dim=lev_coord_name)

  # # Find lowest height exceeding minimum threshold
  # k_coord = xr.DataArray(np.arange(nlev),coords={lev_coord_name:z[lev_coord_name]})
  # kbot_ind = xr.where( z>=z_min, k_coord, -1).max(dim=lev_coord_name)
  # kbot_ind.load() # dask array can't be used in isel() below - so we need to load here
  # if (kbot_ind == -1).any():
  #   raise ValueError(f'ERROR: Could not find model level {z_min} m above surface')

  # # Check that there weren't problems finding the bottom level
  # if np.any(kbot_ind.values==-1) : 
  #   exit(f'ERROR: could not find model level {z_min} m above the surface ')

  # # Define temperature and pressure for "bottom" level
  # tbot = ds_data['T'].isel({lev_coord_name:kbot_ind})
  # pbot = pressure.isel({lev_coord_name:kbot_ind})

  #-----------------------------------------------------------------------------
  # A more performance friendly alternative is to ignore concerns about being
  # too close to the surface and just use the layer closest to the surface

  tbot = ds_data['T'].isel({lev_coord_name:nlev-1})
  pbot = ds_data[pressure_var_name].isel({lev_coord_name:nlev-1})

  #-----------------------------------------------------------------------------

  alpha = lapse_rate*Rdair/gravit                                               # pg 8 eq 6
  
  # provisional extrapolated surface temperature
  Tstar = tbot + alpha*tbot*( ds_data['PS']/pbot - 1.)                          # pg 8 eq 5

  #-----------------------------------------------------------------------------
  # NOTE - The adjustments below originally intended for interpolating data to
  # the mean sea level pressure, and have often been used for initial condition
  # generation without incident. However, tropical cyclone simulations with
  # SCREAM in 2024 revealed that these adjustments can lead to rare edge cases
  # that produce unreasonable values near topography. Disabling the calculations
  # altogether seemed to fix the issue, but they remain here to revisit late.

  # T0 = Tstar + lapse_rate*ds_data['PHIS']/gravit                              # pg 9 eq 13
  
  # # calculate alternate surface geopotential to avoid errors when dividing
  # topo_phis_temp = ds_topo['PHIS']
  # topo_phis_temp = topo_phis_temp.where( topo_phis_temp>topo_min_value, topo_min_value )

  # # The next few lines provide parameter adjustments to deal with  
  # # very high (T_ref1) or low (T_ref2) temperatures 

  # # inhibit low pressure under elevated hot terrain                             pg 9 eq 14.1
  # condition = np.logical_and( Tstar <= T_ref1, T0 > T_ref1 )
  # condition = np.logical_and( condition, ds_topo['PHIS']>topo_min_value )
  # alpha = xr.where(condition, Rdair/topo_phis_temp*(T_ref1-Tstar) , alpha)

  # # inhibit low pressure under elevated hot terrain                             pg 9 eq 14.2
  # condition = np.logical_and( Tstar > T_ref1,  T0 > T_ref1 )
  # condition = np.logical_and( condition, ds_topo['PHIS']>topo_min_value )
  # alpha.values = xr.where(condition, 0, alpha)
  # Tstar.values = xr.where(condition, (T_ref1+Tstar)*0.5 ,Tstar)

  # # inhibit unduly high pressure below elevated cold terrain                    pg 9 eq 14.3
  # condition = ( Tstar < T_ref2 )
  # condition = np.logical_and( condition, ds_topo['PHIS']>topo_min_value )
  # Tstar.values = xr.where(condition, (T_ref2+Tstar)*0.5 ,Tstar)

  # # Calculate new surface pressure                                              pg 9 eq 12
  # del_phis = ds_data['PHIS'] - ds_topo['PHIS']
  # *__, del_phis = xr.broadcast(ds_data['PS'], del_phis)
  # beta = del_phis/(Rdair*Tstar)
  # temp = beta*(1. - 0.5*alpha*beta + (1./3.)*(alpha*beta)**2. )
  # ps_new = ds_data['PS'] * np.exp( temp )
  #-----------------------------------------------------------------------------

  # Calculate new surface pressure                                              pg 9 eq 11
  del_phis = ds_data['PHIS'] - ds_topo['PHIS']
  *__, del_phis = xr.broadcast(ds_data['PS'], del_phis)
  ps_new = ds_data['PS'] * np.power( (1.+alpha*del_phis/(Rdair*Tstar)) , 1./alpha )

  # save attributes to restore later
  ps_attrs = ds_data['PS'].attrs

  if debug: ps_old = ds_data['PS'].copy(deep=True)

  # Only update PHIS if phis difference is not negligible
  ds_data['PS'] = xr.where( np.abs(del_phis)>phis_threshold, ps_new, ds_data['PS'])

  # restore attributes
  ds_data['PS'].attrs = ps_attrs

  # change the dimension name back if it was changed above
  if rename_ncol: ds_data = ds_data.rename({'ncol':'ncol_d'})

  if debug:
    chk_finite(ds_data['PS'],name='ps_new')
    print(f'{verbose_indent}After Adjustment:')
    print_stat(ds_data['PS'],name='PS (new)')
    print_stat(ds_data['PS']-ps_old, name='PS diff')

  return ds_data

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
def adjust_surface_temperature( ds_data, ds_topo, debug=False,
                                verbose=None, verbose_indent='' ):
  """ 
  Adjust the surface temperature based on surace height difference 
  and assumed standard atmosphere lapse rate 
    ds        xarray dataset containing surface temperature and 
              surface geopotential on Model-Coordinateel grid 
    ds_topo   xarray dataset containing smoothed model topography 
              (i.e. target topo)
  """
  if verbose is None : verbose = verbose_default
  if verbose: print(f'\n{verbose_indent}Adjusting surface temperature...')
  if debug: print(f'{verbose_indent}adjust_surface_temperature: DEBUG MODE ENABLED')

  # Make sure to use PHIS_d if file contains both
  if 'PHIS_d' in ds_topo.variables : 
    if 'ncol' in ds_topo.variables: ds_topo = ds_topo.drop(['ncol'])
    if 'PHIS' in ds_topo.data_vars: ds_topo = ds_topo.drop(['PHIS'])
    ds_topo = ds_topo.rename({'PHIS_d':'PHIS','ncol_d':'ncol'})

  # Check for required variables in input datasets
  if 'TS'   not in ds_data.variables : 
    raise KeyError('sfc temperature (TS) variable is missing from ds_data')
  if 'PHIS' not in ds_data.variables : 
    raise KeyError(f'sfc geopotential (PHIS) variable is missing from ds_data')
  if 'PHIS' not in ds_topo.variables : 
    raise KeyError(f'sfc geopotential (PHIS) variable is missing from ds_topo')
  if ds_data.sizes['ncol'] != ds_topo.sizes['ncol'] : 
    topo_ncol = ds_topo.sizes['ncol']
    data_ncol = ds_data.sizes['ncol']
    raise IndexError(f'dimensions of input datasets do not match: data_ncol={data_ncol} / topo_ncol={topo_ncol} ')

  if debug :
    # Debugging print statements
    print(f'{verbose_indent}Before Adjustment:')
    print_stat(ds_data['PHIS'],name='PHIS (old)')
    print_stat(ds_topo['PHIS'],name='PHIS (new)')
    print_stat(ds_data['TS'],name='TS (old)')

  # save attributes to restore later
  ts_attrs = ds_data['TS'].attrs

  ds_data['TS'].values = ds_data['TS'] - ( ds_data['PHIS'] - ds_topo['PHIS'] )*lapse_rate/gravit

  # restore attributes
  ds_data['TS'].attrs = ts_attrs

  if debug :
    # Debugging print statements
    print(f'{verbose_indent}After Adjustment:')
    print_stat(ds_data['TS'],name='TS (new)')

  return ds_data

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# @numba.njit()
# def pressure_interpolation_numba( ntime, ncol, p_mid_new, p_mid_old, T_old, T_new ):
#   # simple interpolation - does not extrapolate
#   for t in range(ntime):
#     for i in range(ncol):
#       T_new[t,:,i] = np.interp( p_mid_new[t,:,i], p_mid_old[t,:,i], T_old[t,:,i] )
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# @numba.njit()
# def adjust_temperature_numba( ntime, ncol, p_mid_new, p_mid_old, T_old, T_new ):
#   """
#   adjust temperature for change in surface pressure using the std atmosphere lapse rate
#   """
#   for t in range(ntime):
#     for i in range(ncol):
#       # pressure thickness
#       dp = p_mid_new[t,:,i] - p_mid_old[t,:,i]
#       # use ideal gass law to get density
#       rho = p_mid_old[t,:,i] / ( Rdair * T_old[t,:,i] )
#       # use hydrostatic euqation to convert dp to dz
#       dz = -1 * dp / ( rho * gravit )
#       # calculate new temperature value
#       T_new[t,:,i] = T_old[t,:,i] + std_lapse*dz

#   return T_new
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def adjust_temperature_eam( ds_data, ps_old, debug=False,
                            verbose=None, verbose_indent='' ):
  """ 
  Adjust the temperature profile based on surace pressure difference 
  created by adjust_surface_pressure(). Input datasets must
  include the following variables:
    time                  time coordinate
    ncol                  column index coordinate
    hyam                  hybrid vertical coordinate A coefficient
    hybm                  hybrid vertical coordinate B coefficient
    PS                    input old surface pressure      [Pa]
    T                     temperature on level centers    [K]
  additionally, the previous sfc pressure "ps_old" must be provided
  """
  if verbose is None : verbose = verbose_default
  if verbose: print(f'\n{verbose_indent}Adjusting temperature profile (EAM)...')
  if debug: print(f'{verbose_indent}adjust_temperature_eam: DEBUG MODE ENABLED')

  rename_ncol = False
  if 'ncol_d' in ds_data.dims: ds_data = ds_data.rename({'ncol_d':'ncol'}) ; rename_ncol=True
  if 'ncol_d' in ps_old.dims : ps_old = ps_old.rename({'ncol_d':'ncol'})

  # Check for required variables in input datasets
  for var in ['time','ncol'] :
    if var not in ds_data.dims : raise KeyError(f'{var} is missing from ds_data')
  for var in ['PS','T','hyam','hybm','P0'] :
    if var not in ds_data.variables : raise KeyError(f'{var} is missing from ds_data')

  # calculate pressure profile for each column
  p_mid_old = ds_data['hyam']*ds_data['P0'] + ds_data['hybm']*ps_old
  p_mid_new = ds_data['hyam']*ds_data['P0'] + ds_data['hybm']*ds_data['PS']

  p_mid_old = p_mid_old.transpose('time','lev','ncol')
  p_mid_new = p_mid_new.transpose('time','lev','ncol')

  T_old = ds_data['T'].copy(deep=True)

  # T_adj = adjust_temperature_numba( len(ds_data.time), len(ds_data.ncol),
  #                                   p_mid_new.values, p_mid_old.values,
  #                                   T_old.values, ds_data['T'].values)
  # ds_data['T'] = ( ds_data['T'].dims, T_adj )

  # pressure thickness
  dp = p_mid_new - p_mid_old
  
  # use ideal gass law to get density
  rho = p_mid_old / ( Rdair * T_old )
  
  # use hydrostatic euqation to convert dp to dz
  dz = -1 * dp / ( rho * gravit )

  # calculate new temperature value
  ds_data['T'] = T_old + lapse_rate*dz

  # change the dimension name back if it was changed above
  if rename_ncol: ds_data = ds_data.rename({'ncol':'ncol_d'})

  return ds_data

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def remove_supersaturation( ds, hybrid_lev=False, pressure_var_name='plev',
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
  qv_sat.values = xr.where(qv_sat.values>=0.0,qv_sat,1.0)

  # Calculate relative humidity for limiter
  rh = ds['Q'] / qv_sat

  if debug:
    print(); print_stat(rh,name='rh in remove_supersaturation')

  # save attributes to restore later
  tmp_attrs = ds['Q'].attrs

  # Apply limiter conditions
  ds['Q'] = xr.where(rh.values>1.,qv_sat,ds['Q'])
  ds['Q'] = xr.where(rh.values<0.,qv_min,ds['Q'])
  
  # restore attributes
  ds['Q'].attrs = tmp_attrs

  if debug:
    print(); print_stat(ds['Q'],name='qv in remove_supersaturation after adjustment')

  return ds

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def adjust_cld_wtr( ds, verbose=None, verbose_indent='' ):
  """
  Adjust cloud water to remove negative values
  """
  if verbose is None : verbose = verbose_default
  if verbose: print(f'\n{verbose_indent}Adjusting cloud water...')

  for var in ['CLDLIQ','CLDICE']:
    if var in ds.data_vars: ds[var].values = xr.where( ds[var].values>=0, ds[var], 0. )

  return ds

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def adjust_cloud_fraction( ds, frac_var_name='FRAC', verbose=None, verbose_indent=''):
  """
  Adjust cloud fraction to remove values outside of [0,1]
  """
  if verbose is None : verbose = verbose_default
  if verbose: print(f'\n{verbose_indent}Adjusting cloud fraction...')

  ds[frac_var_name].values = xr.where(ds[frac_var_name]>=0, ds[frac_var_name], 0. )
  ds[frac_var_name].values = xr.where(ds[frac_var_name]<=1, ds[frac_var_name], 1. )

  return ds

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def apply_random_perturbations( ds, var_list=None, seed=None,
                                verbose=None, verbose_indent='' ):
  """
  Apply random perturbations to the final remapped state variables
  """
  if verbose is None : verbose = verbose_default
  if verbose: print(f'\n{verbose_indent}Applying random perturbation...')

  if var_list is None:
    raise ValueError(f'var_list cannot be None')

  for var in var_list:
    if var not in ds.variables:
      raise KeyError(f'{var} is missing from data')

  if seed is None:
    seed = int(datetime.datetime.utcnow().strftime('%s'))
    seed = seed*hash(os.getenv('USER'))
    seed = seed*hash(' '.join(os.listdir()))
    seed = np.abs(seed)

  # initialize RNG
  rng = np.random.default_rng(seed)

  # apply perturbations
  for var in var_list:
    # use "small" perturbations => 1% of std-dev
    ds[var] = ds[var] + rng.standard_normal( ds[var].shape ) * ds[var].std().values * 0.01

  return ds

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
def get_pressure_from_hybrid( ds, a_coeff_name='hyam', b_coeff_name='hybm' ):
  """
  Calculate 3D pressure field from hybrid vertical coordinates
  following the formulation for CESM/E3SM
  """
  pressure = ds[a_coeff_name] * ds['P0'] + ds[b_coeff_name] * ds['PS']

  # Make sure dimensions are in correct order for mid-point levels
  if a_coeff_name=='hyam' and all(d in ds.dims for d in ['time','lev','ncol']):
    pressure = pressure.transpose('time','lev','ncol')
  
  # Make sure dimensions are in correct order for interface levels
  if a_coeff_name=='hyai' and all(d in ds.dims for d in ['time','ilev','ncol']):
    pressure = pressure.transpose('time','ilev','ncol')
  
  return pressure

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
