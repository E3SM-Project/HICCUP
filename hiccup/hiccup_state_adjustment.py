import xarray as xr, numpy as np, os, datetime
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
z_min = 150.            # min distance [m] from sfc to minimize effects radiation

default_verbose = False

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
# Routines for checking for any number of invalid values
#-------------------------------------------------------------------------------
def chk_finite(x,name=None):
  """
  check the input data for inf values
  input data should be a xarray DataArray 
  """
  inf_cnt = x.where( xr.ufuncs.isinf(x) ).count().values
  nan_cnt = x.where( xr.ufuncs.isnan(x) ).count().values
  if inf_cnt>0 or nan_cnt>0: 
    err_msg = '  '
    if name is not None: err_msg += f'{name}:  '
    err_msg += f'invalid values found! {inf_cnt} infs  /  {nan_cnt} nans  '
    raise ValueError(err_msg)
  return
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
                             lev_coord_name='lev', debug=False, verbose=None ):
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
  if verbose is None : verbose = default_verbose
  if verbose: print('\nAdjusting surface pressure...')
  if debug: print('adjust_surface_pressure: DEBUG MODE ENABLED')

  # define minimum threshold to use when dividing by topo height
  topo_min_value = 10.

  # Make sure to use PHIS_d if file contains both
  if 'PHIS_d' in ds_topo.variables : 
    if 'ncol' in ds_topo.variables: ds_topo = ds_topo.drop(['ncol'])
    if 'PHIS' in ds_topo.data_vars: ds_topo = ds_topo.drop(['PHIS'])
    ds_topo = ds_topo.rename({'PHIS_d':'PHIS','ncol_d':'ncol'})

  if 'ncol_d' in ds_data.dims :
    ds_data = ds_data.rename({'ncol_d':'ncol'})

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

  if debug :
    # Debugging print statements
    print('Before Adjustment:')
    print_stat(ds_data['PS'],name='PS (old)')

  nlev = len(ds_data[lev_coord_name])
        
  # Make 3D pressure variable with surface pressure field added at the bottom
  ps_lev_coord = ds_data[pressure_var_name][lev_coord_name]
  ps_lev_coord = ps_lev_coord.max().values + ps_lev_coord.diff(lev_coord_name).max().values
  ps_tmp = ds_data['PS'].expand_dims({lev_coord_name:[ps_lev_coord]},axis=-1)
  pressure = ds_data[pressure_var_name]
  if 'time' not in pressure.dims : pressure = pressure.expand_dims(time=len(ps_tmp['time']),axis=0)
  if 'ncol' not in pressure.dims : pressure = pressure.expand_dims(ncol=len(ps_tmp['ncol']),axis=2)
  # If ps_tmp has extra lat/lon coords they will cause an error, so just drop them
  if 'lat' in  ps_tmp.coords : ps_tmp = ps_tmp.drop('lat')
  if 'lon' in  ps_tmp.coords : ps_tmp = ps_tmp.drop('lon')
  pressure_with_ps = xr.concat( [ pressure, ps_tmp ], dim=lev_coord_name )

  # calculate pressure thickness
  dp = pressure_with_ps.isel({lev_coord_name:slice(None,None,-1)}).diff(dim=lev_coord_name)
  dp = dp.isel({lev_coord_name:slice(None,None,-1)})
  dp = dp*-1

  # calculate dz from hydrostatic formula
  dz = dp / ( gravit * pressure / (Rdair * ds_data['T']) )

  # calculate z by integrating hydrostatic equation
  # And populate k index array for finding the bottom level
  z = dz.copy(deep=True)
  k_ind = xr.full_like(z,-1,dtype=int)
  # if we don't load these we get a performance bottleneck - not sure why...
  z.load(); k_ind.load()
  for k in range(nlev-1,0-1,-1) :
    k_ind[:,k,:] = np.full_like(k_ind[:,k,:],k)
    if k <= nlev-2 : z[:,k,:] = z[:,k,:] + dz[:,k,:]

  # Find the lowest height that exceeds the minimum
  kbot_ind = xr.where(z>=z_min,k_ind,-1).max(dim=lev_coord_name)

  # Check that there weren't problems finding the bottom level
  if np.any(kbot_ind.values==-1) : 
    exit(f'ERROR: could not find model level {z_min} m above the surface ')

  # Define temperature and pressure for "bottom" level
  tbot = ds_data['T'].isel({lev_coord_name:kbot_ind})
  pbot = pressure.isel({lev_coord_name:kbot_ind})
  
  alpha = lapse*Rdair/gravit                                                    # pg 8 eq 6
  
  # provisional extrapolated surface temperature
  Tstar = tbot + alpha*tbot*( ds_data['PS']/pbot - 1.)                          # pg 8 eq 5
  T0    = Tstar + lapse*ds_data['PHIS']/gravit                                  # pg 9 eq 13
  
  # calculate alternate surface geopotential to avoid errors when dividing
  topo_phis_temp = ds_topo['PHIS']
  topo_phis_temp = topo_phis_temp.where( topo_phis_temp>topo_min_value, topo_min_value )

  # The next few lines provide parameter adjustments to deal with  
  # very high (T_ref1) or low (T_ref2) temperatures 

  # inhibit low pressure under elevated hot terrain                              pg 9 eq 14.1
  condition = np.logical_and( Tstar <= T_ref1, T0 > T_ref1 )
  condition = np.logical_and( condition, ds_topo['PHIS']>topo_min_value )
  alpha = xr.where(condition, Rdair/topo_phis_temp*(T_ref1-Tstar) , alpha)

  # inhibit low pressure under elevated hot terrain                              pg 9 eq 14.2
  condition = np.logical_and( Tstar > T_ref1,  T0 > T_ref1 )
  condition = np.logical_and( condition, ds_topo['PHIS']>topo_min_value )
  alpha.values = xr.where(condition, 0, alpha)
  Tstar.values = xr.where(condition, (T_ref1+Tstar)*0.5 ,Tstar)

  # inhibit unduly high pressure below elevated cold terrain                     pg 9 eq 14.3
  condition = ( Tstar < T_ref2 )
  condition = np.logical_and( condition, ds_topo['PHIS']>topo_min_value )
  Tstar.values = xr.where(condition, (T_ref2+Tstar)*0.5 ,Tstar)
  del_phis = ds_data['PHIS'] - ds_topo['PHIS']

  # Calculate new surface pressure                                               pg 9 eq 12
  *__, del_phis = xr.broadcast(ds_data['PS'], del_phis)
  beta = del_phis/(Rdair*Tstar)
  temp = beta*(1. - 0.5*alpha*beta + (1./3.)*(alpha*beta)**2. )
  ps_new = ds_data['PS'] * np.exp( temp )

  # save attributes to restore later
  ps_attrs = ds_data['PS'].attrs

  # Only update PHIS if phis difference is not negligible
  ds_data['PS'] = xr.where( np.abs(del_phis) > phis_threshold, ps_new, ds_data['PS'])

  # restore attributes
  ds_data['PS'].attrs = ps_attrs

  if debug :
    print()
    print_stat(alpha,name='alpha')
    print()
    chk_finite(xr.DataArray( np.exp( temp.values ) ),name='etmp_da')
    chk_finite(ds_data['PS'],name='ps_new')
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
def adjust_surface_temperature( ds_data, ds_topo, debug=False, verbose=None ):
  """ 
  Adjust the surface temperature based on surace height difference 
  and assumed standard atmosphere lapse rate 
    ds        xarray dataset containing surface temperature and 
              surface geopotential on Model-Coordinateel grid 
    ds_topo   xarray dataset containing smoothed model topography 
              (i.e. target topo)
  """
  if verbose is None : verbose = default_verbose
  if verbose: print('\nAdjusting surface temperature...')
  if debug: print('adjust_surface_temperature: DEBUG MODE ENABLED')

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
    print('Before Adjustment:')
    print_stat(ds_data['PHIS'],name='PHIS (old)')
    print_stat(ds_topo['PHIS'],name='PHIS (new)')
    print_stat(ds_data['TS'],name='TS (old)')

  # save attributes to restore later
  ts_attrs = ds_data['TS'].attrs

  ds_data['TS'].values = ds_data['TS'] - ( ds_data['PHIS'] - ds_topo['PHIS'] )*lapse/gravit

  # restore attributes
  ds_data['TS'].attrs = ts_attrs

  if debug :
    # Debugging print statements
    print('After Adjustment:')
    print_stat(ds_data['TS'],name='TS (new)')

  return 

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def remove_supersaturation( ds, hybrid_lev=False, pressure_var_name='plev', 
                            debug=False, verbose=None ):
  """
  Adjust the surface temperature based on new surace height assumed lapse rate 
    ncol            # columns
    qv              specific humidity
    temperature     temperature at layer mid-points [k]
    pressure        pressure at layer mid-points    (convert to hPa for qv_sat calculation)
  """
  if verbose is None : verbose = default_verbose
  if verbose: print('\nRemoving super saturated data points...')
  if debug: print('remove_supersaturation: DEBUG MODE ENABLED')

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

  return

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def adjust_cld_wtr( ds, verbose=None ):
  """
  Adjust cloud water to remove negative values
  """
  if verbose is None : verbose = default_verbose
  if verbose: print('\nAdjusting cloud water...')

  for var in ['CLDLIQ','CLDICE']:
    if var in ds.data_vars: ds[var].values = xr.where( ds[var].values>=0, ds[var], 0. )

  return

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def adjust_cloud_fraction( ds, frac_var_name='FRAC', verbose=None):
  """
  Adjust cloud fraction to remove values outside of [0,1]
  """
  if verbose is None : verbose = default_verbose
  if verbose: print('\nAdjusting cloud fraction...')

  ds[frac_var_name].values = xr.where(ds[frac_var_name]>=0, ds[frac_var_name], 0. )
  ds[frac_var_name].values = xr.where(ds[frac_var_name]<=1, ds[frac_var_name], 1. )
  return

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def apply_random_perturbations( ds, var_list=None, seed=None, verbose=None):
  """
  Apply random perturbations to the final remapped state variables
  """
  if verbose is None : verbose = default_verbose
  if verbose: print('\nApplying random perturbation...')

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

  return

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
