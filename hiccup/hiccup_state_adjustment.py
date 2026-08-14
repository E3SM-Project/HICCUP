import xarray as xr, numpy as np, os, datetime, copy
#-------------------------------------------------------------------------------
from hiccup.hiccup_constants import std_lapse
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
from hiccup.hiccup_constants import rearth
from hiccup.hiccup_utilities import print_stat
from hiccup.hiccup_utilities import chk_finite

T_ref1    = 290.5       # reference temperature for sfc adjustments
T_ref2    = 255.0       # reference temperature for sfc adjustments

phis_threshold = 1e-3   # threshold for determining if 2 phis values are different
z_min = 150.            # min distance [m] from sfc to minimize effects radiation

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
                             lev_coord_name='lev', hybrid_lev=False,
                             debug=False, verbose=None, verbose_indent='' ):
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

  if hybrid_lev:
    pbot = ds_data['hyam'].isel({lev_coord_name:nlev-1}) * ds_data['P0'] \
          +ds_data['hybm'].isel({lev_coord_name:nlev-1}) * ds_data['PS']
  else:
    pbot = ds_data[pressure_var_name].isel({lev_coord_name:nlev-1})

  #-----------------------------------------------------------------------------

  alpha = std_lapse*Rdair/gravit                                               # pg 8 eq 6
  
  # provisional extrapolated surface temperature
  Tstar = tbot + alpha*tbot*( ds_data['PS']/pbot - 1.)                          # pg 8 eq 5

  #-----------------------------------------------------------------------------
  # NOTE - The adjustments below originally intended for interpolating data to
  # the mean sea level pressure, and have often been used for initial condition
  # generation without incident. However, tropical cyclone simulations with
  # SCREAM in 2024 revealed that these adjustments can lead to rare edge cases
  # that produce unreasonable values near topography. Disabling the calculations
  # altogether seemed to fix the issue, but they remain here to revisit late.

  # T0 = Tstar + std_lapse*ds_data['PHIS']/gravit                              # pg 9 eq 13
  
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

  ds_data['TS'].values = ds_data['TS'] - ( ds_data['PHIS'] - ds_topo['PHIS'] )*std_lapse/gravit

  # restore attributes
  ds_data['TS'].attrs = ts_attrs

  if debug :
    # Debugging print statements
    print(f'{verbose_indent}After Adjustment:')
    print_stat(ds_data['TS'],name='TS (new)')

  return ds_data

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

  # pressure thickness
  dp = p_mid_new - p_mid_old
  
  # use ideal gass law to get density
  rho = p_mid_old / ( Rdair * T_old )
  
  # use hydrostatic euqation to convert dp to dz
  dz = -1 * dp / ( rho * gravit )

  # calculate new temperature value
  ds_data['T'] = T_old + std_lapse*dz

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
def build_gaussian_smoother( lat, lon, corr_length_km, ncol_name='ncol',
                             n_sigma=3. ):
  """
  Build a sparse, row-normalized Gaussian smoothing matrix keyed to a physical
  correlation length, for low-pass filtering fields on an unstructured (ncol)
  grid. Returns a scipy sparse matrix S of shape (ncol,ncol) such that S @ f
  smooths the field f. Smoothing IID Gaussian noise this way yields a field
  whose ~1/e spatial correlation length is approximately corr_length_km.

  lat,lon        : column center latitude/longitude [degrees]
  corr_length_km : approximate 1/e spatial correlation length [km]
  n_sigma        : neighbor search cutoff in units of the kernel std-dev
  """
  from scipy.spatial import cKDTree
  from scipy.sparse import coo_matrix, diags

  if corr_length_km <= 0:
    raise ValueError(f'corr_length_km must be positive, got {corr_length_km}')
  if n_sigma <= 0:
    raise ValueError(f'n_sigma must be positive, got {n_sigma}')

  # convert the requested correlation length into the Gaussian kernel std-dev.
  # smoothing white noise with a kernel of std sigma yields an autocorrelation
  # that reaches 1/e near 2*sigma, so use sigma = corr_length_km/2 to make the
  # interface argument behave like the resulting correlation length.
  sigma_km = corr_length_km / 2.

  # convert lat/lon to xyz on the unit sphere so the KD-tree handles the poles
  # and the antimeridian seam correctly
  lat_r = np.deg2rad(np.asarray(lat,dtype=np.float64))
  lon_r = np.deg2rad(np.asarray(lon,dtype=np.float64))
  xyz = np.column_stack([ np.cos(lat_r)*np.cos(lon_r),
                          np.cos(lat_r)*np.sin(lon_r),
                          np.sin(lat_r) ])
  tree = cKDTree(xyz)

  # neighbor search cutoff as a chord length on the unit sphere
  rearth_km = rearth/1e3
  theta_cut = min( n_sigma*sigma_km/rearth_km, np.pi )
  chord_cut = 2.*np.sin(theta_cut/2.)

  # sparse chord-distance matrix between all columns within the cutoff, then
  # convert chord distance to great-circle distance [km] for the kernel weight
  dmat = tree.sparse_distance_matrix(tree, max_distance=chord_cut,
                                     output_type='coo_matrix')
  chord = np.clip(dmat.data, 0., 2.)
  gc_km = rearth_km * 2.*np.arcsin(chord/2.)
  weights = np.exp( -0.5*(gc_km/sigma_km)**2 )

  S = coo_matrix((weights,(dmat.row,dmat.col)), shape=dmat.shape).tocsr()

  # row-normalize so the filter preserves the field mean (partition of unity)
  row_sum = np.asarray(S.sum(axis=1)).ravel()
  row_sum[row_sum==0.] = 1.
  S = diags(1./row_sum) @ S

  return S
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def apply_correlated_perturbations( ds, var_list=None, corr_length_km=1000.,
                                    seed=None, smoother=None, ncol_name='ncol',
                                    n_sigma=3., verbose=None, verbose_indent='' ):
  """
  Apply spatially correlated random perturbations to the final remapped state
  variables on an unstructured (ncol) grid. As with apply_random_perturbations
  the magnitude is 1% of each variable's standard deviation, but here the noise
  is low-pass filtered with a Gaussian kernel so that it "looks" spatially
  coherent (synoptic-scale) rather than grid-point static. The smoothing matrix
  is built once from the grid coordinates and reused across all variables and
  vertical levels.

  corr_length_km : approximate 1/e spatial correlation length of the
                   perturbations [km] (default 1000)
  smoother       : optional prebuilt matrix from build_gaussian_smoother() to
                   reuse across calls (e.g. when generating many ensemble
                   members); built from the grid coordinates when None
  """
  if verbose is None : verbose = verbose_default
  if verbose: print(f'\n{verbose_indent}Applying spatially correlated random '
                    f'perturbation (corr_length_km={corr_length_km:g})...')

  if var_list is None:
    raise ValueError(f'var_list cannot be None')

  for var in var_list:
    if var not in ds.variables:
      raise KeyError(f'{var} is missing from data')

  if ncol_name not in ds.dims:
    raise ValueError(f'{ncol_name} dimension not found in data; correlated '
                     f'perturbations require an unstructured grid')

  if seed is None:
    seed = int(datetime.datetime.utcnow().strftime('%s'))
    seed = seed*hash(os.getenv('USER'))
    seed = seed*hash(' '.join(os.listdir()))
    seed = np.abs(seed)

  # initialize RNG
  rng = np.random.default_rng(seed)

  # build the smoothing operator once from the grid coordinates (unless a
  # prebuilt smoother was provided for reuse across ensemble members)
  ncol = ds.sizes[ncol_name]
  if smoother is None:
    for coord in ['lat','lon']:
      if coord not in ds.variables:
        raise KeyError(f'{coord} is required for correlated perturbations but '
                       f'is missing from data')
    smoother = build_gaussian_smoother( ds['lat'].values, ds['lon'].values,
                                        corr_length_km, ncol_name=ncol_name,
                                        n_sigma=n_sigma )
  if smoother.shape[0]!=ncol or smoother.shape[1]!=ncol:
    raise ValueError(f'smoother shape {smoother.shape} does not match '
                     f'{ncol_name} size {ncol}')
  S = smoother

  # apply perturbations
  for var in var_list:
    # generate an independent smoothed noise field for each non-ncol slice
    # (e.g. each vertical level), operating on the ncol axis
    if ncol_name not in ds[var].dims:
      raise ValueError(f'variable {var!r} does not have dimension {ncol_name!r}; '
                       f'all variables in var_list must contain the {ncol_name!r} dimension')
    axis = ds[var].dims.index(ncol_name)
    var_data = np.moveaxis( ds[var].values, axis, -1 )
    lead_shape = var_data.shape[:-1]
    n_lead = int(np.prod(lead_shape)) if lead_shape else 1
    iid = rng.standard_normal((n_lead,ncol))
    noise = (S @ iid.T).T                          # (n_lead,ncol)
    # smoothing shrinks variance, so rescale each slice back to unit std
    noise_std = noise.std(axis=1,keepdims=True)
    noise_std[noise_std==0.] = 1.
    noise = noise/noise_std
    # scale to 1% of the variable std-dev, matching apply_random_perturbations
    noise = noise.reshape((*lead_shape,ncol)) * ds[var].std().values * 0.01
    ds[var].values = ds[var].values + np.moveaxis( noise, -1, axis )

  return ds

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def create_perturbed_file( input_file, output_file, var_list=['T','PS','U','V'],
                           spatially_correlated=True, corr_length_km=1000.,
                           seed=None, smoother=None, ncol_name='ncol', n_sigma=3.,
                           clobber=False, verbose=None, verbose_indent='' ):
  """
  Write a copy of input_file with random perturbations applied to var_list.

  The file is copied so that all non-perturbed variables and metadata are
  preserved byte-for-byte; only the perturbed variables are overwritten in
  place, which keeps each variable's dtype, chunking, compression, and fill
  value intact. Intended to be called repeatedly (e.g. seed=member_index in a
  loop) to build a perturbation ensemble from a single final IC file.

  spatially_correlated : if True, low-pass filter the perturbations so they are
                         spatially coherent (synoptic-scale) rather than
                         grid-point noise
  corr_length_km       : approximate 1/e spatial correlation length [km]
  smoother             : optional prebuilt matrix from build_gaussian_smoother()
                         to avoid rebuilding it on each call; the (possibly
                         newly built) smoother is returned so it can be reused

  Returns the smoothing matrix when spatially_correlated is True (or when a
  smoother was supplied) so that it can be passed back in to avoid rebuilding
  it for each ensemble member; returns None when spatially_correlated is False
  and no smoother was provided.
  """
  import shutil
  import netCDF4
  if verbose is None : verbose = verbose_default
  if verbose: print(f'\n{verbose_indent}Creating perturbed file => {output_file}')

  if os.path.abspath(input_file)==os.path.abspath(output_file):
    raise ValueError('input_file and output_file must be different')
  if not os.path.exists(input_file):
    raise OSError(f'input_file does not exist => {input_file}')
  if os.path.exists(output_file) and not clobber:
    raise OSError(f'output_file already exists (use clobber=True) => {output_file}')

  # copy the file so all non-perturbed variables and metadata are preserved
  # exactly; only the perturbed variables will be overwritten below
  out_dir = os.path.dirname(output_file)
  if out_dir!='' and not os.path.exists(out_dir): os.makedirs(out_dir)
  shutil.copy2(input_file, output_file)

  # load only the variables to be perturbed (plus coordinates) into memory
  with xr.open_dataset(input_file) as ds_in:
    present = [v for v in var_list if v in ds_in.variables]
    missing = [v for v in var_list if v not in ds_in.variables]
    if missing:
      print(f'{verbose_indent}WARNING: skipping variables not found in '
            f'{input_file} => {missing}')
    if not present:
      raise KeyError(f'none of var_list found in {input_file} => {var_list}')
    load_vars = list(present)
    if spatially_correlated:
      for coord in ['lat','lon']:
        if coord not in ds_in.variables:
          raise KeyError(f'{coord} is required for correlated perturbations but '
                         f'is missing from {input_file}')
        load_vars.append(coord)
    ds = ds_in[load_vars].load()

  # apply the perturbations in memory using the tested numeric routines
  if spatially_correlated:
    if smoother is None:
      smoother = build_gaussian_smoother( ds['lat'].values, ds['lon'].values,
                                          corr_length_km, ncol_name=ncol_name,
                                          n_sigma=n_sigma )
    ds = apply_correlated_perturbations( ds, var_list=present,
                                         corr_length_km=corr_length_km, seed=seed,
                                         smoother=smoother, ncol_name=ncol_name,
                                         n_sigma=n_sigma, verbose=False,
                                         verbose_indent=verbose_indent )
  else:
    ds = apply_random_perturbations( ds, var_list=present, seed=seed,
                                     verbose=False, verbose_indent=verbose_indent )

  # write the perturbed variable values back into the copy in place
  with netCDF4.Dataset(output_file, 'a') as nc:
    for var in present:
      nc.variables[var][:] = ds[var].values

  return smoother if spatially_correlated else None
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
