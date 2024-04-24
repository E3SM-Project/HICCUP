import os, copy, glob, xarray as xr, numpy as np
plot_backend = 'ngl' # mpl / ngl
if plot_backend == 'mpl': import matplotlib.pyplot as plt
if plot_backend == 'ngl': import ngl
#-------------------------------------------------------------------------------
class tclr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#-------------------------------------------------------------------------------
def print_stat(x,name='(no name)',fmt='f',stat='naxh',indent=''):
  """ Print min, avg, max, and std deviation of input """
  if fmt=='f' : fmt = '%.4f'
  if fmt=='e' : fmt = '%e'
  name_len = 12
  msg = f'{indent}{name:{name_len}}'
  for c in list(stat):
    if c=='h' : msg += '   shp: '+str(x.shape)
    if c=='a' : msg += '   avg: '+fmt%x.mean()
    if c=='n' : msg += '   min: '+fmt%x.min()
    if c=='x' : msg += '   max: '+fmt%x.max()
    if c=='s' : msg += '   std: '+fmt%x.std()
  print(msg); return msg
#-------------------------------------------------------------------------------
def trim_png(fig_file,verbose=True):
   """ crop white space from png file """
   fig_file = fig_file+".png"
   fig_file = fig_file.replace(os.getenv('HOME')+'/Research/E3SM/','')
   if os.path.isfile(fig_file) :
      cmd = "convert -trim +repage "+fig_file+"   "+fig_file
      os.system(cmd)
      if verbose: print("\n"+fig_file+"\n")
   else:
      raise FileNotFoundError(f'\ntrim_png(): {fig_file} does not exist?!\n')
#-------------------------------------------------------------------------------
case,case_name,case_root,case_sub,clr,dsh = [],[],[],[],[],[]
def add_case(case_in,name=None,root=None,sub='run',d=0,c='black',init_date=None):
  # tmp_name = case_in if n is None else n
  # tmp_name = f'E3SM {init_date}'
  case.append(case_in); case_name.append(name)
  case_root.append(root); case_sub.append(sub);
  dsh.append(d) ; clr.append(c)
#-------------------------------------------------------------------------------
sim_var_list,obs_var_list,obs_lev_list,var_str_list = [],[],[],[]
def add_var(sim_var,obs_var,var_str,obs_lev=None):
  sim_var_list.append(sim_var); obs_var_list.append(obs_var);
  obs_lev_list.append(obs_lev); var_str_list.append(var_str)
#-------------------------------------------------------------------------------

fig_file,fig_type = 'fx_skill','png'

init_date = '2020-01-01'
htype = 'h1'

# path to observation data
obs_root = '/global/cfs/projectdirs/m3312/whannah/HICCUP/E3SM_tutorial'
obs_path = f'{obs_root}/ERA5_validation.*.2020-01-*.remap_ne30pg2.nc'

# specify temporal length and frequency of data
spd = 8; time1,time2 = 0,spd*5 # daily files with 3-hourly data

# specify temporal length of reference "climate" (acc only)
ref_time1,ref_time2 = time1,time2 

# list of metrics to calculate
metric_list = ['acc','rmse','mean']

#-------------------------------------------------------------------------------
# Build list of hindcast cases

tmp_root = '/global/homes/w/whannah/E3SM/scratch_pm-cpu'
add_case('E3SM.2024-E3SM-tutorial-hindcast-11.2020-01-01',name='E3SM 11',root=tmp_root,c='blue')
add_case('E3SM.2024-E3SM-tutorial-hindcast-12.2020-01-01',name='E3SM 12',root=tmp_root,c='blue')
add_case('E3SM.2024-E3SM-tutorial-hindcast-13.2020-01-01',name='E3SM 13',root=tmp_root,c='blue')
add_case('E3SM.2024-E3SM-tutorial-hindcast-14.2020-01-01',name='E3SM 14',root=tmp_root,c='blue')
add_case('E3SM.2024-E3SM-tutorial-hindcast-15.2020-01-01',name='E3SM 15',root=tmp_root,c='blue')

#-------------------------------------------------------------------------------
# build list of variables

add_var(sim_var='Z500',obs_var='z',obs_lev=500,var_str='Z500')
add_var(sim_var='T850',obs_var='t',obs_lev=850,var_str='T850')
add_var(sim_var='Q850',obs_var='q',obs_lev=850,var_str='Q850')
add_var(sim_var='U850',obs_var='u',obs_lev=850,var_str='U850')
add_var(sim_var='U200',obs_var='u',obs_lev=200,var_str='U200')
add_var(sim_var='PS',  obs_var='sp',           var_str='PS')

#-------------------------------------------------------------------------------
# specify regional subset via lat/lon bounds

# xlat,xlon,dy,dx =  60,360-45,2,2;

if 'xlat' in locals(): lat1,lat2,lon1,lon2 = xlat-dy/2,xlat+dy/2,xlon-dx/2,xlon+dx/2

if 'lat1' in locals():
  print(f'\n{tclr.RED}  NOTE - regional subset is being applied{tclr.END} - lat: {lat1}:{lat2}  lon: {lon1}:{lon2}')

#-------------------------------------------------------------------------------
num_met = len(metric_list)
num_var = len(sim_var_list)
num_case = len(case)
#---------------------------------------------------------------------------------------------------
# Set up plotting stuff
if plot_backend == 'ngl':
  wks = ngl.open_wks(fig_type,fig_file)
  plot = [None]*num_var*num_met
  res = ngl.Resources()
  res.nglDraw                      = False
  res.nglFrame                     = False
  res.tmXTOn                       = False
  res.tmXBMajorOutwardLengthF      = 0.
  res.tmXBMinorOutwardLengthF      = 0.
  res.tmYLMajorOutwardLengthF      = 0.
  res.tmYLMinorOutwardLengthF      = 0.
  res.tmYLLabelFontHeightF         = 0.015
  res.tmXBLabelFontHeightF         = 0.015
  res.tiXAxisFontHeightF           = 0.015
  res.tiYAxisFontHeightF           = 0.015
  res.tmXBMinorOn                  = False
  res.tmYLMinorOn                  = False
  res.xyLineThicknessF             = 6
  res.tiXAxisString                = '[hours]'
  res.xyLineColors                 = clr
  res.xyDashPatterns               = dsh
#---------------------------------------------------------------------------------------------------
def set_subtitles(wks, plot, left_string='', center_string='', right_string='', font_height=0.01):
  ttres         = ngl.Resources()
  ttres.nglDraw = False
  # Use plot extent to call ngl.text(), otherwise you will see this error:
  # GKS ERROR NUMBER   51 ISSUED FROM SUBROUTINE GSVP  : --RECTANGLE DEFINITION IS INVALID
  strx = ngl.get_float(plot,'trXMinF')
  stry = ngl.get_float(plot,'trYMinF')
  ttres.txFontHeightF = font_height
  # Set annotation resources to describe how close text is to be attached to plot
  amres = ngl.Resources()
  if not hasattr(ttres,'amOrthogonalPosF'):
    amres.amOrthogonalPosF = -0.52   # Top of plot plus a little extra to stay off the border
  else:
    amres.amOrthogonalPosF = ttres.amOrthogonalPosF
  # Add left string
  amres.amJust,amres.amParallelPosF = 'BottomLeft', -0.5   # Left-justified
  tx_id_l   = ngl.text(wks, plot, left_string, strx, stry, ttres)
  anno_id_l = ngl.add_annotation(plot, tx_id_l, amres)
  # Add center string
  amres.amJust,amres.amParallelPosF = 'BottomCenter', 0.0   # Centered
  tx_id_c   = ngl.text(wks, plot, center_string, strx, stry, ttres)
  anno_id_c = ngl.add_annotation(plot, tx_id_c, amres)
  # Add right string
  amres.amJust,amres.amParallelPosF = 'BottomRight', 0.5   # Right-justified
  tx_id_r   = ngl.text(wks, plot, right_string, strx, stry, ttres)
  anno_id_r = ngl.add_annotation(plot, tx_id_r, amres)
  return
#---------------------------------------------------------------------------------------------------
for v,mvar in enumerate(sim_var_list):
  print(f'\n  var: {tclr.MAGENTA}{mvar}{tclr.END}')
  if 'obs_mean' in locals(): del obs_mean
  #---------------------------------------------------------------------------
  # Load Obs data
  print(f'    case: {tclr.CYAN}ERA5{tclr.END}')
  ds_obs = xr.open_mfdataset(obs_path)#.rename({'latitude':'lat','longitude':'lon'})
  data_an = ds_obs[obs_var_list[v]].isel(time=slice(time1,time2))
  # obs_time = data_an.time # this needs to match length of model data
  # area = ds_obs['area'].isel(time=0)
  if 'level' in data_an.coords: data_an = data_an.sel({'level':obs_lev_list[v]})
  #-----------------------------------------------------------------------------
  # unit conversions
  if obs_var_list[v]=='z': data_an = data_an/9.81
  if obs_var_list[v]=='q': data_an = data_an*1e3
  #-----------------------------------------------------------------------------
  # Assign coordinates
  data_an = data_an.assign_coords({'lat':ds_obs['lat'],'lon':ds_obs['lon']})
  #-----------------------------------------------------------------------------
  # define reference state
  data_rf = data_an.isel(time=slice(time1,ref_time2))
  #-----------------------------------------------------------------------------
  # regional subset
  mask = xr.DataArray( np.ones(len(ds_obs['lat']),dtype=bool), coords=ds_obs['lat'].coords )
  if 'lat1' in locals(): mask = mask & (ds_obs['lat']>=lat1) & (ds_obs['lat']<=lat2)
  if 'lon1' in locals(): mask = mask & (ds_obs['lon']>=lon1) & (ds_obs['lon']<=lon2)
  data_an  = data_an.where(mask,drop=True)
  data_rf = data_rf.where(mask,drop=True)
  #-----------------------------------------------------------------------------
  # load to avoid dealing with dask arrays
  data_an.load()
  #-----------------------------------------------------------------------------
  fx_acc_list = []
  fx_rmse_list = []
  fx_mean_list = []
  time_list = []
  for c in range(num_case):
    print(f'    case: {tclr.CYAN}{case[c]}{tclr.END}')
    #---------------------------------------------------------------------------
    # read the simulation data
    ds = xr.open_mfdataset( f'{case_root[c]}/{case[c]}/{case_sub[c]}/{case[c]}.eam.h1.*.nc' ) 
    data_fc = ds[mvar].isel(time=slice(time1,time2))
    #---------------------------------------------------------------------------
    # unit conversions
    if mvar=='TS': data_fc = data_fc + 273.15
    if mvar[0]=='Q': data_fc = data_fc*1e3
    #---------------------------------------------------------------------------
    # regional subset
    data_fc  = data_fc.where(mask,drop=True)
    #---------------------------------------------------------------------------
    # # resample daily
    # data_fc = data_fc.resample(time='3H').mean(dim='time')
    # data_fc = data_fc.isel(time=slice(time1,time2))
    #---------------------------------------------------------------------------
    # deal with potential time mismatch
    if len(data_an.time)>len(data_fc.time): data_an = data_an.isel(time=slice(0,len(data_fc.time)))
    if len(data_fc.time)>len(data_an.time): data_fc = data_fc.isel(time=slice(0,len(data_an.time)))
    data_fc['time'] = data_an['time']
    #---------------------------------------------------------------------------
    # load to avoid dealing with dask arrays
    data_fc.load()
    #---------------------------------------------------------------------------
    # define time coordinate used for plotting
    days_from_zero = ( data_fc['time'] - data_fc['time'][0] ).astype('float')
    days_from_zero = days_from_zero / 3600e9 # convert from nanoseconds to hours
    time_list.append( days_from_zero )
    #---------------------------------------------------------------------------
    for m in range(num_met):
      #-------------------------------------------------------------------------
      if metric_list[m]=='acc':
        if len(data_fc.ncol.values)==1:
          raise ValueError('ACC forecast metric is problematic for a single spatial column')
        # calculate anomaly Correlation Coefficient (ACC) using anomalies from time mean ref state
        data_fc_anomaly = data_fc - data_rf.mean(dim=['time']).values
        data_an_anomaly = data_an - data_rf.mean(dim=['time']).values
        acc_numerator = ( data_fc_anomaly * data_an_anomaly ).sum(dim='ncol').values
        denom_fc = np.sqrt( (data_fc_anomaly**2).sum(dim='ncol').values )
        denom_an = np.sqrt( (data_an_anomaly**2).sum(dim='ncol').values )
        acc = acc_numerator / ( denom_fc * denom_an )
        fx_acc_list.append( acc )
      #-------------------------------------------------------------------------
      # calculate root-mean-square-error from reanalysis
      if metric_list[m]=='rmse':
        rmse = np.sqrt( np.mean( np.square( data_fc.values - data_an.values ), axis=1 ) )
        fx_rmse_list.append( rmse )
      #-------------------------------------------------------------------------
      # calculate spatial mean
      if metric_list[m]=='mean':
        mean = np.mean( data_fc.values, axis=1 )
        fx_mean_list.append( mean )
        if 'obs_mean' not in locals(): obs_mean = np.mean( data_an.values, axis=1 )
  #-----------------------------------------------------------------------------
  # print stats of all skill metrics for this variable
  for m in range(num_met):
    print()
    if metric_list[m]=='mean':  print_stat(obs_mean,       name=f'Mean - {"ERA5":10}',      stat='naxh',indent=' '*4)
    for c in range(num_case):
      if metric_list[m]=='acc': print_stat(fx_acc_list[c], name=f'ACC  - {case_name[c]:10}',stat='naxh',indent=' '*4)
      if metric_list[m]=='rmse':print_stat(fx_rmse_list[c],name=f'RMSE - {case_name[c]:10}',stat='naxh',indent=' '*4)
      if metric_list[m]=='mean':print_stat(fx_mean_list[c],name=f'Mean - {case_name[c]:10}',stat='naxh',indent=' '*4)
  #-----------------------------------------------------------------------------
  for m in range(num_met):
    ip = v*num_met+m
    #---------------------------------------------------------------------------
    if metric_list[m]=='acc':
      tres = copy.deepcopy(res)
      tres.trYMaxF = 1
      tres.trYMinF = np.min(fx_acc_list) - np.std(fx_acc_list)
      tres.tiYAxisString = f'ACC'
      plot[ip] = ngl.xy(wks,np.stack(time_list),np.stack(fx_acc_list),tres)
      set_subtitles(wks, plot[ip], mvar, '', 'ACC', font_height=0.01)
    #---------------------------------------------------------------------------
    if metric_list[m]=='rmse':
      tres = copy.deepcopy(res)
      tres.tiYAxisString = f'RMSE'
      plot[ip] = ngl.xy(wks,np.stack(time_list),np.stack(fx_rmse_list),tres) 
      set_subtitles(wks, plot[ip], mvar, '', 'RMSE', font_height=0.01)
    #---------------------------------------------------------------------------
    if metric_list[m]=='mean':
      tres = copy.deepcopy(res)
      tres.tiYAxisString = f'{mvar}'
      plot[ip] = ngl.xy(wks,np.stack(time_list),np.stack(fx_mean_list),tres)
      tres.xyLineColor = 'black'
      tres.xyMonoLineColor = True
      ngl.overlay(plot[ip], ngl.xy(wks,time_list[0].values,obs_mean,tres))
      set_subtitles(wks, plot[ip], mvar, '', 'Mean', font_height=0.01)
  #-----------------------------------------------------------------------------
  # # Add legend
  # lgres = ngl.Resources()
  # lgres.vpWidthF, lgres.vpHeightF  = 0.05, 0.1
  # lgres.lgLabelFontHeightF = 0.008
  # lgres.lgLineThicknessF   = res.xyLineThicknessF
  # lgres.lgMonoLineColor    = False
  # lgres.lgLineColors       = clr
  # lgres.lgDashIndexes      = dsh
  # lgres.lgLabelJust    = 'CenterLeft'
  # pid = ngl.legend_ndc(wks, len(case_name), case_name, 0.55, 0.9, lgres)
#---------------------------------------------------------------------------------------------------
# Finalize plot
if plot_backend == 'ngl': 
  pres = ngl.Resources()
  pres.nglPanelYWhiteSpacePercent = 5
  pres.nglPanelXWhiteSpacePercent = 5
  layout = [num_var,num_met]
  ngl.panel(wks,plot,layout,pres)
#---------------------------------------------------------------------------------------------------
trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
