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

init_date = '2020-01-01'

obs_root = '/global/cfs/projectdirs/m3312/whannah/HICCUP/E3SM_tutorial'
obs_path = f'{obs_root}/ERA5_validation.*.2020-01-*.remap_ne30pg2.nc'

#-------------------------------------------------------------------------------
# Build list of hindcast cases to load
# tmp_root = '/global/cfs/projectdirs/m3312/whannah/HICCUP/E3SM_tutorial'
# add_case('E3SM.2024-E3SM-tutorial-hindcast.2020-01-01',root=tmp_root,c='red',init_date='2020-01-01')

tmp_root = '/global/homes/w/whannah/E3SM/scratch_pm-cpu'
add_case('E3SM.2024-E3SM-tutorial-hindcast-01.2020-01-01',name='E3SM 01',root=tmp_root,c='red')
add_case('E3SM.2024-E3SM-tutorial-hindcast-02.2020-01-01',name='E3SM 02',root=tmp_root,c='red')
add_case('E3SM.2024-E3SM-tutorial-hindcast-03.2020-01-01',name='E3SM 03',root=tmp_root,c='red')
add_case('E3SM.2024-E3SM-tutorial-hindcast-04.2020-01-01',name='E3SM 04',root=tmp_root,c='red')
add_case('E3SM.2024-E3SM-tutorial-hindcast-05.2020-01-01',name='E3SM 05',root=tmp_root,c='red')

add_case('E3SM.2024-E3SM-tutorial-hindcast-11.2020-01-01',name='E3SM 11',root=tmp_root,c='blue')
add_case('E3SM.2024-E3SM-tutorial-hindcast-12.2020-01-01',name='E3SM 12',root=tmp_root,c='blue')
add_case('E3SM.2024-E3SM-tutorial-hindcast-13.2020-01-01',name='E3SM 13',root=tmp_root,c='blue')
add_case('E3SM.2024-E3SM-tutorial-hindcast-14.2020-01-01',name='E3SM 14',root=tmp_root,c='blue')
# add_case('E3SM.2024-E3SM-tutorial-hindcast-15.2020-01-01',name='E3SM 15',root=tmp_root,c='blue')

#-------------------------------------------------------------------------------
# build list of variables
add_var(sim_var='Z500',obs_var='z',obs_lev=500,var_str='Z500')
# add_var(sim_var='T850',obs_var='t',obs_lev=850,var_str='T850')
# add_var(sim_var='Q850',obs_var='q',obs_lev=850,var_str='Q850')
# add_var(sim_var='U850',obs_var='u',obs_lev=850,var_str='U850')
# add_var(sim_var='U200',obs_var='u',obs_lev=200,var_str='U200')
add_var(sim_var='PS',obs_var='sp',var_str='PS')

#-------------------------------------------------------------------------------

htype = 'h1'

spd = 8; time1,time2 = 0,spd*20 # daily files with 3-hourly data
time2_ref = time2 # use entire record as reference climate (acc only)

fig_file,fig_type = 'fx_skill','png'

num_var = len(sim_var_list)
num_case = len(case)

#---------------------------------------------------------------------------------------------------
# Set up plotting stuff
if plot_backend == 'ngl':
  wks = ngl.open_wks(fig_type,fig_file)
  plot = [None]*num_var*2
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

   ### Use plot extent to call ngl.text(), otherwise you will see this error:
   ### GKS ERROR NUMBER   51 ISSUED FROM SUBROUTINE GSVP  : --RECTANGLE DEFINITION IS INVALID
   strx = ngl.get_float(plot,'trXMinF')
   stry = ngl.get_float(plot,'trYMinF')
   ttres.txFontHeightF = font_height

   ### Set annotation resources to describe how close text is to be attached to plot
   amres = ngl.Resources()
   if not hasattr(ttres,'amOrthogonalPosF'):
      amres.amOrthogonalPosF = -0.52   # Top of plot plus a little extra to stay off the border
   else:
      amres.amOrthogonalPosF = ttres.amOrthogonalPosF

   ### Add left string
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
  #---------------------------------------------------------------------------
  # Load Obs data
  print(f'    case: {tclr.CYAN}ERA5{tclr.END}')
  ds_obs = xr.open_mfdataset(obs_path)#.rename({'latitude':'lat','longitude':'lon'})
  data_an = ds_obs[obs_var_list[v]].isel(time=slice(time1,time2))
  obs_time = data_an.time # this needs to match length of model data
  # area = ds_obs['area'].isel(time=0)
  if 'level' in data_an.coords: data_an = data_an.sel({'level':obs_lev_list[v]})
  #-----------------------------------------------------------------------------
  # unit conversions
  if obs_var_list[v]=='z': data_an = data_an/9.81
  if obs_var_list[v]=='q': data_an = data_an*1e3
  #-----------------------------------------------------------------------------
  # # Assign coordinates
  # data_an = data_an.assign_coords({'lat':ds_obs['lat'],'lon':ds_obs['lon']})
  #-----------------------------------------------------------------------------
  # define reference state
  data_ref = data_an.isel(time=slice(time1,time2_ref))
  #-----------------------------------------------------------------------------
  # # regional subset
  # if 'lat1' in locals(): data_an = data_an.sel(lat=slice(lat1,lat2))
  # if 'lon1' in locals(): data_an = data_an.sel(lon=slice(lon1,lon2))
  #-----------------------------------------------------------------------------
  # load to avoid dealing with dask arrays
  data_an.load()
  #-----------------------------------------------------------------------------
  fx_acc_list = []
  fx_rmse_list = []
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
    # print_stat(data_an,name=f'data_an - {mvar}',stat='naxh',indent=' '*6)
    # print_stat(data_fc,name=f'data_fc - {mvar}',stat='naxh',indent=' '*6)
    # exit()
    #---------------------------------------------------------------------------
    # # regional subset
    # if 'lat1' in locals(): data_fc = data_fc.sel(lat=slice(lat1,lat2))
    # if 'lon1' in locals(): data_fc = data_fc.sel(lon=slice(lon1,lon2))
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
    # calculate anomaly Correlation Coefficient (ACC) using anomalies from time mean ref state
    data_fc_anomaly = data_fc - data_ref.mean(dim=['time']).values
    data_an_anomaly = data_an - data_ref.mean(dim=['time']).values
    acc_numerator = ( data_fc_anomaly * data_an_anomaly ).sum(dim='ncol').values
    denom_fc = np.sqrt( (data_fc_anomaly**2).sum(dim='ncol').values )
    denom_an = np.sqrt( (data_an_anomaly**2).sum(dim='ncol').values )
    acc = acc_numerator / ( denom_fc * denom_an )
    fx_acc_list.append( acc )
    #---------------------------------------------------------------------------
    # calculate root-mean-square-error from reanalysis
    rmse = np.sqrt( np.mean( np.square( data_fc.values - data_an.values ), axis=1 ) )
    fx_rmse_list.append( rmse )
  #-----------------------------------------------------------------------------
  # # print stats of all skill metrics for this variable
  # for n in range(2):
  #   print()
  #   for c in range(num_case): 
  #     if n==0: print_stat(fx_acc_list[c],name=f'ACC  - {case_name[c]}',stat='nxh',indent=' '*4)
  #     if n==1: print_stat(fx_rmse_list[c],name=f'RMSE - {case_name[c]}',stat='nxh',indent=' '*4)
  #-----------------------------------------------------------------------------
  # plot ACC metric
  tres = copy.deepcopy(res)
  tres.trYMaxF = 1
  tres.trYMinF = np.min(fx_acc_list) - np.std(fx_acc_list)
  tres.tiYAxisString = f'ACC'
  plot[v*2+0] = ngl.xy(wks,np.stack(time_list),np.stack(fx_acc_list),tres)
  set_subtitles(wks, plot[v*2+0], mvar, '', 'ACC', font_height=0.015)
  #-----------------------------------------------------------------------------
  # plot RMSE
  tres = copy.deepcopy(res)
  tres.tiYAxisString = f'RMSE'
  plot[v*2+1] = ngl.xy(wks,np.stack(time_list),np.stack(fx_rmse_list),tres) 
  set_subtitles(wks, plot[v*2+1], mvar, '', 'RMSE', font_height=0.015)
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
  layout = [num_var,2]
  ngl.panel(wks,plot,layout,pres)
#---------------------------------------------------------------------------------------------------
trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
