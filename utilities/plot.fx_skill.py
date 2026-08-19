import os, xarray as xr, numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
def mpl_linestyle(d):
  """ map an NGL dash pattern index to a matplotlib linestyle """
  return {0:'solid', 1:(0,(6,3)), 2:(0,(2,2)), 3:(0,(6,3,2,3)), 4:(0,(1,2))}.get(d,'solid')
#-------------------------------------------------------------------------------
case,case_name,case_root,case_sub,clr,dsh = [],[],[],[],[],[]
def add_case(case_in,name=None,root=None,sub='run',d=0,c='black',init_date=None):
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
fig,axs = plt.subplots(num_var,num_met,figsize=(4*num_met,2.5*num_var),squeeze=False)
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
    ax = axs[v,m]
    #---------------------------------------------------------------------------
    if metric_list[m]=='acc':
      for c in range(num_case):
        ax.plot(time_list[c].values,fx_acc_list[c],color=clr[c],
                linestyle=mpl_linestyle(dsh[c]),linewidth=2)
      ax.set_ylim( np.min(fx_acc_list)-np.std(fx_acc_list), 1 )
      ax.set_ylabel('ACC')
      ax.set_title(mvar,loc='left'); ax.set_title('ACC',loc='right')
    #---------------------------------------------------------------------------
    if metric_list[m]=='rmse':
      for c in range(num_case):
        ax.plot(time_list[c].values,fx_rmse_list[c],color=clr[c],
                linestyle=mpl_linestyle(dsh[c]),linewidth=2)
      ax.set_ylabel('RMSE')
      ax.set_title(mvar,loc='left'); ax.set_title('RMSE',loc='right')
    #---------------------------------------------------------------------------
    if metric_list[m]=='mean':
      for c in range(num_case):
        ax.plot(time_list[c].values,fx_mean_list[c],color=clr[c],
                linestyle=mpl_linestyle(dsh[c]),linewidth=2)
      ax.plot(time_list[0].values,obs_mean,color='black',linewidth=2)
      ax.set_ylabel(mvar)
      ax.set_title(mvar,loc='left'); ax.set_title('Mean',loc='right')
    #---------------------------------------------------------------------------
    ax.set_xlabel('[hours]')
#---------------------------------------------------------------------------------------------------
# Add legend
handles = [ Line2D([0],[0],color=clr[c],linestyle=mpl_linestyle(dsh[c]),lw=2,label=case_name[c])
           for c in range(num_case) ]
handles.append( Line2D([0],[0],color='black',lw=2,label='ERA5') )
fig.legend(handles=handles,loc='upper right',fontsize=8)
#---------------------------------------------------------------------------------------------------
# Finalize plot
fig.tight_layout()
fig.savefig(f'{fig_file}.{fig_type}',dpi=150,bbox_inches='tight')
plt.close(fig)
print(f'\n{fig_file}.{fig_type}\n')
#---------------------------------------------------------------------------------------------------
