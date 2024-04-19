import os, ngl, copy, glob, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
#-------------------------------------------------------------------------------
he.default_data_dir=os.getenv('HOME')+'/E3SM/scratch'
he.default_data_sub='run'
#-------------------------------------------------------------------------------
name,case,case_dir,case_sub,clr,dsh,mrk = [],[],[],[],[],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,m=1,c='black'):
  global name,case,case_dir,case_sub
  tmp_name = case_in if n is None else n
  case.append(case_in); name.append(tmp_name)
  case_dir.append(p); case_sub.append(s);
  dsh.append(d) ; clr.append(c) ; mrk.append(m)
#-------------------------------------------------------------------------------
var = []
var_str = []
var_unit = []
file_type_list = []
obs_var_list = []
obs_file_list = []
lev_list = []
def add_var(var_name,obs_var=None,obs_file=None,name='',unit='',lev=None):
   var.append(var_name)
   obs_var_list.append(obs_var)
   obs_file_list.append(obs_file)
   var_str.append(name)
   var_unit.append(unit)
   lev_list.append(lev)
#-------------------------------------------------------------------------------

init_date = '2005-06-01'

add_case('E3SM.2022-HICCUP-SST-TEST-00.F2010.ne30pg2_oECv3',n='CTL',c='red')
add_case('E3SM.2022-HICCUP-SST-TEST-01.F2010.ne30pg2_oECv3',n='EXP',c='blue')


# grid = '180x360'
grid = '90x180'  # use this for ne30pg2!

obs_path = os.getenv('HOME')+f'/HICCUP/data_scratch/ERA5_validation.*.{init_date}.remap_{grid}.nc'

# list of variables to plot
var,lev = [],[]

htype = 'h1'
var.append(('TS','TS')); lev.append(-999)

add_var('LW_flux_up@tom',obs_var='LW_flux_up_at_model_top',obs_file='CERES.LW_flux_up_at_model_top.AVERAGE.ne30pg2.20200126.nc',name='TOA LW up',unit='')


# var.append(('Z500','z')); lev.append(500)

# single time index to load (no averaging) - time max used for reference "climate" calculation
spd = 8
time1,time2 = 0,spd*60
# time2_ref = spd*1  # use first day as reference state (acc only)
time2_ref = time2  # use entire record as reference climate (acc only)

create_legend = False

verbose = True

# output figure type and name


method = 'acc'  # acc / rmse

fig_type,fig_file = os.getenv('HOME')+f'/Research/E3SM/figs_hindcast/hindcast.fx_skill.v1.{method}','png'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
class tcolor: ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'

num_var = len(var)
num_case = len(case)

if 'name' not in vars(): name = case

if 'clr' not in vars(): clr = ['black']*num_case
if 'dsh' not in vars(): dsh = np.zeros(num_case)

# create the plot workstation
wks = ngl.open_wks(fig_type,fig_file)
plot = []
# plot = [None]*(num_var*num_case)

# set oup the plot resources
res = hs.res_xy()
res.tiXAxisString = '[hours]'
res.xyLineThicknessF = 6
# res.xyDashPatterns = [0,1,2,3,4,5,6]
res.xyLineColors = clr
res.xyDashPatterns = dsh

# Open obs dataset
# ds_obs = xr.open_mfdataset(obs_path,combine='by_coords')
# ds_obs = ds_obs.rename({'latitude':'lat','longitude':'lon'})
# ds_obs = ds_obs.isel(time=slice(time1,time2))

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# for v,(mvar,ovar) in enumerate(var.items()):
for v,(mvar,ovar) in enumerate(var):

  print('\n  var: '+clr.MAGENTA+mvar+clr.END+f'   lev={lev[v]}')

  obs_file = os.getenv('HOME')+f'/HICCUP/data_scratch/ERA5_validation.{ovar.upper()}.{init_date}.remap_{grid}.nc'
  print(f'  obs_file: {obs_file}')
  fx_metric_list = []
  time_list = []
  for c in range(num_case):
    #---------------------------------------------------------------------------
    # Load Obs data
    ds_obs = xr.open_dataset(obs_file).rename({'latitude':'lat','longitude':'lon'})
    # ds_obs = ds_obs.isel(time=slice(time1,time2))
    obs_time = ds_obs['time'].isel(time=slice(time1,time2)) # this needs to match length of model data
    data_an = ds_obs[ovar_tmp].isel(time=slice(time1,time2))

    # Select the level - be careful as some fields do not have level coord
    if 'level' in data_an.coords: data_an = data_an.sel({'level':lev[v]})

    if ovar=='z': data_an = data_an/hc.g
    if ovar=='q': data_an = data_an*1e3

    # Assign coordinates
    data_an = data_an.assign_coords({'lat':ds_obs['lat'],'lon':ds_obs['lon']})
    if 'lat1' in locals(): data_an = data_an.sel(lat=slice(lat1,lat2))
    if 'lon1' in locals(): data_an = data_an.sel(lon=slice(lon1,lon2))

    # define reference state
    data_ref = data_an.isel(time=slice(time1,time2_ref))

    #---------------------------------------------------------------------------
    # read the data
    #---------------------------------------------------------------------------
    data_fc = case_obj.load_data(mvar, htype=htype,lev=lev[v],use_remap=True,remap_str=f'remap_{grid}')

    if mvar=='TS': data_fc = data_fc + 273.15

    if 'lat1' in locals(): data_fc = data_fc.sel(lat=slice(lat1,lat2))
    if 'lon1' in locals(): data_fc = data_fc.sel(lon=slice(lon1,lon2))

    if 'plev' in data_fc.coords: data_fc = data_fc.sel({'plev':lev[v]*1e2})

    # print(); print(data_fc); exit()

    data_fc = data_fc.resample(time='3H').mean(dim='time')
    data_fc = data_fc.isel(time=slice(time1,time2))

    # method for dealing with time mismatch due to short runs
    if len(data_an['time'])>len(data_fc['time']):
      data_an = data_an.isel(time=slice(0,len(data_fc['time'])))

    # load to avoid dealing with dask arrays
    data_fc.load()
    data_an.load()

    #---------------------------------------------------------------------------
    # Calculate forecast metric and add to list
    #---------------------------------------------------------------------------
    # Make time coordinate consistent
    data_fc['time'] = data_an['time'] # obs_time

    # calculate nomaly Correlation Coefficient (ACC) using anomalies from time mean ref state
    if method=='acc':
      data_fc = data_fc - data_ref.mean(dim=['time']).values
      data_an = data_an - data_ref.mean(dim=['time']).values
      acc_numerator = ( data_fc * data_an ).sum(dim=['lat','lon']).values
      denom_fc = np.sqrt( (data_fc**2).sum(dim=['lat','lon']).values )
      denom_an = np.sqrt( (data_an**2).sum(dim=['lat','lon']).values )
      acc = acc_numerator / ( denom_fc * denom_an )

      fx_metric_list.append( acc )

    # calculate root-mean-square-error from reanalysis
    if method=='rmse':
      rmse = np.sqrt( np.square( data_fc - data_an ).mean(dim=['lat','lon']) )
      fx_metric_list.append( rmse )

    if verbose: hc.print_stat(fx_metric_list[-1],stat='nxh',compact=True,indent='    ',name=f'{method}')

    # define time coordinate used for plotting
    days_from_zero = ( data_fc['time'] - data_fc['time'][0] ).astype('float')
    days_from_zero = days_from_zero / 3600e9 # convert from nanoseconds to hours
    time_list.append( days_from_zero )

  #-----------------------------------------------------------------------------
  # Create plot
  #-----------------------------------------------------------------------------
  if method=='acc':
    res.trYMaxF = 1
    res.trYMinF = np.min(fx_metric_list) - np.std(fx_metric_list)
    res.tiYAxisString = f'Skill'
  if method=='rsme':
    res.tiYAxisString = f'RMSE'

  # ip = v*num_case+c
  ip = c*num_var+v

  plot.append( ngl.xy(wks,np.stack(time_list),np.stack(fx_metric_list),res) )

  cstr = mvar
  if 'lat1' in locals() and 'lat2' in locals():
    lat1_str = f'{lat1}N' if lat1>=0 else f'{(lat1*-1)}S'
    lat2_str = f'{lat2}N' if lat2>=0 else f'{(lat2*-1)}S'
    cstr += f' ({lat1_str}:{lat2_str})'
  hs.set_subtitles(wks, plot[len(plot)-1], '', cstr, '', font_height=0.015)

  # xyres = hs.res_xy()
  # plot.append( ngl.xy(wks,time.values,max_pt.values,xyres) )
  #-----------------------------------------------------------------------------
  # Ad legend
  #-----------------------------------------------------------------------------
  if create_legend:
    lgres = ngl.Resources()
    # lgres.vpWidthF, lgres.vpHeightF  = 0.1, 0.12
    # lgres.lgLabelFontHeightF = 0.015
    # lgres.lgLineThicknessF   = 16
    lgres.vpWidthF, lgres.vpHeightF  = 0.05, 0.1
    lgres.lgLabelFontHeightF = 0.008
    # lgres.lgLineThicknessF   = 4
    lgres.lgLineThicknessF   = res.xyLineThicknessF
    lgres.lgMonoLineColor    = False
    # lgres.lgMonoDashIndex    = True
    lgres.lgLineColors       = clr
    lgres.lgDashIndexes      = dsh
    lgres.lgLabelJust    = 'CenterLeft'

    # Set legend position
    lx,ly = 0.2, 0.4
    if len(var)==2: lx,ly = 0.2,0.9  # 1x2
    # if len(var)==2: lx,ly = 0.55,0.9  # 2x1
    if len(var)==4: lx,ly = 0.3, 0.9  # 2x2
    if len(var)==6: lx,ly = 0.2,0.8  # 2x3
    if len(var)==8: lx,ly = 0.14,0.64  # 2x4
    # if len(var)==6: lx,ly = 0.16,0.8  # 3x2

    pid = ngl.legend_ndc(wks, len(name), name, lx, ly, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
pres = ngl.Resources()
pres.nglPanelYWhiteSpacePercent = 5
pres.nglPanelXWhiteSpacePercent = 5
# layout = [len(plot),1]
if len(plot)<=3:
  layout = [1,len(plot)]
else:
  layout = [2,np.ceil(len(plot)/2)]
# layout = [np.ceil(len(plot)/2),2]
# layout = [num_case,len(plot)/num_case]
ngl.panel(wks,plot[0:len(plot)],layout,pres)

hc.trim_png(fig_file)