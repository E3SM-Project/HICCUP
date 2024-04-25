#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
import os, datetime, subprocess as sp, datetime
from shutil import copy2
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd): print('\n'+clr.GREEN+cmd+clr.END) ; os.system(cmd); return
#---------------------------------------------------------------------------------------------------
newcase,config,build,submit = False,False,False,False

acct = 'e3sm'
reservation = 'e3sm_testrun'
src_dir = # PATH TO E3SM SOURCE CODE

newcase      = True
config       = True
build        = True
submit       = True

stop_opt,stop_n,resub,walltime = 'ndays',5,0,'0:30:00'

init_date = datetime.datetime.strptime('2023-09-08 00', '%Y-%m-%d %H')

case = '.'.join(['E3SM','2024-E3SM-tutorial-hindcast',init_date.strftime('%Y-%m-%d')])

# init_scratch  = '/global/cfs/projectdirs/e3sm/www/Tutorials/2024/practicum/day_4/atm_breakout'
init_scratch  = # PATH TO INITIAL CONDITON FILES
init_file_atm = f'{init_scratch}/HICCUP.atm_era5.{init_date.strftime("%Y-%m-%d")}.ne30np4.L80.nc'
init_file_sst = f'{init_scratch}/HICCUP.sst_noaa.{init_date.strftime("%Y-%m-%d")}.nc'
scratch_root = f'{os.getenv("SCRATCH")}/e3sm_scratch/pm-cpu'
case_root = f'{scratch_root}/{case}'
#---------------------------------------------------------------------------------------------------
print(f'\n  case : {case}\n')
#---------------------------------------------------------------------------------------------------
if newcase :
   cmd = f'{src_dir}/cime/scripts/create_newcase'
   cmd += f' --case {case}'
   cmd += f' --output-root {case_root} --script-root {case_root}/case_scripts '
   cmd += f' --handle-preexisting-dirs u '
   cmd += f' --compset F2010 --res ne30pg2_oECv3 '
   cmd += f' --mach pm-cpu --pecount 256x1 '
   cmd += f' --project {acct} '
   cmd += f' --walltime {walltime} '
   run_cmd(cmd)
#---------------------------------------------------------------------------------------------------
os.chdir(f'{case_root}/case_scripts')
if config : run_cmd('./case.setup --reset')
if build : run_cmd('./case.build')
#---------------------------------------------------------------------------------------------------
if submit : 
   #----------------------------------------------------------------------------
   # Namelist options
   nfile = 'user_nl_eam'
   file = open(nfile,'w') 
   file.write(f' ncdata = \'{init_file_atm}\' \n')
   # Specify history output frequency and variables - daily files with 3-hourly data
   file.write(' nhtfrq    = 0,-3 \n')
   file.write(' mfilt     = 1,8 \n')
   file.write(" fincl2 = 'PS','TS','PSL'")                     # sfc temperature and pressure
   file.write(          ",'PRECT','TMQ'")                      # precipitation 
   file.write(          ",'TGCLDLWP','TGCLDIWP'")              # liq & ice water path
   file.write(          ",'LHFLX','SHFLX'")                    # surface fluxes
   file.write(          ",'FSNT','FLNT','FLUT'")               # Net TOM rad fluxes
   file.write(          ",'FLNS','FSNS'")                      # Net sfc rad fluxes
   file.write(          ",'TBOT','QBOT','UBOT','VBOT'")        # lowest model level
   file.write(          ",'T850','Q850','U850','V850','Z850'") # 850mb
   file.write(          ",'T500','Q500','U500','V500','Z500'") # 500mb
   file.write(          ",'T200','Q200','U200','V200','Z200'") # 200mb
   file.close()
   #----------------------------------------------------------------------------
   # Specify start date and SST file for hindcast
   sst_yr = int(init_date.strftime('%Y'))
   os.system(f'./xmlchange --file env_run.xml  RUN_STARTDATE={init_date.strftime("%Y-%m-%d")}')
   os.system(f'./xmlchange --file env_run.xml  SSTICE_DATA_FILENAME={init_file_sst}')
   os.system(f'./xmlchange --file env_run.xml  SSTICE_YEAR_ALIGN={sst_yr}')
   os.system(f'./xmlchange --file env_run.xml  SSTICE_YEAR_START={sst_yr}')
   os.system(f'./xmlchange --file env_run.xml  SSTICE_YEAR_END={sst_yr+1}')
   #----------------------------------------------------------------------------
   # Set some run-time stuff
   if 'stop_opt' in locals(): run_cmd(f'./xmlchange STOP_OPTION={stop_opt}')
   if 'stop_n'   in locals(): run_cmd(f'./xmlchange STOP_N={stop_n}')
   if 'resub'    in locals(): run_cmd(f'./xmlchange RESUBMIT={resub}')
   if 'queue'    in locals(): run_cmd(f'./xmlchange JOB_QUEUE={queue}')
   if 'walltime' in locals(): run_cmd(f'./xmlchange JOB_WALLCLOCK_TIME={walltime}')
   run_cmd(f'./xmlchange CHARGE_ACCOUNT={acct},PROJECT={acct}')
   #----------------------------------------------------------------------------
   # Submit the run
   run_cmd(f'./case.submit -a=\"--reservation={reservation}\"')
#---------------------------------------------------------------------------------------------------
# Print the case name again
print(f'\n  case : {case}\n') 
#---------------------------------------------------------------------------------------------------
