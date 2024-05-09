#!/usr/bin/env python3
#---------------------------------------------------------------------------------------------------
import os, datetime, subprocess as sp, datetime
from shutil import copy2
# Make directories created by this script world-readable
os.umask(18) # 18 in decimal is equal to 022 in octal 
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd): print('\n'+clr.GREEN+cmd+clr.END) ; os.system(cmd); return
#---------------------------------------------------------------------------------------------------
# THINGS THAT MAY NEED TO BE EDITED
acct             = 'ntrain6'
reservation      = 'e3sm_day4'
src_dir          = os.getenv('HOME')+'/E3SM/E3SM_SRC0' # PATH TO E3SM SOURCE CODE
init_date,sst_yr = '2023-09-08',2023
init_file_atm    = f'{os.getenv("SCRATCH")}/HICCUP/HICCUP.atm_era5.{init_date}.ne30np4.L80.nc'
init_file_sst    = f'{os.getenv("SCRATCH")}/HICCUP/HICCUP.sst_noaa.{init_date}.nc'
# NO NEED TO EDIT BELOW
#---------------------------------------------------------------------------------------------------
newcase,config,build,submit = False,False,False,False

# newcase      = True
# config       = True
# build        = True
submit       = True

stop_opt,stop_n,resub,walltime = 'ndays',5,0,'0:30:00'

case = '.'.join(['E3SM','2024-E3SM-tutorial-hindcast',init_date])

scratch_root = f'{os.getenv("SCRATCH")}/e3sm_scratch/pm-cpu'
case_root = f'{scratch_root}/{case}'
#---------------------------------------------------------------------------------------------------
print(f'\n  case : {case}\n')
#---------------------------------------------------------------------------------------------------
if newcase :
   cmd = f'{src_dir}/cime/scripts/create_newcase'
   cmd += f' --case {case}'
   cmd += f' --output-root {scratch_root} --script-root {case_root}/case_scripts '
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
   file.write(          ",'PHIS'") # sfc geopotential for TC tracking
   file.close()
   #----------------------------------------------------------------------------
   # Specify start date and SST file for hindcast
   os.system(f'./xmlchange --file env_run.xml  RUN_STARTDATE={init_date}')
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
   if 'reservation' in locals():
      run_cmd(f'./case.submit --batch-args=\"--reservation={reservation}\"')
   else:
      run_cmd(f'./case.submit')
#---------------------------------------------------------------------------------------------------
# Print the case name again
print(f'\n  case : {case}\n') 
#---------------------------------------------------------------------------------------------------
