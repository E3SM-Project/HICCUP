HINDCAST INITIAL CONDITION CREATION UTILITY/PROCESSOR (HICCUP)

This is a tool for creating E3SM initial condition files from reanalysis with 
a focus on simplicity and portability.

The tool is used by editing and running:
  create_initial_condition.py

--------------------------------------------------------------------------------

TABLE OF CONTENTS
  - [SETUP NOTES](#setup-notes)
  - [OBTAINING INPUT DATA](#obtaining-input-data)
  - [VERTICAL GRID FILES](#)
  - [SST AND SEA ICE INITIAL CONDITONS](#)
  - [LAND MODEL INITIAL CONDITONS](#)
  - [RUNNING A HINDCAST](#)
  - [PLOTTING THE OUTPUT INITIAL CONDITION FILE(S)](#)
  - [HINDCAST ANALYSIS AND VALIDATION](#)
  - [DATA FOR TESTING AND DEVELOPMENT](#)

--------------------------------------------------------------------------------

# SETUP NOTES

Dependencies:
  NCO
  TempestRemap
  Python modules:
    xarray  - primary data manipulation tool
    pandas  - helpful for handling time coordinate
    netcdf4 - needed for writing netcdf4 files with xarray
    hdf5    - needed for netcdf4 format - important for fine grids like ne1024
    pynio   - needed for reading grib files in the case of CFS data
    scipy   - needed to fill missing SST data around poles
    cdsapi  - for obtaining ECMWF data
    ftplib  - for obtaining NOAA sst/ice data

It is convenient to create a conda env that includes all these dependencies:
  conda create --name hiccup_env -c conda-forge xarray dask pandas scipy netcdf4 pynio hdf5 cdsapi tempest-remap nco   

After creating the environment it can be activated via:
  source activate hiccup_env

TempestRemap and NCO may already be locally available if you are working on 
a machine at a super-computing center. They can also be installed manually.

To install NCO manually:
  TempestRemap can be easily installed as part of the conda environment above,
  but it can also be installed manually. If using Mac OSX then we recommend 
  using homebrew to install NCO (see https://brew.sh/)
  Otherwise installation information can be found at http://nco.sourceforge.net/

To install TempestRemap manually:
  TempestRemap can be easily installed as part of the conda environment above,
  But it can aslo be downloaded and built from a public github repository
  (https://github.com/ClimateGlobalChange/tempestremap)
  using the following commands:
    git clone https://github.com/paullric/tempestgecore.git
    <edit the Makefile to customize the NetCDF paths>
    make -f Makefile.gmake all

--------------------------------------------------------------------------------

# OBTAINING INPUT DATA

Currently, ERA5 + NOAA SST/ice is the only supported input data option.
To aquire new ERA5 data, be sure "cdsapi" is in your conda environment
and set up your ECMWF API key in ~/.ecmwfapirc,then edit and run:
  get_hindcast_data.ERA5.py

To aquire NOAA OI daily SST and sea ice data, edit and run:
  get_hindcast_data.NOAA_SSTICE.py

--------------------------------------------------------------------------------

# VERTICAL GRID FILES

The current E3SM vertical grid was created through an iterative process 
involving numerous, undocumented, subjective decisions mainly by Phil Rasch 
and Po-Lun Ma who did not document the process, so there is no recipe to 
recreate the grid from scratch. 

A vertical grid file for the L72 grid is included in the HICCUP repository.
  files_vert/vert_coord_L72.nc

To create a new vertical coordinate file it must be extracted from a 
pre-existing model data file as follows:

  1. Dump the vertical grid data into a text file using ncdump:
     ncdump -v P0,hyam,hybm,hyai,hybi,lev,ilev <history_file> > vert_coord.txt

  2. manually edit the file to remove extra header info,
     but keep the general CDL format created by ncdump

  3. Generate a new netcdf file from the edited text file using ncgen:
     ncgen vert_coord.txt -o vert_coord.nc

--------------------------------------------------------------------------------

# SST AND SEA ICE INITIAL CONDITONS

HICCUP can also generate a data file with SST and sea ice data. NOAA OI data is
typically used for this, but HICCUP also currently supports using ERA5 data. 
The NOAA data comes in yearly files of daily averages, and the default behavior 
of the sstice_slice_and_remap() method is to find the time that matches the 
provided atmosphere initial condition file. We plan to implement other methods
for handling the time of SST data, such as producing a weekly average centered 
on the initialization day. 

--------------------------------------------------------------------------------

# LAND MODEL INITIAL CONDITONS

HICCUP does not currently support the generation of land model initial condition
files. This might be possible with the data available from ERA5, but the current
recommendation is to spin up the land model for 5-10 years leading up to the 
desired initialization date using the standalone land model forced by the data
atmosphere component. 

--------------------------------------------------------------------------------

# RUNNING A HINDCAST

After using HICCUP to generate the atmosphere and SST/ice files, an E3SM 
hindcast can be run by following the steps to run a typical "F-compset" run, 
using the compsets such as FC5AV1C-L. The initialization data needs to be copied 
to the scratch space of the machine to ensure they are accessible to the compute 
nodes. 

The atmospheric initial condition file is specified by editing the "user_nl_cam"
file found in the case directory to include:
  
  ncdata = < path to hiccup atmos initial condition file >

The SST file and start date values also need to be specified by modifying the 
env_run.xml file in the case directory. The preferred method for doing this is 
to use the xmlchange command from the case directory as in the example below:

  ./xmlchange SSTICE_DATA_FILENAME=<path to SST file>
  ./xmlchange RUN_STARTDATE=2016-08-01
  ./xmlchange SSTICE_YEAR_ALIGN=2016
  ./xmlchange SSTICE_YEAR_START=2016
  ./xmlchange SSTICE_YEAR_END=2017

If using a python script to run the hindcast, here's a snippet of code that 
does the modifications described above:

  ################################################
  # python code to setup hindcast files
  ################################################
  iyr,imn,idy = 2016,8,1
  init_date = f'{iyr}-{imn:02d}-{idy:02d}'
  init_file_atm = f'<path-to-hiccup-data>/HICCUP.atm_era5.{init_date}.ne30np4.L72.nc'
  init_file_sst = f'<path-to-hiccup-data>/HICCUP.sst_noaa.{init_date}.nc'
  os.chdir(case_directory)
  os.system(f'./xmlchange RUN_STARTDATE={init_date}')
  os.system(f'./xmlchange SSTICE_DATA_FILENAME={init_file_sst}')
  os.system(f'./xmlchange SSTICE_YEAR_ALIGN={iyr}')
  os.system(f'./xmlchange SSTICE_YEAR_START={iyr}')
  os.system(f'./xmlchange SSTICE_YEAR_END={iyr+1}')

  file = open('user_nl_cam','a') 
  file.write(f' ncdata = \'{init_file_atm}\'\n')
  file.close()
  ################################################
  ################################################

--------------------------------------------------------------------------------

# PLOTTING THE OUTPUT INITIAL CONDITION FILE(S)

A plotting script is also included (plot.sanity_check.py), but it requires 
PyNGL (https://www.pyngl.ucar.edu/)to be installed in the python environment.
This was done becase PyNGL has excellent support for plotting data on 
unstructured grids. In the future we may add a plotting script for MatPlotLib.

--------------------------------------------------------------------------------

# HINDCAST ANALYSIS AND VALIDATION

The task of analyzing the hindcast output data is up to user for now, although 
we may include some simple skill/error metrics in the future. For now, we have 
included a few simple scripts for obtaining and remapping ERA5 validation data.
  get_validation_data.ERA5.py
  remap.validation_data.ERA5.py

These scripts are configured to obtain a set of atmospheric fields on common 
pressure levels, like U200 and Z500, that are typically used for calculating 
forecast skill. The remap script is configured to put the data on a relatively 
coarse 2 degree grid in order to simplify the calculation of global metrics. 

--------------------------------------------------------------------------------

# DATA FOR TESTING AND DEVELOPMENT

a low-resolution version of ERA5 pressure level data is included in this repo:
  HICCUP_TEST.ERA5.atm.low-res.nc
  HICCUP_TEST.ERA5.sfc.low-res.nc
To aquire new test data, use the get_ERA5_data.py script, follwed by unpacking 
the data with the following command:
  ncpdq -U HICCUP_TEST.ERA5.atm.nc HICCUP_TEST.ERA5.atm.upack.nc
  ncpdq -U HICCUP_TEST.ERA5.sfc.nc HICCUP_TEST.ERA5.sfc.upack.nc
followed by running (check to make sure file names match):
  remap.test_data.ERA5.py
which is a simple script for reducing the resolution of the test data

--------------------------------------------------------------------------------
