## HINDCAST INITIAL CONDITION CREATION UTILITY/PROCESSOR (HICCUP)

This is a tool for creating [E3SM](https://e3sm.org/) initial condition files from reanalysis with 
a focus on simplicity and portability.

The tool is used by editing and running one of the template scripts:

  `template_scripts/create_initial_condition_from_obs.py`

--------------------------------------------------------------------------------

### TABLE OF CONTENTS
  - [Setup Notes](#setup-notes)
  - [Obtaining Input data](#obtaining-input-data)
  - [Vertical Grid Files](#vertical-grid-files)
  - [SST and Sea Ice Initial Conditions](#sst-and-sea-ice-initial-conditions)
  - [Land Model Initial Conditions](#land-model-initial-conditions)
  - [Generating HICCUP Initial Conditions](#generating-hiccup-initial-conditions)
  - [Running a Hindcast](#running-a-hindcast)
  - [Plotting Initial Condition Data](#plotting-initial-condition-data)
  - [Hindcast Analysis and Verification](#hindcast-analysis-and-validation)
  - [Data for Testing and Development](#data-for-testing-and-development)
  - [Development Plans](#development-plans)

--------------------------------------------------------------------------------

### Setup Notes

Dependencies:
  * [NCO](http://nco.sourceforge.net/)
  * [TempestRemap](https://github.com/ClimateGlobalChange/tempestremap)
  * Python modules:
    * [xarray](http://xarray.pydata.org/en/stable/) - *primary data manipulation tool*
    * [pandas](https://pandas.pydata.org/) - *helpful for handling time coordinate*
    * [netcdf4](https://unidata.github.io/netcdf4-python/) - *needed for writing netcdf4 files with xarray*
    * [hdf5](https://www.h5py.org/) - *needed for netcdf4 format - important for fine grids like ne1024*
    * [pynio](https://www.pyngl.ucar.edu/Nio.shtml) - *needed for reading grib files in the case of CFS data*
    * [scipy](https://www.scipy.org/) - *needed to fill missing SST data around poles*
    * [cdsapi](https://pypi.org/project/cdsapi/) - *for obtaining ECMWF data*
    * [ftplib](https://docs.python.org/3/library/ftplib.html) - *for obtaining NOAA sst/ice data*

It is convenient to create a conda env that includes all these dependencies:
  ```
  conda create --name hiccup_env -c conda-forge xarray dask pandas numpy scipy netcdf4 pynio hdf5 cdsapi tempest-remap nco   
  ```

After creating the environment it can be activated via:

  `source activate hiccup_env`

You can (optionally) install HICCUP into your python environment by running `setup.py`, for example:
```
   ./setup.py install
```
which will allow you to import `hiccup` from any directory.

TempestRemap and NCO may already be locally available if you are working on a machine at a super-computing center. They can also be installed manually, but we recommend including them in the hiccup conda environment to avoid conflicts.

The default paths for things like grid files, mapping files, and output data is set to local directories. However, when working on a machine at a super-computering center, like NERSC or OLCF, it is useful to avoid filling up ones home directory with this data, especially for high resolution output data. We recommend creating a folder on scratch space and using this to set file path variables when calling create_hiccup_data().

--------------------------------------------------------------------------------

### Obtaining Input Data

Currently, ERA5 + NOAA SST/ice is the only supported input data option.
To aquire new ERA5 data, be sure "cdsapi" is in your conda environment
and set up your ECMWF API key in `~/.ecmwfapirc`,then edit and run:

  `get_hindcast_data.ERA5.py`

To aquire NOAA OI daily SST and sea ice data, edit and run:

  `get_hindcast_data.NOAA_SSTICE.py`

--------------------------------------------------------------------------------

### Vertical Grid Files

The current E3SM vertical grid was created through an iterative process 
involving numerous, undocumented, subjective decisions mainly by Phil Rasch 
and Po-Lun Ma who did not document the process, so there is no recipe to 
recreate the grid from scratch. 

A vertical grid file for the L72 grid is included in the HICCUP repository.
  
  `files_vert/vert_coord_L72.nc`

To create a new vertical coordinate file it must be extracted from a 
pre-existing model data file as follows:

  1. Dump the vertical grid data into a text file using ncdump:
     
     `ncdump -v P0,hyam,hybm,hyai,hybi,lev,ilev <history_file> > vert_coord.txt`

  2. manually edit the file to remove extra header info,
     but keep the general CDL format created by ncdump

  3. Generate a new netcdf file from the edited text file using ncgen:
     
     `ncgen vert_coord.txt -o vert_coord.nc`

--------------------------------------------------------------------------------

### SST and Sea Ice Initial Conditions

HICCUP can also generate a data file with SST and sea ice data. NOAA OI data is
typically used for this, but HICCUP also currently supports using ERA5 data. 
The NOAA data comes in yearly files of daily averages, and the default behavior 
of the sstice_slice_and_remap() method is to find the time that matches the 
provided atmosphere initial condition file. We plan to implement other methods
for handling the time of SST data, such as producing a weekly average centered 
on the initialization day. 

--------------------------------------------------------------------------------

### Land Model Initial Conditions

HICCUP does not currently support the generation of land model initial condition
files. This might be possible with the data available from ERA5, but the current
recommendation is to spin up the land model for 5-10 years leading up to the 
desired initialization date using the standalone land model forced by the data
atmosphere component. 

--------------------------------------------------------------------------------

### Generating HICCUP Initial Conditions

After the input data is aquired, HICCUP can be used to generate initial conditions 
by editing and running the `template_scripts/create_initial_condition_from_obs.py` 
script. This script controls the workflow for generating the atmosphere initial 
condition as well as the SST/sea-ice data file.

The HICCUP workflow centers on a "hiccup_data" object that carries the information 
needed for processing the data as well as class methods for processing the data. 
There is also a python dictionary of temporary file names that are used to store 
the data for each variable during processing. This appraoch of separating the 
data variables may seem odd, but it is necessary for very large datasets, so it 
was adopted to avoid supporting multiple workflows. Currently, this dict of files
and the final output file is separate from the hiccup_data object, but we are 
considering putting these into the hiccup_data object to simplify the workflow.

HICCUP is designed to be as modular as possible, but the order in which the input 
data are processed is very important. The most important part of this is the 
regridding and surface adjustment sections. The process must start with the 
horizontal regridding, which alters the surface topography and requires an 
adjustment of surface temperature and pressure. The variable renaming and 
adjustment of time and date information is also done after the horizontal 
regridding. The vertical regridding is the last step in this process because it 
must follow the surface adjustment.

--------------------------------------------------------------------------------

### Running a Hindcast

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

  ```
  ./xmlchange SSTICE_DATA_FILENAME=<path to SST file>
  ./xmlchange RUN_STARTDATE=2016-08-01
  ./xmlchange SSTICE_YEAR_ALIGN=2016
  ./xmlchange SSTICE_YEAR_START=2016
  ./xmlchange SSTICE_YEAR_END=2017
  ```

If using a python script to run the hindcast, here's a snippet of code that 
does the modifications described above:

  ```python
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
  ```

--------------------------------------------------------------------------------

### Plotting Initial Condition Data

A plotting script is also included (plot.sanity_check.py), but it requires 
PyNGL (https://www.pyngl.ucar.edu/) to be installed in the python environment.
This was done becase PyNGL has excellent support for plotting data on 
unstructured grids. However, PyNGL has been put into "maintenance mode", so in 
the future we need to change these scripts to use MatPlotLib https://matplotlib.org/ 
or GeoCAT (https://geocat.ucar.edu/).

--------------------------------------------------------------------------------

### Hindcast Analysis and Validation

The task of analyzing the hindcast output data is up to user for now, although 
we may include some simple skill/error metrics in the future. For now, we have 
included a few simple scripts for obtaining and remapping ERA5 validation data.

  `get_validation_data.ERA5.py`

  `remap.validation_data.ERA5.py`


These scripts are configured to obtain a set of atmospheric fields on common 
pressure levels, like U200 and Z500, that are typically used for calculating 
forecast skill. The remap script is configured to put the data on a relatively 
coarse 2 degree grid in order to simplify the calculation of global metrics. 

--------------------------------------------------------------------------------

### Data for Testing and Development

a low-resolution version of ERA5 pressure level data is included in this repo:
  
  `HICCUP_TEST.ERA5.atm.low-res.nc`
  
  `HICCUP_TEST.ERA5.sfc.low-res.nc`

To aquire new test data, use the get_ERA5_data.py script, follwed by unpacking 
the data with the following command:
  ```
  ncpdq -U HICCUP_TEST.ERA5.atm.nc HICCUP_TEST.ERA5.atm.upack.nc
  ncpdq -U HICCUP_TEST.ERA5.sfc.nc HICCUP_TEST.ERA5.sfc.upack.nc
  ```
followed by running (check to make sure file names match):
  
  `remap.test_data.ERA5.py`
  
which is a simple script for reducing the resolution of the test data

--------------------------------------------------------------------------------

### Development Plans

Current plans for HICCUP enhancments and fixes:

- **Detect whether input files are packed** - The ERA5 data from CDS come "packed" and the NCO unpacking command takes quite a while, even for small files that are already unpacked. The current workflow requires the user to know the state of the input files, so a method for automatically checking whether unpacking needs to be done would be very helpful. This would simplify the workflow a bit because we could delete the flag and line for this unpacking step and just have it done when the hiccup_data object is created.
- **Add run scripts for hindcasts and land spin up** - I have a script for this, but it's very specific to a certain machine. A more general and simplified script for this would be helpful. 
- **Validation data framework?** - Currently there are few scripts for obtaining ERA5 validation data, but I think this could be improved and maybe add a whole separate workflow.
- **Fix support for using ERA5 SST and sea-ice**
- **Add support for CFSR / GFS / MERRA**


--------------------------------------------------------------------------------

