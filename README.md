## HINDCAST INITIAL CONDITION CREATION UTILITY/PROCESSOR (HICCUP)

This is a tool for creating [E3SM](https://e3sm.org/) initial condition files from reanalysis with 
a focus on modularity and portability.

The tool is used by selecting one of the template scripts, such as:

  `template_scripts/create_initial_condition_from_obs.py`

Make a copy of this in the `user_scripts` directory, which will be ignored by git. This template
script copy then needs to be edited to update paths and configure the list of tasks to fit the
user's needs.The new user script can then be excuted after loading a suitable conda environment
(see Setup Notes).

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
    * [scipy](https://www.scipy.org/) - *needed to fill missing SST data around poles*
    * [cdsapi](https://pypi.org/project/cdsapi/) - *for obtaining ECMWF data*
    * [ftplib](https://docs.python.org/3/library/ftplib.html) - *for obtaining NOAA sst/ice data*

It is convenient to create a conda env that includes all these dependencies:
  ```
  conda create --name hiccup_env -c conda-forge xarray dask pandas numpy scipy netcdf4 hdf5 cdsapi tempest-remap "nco>=5.3.1" 
  ```

After creating the environment it can be activated via:

  `source activate hiccup_env`

Alternatively, you can use the E3SM unified environment

On Perlmutter:

`source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh`

On Chrysalis:

`source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh`

You can then install HICCUP into your python environment by running:
```
   pip install ./
```
which will allow you to import `hiccup` from any directory.

Note that if the user is making changes to the underlying hiccup code, the `./setup.py install` process needs to be repeated for any changes to take effect.

TempestRemap and NCO may already be locally available if you are working on a machine at a super-computing center. They can also be installed manually, but we recommend including them in the hiccup conda environment to avoid conflicts.

The default paths for things like grid files, mapping files, and output data is set to local directories. However, when working on a machine at a super-computering center, like NERSC or OLCF, it is useful to avoid filling up ones home directory with this data, especially for high resolution output data. We recommend creating a folder on scratch space and using this to set file path variables when calling create_hiccup_data().

--------------------------------------------------------------------------------

### Obtaining Input Data

Currently. ERA5 + NOAA SST/ice is the preferred input data option.
To aquire new ERA5 data, be sure "cdsapi" is in your conda environment
and you've sset up your CDS API key in `~/.cdsapirc`.

You can then use the `get_hindcast_data.ERA5.py` tool to obtain a single pair of 
ERA5 pressure level and surface data files with

  `python get_hindcast_data.ERA5.py --start-date=<yyyymmdd> --output-root=<path>`

Alternatively, you can obtain ERA5 files over a range of dates with a specified
hourly frequency with

  `python get_hindcast_data.ERA5.py --start-date=<yyyymmdd> --final-date=<yyyymmdd> --start-hour=<hh> --final-hour=<hh> --data-freq=3h --output-root=<path>`

Note that while the `--output-root` argument is optional, it is recommended to 
make sure this points to a location on a scratch disk with sufficient space 
for large data files.

Similarly, 0.25 degree NOAA OI daily SST and sea ice data can be obtained in
yearly files by using the `get_hindcast_data.NOAA_SSTICE.py` tool with command
line arguments to specify a year, or range of years as follows:

  `python get_hindcast_data.NOAA_SSTICE.py --start-year=<yyyy> --final-year=<yyyy> --output-root=<path>`

For a single year, omit the `--final-year` argument.

--------------------------------------------------------------------------------

### Vertical Grid Files

The current E3SM vertical grid was created through an iterative process 
involving numerous, undocumented, subjective decisions mainly by Phil Rasch 
and Po-Lun Ma who did not document the process, so there is no recipe to 
recreate the grid from scratch. 

A vertical grid file for the L80 grid used by E3SMv3 atmosphere is included in the HICCUP repository.
  
  `files_vert/L80_for_E3SMv3.nc`

In addition to other atmosphere vertical grids for other E3SM configuration.

To create a new vertical coordinate file it must be extracted from a 
pre-existing model data file as follows:

  1. Dump the vertical grid data into a text file using ncdump:
     
     `ncdump -v P0,hyam,hybm,hyai,hybi,lev,ilev <history_file> > vert_coord.txt`

  2. manually edit the file to remove extra header info,
     but keep the general CDL format created by ncdump

  3. Generate a new netcdf file from the edited text file using ncgen:
     
     `ncgen vert_coord.txt -o vert_coord.nc`

Vertical grid information can also be procedurally constructed. Future HICCUP updates will bring in template scripts for modifying existing grids and generating new vertical grids from scratch.

--------------------------------------------------------------------------------

### SST and Sea Ice Initial Conditions

HICCUP can generate a file with SST and sea ice data that matches the format that E3SM expects when running a hindcast with prescribed ocean/ice conditions. NOAA OI data is currently the only supported option for this, but HICCUP can easily use ERA5 data. 

Several options are implemented in the `sstice_slice_and_remap()` routine to control how the time coordinate of this data is handled:
```
time_slice_method='match_atmos'   match the time coordinate with the date of the atmospheric initial condition
time_slice_method='initial'       use the first time index of the SST and sea ice data
time_slice_method='use_all'       remap all times provided for transient SSTs
```

For more information on the difference between these approaches see this wiki page => [Fixed vs. Transient SST](https://github.com/E3SM-Project/HICCUP/wiki/Fixed-vs.-Transient-SST)

The first two methods will yield a simulation with SST and sea conditions that are "fixed" at the time of initialization, while the third option provides a simple way to produce a simulation with transient SST and sea ice conditions. 

We plan to implement other methods in the future for handling the time of SST data 
to be more flexible for hig-res runs, such as specifying a specific window of 
SST/ice data to remap and include in the file output file. 

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

See the [wiki article here](https://github.com/E3SM-Project/HICCUP/wiki/Running-a-Hindcast-with-HICCUP-Initial-Conditions)

--------------------------------------------------------------------------------

### Plotting Initial Condition Data

#### ncvis

The [ncvis](https://github.com/SEATStandards/ncvis) tool is a great way to visualize unstructured data, and can help provide a sanity check of HICCUP generated initial condition data.

#### PyNGL

A plotting script is also included (plot.sanity_check.py), but it requires 
PyNGL (https://www.pyngl.ucar.edu/) to be installed in the python environment.
This was done becase PyNGL has excellent support for plotting data on 
unstructured grids. However, PyNGL has been put into "maintenance mode", so in 
the future we need to change these scripts to use an alternative.

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

### Testing

For simple testing of HICCUP functionality the repo includes low-resolution test data from ERA5 and NOAA in the `test_data` folder. These files are used by the `test_scripts/test.*` to exercise the typical HICCUP workflow for generating model input data from observation data and reanalysis. There are also remapping scripts that can be used to regenerate the low-res test data.

To run all unit tests simply type `python test_scripts/unit_test_all.py`. 

--------------------------------------------------------------------------------

### Development Plans

Below are some ideas for future HICCUP enhancments:

- **Fix EAMxx/SCREAM support** - The current workflow is tailored for EAM, and the user needs to use the conversion script in the SCREAM_utils folder to rename and reshape the data in the initial condition file. This step can be completely comitted by adding some special logic to the EAMxx HICCUP data class to do this renaming and reshaping. 
- **Detect whether input files are packed** - The ERA5 data from CDS come "packed" and the NCO unpacking command takes quite a while, even for small files that are already unpacked. The current workflow requires the user to know the state of the input files, so a method for automatically checking whether unpacking needs to be done would be very helpful. This would simplify the workflow a bit because we could delete the flag and line for this unpacking step and just have it done when the hiccup_data object is created.
- **Add simple run script templates for hindcasts and land spin-up** - I have many scripts for these things that are much simpler than the standard [monolithic] E3SM run script for production coupled runs, but they are specific to individual machines and file systems. A more general and simplified script for this would be helpful. 
- **Validation data framework?** - Currently there are few python scripts for obtaining ERA5 validation data, but this could be improved by polishing a separate workflow for this, perhaps including an automated calculation of forecast error.
- **Fix support for using ERA5 SST and sea-ice** - I forget what the issue was here, but this has been requested a few times and I think it would be a valuable feature to have.
- **Add support for CFSR / GFS / MERRA / JRA55** - This has proven difficult due to the ways these datasets are organized. ERA5 offers a lot of flexibilty to facilitate an automated workflow, but other datasets have a single format that must be accomodated. For example, if files are only offered as one variable per file with multiple time steps then a user who needs a single initial condition file at 00Z will have to download orders of magnitude more data than they, and the HICCUP back-end will require special exceptions for how to load each dataset and how the input arguments are sturctured, which will also require many more specialized checks to ensure the data is self-consistent, which seems error-prone. 


--------------------------------------------------------------------------------

