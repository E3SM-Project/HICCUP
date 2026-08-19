## HINDCAST INITIAL CONDITION CREATION UTILITY/PROCESSOR (HICCUP)

This is a tool for creating [E3SM](https://e3sm.org/) initial condition files from reanalysis with 
a focus on modularity and portability.

The tool is used by selecting one of the template scripts, such as:

  `template_hiccup_scripts/create_EAM_IC_from_ERA5-NOAA.py`

Make a copy of this in the `user_scripts` directory, which will be ignored by git. This template
script copy then needs to be edited to update paths and configure the list of tasks to fit the
user's needs. The new user script can then be executed after loading a suitable conda environment
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
  - [Utility Scripts](#utility-scripts)
  - [Hindcast Analysis and Verification](#hindcast-analysis-and-validation)
  - [Testing](#testing)
  - [Development Plans](#development-plans)

--------------------------------------------------------------------------------

### Setup Notes

Dependencies:
  * [NCO](http://nco.sourceforge.net/)
  * [TempestRemap](https://github.com/ClimateGlobalChange/tempestremap)
  * Python modules:
    * [xarray](http://xarray.pydata.org/en/stable/) - *primary data manipulation tool*
    * [numpy](https://numpy.org/) - *numerical operators*
    * [dask](https://www.dask.org/) - *necessary for large grids*
    * [pandas](https://pandas.pydata.org/) - *helpful for handling time coordinate*
    * [netcdf4](https://unidata.github.io/netcdf4-python/) - *needed for writing netcdf4 files with xarray*
    * [hdf5](https://www.h5py.org/) - *needed for netcdf4 format - important for fine grids like ne1024*
    * [scipy](https://www.scipy.org/) - *needed to fill missing SST data around poles*
    * [cdsapi](https://pypi.org/project/cdsapi/) - *for obtaining ECMWF data*

It is convenient to create a conda env that includes all these dependencies:
  ```
  conda env create -f environment.yml
  ```

The plotting scripts in `utilities/` (e.g. `plot.sanity_check.py`) additionally require [matplotlib](https://matplotlib.org/) and [cartopy](https://scitools.org.uk/cartopy/) to create map plots. If you want to make plots, use `conda install -n hiccup_env matplotlib cartopy` to add these libraries to conda env created above.

After creating the environment it can be activated via:

  `conda activate hiccup_env`

Finally, install HICCUP into your python environment by running:
```
pip install -e ./
```
which will allow you to import `hiccup` from any directory.

#### Using the E3SM Unified Env

Alternatively, you can use the E3SM unified environment

On Perlmutter:

`source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh`

On Chrysalis:

`source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh`

Don't forget to install HICCUP into your environment via:
```
pip install -e ./
```

TempestRemap and NCO may already be locally available if you are working on a machine at a HPC center. They can also be installed manually, but we recommend including them in the hiccup conda environment to avoid conflicts.

The default paths for things like grid files, mapping files, and output data are set to local directories. However, when working on a machine at a HPC center, like NERSC or OLCF, it is useful to avoid filling up one's home directory with this data, especially for high resolution output data. We recommend creating a folder on scratch space and using this to set file path variables when calling create_hiccup_data().

--------------------------------------------------------------------------------

### Obtaining Input Data

Currently, ERA5 + NOAA SST/ice is the preferred input data option.
To obtain new ERA5 data, be sure "cdsapi" is in your conda environment
and you've set up your CDS API key in `~/.cdsapirc`.

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

A vertical grid file for the L80 grid used by E3SMv3 atmosphere is included in the HICCUP repository:
  
  `files_vert/L80_for_E3SMv3.nc`

Other atmosphere vertical grids are included for other configurations, or are experimental. Use at own risk!

#### Creating New Vertical Grid files

Vertical grid files can be procedurally constructed, but HICCUP does not contain examples of this. Future updates may bring in template scripts for generating new vertical grids or modifying old ones. If any users are interested in this they can reach out to Walter Hannah for examples.

Vertical coordinate information can be extracted from a pre-existing data file to be used for new HICCUP initial conditions as follows:

  1. Dump the vertical grid data into a text file using ncdump:
     
     `ncdump -v P0,hyam,hybm,hyai,hybi,lev,ilev <history_file> > vert_coord.txt`

  2. manually edit the file to remove extra header info,
     but keep the general CDL format created by ncdump

  3. Generate a new netcdf file from the edited text file using ncgen:
     
     `ncgen vert_coord.txt -o vert_coord.nc`

--------------------------------------------------------------------------------

### SST and Sea Ice Initial Conditions

HICCUP can generate a file with SST and sea ice data that matches the format that E3SM expects when running a hindcast with prescribed ocean/ice conditions. NOAA OI data is currently the only supported option for this.

Several options are implemented in the `sstice_slice_and_remap()` routine to control how the time coordinate of this data is handled:
```
time_slice_method='match_atmos'   match the time coordinate with the date of the atmospheric initial condition
time_slice_method='initial'       use the first time index of the SST and sea ice data
time_slice_method='use_all'       remap all times provided for transient SSTs
```

For more information on the difference between these approaches see this wiki page => [Fixed vs. Transient SST](https://github.com/E3SM-Project/HICCUP/wiki/Fixed-vs.-Transient-SST)

The first two methods will yield a simulation with SST and sea ice conditions that are "fixed" at the time of initialization, while the third option provides a simple way to produce a simulation with transient SST and sea ice conditions. 

We plan to implement other methods in the future for handling the time of SST data 
to be more flexible for high-res runs, such as specifying a specific window of 
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

After the input data is obtained, HICCUP can be used to generate initial conditions 
by editing and running one of the scripts in `template_hiccup_scripts`, such as 
`template_hiccup_scripts/create_EAM_IC_from_ERA5-NOAA.py`. This script orchestrates
the workflow for generating the atmosphere initial condition as well as the
SST/sea-ice data file.

The HICCUP workflow centers on a "hiccup_data" object that carries the information 
needed for processing the data as well as class methods for processing the data. 
There is also a python dictionary of temporary file names that are used to store 
the data for each variable during processing. This approach of separating the 
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

#### ncview

The [ncview](https://cirrus.ucsd.edu/ncview/) tool works well for visualizing atmospheric data, but it is limited to rectilinear grids, such as the common equiangular lat/lon grids. Therefore, in order to use `ncview` the user must remap the unstructured HICCUP initial condition data onto a suitable grid.

##### How to remap HICCUP data for visualization

The `ncremap` tool within the [NCO](http://nco.sourceforge.net/) library is capable of creating these grid files and remapping the data. The documentation ([https://nco.sourceforge.net/nco.html](https://nco.sourceforge.net/nco.html)) is a good resource for understanding how to do this.

Here's an example of what the remapping task might look like:

```shell
NE=30
SRC_GRID=ne${NE}pg2

DST_NY=180
DST_NX=360
DST_GRID=${DST_NY}x${DST_NX}

GRID_ROOT=???
MAPS_ROOT=???

SRC_GRID_FILE=${GRID_ROOT}/${SRC_GRID}_scrip.nc
DST_GRID_FILE=${GRID_ROOT}/${DST_GRID}_scrip.nc
MAP_FILE=${MAPS_ROOT}/map_${SRC_GRID}_to_${DST_GRID}_traave.nc

# generate model grid file
GenerateCSMesh --alt --res ${NE} --file ${GRID_ROOT}/ne${NE}.g
GenerateVolumetricMesh --in ${GRID_ROOT}/ne${NE}.g --out ${GRID_ROOT}/ne${NE}pg2.g --np 2 --uniform
ConvertMeshToSCRIP --in ${GRID_ROOT}/ne${NE}pg2.g --out ${GRID_ROOT}/ne${NE}pg2_scrip.nc

# generate destination (lat/lon) grid file
ncremap -g ${DST_GRID_FILE} -G ttl="Equi-Angular grid, dimensions ${DST_GRID}, cell edges on Poles/Equator and Prime Meridian/Date Line"#latlon=${DST_NY},${DST_NX}#lat_typ=uni#lon_typ=grn_wst

# generate map file
ncremap -6 --alg_typ=traave --grd_src=$SRC_GRID_FILE --grd_dst=$DST_GRID_FILE --map=$MAP_FILE

# apply the map to a HICCUP file
ncremap -m $MAP_FILE hiccup_ic.nc hiccup_ic_${DST_GRID}.nc
```

#### ncvis

The [ncvis](https://github.com/SEATStandards/ncvis) tool is a great way to visualize unstructured data, and can help provide a sanity check of HICCUP generated initial condition data. Although, some users have reported difficulties installing this into their conda env.

#### Included Plotting Scripts

A plotting script is included in the repo (`utilities/plot.sanity_check.py`) that allows for customizable plots of the final HICCUP initial condition data.

--------------------------------------------------------------------------------

### Utility Scripts

Various useful scripts are included in `utilities/` that can help identify issues in the initial condition files that HICCUP produces.

* `utilities/chk.data.py <file>` - This script checks for invalid values, such as `Inf` and `NaN`.

* `utilities/chk.stats.py  <file>` - This script is currently outdated and needs an update, but the idea was to iterate over a fixed list of important variables and print summary statistics to check for unreasonably high or low values.

* `utilities/find_bad_values.py --file=<file>` - Similar to `chk.stats.py`, this script was created when debugging a curious issue where the surface adjustment step was creating odd patterns over ocean regions. It was so helpful that I decided to leave it in the repo.

--------------------------------------------------------------------------------

### Hindcast Analysis and Validation

The task of analyzing the hindcast output data is up to user, although 
HICCUP may incorporate some simple skill/error metrics in the future. For now, we have 
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

Below are issues that will be addressed by future development:

- **Replace NCO commands with python for better scalability** - As we move towards routinely running global cloud resolving simulations we have found that various parts of HICCUP do not scale well. These issues have been addressed in an iterative fashion, but as more of this type of issue comes up it seems that a redesign might be helpful to make HICCUP more robust. NCO will still be required for many tasks, but I think many things could be streamlined and improved if they were done in python/xarray. This would also allow more careful control of the memory footprint for very large grids.
- **Need to expand RRM wiki pages** - There is currently a place-holder page for this, but it needs to be fleshed out to describe the unique challenges that come with RRM cases.
- **Add simple template run scripts for E3SM hindcasts and land spin-up** - I have many scripts for these things that are much simpler than the standard [monolithic] E3SM run script for production coupled runs, but they are specific to individual machines and file systems. A more general and simplified script for this could be helpful for new users.
- **Fix support for using ERA5 SST and sea-ice** - I forget what the issue was here, but this has been requested a few times and I think it would be a valuable feature to have.
- **Add support for ERA5 model level data** - I haven't looked into this enough, but the few times I tried to get this working revealed some limitations related to the vertical interpolation. This has been requested many times, but it does not appear to be a trivial effort to implement.
- **Add support for CFSR / GFS / MERRA / JRA55** - This has proven difficult due to the ways these datasets are organized. ERA5 offers a lot of flexibility to facilitate an automated workflow, but other datasets have a single format that must be accommodated. For example, if files are only offered as one variable per file with multiple time steps then a user who needs a single initial condition file at 00Z will have to download orders of magnitude more data than they need, and the HICCUP back-end will require special exceptions for how to load each dataset and how the input arguments are structured, which will also require many more specialized checks to ensure the data is self-consistent, which seems error-prone. 


--------------------------------------------------------------------------------

