HINDCAST INITIAL CONDITION CREATION UTILITY/PROCESSOR (HICCUP)

This is a tool for creating E3SM initial condition files from reanalysis with 
a focus on simplicity and portability.

The initial design discussion can be found here:
https://confluence.exascaleproject.org/display/ADSE15/Creating+Hindcast+Initial+Conditions

The tool is used by editing and running:
  create_initial_condition.py

--------------------------------------------------------------------------------

SETUP NOTES

Dependencies:
  NCO
  TempestRemap
  Python modules:
    xarray
    pandas
    numpy
    scipy
    netcdf4
    cdsapi (for ECMWF data)
    ftplib (for NOAA sst/ice data)

It is useful to create a conda environment that includes the python dependencies.
This can be created with the following command:

  conda create --name hiccup_env -c conda-forge xarray pandas scipy netcdf4 tempest-remap nco cdsapi 

After creating the env it can be activated via:
  source activate hiccup_env
  
To exit the env simply type: conda deactivate

TempestRemap and NCO may be locally available if you are working with a machine
at a super-computing center (such as NERSC). 

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

OBTAINING INPUT DATA

Currently, ERA5 realanysis and NOAA SST is the only supported input data option.
To aquire new ERA5 data, be sure to use conda to install the "cdsapi" module 
and set up your ECMWF API key in ~/.ecmwfapirc,then edit and run:
  get_ERA5_data.py

To aquire NOAA OI daily SST and sea ice data, edit and run:
  get_NOAA_SST+ICE_data.py

--------------------------------------------------------------------------------

PLOTTING THE OUTPUT INITIAL CONDITION FILE(S)

A plotting script is also included (plot.sanity_check.py), but it requires PyNGL
to be installed in the python environment. This was done becase PyNGL has 
excellent support for plotting data on unstructured grids. In the future we hope 
to add another plotting script that uses MatPlotLib.

--------------------------------------------------------------------------------

VERTICAL GRID FILES

The current E3SM vertical grid was created through an iterative process 
involving numerous, undocumented, subjective decisions mainly by Phil Rasch 
and Po-Lun Ma who did not document the process, so there is no recipe to 
recreate the grid from scratch. To create the vertical coordinate file it is 
easiest to extract it from a pre-existing model data file as follows:

  1. Dump the vertical grid data into a text file using ncdump:
     ncdump -v P0,hyam,hybm,hyai,hybi,lev,ilev <history_file> > vert_coord.txt

  2. manually edit the file to remove extra header info,
     but keep the general CDL format created by ncdump

  3. Generate a new netcdf file from the edited text file using ncgen:
     ncgen vert_coord.txt -o vert_coord.nc

--------------------------------------------------------------------------------

SST AND SEA ICE INITIAL CONDITONS

HICCUP can also generate a data file with SST and sea ice data. NOAA OI data is
typically used for this, but HICCUP also currently supports using ERA5 data. 
The NOAA data comes in yearly files of daily averages, and the default behavior 
of the sstice_slice_and_remap() method is to find the time that matches the 
provided atmosphere initial condition file. We plan to implement other methods
for handling the time of SST data, such as producing a weekly average centered 
on the initialization day. 

--------------------------------------------------------------------------------

DATA FOR TESTING AND DEVELOPMENT

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