Hindcast Initial Condition Creation Utility/Processor (HICCUP)

This is a tool for creating E3SM initial condition files from reanalysis with 
a focus on simplicity and portability.

The initial design discussion can be found here:
https://confluence.exascaleproject.org/display/ADSE15/Creating+Hindcast+Initial+Conditions

The tool is used by editing and running:
  create_initial_condition.py

Currently, ERA5 realanysis is the only supported option for input data.
To aquire new ERA5 data, be sure to use conda to install the "cdsapi" module 
and set up your ECMWF API key in ~/.ecmwfapirc,then edit and run:
  get_ERA5_data.py

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

DATA FOR TESTING AND DEVELOPMENT
a low-resolution version of ERA5 pressure level data is included in this repo:
  HICCUP_TEST.ERA5.atm.low-res.nc
  HICCUP_TEST.ERA5.sfc.low-res.nc
To aquire new test data, use the get_ERA5_data.py script, followed by running:
  remap.ERA5_test_data.py
which is a simple script for reducing the resolution of the test data