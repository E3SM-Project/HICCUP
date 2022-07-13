#!/usr/bin/env python
#===================================================================================================
#  Nov, 2017 - Walter Hannah - Lawrence Livermore National Lab
#  This script prints the tail end of the most recently modified HICCUP output log
#===================================================================================================
import sys, os, glob, errno
import subprocess as sp
import xarray as xr
import numpy as np
from optparse import OptionParser

help_header = 'usage: ./%prog [file] [file] ...\n'
help_header += '\nThis script will check the data of the input files for inf and nan values'
parser = OptionParser(usage=help_header)
parser.add_option('--vars',dest='vars',default='all',help='comma separated list of variables to check')
(opts, args) = parser.parse_args()

# Set up terminal colors
class bcolor:
    ENDC     = '\033[0m';   BLACK    = '\033[30m'
    RED      = '\033[31m';  GREEN    = '\033[32m'
    YELLOW   = '\033[33m';  BLUE     = '\033[34m'
    MAGENTA  = '\033[35m';  CYAN     = '\033[36m'
    WHITE    = '\033[37m';  BOLD     = '\033[1m'

indent = ' '*4

#-------------------------------------------------------------------------------
# Simple routine for chcecking variable values - useful for debugging
#-------------------------------------------------------------------------------
def print_stat(x,name='(no name)',unit='',fmt='f',stat='naxh',indent='  '):
  """ 
  Simple routine for printing various statistics or properites of a variable.
  The characters of the "stat" string variable are used to specify the order 
  and type of quantities to calculate
    n   minimum value
    a   average across all dimensions
    x   maximum value
    s   standard deviation
    h   shape
  """
  if fmt=='f' : fmt = '%f'
  if fmt=='e' : fmt = '%e'
  if unit!='' : unit = '['+str(unit)+']'
  print('\n'+indent+name+' '+unit)
  for c in list(stat):
      if c=='n' : print(indent+'min: '+fmt%x.min() )
      if c=='a' : print(indent+'avg: '+fmt%x.mean())
      if c=='x' : print(indent+'max: '+fmt%x.max() )
      if c=='s' : print(indent+'std: '+fmt%x.std() )
      if c=='h' : print(indent+'shp: '+str(x.shape) )
#===================================================================================================
#===================================================================================================

if len(args)<1: exit('No file name provided!')

print('\nChecking files for invalid values...')

# This method of handling file name arguments allows for wildcards
files = []
for arg in args:
    files.extend(glob.glob(arg))

for file_name in files :

    print('\n  '+bcolor.BOLD+file_name+bcolor.ENDC)

    ds = xr.open_dataset(str(file_name))
    # print(ds)
    # continue

    if opts.vars=='all':
      var_list = ds.data_vars
    else:
      var_list = opts.vars.split(',')

    for var in var_list: 

        # Print stats, but skip time related and other variables
        if var not in ['lat_vertices','lon_vertices','time_bnds','area'] and ds[var].dims!=('time',):
            print_stat(ds[var],name=var,indent=indent,stat='nxh')

        inf_cnt = ds[var].where( xr.ufuncs.isinf(ds[var]) ).count().values
        nan_cnt = ds[var].where( xr.ufuncs.isnan(ds[var]) ).count().values

        if inf_cnt>0: print(f'{indent}{var}: inf values found! ({inf_cnt})')
        if nan_cnt>0: print(f'{indent}{var}: inf values found! ({nan_cnt})')

print('\ndone.\n')

    

#===================================================================================================
#===================================================================================================
