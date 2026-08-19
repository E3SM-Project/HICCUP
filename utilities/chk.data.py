#!/usr/bin/env python
#===================================================================================================
#  This script is for inspecting the input data for HICCUP and check for invalid values
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
def print_stat(x,name='(no name)',unit='',fmt='f',stat='naxh',indent='  ',compact=True):
  """ Print min, avg, max, and std deviation of input """
  if fmt=='f' : fmt = '%20.10f'
  if fmt=='e' : fmt = '%e'
  if unit!='' : unit = f'[{unit}]'
  name_len = 25 if compact else len(name)
  msg = ''
  line = f'{indent}{name:{name_len}} {unit}'
  # if not compact: print(line)
  if not compact: msg += line+'\n'
  for c in list(stat):
    if not compact: line = indent
    if c=='h' : line += '   shp: '+str(x.shape)
    if c=='a' : line += '   avg: '+fmt%x.mean()
    if c=='n' : line += '   min: '+fmt%x.min()
    if c=='x' : line += '   max: '+fmt%x.max()
    if c=='s' : line += '   std: '+fmt%x.std()
    # if not compact: print(line)
    if not compact: msg += line+'\n'
  # if compact: print(line)
  if compact: msg += line#+'\n'
  print(msg)
  return msg
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
        if var not in ['lat_vertices','lon_vertices','time_bnds','area','hyam','hybm','hyai','hybi'] and ds[var].dims!=('time',):
          print_stat(ds[var],name=var,indent=indent,stat='nxh')

          inf_cnt = ds[var].where( np.isinf(ds[var]) ).count().values
          nan_cnt = ds[var].where( np.isnan(ds[var]) ).count().values

          if inf_cnt>0: print(f'{indent}{var}: inf values found! ({inf_cnt})')
          if nan_cnt>0: print(f'{indent}{var}: nan values found! ({nan_cnt})')

print('\ndone.\n')

    

#===================================================================================================
#===================================================================================================
