#!/usr/bin/env python
import os, datetime, xarray as xr
#-------------------------------------------------------------------------------
def print_stat(x,name='(no name)',unit='',fmt='f',stat='naxh',indent='',compact=False):
   """ Print min, avg, max, and std deviation of input """
   if fmt=='f' : fmt = '%.4f'
   if fmt=='e' : fmt = '%+8.4e'
   if unit!='' : unit = f'[{unit}]'
   name_len = 12 if compact else len(name)
   msg = ''
   line = f'{indent}{name:{name_len}} {unit}'
   # if not compact: print(line)
   if not compact: msg += line+'\n'
   for c in list(stat):
      if not compact: line = indent
      if c=='h' : line += '     shp: '+str(x.shape)
      if c=='a' : line += '     avg: '+fmt%x.mean()
      if c=='n' : line += '     min: '+fmt%x.min()
      if c=='x' : line += '     max: '+fmt%x.max()
      if c=='s' : line += '     std: '+fmt%x.std()
      # if not compact: print(line)
      if not compact: msg += line+'\n'
   # if compact: print(line)
   if compact: msg += line#+'\n'
   print(msg)
   return msg
#-------------------------------------------------------------------------------
from optparse import OptionParser
parser = OptionParser()
# parser.add_option('--no-indent',action='store_false', dest='indent_flag', default=True,help='do not indent long variables')
# parser.add_option('--params', dest='params', default=None,help='Comma separated list of params')
# parser.add_option("--all",action="store_true", dest="use_all_logs", default=False,help="check all available logs")
# parser.add_option("--nstep",action="store_true", dest="nstep_only", default=False,help="only print nstep values")
# parser.add_option("--partial",action="store_true", dest="allow_partial_match", default=False,help="allow partial matches of input search strings")
# parser.add_option('-n',dest='num_line',default='10',help='sets number of lines to print')
# parser.add_option('--alt',dest='alt_search_str',default='E3SM',help='Sets search string (i.e. "E3SM") to use when searching case name')
(opts, args) = parser.parse_args()
if len(args) < 1 : print('\nERROR: no file provided!\n'); exit()
#-------------------------------------------------------------------------------

file_names = args

eam_var_list = []
scream_var_list = []

eam_var_list.append('T');      scream_var_list.append('T_mid')
eam_var_list.append('Q');      scream_var_list.append('qv')
eam_var_list.append('CLDLIQ'); scream_var_list.append('qc')
eam_var_list.append('CLDICE'); scream_var_list.append('qi')
eam_var_list.append('NUMLIQ'); scream_var_list.append('nc')
eam_var_list.append('NUMICE'); scream_var_list.append('ni')



for f in file_names:
   print()
   print(f)

   ds = xr.open_dataset(f)

   # for key in ds.variables.keys(): print(key)

   if 'NUMLIQ' in ds.variables.keys():
      tmp_var_list = eam_var_list
   else:
      tmp_var_list = scream_var_list

   for key in tmp_var_list:
      print_stat(ds[key],name=key,compact=True,indent=' '*4,stat='nx',fmt='e')
   

   ds.close()

#-------------------------------------------------------------------------------
