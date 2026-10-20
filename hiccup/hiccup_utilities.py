# This module contains some useful utility routines used in hiccup_data_class.py
import os, sys, re, shutil, shlex, subprocess as sp
import numpy as np
verbose_default = False # local verbosity default
# ------------------------------------------------------------------------------
# class for terminal colors
class tcolor:
    """ simple class for coloring terminal text """
    ENDC, BLACK, RED     = '\033[0m','\033[30m','\033[31m'
    GREEN, YELLOW, BLUE  = '\033[32m','\033[33m','\033[34m'
    MAGENTA, CYAN, WHITE = '\033[35m','\033[36m','\033[37m'
# ------------------------------------------------------------------------------
# Common method for printing and running commands
def run_cmd(cmd,verbose=None,prepend_line=True,use_color=True,shell=False,prefix='  ',suffix=''):
    """
    Method to encapsulate running system commands and checking for failures
    cmd can be a string or an argv list - use a list when arguments
    (e.g. file paths) may contain whitespace
    """
    if prepend_line : prefix = '\n'+prefix
    if verbose is None : verbose = verbose_default
    is_argv = isinstance(cmd,(list,tuple))
    cmd_str = shlex.join(cmd) if is_argv else cmd
    msg = f'{prefix}{cmd_str}{suffix}'
    if use_color : msg = tcolor.GREEN + msg + tcolor.ENDC
    if verbose : print(msg)
    if shell:
        sp.check_call(cmd_str,shell=True)
    else:
        sp.check_call(list(cmd) if is_argv else cmd.split())
    return
# ------------------------------------------------------------------------------
# Method for checking if required software is installed
def check_dependency(cmd):
    """ 
    Check for required system commands 
    """
    if shutil.which(cmd) is None : raise OSError(f'{cmd} is not in system path')
    return
# ------------------------------------------------------------------------------
# Methods for parsing and comparing version strings for required software
def suffix_as_tuple(suffix):
    """
    """
    order = ['alpha', 'beta', '']
    suffix_text = "".join(c for c in suffix if not c.isdigit())
    suffix_num = "".join(c for c in suffix if c.isdigit())
    assert(suffix_text in order)
    return (order.index(suffix_text), int(suffix_num) if suffix_num else -1)
def parse_version(version=None):
    """
    parse a version string into a tuple of values, 
    plus a suffix described alpha or beta modifiers
    """
    version_list = version.split('-')
    main_version = tuple(int(n) for n in version_list[0].split("."))
    suffix = version_list[-1] if len(version_list) == 2 else ""
    return main_version, suffix_as_tuple(suffix)
def compare_version(version, required_version=None):
    """
    use tuple version of parsed version string to
    return True if version >= required_version 
    """
    if version == required_version: return True
    v0, suffix0 = parse_version(version)
    v1, suffix1 = parse_version(required_version)
    return (v0 > v1) or ((v0 == v1) and suffix0 >= suffix1)
# ------------------------------------------------------------------------------
# Check version of NCO - and fail if not recent enough
def check_nco_version():
    """
    HICCUP requires features available after a specific NCO version:
    - vertical interpolation bug fix added in 4.9.2-alpha09
    - netcdf string variable handling adding in 5.3.1
    This method parses the version string to check if the version is correct.
    I'm not sure how to handle the "alpha" part of the version string...
    Note - ncks reports the version information through STDERR instead of STDOUT
    """
    msg,err = sp.Popen(['ncks','--version'],stdout=sp.PIPE,stderr=sp.PIPE
                      ,universal_newlines=True).communicate()
    # grab the second line of the version string
    version_str = err.split('\n',1)[1]
    # grab the characters that come after "version" and remove newline character
    version_str = version_str.split('version ',1)[1].replace('\n','')
    min_version = '5.3.1'
    if not compare_version(version_str, required_version=min_version): 
        # current version is not valid, so exit
        err_msg = f'NCO version {version_str} is too old.'
        err_msg += f'\nHICCUP requires NCO version {min_version} or higher'
        raise EnvironmentError(err_msg)
    return
# ------------------------------------------------------------------------------
# Simple routine to check statistics of an array - useful for debugging
def print_stat(x,name='(no name)',unit='',fmt='f',stat='naxh',indent=' '*2,compact=True):
   """
   By default print the min, avg, max, and std dev of input data
   """
   if fmt=='f' : fmt = '%.4f'
   if fmt=='e' : fmt = '%e'
   if unit!='' : unit = f'[{unit}]'
   name_len = 20 if compact else len(name)
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
      if not compact: msg += line+'\n'
   if compact: msg += line
   print(msg)
   return msg
# ------------------------------------------------------------------------------
# Method to check for any number of invalid values
def chk_finite(x,name=None):
  """
  check the input data for inf values
  input data should be a xarray DataArray 
  """
  inf_cnt = np.isinf(x).sum().values
  nan_cnt = np.isnan(x).sum().values
  if inf_cnt>0 or nan_cnt>0: 
    err_msg = '  '
    if name is not None: err_msg += f'{name}:  '
    err_msg += f'invalid values found! {inf_cnt} infs  /  {nan_cnt} nans  '
    raise ValueError(err_msg)
  return
# ------------------------------------------------------------------------------
def inherit_group_from_dir(file_path):
    """Set the group of a file to match the group of its parent directory."""
    parent_dir = os.path.dirname(os.path.abspath(file_path))
    dir_gid = os.stat(parent_dir).st_gid
    os.chown(file_path, -1, dir_gid)
# ------------------------------------------------------------------------------
# Methods for detecting and fixing NaN _FillValue metadata
# (NaN is not a valid _FillValue for ncremap, so it needs to be replaced
# before remapping)
def get_nan_FillValue_vars(file_path):
    """
    Return dict of {variable name : ncatted type code} for all variables
    in file_path that have NaN as the _FillValue (read-only)
    """
    check_dependency('ncdump')
    result = sp.run(['ncdump','-h',file_path], capture_output=True, check=True, text=True)
    # ncdump prints a float NaN as "NaNf" and a double NaN as "NaN"
    nan_vars = {}
    for m in re.finditer(r'^\s*(\S+):_FillValue = NaN(f?) ;', result.stdout, re.MULTILINE):
        nan_vars[m.group(1)] = 'f' if m.group(2)=='f' else 'd'
    return nan_vars
def file_has_nan_FillValue(file_path):
    """Return True if any variable in file_path has NaN as the _FillValue."""
    return len(get_nan_FillValue_vars(file_path)) > 0
def fix_nan_FillValue(file_path, fill_value=1.0e36, verbose=None):
    """
    Replace any NaN _FillValue in file_path with fill_value, modifying the
    file in place. Returns True if any variables were modified.
    """
    if verbose is None : verbose = verbose_default
    check_dependency('ncatted')
    nan_vars = get_nan_FillValue_vars(file_path)
    if len(nan_vars)==0: return False
    # only modify variables with a NaN _FillValue, keeping the variable's type
    # build an argv list so that paths containing whitespace are preserved
    att_opts = []
    for var,typ in nan_vars.items(): att_opts += ['-a',f'_FillValue,{var},m,{typ},{fill_value}']
    run_cmd(['ncatted','-O',*att_opts,file_path], verbose, prepend_line=False)
    return True
# ------------------------------------------------------------------------------
