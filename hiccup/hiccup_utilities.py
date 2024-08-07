# This module contains some useful utility routines used in hiccup_data_class.py
import os, sys, re, shutil, subprocess as sp
# from time import perf_counter
# timer_msg_all = []
# ------------------------------------------------------------------------------
# Global verbosity default
hiccup_verbose = False
verbose_indent = ''
# ------------------------------------------------------------------------------
# class for terminal colors
# ------------------------------------------------------------------------------
class tcolor:
    """ simple class for coloring terminal text """
    ENDC, BLACK, RED     = '\033[0m','\033[30m','\033[31m'
    GREEN, YELLOW, BLUE  = '\033[32m','\033[33m','\033[34m'
    MAGENTA, CYAN, WHITE = '\033[35m','\033[36m','\033[37m'
# ------------------------------------------------------------------------------
# Common method for printing and running commands
# ------------------------------------------------------------------------------
def run_cmd(cmd,verbose=None,prepend_line=True,use_color=True,shell=False,prefix='  ',suffix=''):
    """
    Method to encapsulate running system commands and checking for failures
    """
    if prepend_line : prefix = '\n'+prefix
    if verbose is None : verbose = hiccup_verbose
    msg = f'{prefix}{cmd}{suffix}'
    if use_color : msg = tcolor.GREEN + msg + tcolor.ENDC
    if verbose : print(msg)
    if shell:
        sp.check_call(cmd,shell=True)
    else:
        sp.check_call(cmd.split())
    return
# ------------------------------------------------------------------------------
# Method for checking if required software is installed
# ------------------------------------------------------------------------------
def check_dependency(cmd):
    """ 
    Check for required system commands 
    """
    if shutil.which(cmd) is None : raise OSError(f'{cmd} is not in system path')
    return
# ------------------------------------------------------------------------------
# Methods for parsing and comparing version strings for required software
# ------------------------------------------------------------------------------
def suffix_as_tuple(suffix):
    """
    """
    order = ['alpha', 'beta', '']
    suffix_text = "".join(c for c in suffix if not c.isdigit())
    suffix_num = "".join(c for c in suffix if c.isdigit())
    assert(suffix_text in order)
    return (order.index(suffix_text), int(suffix_num) if suffix_num else -1)
def parse_version(version='4.9.2-alpha'):
    """
    parse a version string into a tuple of values, 
    plus a suffix described alpha or beta modifiers
    """
    version_list = version.split('-')
    main_version = tuple(int(n) for n in version_list[0].split("."))
    suffix = version_list[-1] if len(version_list) == 2 else ""
    return main_version, suffix_as_tuple(suffix)
def compare_version(version, required_version='4.9.2-alpha'):
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
# ------------------------------------------------------------------------------
def check_nco_version():
    """
    NCO needs to include a vertical interpolation bug fix added in 4.9.2-alpha09
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
    min_version = '4.9.2-alpha9'
    if not compare_version(version_str, required_version='4.9.2-alpha'): 
        # current version is not valid, so exit
        err_msg = f'NCO version {version_str} is too old.'
        err_msg += f'\nHICCUP requires NCO version {min_version} or higher'
        raise EnvironmentError(err_msg)
    return
# # ------------------------------------------------------------------------------
# # Print individual timer information
# # ------------------------------------------------------------------------------
# def print_timer(timer_start,use_color=True,prefix='\n',caller=None,print_msg=True):
#     """
#     Print the final timer result based on input start time
#     Also update timer_msg_all for use in print_timer_summary
#     """
#     # if caller is not provider get name of parent routine
#     if caller is None: caller = sys._getframe(1).f_code.co_name
#     # calculate elapsed time
#     etime = perf_counter()-timer_start
#     time_str = f'{etime:10.1f} sec'
#     # add minutes if longer than 60 sec or 2 hours
#     if etime>60       : time_str += f' ({(etime/60):4.1f} min)'
#     # if etime>(2*3600) : time_str += f' ({(etime/3600):.1f} hr)'
#     # create the timer result message
#     msg = f'{caller:40} elapsed time: {time_str}'
#     # add message to list of messages for print_timer_summary
#     timer_msg_all.append(msg)
#     # Apply color
#     if use_color : msg = tcolor.YELLOW + msg + tcolor.ENDC
#     # print the message
#     print(prefix+msg)
#     return
# # ------------------------------------------------------------------------------
# # Print a summary of timer information
# # ------------------------------------------------------------------------------
# def print_timer_summary():
#     """
#     Print timer summary based on information compiled by print_timer()
#     """
#     # Add timer info for all if timer_start_total was set
#     if timer_start_total is not None: 
#         print_timer(timer_start_total,caller=f'Total',print_msg=False)
#     if do_timers:
#         print(verbose_indent+'\nHICCUP Timer results:')
#         for msg in timer_msg_all:
#             print(verbose_indent+f'  {msg}')
#     return
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
