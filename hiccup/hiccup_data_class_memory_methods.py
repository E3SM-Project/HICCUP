# ------------------------------------------------------------------------------
# Methods for monitoring memory usage
# ------------------------------------------------------------------------------
from hiccup.hiccup_data_class_common import *
from hiccup.hiccup_utilities import tcolor
# ------------------------------------------------------------------------------
def get_mem_usage():
    """
    get current RSS and VMS memory usage
    """
    mem_rss_GB = psutil.Process().memory_info().rss * 1e-9
    mem_vms_GB = psutil.Process().memory_info().vms * 1e-9
    return (mem_rss_GB, mem_vms_GB)
# ------------------------------------------------------------------------------
def print_mem_usage(indent=None,msg=None,use_color=True,print_msg=True):
    """
    construct an informative print statement to describe current memory usage
    """
    if indent is None: indent = '  '
    # if no text is provided then use name of parent routine calling this method
    if msg is None: msg = sys._getframe(1).f_code.co_name
    (mem_rss_GB, mem_vms_GB) = get_mem_usage()
    print_msg = f'{indent}Memory: {mem_rss_GB:8.2f} / {mem_vms_GB:8.2f}  GB  (RSS/VMS)'
    if use_color: print_msg = tcolor.CYAN + print_msg + tcolor.ENDC
    print_msg += f'  {msg}'
    if print_msg: print(print_msg)
    return (print_msg, mem_rss_GB, mem_vms_GB)
# ------------------------------------------------------------------------------