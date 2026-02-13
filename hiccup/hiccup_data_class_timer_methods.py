from hiccup.hiccup_data_class_common import __all__
from hiccup.hiccup_utilities import tcolor
timer_start_total = None
timer_msg_all = []
# ------------------------------------------------------------------------------
def print_timer(timer_start,indent=None,use_color=True,caller=None,print_msg=True):
    """
    Print the final timer result based on input start time
    """
    if indent is None: indent = '  '
    # if caller is not provider get name of parent routine
    if caller is None: caller = sys._getframe(1).f_code.co_name
    # calculate elapsed time
    etime = perf_counter()-timer_start
    time_str = f'{etime:10.1f} sec'
    # add minutes if longer than 60 sec or 2 hours
    if etime>60       : time_str += f' ({(etime/60):4.1f} min)'
    # if etime>(2*3600) : time_str += f' ({(etime/3600):.1f} hr)'
    # create the timer result message
    msg = f'{caller:40} elapsed time: {time_str}'
    # Apply color
    if use_color : msg = tcolor.YELLOW + msg + tcolor.ENDC
    # print the message
    if print_msg: print(f'\n{msg}')
    return msg
# ------------------------------------------------------------------------------
def print_timer_summary(timer_start_total=None,timer_msg_all=None):
    """
    Print timer summary based on information compiled by print_timer()
    """
    # Add timer info for all if timer_start_total was set
    if timer_start_total is not None:
        print_timer(timer_start_total,caller=f'Total',print_msg=False)
    if timer_msg_all is not None:
        print(f'\nHICCUP timer results:')
        for msg in timer_msg_all:
            print(f'  {msg}')
    return
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
