from hiccup.hiccup_data_class_common import *
from hiccup.hiccup_utilities import tcolor
timer_start_total = None
timer_msg_all = []
# ------------------------------------------------------------------------------
def print_timer(timer_start,indent=None,use_color=True,caller=None,print_msg=True):
    """
    Print the final timer result based on input start time
    """
    if indent is None: indent = '' # no indent by default
    # if caller is not provider get name of parent routine
    if caller is None: caller = sys._getframe(1).f_code.co_name
    # calculate elapsed time
    etime = perf_counter()-timer_start
    time_str = f'{etime:10.1f} sec'
    # add minutes if longer than 60 sec or 2 hours
    if etime>60       : time_str += f' ({(etime/60):4.1f} min)'
    # if etime>(2*3600) : time_str += f' ({(etime/3600):.1f} hr)'
    # create the timer result message
    msg = f'{indent}{caller:40} elapsed time: {time_str}'
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
    msg_list = list(timer_msg_all) if timer_msg_all is not None else []
    # Add timer info for all if timer_start_total was set
    if timer_start_total is not None:
        msg_list.append( print_timer(timer_start_total,caller=f'Total',print_msg=False) )
    if len(msg_list)==0: return
    # use hours instead of minutes if any line is longer than 120 minutes
    time_pattern = re.compile(r'(\d+\.\d+) sec(?: \(\s*\d+\.\d+ min\))?')
    etime_list = [ float(m.group(1)) for m in map(time_pattern.search,msg_list) if m ]
    use_hours = any( etime>(120*60) for etime in etime_list )
    def convert_time_str(match):
        etime = float(match.group(1))
        time_str = f'{etime:10.1f} sec'.strip()
        if use_hours:
            if etime>60 : time_str += f' ({(etime/3600):4.1f} hr)'
        else:
            if etime>60 : time_str += f' ({(etime/60):4.1f} min)'
        return time_str
    print(f'\nHICCUP timer results:')
    for msg in msg_list:
        print(f'  {time_pattern.sub(convert_time_str,msg)}')
    return
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
