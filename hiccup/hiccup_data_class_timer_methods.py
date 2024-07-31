import sys
from time import perf_counter
from hiccup.hiccup_utilities import tcolor
from hiccup.hiccup_utilities import verbose_indent
timer_msg_all = []
timer_start_total = None
# ------------------------------------------------------------------------------
# Print individual timer information
# ------------------------------------------------------------------------------
def print_timer(self,timer_start,use_color=True,prefix='\n',caller=None,print_msg=True):
    """
    Print the final timer result based on input start time
    Also update timer_msg_all for use in print_timer_summary
    """
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
    # add message to list of messages for print_timer_summary
    timer_msg_all.append(msg)
    # Apply color
    if use_color : msg = tcolor.YELLOW + msg + tcolor.ENDC
    # print the message
    print(prefix+msg)
    return
# ------------------------------------------------------------------------------
# Print a summary of timer information
# ------------------------------------------------------------------------------
def print_timer_summary(self,):
    """
    Print timer summary based on information compiled by print_timer()
    """
    # Add timer info for all if timer_start_total was set
    if self.timer_start_total is not None: 
        self.print_timer(self.timer_start_total,caller=f'Total',print_msg=False)
    if self.do_timers:
        print(verbose_indent+'\nHICCUP Timer results:')
        for msg in self.timer_msg_all:
            print(verbose_indent+f'  {msg}')
    return
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
