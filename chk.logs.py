#!/usr/bin/env python
#===================================================================================================
#  Nov, 2017 - Walter Hannah - Lawrence Livermore National Lab
#  This script prints the tail end of the most recently modified HICCUP output log
#===================================================================================================
import sys, os, glob
import subprocess as sp
from optparse import OptionParser

help_header = 'usage: ./%prog [-n] [-l]\n'
help_header += '\nThis script will look for log files in the current directory (./logs*/*) '
help_header += '\nand print the tail end of the most recently modified file'
parser = OptionParser(usage=help_header)
parser.add_option('-l',dest='num_logs',default=1,help='sets number of logs to print')
parser.add_option('-n',dest='num_line',default=10,help='sets number of lines to print')
(opts, args) = parser.parse_args()

# Set up terminal colors
class bcolor:
    ENDC     = '\033[0m';   BLACK    = '\033[30m'
    RED      = '\033[31m';  GREEN    = '\033[32m'
    YELLOW   = '\033[33m';  BLUE     = '\033[34m'
    MAGENTA  = '\033[35m';  CYAN     = '\033[36m'
    WHITE    = '\033[37m';  BOLD     = '\033[1m'

#===================================================================================================
#===================================================================================================

log_path = './logs*/*'

# get list of log files sorted by modification time
logs = sorted( glob.glob( log_path ), key=lambda file: os.path.getmtime(file) )

#===================================================================================================
# Loop through log files
#===================================================================================================
cnt = 0
nmsg = int(opts.num_logs)
nlog = len(logs)
for l in reversed(range(nlog-nmsg,nlog)) :

    log_file = logs[l]

    # print the log file name
    print('\n'+bcolor.BOLD+log_file+bcolor.ENDC)

    # create the command to look at the tail of the file
    cmd = f'tail {log_file} -n {opts.num_line}'

    # indent the output
    cmd += '| sed -e \'s/^/    /\' '

    # run the command
    sp.call(cmd,shell=True)

#===================================================================================================
#===================================================================================================
