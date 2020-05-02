#!/usr/bin/env python 
#
#  file: run_area_cmdline1.py
#
#  This python script runs the C++ code area_cmdline for a sequence
#   of values that are hardcoded in the script.  The output is to
#   the screen.
#
#  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
#
#  Revision history:
#      28-Dec-2010  original version
#
#  Notes:  
#   * run with: "python run_area_cmdline1.py" or just "./run_area_cmdline1.py"
#   * should check the return code to make sure it worked!
#   * could use Popen instead of call
#   * pass a full string to the shell with the first method below
#      or pass a list of arguments without "shell=True".
#  
#*************************************************************************
from subprocess import call    # make the "call" function available

import numpy as np

# Define the inputs at the top, so they are easy to find and change
value_list = np.arange(5,26,5)
again_list = np.arange(1,6,1)
# Ok, let's do it . . .

for i in range(value_list.size):    # don't forget the colon!
  my_command = "./area_cmdline.x " + str(value_list[i])+" "+str(again_list[i])  # convert radius to a string
  retcode = call(my_command, shell=True)    # pass "my_command" to be executed

#*************************************************************************
