import subprocess as sp
import shlex
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#Usage: Put this script in the same directory as the Fortran (*.f95) files.
# Then, Just change the ps_name to the problem set intended and run the script!
# If you like to change how graphs is plot, modify the code in the end.

def run_fortran(instr):
    proc = sp.Popen(shlex.split(instr), stdout=sp.PIPE, stderr=sp.PIPE)
    while proc.poll() is None:  # print stdout realtime/ for debugging
        l = proc.stdout.readline()  # This blocks until it receives a newline.
        print(l)
    proc.stdout.read()
    proc.wait()
    stdout = proc.stdout.read().decode()
    stderr = proc.stderr.read().decode()
    if proc.returncode:
        raise ValueError(stdout + stderr)


# compile the code on command line
ps_name = 'ps4'
ffile_name = "{}.f95".format(ps_name)
f_path = '/Users/chek_choi/Downloads/fortran/'   # It is also bin path
fullname = os.path.join(f_path, ffile_name)
program_name = os.path.join(f_path, ps_name)
in_str = "gfortran -o {} {}".format(program_name, fullname)
run_fortran(in_str)
# ================================================================
# run the compiled program
in_str = '{}'.format(program_name)
run_fortran(in_str)
# get the arrays stored in specific text file
#TODO: MAKE IT MORE GENERIC, SAVE IT TO A SPECIFIC FOLDER AND LOOP THROUGH ALL TXT FILE IN FOLDER
f_list = []
for c in ['VALUEFUN', 'POLICYFUN', 'STATDIST', 'AGRID']:
    f = os.path.join(f_path, c)
    vf = np.transpose(np.genfromtxt(f))
    f_list.append(vf)
print(f_list)
# use mpl to plot graphs
plt.plot(f_list[0], f_list[3])