import subprocess as sp
import shlex
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#
# #Usage: Put this script in the same directory as the Fortran (*.f95) files.
# # Then, Just change the ps_name to the problem set intended and run the script!
# # If you like to change how graphs is plot, modify the code in the end.
#
# def run_fortran(instr):
#     proc = sp.Popen(shlex.split(instr), stdout=sp.PIPE, stderr=sp.PIPE)
#     while proc.poll() is None:  # print stdout realtime/ for debugging
#         l = proc.stdout.readline()  # This blocks until it receives a newline.
#         print(l)
#     proc.stdout.read()
#     proc.wait()
#     stdout = proc.stdout.read().decode()
#     stderr = proc.stderr.read().decode()
#     if proc.returncode:
#         raise ValueError(stdout + stderr)
#
#
# # compile the code on command line
# ps_name = 'ps2' # change for whatever file
# ffile_name = "{}.f95".format(ps_name)
f_path = '/Users/chek_choi/Downloads/fortran/'   # It is also bin path
# fullname = os.path.join(f_path, ffile_name)
# program_name = os.path.join(f_path, ps_name)
# in_str = "gfortran -o {} {}".format(program_name, fullname)
# run_fortran(in_str)
# # ================================================================
# # run the compiled program
# in_str = '{}'.format(program_name)
# run_fortran(in_str)
# get the arrays stored in specific text file
#TODO: MAKE IT MORE GENERIC, SAVE IT TO A SPECIFIC FOLDER AND LOOP THROUGH ALL TXT FILE IN FOLDER
arr_li = ['VFUNC', 'PFUNC', 'STATDIST', 'AGRID', 'LORENZ']
f_dict = dict()
for c in arr_li:
    f = os.path.join(f_path, c)
    try:
        vf = np.transpose(np.genfromtxt(f))
    except OSError:
        vf = None
    f_dict[c] = vf
# # use mpl to plot graphs
# # LORENZ
# ozspace = np.linspace(0,1,700)
# plt.plot(f_list[4], ozspace, 'r')
# plt.plot(ozspace, ozspace, 'b')
# # Policy function
#
# # STAT DIST
# plt.hist(f_list[2][0], color='b', range=(-2, 5))
# plt.hist(f_list[2][1], color='r', range=(-2, 5))
# plt.show()

# ps2
# df1 = [f_dict['AGRID'][int(i-1)] for i in f_dict['PFUNC'][0]]
# df2 = [f_dict['AGRID'][int(i-1)] for i in f_dict['PFUNC'][1]]
# plt.plot( f_dict['AGRID'],(df1-f_dict['AGRID']), 'r')
# plt.plot( f_dict['AGRID'],(df2-f_dict['AGRID']), 'b')
# plt.show()