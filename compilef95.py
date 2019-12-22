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
arr_li = ['VFUNC', 'PFUNC', 'STATDIST', 'AGRID', 'LORENZ','LAMBDA',
          'VFUND', 'PFUND', 'Q', 'CONSUM_EQ', 'CAPPATH', 'LABPATH', 'RPATH', 'WPATH']
f_dict = dict()
ps_spec_suffix = ""
for c in arr_li:
    c = c+ps_spec_suffix
    f = os.path.join(f_path, c)
    try:
        vf = np.transpose(np.genfromtxt(f))
    except OSError:
        vf = None
    f_dict[c] = vf

#PS6 - GRAPH
time = np.linspace(0,30,30)
age = np.linspace(0,66,66)
plt.plot(age,f_dict['CONSUM_EQ'][0], 'b')
plt.plot(age,f_dict['CONSUM_EQ'][1], 'r')
plt.show()
# required_graph = ['CAPPATH', 'LABPATH',  'RPATH', 'WPATH']
# fig, axs = plt.subplots(2, 2)
# order = [axs[0,0],axs[0,1],axs[1,0],axs[1,1]]
# print(axs)
# def_x = time
# colors = ['r', 'b', 'g', 'pink', 'c']
# for i, item in enumerate(required_graph):
#     g = f_dict[item]
#     current_ax = order[i]
#     print(current_ax)
#     current_ax.set_title(item)
#     current_ax.plot(def_x, f_dict[item], colors[1])
# plt.show()

# ps5 -graphs
# fig, axs = plt.subplots(1, 2)
# required_graph = ['VFUNC', 'PFUNC']
# def_x = f_dict['AGRID']
# colors = ['r', 'b', 'g', 'pink', 'c']
# for i, item in enumerate(required_graph):
#     g = f_dict[item]
#     current_ax = axs[i]
#     print(current_ax)
#     current_ax.set_title(item)
#     for di, dim in enumerate(g):
#         current_ax.plot(def_x, f_dict[item][di], colors[di])
# plt.show()
# plt.plot( f_dict['AGRID'],f_dict['VFUNC'][0], 'r')
# plt.plot( f_dict['AGRID'],f_dict['VFUNC'][1], 'b')
# plt.show()

#PS4B - GRAPHS
# required_graph = ['VFUNC', 'PFUNC',  'VFUND', 'PFUND']
# fig, axs = plt.subplots(2, 2)
# order = [axs[0,0],axs[0,1],axs[1,0],axs[1,1]]
# print(axs)
# def_x = f_dict['AGRID']
# colors = ['r', 'b', 'g', 'pink', 'c']
# for i, item in enumerate(required_graph):
#     g = f_dict[item]
#     current_ax = order[i]
#     print(current_ax)
#     current_ax.set_title(item)
#     for di, dim in enumerate(g):
#         current_ax.plot(def_x, f_dict[item][di], colors[di])
# plt.show()
# plt.plot( f_dict['AGRID'],f_dict['CONSUM_EQ'][2], 'r')
# plt.plot( f_dict['AGRID'],f_dict['CONSUM_EQ'][3], 'b')
# plt.show()

# ps2 - PLOT GRAPHS
# df1 = [f_dict['AGRID'][int(i-1)] for i in f_dict['PFUNC'][0]]
# df2 = [f_dict['AGRID'][int(i-1)] for i in f_dict['PFUNC'][1]]
# plt.plot( f_dict['AGRID'],(df1-f_dict['AGRID']), 'r')
# plt.plot( f_dict['AGRID'],(df2-f_dict['AGRID']), 'b')
# plt.show()

# PS4 - PLOT GRAPHS
# LORENZ
# print(f_dict['LORENZ'])
# ozspace = np.linspace(0,1,100)
# plt.plot(ozspace, f_dict['LORENZ'], 'r')
# plt.plot(ozspace, ozspace, 'b')
# plt.plot( f_dict['AGRID'],f_dict['PFUNC'][0], 'r')
# plt.plot( f_dict['AGRID'],f_dict['PFUNC'][1], 'b')
# plt.show()