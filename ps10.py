import numpy as np
import scipy.integrate as intgr
import scipy.optimize as opt
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from ps10_func import Ps10


# np.random.seed(4200)
# # SetUp
# smm = Ps10() # no change of para
# ###############################################################################
# # Question 2
# # plt.plot(range(smm.t), smm.true_data)
# # plt.show()
# ###############################################################################
# # Question 3
# smm.td()
# smm.moments()
# # np.random.set_state(smm.randstate)
# # print( np.array_equal(np.random.get_state(), smm.randstate))
# # rand_var = smm.sigmaz*np.random.randn(smm.h,smm.t) + smm.mean # in q3, use true vals
# ###############################################################################
# # Setup - mainly for plotting
# rho_range = np.linspace(0.35, 0.65, 100)
# sig_range = np.linspace(0.8, 1.2, 100)
# ###############################################################################
# # Question 4a
# w = np.identity(smm.l_len)
# # smm.plot_j(w)
# init_guess = (0.5,1)
# xfunc = lambda b: smm.j(b, w)
# b_hat = opt.fmin(xfunc, init_guess)
# print('bhat',b_hat)
# ###############################################################################
# # Question 4b
# init_guess = b_hat
# w2 = smm.cal_shat(b_hat)
# xfunc = lambda b: smm.j(b, w2)
# b_hat2 = opt.fmin(xfunc, init_guess)
# print('bhat2',b_hat2)
# ###############################################################################
# # Question 4c
# small_s = 0.001
# delta_g = []
# for i in range(2):
#     b_changed = b_hat
#     b_changed[i] = b_hat[i] - small_s
#     smm.y(b_hat)
#     mom_hat = smm.my
#     smm.y(b_changed)
#     mom_changed = smm.my
#     delta = (mom_hat-mom_changed)/small_s
#     delta_g.append(delta)
# delta_g = np.array(delta_g)
# print("delta_g", delta_g)
# # standard error
# # se = np.matmul(np.transpose(delta_g),w2)
# # print(se.shape, delta_g.shape)
# # se = 1/np.matmul(se, delta_g)/smm.t
# # se = np.sqrt(se)
# # print("standard error",se)
# ###############################################################################
# # Question 4d
# j_test = smm.t*smm.h/(1+smm.h)*smm.j(b_hat2, w2)
# print(j_test)  # test if far from zero
# # ###############################################################################
# # # Question 5
# # #smm = Ps10(l=(1,2))  # use the second and third moment instead
# # ###############################################################################
# # Question 6
rho1_boot = []
rho2_boot = []
for i in range(100):
    np.random.seed(i)
    smm = Ps10(l=(0, 1, 2))
    smm.td()
    smm.moments()
    b_hat = smm.b_hat()
    rho1_boot.append(b_hat[0])
    b_hat2 = smm.b_hat_two(b_hat)
    rho2_boot.append(b_hat2[0])
# Histogram
bs = np.linspace(0.35,0.75, 20)
plt.hist(rho2_boot, bs, alpha=0.5, color='r')
plt.hist(rho1_boot, bs,alpha=0.5, color='b')
plt.show()
