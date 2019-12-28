import numpy as np
import scipy.integrate as intgr
import scipy.optimize as opt
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from ps10_func import Ps10


np.random.seed(4200)
# SetUp
smm = Ps10() # no change of para
###############################################################################
# Question 2
# plt.plot(range(smm.t), smm.true_data)
# plt.show()
###############################################################################
# Question 3
smm.moments()
# np.random.set_state(smm.randstate)
# print( np.array_equal(np.random.get_state(), smm.randstate))
# rand_var = smm.sigmaz*np.random.randn(smm.h,smm.t) + smm.mean # in q3, use true vals
###############################################################################
# Setup - mainly for plotting
rho_range = np.linspace(0.35, 0.65, 100)
sig_range = np.linspace(0.8, 1.2, 100)
###############################################################################
# Question 4a
w = np.identity(2)
smm.plot_j(w)
init_guess = (0.5,1)
xfunc = lambda b: smm.j(b, w)
b_hat = opt.fmin(xfunc, init_guess)
print(b_hat)
###############################################################################
# Question 4b
ref_mu = smm.ref_moments()
for i in range(4):
    s_vec = smm.mx - ref_mu[smm.l]
    s = np.matmul(s_vec, np.transpose(s_vec))
    print()
w2 = np.invert(s)
init_guess = (0.5,1)
xfunc = lambda b: smm.j(b,w2)
b_hat2 = opt.fmin(xfunc, init_guess)
print(b_hat2)
###############################################################################
# Question 4c
small_s = 0.001
delta_g = []
for i in range(2):
    b_changed = b_hat
    b_changed[i] = b_hat[i] - small_s
    smm.y(b_hat)
    mom_hat = smm.my
    smm.y(b_changed)
    mom_changed = smm.my
    delta = (mom_hat-mom_changed)/small_s
    delta_g.append(delta)
print(delta_g)
# standard error
se = np.matmul(np.transpose(delta_g),w2)
se = 1/np.matmul(se, delta_g)/smm.t
se = np.sqrt(se)
print(se)
###############################################################################
# Question 4d
j_test = smm.t*smm.h/(1+smm.h)*smm.j(b_hat2, w2)
print(j_test)  # test if far from zero
###############################################################################
# Question 5
#smm = Ps10(l=(1,2))  # use the second and third moment instead
###############################################################################
# Question 6
rho1_boot = []
rho2_boot = []
for i in range(200):
    smm = Ps10(l=(0, 1, 2))
# Histogram
probs = np.linspace(0.1,1, 10)
plt.hist(probs, rho1_boot, 'b')
plt.hist(probs, rho2_boot, 'r')
plt.show()
