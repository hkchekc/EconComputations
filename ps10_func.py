import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


class Ps10:

    def __init__(self,t=200, sigmaz=2,rhoz=0.5,xz=0, h=10, l=(0,1,2)):
        self.t = t
        self.sigmaz = sigmaz
        self.rhoz = rhoz
        self.xz = xz
        self.randstate = np.random.get_state()
        self.h =h
        self.l_len=len(l)
        self.l = np.s_[0:3]
        self.mx = None
        self.my = None
        self.true_data = None
        self.mean = None

    def td(self):
        x_old = 0.
        x = []
        for i in range(self.t):
            epsilon = np.random.randn(1)
            x_new = self.rhoz*x_old + epsilon
            x.append(x_new)
            x_old = x_new
        self.true_data = x

    def moments(self, da=None): # based on the calculation of question 1
        td_bool =False
        if da is None:
            da = self.true_data
            td_bool = True
        mean = np.mean(da)
        vari = np.var(da)
        if td_bool:
            first_ar_co = np.mean([(da[i]-mean)*
                                    (da[i-1]-mean) for i in range(1, self.t)])
        else:
            first_ar_co = 0
            for hi in range(self.h):
                tmp_da = da[hi]
                ar_co_h = np.mean([(tmp_da[i]-mean)*
                                    (tmp_da[i-1]-mean) for i in range(1, self.t)])
                first_ar_co += ar_co_h
            first_ar_co /= self.h
        m_vec = [mean, vari, first_ar_co]
        if td_bool:
            self.mean = m_vec[0]
            self.mx = np.array(m_vec)[self.l]
        else:
            self.my = np.array(m_vec)[self.l]

    def y(self, b):
        y_arr=[]
        rho = b[0]
        sigma = b[1]
        rvar = sigma*np.random.randn(self.h, self.t) + self.mean
        for hi in range(self.h):
            rv_h = rvar[hi]
            y_h = [float(rv_h[0])]
            for ti in range(1, self.t):
                tmp = rv_h[ti]+rho*rv_h[ti-1]
                y_h.append(tmp)
            y_arr.append(y_h)
        self.moments(da=y_arr)


    def j(self,b, w):
        self.y(b)
        m_diff = self.mx - self.my
        j_val = np.matmul(m_diff, w)
        j_val = np.matmul(j_val, m_diff.reshape(self.l_len, 1))
        return j_val

    def cal_shat(self, b):
        self.y(b)
        tmp_diff = self.true_data - self.mean
        tmp_ac = np.append(0, tmp_diff[1:]*tmp_diff[:self.t-1])
        tmp1_arr = np.vstack((tmp_diff.reshape(self.t), tmp_diff.reshape(self.t)**2, tmp_ac))[self.l]
        gamma0 = np.matmul(tmp1_arr, np.transpose(tmp1_arr))
        gamma_arr = []
        for i in range(4):
            gamma = np.zeros((self.l_len, self.l_len))
            for tidx in range(i+1, self.t):
                tmp1 = tmp1_arr[:,tidx]
                tmp2 = tmp1_arr[:,tidx-i]
                gamma += np.matmul(tmp1, np.transpose(tmp2))
            gamma_arr.append(gamma)
        gamma_arr = np.array(gamma_arr)
        shat = np.zeros((self.l_len, self.l_len)) + gamma0
        for j in range(4):
            shat += (1 - j - 1 / (4 + 1)) * (gamma_arr[j] + np.transpose(gamma_arr[j]))
        w2 = np.linalg.inv(shat)
        return w2

    def ref_moments(self): # based on the calculation of question 1
        mean = 0
        vari = self.sigmaz**2/(1-self.rhoz**2)
        first_ar_co = self.rhoz*vari
        m_vec = [mean, vari, first_ar_co]
        return  m_vec

    @staticmethod
    def pop_b_arr(rho_range, sig_range):
        b_arr = []
        for ri, rh in enumerate(rho_range):
            for si, sig in enumerate(sig_range):
                b_arr.append((rh, sig))
        return b_arr

    def plot_j(self, w):
        rho_range = np.linspace(0.35, 0.65, 100)
        sig_range = np.linspace(0.8, 1.2, 100)
        b_arr = self.pop_b_arr(rho_range, sig_range)
        j_arr = []
        for b in b_arr:
            j_arr.append(self.j(b, w))
        j_arr = np.array(j_arr).reshape(100, 100)
        # plot 3d
        print(j_arr.shape)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.contour3D(rho_range, sig_range, j_arr, 100, cmap='viridis')
        ax.set_xlabel('rho')
        ax.set_ylabel('sigma')
        ax.set_zlabel('J value')
        plt.show()

    def b_hat(self):
        w = np.identity(self.l_len)
        init_guess = (0.5, 1)
        xfunc = lambda b: self.j(b, w)
        b_hat = opt.fmin(xfunc, init_guess)
        return b_hat

    def b_hat_two(self, b_hat):
        init_guess = b_hat
        w2 = self.cal_shat(b_hat)
        xfunc = lambda b: self.j(b, w2)
        b_hat2 = opt.fmin(xfunc, init_guess)
        return b_hat2


