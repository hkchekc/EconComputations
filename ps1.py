import numpy as np
import math
import collections
from sklearn.impute import SimpleImputer
import pandas as pd

# fname = "/Users/chek_choi/Downloads/J269830.xlsx" # wealth
# fname = "/Users/chek_choi/Downloads/J269946.xlsx"  # income
fname = "/Users/chek_choi/Downloads/J269947.xlsx"  # earnings

# get downloaded data
data = pd.read_excel(fname)
# data[data<0] = 1 # for income
data = data[:7328] # FOR EARNINGS
data[data> 9999990] = np.nan # for income and earning
print(data['w1989'][0:4])

data = data.dropna()
da_len = len(data['w1989'])
# define quantile
no_bins = 5
quan89 = []
quan94 = []
p = 1/no_bins
for i in range(no_bins):
    q = i*p
    quan89.append(np.quantile(data['w1989'], q))
    quan94.append(np.quantile(data['w1994'], q))

# mark the arrays in [q89] and [q94] matrices
mask89 = np.zeros(len(data['w1989']))
mask94 = np.zeros(len(data['w1994']))
tmp = mask89
for i in range(no_bins):
    tmp1 = np.ma.masked_greater(data['w1989'].to_numpy(), quan89[i]).mask.astype(int)  # mask 1989
    tmp2 = np.ma.masked_greater(data['w1994'].to_numpy(), quan94[i]).mask.astype(int)  # mask 1994
    mask89 = [x+y for x,y in zip(tmp1, mask89)]
    mask94 = [(x+y) for x, y in zip(tmp2, mask94)]
print(mask89)
print(mask94)
print(quan89)
print(quan94)
# add up the mask and put in transition matrix
trans_mat = np.zeros((no_bins, no_bins))
tmp = np.zeros(shape=(len(data['w1994'])))
for i in range(no_bins): # loop for 89
    com89 = np.equal(np.array(mask89), np.ones(da_len)*(i+1))
    for j in range(no_bins):  # loop for 94
        com94 = np.equal(np.array(mask94), np.ones(da_len)*(j+1))
        count = np.where(np.logical_and(com89,com94), 1, 0)
        # print(count94, 'cooou')
        trans_mat[i, j] = sum(count)
trans_mat /= sum(trans_mat)
print(sum(trans_mat), "sum")

trans_mat = pd.DataFrame(trans_mat)
print(trans_mat)

