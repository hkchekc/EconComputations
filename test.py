import math
import numpy as np
import itertools

com = itertools.combinations([1,2,3,4,5], 2)
a = np.array([1,2,3,4,5,6,7])
for item in com:
    li = list(item)
    sli = slice(li)
    nli = np.s_[li]
    print(sli)
    print(a[nli])
