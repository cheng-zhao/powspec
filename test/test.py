import proplot as pplt
import numpy as np
import sys
print(sys.path)
sys.path.append("/global/u1/d/dforero/codes/powspec_py/powspec/")
from pypowspec import *


pk = test()

fig, ax = pplt.subplots(nrows=1, ncols=3, share = 0)
for i in range(1,4):
    ax[i-1].semilogx(pk[:,0], pk[:,0] * pk[:,i])
fig.savefig('test/test.png', dpi=300)