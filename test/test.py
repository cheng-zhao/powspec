import proplot as pplt
import numpy as np
import sys
print(sys.path)
sys.path.append("/global/u1/d/dforero/codes/powspec_py/powspec/")
from pypowspec import compute_auto_powspec

seed = 42
np.random.seed(seed)
nobj = int(1e4)
data = 1000. * np.random.random((nobj, 3)).astype(np.double)
w = np.ones(nobj)
pk = compute_auto_powspec(data[:,0], data[:,1], data[:,2], w)
print("Plotting", flush=True)
fig, ax = pplt.subplots(nrows=1, ncols=3, share = 0)
for i in range(1,4):
    ax[i-1].semilogx(pk[:,0], pk[:,0] * pk[:,i])
fig.savefig('test/test.png', dpi=300)