import proplot as pplt
import numpy as np
import sys
sys.path.append("/global/u1/d/dforero/codes/powspec_py/powspec/")
from pypowspec import compute_auto_powspec

try:
    import pandas as pd
    test_fname = "/global/project/projectdirs/desi/mocks/UNIT/HOD_Shadab/HOD_boxes/redshift0.9873/UNIT_DESI_Shadab_HOD_snap97_ELG_v0.txt"
    data = pd.read_csv(test_fname, usecols = (0,1,3), engine='c', delim_whitespace=True, names = ['x', 'y', 'zrsd']).values
except:
    print("WARNING: Read catalog failed, testing with uniform random.", flush=True)
    seed = 42
    np.random.seed(seed)
    nobj = int(1e4)
    data = 1000. * np.random.random((nobj, 3)).astype(np.double)


nobj = data.shape[0]
w = np.ones(nobj)
pk = compute_auto_powspec(data[:,0], data[:,1], data[:,2], w, output_file = "test/test.powspec")
print(pk)
print("Plotting", flush=True)
fig, ax = pplt.subplots(nrows=1, ncols=3, share = 0)
for i in range(3):
    ax[i-1].semilogx(pk['k'], pk['k'] * pk['multipoles'][:,i])
fig.savefig('test/test.png', dpi=300)