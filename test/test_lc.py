import matplotlib
matplotlib.use("Agg")
import proplot as pplt
import numpy as np
import sys
sys.path.append("/global/u1/d/dforero/codes/powspec_py/powspec/")
from pypowspec import compute_auto_lc, compute_cross_lc


import pandas as pd
dat_fname = "/global/project/projectdirs/desi/mocks/UNIT/HOD_Shadab/multiple_snapshot_lightcone/UNIT_lightcone_multibox_ELG_footprint_nz_NGC.dat"
ran_fname = "/global/project/projectdirs/desi/mocks/UNIT/HOD_Shadab/multiple_snapshot_lightcone/UNIT_lightcone_multibox_ELG_footprint_nz_1xdata_5.ran_NGC.dat"
data = pd.read_csv(dat_fname, usecols = (0,1,3,4), engine='c', delim_whitespace=True, names = ['ra', 'dec', 'zrsd', 'nz']).values.astype(np.float32)
rand = pd.read_csv(ran_fname, usecols = (0,1,3,4), engine='c', delim_whitespace=True, names = ['ra', 'dec', 'zrsd', 'nz']).values.astype(np.float32)

data = data[(data[:,2] > 0.7) & (data[:,2] < 1.)]
rand = rand[(rand[:,2] > 0.7) & (rand[:,2] < 1.)]

fkp_data = 1. / (1 + 5000 * data[:,3]).astype(np.float32)
fkp_rand = 1. / (1 + 5000 * rand[:,3]).astype(np.float32)
nobj = data.shape[0]
wdata = np.ones(nobj).astype(np.float32)

wrand = np.ones(rand.shape[0]).astype(np.float32)



pk = compute_auto_lc(data[:,0], data[:,1], data[:,2], wdata, fkp_data, data[:,3],
                    rand[:,0], rand[:,1], rand[:,2], wrand, fkp_rand, rand[:,3],
                    powspec_conf_file = "test/powspec_lc.conf",
                    output_file = "test/lc_auto_test.powspec")

print("Plotting", flush=True)
fig, ax = pplt.subplots(nrows=1, ncols=3, share = 0)
for i in range(3):
    ax[i].plot(pk['k'], pk['k'] * pk['multipoles'][:,i], label='auto')
    


pk = compute_cross_lc(data[:,0], data[:,1], data[:,2], wdata, fkp_data, data[:,3],
                    rand[:,0], rand[:,1], rand[:,2], wrand, fkp_rand, rand[:,3],
                    data[:,0], data[:,1], data[:,2], wdata, fkp_data, data[:,3],
                    rand[:,0], rand[:,1], rand[:,2], wrand, fkp_rand, rand[:,3],
                    powspec_conf_file = "test/powspec_lc_cross.conf",
                    output_auto = ["test/lc_auto_test_1.powspec","test/lc_auto_test_1.powspec"],
                    output_cross = "test/lc_cross_test.powspec")
for i in range(3):
    for j in range(2):
        ax[i].plot(pk['k'], pk['k'] * pk['auto_multipoles'][j,:,i], label=j+1)
    ax[i].plot(pk['k'], pk['k'] * pk['cross_multipoles'][:,i], label="1x2")

try:
    res = np.loadtxt("test/lc_auto_ref.powspec")
    for i in range(3):
        ax[i].plot(res[:,0], res[:,0]*res[:,5+i], label="auto ref", ls='--')
    res = np.loadtxt("test/lc_auto_ref_1.powspec")
    for i in range(3):
        ax[i].plot(res[:,0], res[:,0]*res[:,5+i], label="1 ref", ls='--')
    res = np.loadtxt("test/lc_auto_ref_2.powspec")
    for i in range(3):
        ax[i].plot(res[:,0], res[:,0]*res[:,5+i], label="2 ref", ls='--')
    res = np.loadtxt("test/lc_cross_ref.powspec")
    for i in range(3):
        ax[i].plot(res[:,0], res[:,0]*res[:,5+i], label="1x2 ref", ls='--')
    
except OSError as e:
    print(e)
    pass
    


[a.legend(loc='top') for a in ax]
fig.savefig('test/test_lc.png', dpi=300)
