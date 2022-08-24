import matplotlib
matplotlib.use("Agg")
import proplot as pplt
import numpy as np
import sys
sys.path.append("/global/u1/d/dforero/codes/powspec_py/powspec/")
from pypowspec import compute_auto_box_rand, compute_cross_box_rand

try:
    import pandas as pd
    test_fname = "/global/project/projectdirs/desi/mocks/UNIT/HOD_Shadab/HOD_boxes/redshift0.9873/UNIT_DESI_Shadab_HOD_snap97_ELG_v0.txt"
    data = pd.read_csv(test_fname, usecols = (0,1,3), engine='c', delim_whitespace=True, names = ['x', 'y', 'zrsd']).values
    nobj = data.shape[0]
except:
    print("WARNING: Read catalog failed, testing with uniform random.", flush=True)
    seed = 42
    np.random.seed(seed)
    nobj = int(1e4)
    data = 1000. * np.random.random((nobj, 3)).astype(np.double)

np.random.seed(0)
rand = (data.max() - data.min()) * np.random.random((50 * nobj, 3)).astype(np.double)
wdata = np.ones(data.shape[0])
wrand = np.ones(rand.shape[0])

pk = compute_cross_box_rand(data[:,0], data[:,1], data[:,2], wdata, 
                       rand[:,0], rand[:,1], rand[:,2], wrand, 
                       data[:,0], data[:,1], data[:,2], wdata, 
                       rand[:,0], rand[:,1], rand[:,2], wrand,
                       powspec_conf_file = "test/powspec_cross.conf",
                       output_auto = ["test/box_auto_test_rand_1.powspec", "test/box_auto_test_rand_2.powspec"],
                       output_cross = "test/box_cross_test_rand.powspec")

print("Plotting", flush=True)
fig, ax = pplt.subplots(nrows=1, ncols=3, share = 0)

for i in range(3):
    for j in range(2):
        ax[i].semilogx(pk['k'], pk['k'] * pk['auto_multipoles'][j,:,i], label=j+1)
    ax[i].semilogx(pk['k'], pk['k'] * pk['cross_multipoles'][:,i], label="1x2")
fig.savefig('test/test_box_rand.png', dpi=300)

pk = compute_auto_box_rand(data[:,0], data[:,1], data[:,2], wdata,
                      rand[:,0], rand[:,1], rand[:,2], wrand, 
                      powspec_conf_file = "test/powspec_auto.conf",
                      output_file = "test/box_auto_test_rand.powspec")
print("Plotting", flush=True)

for i in range(3):
    ax[i].semilogx(pk['k'], pk['k'] * pk['multipoles'][:,i], label='auto')
fig.savefig('test/test_box_rand.png', dpi=300)


try:
    res = np.loadtxt("test/box_auto_ref.powspec")
    for i in range(3):
        ax[i].semilogx(res[:,0], res[:,0]*res[:,5+i], label="auto ref", ls='--')
    res = np.loadtxt("test/box_auto_ref_1.powspec")
    for i in range(3):
        ax[i].semilogx(res[:,0], res[:,0]*res[:,5+i], label="1 ref", ls='--')
    res = np.loadtxt("test/box_auto_ref_2.powspec")
    for i in range(3):
        ax[i].semilogx(res[:,0], res[:,0]*res[:,5+i], label="2 ref", ls='--')
    res = np.loadtxt("test/box_cross_ref.powspec")
    for i in range(3):
        ax[i].semilogx(res[:,0], res[:,0]*res[:,5+i], label="1x2 ref", ls='--')

    
except OSError as e:
    print(e)
    pass

[a.legend(loc='top') for a in ax]
fig.savefig('test/test_box_rand.png', dpi=300)
