import matplotlib
matplotlib.use("Agg")
import proplot as pplt
import numpy as np
import sys
sys.path.append("/global/u1/d/dforero/codes/powspec_py/powspec/")
from pypowspec import compute_auto_box_rand, compute_cross_box_rand, compute_auto_box


def hankel_sum(k, P, smooth_a, sarr, ell):
    P = P.copy()
    #P *= k ** 3 * np.exp(-(k ** 2) * smooth_a ** 2) * 0.5 / np.pi ** 2
    #delta = np.log(k[1]) - np.log(k[0])
    P *= k ** 2 * np.exp(-(k ** 2) * smooth_a ** 2) * 0.5 / np.pi ** 2
    delta = k[1] - k[0]

    ks = sarr[:, None] * k[None, :]
    j0 = np.sin(ks) / ks
    # j0 = np.sinc(ks)
    j2 = -(3.0 / ks ** 2 - 1.0) * j0 + 3.0 * np.cos(ks) / ks ** 2
    j4 = (5 * (2 * ks ** 2 - 21.0) * np.cos(ks)) / ks ** 4 + (
        (ks ** 4 - 45 * ks ** 2 + 105.0) * np.sin(ks)
    ) / ks ** 5

    if ell == 0:
        xi = (P[None, :] * j0 * delta).sum(axis=1)
    elif ell == 2:
     xi = (P[None, :] * j2 * delta).sum(axis=1)
    elif ell == 4:
        xi = (P[None, :] * j4 * delta).sum(axis=1)
    else:
        raise NotImplementedError

    return sarr, xi#, xi2, xi4


try:
    import pandas as pd
    data_fname = "/global/cfs/cdirs/desi/users/UNIT-BAO-RSD-challenge/Reconstruction/Results_stage2/PyRec_FD/UNITSIM/LRG/LRG-wpmax-v3-snap103-redshift0.74_dens0_s10.0.dat.npy"
    data = np.load(data_fname).astype(np.double)
    rand_fname = "/global/cfs/cdirs/desi/users/UNIT-BAO-RSD-challenge/Reconstruction/Results_stage2/PyRec_FD/UNITSIM/LRG/LRG-wpmax-v3-snap103-redshift0.74_dens0_s10.0_seed7.ran.npz"
    rand = np.load(rand_fname)['sym'].astype(np.double)
    nobj = data.shape[0]
    ref_fn = "/global/cfs/cdirs/desi/users/UNIT-BAO-RSD-challenge/Reconstruction/Results_stage2/PyRec_FD/UNITSIM/LRG/LRG-wpmax-v3-snap103-redshift0.74_dens0_s10.0_seed7.sym_powspec.dat"
except:
    print("WARNING: Read catalog failed, testing with uniform random.", flush=True)
    seed = 42
    np.random.seed(seed)
    nobj = int(1e4)
    data = 1000. * np.random.random((nobj, 3)).astype(np.double)
    np.random.seed(0)
    rand = (data.max() - data.min()) * np.random.random((50 * nobj, 3)).astype(np.double)

kmin = 5e-3
kmax = 1
wdata = np.ones(data.shape[0])
wrand = np.ones(rand.shape[0])

pk = compute_cross_box_rand(data[:,0], data[:,1], data[:,2], wdata, 
                       rand[:,0], rand[:,1], rand[:,2], wrand, 
                       data[:,0], data[:,1], data[:,2], wdata, 
                       rand[:,0], rand[:,1], rand[:,2], wrand,
                       powspec_conf_file = "test/powspec_cross.conf",
                       output_auto = ["test/box_auto_test_recon_1.powspec", "test/box_auto_test_recon_2.powspec"],
                       output_cross = "test/box_cross_test_recon.powspec")
kmask = (pk['k'] > kmin) & (pk['k'] < kmax)
print("Plotting", flush=True)
fig, ax = pplt.subplots(nrows=2, ncols=3, share = 0)

for i in range(3):
    for j in range(2):
        ax[0,i].semilogx(pk['k'], pk['k'] * pk['auto_multipoles'][j,:,i], label=j+1)
        s, xi = hankel_sum(pk['k'][kmask],  pk['auto_multipoles'][j,:,i][kmask], smooth_a = 3, sarr= np.arange(1e-5, 200, 1), ell = 2 * i)
        ax[1,i].plot(s, s**2*xi)
    ax[0,i].semilogx(pk['k'], pk['k'] * pk['cross_multipoles'][:,i], label="1x2")
    s, xi = hankel_sum(pk['k'][kmask],  pk['cross_multipoles'][:,i][kmask], smooth_a = 3, sarr= np.arange(1e-5, 200, 1), ell = 2 * i)
    ax[1,i].plot(s, s**2*xi)
fig.savefig('test/test_box_recon.png', dpi=300)

pk = compute_auto_box_rand(data[:,0], data[:,1], data[:,2], wdata,
                      rand[:,0], rand[:,1], rand[:,2], wrand, 
                      powspec_conf_file = "test/powspec_auto.conf",
                      output_file = "test/box_auto_test_recon.powspec")
print("Plotting", flush=True)

for i in range(3):
    ax[0,i].semilogx(pk['k'], pk['k'] * pk['multipoles'][:,i], label='auto')
    s, xi = hankel_sum(pk['k'][kmask],  pk['multipoles'][:,i][kmask], smooth_a = 3, sarr= np.arange(1e-5, 200, 1), ell = 2 * i)
    ax[1,i].plot(s, s**2*xi)
fig.savefig('test/test_box_recon.png', dpi=300)


try:
    res = np.loadtxt(ref_fn)
    for i in range(3):
        ax[0,i].semilogx(res[:,0], res[:,0]*res[:,1+i], label="auto ref", ls='--')
        s, xi = hankel_sum(res[:,0], res[:,1+i], smooth_a = 3, sarr= np.arange(1e-5, 200, 1), ell = 2 * i)
        ax[1,i].plot(s, s**2*xi, ls='--', label='auto ref')
    
except OSError as e:
    print(e)
    pass

[a.legend(loc='top') for a in ax]
fig.savefig('test/test_box_recon.png', dpi=300)
