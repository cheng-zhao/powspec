import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/users/d/dforeros/codes/py_powspec/powspec")
from pypowspec import compute_auto_box_mesh

data = np.load("/home/users/d/dforeros/projects/desi-patchy/run/collapsed.npy").T#1000. * np.random.random(size = (100000, 3))
def ngp(grid_size, box_size, data):
    i = ((data[:,0] / box_size) * grid_size).astype(int)
    j = ((data[:,1] / box_size) * grid_size).astype(int)
    k = ((data[:,2] / box_size) * grid_size).astype(int)
    
    grid = np.zeros((grid_size,)*3)
    for ii in range(data.shape[0]):
        grid[i[ii], j[ii], k[ii]] += 1
    #grid[(i,j,k)] += 1
    return grid
grid = ngp(512, 2000., data)


print(grid.sum())
pk = compute_auto_box_mesh(grid, "test/powspec_auto.conf")

print(pk['multipoles'].shape)
plt.plot(pk['k'], pk['k'] * pk['multipoles'][:,0])
plt.savefig("test/test_mesh.png")
