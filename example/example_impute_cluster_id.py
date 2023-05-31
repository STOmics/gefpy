import h5py
import numpy as np

cluster = {}
with open("../test_data/shanghaiLab_mouseBrain_210/FP200000442TL_A2.5.cgef_leiden.txt") as f:
    for line in f:
        if line.startswith(','):
            continue
        fields = line.strip().split(',')
        cluster[int(fields[1])] = int(fields[2]) + 1

h5f = h5py.File("../test_data/shanghaiLab_mouseBrain_210/FP200000442TL_A2.6.cluster.cgef", 'r+')
cell_names = np.bitwise_or(np.left_shift(h5f['cellBin']['cell']['x'].astype('uint64'), 32), h5f['cellBin']['cell']['y'])

celltid = np.zeros(h5f['cellBin']['cell'].shape, dtype='uint16')
n = 0
for cell_name in cell_names:
    if cell_name in cluster:
        celltid[n] = cluster[cell_name]
    n += 1

# h5f['cellBin']['cell']['cellTypeID'] = celltid
h5f['cellBin']['cell']['clusterID'] = celltid