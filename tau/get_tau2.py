import numpy as np
from netCDF4 import Dataset

nlat = 31
nlon = 30
year0 = 1961
yearf = 1994
gridName = "%dx%d" % (nlat, nlon)
periodName = "%d%d" % (year0, yearf)

sstFile = "ersst.%s_%s.nc"  % (gridName, periodName)
stressFile = "pac.%s_%s.nc" % (gridName, periodName)

# Read sst
dset = Dataset(sstFile, "r")
sst = dset.variables["sst"][:]
nt = sst.shape[0]
dset.close()

# Remove means
sstm = np.mean(sst, 0)
Wum = np.mean(Wu, 0)
Wvm = np.mean(Wv, 0)
for t in np.arange(nt):
    sst[t] -= sstm
    Wu[t] -= Wum
    Wv[t] -= Wvm

