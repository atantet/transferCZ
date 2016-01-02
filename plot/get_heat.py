import numpy as np

initDir = '../init/'
nlat = 31
nlon = 30
year0 = 1961
yearf = 1994
gridName = "%dx%d" % (nlat, nlon)
periodName = "%d%d" % (year0, yearf)
postfix = '%s_%s' % (gridName, periodName)

# Case definition
#mu = 2.7
#mu = 3.0
#mu = 3.3
mu = 7.0
#simType = '_short'
simType = '_deter'

caseDir = '../'
caseSubDir = 'mu%d%s/' % (int(mu * 10), simType)
casePath = '%s/%s/' % (caseDir, caseSubDir)

# Dataset definition
dsetFile = 'fort.49'

hCol = 3   #  h
hName = 'Thermocline depth'
hFile = 'h'
TCol = 4    # T
TName = 'SST'
TFile = 'T'
# Read dataset
dset = np.loadtxt('%s/%s' % (casePath, dsetFile))
time = dset[:, 0]
nt = time.shape[0] / (nlat * nlon)
time = time.reshape(nt, nlon, nlat)[:, 0, 0]
# X = dset[:, 1].reshape(nt, nlon, nlat)[0].T
# Y = dset[:, 2].reshape(nt, nlon, nlat)[0].T
h = dset[:, hCol].reshape(nt, nlon, nlat)
T = dset[:, TCol].reshape(nt, nlon, nlat)

# Heat content
hc = h * T
hc = hc.swapaxes(1, 2)

# Save
np.savetxt('%s/heat_%s.txt' % (casePath, postfix), hc.reshape(nt, nlat*nlon))
