import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors

initDir = '../init/'
nlat = 31
nlon = 30
year0 = 1961
yearf = 1994
gridName = "%dx%d" % (nlat, nlon)
periodName = "%d%d" % (year0, yearf)
postfix = '%s_%s' % (gridName, periodName)

L = 1.5e7
c0 = 2
timeDim = L / c0 / (60. * 60. * 24)

# Case definition
mu = 2.7
#mu = 3.0
#mu = 3.3
simLength = '_veryshort'
caseDir = '../'
caseSubDir = 'mu%d%s/' % (int(mu * 10), simLength)
casePath = '%s/%s/' % (caseDir, caseSubDir)

pltDir = caseSubDir
os.system('mkdir %s 2> /dev/null' % pltDir)

# Dataset definition
dsetFile = 'fort.41'

# fieldCol = 2   # w1
# fieldName = 'Upwelling'
# fieldFile = 'w1'
# fieldCol = 3   #  h
# fieldName = 'Thermocline depth'
# fieldFile = 'h'
fieldCol = 4    # T - T0
fieldName = 'SST Anomaly'
fieldFile = 'T'
# fieldCol = 5    # u_A
# fieldName = 'Zonal wind'
# fieldFile = 'u_A'
# fieldCol = 6    # wind
# fieldName = 'Total wind-stress'
# fieldFile = 'wind'

# # Read cartesian coordinates
lon = np.loadtxt('%s/lon_%s.txt' % (initDir, gridName))

# Read dataset
dset = np.loadtxt('%s/%s' % (casePath, dsetFile))
timeND = dset[:, 1]
nt = timeND.shape[0] / nlon
time = timeND.reshape(nt, nlon)[:, 0] * timeDim
# X = dset[:, 0].reshape(nt, nlon, nlat)[0].T
var = dset[:, fieldCol].reshape(nt, nlon)
(X, T) = np.meshgrid(lon, time)

# Plot stats
linewidth = 2

fig = plt.figure()
profile = var.mean(0)
plt.plot(lon, profile, linewidth=linewidth)
plt.title('Mean of equatorial %s' % fieldName)
fig.savefig('%s/eq%sMean_%s.png' % (pltDir, fieldFile, postfix),
            bbox_inches='tight')

fig = plt.figure()
profile = var.std(0)
plt.plot(lon, profile, linewidth=linewidth)
plt.title('Std of equatorial %s' % fieldName)
fig.savefig('%s/eq%sStd_%s.png' % (pltDir, fieldFile, postfix),
            bbox_inches='tight')

# Plot Hovmoller diagram
nlev = 20
fig = plt.figure()
vmax = np.max(np.abs(var))
vmin = -vmax
#vmax = np.max(var)
#vmin = np.min(var)
levels = np.linspace(vmin, vmax, nlev)
cs = plt.contourf(X, T, var, levels, cmap=cm.RdBu_r)
plt.title('Hovmoller diagram of %s' % fieldName)
plt.colorbar()
#plt.colorbar(orientation='horizontal')
fig.savefig('%s/eq%sHov_%s.png' % (pltDir, fieldFile, postfix),
            bbox_inches='tight')
