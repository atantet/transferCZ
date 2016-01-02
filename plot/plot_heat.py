import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.basemap import Basemap

initDir = '../init/'
nlat = 31
nlon = 30
year0 = 1961
yearf = 1994
gridName = "%dx%d" % (nlat, nlon)
periodName = "%d%d" % (year0, yearf)
postfix = '%s_%s' % (gridName, periodName)

# Dimensions
L = 1.5e7
c0 = 2
timeDim = L / c0 / (60. * 60. * 24)

# Case definition
#mu = 2.7
#mu = 3.0
# mu = 3.3
mu = 7.0
# simType = '_short'
simType = '_deter'
caseDir = '../'
caseSubDir = 'mu%d%s/' % (int(mu * 10), simType)
casePath = '%s/%s/' % (caseDir, caseSubDir)

pltDir = caseSubDir
os.system('mkdir %s 2> /dev/null' % pltDir)

# Dataset definition
fieldName = 'Heat Content'
fieldFile = 'heat'
caseFile = '%s_%s.txt' % (fieldFile, postfix)

# # Read cartesian coordinates
lat = np.loadtxt('%s/lat_%s.txt' % (initDir, gridName))
lon = np.loadtxt('%s/lon_%s.txt' % (initDir, gridName))
(LON, LAT) = np.meshgrid(lon, lat)
llcrnrlon=lon.min()
llcrnrlat=lat.min()
urcrnrlon=lon.max()
urcrnrlat=lat.max()
map = Basemap(projection='merc', llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, resolution='c')
(X, Y) = map(LON, LAT)

# Read dataset
print 'Reading dataset %s/%s...' % (casePath, caseFile)
var = np.loadtxt('%s/%s' % (casePath, caseFile))
nt = var.shape[0]
var = var.reshape(nt, nlat, nlon)
# Read time
dsetFile = 'fort.41'
dset = np.loadtxt('%s/%s' % (casePath, dsetFile))
timeND = dset[:, 1]
nt = timeND.shape[0] / nlon
time = timeND.reshape(nt, nlon)[:, 0] * timeDim

# Plot field statistics
nlev = 20
fig = plt.figure()
field = hc.mean(0)
vmax = np.max(field)
vmin = np.min(field)
levels = np.linspace(vmin, vmax, nlev)
cs = map.contourf(X, Y, field, levels, cmap=cm.RdBu_r)
map.drawcoastlines()
map.drawparallels(np.arange(0, 81.,10.))
map.drawmeridians(np.arange(-180.,181.,30.))
plt.title('Mean of %s' % fieldName)
plt.colorbar(orientation='horizontal')
fig.savefig('%s/%sMean_%s.png' % (pltDir, fieldFile, postfix),
            bbox_inches='tight')

fig = plt.figure()
field = var.std(0)
vmax = np.max(field)
vmin = np.min(field)
levels = np.linspace(vmin, vmax, nlev)
cs = map.contourf(X, Y, field, levels, cmap=cm.RdBu_r)
map.drawcoastlines()
map.drawparallels(np.arange(0, 81.,10.))
map.drawmeridians(np.arange(-180.,181.,30.))
plt.title('Std of %s' % fieldName)
plt.colorbar(orientation='horizontal')
fig.savefig('%s/%sStd_%s.png' % (pltDir, fieldFile, postfix),
            bbox_inches='tight')

# Plot total heat content
totalHc = var.mean(1).mean(1)

linewidth = 2
fig = plt.figure()
plt.plot(time[:200], totalHc[:200], linewidth=linewidth)
plt.title('Time-series of the total %s' % fieldName)
fig.savefig('%s/total%s_%s.png' % (pltDir, fieldName, postfix),
            bbox_inches='tight')
