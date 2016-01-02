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

# Case definition
#mu = 2.7
#mu = 3.0
#mu = 3.3
#mu = 3.5
mu = 4.0
#mu = 4.5
#mu = 5.0
#mu = 5.5
#mu = 6.0
#mu = 7.0
#ampMean = 2
#ampMean = 2.5
ampMean = 3.0
# simType = '_short'
#simType = '_veryshort'
#simType = '_deter'
simType = '_amp%02d_deter' % (int(ampMean * 10),)
caseDir = '../'
caseSubDir = 'mu%d%s/' % (int(mu * 10), simType)
casePath = '%s/%s/' % (caseDir, caseSubDir)

pltDir = caseSubDir
os.system('mkdir %s 2> /dev/null' % pltDir)

# Dataset definition
dsetFile = 'fort.49'

# fieldDef = (3, 'Thermocline depth', 'h') #  h
# fieldDef = (4, 'SST', 'T') # T
# fieldDef = (5, 'Zonal Wind', 'u_A')
# fieldDef = (6, 'Radiative equilibrium temperature', 'T0')
fieldDef = (7, 'Total wind-stress', 'wind')
# fieldDef = (8, 'External wind-stress', 'taux')

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
dset = np.loadtxt('%s/%s' % (casePath, dsetFile))
time = dset[:, 0]
nt = time.shape[0] / (nlat * nlon)
time = time.reshape(nt, nlon, nlat)[:, 0, 0]
# X = dset[:, 1].reshape(nt, nlon, nlat)[0].T
# Y = dset[:, 2].reshape(nt, nlon, nlat)[0].T
var = dset[:, fieldDef[0]].reshape(nt, nlon, nlat)
var = var.swapaxes(1, 2)

# Remove spinup
spinup = 12 * 6
var = var[spinup:]
time = var[spinup:]
nt = time.shape[0]

nlev = 20
fig = plt.figure()
field = var.mean(0)
vmax = np.max(field)
vmin = np.min(field)
levels = np.linspace(vmin, vmax, nlev)
cs = map.contourf(X, Y, field, levels, cmap=cm.RdBu_r)
map.drawcoastlines()
map.drawparallels(np.arange(0, 81.,10.))
map.drawmeridians(np.arange(-180.,181.,30.))
plt.title('Mean of %s' % fieldDef[1])
plt.colorbar(orientation='horizontal')
fig.savefig('%s/%sMean_%s.png' % (pltDir, fieldDef[2], postfix),
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
plt.title('Std of %s' % fieldDef[1])
plt.colorbar(orientation='horizontal')
fig.savefig('%s/%sStd_%s.png' % (pltDir, fieldDef[2], postfix),
            bbox_inches='tight')
