import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.basemap import Basemap, addcyclic
from scipy.io import FortranFile

# Amplification factor for the mean wind stress
ampMean = 3.0

initDir = '../init/'
nlat = 31
nlon = 30
year0 = 1961
yearf = 1994
gridName = "%dx%d" % (nlat, nlon)
periodName = "%d%d" % (year0, yearf)
postfix = '%s_%s' % (gridName, periodName)

sstFile = "ersst.%s.nc"  % postfix
psFile = "pac.%s.nc" % postfix

T0 = 30.
rhoAir = 1.2
CD = 1.2e-3 # Large & Pond 1981 for w < 11 m/s
L = 1.5e7
c0 = 2.
rho = 1024
H = 200.
toStress = rhoAir * CD
toND = L / (c0**2 * rho * H) # = F0 / tau0

# Read cartesian coordinates
x = np.loadtxt('%s/x_%s.txt' % (initDir, gridName))
y = np.loadtxt('%s/y_%s.txt' % (initDir, gridName))
(X2, Y2) = np.meshgrid(x, y)

# Read sst
dset = Dataset(sstFile, "r")
sst = dset.variables["sst"][:]
lat = dset.variables["lat"][:]
lon = dset.variables["lon"][:]
(LON, LAT) = np.meshgrid(lon, lat)
nt = sst.shape[0]
dset.close()

# Map definition
llcrnrlon = lon.min()
llcrnrlat = lat.min()
urcrnrlon = lon.max()
urcrnrlat = lat.max()
nlev = 10
map = Basemap(projection='merc', llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, resolution='c')
(x, y) = map(LON, LAT)

# Read zonal pseudo wind-stress
dset = Dataset(psFile, "r")
Wu = dset.variables["Wu"][:]
dset.close()

# Get mask
N = nlat * nlon
sst = sst.reshape(nt, N)
Wu = Wu.reshape(nt, N)
mask = np.any(sst.mask, 0) | np.any(Wu.mask, 0)
sstMasked = np.array(sst[:, ~mask])
WuMasked = np.array(Wu[:, ~mask])
lonlMasked = LON.flatten()[~mask]
latlMasked = LAT.flatten()[~mask]
nValid = N - mask.sum()

# Remove radiative equilibrium temperature : Ta = T - T0
sstMasked -= T0

# Remove mean, remove seasonal cycle, put back mean
print 'Getting anomalies...'
sstMean = sstMasked.mean(0)
WuMean = WuMasked.mean(0)
sstMaskedAnom = sstMasked - np.tile(np.expand_dims(sstMean, 0), (nt, 1))
WuMaskedAnom = WuMasked - np.tile(np.expand_dims(WuMean, 0), (nt, 1))
ssta = np.copy(sstMaskedAnom,)
Wua = np.copy(WuMaskedAnom)
for k in np.arange(12):
    ssta[k::12] -= np.tile(np.expand_dims(sstMaskedAnom[k::12].mean(0), 0),
                           (sstMaskedAnom[k::12].shape[0], 1))
    Wua[k::12] -= np.tile(np.expand_dims(WuMaskedAnom[k::12].mean(0), 0),
                          (WuMaskedAnom[k::12].shape[0], 1))
ssta += np.tile(np.expand_dims(sstMean, 0), (nt, 1))
Wua += np.tile(np.expand_dims(WuMean, 0), (nt, 1))

# Regressions
print 'Getting wind stress residuals...'
WuResDim = np.copy(Wua)
for ij in np.arange(nValid):
    X = np.matrix(ssta[:, ij]).T
    Y = np.matrix(Wua[:, ij]).T
    A = (X.T * X)**(-1) * (X.T * Y)
    WuResDim[:, ij] = np.squeeze(np.array(Y - X * A))
        
# Adimensionalize ! F0=L tau0/(co^2 rho H)  
print 'Get adimensional residual wind stress...'
WuRes = WuResDim * toStress * toND

# Decompose residual wind-stress
WuResMean = WuRes.mean(0)
WuResAnom = WuRes - np.tile(np.expand_dims(WuResMean, 0), (nt, 1))

# Plot WuResMean
fig = plt.figure()
field = np.ma.masked_all((nlat*nlon,))
field[~mask] = WuResMean
field = field.reshape(nlat, nlon)
vmax = np.max(field)
vmin = np.min(field)
levels = np.linspace(vmin, vmax, nlev)
cs = map.contourf(x, y, field, levels, cmap=cm.RdBu_r)
plt.title('Mean of the residual wind-stress')
map.drawcoastlines()
# draw parallels and meridians.
map.drawparallels(np.arange(0, 81.,10.))
map.drawmeridians(np.arange(-180.,181.,30.))
plt.colorbar(orientation='horizontal')
fig.savefig('tauResMean_%s.png' % postfix, bbox_inches='tight') 

# Plot std of WuRes
fig = plt.figure()
field = np.ma.masked_all((nlat*nlon,))
field[~mask] = WuResAnom.std(0)
field = field.reshape(nlat, nlon)
vmax = np.max(field)
vmin = 0.
levels = np.linspace(vmin, vmax, nlev)
cs = map.contourf(x, y, field, levels, cmap=cm.hot_r)
plt.title('Std of residual wind-stress')
map.drawcoastlines()
# draw parallels and meridians.
map.drawparallels(np.arange(0, 81.,10.))
map.drawmeridians(np.arange(-180.,181.,30.))
plt.colorbar(orientation='horizontal')
fig.savefig('tauResStd_%s.png' % postfix, bbox_inches='tight') 

# EOF Decomposition
print 'EOF decomposition...'
# Get covariance matrix
sim = np.cov(WuResAnom, rowvar=False)
# Eigenvector decomposition
(w, v) = np.linalg.eigh(sim)
# Get principal component
pc = np.dot(WuResAnom, v)
isort = np.argsort(w)
w = w[isort][::-1]
v = v[:, isort][:, ::-1]
pc = pc[:, isort][:, ::-1]
wn = w / w.sum()

# Plot first EOFs
nEOF = 1
#nEOF = 3
print 'First %d EOFs explained variance: ' % nEOF, (wn[:nEOF] * 100).astype(int)
print 'First %d EOFs cumulated explained variance: ' % nEOF, (wn[:nEOF].cumsum() * 100).astype(int)
for k in np.arange(nEOF):
    fig = plt.figure()
    eof = np.ma.masked_all((nlat*nlon,))
    eof[~mask] = v[:, k]
    eof = eof.reshape(nlat, nlon)
    vmax = np.max(np.abs(eof))
    vmin = -vmax
    levels = np.linspace(vmin, vmax, nlev)
    cs = map.contourf(x, y, eof, levels, cmap=cm.RdBu_r)
    plt.title('EOF #' + str(k) + " explaining %2d%% of variance" % (wn[k] * 100,))
    map.drawcoastlines()
    # draw parallels and meridians.
    map.drawparallels(np.arange(0, 81.,10.))
    map.drawmeridians(np.arange(-180.,181.,30.))
    plt.colorbar(orientation='horizontal')
    fig.savefig('tau_eof%d_%s.png' % (k, postfix), bbox_inches='tight')

# Get Periodograms
T = 100000
sampPeriod = 1.
window = np.hamming(nt)
# Get nearest larger power of 2
if np.log2(nt) != int(np.log2(nt)):
    nfft = 2**(int(np.log2(nt)) + 1)
else:
    nfft = nt
# Get frequencies and shift zero frequency to center
freq = np.fft.fftfreq(nfft, d=sampPeriod)
freq = np.fft.fftshift(freq)
freqYear = freq * 12

windowrn = np.hamming(T)
if np.log2(T) != int(np.log2(T)):
    nfftrn = 2**(int(np.log2(T)) + 1)
else:
    nfftrn = T
# Get frequencies and shift zero frequency to center
#freqrn = np.fft.fftfreq(nfftrn, d=sampPeriod)
freqrn = np.fft.fftfreq(nfftrn, d=sampPeriod)
freqrn = np.fft.fftshift(freqrn)
freqrnYear = freqrn * 12
    
nRAVG = 5
nRAVGrn = int(nRAVG * nfftrn * 1. / nfft * 0.1)

# Get NINO4
print 'Getting NINO4 index...'
nino4slat = -5.
nino4nlat = 5.
nino4wlon = 160.
nino4elon = 210.
nino4 = WuResAnom[:, (lonlMasked >= nino4wlon) & (lonlMasked <= nino4elon)
              & (latlMasked >= nino4slat) & (latlMasked <= nino4nlat)].mean(1) / toND

# Get periodogram of zonal wind stress averaged over nino4
ts = nino4# / nt
# Apply window
tsWindowed = ts * window
# Fourier transform and shift zero frequency to center
fts = np.fft.fft(tsWindowed, nfft, 0)
fts = np.fft.fftshift(fts)
# Get periodogram
perio = np.abs(fts / nt)#**2

# Apply running average
perioRAVG = perio.copy()
for k in np.arange(nRAVG/2, nfft-nRAVG/2):
    perioRAVG[k] = perio[k-nRAVG/2:k+nRAVG/2 + 1].mean() / nRAVG

# Plot
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(freqYear, np.log10(perioRAVG))
#    ax.set_xscale('log')
#    ax.set_yscale('log')
ax.set_xlim(0, 4)
#ax.set_ylim(0, vmax)
plt.title('Periodogram of residual wind stress averaged over Nino4')
fig.savefig('nino4_perio_%s.png' % postfix, bbox_inches='tight')

# Get red-noise parameters
print 'Getting red-noise parameters of principal components...'
rn = np.empty((2, nEOF))
for k in np.arange(nEOF):
    ts = pc[:, k]
    # Get autocorrelation
    rn[0, k] = np.corrcoef(ts[1:], ts[:-1])[0, 1]
    rn[1, k] = ts.std() * np.sqrt(1. - rn[0, k]**2)

# Generate the red noise time series and plot the FFTs
print 'Generating red noise principal components...'
pcrn = np.empty((T, nEOF))
for k in np.arange(nEOF):
    pcrn[0, k] = np.random.normal(0, rn[1, k])
    for t in np.arange(1, T):
        pcrn[t, k] = rn[0, k] * pcrn[t-1, k] + np.random.normal(0, rn[1, k])

# Plot FFTs
print 'Plotting periodograms...'
for k in np.arange(nEOF):
    ts = pc[:, k]
    # FFT
    # Apply window
    tsWindowed = ts * window
    # Fourier transform and shift zero frequency to center
    fts = np.fft.fft(tsWindowed, nfft, 0)
    fts = np.fft.fftshift(fts)
    # Get periodogram
    perio = np.abs(fts / nt)**2
    # Apply running average
    perioRAVG = perio.copy()
    for i in np.arange(nRAVG/2, nfft-nRAVG/2):
        perioRAVG[i] = perio[i-nRAVG/2:i+nRAVG/2 + 1].mean() / nRAVG

    tsrn = pcrn[:, k]
    # FFT
    # Apply window
    tsrnWindowed = tsrn * windowrn
    # Fourier transform and shift zero frequency to center
    ftsrn = np.fft.fft(tsrnWindowed, nfftrn, 0)
    ftsrn = np.fft.fftshift(ftsrn)
    # Get periodogram
    periorn = np.abs(ftsrn / T)**2
    # Apply running average
    periornRAVG = periorn.copy()
    for i in np.arange(nRAVGrn/2, nfftrn-nRAVGrn/2):
        periornRAVG[i] = periorn[i-nRAVGrn/2:i+nRAVGrn/2 + 1].mean() / nRAVGrn


    Prn = rn[1, k] * (1. - rn[0, k]**2) / (1. - 2.*rn[0, k]*np.cos(freq*2*np.pi) + rn[0, k]**2) / nt**2

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(freqYear, np.log10(perioRAVG), '-k', linewidth=2)
    ax.plot(freqrnYear, np.log10(periornRAVG), ':k', linewidth=2)
    ax.plot(freqYear, np.log10(Prn), '--k')
    #    ax.set_xscale('log')
#    ax.set_yscale('log')
    ax.set_xlim(0, 4)
    #ax.set_ylim(0, vmax)
    plt.title('Periodogram of principal component %d' % k)
    plt.legend(('Periodogram of pc %d' % k,
                'Periodogram of red-noise fitted to pc %d' % k,
                'Theoretical periodogram of red-noise fitted to pc %d' % k), loc='lower left')
    fig.savefig('perio_pc%d_%s.png' % (k, postfix),
                bbox_inches='tight')



# Write red noise time series and its derivative w.r.t x and y
# Use a centered order 1 skim
print 'Writing adimensional residual stochastic wind stress...'
X3 = np.tile(np.expand_dims(X2, 0), (T, 1, 1))
dX3Center = X3[:, :, 2:] - X3[:, :, :-2]
dX3Left = X3[:, :, 1] - X3[:, :, 0]
dX3Right = X3[:, :, -1] - X3[:, :, -2]
Y3 = np.tile(np.expand_dims(Y2, 0), (T, 1, 1))
dY3Center = Y3[:, 2:] - Y3[:, :-2]
dY3Left = Y3[:, 1] - Y3[:, 0]
dY3Right = Y3[:, -1] - Y3[:, -2]

# Write mean tau
WuResMeanValid = np.zeros((nlat * nlon,))
WuResMeanValid[~mask] = WuResMean
WuResMeanValid = WuResMeanValid.reshape(nlat, nlon)
tau = np.tile(np.expand_dims(WuResMeanValid, 0), (T, 1, 1)) * ampMean
dtaudx = np.empty(tau.shape)
dtaudx[:, :, 1:-1] = (tau[:, :, 2:] - tau[:, :, :-2]) / dX3Center
dtaudx[:, :, 0] = (tau[:, :, 1] - tau[:, :, 0]) / dX3Left
dtaudx[:, :, -1] = (tau[:, :, -1] - tau[:, :, -2]) / dX3Right
dtaudy = np.empty(tau.shape)
dtaudy[:, 1:-1] = (tau[:, 2:] - tau[:, :-2]) / dY3Center
dtaudy[:, 0] = (tau[:, 1] - tau[:, 0]) / dY3Left
dtaudy[:, -1] = (tau[:, -1] - tau[:, -2]) / dY3Right
f = FortranFile('%s/tau_mean_ampMean%02d%s.bin' % (initDir, int(ampMean * 10), postfix), 'w')
for t in np.arange(T):
    f.write_record(tau[t])
f.close()
f = FortranFile('%s/dtaudx_mean_ampMean%02d%s.bin' % (initDir, int(ampMean * 10), postfix), 'w')
for t in np.arange(T):
    f.write_record(dtaudx[t])
f.close()
f = FortranFile('%s/dtaudy_mean_ampMean%02d%s.bin' % (initDir, int(ampMean * 10), postfix), 'w')
for t in np.arange(T):
    f.write_record(dtaudy[t])
f.close()

# Write tau with variations
for k in np.arange(nEOF):
    print 'Writing tau with %d eofs...' % (k+1,)
    eof = np.zeros((nlat*nlon,))
    eof[~mask] = v[:, k]
    eof = eof.reshape(nlat, nlon)
    tau += np.tile(np.expand_dims(np.expand_dims(pcrn[:, k], 1), 2),
                   (1, nlat, nlon)) \
        * np.tile(np.expand_dims(eof, 0), (T, 1, 1))
    # Normalize
    tau =  tau / np.tile(tau.std(0).mean(), (T, 1, 1))
    dtaudx = np.empty(tau.shape)
    dtaudx[:, :, 1:-1] = (tau[:, :, 2:] - tau[:, :, :-2]) / dX3Center
    dtaudx[:, :, 0] = (tau[:, :, 1] - tau[:, :, 0]) / dX3Left
    dtaudx[:, :, -1] = (tau[:, :, -1] - tau[:, :, -2]) / dX3Right
    dtaudy = np.empty(tau.shape)
    dtaudy[:, 1:-1] = (tau[:, 2:] - tau[:, :-2]) / dY3Center
    dtaudy[:, 0] = (tau[:, 1] - tau[:, 0]) / dY3Left
    dtaudy[:, -1] = (tau[:, -1] - tau[:, -2]) / dY3Right
    f = FortranFile('%s/tauNorm_%deofs_ampMean%02d%s.bin' \
                    % (initDir, k+1, int(ampMean * 10), postfix), 'w')
    for t in np.arange(T):
        f.write_record(tau[t])
    f.close()
    f = FortranFile('%s/dtaudxNorm_%deofs_ampMean%02d%s.bin' \
                    % (initDir, k+1, int(ampMean * 10), postfix), 'w')
    for t in np.arange(T):
        f.write_record(dtaudx[t])
    f.close()
    f = FortranFile('%s/dtaudyNorm_%deofs_ampMean%02d%s.bin' % \
                    (initDir, k+1, int(ampMean * 10), postfix), 'w')
    for t in np.arange(T):
        f.write_record(dtaudy[t])
    f.close()

    # # Plot
    # fig = plt.figure()
    # field = tau.mean(0)
    # vmax = np.max(field)
    # vmin = np.min(field)
    # levels = np.linspace(vmin, vmax, nlev)
    # cs = map.contourf(x, y, field, levels, cmap=cm.RdBu_r)
    # plt.title('RN Mean with %d EOFs' % (k,))
    # map.drawcoastlines()
    # # draw parallels and meridians.
    # map.drawparallels(np.arange(0, 81.,10.))
    # map.drawmeridians(np.arange(-180.,181.,30.))
    # plt.colorbar(orientation='horizontal')

    # fig = plt.figure()
    # field = tau.std(0)
    # vmax = np.max(field)
    # vmin = np.min(field)
    # levels = np.linspace(vmin, vmax, nlev)
    # cs = map.contourf(x, y, field, levels, cmap=cm.RdBu_r)
    # plt.title('RN std with %d EOFs' % (k,))
    # map.drawcoastlines()
    # # draw parallels and meridians.
    # map.drawparallels(np.arange(0, 81.,10.))
    # map.drawmeridians(np.arange(-180.,181.,30.))
    # plt.colorbar(orientation='horizontal')
