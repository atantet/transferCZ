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
dstPostfix = ''
noiseType = 'White'

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

# time
T = 600 # duration in years
seed = 4
nBlocks = 1
dt = 0.060 # adimensional time-step
adim2Sec = L / c0 # Conversion from adimensional to seconds
sec2Years =  1. / (60*60*24*365)
adim2Years = adim2Sec * sec2Years

# Read cartesian coordinates
x = np.loadtxt('%s/lon_%s.txt' % (initDir, gridName))
y = np.loadtxt('%s/lat_%s.txt' % (initDir, gridName))
(X2, Y2) = np.meshgrid(x, y)

# Read sst
dset = Dataset(sstFile, "r")
sst = dset.variables["sst"][:]
lat = dset.variables["lat"][:]
lon = dset.variables["lon"][:]
(LON, LAT) = np.meshgrid(lon, lat)
ntObs = sst.shape[0]
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
sst = sst.reshape(ntObs, N)
Wu = Wu.reshape(ntObs, N)
mask = np.any(sst.mask, 0) | np.any(Wu.mask, 0)
sstMasked = np.array(sst[:, ~mask])
WuMasked = np.array(Wu[:, ~mask])
lonlMasked = LON.flatten()[~mask]
latlMasked = LAT.flatten()[~mask]
nValid = N - mask.sum()

# Remove mean, remove seasonal cycle, put back mean
print 'Getting anomalies...'
sstMean = sstMasked.mean(0)
WuMean = WuMasked.mean(0)
sstMaskedAnom = sstMasked - np.tile(np.expand_dims(sstMean, 0), (ntObs, 1))
WuMaskedAnom = WuMasked - np.tile(np.expand_dims(WuMean, 0), (ntObs, 1))
ssta = np.copy(sstMaskedAnom,)
Wua = np.copy(WuMaskedAnom)
for k in np.arange(12):
    ssta[k::12] -= np.tile(np.expand_dims(sstMaskedAnom[k::12].mean(0), 0),
                           (sstMaskedAnom[k::12].shape[0], 1))
    Wua[k::12] -= np.tile(np.expand_dims(WuMaskedAnom[k::12].mean(0), 0),
                          (WuMaskedAnom[k::12].shape[0], 1))

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

# Plot std of WuRes
fig = plt.figure()
field = np.ma.masked_all((nlat*nlon,))
field[~mask] = WuRes.std(0)
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
fig.savefig('tauResStd%s.png' % dstPostfix, bbox_inches='tight') 

# EOF Decomposition
print 'EOF decomposition...'
# Get covariance matrix
sim = np.cov(WuRes, rowvar=False)
# Eigenvector decomposition
(w, v) = np.linalg.eigh(sim)
# Get principal component
pc = np.dot(WuRes, v)
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
    fig.savefig('tau_eof%d%s.png' % (k, dstPostfix), bbox_inches='tight')

# Get Periodograms
sampPeriod = 1.
window = np.hamming(ntObs)
# Get nearest larger power of 2
if np.log2(ntObs) != int(np.log2(ntObs)):
    nfft = 2**(int(np.log2(ntObs)) + 1)
else:
    nfft = ntObs
# Get frequencies and shift zero frequency to centObser
freq = np.fft.fftfreq(nfft, d=sampPeriod)
freq = np.fft.fftshift(freq)
freqYear = freq * 12

nt =  int(T / (dt * adim2Years))
windowrn = np.hamming(nt)
if np.log2(nt) != int(np.log2(nt)):
    nfftrn = 2**(int(np.log2(nt)) + 1)
else:
    nfftrn = nt
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
nino4 = WuRes[:, (lonlMasked >= nino4wlon) & (lonlMasked <= nino4elon)
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
fig.savefig('nino4_perio%s.png' % dstPostfix, bbox_inches='tight')

# Get white noise parameters (std)
print 'Getting white-noise parameters of principal components...'
rn = np.empty((1, nEOF))
for k in np.arange(nEOF):
    rn[0, k] = pc[:, k].std()

# # Plot FFTs
# print 'Plotting periodograms...'
# for k in np.arange(nEOF):
#     ts = pc[:, k]
#     # FFT
#     # Apply window
#     tsWindowed = ts * window
#     # Fourier transform and shift zero frequency to center
#     fts = np.fft.fft(tsWindowed, nfft, 0)
#     fts = np.fft.fftshift(fts)
#     # Get periodogram
#     perio = np.abs(fts / nt)**2
#     # Apply running average
#     perioRAVG = perio.copy()
#     for i in np.arange(nRAVG/2, nfft-nRAVG/2):
#         perioRAVG[i] = perio[i-nRAVG/2:i+nRAVG/2 + 1].mean() / nRAVG

#     tsrn = pcwn[:, k]
#     # FFT
#     # Apply window
#     tsrnWindowed = tsrn * windowrn
#     # Fourier transform and shift zero frequency to center
#     ftsrn = np.fft.fft(tsrnWindowed, nfftrn, 0)
#     ftsrn = np.fft.fftshift(ftsrn)
#     # Get periodogram
#     periorn = np.abs(ftsrn / nt)**2
#     # Apply running average
#     periornRAVG = periorn.copy()
#     for i in np.arange(nRAVGrn/2, nfftrn-nRAVGrn/2):
#         periornRAVG[i] = periorn[i-nRAVGrn/2:i+nRAVGrn/2 + 1].mean() / nRAVGrn


#     Pwn = 1. / nt**2

#     # Plot
#     fig = plt.figure()
#     ax = fig.add_subplot(1,1,1)
#     ax.plot(freqYear, np.log10(perioRAVG), '-k', linewidth=2)
#     ax.plot(freqrnYear, np.log10(periornRAVG), ':k', linewidth=2)
# #    ax.plot(freqYear, np.log10(Pwn), '--k')
#     #    ax.set_xscale('log')
# #    ax.set_yscale('log')
#     ax.set_xlim(0, 4)
#     #ax.set_ylim(0, vmax)
#     plt.title('Periodogram of principal component %d' % k)
#     plt.legend(('Periodogram of pc %d' % k,
#                 'Periodogram of white-noise fitted to pc %d' % k),
# #                'Theoretical periodogram of white-noise fitted to pc %d' % k),
#                 loc='lower left')
#     fig.savefig('perio_pc%d%s.png' % (k, dstPostfix),
#                 bbox_inches='tight')



ntBlock = int(nt / nBlocks)
ftau = []
fdtaudx = []
fdtaudy = []
# Use a centBlockered order 1 skim
X3 = np.tile(np.expand_dims(X2, 0), (ntBlock, 1, 1))
dX3CentBlocker = X3[:, :, 2:] - X3[:, :, :-2]
dX3Left = X3[:, :, 1] - X3[:, :, 0]
dX3Right = X3[:, :, -1] - X3[:, :, -2]
Y3 = np.tile(np.expand_dims(Y2, 0), (ntBlock, 1, 1))
dY3CentBlocker = Y3[:, 2:] - Y3[:, :-2]
dY3Left = Y3[:, 1] - Y3[:, 0]
dY3Right = Y3[:, -1] - Y3[:, -2]

# Set the seed
np.random.seed(seed)
for k in np.arange(nEOF):
    ftau.append(FortranFile('%s/tau%s_%deofs_amp%02d_seed%d%s.bin' \
                            % (initDir, noiseType, k+1, int(ampMean * 10),
                               seed, dstPostfix), 'w'))
    fdtaudx.append(FortranFile('%s/dtaudx%s_%deofs_amp%02d_seed%d%s.bin' \
                              % (initDir, noiseType, k+1, int(ampMean * 10),
                                 seed, dstPostfix), 'w'))
    fdtaudy.append(FortranFile('%s/dtaudy%s_%deofs_amp%02d_seed%d%s.bin' % \
                              (initDir, noiseType, k+1, int(ampMean * 10),
                               seed, dstPostfix), 'w'))

for block in np.arange(nBlocks):
    print 'Saving block ', (block+1), ' of ', nBlocks
    pcwn = np.empty((ntBlock, nEOF))

    # Write-noise time series and its derivative w.r.t x and y
    # Write tau with variations
    tauNoNorm = np.zeros((ntBlock, nlat, nlon))
    for k in np.arange(nEOF):
        pcwn[:, k] = np.random.normal(0, rn[0, k], size=ntBlock)
        eof = np.zeros((nlat*nlon,))
        eof[~mask] = v[:, k]
        eof = eof.reshape(nlat, nlon)
        tauNoNorm += np.tile(np.expand_dims(np.expand_dims(pcwn[:, k], 1), 2),
                             (1, nlat, nlon)) \
                             * np.tile(np.expand_dims(eof, 0), (ntBlock, 1, 1))
        # Normalize
        tau =  tauNoNorm / np.tile(tauNoNorm.std(0).mean(), (ntBlock, 1, 1))
        dtaudx = np.empty(tau.shape)
        dtaudx[:, :, 1:-1] = (tau[:, :, 2:] - tau[:, :, :-2]) / dX3CentBlocker
        dtaudx[:, :, 0] = (tau[:, :, 1] - tau[:, :, 0]) / dX3Left
        dtaudx[:, :, -1] = (tau[:, :, -1] - tau[:, :, -2]) / dX3Right
        dtaudy = np.empty(tau.shape)
        dtaudy[:, 1:-1] = (tau[:, 2:] - tau[:, :-2]) / dY3CentBlocker
        dtaudy[:, 0] = (tau[:, 1] - tau[:, 0]) / dY3Left
        dtaudy[:, -1] = (tau[:, -1] - tau[:, -2]) / dY3Right
        print 'Writing tau with %d eofs...' % (k+1,)
        for t in np.arange(ntBlock):
            ftau[k].write_record(tau[t])
            fdtaudx[k].write_record(dtaudx[t])
            fdtaudy[k].write_record(dtaudy[t])

for k in np.arange(nEOF):
    ftau[k].close()
    fdtaudx[k].close()
    fdtaudy[k].close()
