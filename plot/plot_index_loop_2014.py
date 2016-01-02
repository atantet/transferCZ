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

# Dimensions
L = 1.5e7
c0 = 2
timeDim = L / c0 / (60. * 60. * 24)

# Case definition
muRng = np.empty((17,))
muRng[:3] = np.arange(2., 3.1, 0.5)
muRng[3:15] = np.arange(3.1, 4.21, 0.1)
muRng[15:] = np.arange(4.5, 5.1, 0.5)

#ampMean = 2
#ampMean = 2.5
ampMean = 3.0
neof = 1
# simType = '_short'
#simType = '_veryshort'
#simType = '_deter'
#simType = '_amp%02d' % (int(ampMean * 10),)
simType = '_amp%02d_%deof' % (int(ampMean * 10), neof)
#simType = '_amp%02d_deter' % (int(ampMean * 10),)
#simType = ''
caseDir = '../'

# Dataset definition
dsetFile = 'fort.49'

field_h = (3, 'Thermocline depth', 'h')
field_T = (4, 'SST', 'T')
field_u_A = (5, 'Zonal Wind', 'u_A')
field_T0 = (6, 'Radiative equilibrium temperature', 'T0')
field_wind = (7, 'Total wind-stress', 'wind')
field_taux = (8, 'External wind-stress', 'taux')

# Indices definition
# Nino3
nino3Def = ([210., -5., 270., 5.], 'NINO3', 'nino3')
# Nino4
nino4Def = ([160., -5., 210., 5.], 'NINO4', 'nino4')

# Field and index choice
#fieldDef = field_T
#indexDef = nino3Def
fieldDef = field_h
indexDef = nino4Def

# # Read cartesian coordinates
lat = np.loadtxt('%s/lat_%s.txt' % (initDir, gridName))
lon = np.loadtxt('%s/lon_%s.txt' % (initDir, gridName))
(LON, LAT) = np.meshgrid(lon, lat)

for k in np.arange(muRng.shape[0]):
    mu = muRng[k]
    caseSubDir = 'mu%d%s/' % (int(mu * 10), simType)
    casePath = '%s/%s/' % (caseDir, caseSubDir)
    pltDir = caseSubDir
    os.system('mkdir %s 2> /dev/null' % pltDir)

    # Read dataset
    print 'Reading dataset %s/%s' % (casePath, dsetFile)
    dset = np.loadtxt('%s/%s' % (casePath, dsetFile))
    timeND = dset[:, 0]
    nt = timeND.shape[0] / (nlat * nlon)
    timeFull = timeND.reshape(nt, nlon, nlat)[:, 0, 0] * timeDim / 365
    var = dset[:, fieldDef[0]].reshape(nt, nlon, nlat)
    var = var.swapaxes(1, 2)

    # Index calculation
    print 'Getting %s index...' % (nino3Def[1],)
    ilatIndex = (lat >= indexDef[0][1]) & (lat <= indexDef[0][3])
    ilonIndex = (lon >= indexDef[0][0]) & (lon <= indexDef[0][2])
    indexFull = var[:, ilatIndex][:, :, ilonIndex].mean(1).mean(1)
    # Save
    np.savetxt('%s/%s/%s%s_%s.txt'\
               % (caseDir, caseSubDir, indexDef[2], fieldDef[2], postfix),
               indexFull)

    # Remove spinup
    spinup = 12 * 18
    #spinup = 0
    index = indexFull[spinup:]
    time = timeFull[spinup:]
    nt = time.shape[0]
    
    # Plot time-series
    linewidth = 2
    fig = plt.figure()
    plt.plot(time, index, linewidth=linewidth)
    plt.title('Time-series of %s averaged over %s' % (fieldDef[1], indexDef[1]))
    fig.savefig('%s/%s%s_%s.png' % (pltDir, indexDef[2], fieldDef[2], postfix),
                bbox_inches='tight')
    
    
    # Get periodogram of zonal wind stress averaged over index
    sampPeriod = 1.
    nRAVG = 1
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
    
    ts = index - index.mean(0)
    # Apply window
    tsWindowed = ts * window
    # Fourier transform and shift zero frequency to center
    fts = np.fft.fft(tsWindowed, nfft, 0)
    fts = np.fft.fftshift(fts)
    # Get periodogram
    perio = np.abs(fts / nt)**2
    
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
    plt.title('Periodogram of %s averaged over %s' % (fieldDef[1], indexDef[1]))
    fig.savefig('%s/%s%sPerio_%s.png' % (pltDir, indexDef[2], fieldDef[2], postfix),
                bbox_inches='tight')
