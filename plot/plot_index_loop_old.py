import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors

initDir = '../init/'
nlat = 31
nlon = 30

# Dimensions
L = 1.5e7
c0 = 2
timeDim = L / c0 / (60. * 60. * 24)

# Case definition
epsRng =  np.array([0.])
muRng = np.arange(3.6, 4.1, 0.1)

neof = 1

outputDir = '../output/'
dsetFile = 'fort.49'
prefix = 'zc'
simType = '_%deof' % (neof,)
indicesDir = '../observables/'

field_h = (1, 'Thermocline depth', 'h')
field_T = (2, 'SST', 'T')
field_wind = (3, 'Total wind-stress', 'wind')
#field_u_A = (5, 'Zonal Wind', 'u_A')
#field_T0 = (6, 'Radiative equilibrium temperature', 'T0')
#field_taux = (8, 'External wind-stress', 'taux')

# Indices definition
# Nino3
nino3Def = ([210., -5., 270., 5.], 'NINO3', 'nino3')
# Nino4
nino4Def = ([160., -5., 210., 5.], 'NINO4', 'nino4')

# Field and index choice
fieldDef = field_T
indexDef = nino3Def
#fieldDef = field_h
#indexDef = nino4Def

# # Read cartesian coordinates
lat = np.loadtxt('%s/lat_%s.txt' % (initDir, gridName))
lon = np.loadtxt('%s/lon_%s.txt' % (initDir, gridName))
(LON, LAT) = np.meshgrid(lon, lat)

for eps in epsRng:
    for mu in muRng:
        resDir = '%s%s_mu%02d_eps%02d/' % (prefix, simType, mu * 10, eps * 100)
        resPath = '%s/%s/' % (outputDir, resDir)
        indicesPath = '%s/%s' % (indicesDir, resDir)
        pltDir = resDir
        os.system('mkdir %s %s 2> /dev/null' % (pltDir, indicesPath))

        # # Read dataset
        # print 'Reading dataset %s/%s' % (resPath, dsetFile)
        # dset = np.loadtxt('%s/%s' % (resPath, dsetFile))
        # timeND = dset[:, 0]
        # var = dset[:, fieldDef[0]]
        # del dset
        # # Remove possible extra lines if incomplete
        # nt = timeND.shape[0] / (nlat * nlon)
        # timeND = timeND[:nt*nlat*nlon]
        # var = var[:nt*nlat*nlon]
        # timeFull = timeND.reshape(nt, nlon, nlat)[:, 0, 0] * timeDim / 365
        # var = var.reshape(nt, nlon, nlat).swapaxes(1, 2)

        # Index calculation
        print 'Getting %s index...' % (nino3Def[1],)
        ilatIndex = (lat >= indexDef[0][1]) & (lat <= indexDef[0][3])
        ilonIndex = (lon >= indexDef[0][0]) & (lon <= indexDef[0][2])
        indexFull = var[:, ilatIndex][:, :, ilonIndex].mean(1).mean(1)
        # Save
        np.savetxt('%s/%s%s.txt'\
                   % (indicesPath, indexDef[2], fieldDef[2]), indexFull)

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
        fig.savefig('%s/%s%s.png' % (pltDir, indexDef[2], fieldDef[2]),
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
        for iavg in np.arange(nRAVG/2, nfft-nRAVG/2):
            perioRAVG[iavg] = perio[iavg-nRAVG/2:iavg+nRAVG/2 + 1].mean() / nRAVG
        
        # Plot
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(freqYear, np.log10(perioRAVG))
        #    ax.set_xscale('log')
        #    ax.set_yscale('log')
        ax.set_xlim(0, 4)
        #ax.set_ylim(0, vmax)
        plt.title('Periodogram of %s averaged over %s' % (fieldDef[1], indexDef[1]))
        fig.savefig('%s/%s%sPerio.png' % (pltDir, indexDef[2], fieldDef[2]),
                    bbox_inches='tight')
