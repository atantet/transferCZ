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
H = 200
tau_0 = 1.0922666667e-2
delta_T = 1.
sampFreq = 0.35 / 0.06 * 12 # (in year^-1)

# Case definition
muRng = np.array([2.1, 2.5, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95,
                  3., 3.1, 3.3, 3.7])
#amu0Rng = np.arange(0.1, 0.51, 0.1)
#amu0Rng = np.array([0.2])
amu0Rng = np.array([0.75, 1.])
#epsRng = np.arange(0., 0.11, 0.05)
epsRng = np.array([0.])
sNoise = 'without noise'

# spinup = int(sampFreq * 10)
# timeWin = np.array([0, -1])
spinup = 0
timeWin = np.array([0, int(sampFreq * 1000)])

neof = 1

prefix = 'zc'
simType = '_%deof_seasonal' % (neof,)
indicesDir = '../observables/'

field_h = (1, 'Thermocline depth', r'$h$', 'm', H, r'$h^2$')
field_T = (2, 'SST', 'T', r'$^\circ C$', delta_T, r'$(^\circ C)^2$')
#field_u_A = (3, 'Wind stress due to coupling', r'$u_A$', 'm/s', tau_0)
#field_taux = (4, 'External wind-stress', 'taux', 'm/s', tau_0)

# Indices definition
# Nino3
#nino3Def = ('NINO3', 'nino3')
nino3Def = ('Eastern', 'nino3')
# Nino4
#nino4Def = ('NINO4', 'nino4')
nino4Def = ('Western', 'nino4')

# Field and index choice
fieldDef = field_T
indexDef = nino3Def
#fieldDef = field_h
#indexDef = nino4Def
#fieldDef = field_u_A
#indexDef = nino4Def
#fieldDef = field_taux
#indexDef = nino4Def

fs_default = 'x-large'
fs_latex = 'xx-large'
fs_xlabel = fs_default
fs_ylabel = fs_default
fs_xticklabels = fs_default
fs_yticklabels = fs_default
fs_legend_title = fs_default
fs_legend_labels = fs_default
fs_cbar_label = fs_default
#            figFormat = 'eps'
figFormat = 'png'
dpi = 300

for eps in epsRng:
    for amu0 in amu0Rng:
        for mu in muRng:
            postfix = '_mu%04d_amu0%04d_eps%04d' \
                      % (np.round(mu * 1000, 1), np.round(amu0 * 1000, 1),
                         np.round(eps * 1000, 1))
            resDir = '%s%s%s/' % (prefix, simType, postfix)
            indicesPath = '%s/%s' % (indicesDir, resDir)
            pltDir = resDir
            os.system('mkdir %s %s 2> /dev/null' % (pltDir, indicesPath))

            # Read dataset
            indexData = np.loadtxt('%s/%s.txt' % (indicesPath, indexDef[1]))
            timeND = indexData[:, 0]
            timeFull = timeND * timeDim / 365
            indexFull = indexData[:, fieldDef[0]] * fieldDef[4]

            # Remove spinup
            index = indexFull[spinup:]
            time = timeFull[spinup:]
            nt = time.shape[0]
    
            # Plot time-series
            linewidth = 2
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(time[timeWin[0]:timeWin[1]], index[timeWin[0]:timeWin[1]],
                    linewidth=linewidth)
            #        plt.title('Time-series of %s averaged over %s\nfor mu = %.4f and eps = %.4f' % (fieldDef[1], indexDef[0], mu, eps))
            ax.set_xlabel('years', fontsize=fs_latex)
            ax.set_ylabel('%s %s (%s)' \
                          % (indexDef[0], fieldDef[1], fieldDef[3]),
                          fontsize=fs_latex)
#            ax.set_xlim(0, 250)
#            ax.set_ylim(29.70, 30.)
            plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
            plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
            fig.savefig('%s/%s%s%s.png' % (pltDir, indexDef[1], fieldDef[2],
                                           postfix), bbox_inches='tight')
        
            
            # Get periodogram of zonal wind stress averaged over index
            nRAVG = 1
            window = np.hamming(nt)
            # Get nearest larger power of 2
            if np.log2(nt) != int(np.log2(nt)):
                nfft = 2**(int(np.log2(nt)) + 1)
            else:
                nfft = nt
            # Get frequencies and shift zero frequency to center
            freq = np.fft.fftfreq(nfft, d=1./sampFreq)
            freq = np.fft.fftshift(freq)

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
                perioRAVG[iavg] = perio[iavg-nRAVG/2:iavg+nRAVG/2 + 1].mean()\
                              / nRAVG
        
            # Plot
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.plot(freq, np.log10(perioRAVG), '-k')
            #    ax.set_xscale('log')
            #    ax.set_yscale('log')
            ax.set_xlim(0, 4)
            #ax.set_ylim(0, vmax)
            ax.set_xlabel(r'years$^{-1}$', fontsize=fs_latex)
            ax.set_ylabel(fieldDef[5], fontsize=fs_latex)
            plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
            plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            ax.text(xlim[0] + 0.2 * (xlim[1] - xlim[0]),
                    ylim[0] + 0.1 * (ylim[1] - ylim[0]),
                    r'$\mu = %.3f$ %s' % (mu, sNoise), fontsize=32)
            #        plt.title('Periodogram of %s averaged over %s\nfor mu = %.1f and eps = %.2f' % (fieldDef[1], indexDef[0], mu, eps))
            fig.savefig('%s/%s%sPerio%s.%s' \
                        % (pltDir, indexDef[1], fieldDef[2],
                           postfix, figFormat),
                        bbox_inches='tight', dpi=dpi)
