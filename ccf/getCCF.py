import os
import numpy as np
import matplotlib.pyplot as plt
import atmath, atplot

# Define the observable
neof = 1
prefix = 'zc'
simType = '_%deof' % (neof,)
indicesDir = '../data/observables/'
resDir = '../results/'
plotDir = '../results/plot/'

# Case definition
timeFreq = 0.35 / 0.060 # (integration frequency in months
timeScaleConversion = 1. / 12
spinupYear = 100
sampFreq = timeFreq / timeScaleConversion
lagMax = 100
#lagMax = 50
tapeWindow = 100 * sampFreq
yminPerio = 1.e-8
ymaxPerio = 1.

muRng = np.array([2.5, 2.8, 2.9, 3.5])
#muRng = np.array([3.5])
epsRng = np.array([0.05])
seedRng = np.arange(10)

nSeeds = seedRng.shape[0]
spinup = int(spinupYear * sampFreq)

# Dimensions
L = 1.5e7
c0 = 2
timeDim = L / c0 / (60. * 60. * 24)
H = 200
tau_0 = 1.0922666667e-2
delta_T = 1.

# Variables definition
field_h = (1, 'Thermocline depth', 'h', 'm', H)
field_T = (2, 'SST', 'T', r'$^\circ C$', delta_T)
field_u_A = (3, 'Wind stress due to coupling', 'u_A', 'm/s', tau_0)
field_taux = (4, 'External wind-stress', 'taux', 'm/s', tau_0)

# Indices definition
# Nino3
#nino3Def = ('NINO3', 'nino3')
nino3Def = ('Eastern', 'nino3')
# Nino4
#nino4Def = ('NINO4', 'nino4')
nino4Def = ('Western', 'nino4')

#  Field and index choice
fieldDef1 = field_T
indexDef1 = nino3Def
fieldDef2 = fieldDef1
indexDef2 = indexDef1

lagMaxSample = int(np.round(lagMax * sampFreq))
lags = np.arange(0, lagMax + 0.999 / sampFreq, 1. / sampFreq)
for eps in epsRng:
    print 'eps = %f' % eps
    for mu in muRng:
        print 'mu = %f' % mu
        outDir = '%s%s_mu%04d_eps%04d/' \
                 % (prefix, simType,
                    np.round(mu * 1000, 1), np.round(eps * 1000, 1))
        # Create directories
        os.system('mkdir %s 2> /dev/null' % outDir)
        postfix = '_%s_%s_%s_%s_mu%04d_eps%04d' \
                  % (fieldDef1[2], indexDef1[1],
                     fieldDef2[2], indexDef2[1],
                     np.round(mu * 1000, 1), np.round(eps * 1000, 1))
        ccf = np.zeros((lags.shape[0],))
        for s in np.arange(nSeeds):
            seed = seedRng[s]
            print 'for seed', seed
            caseDir = '%s%s_mu%04d_eps%04d_seed%d/' \
                      % (prefix, simType, np.round(mu * 1000, 1),
                         np.round(eps * 1000, 1), seed)
            indicesPath = '%s/%s/' % (indicesDir, caseDir)

            indexPath1 = '%s/%s.txt' % (indicesPath, indexDef1[1])
            indexPath2 = '%s/%s.txt' % (indicesPath, indexDef2[1])

            # Read datasets
            indexData1 = np.loadtxt(indexPath1)
            timeFull1 = indexData1[spinup:, 0]
            observable1 = indexData1[spinup:, fieldDef1[0]]*fieldDef1[4]

            indexData2 = np.loadtxt(indexPath2)
            timeFull2 = indexData2[spinup:, 0]
            observable2 = indexData2[spinup:, fieldDef2[0]]*fieldDef2[4]

            nt = np.min([timeFull1.shape[0], timeFull2.shape[0]])
            time = timeFull1[:nt]
            observable1 = observable1[:nt]
            observable2 = observable2[:nt]

            # Get ccf averaged over seeds (should add weights based on length)
            ccf += atmath.ccf(observable1, observable2, lagMax=lagMax,
                              sampFreq=sampFreq)[lagMaxSample:]

            # Get common frequencies
            if s == 0:
                nTapes = int(nt / tapeWindow)
                freq = atmath.getFreqPow2(tapeWindow, sampFreq=sampFreq)
                nfft = freq.shape[0]
                perio = np.zeros((nfft,))
                perioSTD = np.zeros((nfft,))
                
            # Get perio averaged over seeds (should add weights based on length)
            (freq, perioSeed, perioSTDSeed) \
                = atmath.getPerio(observable1, freq=freq, sampFreq=sampFreq,
                                  tapeWindow=tapeWindow)
            perio += perioSeed
            perioSTD += perioSTDSeed**2 * nTapes
            
        ccf /= nSeeds
        perio /= nSeeds
        perioSTD = np.sqrt(perioSTD / (nSeeds * nTapes))

        # Save results
        np.savetxt('%s/%s/ccf%s_nSeeds%d_lagMax%dyr.txt'\
                   % (resDir, outDir, postfix, nSeeds, lagMax), ccf)
        np.savetxt('%s/%s/perio%s_nSeeds%d.txt' % (resDir, outDir, postfix, nSeeds), perio)
        np.savetxt('%s/%s/perioSTD%s_nSeeds%d.txt' % (resDir, outDir, postfix, nSeeds), perioSTD)
        np.savetxt('%s/%s/freq%s_nSeeds%d.txt' % (resDir, outDir, postfix, nSeeds), freq)
        
        # Plot CCF
        print 'Plotting correlation function...'
        (fig, ax) = atplot.plotCCF(ccf, lags)
        plt.savefig('%s/%s/ccf%s_nSeeds%d_lagMax%dyr.%s'\
                    % (plotDir, outDir, postfix, nSeeds, lagMax, atplot.figFormat),
                    dpi=atplot.dpi, bbox_inches=atplot.bbox_inches)

        # Plot perio
        print 'Plotting periodogram...'
        angFreq = freq * 2 * np.pi
        (fig, ax) = atplot.plotPerio(perio, perioSTD=perioSTD, freq=angFreq,
                                     absUnit='yr$^{-1}$', yscale='log',
                                     xlim=(0, 4*2*np.pi),
                                     ylim=(yminPerio, ymaxPerio))
        fig.savefig('%s/%s/perio%s_nSeeds%d.%s'\
                    % (plotDir, outDir, postfix, nSeeds, atplot.figFormat),
                    dpi=atplot.dpi, bbox_inches=atplot.bbox_inches)

