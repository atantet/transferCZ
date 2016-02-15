import os
import numpy as np
import matplotlib.pyplot as plt
import pylibconfig2
import ergoPlot, ergoStat

configFile = '../cfg/transferCZ.cfg'
cfg = pylibconfig2.Config()
cfg.read_file(configFile)

timeFreq = 0.35 / 0.060 # integration frequency in months
timeScaleConversion = 1. / 12
sampFreq = timeFreq / timeScaleConversion
dim = len(cfg.caseDef.indicesName)

field_h = (1, 'Thermocline depth', 'h', 'm', cfg.dimension.H)
field_T = (2, 'SST', 'T', r'$^\circ C$', cfg.dimension.delta_T)
field_u_A = (3, 'Wind stress due to coupling', 'u_A', 'm/s', cfg.dimension.tau_0)
field_taux = (4, 'External wind-stress', 'taux', 'm/s', cfg.dimension.tau_0)

nino3 = ('Eastern', 'nino3')
nino4 = ('Western', 'nino4')
indicesName = []
fieldsDef = []
for d in np.arange(dim):
    if cfg.stat.indicesName[d] == 'nino3':
        indicesName.append(nino3)
    elif cfg.stat.indicesName[d] == 'nino4':
        indicesName.append(nino4)
    if cfg.stat.fieldsName[d] == 'T':
        fieldsDef.append(field_T)
    if cfg.stat.fieldsName[d] == 'h':
        fieldsDef.append(field_h)

muRng = np.array([2.5, 2.8, 2.9, 3.5])
#muRng = np.array([2.9])
epsRng = np.array([0.05])
spinup = int(cfg.simulation.spinupMonth / 12 * sampFreq)
nSeeds = len(cfg.caseDef.seedRng)


lagMaxSample = int(np.round(cfg.stat.lagMax * sampFreq))
lags = np.arange(-cfg.stat.lagMax, cfg.stat.lagMax + 0.999 / sampFreq, 1. / sampFreq)
for eps in epsRng:
    print 'eps = %f' % eps
    for mu in muRng:
        print 'mu = %f' % mu
        outDir = '%s%s_mu%04d_eps%04d' \
                 % (cfg.caseDef.prefix, cfg.caseDef.simType,
                    np.round(mu * 1000, 1), np.round(eps * 1000, 1))
        # Create directories
        postfix = '_%s_%s_%s_%s_mu%04d_eps%04d' \
                  % (fieldsDef[0][2], indicesName[0][1],
                     fieldsDef[1][2], indicesName[1][1],
                     np.round(mu * 1000, 1), np.round(eps * 1000, 1))
        corrSample = np.zeros((lags.shape[0],))
        for s in np.arange(nSeeds):
            seed = cfg.caseDef.seedRng[s]
            print 'for seed', seed
            caseDir = '%s_seed%d/' \
                      % (outDir, seed)
            indicesPath = '%s/%s/' % (cfg.general.indicesDir, caseDir)

            indexPath1 = '%s/%s.txt' % (indicesPath, indicesName[0][1])
            indexPath2 = '%s/%s.txt' % (indicesPath, indicesName[1][1])

            # Read datasets
            indexData1 = np.loadtxt(indexPath1)
            timeFull1 = indexData1[spinup:, 0]
            observable1 = indexData1[spinup:, fieldsDef[0][0]]*fieldsDef[0][4]

            indexData2 = np.loadtxt(indexPath2)
            timeFull2 = indexData2[spinup:, 0]
            observable2 = indexData2[spinup:, fieldsDef[1][0]]*fieldsDef[1][4]

            nt = np.min([timeFull1.shape[0], timeFull2.shape[0]])
            time = timeFull1[:nt]
            observable1 = observable1[:nt]
            observable2 = observable2[:nt]

            # Get corrSample averaged over seeds (should add weights based on length)
            corrSample += ergoStat.ccf(observable1, observable2, lagMax=cfg.stat.lagMax,
                                sampFreq=sampFreq)

            # Get common frequencies
            if s == 0:
                nChunks = int(nt / (cfg.stat.chunkWidth * sampFreq))
                freq = ergoStat.getFreqPow2(cfg.stat.chunkWidth, sampFreq=sampFreq)
                nfft = freq.shape[0]
                powerSample = np.zeros((nfft,))
                powerSampleSTD = np.zeros((nfft,))
                
            # Get powerSample averaged over seeds (should add weights based on length)
            (freq, powerSampleSeed, powerSampleSTDSeed) \
                = ergoStat.getPerio(observable1, observable2, freq=freq, sampFreq=sampFreq,
                                    chunkWidth=cfg.stat.chunkWidth, norm=False)
            powerSample += powerSampleSeed
            powerSampleSTD += powerSampleSTDSeed**2 * nChunks
            
        corrSample /= nSeeds
        powerSample /= nSeeds
        powerSampleSTD = np.sqrt(powerSampleSTD / (nSeeds * nChunks))

        # Save results
        np.savetxt('%s/%s/corrSample%s_nSeeds%d_lagMax%dyr.txt'\
                   % (cfg.general.resDir, outDir, postfix, nSeeds, cfg.stat.lagMax), corrSample)
        np.savetxt('%s/%s/lags%s_nSeeds%d_lagMax%dyr.txt' \
                   % (cfg.general.resDir, outDir, postfix, nSeeds, cfg.stat.lagMax), lags)
        np.savetxt('%s/%s/powerSample%s_nSeeds%d_chunk%dyr.txt' \
                   % (cfg.general.resDir, outDir, postfix, nSeeds, cfg.stat.chunkWidth), powerSample)
        np.savetxt('%s/%s/powerSampleSTD%s_nSeeds%d_chunk%dyr.txt' \
                   % (cfg.general.resDir, outDir, postfix, nSeeds, cfg.stat.chunkWidth), powerSampleSTD)
        np.savetxt('%s/%s/freq%s_nSeeds%d_chunk%dyr.txt' \
                   % (cfg.general.resDir, outDir, postfix, nSeeds, cfg.stat.chunkWidth), freq)
        
        # Plot CORRSAMPLE
        print 'Plotting correlation function...'
        (fig, ax) = ergoPlot.plotCCF(corrSample, lags, plotPositive=True)
        plt.savefig('%s/%s/corrSample%s_nSeeds%d_lagMax%dyr.%s'\
                    % (cfg.general.plotDir, outDir, postfix, nSeeds, cfg.stat.lagMax, ergoPlot.figFormat),
                    dpi=ergoPlot.dpi, bbox_inches=ergoPlot.bbox_inches)

        # Plot powerSample
        print 'Plotting powerSampledogram...'
        angFreq = freq * 2 * np.pi
        (fig, ax) = ergoPlot.plotPerio(powerSample, perioSTD=powerSampleSTD, freq=angFreq, 
                                       plotPositive=True,
                                       absUnit='', yscale='log',
                                       xlim=(0, cfg.stat.angFreqMax),
                                       ylim=(cfg.stat.yminPower, cfg.stat.ymaxPower))
        fig.savefig('%s/%s/powerSample%s_nSeeds%d_chunk%dyr.%s'\
                    % (cfg.general.plotDir, outDir, postfix, nSeeds, cfg.stat.chunkWidth, ergoPlot.figFormat),
                    dpi=ergoPlot.dpi, bbox_inches=ergoPlot.bbox_inches)

