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
sampFreq = timeFreq / timeScaleConversion # sampling frequency in years
dim = len(cfg.caseDef.indicesName)

field_h = (1, 'Thermocline depth', 'h', 'm', cfg.dimension.H)
field_T = (2, 'SST', 'T', r'$^\circ C$', cfg.dimension.delta_T)
field_u_A = (3, 'Wind stress due to coupling', 'u_A', 'm/s',
             cfg.dimension.tau_0)
field_taux = (4, 'External wind-stress', 'taux', 'm/s',
              cfg.dimension.tau_0)

nino3 = ('Eastern', 'nino3')
nino4 = ('Western', 'nino4')
indicesName = []
fieldsDef = []
for d in np.arange(dim):
    if cfg.stat.indicesName[d] == 'nino3':
        indicesName.append(nino3)
    elif cfg.stat.indicesName[d] == 'nino4':
        indicesName.append(nino4)
    else:
        raise ValueError('Index name %d not recognized' % d)

    if cfg.stat.fieldsName[d] == 'T':
        fieldsDef.append(field_T)
    elif cfg.stat.fieldsName[d] == 'h':
        fieldsDef.append(field_h)
    else:
        raise ValueError('Field name %d not recognized' % d)

spinup = int(cfg.simulation.spinupMonth / 12 * sampFreq)
nSeeds = len(cfg.caseDef.seedRng)


lagMaxSample = int(cfg.stat.lagMax * sampFreq + 0.1)
lags = np.arange(-cfg.stat.lagMax, cfg.stat.lagMax + 0.999 / sampFreq,
                 1. / sampFreq)
nLags = lags.shape[0]

eps = cfg.caseDef.eps
print 'eps = %f' % eps
mu = cfg.caseDef.mu
print 'mu = %f' % mu
outDir = '%s%s_mu%04d_eps%04d' \
         % (cfg.caseDef.prefix, cfg.caseDef.simType,
            int(mu * 1000 + 0.1), int(eps * 1000 + 0.1))
postfix = '_%s_%s_%s_%s_mu%04d_eps%04d' \
          % (fieldsDef[0][2], indicesName[0][1],
             fieldsDef[1][2], indicesName[1][1],
             int(mu * 1000 + 0.1), int(eps * 1000 + 0.1))
dstPostfix = "%s_nSeeds%d" % (postfix, nSeeds)
corrSample = np.zeros((nLags,))
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
    # (do not normalize here, because we summup the seeds)
    corrSample += ergoStat.ccf(observable1, observable2,
                               lagMax=cfg.stat.lagMax,
                               sampFreq=sampFreq, norm=False)

    # Get common frequencies
    if s == 0:
        nChunks = int(nt / (cfg.stat.chunkWidth * sampFreq))
        freq = ergoStat.getFreqPow2(cfg.stat.chunkWidth,
                                    sampFreq=sampFreq)
        nfft = freq.shape[0]
        powerSample = np.zeros((nfft,))
        powerSampleSTD = np.zeros((nfft,))
                
    # Get powerSample averaged over seeds
    # (should add weights based on length)
    (freq, powerSampleSeed, powerSampleSTDSeed) \
        = ergoStat.getPerio(observable1, observable2,
                            freq=freq, sampFreq=sampFreq,
                            chunkWidth=cfg.stat.chunkWidth, norm=False)
    powerSample += powerSampleSeed
    powerSampleSTD += powerSampleSTDSeed**2 * nChunks
            
corrSample /= nSeeds
powerSample /= nSeeds
powerSampleSTD = np.sqrt(powerSampleSTD / (nSeeds * nChunks))
if cfg.stat.norm:
    cov = corrSample[(lags.shape[0] - 1) / 2]
    corrSample /= cov
    powerSample /= cov
    powerSampleSTD /= cov

# Save results
np.savetxt('%s/correlation/corrSample%s_lagMax%dyr.txt'\
           % (cfg.general.resDir, dstPostfix,
              cfg.stat.lagMax), corrSample)
np.savetxt('%s/correlation/lags%s_lagMax%dyr.txt' \
           % (cfg.general.resDir, dstPostfix,
              cfg.stat.lagMax), lags)
np.savetxt('%s/power/powerSample%s_chunk%dyr.txt' \
           % (cfg.general.resDir, dstPostfix,
              cfg.stat.chunkWidth), powerSample)
np.savetxt('%s/power/powerSampleSTD%s_chunk%dyr.txt' \
           % (cfg.general.resDir, dstPostfix,
              cfg.stat.chunkWidth), powerSampleSTD)
np.savetxt('%s/power/freq%s_chunk%dyr.txt' \
           % (cfg.general.resDir, dstPostfix,
              cfg.stat.chunkWidth), freq)
        
# Plot corrSample
print 'Plotting correlation function...'
(fig, ax) = ergoPlot.plotCCF(corrSample, lags, plotPositive=True)
plt.savefig('%s/correlation/corrSample%s_lagMax%dyr.%s'\
            % (cfg.general.plotDir, dstPostfix,
               cfg.stat.lagMax, ergoPlot.figFormat),
            dpi=ergoPlot.dpi, bbox_inches=ergoPlot.bbox_inches)

# Plot powerSample
print 'Plotting powerSampledogram...'
angFreq = freq * 2 * np.pi
(fig, ax) = ergoPlot.plotPerio(powerSample, perioSTD=powerSampleSTD,
                               freq=angFreq,  plotPositive=True,
                               absUnit='', yscale='log',
                               xlim=(0, cfg.stat.angFreqMax),
                               ylim=(cfg.stat.powerMin, cfg.stat.powerMax))
fig.savefig('%s/power/powerSample%s_chunk%dyr.%s'\
            % (cfg.general.plotDir, dstPostfix,
               cfg.stat.chunkWidth, ergoPlot.figFormat),
            dpi=ergoPlot.dpi, bbox_inches=ergoPlot.bbox_inches)

