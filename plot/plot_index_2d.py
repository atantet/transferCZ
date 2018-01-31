import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import pylibconfig2
import ergoPlot

configFile = '../cfg/transferCZ.cfg'
cfg = pylibconfig2.Config()
cfg.read_file(configFile)

timeFreq = 0.35 / 0.06 # integration frequency in months
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

# Dimensions
timeDim = cfg.units.L / cfg.units.c0 / (60. * 60. * 24)

# Case definition
#muRng = np.arange(2.1, 4., 0.2)
#muRng = np.array([2.5, 2.9, 3.5])
muRng = np.array([2.8])
#epsRng = np.arange(0., 0.41, 0.05)
#epsRng = np.array([0.])
#sNoise = 'without noise'
epsRng = np.array([0.05])
sNoise = 'with noise'

#spinup = int(cfg.simulation.spinupMonth / 12 * sampFreq)
timeWin = np.array([0, int(sampFreq * 50)])
nSeeds = len(cfg.caseDef.seedRng)

for eps in epsRng:
    for mu in muRng:
        srcPostfix = '_mu%04d_eps%04d' % (mu * 1000, eps * 1000)
        caseDirPrefix = '%s%s%s' % (cfg.caseDef.prefix, cfg.caseDef.simType,
                                    srcPostfix)
        pltDir = '%s/' % caseDirPrefix
        os.system('mkdir %s 2> /dev/null' % (pltDir,))
        print 'Treating case (eps, mu) = (', eps, ',', mu, ')'

        indexFull1 = np.array([])
        timeFull = np.array([])
        indexFull2 = np.array([])
        for s in np.arange(nSeeds):
            seed = cfg.caseDef.seedRng[s]
            obsPath = '%s/%s_seed%d' % (cfg.general.indicesDir,
                                        caseDirPrefix, seed)
            print 'Loading data for seed ', seed

            # Read first component
            data = np.loadtxt("%s/%s.txt" % (obsPath, indicesName[0][1]))
            timeND = data[:, 0]
            timeSeed = timeND * timeDim / 365
            if seed > 0:
                timeSeed += timeFull[-1]
            timeFull = np.concatenate((timeFull, timeSeed))
            indexSeed = data[:, fieldsDef[0][0]] * fieldsDef[0][4]
            indexFull1 = np.concatenate((indexFull1, indexSeed))

            # Read second component
            data = np.loadtxt("%s/%s.txt" % (obsPath, indicesName[1][1]))
            indexSeed = data[:, fieldsDef[1][0]] * fieldsDef[1][4]
            indexFull2 = np.concatenate((indexFull2, indexSeed))
            
        nt = timeFull.shape[0]
    
        # Plot time-series
        linewidth = 2
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(indexFull1[timeWin[0]:timeWin[1]],
                indexFull2[timeWin[0]:timeWin[1]],
                linewidth=linewidth)
#        plt.title('Time-series of %s averaged over %s\nfor mu = %.1f and eps = %.2f' % (fieldsDef[1], indexDef[0], mu, eps))
        ax.set_xlabel('%s %s (%s)' % (indicesName[0][0], fieldsDef[0][1],
                                      fieldsDef[0][3]),
                      fontsize=ergoPlot.fs_latex)
        ax.set_ylabel('%s %s (%s)' % (indicesName[1][0], fieldsDef[1][1],
                                      fieldsDef[1][3]),
                      fontsize=ergoPlot.fs_latex)
        plt.setp(ax.get_xticklabels(), fontsize=ergoPlot.fs_xticklabels)
        plt.setp(ax.get_yticklabels(), fontsize=ergoPlot.fs_yticklabels)
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.text(xlim[0] + 0.2 * (xlim[1] - xlim[0]),
                ylim[0] + 0.9 * (ylim[1] - ylim[0]),
                r'$\mu = %.3f$ %s' % (mu, sNoise), fontsize=32)
        figName = '%s/%s%s%s%s%s.%s' \
                    % (pltDir, indicesName[0][1], fieldsDef[0][2],
                       indicesName[1][1], fieldsDef[1][2],
                       srcPostfix, figFormat)
        fig.savefig(figName, bbox_inches=ergoPlot.bbox_inches,
                    dpi=ergoPlot.dpi)
    
