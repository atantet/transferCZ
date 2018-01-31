import os
import numpy as np
import matplotlib.pyplot as plt
import pylibconfig2
import ergoPlot

configFile = '../cfg/transferCZ.cfg'
cfg = pylibconfig2.Config()
cfg.read_file(configFile)
fileFormat = cfg.general.fileFormat

# Transition lag
if (hasattr(cfg.stat, 'tauDimPlot')):
    tauDim = cfg.stat.tauDimPlot
else:
    tauDim = cfg.transfer.tauRng[0]

timeScaleConversion = 1. / 12
dimObs = len(cfg.caseDef.indicesName)
nSeeds = len(cfg.caseDef.seedRng)

nev = cfg.spectrum.nev
evPlot = np.array([1, 3, 5])
plotForward = False
#plotForward = True
#plotBackward = False
plotBackward = True
ampMin = 0.
ampMax = 0.07
nlevAmp = int((ampMax - ampMin) * 100 + 0.1) + 1
csfilter = 0.5
csfmt = '%1.2f'

field_h = (1, 'H', 'h', 'm')
field_T = (2, 'SST', 'T', r'$^\circ C$')
field_u_A = (3, 'wind stress due to coupling', 'u_A', 'm/s')
field_taux = (4, 'external wind-stress', 'taux', 'm/s')

nino3 = ('E', 'nino3')
nino4 = ('W', 'nino4')
indicesName = []
fieldsDef = []
for d in np.arange(dimObs):
    if cfg.caseDef.indicesName[d] == 'nino3':
        indicesName.append(nino3)
    elif cfg.caseDef.indicesName[d] == 'nino4':
        indicesName.append(nino4)
    if cfg.caseDef.fieldsName[d] == 'T':
        fieldsDef.append(field_T)
    if cfg.caseDef.fieldsName[d] == 'h':
        fieldsDef.append(field_h)


compName0 = '%s%s' % (indicesName[0][0], fieldsDef[0][1])
ev_xlabel = '%s (%s)' % (compName0, fieldsDef[0][3])
compName1 = '%s%s' % (indicesName[1][0], fieldsDef[1][1])
ev_ylabel = '%s (%s)' % (compName1, fieldsDef[1][3])

srcPostfix = "%s%s_mu%04d_eps%04d" % (cfg.caseDef.prefix, cfg.caseDef.simType,
                                      np.round(cfg.caseDef.mu * 1000, 1),
                                      np.round(cfg.caseDef.eps * 1000, 1))
obsName = ''
gridPostfix = ''
N = 1
for d in np.arange(dimObs):
    obsName += '_%s_%s' % (fieldsDef[d][2], indicesName[d][1])
    N *= cfg.grid.nx[d]
    gridPostfix = "%s_n%dl%dh%d" % (gridPostfix, cfg.grid.nx[d],
                                    cfg.grid.nSTDLow[d], cfg.grid.nSTDHigh[d])
cpyBuffer = gridPostfix
gridPostfix = '_%s%s%s' % (srcPostfix, obsName, cpyBuffer)

# Read grid
gridFile = '%s/grid/grid%s.txt' % (cfg.general.resDir, gridPostfix)
coord = ergoPlot.readGrid(gridFile, dimObs)

# Coordinate matrices read in 'ij' indexing (not 'xy')!
if dimObs == 1:
    X = coord[0]
elif dimObs == 2:
    X, Y = np.meshgrid(coord[0], coord[1], indexing='ij')
    coord = (X.flatten(), Y.flatten())
elif dimObs == 3:
    X, Y, Z = np.meshgrid(coord[0], coord[1], coord[2], indexing='ij')
    coord = (X.flatten(), Y.flatten(), Z.flatten())

tau = tauDim * timeScaleConversion
dstPostfix = gridPostfix
dstPostfixTau = "%s_tau%03d" % (gridPostfix, int(tauDim * 1000 + 0.1))
specDir = '%s/spectrum/' % cfg.general.plotDir

# File names
eigValForwardFile = '%s/eigval/eigValForward_nev%d%s.%s' \
                    % (cfg.general.specDir, nev, dstPostfixTau, fileFormat)
eigVecForwardFile = '%s/eigvec/eigVecForward_nev%d%s.%s' \
                    % (cfg.general.specDir, nev, dstPostfixTau, fileFormat)
eigValBackwardFile = '%s/eigval/eigValBackward_nev%d%s.%s' \
                    % (cfg.general.specDir, nev, dstPostfixTau, fileFormat)
eigVecBackwardFile = '%s/eigvec/eigVecBackward_nev%d%s.%s' \
                    % (cfg.general.specDir, nev, dstPostfixTau, fileFormat)
statDistFile = '%s/transfer/initDist/initDist%s.%s' \
               % (cfg.general.resDir, dstPostfix, fileFormat)
maskFile = '%s/transfer/mask/mask%s.%s' \
           % (cfg.general.resDir, dstPostfix, fileFormat)

# Read stationary distribution
if statDistFile is not None:
    if fileFormat == 'bin':
        statDist = np.fromfile(statDistFile, float)
    else:
        statDist = np.loadtxt(statDistFile, float)
else:
    statDist = None

# Read mask
if maskFile is not None:
    if fileFormat == 'bin':
        mask = np.fromfile(maskFile, np.int32)
    else:
        mask = np.loadtxt(maskFile, np.int32)
else:
    mask = np.arange(N)
NFilled = np.max(mask[mask < N]) + 1

# Read transfer operator spectrum from file and create a bi-orthonormal basis
# of eigenvectors and backward eigenvectors:
print 'Readig spectrum for tauDim = %.3f...' % tauDim
(eigValForward, eigValBackward, eigVecForward, eigVecBackward) \
    = ergoPlot.readSpectrum(eigValForwardFile, eigValBackwardFile,
                            eigVecForwardFile, eigVecBackwardFile,
                            makeBiorthonormal=~cfg.spectrum.makeBiorthonormal,
                            fileFormat=fileFormat) 

print 'Getting conditionning of eigenvectors...'
eigenCondition = ergoPlot.getEigenCondition(eigVecForward, eigVecBackward)

# Get generator eigenvalues
eigValGen = (np.log(np.abs(eigValForward)) + np.angle(eigValForward)*1j) / tau


# Plot eigenvectors of transfer operator
alpha = 0.05
os.system('mkdir %s/spectrum/eigvec 2> /dev/null' % cfg.general.plotDir)
os.system('mkdir %s/spectrum/reconstruction 2> /dev/null' \
          % cfg.general.plotDir)
for ev in evPlot:
    if ev == 0:
        positive = True
        cmap = cm.hot_r
    else:
        positive = False
        cmap = cm.RdBu_r
        
    if plotForward:
        print 'Plotting polar eigenvector %d...' % (ev + 1,)
        fig, = ergoPlot.plotEigVecPolar(X, Y, eigVecForward[ev],
                                        mask=mask, combine=True,
                                        xlabel=ev_xlabel, ylabel=ev_ylabel,
                                        alpha=alpha, cmap=cmap, csfmt=csfmt,
                                        ampMin=ampMin, ampMax=ampMax,
                                        nlevAmp=nlevAmp, csfilter=csfilter)
        dstFile = '%s/eigvec/eigvecForwardPolar_nev%d_ev%03d%s.%s' \
                  % (specDir, nev, ev + 1, dstPostfixTau, ergoPlot.figFormat)
        fig.savefig(dstFile, bbox_inches=ergoPlot.bbox_inches,
                    dpi=ergoPlot.dpi)
    if plotBackward:
        print 'Plotting polar backward eigenvector %d...' % (ev + 1,)
        fig, = ergoPlot.plotEigVecPolar(X, Y, eigVecBackward[ev],
                                        mask=mask, combine=True,
                                        xlabel=ev_xlabel, ylabel=ev_ylabel,
                                        alpha=alpha, cmap=cmap, csfmt=csfmt,
                                        ampMin=ampMin, ampMax=ampMax,
                                        nlevAmp=nlevAmp, csfilter=csfilter)
        dstFile = '%s/eigvec/eigvecBackwardPolar_nev%d_ev%03d%s.%s' \
                  % (specDir, nev, ev + 1, dstPostfixTau, ergoPlot.figFormat)
        fig.savefig(dstFile, bbox_inches=ergoPlot.bbox_inches,
                    dpi=ergoPlot.dpi)
