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
plotForward = False
#plotForward = True
#plotBackward = False
plotBackward = True
plotImag = False
#plotImag = True
plotPolar = True
#plotPolar = False
xmineigVal = -cfg.stat.rateMax
ymineigVal = -cfg.stat.angFreqMax
xlimEig = [xmineigVal, -xmineigVal/100]
ylimEig = [ymineigVal, -ymineigVal]
zlimEig = [cfg.stat.powerMin, cfg.stat.powerMax]
xticks = None
yticksPos = np.arange(0, ylimEig[1], 5.)
yticksNeg = np.arange(0, ylimEig[0], -5.)[::-1]
yticks = np.concatenate((yticksNeg, yticksPos))
zticks = np.logspace(np.log10(zlimEig[0]), np.log10(zlimEig[1]),
                    int(np.round(np.log10(zlimEig[1]/zlimEig[0]) + 1)))
zticks = np.logspace(np.log10(zlimEig[0]), np.log10(zlimEig[1]),
                     int(np.round(np.log10(zlimEig[1]/zlimEig[0])/2 + 1)))

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


compNames = ('%s%s' % (indicesName[0][0], fieldsDef[0][1]),
             '%s%s' % (indicesName[1][0], fieldsDef[1][1]))
evlabels = ('%s (%s)' % (compNames[0], fieldsDef[0][3]),
            '%s (%s)' % (compNames[1], fieldsDef[1][3]))

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

# Plot eigenvectors of transfer operator
alpha = 0.05
os.system('mkdir %s/spectrum/eigvec 2> /dev/null' % cfg.general.plotDir)

coordCut = 1
coordPlot = (coordCut + 1) % dimObs
cutVal = (coord[coordCut][mask < N] * statDist).sum()
idCut0 = np.argmin(np.abs(coord[coordCut] - cutVal))
idCut = np.abs(coord[1] - coord[1][idCut0]) < 1.e-6
xx = coord[coordPlot][idCut]

# Figure for the phase
figPhase = plt.figure()
axPhase = figPhase.add_subplot(111)
dstFilePhase = '%s/eigvec/eigvecBackwardPhaseCut_nev%d%s.%s' \
               % (specDir, nev, dstPostfixTau, ergoPlot.figFormat)
# Figure for the amp
figAmp = plt.figure()
axAmp = figAmp.add_subplot(111)
dstFileAmp = '%s/eigvec/eigvecBackwardAmpCut_nev%d%s.%s' \
             % (specDir, nev, dstPostfixTau, ergoPlot.figFormat)

# Plot
idPlot = np.array([1, 3, 5])
locPhase = 'upper right'
locAmp = 'upper right'
for ev in idPlot:
    v = np.zeros((N,), complex)
    v[mask < N] = eigVecBackward[ev]
    amp = np.abs(v[idCut])
    phase = np.angle(v[idCut])
    # Plot phase
    axPhase.plot(xx, phase, label=r'$\mathrm{arg} \psi_%d$' % ev)
    # Plot amplitude
    axAmp.plot(xx, amp, label=r'$|\psi_%d$|' % ev)
    
# Configure plot of phase and save
axPhase.legend(loc=locPhase)
axPhase.set_xlabel(evlabels[coordPlot], fontsize=ergoPlot.fs_latex)
axPhase.set_ylabel(r'$\mathrm{arg} \psi$', fontsize=ergoPlot.fs_latex)
plt.setp(axPhase.get_xticklabels(), fontsize=ergoPlot.fs_xticklabels)
plt.setp(axPhase.get_yticklabels(), fontsize=ergoPlot.fs_yticklabels)
figPhase.savefig(dstFilePhase, bbox_inches=ergoPlot.bbox_inches,
                 dpi=ergoPlot.dpi)
# Configure plot of amplitude and save
axAmp.legend(loc=locAmp)
axAmp.set_xlabel(evlabels[coordPlot], fontsize=ergoPlot.fs_latex)
axAmp.set_ylabel(r'$|\psi|$', fontsize=ergoPlot.fs_latex)
plt.setp(axAmp.get_xticklabels(), fontsize=ergoPlot.fs_xticklabels)
plt.setp(axAmp.get_yticklabels(), fontsize=ergoPlot.fs_yticklabels)
figAmp.savefig(dstFileAmp, bbox_inches=ergoPlot.bbox_inches,
               dpi=ergoPlot.dpi)

    
