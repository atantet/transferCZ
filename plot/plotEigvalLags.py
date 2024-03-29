import os
import numpy as np
import matplotlib.pyplot as plt
import pylibconfig2
from ergopack import ergoplot
from matplotlib import cm, rcParams
import os
import pandas as pd

configFile = os.path.join('..', 'cfg', 'transferCZ.cfg')
cfg = pylibconfig2.Config()
cfg.read_file(configFile)
fileFormat = cfg.general.fileFormat

# Transition lag
timeScaleConversion = 1. / 12
dimObs = len(cfg.caseDef.indicesName)
nSeeds = len(cfg.caseDef.seedRng)

nev = cfg.spectrum.nev
nevPlot = 10
plotForward = False
#plotForward = True
#plotBackward = False
plotBackward = True
plotImag = False
#plotImag = True
plotPolar = True
colors = rcParams['axes.prop_cycle'].by_key()['color']
#plotPolar = False
ampMin = 0.
ampMax = 0.07
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
fileName = 'grid{}.txt'.format(gridPostfix)
gridFilePath = os.path.join(cfg.general.resDir, 'grid', fileName)
coord = ergoplot.readGrid(gridFilePath, dimObs)

# Coordinate matrices read in 'ij' indexing (not 'xy')!
if dimObs == 1:
    X = coord[0]
elif dimObs == 2:
    X, Y = np.meshgrid(coord[0], coord[1], indexing='ij')
    coord = (X.flatten(), Y.flatten())
elif dimObs == 3:
    X, Y, Z = np.meshgrid(coord[0], coord[1], coord[2], indexing='ij')
    coord = (X.flatten(), Y.flatten(), Z.flatten())
dstPostfix = gridPostfix

iTauDimRng = np.array(cfg.transfer.tauDimRng).astype(int)
idx = pd.Index(iTauDimRng, name='tauDim')
cols = pd.Index(range(nev), name='iev')
eigValTau = pd.DataFrame(index=idx, columns=cols, dtype=complex)
eigCondTau = pd.DataFrame(index=idx, columns=cols)
tauRng = iTauDimRng * timeScaleConversion
for ktau, tauDim in enumerate(cfg.transfer.tauDimRng):
    iTauDim = int(tauDim)
    tau = tauRng[ktau]
    dstPostfixTau = "%s_tau%03d" % (gridPostfix, int(tauDim * 1000 + 0.1))
    specDir = os.path.join(cfg.general.plotDir, 'spectrum')

    # File names
    fileName = 'eigValForward_nev{:d}{}.{}'.format(nev, dstPostfixTau, fileFormat)
    eigValForwardFile = os.path.join(cfg.general.specDir, 'eigval', fileName)
    fileName = 'eigVecForward_nev{:d}{}.{}'.format(nev, dstPostfixTau, fileFormat)
    eigVecForwardFile = os.path.join(cfg.general.specDir, 'eigvec', fileName)
    fileName = 'eigValBackward_nev{:d}{}.{}'.format(nev, dstPostfixTau, fileFormat)
    eigValBackwardFile = os.path.join(cfg.general.specDir, 'eigval', fileName)
    fileName = 'eigVecBackward_nev{:d}{}.{}'.format(nev, dstPostfixTau, fileFormat)
    eigVecBackwardFile = os.path.join(cfg.general.specDir, 'eigvec', fileName)
    fileName = 'initDist{}.{}'.format(dstPostfix, fileFormat)
    statDistFile = os.path.join(cfg.general.resDir, 'transfer', 'initDist', fileName)
    fileName = 'mask{}.{}'.format(dstPostfix, fileFormat)
    maskFile = os.path.join(cfg.general.resDir, 'transfer', 'mask', fileName)

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
    print('Readig spectrum for tauDim = {:.3f}...'.format(tauDim))
    (eigValForward, eigValBackward, eigVecForward, eigVecBackward) \
        = ergoplot.readSpectrum(eigValForwardFile, eigValBackwardFile,
                                eigVecForwardFile, eigVecBackwardFile,
                                makeBiorthonormal=~cfg.spectrum.makeBiorthonormal,
                                fileFormat=fileFormat) 

    print('Getting conditionning of eigenvectors...')
    eigenCondition = ergoplot.getEigenCondition(eigVecForward, eigVecBackward)
    eigCondTau.loc[iTauDim, :nevPlot-1] = eigenCondition[:nevPlot]
        

    # Get generator eigenvalues
    eigValTau.loc[iTauDim, :nevPlot-1] = ((np.log(np.abs(eigValForward)) + \
                                           np.angle(eigValForward)*1j) / tau)[:nevPlot]

##
fig = plt.figure()
ax = fig.add_subplot(111)
lw = 2
ls = '-'
for iev in range(nev):
    ax.plot(tauRng, eigValTau.loc[:, iev].real, linewidth=lw, linestyle=ls,
            color=colors[iev%len(colors)])
ax.set_xlabel(r'$\tau$', fontsize=ergoplot.fs_latex)
ax.set_ylabel(r'$\mathrm{Re}(\lambda_k)(\tau)$', fontsize=ergoplot.fs_latex)
ax.set_ylim(-1., 0.)

fig = plt.figure()
ax = fig.add_subplot(111)
lw = 2
ls = '-'
for iev in range(nev):
    ax.plot(tauRng, eigCondTau.loc[:, iev].real, linewidth=lw, linestyle=ls,
            color=colors[iev%len(colors)])
ax.set_xlabel(r'$\tau$', fontsize=ergoplot.fs_latex)
ax.set_ylabel(r'$\sigma(\lambda_k)(\tau)$', fontsize=ergoplot.fs_latex)
ax.set_ylim(1., 5.)

plt.show(block=False)
##
