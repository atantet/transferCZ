import os
import numpy as np
import matplotlib.pyplot as plt
import pylibconfig2
from ergopack import ergoplot

configFile = '../cfg/transferCZ.cfg'
cfg = pylibconfig2.Config()
cfg.read_file(configFile)
fileFormat = cfg.general.fileFormat
indicesName = cfg.stat.indicesName
makeBiorthonormal = ~cfg.spectrum.makeBiorthonormal

nevPlot = 0
plotForward = True
plotBackward = True

# Transition lag
if (hasattr(cfg.stat, 'tauDimPlot')):
    tauDim = cfg.stat.tauDimPlot
else:
    tauDim = cfg.transfer.tauDimRng[0]
unfold = False
tauDimUnfold = None
if hasattr(cfg.stat, 'unfold'):
    if cfg.stat.unfold:
        unfold = True
        if hasattr(cfg.stat, 'tauDimUnfold'):
            tauDimUnfold = cfg.stat.tauDimUnfold
        else:
            tauDimUnfold = tauDim + 1

timeScaleConversion = 1. / 12
dimObs = len(indicesName)

field_h = (1, 'thermocline depth', 'h', 'm')
field_T = (2, 'SST', 'T', r'$^\circ C$')
field_u_A = (3, 'wind stress due to coupling', 'u_A', 'm/s')
field_taux = (4, 'external wind-stress', 'taux', 'm/s')

nino3 = ('Eastern', 'nino3')
nino4 = ('Western', 'nino4')
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

realLabel = r'$\Re(\lambda_k)$'
imagLabel = r'$\Im(\lambda_k)$'

ampMin = None
ampMax = None
nev = cfg.spectrum.nev
# plotBackward = False
plotBackward = True
# plotImag = False
# plotImag = True
plotPolar = True
xmineigVal = -cfg.stat.rateMax
ymineigVal = -cfg.stat.angFreqMax
plotImag = False
# plotImag = True
xlimEig = [xmineigVal, -xmineigVal/100]
ylimEig = [ymineigVal, -ymineigVal]
xticks = None
yticksPos = np.arange(0, ylimEig[1], 5.)
yticksNeg = np.arange(0, ylimEig[0], -5.)[::-1]
yticks = np.concatenate((yticksNeg, yticksPos))

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
coord = ergoplot.readGrid(gridFile, dimObs)
if dimObs == 1:
    X = coord[0]
elif dimObs == 2:
    X, Y = np.meshgrid(coord[0], coord[1], indexing='ij')
    coord = (X.flatten(), Y.flatten())
elif dimObs == 3:
    X, Y, Z = np.meshgrid(coord[0], coord[1], coord[2], indexing='ij')
    coord = (X.flatten(), Y.flatten(), Z.flatten())

seedRngRg = [[0, 1, 2, 3, 4, 5, 6, 7, 8],
             [0, 1, 2, 3, 4, 5, 6, 7, 9],
             [0, 1, 2, 3, 4, 5, 6, 8, 9],
             [0, 1, 2, 3, 4, 5, 7, 8, 9],
             [0, 1, 2, 3, 4, 6, 7, 8, 9],
             [0, 1, 2, 3, 5, 6, 7, 8, 9],
             [0, 1, 2, 4, 5, 6, 7, 8, 9],
             [0, 1, 3, 4, 5, 6, 7, 8, 9],
             [0, 2, 3, 4, 5, 6, 7, 8, 9],
             [1, 2, 3, 4, 5, 6, 7, 8, 9]]

fig, ax = plt.subplots()
for kk, seedRng in enumerate(seedRngRg):
    nSeeds = len(seedRng)
    seedPostfix = '_seeds' + ''.join(str(s) for s in seedRng)
    dstPostfix = '{}{}'.format(gridPostfix, seedPostfix)
    tau = tauDim * timeScaleConversion
    dstPostfixTau = "%s_tau%03d" % (dstPostfix, int(tauDim * 1000 + 0.1))
    specDir = os.path.join(cfg.general.plotDir, 'spectrum')

    # File names
    fileName = 'eigvalForward_nev{:d}{}.{}'.format(
        nev, dstPostfixTau, fileFormat)
    eigValForwardFile = os.path.join(cfg.general.specDir, 'eigval', fileName)
    fileName = 'eigvecForward_nev{:d}{}.{}'.format(
        nev, dstPostfixTau, fileFormat)
    eigVecForwardFile = os.path.join(cfg.general.specDir, 'eigvec', fileName)
    fileName = 'eigvalBackward_nev{:d}{}.{}'.format(
        nev, dstPostfixTau, fileFormat)
    eigValBackwardFile = os.path.join(cfg.general.specDir, 'eigval', fileName)
    fileName = 'eigvecBackward_nev{:d}{}.{}'.format(
        nev, dstPostfixTau, fileFormat)
    eigVecBackwardFile = os.path.join(cfg.general.specDir, 'eigvec', fileName)
    fileName = 'initDist{}.{}'.format(dstPostfix, fileFormat)
    statDistFile = os.path.join(cfg.general.resDir, 'transfer', 'initDist',
                                fileName)
    fileName = 'mask{}.{}'.format(dstPostfix, fileFormat)
    maskFile = os.path.join(cfg.general.resDir, 'transfer', 'mask', fileName)

    if unfold:
        tauUnfold = tauDimUnfold * timeScaleConversion
        dstPostfixTauUnfold = "{}_tau{:03d}".format(
            dstPostfix, int(tauDimUnfold * 1000 + 0.1))
        fileName = 'eigvalForward_nev{:d}{}.{}'.format(
            nev, dstPostfixTauUnfold, fileFormat)
        eigValForwardFileUnfold = os.path.join(
            cfg.general.specDir, 'eigval', fileName)
        fileName = 'eigvalBackward_nev{:d}{}.{}'.format(
            nev, dstPostfixTauUnfold, fileFormat)
        eigValBackwardFileUnfold = os.path.join(
            cfg.general.specDir, 'eigval', fileName)
        fileName = 'eigvecForward_nev{:d}{}.{}'.format(
            nev, dstPostfixTauUnfold, fileFormat)
        eigVecForwardFileUnfold = os.path.join(
            cfg.general.specDir, 'eigvec', fileName)
        fileName = 'eigvecBackward_nev{:d}{}.{}'.format(
            nev, dstPostfixTauUnfold, fileFormat)
        eigVecBackwardFileUnfold = os.path.join(
            cfg.general.specDir, 'eigvec', fileName)

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

    # Read transfer operator spectrum from file and create a bi-orthonormal
    # basis of eigenvectors and backward eigenvectors:
    print('Readig spectrum for tauDim = {:.3f}...'.format(tauDim))
    (eigValForward, eigValBackward, eigVecForward, eigVecBackward) \
        = ergoplot.readSpectrum(eigValForwardFile, eigValBackwardFile,
                                eigVecForwardFile, eigVecBackwardFile,
                                makeBiorthonormal=makeBiorthonormal,
                                fileFormat=fileFormat)
    eigenCondition = ergoplot.getEigenCondition(eigVecForward, eigVecBackward,
                                                statDist)

    eigValGenOrig = ergoplot.eig2Generator(eigValForward, tau)

    # Read unfolding eigenvalues
    # eigValForwardUnfold = None
    st = statDist.copy()
    st2 = np.concatenate((st, st))
    alpha = 0.05
    if unfold:
        print('Readig spectrum for tauDimUnfold = {:.3f} to unfold...'.format(
            tauDimUnfold))
        (eigValForwardUnfold, eigValBackwardUnfold,
         eigVecForwardUnfold, eigVecBackwardUnfold) \
            = ergoplot.readSpectrum(eigValForwardFileUnfold,
                                    eigValBackwardFileUnfold,
                                    eigVecForwardFileUnfold,
                                    eigVecBackwardFileUnfold,
                                    makeBiorthonormal=makeBiorthonormal,
                                    fileFormat=fileFormat)
        eigenConditionUnfold = ergoplot.getEigenCondition(
            eigVecForwardUnfold, eigVecBackwardUnfold, statDist)

        # Filter out based on eigen condition number
        maskCondition = eigenCondition > cfg.stat.maxCondition
        eigValForward = np.ma.masked_array(eigValForward, mask=maskCondition)
        maskCondition = np.tile(eigenCondition > cfg.stat.maxCondition,
                                (eigVecForward.shape[1], 1)).T
        eigVecForward = np.ma.masked_array(eigVecForward, mask=maskCondition)
        eigVecBackward = np.ma.masked_array(eigVecBackward, mask=maskCondition)

        maskCondition = eigenConditionUnfold > cfg.stat.maxCondition
        eigValForwardUnfold = np.ma.masked_array(
            eigValForwardUnfold, mask=maskCondition)
        maskCondition = np.tile(eigenConditionUnfold > cfg.stat.maxCondition,
                                (eigVecForwardUnfold.shape[1], 1)).T
        eigVecForwardUnfold = np.ma.masked_array(
            eigVecForwardUnfold, mask=maskCondition)
        eigVecBackwardUnfold = np.ma.masked_array(
            eigVecBackwardUnfold, mask=maskCondition)

        # Sort thanks to the distance between the eigenvectors
        print('Sorting...')
        eigVecForwardUnfold \
            = eigVecForwardUnfold[np.argsort(-np.abs(eigValForwardUnfold))]
        valid = np.nonzero(~eigValForward.mask)[0]
        isortUnfold = np.empty((valid.shape[0],), dtype=int)
        for iev, ev in enumerate(valid):
            ratio = (np.tile(eigVecForward[ev], (eigVecForwardUnfold.shape[0], 1))
                     / eigVecForwardUnfold)
            eigVecDist = np.std(ratio, 1) / np.mean(np.abs(ratio), 1)
            isortUnfold[iev] = np.ma.argmin(eigVecDist)
            print('{:03d} <-> {:03d} (dist = {:.12f})'.format(
                ev, isortUnfold[iev], eigVecDist[isortUnfold[iev]]))
        eigValForwardUnfoldSort = eigValForwardUnfold[isortUnfold]

        # Get generator eigenvalues
        print('Converting to generator eigenvalues...')
        eigValForward = np.array(eigValForward[valid])
        eigVecForward = np.array(eigVecForward[valid])
        eigVecBackward = np.array(eigVecBackward[valid])
        eigenCondition = np.array(eigenCondition[valid])
        eigValGen = ergoplot.eig2Generator(
            eigValForward, tau, eigValForwardUnfoldSort, tauUnfold)
    else:
        eigValGen = eigValGenOrig.copy()

    ergoplot.plotEig(eigValGen, xlabel=realLabel, ylabel=imagLabel,
                     xlim=xlimEig, ylim=ylimEig, xticks=xticks,
                     yticks=yticks, markersize=1, marker='o',
                     fig_ax=(fig, ax))
    # plt.scatter(eigValGen.real, eigValGen.imag, s=1)


fileName = 'eigval_bootstrap_nev{:d}{}.png'.format(nev, dstPostfixTau)
filePath = os.path.join(specDir, 'bootstrap', fileName)
fig.savefig(filePath, dpi=ergoplot.dpi, bbox_inches=ergoplot.bbox_inches)

plt.show(block=False)
