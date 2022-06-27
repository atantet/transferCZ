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

nevPlot = 10

# Transition lag
tauDimRng = np.arange(
    cfg.transfer.lag0,
    cfg.transfer.lag0 + (cfg.transfer.nLags - 1) * cfg.transfer.stepLag,
    cfg.transfer.stepLag)
unfold = False
if hasattr(cfg.stat, 'unfold'):
    if cfg.stat.unfold:
        unfold = True

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

nev = cfg.spectrum.nev
xmineigVal = -0.3
xlimEig = [xmineigVal, -xmineigVal/100]

srcPostfix = '{}{}_mu{:04d}_eps{:04d}'.format(
    cfg.caseDef.prefix, cfg.caseDef.simType,
    int(cfg.caseDef.mu * 1000 + 0.1),
    int(cfg.caseDef.eps * 1000 + 0.1))
print(srcPostfix)

# Obs name
obsName = ''.join('_{}_{}'.format(fieldsDef[d][2], indicesName[d][1])
                  for d in np.arange(dimObs))
# Seed postfix
nSeeds = len(cfg.caseDef.seedRng)
seedPostfix = '_seeds' + ''.join(str(s) for s in cfg.caseDef.seedRng)

gridPostfix = ''
N = 1
for d in np.arange(dimObs):
    N *= cfg.grid.nx[d]
    gridPostfix = '{}_n{:d}l{:d}h{:d}'.format(
        gridPostfix, cfg.grid.nx[d], int(cfg.grid.nSTDLow[d]),
        int(cfg.grid.nSTDHigh[d]))
cpyBuffer = gridPostfix
gridPostfix = '_{}{}{}'.format(srcPostfix, obsName, cpyBuffer)

# Read grid
gridFile = '{}/grid/grid{}.txt'.format(cfg.general.resDir, gridPostfix)
coord = ergoplot.readGrid(gridFile, dimObs)
if dimObs == 1:
    X = coord[0]
elif dimObs == 2:
    X, Y = np.meshgrid(coord[0], coord[1], indexing='ij')
    coord = (X.flatten(), Y.flatten())
elif dimObs == 3:
    X, Y, Z = np.meshgrid(coord[0], coord[1], coord[2], indexing='ij')
    coord = (X.flatten(), Y.flatten(), Z.flatten())
dstPostfix = '{}{}'.format(gridPostfix, seedPostfix)

eigValGenRng = np.ma.masked_all((tauDimRng.shape[0], nev))
diffRng = np.ma.masked_all((tauDimRng.shape[0], nev))
for itau, tauDim in enumerate(tauDimRng):
    tauDimUnfold = tauDim + 1
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
    (eigValForward, eigValBackward, eigVecForward, eigVecBackward)\
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
         eigVecForwardUnfold, eigVecBackwardUnfold)\
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
        # print('Sorting...')
        eigVecForwardUnfold\
            = eigVecForwardUnfold[np.argsort(-np.abs(eigValForwardUnfold))]
        valid = np.nonzero(~eigValForward.mask)[0]
        isortUnfold = np.empty((valid.shape[0],), dtype=int)
        for iev, ev in enumerate(valid):
            ratio = (np.tile(eigVecForward[ev], (eigVecForwardUnfold.shape[0], 1))
                     / eigVecForwardUnfold)
            eigVecDist = np.std(ratio, 1) / np.mean(np.abs(ratio), 1)
            isortUnfold[iev] = np.ma.argmin(eigVecDist)
            # print('{:03d} <-> {:03d} (dist = {:.12f})'.format(
            #     ev, isortUnfold[iev], eigVecDist[isortUnfold[iev]]))
        eigValForwardUnfoldSort = eigValForwardUnfold[isortUnfold]

        # Get generator eigenvalues
        # print('Converting to generator eigenvalues...')
        eigValForward = np.array(eigValForward[valid])
        eigVecForward = np.array(eigVecForward[valid])
        eigVecBackward = np.array(eigVecBackward[valid])
        eigenCondition = np.array(eigenCondition[valid])
        eigValGen = ergoplot.eig2Generator(
            eigValForward, tau, eigValForwardUnfoldSort, tauUnfold)
    else:
        eigValGen = eigValGenOrig.copy()

    eigValGen = eigValGen[:nev]
    eigValGenRng[itau, :len(eigValGen)] = eigValGen.real

fig, ax = plt.subplots()
for ev in range(nevPlot):
    ax.plot(tauDimRng, eigValGenRng[:, ev])
ax.set_xlim(*tauDimRng[[0, -1]])
ax.set_ylim(xlimEig)
ax.vlines([15.], *xlimEig, linestyle='--', linewidth=1, colors='k')
ax.set_ylabel(realLabel, fontsize=ergoplot.fs_latex)
ax.set_xlabel(r'$\tau$', fontsize=ergoplot.fs_latex)
fileName = 'eigval_tau_nev{:d}{}.png'.format(nev, dstPostfixTau)
filePath = os.path.join(specDir, 'bootstrap', fileName)
fig.savefig(filePath, dpi=ergoplot.dpi, bbox_inches=ergoplot.bbox_inches)

fig, ax = plt.subplots()
for ev in range(1, nevPlot):
    w = eigValGenRng[:, ev]
    diff = np.abs((w[:-1] - w[1:]) / w[:-1] * 100)
    diffRng[1:len(diff) + 1, ev] = diff
    ax.plot(tauDimRng[1:], diff)
ax.set_xlim(*tauDimRng[[0, -1]])
ylim = [0., 50.]
ax.set_ylim(*ylim)
ax.hlines([1.], *tauDimRng[[0, -1]], linestyle='--', linewidth=1, colors='k')
ax.vlines([100.], *ylim, linestyle='--', linewidth=1, colors='k')
ax.set_ylabel(r'$\Delta \Re(\lambda_k)$ (%)', fontsize=ergoplot.fs_latex)
ax.set_xlabel(r'$n_x$', fontsize=ergoplot.fs_latex)
fileName = 'eigval_taudiff_nev{:d}{}.png'.format(nev, dstPostfixTau)
filePath = os.path.join(specDir, 'bootstrap', fileName)
fig.savefig(filePath, dpi=ergoplot.dpi, bbox_inches=ergoplot.bbox_inches)

plt.show(block=False)
