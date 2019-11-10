import os
import numpy as np
import matplotlib.pyplot as plt
import pylibconfig2
from ergoPack import ergoPlot

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
nSeeds = len(cfg.caseDef.seedRng)

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

obsIdx0 = cfg.stat.idxf
obsIdx1 = cfg.stat.idxg
compName0 = '%s%s' % (indicesName[0][0], fieldsDef[0][1])
ev_xlabel = '%s (%s)' % (compName0, fieldsDef[0][3])
compName1 = '%s%s' % (indicesName[1][0], fieldsDef[1][1])
ev_ylabel = '%s (%s)' % (compName1, fieldsDef[1][3])
corrName = 'C%d%d' % (obsIdx0, obsIdx1)
powerName = 'S%d%d' % (obsIdx0, obsIdx1)
corrLabel = r'$C_{x_%d, x_%d}(t)$' % (obsIdx0 + 1, obsIdx1 + 1)
powerLabel = r'$S_{x_%d, x_%d}(\omega)$' % (obsIdx0 + 1, obsIdx1 + 1)
realLabel = r'$\Re(\lambda_k)$'
imagLabel = r'$\Im(\lambda_k)$'
xlabelCorr = r'$t$'

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
zlimEig = [cfg.stat.powerMin, cfg.stat.powerMax]
xticks = None
yticksPos = np.arange(0, ylimEig[1], 5.)
yticksNeg = np.arange(0, ylimEig[0], -5.)[::-1]
yticks = np.concatenate((yticksNeg, yticksPos))
zticks = np.logspace(np.log10(zlimEig[0]), np.log10(zlimEig[1]),
                     int(np.round(np.log10(zlimEig[1]/zlimEig[0]) + 1)))
zticks = np.logspace(np.log10(zlimEig[0]), np.log10(zlimEig[1]),
                     int(np.round(np.log10(zlimEig[1]/zlimEig[0])/2 + 1)))

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
statDistFile = os.path.join(cfg.general.resDir, 'transfer', 'initDist',
                            fileName)
fileName = 'mask{}.{}'.format(dstPostfix, fileFormat)
maskFile = os.path.join(cfg.general.resDir, 'transfer', 'mask', fileName)

if unfold:
    tauUnfold = tauDimUnfold * timeScaleConversion
    dstPostfixTauUnfold = "{}_tau{:03d}".format(
        gridPostfix, int(tauDimUnfold * 1000 + 0.1))
    fileName = 'eigValForward_nev{:d}{}.{}'.format(
        nev, dstPostfixTauUnfold, fileFormat)
    eigValForwardFileUnfold = os.path.join(
        cfg.general.specDir, 'eigval', fileName)
    fileName = 'eigValBackward_nev{:d}{}.{}'.format(
        nev, dstPostfixTauUnfold, fileFormat)
    eigValBackwardFileUnfold = os.path.join(
        cfg.general.specDir, 'eigval', fileName)
    fileName = 'eigVecForward_nev{:d}{}.{}'.format(
        nev, dstPostfixTauUnfold, fileFormat)
    eigVecForwardFileUnfold = os.path.join(
        cfg.general.specDir, 'eigvec', fileName)
    fileName = 'eigVecBackward_nev{:d}{}.{}'.format(
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

# Read transfer operator spectrum from file and create a bi-orthonormal basis
# of eigenvectors and backward eigenvectors:
print('Readig spectrum for tauDim = {:.3f}...'.format(tauDim))
(eigValForward, eigValBackward, eigVecForward, eigVecBackward) \
    = ergoPlot.readSpectrum(eigValForwardFile, eigValBackwardFile,
                            eigVecForwardFile, eigVecBackwardFile,
                            makeBiorthonormal=makeBiorthonormal,
                            fileFormat=fileFormat)
eigenCondition = ergoPlot.getEigenCondition(eigVecForward, eigVecBackward,
                                            statDist)

eigValGenOrig = ergoPlot.eig2Generator(eigValForward, tau)

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
        = ergoPlot.readSpectrum(eigValForwardFileUnfold,
                                eigValBackwardFileUnfold,
                                eigVecForwardFileUnfold,
                                eigVecBackwardFileUnfold,
                                makeBiorthonormal=makeBiorthonormal,
                                fileFormat=fileFormat)
    eigenConditionUnfold = ergoPlot.getEigenCondition(
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
    eigValGen = ergoPlot.eig2Generator(
        eigValForward, tau, eigValForwardUnfoldSort, tauUnfold)
else:
    eigValGen = eigValGenOrig.copy()


# Define observables
# (f, g) = (T, T)
f = coord[cfg.stat.idxf][mask < N]
g = coord[cfg.stat.idxg][mask < N]
obsIdx0 = 0
obsIdx1 = 0
obsName0 = compName0
corrSamplePostfix = '_%s_%s_%s_%s_mu%04d_eps%04d' \
                    % (fieldsDef[obsIdx0][2], indicesName[obsIdx0][1],
                       fieldsDef[obsIdx1][2], indicesName[obsIdx1][1],
                       np.round(cfg.caseDef.mu * 1000, 1),
                       np.round(cfg.caseDef.eps * 1000, 1))
corrName = 'C%d%d' % (obsIdx0, obsIdx1)
powerName = 'S%d%d' % (obsIdx0, obsIdx1)
corrLabel = r'$C_{\mathrm{%s}}(t)$' % obsName0
powerLabel = r'$S_{\mathrm{%s}}(z)$' % obsName0
realLabel = r'$\Re(\lambda_k)$'
imagLabel = r'$\Im(\lambda_k)$'
xlabelCorr = r'$t$'

# Get covariance
mean_f = (f * statDist).sum()
mean_g = (statDist * np.conjugate(g)).sum()
cfg0 = ((f - mean_f) * statDist * (g - mean_g)).sum()

# Read ccf
fileName = 'corrSample{}_nSeeds{:d}_lagMax{:d}yr.txt'.format(
    corrSamplePostfix, nSeeds, cfg.stat.lagMax)
filePath = os.path.join(cfg.general.resDir, 'correlation', fileName)
corrSample = np.loadtxt(filePath)
fileName = 'lags{}_nSeeds{:d}_lagMax{:d}yr.txt'.format(
    corrSamplePostfix, nSeeds, cfg.stat.lagMax)
filePath = os.path.join(cfg.general.resDir, 'correlation', fileName)
lags = np.loadtxt(filePath)
fileName = 'powerSample{}_nSeeds{:d}_chunk{:d}yr.txt'.format(
    corrSamplePostfix, nSeeds, cfg.stat.chunkWidth)
filePath = os.path.join(cfg.general.resDir, 'power', fileName)
powerSample = np.loadtxt(filePath)
fileName = 'freq{}_nSeeds{:d}_chunk{:d}yr.txt'.format(
    corrSamplePostfix, nSeeds, cfg.stat.chunkWidth)
filePath = os.path.join(cfg.general.resDir, 'power', fileName)
freq = np.loadtxt(filePath)
# Convert to angular frequencies
angFreq = freq * 2*np.pi
powerSample /= 2*np.pi


# Reconstruct correlation and power spectrum
# Get normalized weights
weights = ergoPlot.getSpectralWeights(f, g, eigVecForward, eigVecBackward)

# Remove components with heigh condition number
weights[eigenCondition > cfg.stat.maxCondition] = 0.
# weights[3:] = 0.
# weights = np.abs(weights)
condition = np.empty(eigenCondition.shape, dtype='|U1')
condition[:] = 'k'
condition[eigenCondition > cfg.stat.maxCondition - 0.001] = 'w'
(corrRec, compCorrRec) = ergoPlot.spectralRecCorrelation(lags, eigValGen,
                                                         weights,
                                                         norm=cfg.stat.norm)
(powerRec, compPowerRec) = ergoPlot.spectralRecPower(angFreq, eigValGen,
                                                     weights,
                                                     norm=cfg.stat.norm)

# Plot correlation reconstruction
fig = ergoPlot.plotRecCorrelation(lags, corrSample, corrRec,
                                  plotPositive=True,
                                  ylabel=corrLabel, xlabel=xlabelCorr)
fileName = '{}Rec_lag{:d}_nev{:d}{}.{}'.format(
    corrName, int(cfg.stat.lagMax), nev, dstPostfixTau, ergoPlot.figFormat)
filePath = os.path.join(specDir, 'reconstruction', fileName)
fig.savefig(filePath, dpi=ergoPlot.dpi, bbox_inches=ergoPlot.bbox_inches)

# PLot spectrum, powerSampledogram and spectral reconstruction
w = weights.copy()
if cfg.stat.norm:
    w /= w[1:].sum()
msize = np.zeros((w.shape[0]))
msize[w.real > 0] = np.log10(w[w.real > 0].real)
msize[w.real > 0] = (msize[w.real > 0] + 8) * 10
# msize[w.real > 0] = (msize[w.real > 0] + 6) * 3
msize[msize < 0] = 0.
msize[0] = 12.
fig = ergoPlot.plotEigPowerRec(angFreq, eigValGen, powerSample, powerRec,
                               markersize=msize, condition=condition,
                               xlabel=realLabel, ylabel=imagLabel,
                               zlabel=powerLabel,
                               xlim=xlimEig, ylim=ylimEig, zlim=zlimEig,
                               xticks=xticks, yticks=yticks, zticks=zticks)
fileName = '{}Rec_chunk{:d}_nev{:d}{}.{}'.format(
    powerName, int(cfg.stat.chunkWidth), nev, dstPostfixTau,
    ergoPlot.figFormat)
filePath = os.path.join(specDir, 'reconstruction', fileName)
fig.savefig(filePath, dpi=ergoPlot.dpi, bbox_inches=ergoPlot.bbox_inches)

nw = weights[1:]
nw /= nw.sum()
nw = np.abs(nw)

print(nw / eigenCondition[1:])

plt.show(block=False)
