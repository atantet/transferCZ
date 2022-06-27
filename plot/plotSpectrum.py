import os
import numpy as np
import matplotlib.pyplot as plt
import pylibconfig2
from ergopack import ergoplot
from matplotlib import cm
import os

configFile = os.path.join('..', 'cfg', 'transferCZ.cfg')
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
#nevPlot = 6
nevPlot = 0
plotForward = False
#plotForward = True
#plotBackward = False
plotBackward = True
plotImag = False
#plotImag = True
plotPolar = True
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

# Get generator eigenvalues
eigValGen = (np.log(np.abs(eigValForward)) + np.angle(eigValForward)*1j) / tau


# Plot eigenvectors of transfer operator
alpha = 0.05
plotDir = os.path.join(cfg.general.plotDir, 'spectrum')
os.makedirs(os.path.join(plotDir, 'eigvec'), exist_ok=True)
os.makedirs(os.path.join(plotDir, 'reconstruction'), exist_ok=True)
for ev in np.arange(nevPlot):
    if ev == 0:
        positive = True
        cmap = cm.hot_r
    else:
        positive = False
        cmap = cm.RdBu_r
    if plotPolar:
        if plotForward:
            print('Plotting polar eigenvector {:d}...'.format(ev + 1))
            (figPhase, figAmp) \
                = ergoplot.plotEigVecPolar(X, Y, eigVecForward[ev],
                                           mask=mask, combine=False,
                                           xlabel=ev_xlabel, ylabel=ev_ylabel,
                                           alpha=alpha, cmap=cmap,
                                           ampMin=ampMin, ampMax=ampMax)
            fileName = 'eigvecForwardPhase_nev{:d}_ev{:03d}{}.{}'.format(
                nev, ev + 1, dstPostfixTau, ergoplot.figFormat)
            dstFilePhase = os.path.join(specDir, 'eigvec', fileName)
            fileName = 'eigvecForwardAmp_nev{:d}_ev{:03d}{}.{}'.format(
                nev, ev + 1, dstPostfixTau, ergoplot.figFormat)
            dstFileAmp = os.path.join(specDir, 'eigvec', fileName)
            figPhase.savefig(dstFilePhase, bbox_inches=ergoplot.bbox_inches,
                           dpi=ergoplot.dpi)
            figAmp.savefig(dstFileAmp, bbox_inches=ergoplot.bbox_inches,
                           dpi=ergoplot.dpi)
        if plotBackward:
            print('Plotting polar backward eigenvector {:d}...'.format(ev + 1))
            (figPhase, figAmp) \
                = ergoplot.plotEigVecPolar(X, Y, eigVecBackward[ev],
                                           mask=mask, combine=False,
                                           xlabel=ev_xlabel, ylabel=ev_ylabel,
                                           alpha=alpha, cmap=cmap,
                                           ampMin=ampMin, ampMax=ampMax)
            fileName = 'eigvecBackwardPhase_nev{:d}_ev{:03d}{}.{}'.format(
                nev, ev + 1, dstPostfixTau, ergoplot.figFormat)
            dstFilePhase = os.path.join(specDir, 'eigvec', fileName)
            fileName = 'eigvecBackwardAmp_nev{:d}_ev{:03d}{}.{}'.format(
                nev, ev + 1, dstPostfixTau, ergoplot.figFormat)
            dstFileAmp = os.path.join(specDir, 'eigvec', fileName)
            figPhase.savefig(dstFilePhase, bbox_inches=ergoplot.bbox_inches,
                           dpi=ergoplot.dpi)
            figAmp.savefig(dstFileAmp, bbox_inches=ergoplot.bbox_inches,
                           dpi=ergoplot.dpi)

    else:
        if plotForward:
            print('Plotting real part of forward eigenvector {:d}...'.format(ev + 1))
            fig = ergoplot.plotEigVec(X, Y, eigVecForward[ev].real,
                                      xlabel=ev_xlabel, ylabel=ev_ylabel,
                                      mask=mask, alpha=alpha,
                                      positive=positive, cmap=cmap)
            fileName = 'eigvecForwardReal_nev{:d}_ev{:03d}{}.{}'.format(
                nev, ev + 1, dstPostfixTau, ergoplot.figFormat)
            dstFile = os.path.join(specDir, 'eigvec', fileName)
            fig.savefig(dstFile, bbox_inches=ergoplot.bbox_inches,
                        dpi=ergoplot.dpi)
            if plotImag & (eigValForward[ev].imag != 0):
                print('Plotting imaginary part of forward eigenvector',
                      '{:d}...'.format(ev + 1))
                fig = ergoplot.plotEigVec(X, Y, eigVecForward[ev].imag,
                                          mask=mask, cmap=cmap,
                                          xlabel=ev_xlabel, ylabel=ev_ylabel,
                                          positive=positive, alpha=alpha)
                fileName = 'eigvecForwardImag_nev{:d}_ev{:03d}{}.{}'.format(
                    nev, ev + 1, dstPostfixTau, ergoplot.figFormat)
                dstFile = os.path.join(specDir, 'eigvec', fileName)
                fig.savefig(dstFile, bbox_inches=ergoplot.bbox_inches,
                            dpi=ergoplot.dpi)
    
        # Plot eigenvectors of backward operator
        if plotBackward:
            print('Plotting real part of backward eigenvector {:d}...'.format(ev + 1))
            fig = ergoplot.plotEigVec(X, Y, eigVecBackward[ev].real,
                                      mask=mask, cmap=cmap,
                                      xlabel=ev_xlabel, ylabel=ev_ylabel,
                                      positive=positive, alpha=alpha)
            fileName = 'eigvecBackwardReal_nev{:d}_ev{:03d}{}.{}'.format(
                nev, ev + 1, dstPostfixTau, ergoplot.figFormat)
            dstFile = os.path.join(specDir, 'eigvec', fileName)
            fig.savefig(dstFile, bbox_inches=ergoplot.bbox_inches,
                        dpi=ergoplot.dpi)
            
            if plotImag & (eigValForward[ev].imag != 0):
                print('Plotting imag. part of backward eigenvector',
                      '{:d}...'.format(ev + 1))
                fig = ergoplot.plotEigVec(X, Y, eigVecBackward[ev].imag,
                                          mask=mask, cmap=cmap,
                                          xlabel=ev_xlabel, ylabel=ev_ylabel,
                                          positive=positive, alpha=alpha)
                fileName = 'eigvecBackwardImag_nev{:d}_ev{:03d}{}.{}'.format(
                    nev, ev + 1, dstPostfixTau, ergoplot.figFormat)
                dstFile = os.path.join(specDir, 'eigvec', fileName)
                fig.savefig(dstFile, bbox_inches=ergoplot.bbox_inches,
                            dpi=ergoplot.dpi)

            
# Define observables
print('Reading corrSample and perio')
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
weights = ergoplot.getSpectralWeights(f, g, eigVecForward, eigVecBackward)

# Remove components with heigh condition number
weights[eigenCondition > cfg.stat.maxCondition] = 0.
# weights[3:] = 0.
#weights = np.abs(weights)
condition = np.empty(eigenCondition.shape, dtype='|U1')
condition[:] = 'k'
condition[eigenCondition > cfg.stat.maxCondition - 0.001] = 'w'
(corrRec, compCorrRec) = ergoplot.spectralRecCorrelation(lags, eigValGen,
                                                         weights,
                                                         norm=cfg.stat.norm)
(powerRec, compPowerRec) = ergoplot.spectralRecPower(angFreq, eigValGen,
                                                     weights,
                                                     norm=cfg.stat.norm)

# Plot correlation reconstruction
fig = ergoplot.plotRecCorrelation(lags, corrSample, corrRec,
                                  plotPositive=True,
                                  ylabel=corrLabel, xlabel=xlabelCorr)
fileName = '{}Rec_lag{:d}_nev{:d}{}.{}'.format(
    corrName, int(cfg.stat.lagMax), nev, dstPostfixTau, ergoplot.figFormat)
filePath = os.path.join(specDir, 'reconstruction', fileName)
fig.savefig(filePath, dpi=ergoplot.dpi, bbox_inches=ergoplot.bbox_inches)

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
fig = ergoplot.plotEigPowerRec(angFreq, eigValGen, powerSample, powerRec,
                               markersize=msize, condition=condition,
                               xlabel=realLabel, ylabel=imagLabel,
                               zlabel=powerLabel,
                               xlim=xlimEig, ylim=ylimEig, zlim=zlimEig,
                               xticks=xticks, yticks=yticks, zticks=zticks)
fileName = '{}Rec_chunk{:d}_nev{:d}{}.{}'.format(
    powerName, int(cfg.stat.chunkWidth), nev, dstPostfixTau, ergoplot.figFormat)
filePath = os.path.join(specDir, 'reconstruction', fileName)
fig.savefig(filePath, dpi=ergoplot.dpi, bbox_inches=ergoplot.bbox_inches)

nw = weights[1:]
nw /= nw.sum()
nw = np.abs(nw)

print(nw / eigenCondition[1:])

plt.show(block=False)
