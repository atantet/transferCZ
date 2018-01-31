import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import pylibconfig2
import ergoPlot

configFile = '../cfg/transferCZ.cfg'
cfg = pylibconfig2.Config()
cfg.read_file(configFile)

#tauDim = 1.
tauDim = 3.
#tauDim = 6.

timeScaleConversion = 1. / 12
dim = len(cfg.caseDef.indicesName)
nSeeds = len(cfg.caseDef.seedRng)

nev = cfg.spectrum.nev
nevPlot = 0
#plotBackward = False
plotBackward = True
#plotImag = False
plotImag = True
xminEigVal = -1.2
yminEigVal = -11.5
ev_xlabel = r'$x_1$'
ev_ylabel = r'$x_2$'
xlimEig = [cfg.stat.rateMin, -cfg.stat.rateMin / 100]
ylimEig = [-cfg.stat.angFreqMax, cfg.stat.angFreqMax]
zlimEig = [cfg.stat.yminPower, cfg.stat.ymaxPower]
xticks = None
yticksPos = np.arange(0, ylimEig[1], 2.5)
yticksNeg = np.arange(0, ylimEig[0], -2.5)[::-1]
yticks = np.concatenate((yticksNeg, yticksPos))
zticks = np.logspace(np.log10(zlimEig[0]), np.log10(zlimEig[1]),
                     int(np.round(np.log10(zlimEig[1]/zlimEig[0]) + 1)))
maxCondition = 5


field_h = (1, 'thermocline depth', 'h', 'm')
field_T = (2, 'SST', 'T', r'$^\circ C$')
field_u_A = (3, 'wind stress due to coupling', 'u_A', 'm/s')
field_taux = (4, 'external wind-stress', 'taux', 'm/s')

nino3 = ('Eastern', 'nino3')
nino4 = ('Western', 'nino4')
indicesName = []
fieldsDef = []
for d in np.arange(dim):
    if cfg.caseDef.indicesName[d] == 'nino3':
        indicesName.append(nino3)
    elif cfg.caseDef.indicesName[d] == 'nino4':
        indicesName.append(nino4)
    if cfg.caseDef.fieldsName[d] == 'T':
        fieldsDef.append(field_T)
    if cfg.caseDef.fieldsName[d] == 'h':
        fieldsDef.append(field_h)

srcPostfix = "%s%s_mu%04d_eps%04d" % (cfg.caseDef.prefix, cfg.caseDef.simType,
                                      np.round(cfg.caseDef.mu * 1000, 1),
                                      np.round(cfg.caseDef.eps * 1000, 1))
obsName = ''
gridPostfix = ''
N = 1
for d in np.arange(dim):
    obsName += '_%s_%s' % (fieldsDef[d][2], indicesName[d][1])
    N *= cfg.grid.nx[d]
    gridPostfix = "%s_n%dl%dh%d" % (gridPostfix, cfg.grid.nx[d],
                                    cfg.grid.nSTDLow[d], cfg.grid.nSTDHigh[d])
cpyBuffer = gridPostfix
gridPostfix = '_%s%s%s' % (srcPostfix, obsName, cpyBuffer)

# Read grid
gridFile = '%s/grid/grid%s.txt' % (cfg.general.resDir, gridPostfix)
(X, Y) = ergoPlot.readGrid(gridFile, dim)
coord = (X.flatten(), Y.flatten())

tauConv = tauDim * timeScaleConversion
postfix = "%s_tau%03d" % (gridPostfix, tauDim * 1000)

print 'Readig spectrum...'
EigValForwardFile = '%s/eigval/eigValForward_nev%d%s.txt' % (cfg.general.specDir, nev, postfix)
EigVecForwardFile = '%s/eigvec/eigVecForward_nev%d%s.txt' % (cfg.general.specDir, nev, postfix)
EigValBackwardFile = '%s/eigval/eigValBackward_nev%d%s.txt' \
                    % (cfg.general.specDir, nev, postfix)
EigVecBackwardFile = '%s/eigvec/eigVecBackward_nev%d%s.txt' \
                    % (cfg.general.specDir, nev, postfix)
statDist = np.loadtxt('%s/transfer/initDist/initDist%s.txt' % (cfg.general.resDir, gridPostfix))
(eigValForward, eigVecForward, eigValBackward, eigVecBackward) \
    = ergoPlot.readSpectrum(EigValForwardFile, EigVecForwardFile,
                            EigValBackwardFile, EigVecBackwardFile,
                            statDist, makeBiorthonormal=False)

print 'Getting conditionning of eigenvectors...'
eigenCondition = ergoPlot.getEigenCondition(eigVecForward, eigVecBackward, statDist)

# Get generator eigenvalues
eigValGen = (np.log(np.abs(eigValForward)) + np.angle(eigValForward)*1j) / tauConv

# Plot eigenvectors
alpha = 0.01
for ev in np.arange(nevPlot):
    print 'Plotting real part of eigenvector %d...' % (ev+1,)
    ergoPlot.plot2D(X, Y, eigVecForward[:, ev].real, ev_xlabel, ev_ylabel, alpha)
    plt.savefig('%s/spectrum/eigvec/eigVecForwardReal_nev%d_ev%03d%s.%s' \
                % (cfg.general.plotDir, nev, ev+1, postfix, ergoPlot.figFormat),
                bbox_inches='tight', dpi=ergoPlot.dpi)
    
    if plotImag & (eigValForward[ev].imag != 0):
        print 'Plotting imaginary  part of eigenvector %d...' % (ev+1,)
        ergoPlot.plot2D(X, Y, eigVecForward[:, ev].imag, ev_xlabel, ev_ylabel, alpha)
        plt.savefig('%s/spectrum/eigvec/eigVecForwardImag_nev%d_ev%03d%s.%s' \
                    % (cfg.general.plotDir, nev, ev+1, postfix, ergoPlot.figFormat),
                    bbox_inches='tight', dpi=ergoPlot.dpi)
    
    # Plot eigenvectors of backward operator
    if plotBackward:
        print 'Plotting real part of backward eigenvector %d...' % (ev+1,)
        ergoPlot.plot2D(X, Y, eigVecBackward[:, ev].real, ev_xlabel, ev_ylabel, alpha)
        plt.savefig('%s/spectrum/eigvec/eigVecBackwardReal_nev%d_ev%03d%s.%s' \
                    % (cfg.general.plotDir, nev, ev+1, postfix, ergoPlot.figFormat),
                    bbox_inches='tight', dpi=ergoPlot.dpi)
        
        if plotImag & (eigValForward[ev].imag != 0):
            print 'Plotting imaginary  part of backward eigenvector %d...' % (ev+1,)
            ergoPlot.plot2D(X, Y, eigVecBackward[:, ev].imag, ev_xlabel, ev_ylabel, alpha)
            plt.savefig('%s/spectrum/eigvec/eigVecBackwardImag_nev%d_ev%03d%s.%s' \
                        % (cfg.general.plotDir, nev, ev+1, postfix, ergoPlot.figFormat),
                        bbox_inches='tight', dpi=ergoPlot.dpi)

# Get ccf
#    Get sample cross-correlation
print 'Reading corrSample and perio'
f = X.flatten()
g = f
obsIdx0 = 0
obsIdx1 = 0
corrSamplePath = '../results/%s/' % srcPostfix
corrSamplePostfix = '_%s_%s_%s_%s_mu%04d_eps%04d' \
                    % (fieldsDef[obsIdx0][2], indicesName[obsIdx0][1],
                       fieldsDef[obsIdx1][2], indicesName[obsIdx1][1],
                       np.round(cfg.caseDef.mu * 1000, 1), np.round(cfg.caseDef.eps * 1000, 1))
corrName = 'C%d%d' % (obsIdx0, obsIdx1)
powerName = 'S%d%d' % (obsIdx0, obsIdx1)
corrLabel = r'$C_{x_%d, x_%d}(t)$' % (obsIdx0 + 1, obsIdx1 + 1)
powerLabel = r'$S_{x_%d, x_%d}(\omega)$' % (obsIdx0 + 1, obsIdx1 + 1)
realLabel = r'$\Re(\bar{\lambda}_k)$'
imagLabel = r'$\Im(\bar{\lambda}_k)$'

corrSample = np.loadtxt('%s/corrSample%s_nSeeds%d_lagMax%dyr.txt' \
                        % (corrSamplePath, corrSamplePostfix, nSeeds, cfg.stat.lagMax))
lags = np.loadtxt('%s/lags%s_nSeeds%d_lagMax%dyr.txt' \
                  % (corrSamplePath, corrSamplePostfix, nSeeds, cfg.stat.lagMax))
powerSample = np.loadtxt('%s/powerSample%s_nSeeds%d_chunk%dyr.txt' \
                         % (corrSamplePath, corrSamplePostfix, nSeeds, cfg.stat.chunkWidth))
freq = np.loadtxt('%s/freq%s_nSeeds%d_chunk%dyr.txt' \
                  % (corrSamplePath, corrSamplePostfix, nSeeds, cfg.stat.chunkWidth))

angFreq = freq * 2*np.pi
cfg0 = ((f - (f * statDist).sum()) * statDist * (g - (g * statDist).sum())).sum()
powerSample /= 2 * np.pi * cfg0

# Reconstruct correlation and power spectrum
# Get normalized weights
weights = ergoPlot.getSpectralWeights(f, g, eigVecForward, eigVecBackward, statDist, skipMean=True)
# Remove components with heigh condition number
weights[eigenCondition > maxCondition] = 0.
eigenCondition[eigenCondition > maxCondition] = maxCondition
(corrRec, compCorrRec) = ergoPlot.spectralRecCorrelation(lags, f, g, eigValGen, weights,
                                                         statDist, skipMean=True, norm=True)
(powerRec, compPowerRec) = ergoPlot.spectralRecPower(angFreq, f, g, eigValGen, weights, statDist, norm=True)

# Plot correlation reconstruction
ergoPlot.plotRecCorrelation(lags, corrSample, corrRec, plotPositive=True,
                            ylabel=corrLabel)
plt.savefig('%s/spectrum/reconstruction/%sRec_lag%d_nev%d%s.%s'\
            % (cfg.general.plotDir, corrName, int(cfg.stat.lagMax),
               nev, postfix, ergoPlot.figFormat),
            dpi=ergoPlot.dpi, bbox_inches=ergoPlot.bbox_inches)

# PLot spectrum, powerSampledogram and spectral reconstruction
msize = np.zeros((weights.shape[0]))
msize[weights.real > 0] = np.log(weights[weights.real > 0].real)
msize[weights.real > 0] = (msize[weights.real > 0] + 15) * 10
# msize[weights.real > 0] = (msize[weights.real > 0] + 6) * 3
msize[msize < 0] = 0.
#msize = 20
ergoPlot.plotEigPowerRec(angFreq, eigValGen, powerSample, powerRec,
                         markersize=msize, condition=eigenCondition,
                         xlabel=realLabel, ylabel=imagLabel, zlabel=powerLabel,
                         xlim=xlimEig, ylim=ylimEig, zlim=zlimEig,
                         xticks=xticks, yticks=yticks, zticks=zticks)
plt.savefig('%s/spectrum/reconstruction/%sRec_nev%d%s.%s'\
            % (cfg.general.plotDir, powerName,
               nev, postfix, ergoPlot.figFormat),
            dpi=ergoPlot.dpi, bbox_inches=ergoPlot.bbox_inches)

