import os
import numpy as np
import matplotlib.pyplot as plt
import pylibconfig2
import ergoPlot

configFile = '../cfg/transferCZ.cfg'
cfg = pylibconfig2.Config()
cfg.read_file(configFile)

# Transition lag
if (hasattr(cfg.stat, 'tauDimPlot')):
    tauDim = cfg.stat.tauDimPlot
else:
    tauDim = cfg.transfer.tauDimRng[0]
unfold = False
tauDimUnfold = None
if (hasattr(cfg.stat, 'unfold')) & (hasattr(cfg.stat, 'tauDimUnfold')):
    if cfg.stat.unfold:
        unfold = True
        tauDimUnfold = cfg.stat.tauDimUnfold

timeScaleConversion = 1. / 12
dimObs = len(cfg.caseDef.indicesName)
nSeeds = len(cfg.caseDef.seedRng)

obsIdx0 = cfg.stat.idxf
obsIdx1 = cfg.stat.idxg
ev_xlabel = r'East SST'
ev_ylabel = r'West TC'
corrName = 'C%d%d' % (obsIdx0, obsIdx1)
powerName = 'S%d%d' % (obsIdx0, obsIdx1)
corrLabel = r'$C_{x_%d, x_%d}(t)$' % (obsIdx0 + 1, obsIdx1 + 1)
powerLabel = r'$S_{x_%d, x_%d}(\omega)$' % (obsIdx0 + 1, obsIdx1 + 1)
realLabel = r'$\Re(\lambda_k)$'
imagLabel = r'$\Im(\lambda_k)$'
xlabelCorr = r'$t$'

nevPlot = 0
nev = cfg.spectrum.nev
#plotBackward = False
plotBackward = True
#plotImag = False
#plotImag = True
plotPolar = True
xmineigVal = -cfg.stat.rateMax
ymineigVal = -cfg.stat.angFreqMax
#plotBackward = False
plotBackward = True
plotImag = False
#plotImag = True
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
    X, Y = np.meshgrid(coord[0], coord[1])
    coord = (X.flatten(), Y.flatten())
elif dimObs == 3:
    X, Y, Z = np.meshgrid(coord[0], coord[1], coord[2], indexing='ij')
    coord = (X.T.flatten(), Y.T.flatten(), Z.T.flatten())

tau = tauDim * timeScaleConversion
postfix = "%s_tau%03d" % (gridPostfix, tauDim * 1000)

# File names
eigValForwardFile = '%s/eigval/eigValForward_nev%d%s.txt' % (cfg.general.specDir, nev, postfix)
eigVecForwardFile = '%s/eigvec/eigVecForward_nev%d%s.txt' % (cfg.general.specDir, nev, postfix)
eigValBackwardFile = '%s/eigval/eigValBackward_nev%d%s.txt' \
                    % (cfg.general.specDir, nev, postfix)
eigVecBackwardFile = '%s/eigvec/eigVecBackward_nev%d%s.txt' \
                    % (cfg.general.specDir, nev, postfix)
statDistFile = '%s/transfer/initDist/initDist%s.txt' \
               % (cfg.general.resDir, gridPostfix)
if unfold:
    tauUnfold = tauDimUnfold * timeScaleConversion
    postfixUnfold = "%s_tau%03d" % (gridPostfix, tauDimUnfold * 1000)
    eigValForwardFileUnfold = '%s/eigval/eigvalForward_nev%d%s.txt' \
                              % (cfg.general.specDir, nev, postfixUnfold)
    eigVecForwardFileUnfold = '%s/eigvec/eigvecForward_nev%d%s.txt' \
                              % (cfg.general.specDir, nev, postfixUnfold)
    tauMax = np.max([tau, tauUnfold])

# Read transfer operator spectrum from file and create a bi-orthonormal basis
# of eigenvectors and backward eigenvectors:
print 'Readig spectrum for tauDim = %.3f...' % tauDim
(eigValForward, eigValBackward, statDist, eigVecForward, eigVecBackward) \
    = ergoPlot.readSpectrum(eigValForwardFile, eigValBackwardFile, statDistFile,
                            eigVecForwardFile, eigVecBackwardFile,
                            makeBiorthonormal=~cfg.spectrum.makeBiorthonormal)

print 'Getting conditionning of eigenvectors...'
eigenCondition = ergoPlot.getEigenCondition(eigVecForward, eigVecBackward, statDist)

eigValGenOrig = ergoPlot.eig2Generator(eigValForward, tau)

# Read unfolding eigenvalues
#eigValForwardUnfold = None
st = statDist.copy()
st2 = np.concatenate((st, st))
if unfold:
    print 'Readig spectrum for tauDimUnfold = %.3f to unfold...' % tauDimUnfold
    (eigValForwardUnfold, eigVecForwardUnfold) \
        = ergoPlot.readSpectrum(eigValForwardFileUnfold,
                                eigVecForwardFile=eigVecForwardFileUnfold)
    eigVecForwardUnfold *= 0.1
    eigValGenShift = np.log(np.abs(eigValForward)) / tau \
                     + 1j*(np.mod(np.angle(eigValForward) / tau+np.pi/tauMax,
                                  2*np.pi/tauMax)-np.pi/tauMax)
    
    # Sort thanks to the distance between the eigenvectors
    print 'Sorting...'
    isortUnfold = np.empty((nev,), dtype=int)
    eigVecForwardUnfold \
        = eigVecForwardUnfold[:, np.argsort(-np.abs(eigValForwardUnfold))]
    for ev in np.arange(nev):
        eigVecDist = np.empty((nev,))
        vec = eigVecForward[:, ev].copy()
        B = np.matrix(np.empty((N*2, 1)))
        B[:N, 0] = np.matrix(vec.real).T
        B[N:, 0] = np.matrix(vec.imag).T
        Ba = np.array(B)[:, 0]
        for iev in np.arange(nev):
            # Get their distance to a phase shift
            vecUnfold = eigVecForwardUnfold[:, iev].copy()
            Z = np.matrix(np.empty((N*2, 2)))
            Z[:N, 0] = np.matrix(vecUnfold.real).T
            Z[N:, 0] = np.matrix(vecUnfold.imag).T
            Z[:N, 1] = -np.matrix(vecUnfold.imag).T
            Z[N:, 1] = np.matrix(vecUnfold.real).T
            Zd = np.matrix(np.empty((2, N*2)))
            Zd[0, :N] = np.matrix(vecUnfold.real * st)
            Zd[0, N:] = np.matrix(vecUnfold.imag * st)
            Zd[1, :N] = -np.matrix(vecUnfold.imag * st)
            Zd[1, :N] = np.matrix(vecUnfold.real * st)
            ZdZ = Zd * Z
            detZdZ = ZdZ[0,0]*ZdZ[1,1] - ZdZ[0,1]*ZdZ[1,0]
            if np.abs(detZdZ) < 1.e-10:
                A = np.matrix(np.zeros((2, 1)))
            else:
                ZdZm1 = np.matrix(np.empty((2, 2)))
                ZdZm1[0, 0] = ZdZ[1, 1]
                ZdZm1[0, 1] = -ZdZ[0, 1]
                ZdZm1[1, 0] = -ZdZ[1, 0]
                ZdZm1[1, 1] = ZdZ[0, 0]
                ZdZm1 /= detZdZ
                A = ZdZm1 * (Zd*B)
            ZAa = np.array(Z * A)[:, 0]
            Sres = np.sum((Ba - ZAa)**2 * st2)
            eigVecDist[iev] = Sres
        isortUnfold[ev] = np.argmin(eigVecDist)
#        print eigVecDist
        print '%03d <-> %03d (dist = %.12f)' % (ev, isortUnfold[ev], eigVecDist[isortUnfold[ev]])
        # vecUnfold = eigVecForwardUnfold[:, isortUnfold[ev]]
        # a = np.sum(vecUnfold * st * vecUnfold)**(-1) * np.sum(vecUnfold * st * vec)
        # ergoPlot.plot2D(X, Y, vec.real, ev_xlabel, ev_ylabel)
        # plt.title('ev.real = %d' % ev)
        # ergoPlot.plot2D(X, Y, (a*vecUnfold).real, ev_xlabel, ev_ylabel)
        # plt.title('ev.real unfold = %d' % isortUnfold[ev])
        # if np.abs(eigValForward[ev].imag) > 1.e-7:
        #     ergoPlot.plot2D(X, Y, vec.imag, ev_xlabel, ev_ylabel)
        #     plt.title('ev.imag = %d' % ev)
        #     ergoPlot.plot2D(X, Y, (a*vecUnfold).imag, ev_xlabel, ev_ylabel)
        #     plt.title('ev.imag unfold = %d' % isortUnfold[ev])
    eigValForwardUnfoldSort = eigValForwardUnfold[isortUnfold]

    # Get generator eigenvalues
    print 'Converting to generator eigenvalues...'
    eigValGen = ergoPlot.eig2Generator(eigValForward, tau, eigValForwardUnfoldSort, tauUnfold)
else:
    eigValGen = eigValGenOrig.copy()

# Plot eigenvectors of transfer operator
alpha = 0.01
for ev in np.arange(nevPlot):
    if plotPolar:
        print 'Plotting polar eigenvector %d...' % (ev + 1,)
        #ergoPlot.plot2D(X, Y, eigVecForward[:, ev].real, ev_xlabel, ev_ylabel, alpha)
        ergoPlot.plotEigVecPolar(X, Y, eigVecForward[:, ev], ev_xlabel, ev_ylabel, alpha)
        dstFile = '%s/spectrum/eigvec/eigvecForwardPolar_nev%d_ev%03d%s.%s' \
                  % (cfg.general.plotDir, cfg.spectrum.nev, ev + 1, postfix, ergoPlot.figFormat)
        plt.savefig(dstFile, bbox_inches=ergoPlot.bbox_inches, dpi=ergoPlot.dpi)
        if plotBackward:
            print 'Plotting polar backward eigenvector %d...' % (ev + 1,)
            #ergoPlot.plot2D(X, Y, eigVecForward[:, ev].real, ev_xlabel, ev_ylabel, alpha)
            ergoPlot.plotEigVecPolar(X, Y, eigVecBackward[:, ev], ev_xlabel, ev_ylabel, alpha)
            dstFile = '%s/spectrum/eigvec/eigvecBackwardPolar_nev%d_ev%03d%s.%s' \
                      % (cfg.general.plotDir, cfg.spectrum.nev, ev + 1, postfix, ergoPlot.figFormat)
            plt.savefig(dstFile, bbox_inches=ergoPlot.bbox_inches, dpi=ergoPlot.dpi)

    else:
        print 'Plotting real part of eigenvector %d...' % (ev + 1,)
        #ergoPlot.plot2D(X, Y, eigVecForward[:, ev].real, ev_xlabel, ev_ylabel, alpha)
        ergoPlot.plot2D(X, Y, eigVecForward[:, ev].real, ev_xlabel, ev_ylabel, alpha)
        dstFile = '%s/spectrum/eigvec/eigvecForwardReal_nev%d_ev%03d%s.%s' \
                  % (cfg.general.plotDir, cfg.spectrum.nev, ev + 1, postfix, ergoPlot.figFormat)
        plt.savefig(dstFile, bbox_inches=ergoPlot.bbox_inches, dpi=ergoPlot.dpi)
    
        if plotImag & (eigValForward[ev].imag != 0):
            print 'Plotting imaginary  part of eigenvector %d...' % (ev + 1,)
            #        ergoPlot.plot2D(X, Y, eigVecForward[:, ev].imag, ev_xlabel, ev_ylabel, alpha)
            ergoPlot.plot2D(X, Y, eigVecForward[:, ev].imag, ev_xlabel, ev_ylabel, alpha)
            dstFile = '%s/spectrum/eigvec/eigvecForwardImag_nev%d_ev%03d%s.%s' \
                      % (cfg.general.plotDir, cfg.spectrum.nev, ev + 1, postfix, ergoPlot.figFormat)
            plt.savefig(dstFile, bbox_inches=ergoPlot.bbox_inches, dpi=ergoPlot.dpi)
    
        # Plot eigenvectors of backward operator
        if plotBackward:
            print 'Plotting real part of backward eigenvector %d...' % (ev + 1,)
            ergoPlot.plot2D(X, Y, eigVecBackward[:, ev].real, ev_xlabel, ev_ylabel, alpha)
            dstFile = '%s/spectrum/eigvec/eigvecBackwardReal_nev%d_ev%03d%s.%s' \
                      % (cfg.general.plotDir, cfg.spectrum.nev, ev + 1, postfix, ergoPlot.figFormat)
            plt.savefig(dstFile, bbox_inches=ergoPlot.bbox_inches, dpi=ergoPlot.dpi)
            
            if plotImag & (eigValForward[ev].imag != 0):
                print 'Plotting imaginary  part of backward eigenvector %d...' % (ev + 1,)
                ergoPlot.plot2D(X, Y, eigVecBackward[:, ev].imag, ev_xlabel, ev_ylabel, alpha)
                dstFile = '%s/spectrum/eigvec/eigvecBackwardImag_nev%d_ev%03d%s.%s' \
                          % (cfg.general.plotDir, cfg.spectrum.nev, ev + 1, postfix,
                             ergoPlot.figFormat)
                plt.savefig(dstFile, bbox_inches=ergoPlot.bbox_inches, dpi=ergoPlot.dpi)

            
# Define observables
print 'Reading corrSample and perio'
f = X.flatten()
g = f
corrSamplePath = '../results/%s/' % srcPostfix
corrSamplePostfix = '_%s_%s_%s_%s_mu%04d_eps%04d' \
                    % (fieldsDef[obsIdx0][2], indicesName[obsIdx0][1],
                       fieldsDef[obsIdx1][2], indicesName[obsIdx1][1],
                       np.round(cfg.caseDef.mu * 1000, 1), np.round(cfg.caseDef.eps * 1000, 1))

# Read ccf
corrSample = np.loadtxt('%s/corrSample%s_nSeeds%d_lagMax%dyr.txt' \
                        % (corrSamplePath, corrSamplePostfix, nSeeds, cfg.stat.lagMax))
lags = np.loadtxt('%s/lags%s_nSeeds%d_lagMax%dyr.txt' \
                  % (corrSamplePath, corrSamplePostfix, nSeeds, cfg.stat.lagMax))
powerSample = np.loadtxt('%s/powerSample%s_nSeeds%d_chunk%dyr.txt' \
                         % (corrSamplePath, corrSamplePostfix, nSeeds, cfg.stat.chunkWidth))
freq = np.loadtxt('%s/freq%s_nSeeds%d_chunk%dyr.txt' \
                  % (corrSamplePath, corrSamplePostfix, nSeeds, cfg.stat.chunkWidth))

# Convert to angular frequencies and normalize by covariance
angFreq = freq * 2*np.pi
cfg0 = ((f - (f * statDist).sum()) * statDist * (g - (g * statDist).sum())).sum()
powerSample /= 2 * np.pi * cfg0

# Reconstruct correlation and power spectrum
# Get normalized weights
weights = ergoPlot.getSpectralWeights(f, g, eigVecForward, eigVecBackward,
                                      statDist, skipMean=True)
# Remove components with heigh condition number
weights[eigenCondition > cfg.stat.maxCondition] = 0.
condition = np.empty(eigenCondition.shape, dtype='S1')
condition[:] = 'k'
condition[eigenCondition > cfg.stat.maxCondition - 0.001] = 'w'
(corrRec, compCorrRec) = ergoPlot.spectralRecCorrelation(lags, f, g,
                                                         eigValGen, weights,
                                                         statDist,
                                                         skipMean=True,
                                                         norm=True)
(powerRec, compPowerRec) = ergoPlot.spectralRecPower(angFreq, f, g,
                                                     eigValGen, weights,
                                                     statDist, norm=True)

# Plot correlation reconstruction
ergoPlot.plotRecCorrelation(lags, corrSample, corrRec, plotPositive=True,
                            ylabel=corrLabel, xlabel=xlabelCorr)
plt.savefig('%s/spectrum/reconstruction/%sRec_lag%d_nev%d%s.%s'\
            % (cfg.general.plotDir, corrName, int(cfg.stat.lagMax),
               cfg.spectrum.nev, postfix, ergoPlot.figFormat),
            dpi=ergoPlot.dpi, bbox_inches=ergoPlot.bbox_inches)

# PLot spectrum, powerSampledogram and spectral reconstruction
weights /= cfg0
msize = np.zeros((weights.shape[0]))
msize[weights.real > 0] = np.log10(weights[weights.real > 0].real)
msize[weights.real > 0] = (msize[weights.real > 0] + 8) * 10
# msize[weights.real > 0] = (msize[weights.real > 0] + 6) * 3
msize[msize < 0] = 0.
ergoPlot.plotEigPowerRec(angFreq, eigValGen, powerSample, powerRec,
                         markersize=msize, condition=condition,
                         xlabel=realLabel, ylabel=imagLabel,
                         zlabel=powerLabel,
                         xlim=xlimEig, ylim=ylimEig, zlim=zlimEig,
                         xticks=xticks, yticks=yticks, zticks=zticks)
plt.savefig('%s/spectrum/reconstruction/%sRec_chunk%d_nev%d%s.%s'\
            % (cfg.general.plotDir, powerName, int(cfg.stat.chunkWidth),
               cfg.spectrum.nev, postfix, ergoPlot.figFormat),
            dpi=ergoPlot.dpi, bbox_inches=ergoPlot.bbox_inches)

