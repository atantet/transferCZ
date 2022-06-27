import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import itertools
import pylibconfig2
import ergoplot

configFile = '../cfg/transferCZ.cfg'
cfg = pylibconfig2.Config()
cfg.read_file(configFile)

timeScaleConversion = 1. / 12
dim = len(cfg.caseDef.indicesName)

nev = cfg.spectrum.nev
nevPlot = 10
#plotBackward = False
plotBackward = True
plotImag = False
#plotImag = True
xminEigVal = -1.2
yminEigVal = -11.5
ev_xlabel = r'$x_1$'
ev_ylabel = r'$x_2$'

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
(X, Y) = ergoplot.readGrid(gridFile, dim)
coord = (X.flatten(), Y.flatten())

if hasattr(cfg.transfer, "tauDimRng"):
    tauDimRng = np.array(cfg.transfer.tauDimRng)
    nLags = len(tauDimRng)
elif hasattr(cfg.transfer, "stepLag") \
     and hasattr(cfg.transfer, "lag0") \
     and hasattr(cfg.transfer, "nLags"):
    nLags = cfg.transfer.nLags
    tauDimRng = np.empty((nLags,))
    for lag in np.arange(nLags):
        tauDimRng[lag] = cfg.transfer.lag0 + lag*cfg.transfer.stepLag

tauDimRng = np.arange(1., 80., 1.)
nLags = tauDimRng.shape[0]
 
eigenCondition = np.empty((nLags, cfg.spectrum.nev))
eigVal = np.empty((nLags, cfg.spectrum.nev), dtype=complex)
eigValGen = np.empty((nLags, cfg.spectrum.nev), dtype=complex)
for lag in np.arange(nLags):
    tau = tauDimRng[lag]
    postfix = "%s_tau%03d" % (gridPostfix, tau * 1000)
    tauConv = tau * timeScaleConversion
    
    # Define file names
    EigValForwardFile = '%s/eigval/eigValForward_nev%d%s.txt' \
                        % (cfg.general.specDir, cfg.spectrum.nev, postfix)
    EigVecForwardFile = '%s/eigvec/eigVecForward_nev%d%s.txt' \
                        % (cfg.general.specDir, cfg.spectrum.nev, postfix)
    EigValBackwardFile = '%s/eigval/eigValBackward_nev%d%s.txt' \
                        % (cfg.general.specDir, cfg.spectrum.nev, postfix)
    EigVecBackwardFile = '%s/eigvec/eigVecBackward_nev%d%s.txt' \
                        % (cfg.general.specDir, cfg.spectrum.nev, postfix)

    # Read stationary distribution
    statDist = np.loadtxt('%s/transfer/initDist/initDist%s.txt' \
                          % (cfg.general.resDir, postfix))

    # Read transfer operator spectrum from file and create a bi-orthonormal basis
    # of eigenvectors and adjoint eigenvectors:
    print 'Readig spectrum of tau = ', tau
    (eigValForward, eigVecForward, eigValBackward, eigVecBackward) \
        = ergoplot.readSpectrum(cfg.spectrum.nev, EigValForwardFile, EigVecForwardFile,
                                EigValBackwardFile, EigVecBackwardFile, statDist)
    # Save eigenvalues
    eigVal[lag] = eigValForward.copy()
    # Get generator eigenvalues
    eigValGen[lag] = (np.log(np.abs(eigValForward)) + np.angle(eigValForward)*1j) / tauConv
    # Get condition number
    eigenCondition[lag] = ergoplot.getEigenCondition(eigVecForward, eigVecBackward,
                                                     statDist)

    print 'lag ', lag
    print eigVal[lag]
    print eigValGen[lag]
    print eigenCondition[lag]

    # Smoothen
    # if lag > 0:
    #     tmpEigVal = eigVal[lag].copy()
    #     tmpEigValGen = eigValGen[lag].copy()
    #     tmpEigenCondition = eigenCondition[lag].copy()
    #     for ev in np.arange(cfg.spectrum.nev):
    #         idx = np.argmin(np.abs(tmpEigVal[ev] - eigVal[lag-1]))
    #         eigVal[lag, idx] = tmpEigVal[ev]
    #         eigValGen[lag, idx] = tmpEigValGen[ev]
    #         eigenCondition[lag, idx] = tmpEigenCondition[ev]
    #     print 'sort'
    #     print eigVal[lag]
    #     print eigValGen[lag]
    #     print eigenCondition[lag]

    plt.figure()
    plt.scatter(eigValGen[lag].real, eigValGen[lag].imag)
    plt.xlim(xminEigVal / 4, -xminEigVal / 100)
    plt.ylim(yminEigVal, -yminEigVal)



# Plot condition numbers versus absolute value of eigenvalue
# for different lags (and for each eigenvalue)
fig = plt.figure()
ax = fig.add_subplot(111)
markers = itertools.cycle((',', '+', '.', 'o', '*', '^', 'v', '<', '>', '8', 's'))
colors = itertools.cycle(('r', 'g', 'b', 'c', 'm', 'y', 'k'))
msize = 20
maxCond = 2.5
for ev in np.arange(nevPlot):
    color = colors.next()
    plt.plot(eigValGen[:, ev].real, eigenCondition[:, ev],
             linestyle='-', linewidth=1, color=color)
    while markers.next() != ',':
        continue
    for k in np.arange(nLags):
        marker = markers.next()
        plt.scatter(eigValGen[k, ev].real, eigenCondition[k, ev],
                    color=color, marker=marker, s=msize)
ax.set_xlim(xminEigVal / 2, -xminEigVal / 100)
ax.set_ylim(0., maxCond)

