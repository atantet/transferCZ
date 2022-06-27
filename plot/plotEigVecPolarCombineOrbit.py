import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pylibconfig2
from ergopack import ergoplot

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
# Seed postfix
nSeeds = len(cfg.caseDef.seedRng)
seedPostfix = '_seeds' + ''.join(str(s) for s in cfg.caseDef.seedRng)

nev = cfg.spectrum.nev

# evPlot = np.array([0])
# plotForward = True
# plotBackward = False
# cbar_format = '{:.2e}'

evPlot = np.array([1, 2, 3, 4, 5, 6])
plotForward = False
plotBackward = True
cbar_format = '{:.3e}'

if cfg.caseDef.mu >= 2.9:
    xmin = 24.
    xmax = 29.5
    ymin = 80.
    ymax = 180.
else:
    xmin = 25.
    xmax = 27.
    ymin = 100.
    ymax = 130.
plot_orbit = True
# plot_orbit = False

# ampMin = 0.
# ampMax = 0.07
# nlevAmp = 11
ampMin = None
ampMax = None
nlevAmp = None
orbit_color = '0.5'
orbit_size = 48
orbit_lw = 2


def d_formatter(x, pos=None):
    fmt = '' if x % 1 > 1.e-6 else '{:.0f}'.format(x)
    return fmt


xtick_formatter = d_formatter
ytick_formatter = d_formatter

field_h = (1, 'H', 'h', 'm')
field_T = (2, 'SST', 'T', r'$^\circ$C')
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
coord = ergoplot.readGrid(gridFile, dimObs)

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
dstPostfix = '{}{}'.format(gridPostfix, seedPostfix)
dstPostfixTau = "%s_tau%03d" % (dstPostfix, int(tauDim * 1000 + 0.1))
specDir = '%s/spectrum/' % cfg.general.plotDir

# File names
eigValForwardFile = '%s/eigval/eigvalForward_nev%d%s.%s' \
                    % (cfg.general.specDir, nev, dstPostfixTau, fileFormat)
eigVecForwardFile = '%s/eigvec/eigvecForward_nev%d%s.%s' \
    % (cfg.general.specDir, nev, dstPostfixTau, fileFormat)
eigValBackwardFile = '%s/eigval/eigvalBackward_nev%d%s.%s' \
    % (cfg.general.specDir, nev, dstPostfixTau, fileFormat)
eigVecBackwardFile = '%s/eigvec/eigvecBackward_nev%d%s.%s' \
    % (cfg.general.specDir, nev, dstPostfixTau, fileFormat)
maskFile = '%s/transfer/mask/mask%s.%s' \
           % (cfg.general.resDir, dstPostfix, fileFormat)

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


# Read orbit
units = {'T': cfg.units.delta_T, 'h': cfg.units.H}
if plot_orbit:
    dataDir = 'zc_1eof_mu{:04d}_eps0000_seed0'.format(
        int(cfg.caseDef.mu * 1000 + 0.1))
    indicesDir = cfg.general.indicesDir
    xorbit = []
    for idxName, fieldName in zip(
            cfg.caseDef.indicesName, cfg.caseDef.fieldsName):
        filePath = os.path.join(indicesDir, dataDir, idxName + '.txt')
        xorbit.append(np.expand_dims(
            np.loadtxt(filePath)[:, 1] * units[fieldName], axis=1))
    nt = np.min([xo.shape[0] for xo in xorbit])
    xorbit = np.concatenate([xo[:nt] for xo in xorbit], axis=1)
    it0 = int(xorbit.shape[0] * 0.99)
    xorbit = xorbit[it0:]
    if np.var(xorbit[-1] - xorbit[-2]) < 1.e-8 * np.var(xorbit[-1]):
        plot_orbit_scatter = True
    else:
        plot_orbit_scatter = False


def plot_orbit(xorbit):
    if plot_orbit:
        # Add orbit
        ax = fig.gca()
        if plot_orbit_scatter:
            ax.scatter(xorbit[:, 0], xorbit[:, 1], c=orbit_color, s=orbit_size,
                       zorder=3)
        else:
            ax.plot(xorbit[:, 0], xorbit[:, 1], linestyle='-',
                    linewidth=orbit_lw, color=orbit_color, zorder=3)


# Plot eigenvectors of transfer operator
alpha = 0.05
csfilter = 0.5
os.makedirs(os.path.join(cfg.general.plotDir, 'spectrum', 'eigvec'),
            exist_ok=True)
os.makedirs(os.path.join(cfg.general.plotDir, 'spectrum', 'reconstruction'),
            exist_ok=True)
for ev in evPlot:
    cmap = cm.hot_r if ev == 0 else cm.RdBu_r

    if plotForward:
        print('Plotting polar eigenvector {:d}...'.format(ev + 1))
        fig = ergoplot.plotEigVecPolarCombine(
            X, Y, eigVecForward[ev], mask=mask, xlabel=ev_xlabel,
            ylabel=ev_ylabel, alpha=alpha, cmap=cmap, ampMin=ampMin,
            ampMax=ampMax, nlevAmp=nlevAmp, csfilter=csfilter, xmin=xmin,
            xmax=xmax, ymin=ymin, ymax=ymax, xtick_formatter=xtick_formatter,
            ytick_formatter=ytick_formatter, cbar_format=cbar_format)

        plot_orbit(xorbit)

        dstFile = os.path.join(
            specDir, 'eigvec',
            'eigvecForwardPolar_nev{:d}_ev{:03d}{}.{}'.format(
                nev, ev + 1, dstPostfixTau, ergoplot.figFormat))
        fig.savefig(dstFile, bbox_inches=ergoplot.bbox_inches,
                    dpi=ergoplot.dpi)
    if plotBackward:
        print('Plotting polar backward eigenvector {:d}...'.format(ev + 1))
        fig = ergoplot.plotEigVecPolarCombine(
            X, Y, eigVecBackward[ev], mask=mask, xlabel=ev_xlabel,
            ylabel=ev_ylabel, alpha=alpha, cmap=cmap, ampMin=ampMin,
            ampMax=ampMax, nlevAmp=nlevAmp, csfilter=csfilter,
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
            xtick_formatter=xtick_formatter, ytick_formatter=ytick_formatter,
            cbar_format=cbar_format)

        plot_orbit(xorbit)

        dstFile = os.path.join(
            specDir, 'eigvec',
            'eigvecBackwardPolar_nev{:d}_ev{:03d}{}.{}'.format(
                nev, ev + 1, dstPostfixTau, ergoplot.figFormat))
        fig.savefig(dstFile, bbox_inches=ergoplot.bbox_inches,
                    dpi=ergoplot.dpi)

plt.show(block=False)
