import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pylibconfig2
from ergopack import ergoplot
sys.path.append('../cfg/')
from coupledRO import *
from ergopack.ergoInt import *

#figFormat = ergoplot.figFormat
figFormat = 'png'

direct = 'Forward'
#direct = 'Backward'

configFile = '../cfg/transferCZ.cfg'
cfg = pylibconfig2.Config()
cfg.read_file(configFile)
os.system('mkdir %s/spectrum/ 2> /dev/null' % cfg.general.plotDir)
os.system('mkdir %s/spectrum/eigval/ 2> /dev/null' % cfg.general.plotDir)
os.system('mkdir %s/spectrum/eigvec/ 2> /dev/null' % cfg.general.plotDir)

compNames = (r'$T_E~({}^\circ C)$', r'$h_W~(m)$')
realLabel = r'$\Re(\lambda_k)$'
imagLabel = r'$\Im(\lambda_k)$'
xmineigVal = -cfg.stat.rateMax
ymineigVal = -cfg.stat.angFreqMax
xlimEig = [xmineigVal, -xmineigVal/100]
ylimEig = [ymineigVal, -ymineigVal]
(ev_xlabel, ev_ylabel) = compNames
#xlimEV = [-20., 20.]
#ylimEV = [-20., 20.]


# Transition lag
if (hasattr(cfg.stat, 'tauDimPlot')):
    tauDim = cfg.stat.tauDimPlot
else:
    tauDim = cfg.transfer.tauRng[0]
timeScaleConversion = 1. / 12
tau = tauDim * timeScaleConversion

L = cfg.simulation.LCut + cfg.simulation.spinup
spinup = cfg.simulation.spinup
printStepNum = int(cfg.simulation.printStep / cfg.simulation.dt + 0.1)
caseName = cfg.model.caseName
dim = cfg.model.dim
obsComp = np.array(cfg.observable.components)
dimObs = obsComp.shape[0]
nProc = ''
if (hasattr(cfg.sprinkle, 'nProc')):
    nProc = '_nProc' + str(cfg.sprinkle.nProc)

N = np.prod(np.array(cfg.grid.nx))
gridPostfix = ""
for d in np.arange(dimObs):
    if (hasattr(cfg.grid, 'gridLimitsLow')
        & hasattr(cfg.grid, 'gridLimitsHigh')):
        gridPostfix = "%s_n%dl%dh%d" % (gridPostfix, cfg.grid.nx[d],
                                        cfg.grid.gridLimitsLow[d],
                                        cfg.grid.gridLimitsHigh[d])
    else:
        gridPostfix = "%s_n%dl%dh%d" % (gridPostfix, cfg.grid.nx[d],
                                        cfg.sprinkle.minInitState[d],
                                        cfg.sprinkle.maxInitState[d])
gridFile = '%s/grid/grid_%s%s.txt' \
           % (cfg.general.resDir, caseName, gridPostfix)
coord = ergoplot.readGrid(gridFile, dimObs)
xyz = np.zeros((coord[0].shape[0], 3))
xyz[:, 0] = coord[0]
xyz[:, 1] = coord[1]
diag = diagnose(xyz, p)
X, Y = np.meshgrid(diag['TE'], diag['hW'])
X, Y = X.T, Y.T

# Read the fixed point
dim = cfg.model.dim
fileFormat = cfg.general.fileFormat
initContRng = [0.]
contStepRng = [0.01]
initCont = initContRng[0]
contStep = contStepRng[0]
contAbs = np.sqrt(contStep*contStep)
sign = contStep / contAbs
exp = np.log10(contAbs)
mantis = sign * np.exp(np.log(contAbs) / exp)
contPostfix = "_cont%04d_contStep%de%d" \
    % (int(initCont * 1000 + 0.1), int(mantis*1.01), (int(exp*1.01)))
dstPostfixFP = "%s_eta2%04d_r%04d_gamma%04d%s" \
    % (srcPostfix, int(p["eta2"] * 1000 + 0.1),
       int(p["r"] * 1000 + 0.1),
       int(p["gamma"] * 1000 + 0.1), contPostfix)
fpFileName = '%s/fpState/fpState%s.%s' \
    % (contDir, dstPostfixFP, fileFormat)
stateFP = readFile(fpFileName).reshape(-1, dim+1)
fp = stateFP[:, :dim]
contRngFP = stateFP[:, dim]
dstPostfixPO = "%s_dt%d" \
    % (dstPostfixFP, int(np.round(-np.log10(cfg.simulation.dt)) + 0.1))
poFileName = '%s/poState/poState%s.%s' % (contDir, dstPostfixPO, fileFormat)
statePO = readFile(poFileName)
statePO = statePO.reshape(-1, dim+2)
po = statePO[:, :dim]
TRng = statePO[:, dim+1]
contRngPO = statePO[:, dim]

alpha = 0.0
ss = 4
os.system('mkdir %s/spectrum/eigvec 2> /dev/null' % cfg.general.plotDir)
ampMin = 0.
ampMax = None
# nlevAmp = int((ampMax - ampMin) * 100 + 0.1) + 1
nlevAmp = 11
# csfilter = 0.
csfilter = 1
csfmt = '%1.1f'
lw = 2
msize = 30
# for mu in np.arange(0.6, 1.05, 0.02):
for mu in np.arange(0.6, 1.21, 0.06):
    print('mu: {:.2f}'.format(mu))
    icontFP = np.argmin(np.abs(contRngFP - mu))
    p['mu'] = contRngFP[icontFP]
    diagFP = diagnose(np.expand_dims(fp[icontFP], axis=0), p)
    # diagFP = {list(diagFP.keys())[k]: list(diagFP.values())[k][0]
    #           for k in range(len(diagFP))}

    # Read the periodic orbit
    plotPO = False
    dCont = np.abs(contRngPO - mu)
    icontPO = np.argmin(dCont)
    if np.abs(dCont[icontPO]) < 0.02:
        p['mu'] = contRngPO[icontPO]
        x0 = po[icontPO]
        T = TRng[icontPO]
        ntOrbit = int(np.ceil(T / cfg.simulation.dt) + 0.1)
        dtOrbit = T / ntOrbit
        xt = propagate(x0, fieldRO2D, p, stepRK4, dtOrbit, ntOrbit)
        diagPO = diagnose(xt, p)
        plotPO = True

    # Read grid
    srcPostfixSim = "_%s_mu%03d_eta2%03d_gamma%03d_r%03d_sigmahInf2%03d\
_L%d_spinup%d_dt%d_samp%d" \
    % (caseName, int(mu * 100 + 0.1),
       int(p["eta2"] * 1000 + 0.1), int(p["gamma"] * 1000 + 0.1),
       int(p["r"] * 1000 + 0.1),
       int(p["sigmahInf2"] * 1000 + 0.1), int(L + 0.1),
       int(spinup + 0.1),
       int(-np.round(np.log10(cfg.simulation.dt)) + 0.1),
       printStepNum)
    postfix = "%s_nTraj%d%s" % (srcPostfixSim, cfg.sprinkle.nTraj,
                                gridPostfix)
    postfixTau = "%s_tau%03d" % (postfix, int(tau * 100 + 0.1))

    # Define file names
    eigValFile = '%s/eigval/eigval%s_nev%d%s.%s' \
        % (cfg.general.specDir, direct, cfg.spectrum.nev, postfixTau,
           cfg.general.fileFormat)
    eigVecFile = '%s/eigvec/eigvec%s_nev%d%s.%s' \
        % (cfg.general.specDir, direct, cfg.spectrum.nev, postfixTau,
           cfg.general.fileFormat)
    maskFile = '%s/transfer/mask/mask%s.%s' \
        % (cfg.general.resDir, postfix, cfg.general.fileFormat)

    # Read mask
    if maskFile is not None:
        if cfg.general.fileFormat == 'bin':
            mask = np.fromfile(maskFile, np.int32)
        else:
            mask = np.loadtxt(maskFile, np.int32)
    else:
        mask = np.arange(N)
    NFilled = np.max(mask[mask < N]) + 1

    # Read transfer operator spectrum from file and create a bi-orthonormal basis
    # of eigenvectors and backward eigenvectors:
    if direct == 'Forward':
        (eigVal, eigVec) \
            = ergoplot.readSpectrum(eigValForwardFile=eigValFile,
                                    eigVecForwardFile=eigVecFile)
    if direct == 'Backward':
        (eigVal, eigVec) \
            = ergoplot.readSpectrum(eigValBackwardFile=eigValFile,
                                    eigVecBackwardFile=eigVecFile)
    eigVec = eigVec.T
    nev = eigVal.shape[0]
    eigValGen = np.log(eigVal) / (tau * pdim['tadim2year'])

    # Sort
    isort = np.argsort(-np.abs(eigVal))
    eigVal = eigVal[isort]
    eigVec = eigVec[:, isort]

    # Plot eigenvectors of transfer operator
    for ev in np.arange(1, cfg.spectrum.nEigVecPlot):
        print('\tPlotting eigenvector {:d}...'.format(ev + 1))
        cmap = cm.hot_r if ev == 0 else cm.RdBu_r
        v0 = eigVec[:, 0] / eigVec[:, 0].sum()
        vec = eigVec[:, ev] # / v0
        # vec /= np.max(np.abs(vec))
        ergoplot.plotEigVecPolarCombine(X, Y, vec, mask,
                                        xlabel=ev_xlabel, ylabel=ev_ylabel,
                                        alpha=alpha, cmap=cmap, csfmt=csfmt,
                                        ampMin=ampMin, ampMax=ampMax, nlevAmp=nlevAmp)
        if plotPO:
            plt.plot(diagPO['TE'], diagPO['hW'], linestyle='-',
                     linewidth=lw, color='k')
        plt.scatter(diagFP['TE'], diagFP['hW'], s=msize,
                    c='k', edgecolor='face', marker='o')


        #, xlim=xlimEV, ylim=ylimEV)
        dstFile = '%s/spectrum/eigvec/eigvec%sPolar_ev%03d%s.%s' \
            % (cfg.general.plotDir, direct, ev + 1, postfixTau, figFormat)
        plt.savefig(dstFile, bbox_inches=ergoplot.bbox_inches, dpi=ergoplot.dpi)
