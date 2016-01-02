import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
import atmath

prefix = 'zc_1eof_'
timeFreq = 0.35 / 0.060 # integration frequency in months
timeScaleConversion = 1. / 12

field_h = (1, 'thermocline depth', 'h', 'm')
field_T = (2, 'SST', 'T', r'$^\circ C$')
field_u_A = (3, 'wind stress due to coupling', 'u_A', 'm/s')
field_taux = (4, 'external wind-stress', 'taux', 'm/s')

nino3 = ('Eastern', 'nino3')
nino4 = ('Western', 'nino4')

mu = 3.5
eps = 0.05
indicesName = [nino3, nino4]
fieldsDef = [field_T, field_h]

dim = 2
nx = 100
nSTDLow = [5, 4]
nSTDHigh = [3, 4]

tauDimRng = np.array([1.])

nev = 200
#nevPlot = 0
nevPlot = 4
#plotAdjoint = False
plotAdjoint = True
plotCCF = False
#plotCCF = True
xminEigval = -2.2
yminEigval = -23

os.system('mkdir spectrum/eigval/figs spectrum/eigvec/figs spectrum/ccf 2> /dev/null')

srcPostfix = "%smu%04d_eps%04d" % (prefix, mu*1000, eps*1000)
obsName = ''
N = 1
for d in np.arange(dim):
    obsName += '_%s_%s' % (fieldsDef[d][2], indicesName[d][1])
    N *= nx
    gridPostfix = "n%dl%dh%d" % (nx, nSTDLow[d], nSTDHigh[d])
cpyBuffer = gridPostfix
gridPostfix = '_%s%s_%s' % (srcPostfix, obsName, cpyBuffer)

# Read grid
gridFile = 'grid/grid%s.txt' % gridPostfix
f = open(gridFile, 'r')
bounds = []
coord = []
for k in np.arange(dim):
    bounds.append(np.array(f.readline().split()).astype(float))
    coord.append((bounds[k][1:] + bounds[k][:-1]) / 2)
f.close()
X, Y = np.meshgrid(coord[0], coord[1])

# Plot
levels = 20
fs_default = 'x-large'
fs_latex = 'xx-large'
fs_xlabel = fs_default
fs_ylabel = fs_default
fs_xticklabels = fs_default
fs_yticklabels = fs_default
fs_legend_title = fs_default
fs_legend_labels = fs_default
fs_cbar_label = fs_default
msize = 48
#            figFormat = 'eps'
figFormat = 'png'
dpi = 300
readMap = False
gridXlim = [coord[0].min(), coord[0].max()]
gridYlim = [coord[1].min(), coord[1].max()]


for lag in np.arange(tauDimRng.shape[0]):
    tauDim = tauDimRng[lag]
    tauConv = tauDim * timeScaleConversion
    maxImagRes = np.pi / tauConv
    postfix = "%s_tau%03d" % (gridPostfix, tauDim * 1000)

    print 'Readig spectrum...'
    EigValFile = 'spectrum/eigval/eigval_nev%d%s.txt' % (nev, postfix)
    eigval = np.loadtxt(EigValFile)
    eigval = eigval[:, 0] + eigval[:, 1]*1j
    isort = np.argsort(np.abs(eigval))[::-1]
    eigval = eigval[isort]
    nevSingle = eigval.shape[0]
    if nevPlot > 0:
        EigVecFile = 'spectrum/eigvec/eigvec_nev%d%s.txt' % (nev, postfix)
        eigvec = np.loadtxt(EigVecFile)
        eigvec = eigvec[::2] + eigvec[1::2]*1j
        eigvec = eigvec[isort]

    if plotAdjoint and (nevPlot > 0):
        print 'Readig adjoint spectrum...'
        EigValAdjointFile = 'spectrum/eigval/eigvalAdjoint_nev%d%s.txt' \
                            % (nev, postfix)
        EigVecAdjointFile = 'spectrum/eigvec/eigvecAdjoint_nev%d%s.txt' \
                            % (nev, postfix)
        eigvalAdjoint = np.loadtxt(EigValAdjointFile)
        eigvalAdjoint = eigvalAdjoint[:, 0] + eigvalAdjoint[:, 1]*1j
        eigvecAdjoint = np.loadtxt(EigVecAdjointFile)
        # From the transpose we get the conjugate of the adjoint eigenvectors
        # so we take back the conjugate
        eigvecAdjoint = eigvecAdjoint[::2] - eigvecAdjoint[1::2]*1j
        isort = np.argsort(np.abs(eigvalAdjoint))[::-1]
        eigvalAdjoint = eigvalAdjoint[isort]
        eigvecAdjoint = eigvecAdjoint[isort]
        # eigvecAdjointScale = np.zeros((nevSingle, N), dtype=complex)
        # for k in np.arange(nevSingle):
        #     eigvecAdjointScale[k] = eigvecAdjoint[k] \
        #                         / np.conjugate(np.vdot(eigvecAdjoint[k],
        #                                                eigvec[k]))

    # Get generator eigenvalues
    eigvalGen = np.empty((nevSingle,), dtype=complex)
    ev = 0
    for count in np.arange(eigval.shape[0]):
        eigvalGen[ev] = (np.log(np.abs(eigval[count])) \
                         + np.angle(eigval[count]) * 1j) / tauConv
        ev += 1
        if ev >= nevSingle:
            break
        if eigval[count].imag != 0:
            eigvalGen[ev] = np.conjugate(eigvalGen[ev-1])
            ev +=1
            if ev >= nevSingle:
                break


    # Plot spectrum
    print 'Plotting spectrum slowest rate ', -1. / eigvalGen[1].real
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(eigvalGen[1:].real, eigvalGen[1:].imag,
               c='k', s=msize, marker='o')
    ax.scatter(eigvalGen[0].real, eigvalGen[0].imag,
               c='r', s=msize, marker='o')
    ax.set_xlabel(r'$\Re(\hat{\lambda}^\mathcal{R}_k)$', fontsize=fs_latex)
    ax.set_ylabel(r'$\Im(\hat{\lambda}^\mathcal{R}_k)$', fontsize=fs_latex)
    plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
    plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
    #ax.set_title('%d-time-step spectrum for %s\nSlowest time-scale: %.1f' \
        #    % (tau, srcPostfix, -1. / rate[0]))
    ax.set_xlim(xminEigval, -xminEigval / 100)
    ax.set_ylim(yminEigval, -yminEigval)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    plt.plot([xlim[0], xlim[1]], [maxImagRes, maxImagRes], '--k')
    plt.plot([xlim[0], xlim[1]], [-maxImagRes, -maxImagRes], '--k')
    plt.text(xlim[1] - (xlim[1] - xlim[0])*0.21,
             ylim[0] + (ylim[1] - ylim[0])*0.04,
             r'$\mu = %.2f$' % (mu,),
             fontsize=fs_latex)
    fig.savefig('spectrum/eigval/figs/eigval_nev%d%s.%s' \
                % (nev, postfix, figFormat), bbox_inches='tight', dpi=dpi)
        
    
    # Plot eigenvectors of transfer operator
    tol = 0.
    alpha = 0.01
    for k in np.arange(nevPlot):
        if np.abs(eigval[k].real - 1) < 1.e-3:
            # Plot invariant measure
            print 'Plotting stationary density...'
            statDen = eigvec[0].real
            statDen /= statDen.sum()
            v2Real = statDen.copy()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            vecAlpha = v2Real[v2Real != 0]
            vmax = np.sort(vecAlpha)[(1. - alpha) * vecAlpha.shape[0]]
            v2Real[v2Real > vmax] = vmax
            h = ax.contourf(X, Y, v2Real.reshape(nx, nx), levels,
                            cmap=cm.hot_r, vmax=vmax, vmin=0)
            ax.set_xlim(gridXlim)
            ax.set_ylim(gridYlim)
            cbar = plt.colorbar(h)
            d = 0
            ax.set_xlabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                           fieldsDef[d][3]), fontsize=fs_latex)
            d = 1
            ax.set_ylabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                           fieldsDef[d][3]), fontsize=fs_latex)
            # ax.set_title("Approximation of the invariant measure",
            #              fontsize=fs_default)
            plt.setp(cbar.ax.get_yticklabels(), fontsize=fs_yticklabels)
            plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
            plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
            fig.savefig('spectrum/eigvec/figs/eigvecReal_nev%d_ev%03d%s.%s' \
                        % (nev, 1, postfix, figFormat),
                        bbox_inches='tight', dpi=dpi)
        else:
            print 'Plotting real part of eigenvector %d...' % (k+1,)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            v2Real = eigvec[k].real
            v2Real[statDen != 0] /= statDen[statDen != 0]
            vecAlpha = v2Real[v2Real != 0]
            vmax = np.sort(np.abs(vecAlpha))[(1. - 2*alpha) \
                                             * vecAlpha.shape[0]]
            v2Real[v2Real > vmax] = vmax
            v2Real[v2Real < -vmax] = -vmax
            h = ax.contourf(X, Y, v2Real.reshape(nx, nx), levels,
                            cmap=cm.RdBu_r, vmin=-vmax, vmax=vmax)
            ax.set_xlim(gridXlim)
            ax.set_ylim(gridYlim)
            plt.colorbar(h)
            d = 0
            ax.set_xlabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                           fieldsDef[d][3]), fontsize=fs_latex)
            d = 1
            ax.set_ylabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                           fieldsDef[d][3]), fontsize=fs_latex)
            # ax.set_title("Real part of the eigenvector %d" % (k+1,),
            #              fontsize=fs_default)
            plt.setp(cbar.ax.get_yticklabels(), fontsize=fs_yticklabels)
            plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
            plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
            fig.savefig('spectrum/eigvec/figs/eigvecReal_nev%d_ev%03d%s.%s' \
                        % (nev, k+1, postfix, figFormat),
                        bbox_inches='tight', dpi=dpi)

            if eigval[k].imag != 0:
                print 'Plotting imaginary  part of eigenvector %d...' % (k+1,)
                fig = plt.figure()
                ax = fig.add_subplot(111)
                v2Imag = eigvec[k].imag
                v2Imag[statDen != 0] /= statDen[statDen != 0]
                vecAlpha = v2Imag[v2Imag != 0]
                vmax = np.sort(np.abs(vecAlpha))[(1. - 2*alpha) \
                                                 * vecAlpha.shape[0]]
                v2Imag[v2Imag > vmax] = vmax
                v2Imag[v2Imag < -vmax] = -vmax
                h = ax.contourf(X, Y, v2Imag.reshape(nx, nx), levels,
                                cmap=cm.RdBu_r, vmin=-vmax, vmax=vmax)
                ax.set_xlim(gridXlim)
                ax.set_ylim(gridYlim)
                plt.colorbar(h)
                d = 0
                ax.set_xlabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                               fieldsDef[d][3]), fontsize=fs_latex)
                d = 1
                ax.set_ylabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                               fieldsDef[d][3]), fontsize=fs_latex)
                # ax.set_title("Imaginary part of the eigenvector %d" % (k+1,),
                #              fontsize=fs_default)
                plt.setp(cbar.ax.get_yticklabels(), fontsize=fs_yticklabels)
                plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
                plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
                fig.savefig('spectrum/eigvec/figs/eigvecImag_nev%d_ev%03d%s.%s' % (nev, k, postfix, figFormat), bbox_inches='tight', dpi=dpi)

 
    # Plot eigenvectors of Koopman operator
    if plotAdjoint:
        eigvecAdjointScale = np.zeros((nevSingle, N), dtype=complex)
        for k in np.arange(nevPlot):
            eigvecAdjointScale[k] = eigvecAdjoint[k] \
                                    / np.conjugate(np.vdot(eigvecAdjoint[k],
                                                           eigvec[k]))
            if np.abs(eigval[k].real - 1) < 1.e-3:
                # Plot invariant measure
                print 'Plotting ergodic vector...'
                ergodicVec = np.abs(eigvecAdjoint[0].real)
                ergodicVec /= np.abs(ergodicVec).max()
                fig = plt.figure()
                ax = fig.add_subplot(111)
                alpha = 0.01
                h = ax.contourf(X, Y, ergodicVec.reshape(nx, nx), levels,
                                cmap=cm.hot_r)
                ax.set_xlim(gridXlim)
                ax.set_ylim(gridYlim)
                plt.colorbar(h)
                d = 0
                ax.set_xlabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                               fieldsDef[d][3]), fontsize=fs_latex)
                d = 1
                ax.set_ylabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                               fieldsDef[d][3]), fontsize=fs_latex)
                plt.setp(cbar.ax.get_yticklabels(), fontsize=fs_yticklabels)
                plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
                plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
                # ax.set_title("Approximation of the ergodic vector",
                #              fontsize=fs_default)
                fig.savefig('spectrum/eigvec/figs/eigvecAdjointReal_nev%d_ev%03d%s.%s' % (nev, 1, postfix, figFormat), bbox_inches='tight', dpi=dpi)
            else:
                print 'Plotting real part of Koopman eigenvector %d...' \
                    % (k+1,)
                fig = plt.figure()
                ax = fig.add_subplot(111)
                v2Real = eigvecAdjoint[k].real
                vecAlpha = v2Real[v2Real != 0]
                vmax = np.sort(np.abs(vecAlpha))[(1. - 2*alpha) \
                                                 * vecAlpha.shape[0]]
                v2Real[v2Real > vmax] = vmax
                v2Real[v2Real < -vmax] = -vmax
                h = ax.contourf(X, Y, v2Real.reshape(nx, nx), levels,
                                cmap=cm.RdBu_r, vmin=-vmax, vmax=vmax)
                ax.set_xlim(gridXlim)
                ax.set_ylim(gridYlim)
                plt.colorbar(h)
                d = 0
                ax.set_xlabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                               fieldsDef[d][3]), fontsize=fs_latex)
                d = 1
                ax.set_ylabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                               fieldsDef[d][3]), fontsize=fs_latex)
                plt.setp(cbar.ax.get_yticklabels(), fontsize=fs_yticklabels)
                plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
                plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
                # ax.set_title("Real part of the Koopman eigenvector %d" \
                #              % (k+1,),
                #              fontsize=fs_default)
                fig.savefig('spectrum/eigvec/figs/eigvecAdjointReal_nev%d_ev%03d%s.%s' % (nev, k+1, postfix, figFormat), bbox_inches='tight', dpi=dpi)

                if eigval[k].imag != 0:
                    print 'Plotting imaginary  part of Koopman eigenvector %d...' % (k+1,)
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    v2Imag = eigvecAdjoint[k].imag
                    vecAlpha = v2Imag[v2Imag != 0]
                    vmax = np.sort(np.abs(vecAlpha))[(1. - 2*alpha) \
                                                     * vecAlpha.shape[0]]
                    v2Imag[v2Imag > vmax] = vmax
                    v2Imag[v2Imag < -vmax] = -vmax
                    h = ax.contourf(X, Y, v2Imag.reshape(nx, nx), levels,
                                    cmap=cm.RdBu_r, vmin=-vmax, vmax=vmax)
                    ax.set_xlim(gridXlim)
                    ax.set_ylim(gridYlim)
                    plt.colorbar(h)
                    d = 0
                    ax.set_xlabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                                   fieldsDef[d][3]), fontsize=fs_latex)
                    d = 1
                    ax.set_ylabel(r'%s %s (%s)' % (indicesName[d][0], fieldsDef[d][1],
                                                   fieldsDef[d][3]), fontsize=fs_latex)
                    plt.setp(cbar.ax.get_yticklabels(),
                             fontsize=fs_yticklabels)
                    plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
                    plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
                    # ax.set_title("Imaginary part of the Koopman eigenvector %d" % (k+1,), fontsize=fs_default)
                    fig.savefig('spectrum/eigvec/figs/eigvecAdjointImag_nev%d_ev%03d%s.%s' % (nev, k, postfix, figFormat), bbox_inches='tight', dpi=dpi)

if plotCCF:
    # Get ccf
#    lagMax = 1e3
    lagMax = 1e2
    lagMaxSample = lagMax / (dt * sampling)
    lags = np.arange(0, lagMax+dt*sampling, dt*sampling)
    f = X.flatten()
    g = f
    obsIdx0 = 0
    obsIdx1 = 0
    
#    Get sample cross-correlation
    postfix = '_%s_L%d_spinup%d_dt%d_samp%d' \
              % (caseName, L, spinup, -np.log10(dt), sampling)
    print 'Reading sim'
    sim = np.loadtxt('sim/sim%s.txt' % postfix)
    print 'Getting acf'
    simMax = L / (dt * sampling)
    acf = atmath.ccf(sim[:simMax], sim[:simMax],
                     lagMax=lagMaxSample)[lagMaxSample:]
    np.savetxt('sim/acf%s.txt' % postfix, acf)
    
    # Get reconstructed cross correlation
    acfRec = np.zeros((lags.shape[0],), dtype=complex)
    for ev in np.arange(1, nevSingle):
        acfRec += np.abs(np.exp(eigvalGen[ev]*lags) \
                  * (f * statDen * np.conjugate(eigvecAdjoint[ev])).sum() \
                  * (eigvec[ev] * np.conjugate(g)).sum())
    acfRec /= acfRec[0]
    
    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(lags, acf, linewidth=2)
    ax.plot(lags, acfRec, '--', linewidth=2)
    fig.savefig('spectrum/ccf/ccf_xx_nev%d%s.%s' % (nev, postfix, figFormat),
                bbox_inches='tight', dpi=dpi)
