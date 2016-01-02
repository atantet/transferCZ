import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
import atmath

# Dimensions
L = 1.5e7
c0 = 2
timeDim = L / c0 / (60. * 60. * 24)
H = 200
tau_0 = 1.0922666667e-2
delta_T = 1.
sampFreq = 0.35 / 0.06 * 12 # (in year^-1)
timeDim = L / c0 / (60. * 60. * 24)
spinup = int(sampFreq * 10)

indicesDir = '../observables/'
prefix = 'zc_1eof'
timeFreq = 0.35 / 0.060 # integration frequency in months
timeScaleConversion = 1. / 12

field_h = (1, 'thermocline depth', 'h', 'm', H, r'$h^2$')
field_T = (2, 'SST', 'T', r'$^\circ C$', delta_T, r'$(^\circ C)^2$')

nino3 = ('Eastern', 'nino3')
nino4 = ('Western', 'nino4')

mu = 2.5
#mu = 2.8
#mu = 2.9
#mu = 3.5
eps = 0.05
indicesName = [nino3, nino4]
fieldsDef = [field_T, field_h]

dim = 2
nx = 100
nSTDLow = [5, 4]
nSTDHigh = [3, 4]

tauDimRng = np.array([1.])

nev = 200
plotCCF = False
#plotCCF = True
#xminEigval = -2.2
#yminEigval = -23
xminEigval = -1.2
yminEigval = -11.5

os.system('mkdir spectrum/eigval/figs spectrum/eigvec/figs spectrum/ccf 2> /dev/null')

srcPostfix = "%s_mu%04d_eps%04d" % (prefix, mu*1000, eps*1000)
indicesPath = '%s/%s/' % (indicesDir, srcPostfix)
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
#msize = 48
msize = 32
#            figFormat = 'eps'
figFormat = 'png'
dpi = 300
readMap = False
gridXlim = [coord[0].min(), coord[0].max()]
gridYlim = [coord[1].min(), coord[1].max()]


lag = 0
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


# # Plot spectrum
# print 'Plotting spectrum slowest rate ', -1. / eigvalGen[1].real
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.scatter(eigvalGen[1:].real, eigvalGen[1:].imag,
#            c='k', s=msize, marker='o')
# ax.scatter(eigvalGen[0].real, eigvalGen[0].imag,
#            c='r', s=msize, marker='o')
# ax.set_xlabel(r'$\Re(\hat{\lambda}^\mathcal{R}_k)$', fontsize=fs_latex)
# ax.set_ylabel(r'$\Im(\hat{\lambda}^\mathcal{R}_k)$', fontsize=fs_latex)
# plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
# plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
# #ax.set_title('%d-time-step spectrum for %s\nSlowest time-scale: %.1f' \
#     #    % (tau, srcPostfix, -1. / rate[0]))
# ax.set_xlim(xminEigval, -xminEigval / 100)
# ax.set_ylim(yminEigval, -yminEigval)
# xlim = ax.get_xlim()
# ylim = ax.get_ylim()
# plt.plot([xlim[0], xlim[1]], [maxImagRes, maxImagRes], '--k')
# plt.plot([xlim[0], xlim[1]], [-maxImagRes, -maxImagRes], '--k')
# plt.text(xlim[1] - (xlim[1] - xlim[0])*0.21,
#          ylim[0] + (ylim[1] - ylim[0])*0.04,
#          r'$\mu = %.2f$' % (mu,),
#          fontsize=fs_latex)
# fig.savefig('spectrum/eigval/figs/eigval_nev%d%s.%s' \
#             % (nev, postfix, figFormat), bbox_inches='tight', dpi=dpi)
        
    
print 'Reading sim'
indexDef = nino3
fieldDef = field_T
indexData = np.loadtxt('%s/%s.txt' % (indicesPath, indexDef[1]))
timeND = indexData[:, 0]
timeFull = timeND * timeDim / 365
indexFull = indexData[:, fieldDef[0]] * fieldDef[4]

# Remove spinup
sim = indexFull[spinup:]
time = timeFull[spinup:]
nt = time.shape[0]

#nt = simMax / 5
nRAVG = 1
window = np.hamming(nt)
# Get nearest larger power of 2
if np.log2(nt) != int(np.log2(nt)):
    nfft = 2**(int(np.log2(nt)) + 1)
else:
    nfft = nt
# Get frequencies and shift zero frequency to center
freq = np.fft.fftfreq(nfft, d=1./sampFreq)
freq = np.fft.fftshift(freq)

# Apply window
ts = sim[:nt] 
ts -= ts.mean(0)
tsWindowed = ts * window
# Fourier transform and shift zero frequency to center
fts = np.fft.fft(tsWindowed, nfft, 0)
fts = np.fft.fftshift(fts)
# Get periodogram
perio = np.abs(fts)**2 / nt / nfft
perio /= perio.sum()

# Plot
print 'Plotting periodogram...'
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#ax.plot(freq[nfft/2+1:] * 2*np.pi, perio[nfft/2+1:], 'k.', markersize=5)
ax.plot(freq[nfft/2+1:] * 2*np.pi, perio[nfft/2+1:], 'k-')
#    ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0, 4 * 2*np.pi)
zmin = 1.e-12
zmax = 1.e-3
ax.set_ylim(zmin, zmax)
#   ax.set_xlabel(r'years$^{-1}$', fontsize=fs_latex)
#   ax.set_ylabel(fieldDef[5], fontsize=fs_latex)
ax.set_xlabel(r'$\omega$', fontsize=fs_latex)
ax.set_ylabel(r'$\hat{S}_{x,x}(\omega)$', fontsize=fs_latex)
plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
fig.savefig('spectrum/ccf/perio_xx_nev%d%s.%s' % (nev, postfix, figFormat),
            bbox_inches='tight', dpi=dpi)


# PLot 3d
#zmin = 1.e-8
#zmax = 1.e-2
zmin = -8
zmax = -1
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
iEigVal = (eigvalGen.real >= xminEigval) \
          & (eigvalGen.real <= -xminEigval / 100) \
          & (eigvalGen.imag >= yminEigval) \
          & (eigvalGen.imag <= -yminEigval)
ax.scatter(eigvalGen[iEigVal][1:].real,
           eigvalGen[iEigVal][1:].imag,
           np.ones((eigvalGen[iEigVal][1:].shape[0],)) * zmin,
           c='b', s=msize, marker='o', depthshade=False)
ax.set_xlim3d(xminEigval, -xminEigval / 100)
ax.set_ylim3d(yminEigval, -yminEigval)
ax.set_zlim3d(zmin, zmax)
ax.scatter(eigvalGen[0].real, eigvalGen[0].imag, zmin,
           c='r', s=msize, marker='o')
angFreq = freq * 2*np.pi
iangFreq = (angFreq >= yminEigval) & (angFreq <= -yminEigval)
perioPlot = np.log10(perio)
#perioPlot = perio
perioPlot[perioPlot < zmin] = zmin
ax.plot(np.zeros((iangFreq.sum(),)), angFreq[iangFreq],
        perioPlot[iangFreq], 'k-')
eigvalPlot = eigvalGen[iEigVal]
iFirst = np.abs(eigvalPlot.imag / eigvalPlot.real) \
         > np.abs(eigvalPlot[1].imag / eigvalPlot[1].real) / 10
nPlot = 2*7
if iEigVal.sum() < nPlot:
    nPlot = iEigVal.sum()
eigvalPlottmp = np.empty((nPlot,), dtype=complex)
eigvalPlottmp[:np.min([iFirst.sum(), nPlot])] = eigvalPlot[iFirst][:nPlot]
if np.sum(iFirst) < nPlot:
    eigvalPlottmp[iFirst.sum():] = eigvalPlot[:nPlot-iFirst.sum()]
eigvalPlot = eigvalPlottmp
for k in np.arange(eigvalPlot.shape[0]):
    ax.plot([eigvalPlot[k].real, 0.],
            [eigvalPlot[k].imag, eigvalPlot[k].imag],
            [zmin, zmin],
            '--b')
    ax.plot([0., 0.],
            [eigvalPlot[k].imag, eigvalPlot[k].imag],
            [zmin,
             perioPlot[np.argmin(np.abs(angFreq - eigvalPlot[k].imag))]],
            '--b')
if np.sum(~iFirst) > 1:
    ax.plot([eigvalGen[iEigVal][~iFirst][1].real, 0.],
            [eigvalGen[iEigVal][~iFirst][1].imag,
             eigvalGen[iEigVal][~iFirst][1].imag],
            [zmin, zmin], '--g')
    ax.plot([0, 0.],
            [eigvalGen[iEigVal][~iFirst][1].imag,
             eigvalGen[iEigVal][~iFirst][1].imag],
            [zmin, perioPlot[np.argmin(np.abs(angFreq - eigvalGen[iEigVal][~iFirst][1].imag))]], '--g')
#ax.set_xlabel(r'$\Re(\hat{\lambda}^\mathcal{R}_k)$', fontsize=fs_latex)
#ax.set_ylabel(r'$\Im(\hat{\lambda}^\mathcal{R}_k)$', fontsize=fs_latex)
ax.set_xlabel('\n' + r'$\Re(\bar{\lambda}_k)$', fontsize=fs_default,
              linespacing=1.5)
ax.set_ylabel('\n' + r'$\Im(\bar{\lambda}_k)$', fontsize=fs_default,
              linespacing=1)
ax.set_zlabel(r'$\hat{S}_{x, x}(\omega)$', fontsize=fs_default)
ax.set_xticks(np.arange(0., xminEigval, -0.2))[::-1]
zticks = np.arange(zmin, zmax, 2)
ax.set_zticks(zticks)
zticklabels = []
for k in np.arange(zticks.shape[0]):
    zticklabels.append(r'$10^{%d}$' % zticks[k])
ax.set_zticklabels(zticklabels)
ax.view_init(30, -150)
#plt.setp(ax.get_xticklabels(), fontsize=fs_xticklabels)
#plt.setp(ax.get_yticklabels(), fontsize=fs_yticklabels)
#ax.set_title('%d-time-step spectrum for %s\nSlowest time-scale: %.1f' \
        #    % (tau, srcPostfix, -1. / rate[0]))
#ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.))
#ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.))
#ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.))
#ax.grid(False)
fig.savefig('spectrum/ccf/eigPerio%s.%s' \
            % (postfix, figFormat), bbox_inches='tight', dpi=dpi)
