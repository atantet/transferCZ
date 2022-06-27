import os
import numpy as np
import matplotlib.pyplot as plt
import pylibconfig2
from ergopack import ergoplot

configFile = '../cfg/transferCZ.cfg'
cfg = pylibconfig2.Config()
cfg.read_file(configFile)
fileFormat = cfg.general.fileFormat

# Time dimension in days
second_to_year = 1 / (60 * 60 * 24 * 365)
time_dim = cfg.units.L / cfg.units.c0 * second_to_year

# Cases
cases = {'Deterministic': 0., 'Stochastic': 0.01}

field_h = (1, 'H', 'h', 'm')
field_T = (2, 'SST', 'T', r'$^\circ$C')
field_u_A = (3, 'wind stress due to coupling', 'u_A', 'm/s')
field_taux = (4, 'external wind-stress', 'taux', 'm/s')

nino3 = ('E', 'nino3')
nino4 = ('W', 'nino4')
indicesName = []
fieldsDef = []
dimObs = len(cfg.caseDef.indicesName)
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

units = {'T': cfg.units.delta_T, 'h': cfg.units.H}
ev_xlabel = '%s (%s)' % (compName0, fieldsDef[0][3])
ev_ylabel = '%s (%s)' % (compName1, fieldsDef[1][3])
tss = {}
times = {}
for case, eps in cases.items():
    dataDir = 'zc_1eof_mu{:04d}_eps{:04d}_seed0'.format(
        int(cfg.caseDef.mu * 1000 + 0.1), int(eps * 1000 + 0.1))
    indicesDir = cfg.general.indicesDir
    xorbit = []
    t = []
    for idxName, fieldName in zip(
            cfg.caseDef.indicesName, cfg.caseDef.fieldsName):
        filePath = os.path.join(indicesDir, dataDir, idxName + '.txt')
        data = np.loadtxt(filePath)
        xorbit.append(np.expand_dims(data[:, 1] * units[fieldName], axis=1))
    nt = np.min([xo.shape[0] for xo in xorbit])
    tss[case] = np.concatenate([xo[:nt] for xo in xorbit], axis=1)
    times[case] = data[:nt, 0] * time_dim

srcPostfix = "%s%s_mu%04d_eps%04d" % (cfg.caseDef.prefix, cfg.caseDef.simType,
                                      np.round(cfg.caseDef.mu * 1000, 1),
                                      np.round(cfg.caseDef.eps * 1000, 1))
obsName = ''
for d in np.arange(dimObs):
    obsName += '_%s_%s' % (fieldsDef[d][2], indicesName[d][1])
dstPostfix = '_%s%s' % (srcPostfix, obsName)
fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
itmax = 10000
lss = ['-', '--']
colors = [c['color'] for c in plt.rcParams['axes.prop_cycle']]
for k, case in enumerate(tss):
    time, ts = times[case], tss[case]
    ax.plot(time[:itmax], ts[:itmax, 0], linestyle=lss[k],
            color=colors[0], label=case)
    ax2.plot(time[:itmax], ts[:itmax, 1], linestyle=lss[k],
             color=colors[1], label=case)
ax.set_xlabel('time (y)')
ax.set_ylabel(ev_xlabel)
ax2.set_ylabel(ev_ylabel)
ax.legend(loc='best')

if cfg.caseDef.mu >= 2.9:
    xticks = np.arange(25., 28.6, 0.5)
    time_slices = [slice(nt - 1000, nt), slice(nt-1000, nt)]
    lws = [2, 2]
else:
    xticks = np.arange(25., 27.1, 0.5)
    time_slices = [slice(1000, nt), slice(nt-2000, nt)]
    lws = [1, 2]

fig = plt.figure()
ax = fig.add_subplot(111)
lss = ['-', '--']
colors = [c['color'] for c in plt.rcParams['axes.prop_cycle']]
# colors = ['k', 'k']
for k, case in enumerate(tss):
    time, ts = times[case], tss[case]
    time_slice = time_slices[k]
    ax.plot(ts[time_slice, 0], ts[time_slice, 1], linestyle=lss[k],
            linewidth=lws[k], color=colors[k], label=case)
ax.set_xticks(xticks)
ax.set_xlabel(ev_xlabel, fontsize=ergoplot.fs_xlabel)
ax.set_ylabel(ev_ylabel, fontsize=ergoplot.fs_ylabel)
plt.setp(ax.get_xticklabels(), fontsize=ergoplot.fs_xticklabels)
plt.setp(ax.get_yticklabels(), fontsize=ergoplot.fs_yticklabels)
ax.legend(loc='upper left', fontsize=ergoplot.fs_legend_labels)
series_dir = os.path.join(cfg.general.plotDir, 'series')
os.makedirs(series_dir, exist_ok=True)
dstFile = os.path.join(
    series_dir, 'series{}.{}'.format(dstPostfix, ergoplot.figFormat))
fig.savefig(dstFile, bbox_inches=ergoplot.bbox_inches,
            dpi=ergoplot.dpi)

plt.show(block=False)
