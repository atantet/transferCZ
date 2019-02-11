import numpy as np
import matplotlib.pyplot as plt

idx0 = 2
idx1 = 1
seed = 0
eps = 0.05

mu = 3.5

caseDir = 'zc_1eof_mu%d_eps%04d_seed%d' % (int(mu * 1000 + 0.1), int(eps * 1000 + 0.1), seed)

(ts0, ts1) = (np.loadtxt('%s/nino3.txt' % caseDir)[:, idx0], np.loadtxt('%s/nino4.txt' % caseDir)[:, idx1])
imax = np.min((ts0.shape[0], ts1.shape[0]))
(ts0, ts1) = (ts0[:imax], ts1[:imax])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ts0, ts1)
ax.set_xlabel('T')
ax.set_ylabel('H')
ax.set_title(caseDir)

