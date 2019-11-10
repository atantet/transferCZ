import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

seed = 0
indicesName =  ["nino3", "nino4"];
fieldsName  = ["T", "h"];
fieldIdx = {'T': 2, 'h': 1}

# mu = '2500'
# mu = '2800'
mu = '2900'
# mu = '3500'

# eps = '0050'
eps = '0000'


rootDir = os.path.join('data', 'observables')
dirName = 'zc_1eof_mu{}_eps{}_seed{:d}'.format(mu, eps, seed)

a = []
for idxName, fieldName in zip(indicesName, fieldsName):
    filePath = os.path.join(rootDir, dirName, idxName + '.txt')
    a.append(np.expand_dims(np.loadtxt(filePath)[:, fieldIdx[fieldName]], axis=1))
a = np.concatenate(a, axis=1)

# Save 
dstFilePath = 'zc_1eof_mu{}_eps{}_seed{:d}.txt'.format(mu, eps, seed)
np.savetxt(dstFilePath, a)

s = 5000
plt.figure()
plt.plot(a[-s:, 0], a[-s:, 1])
plt.show(block=False)
