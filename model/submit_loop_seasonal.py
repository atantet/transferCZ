import os
import numpy as np

muRng = np.array([2.1, 2.5, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95,
                  3., 3.1, 3.3, 3.7])
#amu0Rng = np.arange(0.1, 0.51, 0.1)
amu0Rng = np.array([1.5, 2., 3.,])
#epsRng = np.arange(0., 0.11, 0.05)
epsRng = np.array([0.])

prefix = 'zc_1eof_seasonal'
templateName = '_template'
designFile = 'design.par'
jobFile = 'submit.sh'
jobOutDir = 'job_output'
templateDir = '%s%s' % (prefix, templateName)

# Get root directory
cwd = os.getcwd()
for eps in epsRng:
    for amu0 in amu0Rng:
        for mu in muRng:
            runDir = '%s_mu%04d_amu0%04d_eps%04d' \
                % (prefix, np.round(mu*1000, 1),
                   np.round(amu0*1000, 1), np.round(eps*1000, 1))
            # Get back to root directory
            os.chdir(cwd)
            # Remove existing run directory
            os.system('rm -r %s 2> /dev/null' % runDir)
            # Copy template directory to run directory
            os.system('cp -r %s %s' % (templateDir, runDir))
            # Change to run directory
            os.chdir(runDir)
            # Write parameter file anew
            os.system('echo "mu_0= %.3f" > %s' % (mu, designFile))
            os.system('echo "amu0= %.3f" >> %s' % (amu0, designFile))
            os.system('echo "epsn= %.3f" >> %s' % (eps, designFile))
            # Submit job
            os.system('mkdir %s 2> /dev/null' % jobOutDir)
            os.system('sbatch %s' % jobFile)
        

