import os
import numpy as np

#muRng = np.arange(2.1, 4.0, 0.2)
#epsRng = np.arange(0., 0.41, 0.05)
#muRng = np.arange(2.725, 2.9, 0.025)
# epsRng = np.arange(0.05, 0.21, 0.05)
#epsRng = np.array([0.01, 0.025])
#epsRng = np.array([0.])
muRng = np.array([2.5, 2.8, 2.9, 3.5])
epsRng = np.array([0.05])
#seedRng = np.arange(0, 10)
seedRng = np.arange(10, 20)

prefix = 'zc_1eof'
templateName = '_template_20151224'
designFile = 'design.par'
jobFile = 'submit.sh'
jobOutDir = 'job_output'
templateDir = '%s%s' % (prefix, templateName)

# Get root directory
cwd = os.getcwd()
for seed in seedRng:
    for eps in epsRng:
        for mu in muRng:
            runDir = '%s_mu%04d_eps%04d_seed%d' \
                % (prefix, np.round(mu*1000, 1), np.round(eps*1000, 1), seed)
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
            os.system('echo "epsn= %.3f" >> %s' % (eps, designFile))
            os.system('echo "seed= %d" >> %s' % (seed, designFile))
            # Submit job
            os.system('mkdir %s 2> /dev/null' % jobOutDir)
            os.system('sbatch %s' % jobFile)
        

