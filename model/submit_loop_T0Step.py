import os
import numpy as np

#muRng = np.array([2.3, 2.7, 2.8, 2.825, 2.85, 2.875, 2.9, 3., 3.4])
muRng = np.arange(2.5, 2.7, 0.025)
epsRng = np.array([0.])
T0 = 0.5

prefix = 'zc_1eof'
templateName = '_T0Step%03d_template' % (T0*100,)
designFile = 'design.par'
jobFile = 'submit.sh'
jobOutDir = 'job_output'
templateDir = '%s%s' % (prefix, templateName)

# Get root directory
cwd = os.getcwd()
for eps in epsRng:
    for mu in muRng:
        runDir = '%s_T0Step%03d_mu%04d_eps%04d' \
            % (prefix, T0*100, np.round(mu*1000, 1), np.round(eps*1000, 1))
        # Get back to root directory
        os.chdir(cwd)
        # Remove existing run directory
        os.system('rm -r %s 2> /dev/null' % runDir)
        # Copy template directory to run directory
        os.system('cp -r %s %s' % (templateDir, runDir))
        # Change to run directory
        os.chdir(runDir)
        # Write parameter file anew
        os.system('echo "mu_0= %.4f" > %s' % (mu, designFile))
        os.system('echo "epsn= %.4f" >> %s' % (eps, designFile))
        # Submit job
        os.system('mkdir %s 2> /dev/null' % jobOutDir)
        os.system('sbatch %s' % jobFile)
        

