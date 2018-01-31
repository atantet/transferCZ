include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <ergoPack/transferOperator.hpp>
#include <ergoPack/transferSpectrum.hpp>
#include "../cfg/readConfig.hpp"


/** \file spectrumCZ.cpp
 *  \brief Get spectrum of transfer operators for the Cane-Zebiak model.
 *   
 * Get spectrum of transfer operators for the Cane-Zebiak model.
 * For configuration, the file cfg/transferCZ.cfg is parsed using libconfig C++ library.
 * The transition matrices are then read from matrix files in compressed format.
 * The Eigen problem is then defined and solved using ARPACK++.
 * Finally, the results are written to file.
 */


// Configuration
char configFileName[256];
// General
char indicesDir[256];
char resDir[256];
char specDir[256];
char plotDir[256];
// Units
double L;                   //!< Horizontal unit of length
double c0;                  //!< Ocean wave celerity
double timeDim;             //!< Time dimension in days
double sampFreq;            //!< Sampling frequency in months
double H;                   //!< Vertical unit of length
double tau_0;               //!< Wind-stress unit
double delta_T;             //!< Temperature unit
// CaseDef
char prefix[256];
char simType[256];
std::string indicesName[DIM];
fieldDef field_h;   //!< Thermocline depth record definition
fieldDef field_T;   //!< SST record definition
fieldDef fields[DIM];
double mu;
double eps;
gsl_vector_uint *seedRng;
size_t nSeeds;
char obsName[256];
char srcPostfix[256];
// Simulation
double dt;
double spinupMonth;
int dimSeries;
// Grid
bool readGridMem;
size_t N;
gsl_vector_uint *nx;
gsl_vector *nSTDLow;
gsl_vector *nSTDHigh;
char gridPostfix[256];
char gridCFG[256];
char gridFileName[256];
// Transfer
int nLags;
gsl_vector *tauDimRng;
bool stationary;                //!< Whether the problem is stationary or not
// Spectrum
int nev;                        //!< Number of eigenvalues to look for
configAR config;                //!< Configuration data for the eigen problem
bool getForwardEigenvectors;    //!< Whether to get forward eigenvectors
bool getBackwardEigenvectors;   //!< Whether to get backward eigenvectors
bool makeBiorthonormal;         //!< Whether to make eigenvectors biorthonormal


// Main program
int main(int argc, char * argv[])
{
  // Read configuration file
  if (argc < 2)
    {
      std::cout << "Enter path to configuration file:" << std::endl;
      std::cin >> configFileName;
    }
  else
    {
      strcpy(configFileName, argv[1]);
    }
  try
    {
      readConfig(configFileName);
    }
  catch (...)
    {
      std::cerr << "Error reading configuration file" << std::endl;
      return(EXIT_FAILURE);
    }

  
  // Declarations
  // Lags
  double tauDim;

  // Names and Files
  char postfix[256];
  char forwardTransitionFileName[256], initDistFileName[256];

  char EigValForwardFileName[256], EigVecForwardFileName[256],
    EigValBackwardFileName[256], EigVecBackwardFileName[256];

  // Transition Matrix
  transferOperator *transferOp;
  transferSpectrum *transferSpec;
  size_t ev;

  ev = 0;
  // Remove first eigenvalue

  ev++;
  
  // Find first (pair of) eigenvalue(s) increasing the lag
  lag = 0;
  while ((lag < nLags))
    {
      // Read transfer operator

      // Get spectrum
      tauDim = gsl_vector_get(tauDimRng, lag);
      std::cout << "Getting spectrum of transfer operator for lag " << tauDim << std::endl;

      // Get file names
      sprintf(postfix, "%s_tau%03d", gridPostfix, (int) (tauDim * 1000));
      sprintf(forwardTransitionFileName, \
	      "%s/transfer/forwardTransition/forwardTransition%s.coo", resDir, postfix);
      sprintf(initDistFileName, "%s/transfer/initDist/initDist%s.txt",
	      resDir, postfix);
      sprintf(EigValForwardFileName, "%s/eigval/eigvalForward_nev%d%s.txt",
	      specDir, nev, postfix);
      sprintf(EigVecForwardFileName, "%s/eigvec/eigvecForward_nev%d%s.txt",
	      specDir, nev, postfix);
      sprintf(EigValBackwardFileName, "%s/eigval/eigvalBackward_nev%d%s.txt",
	      specDir, nev, postfix);
      sprintf(EigVecBackwardFileName, "%s/eigvec/eigvecBackward_nev%d%s.txt",
	      specDir, nev, postfix);

      
      lag++;
    }

  // Get eigenvalue from previous operator to get imaginary part of generator

  // Save generator eigenvalue

  // Once eigenvalue with largest real part, remove it from operator

  // Find other eigenvalues one by one
  ev = 
  for (size_t ev = 0; ev < (size_t) nev; ev++)
    {
      // Get eigenvalue

      // Get condition number
      
      // Decrease lag until condition number is small enough,
      // removing the component of each recorded eigenvalue

      // Get eigenvalue from previous operator to get imaginary part of generator

      // Save generator eigenvalue

    }

      
  for (size_t lag = 0; lag < (size_t) nLags; lag++)
    {
      tauDim = gsl_vector_get(tauDimRng, lag);
      std::cout << "Getting spectrum of transfer operator for lag " << tauDim << std::endl;

      // Get file names
      sprintf(postfix, "%s_tau%03d", gridPostfix, (int) (tauDim * 1000));
      sprintf(forwardTransitionFileName, \
	      "%s/transfer/forwardTransition/forwardTransition%s.coo", resDir, postfix);
      sprintf(initDistFileName, "%s/transfer/initDist/initDist%s.txt",
	      resDir, postfix);
      sprintf(EigValForwardFileName, "%s/eigval/eigvalForward_nev%d%s.txt",
	      specDir, nev, postfix);
      sprintf(EigVecForwardFileName, "%s/eigvec/eigvecForward_nev%d%s.txt",
	      specDir, nev, postfix);
      sprintf(EigValBackwardFileName, "%s/eigval/eigvalBackward_nev%d%s.txt",
	      specDir, nev, postfix);
      sprintf(EigVecBackwardFileName, "%s/eigvec/eigvecBackward_nev%d%s.txt",
	      specDir, nev, postfix);

      // Read transfer operator
      std::cout << "Reading transfer operator..." << std::endl;
      transferOp = new transferOperator(N, stationary);
      transferOp->scanInitDist(initDistFileName);
      transferOp->scanForwardTransition(forwardTransitionFileName);

      // Get spectrum
      try
	{
	  // Solve eigen value problem with default configuration
	  transferSpec = new transferSpectrum(nev, transferOp, config);

	  if (getForwardEigenvectors)
	    {
	      std::cout << "Solving eigen problem for forward transition matrix..." << std::endl;
	      transferSpec->getSpectrumForward();
	      std::cout << "Found " << transferSpec->EigProbForward.ConvergedEigenvalues()
			<< "/" << nev << " eigenvalues." << std::endl;	      
	    }
	  if (getBackwardEigenvectors)
	    {
	      std::cout << "Solving eigen problem for backward transition matrix..." << std::endl;
	      transferSpec->getSpectrumBackward();
	      std::cout << "Found " << transferSpec->EigProbBackward.ConvergedEigenvalues()
			<< "/" << nev << " eigenvalues." << std::endl;
	    }
	  if (makeBiorthonormal)
	    {
	      std::cout << "Making set of forward and backward eigenvectors biorthonormal..."
			<< std::endl;
	      transferSpec->makeBiorthonormal();
	    }
	}
      catch (std::exception &ex)
	{
	  std::cerr << "Error calculating spectrum: " << ex.what() << std::endl;
	  return EXIT_FAILURE;
	}
  
      // Write spectrum 
      try
	{
	  if (getForwardEigenvectors)
	    {
	      std::cout << "Writing forward eigenvalues and eigenvectors..." << std::endl;
	      transferSpec->writeSpectrumForward(EigValForwardFileName,
						 EigVecForwardFileName);
	    }
	  if (getBackwardEigenvectors)
	    {
	      std::cout << "Writing backward eigenvalues and eigenvectors..." << std::endl;
	      transferSpec->writeSpectrumBackward(EigValBackwardFileName,
						  EigVecBackwardFileName);
	    }
	}
      catch (std::exception &ex)
	{
	  std::cerr << "Error writing spectrum: " << ex.what() << std::endl;
	  return EXIT_FAILURE;
	}

      // Free
      delete transferOp;
      delete transferSpec;
    }

  // Free configuration
  freeConfig();
  
  return 0;
}
