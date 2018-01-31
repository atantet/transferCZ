#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <transferOperator.hpp>
#include <transferSpectrum.hpp>
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
char fileFormat[256];
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
// Spectrum
int nev;                        //!< Number of eigenvalues to look for
configAR config;                //!< Configuration data for the eigen problem
bool getForwardEigenvectors;    //!< Whether to get forward eigenvectors
bool getBackwardEigenvectors;   //!< Whether to get backward eigenvectors
bool makeBiorthonormal;      //!< Whether to make eigenvectors biorthonormal


// Main program
int main(int argc, char * argv[])
{
  // Read configuration file
  if (argc < 2) {
    std::cout << "Enter path to configuration file:" << std::endl;
    std::cin >> configFileName;
  }
  else {
    strcpy(configFileName, argv[1]);
  }
  try {
    readConfig(configFileName);
  }
  catch (...) {
    std::cerr << "Error reading configuration file" << std::endl;
    return(EXIT_FAILURE);
  }

  
  // Declarations
  // Transfer
  double tauDim;
  char dstPostfix[256], dstPostfixTau[256];
  char transitionFileName[256], initDistFileName[256], maskFileName[256];
  transferOperator *transferOp;
  gsl_vector *initDist;
  gsl_vector_uint *mask;

  // Eigenproblem
  char EigValForwardFileName[256], EigVecForwardFileName[256],
    EigValBackwardFileName[256], EigVecBackwardFileName[256];
  transferSpectrum *transferSpec;
  
  sprintf(dstPostfix, "%s", gridPostfix);
  
  // Scan transition matrices and distributions for different lags
  for (size_t lag = 0; lag < (size_t) nLags; lag++) {
    tauDim = gsl_vector_get(tauDimRng, lag);
    std::cout << "Getting spectrum of transfer operator for lag "
	      << tauDim << std::endl;

    // Get file names
    sprintf(dstPostfixTau, "%s_tau%03d", dstPostfix,
	    (int) (tauDim * 1000 + 0.1));
    sprintf(transitionFileName, \
	    "%s/transfer/transition/transition%s.crs%s",
	    resDir, dstPostfixTau, fileFormat);
    sprintf(EigValForwardFileName, "%s/eigval/eigvalForward_nev%d%s.%s",
	    specDir, nev, dstPostfixTau, fileFormat);
    sprintf(EigVecForwardFileName, "%s/eigvec/eigvecForward_nev%d%s.%s",
	    specDir, nev, dstPostfixTau, fileFormat);
    sprintf(EigValBackwardFileName, "%s/eigval/eigvalBackward_nev%d%s.%s",
	    specDir, nev, dstPostfixTau, fileFormat);
    sprintf(EigVecBackwardFileName, "%s/eigvec/eigvecBackward_nev%d%s.%s",
	    specDir, nev, dstPostfixTau, fileFormat);

    // Read transfer operator
    std::cout << "Reading transfer operator..." << std::endl;
    try {
      /** Construct transfer operator without allocating memory
	  to the distributions (only to the mask) ! */
      transferOp = new transferOperator(N);

      // Scan transition matrix (this sets NFilled)
      std::cout << "Scanning transition matrix from "
		<< transitionFileName << std::endl;
      transferOp->scanTransition(transitionFileName, fileFormat);

      // Scan mask for the first lag
      if (lag == 0) {
	sprintf(maskFileName, "%s/transfer/mask/mask%s.%s",
		resDir, dstPostfix, fileFormat);
	std::cout << "Scanning mask from "
		  << maskFileName << std::endl;
	transferOp->scanMask(maskFileName, fileFormat);
	// Save mask
	mask = gsl_vector_uint_alloc(transferOp->getN());
	gsl_vector_uint_memcpy(mask, transferOp->mask);
      
	// Allocate memory to distributions
	transferOp->allocateDist();

	// Scan initial distribution for the first lag
	sprintf(initDistFileName, "%s/transfer/initDist/initDist%s.%s",
		resDir, dstPostfix, fileFormat);
	std::cout << "Scanning initial distribution from "
		  << initDistFileName << std::endl;
	transferOp->scanInitDist(initDistFileName, fileFormat);
	// Save initial distribution
	initDist = gsl_vector_alloc(transferOp->getNFilled());
	gsl_vector_memcpy(initDist, transferOp->initDist);
      }
    }
    catch (std::exception &ex) {
      std::cerr << "Error reading transfer operator: " << ex.what()
		<< std::endl;
      return EXIT_FAILURE;
    }
    // Get spectrum
    try {
      // Solve eigen value problem with default configuration
      transferSpec = new transferSpectrum(nev, transferOp, config);

      if (getForwardEigenvectors) {
	std::cout
	  << "Solving left eigenproblem for transition matrix..."
	  << std::endl;
	transferSpec->getSpectrumForward();
	std::cout << "Found " << transferSpec->getNev()
		  << "/" << nev << " eigenvalues." << std::endl;	      
      }
      if (getBackwardEigenvectors) {
	std::cout
	  << "Solving right eigenproblem for transition matrix..."
	  << std::endl;
	transferSpec->getSpectrumBackward();
	std::cout << "Found " << transferSpec->getNev()
		  << "/" << nev << " eigenvalues." << std::endl;
      }
      if (getForwardEigenvectors && getBackwardEigenvectors
	  && makeBiorthonormal) {
	std::cout << "Making set of forward and backward eigenvectors \
biorthonormal..." << std::endl;
	transferSpec->makeBiorthonormal();
      }
    }
    catch (std::exception &ex) {
      std::cerr << "Error calculating spectrum: " << ex.what() << std::endl;
      return EXIT_FAILURE;
    }
  
    // Write spectrum 
    try {
      if (getForwardEigenvectors) {
	std::cout << "Writing forward eigenvalues and eigenvectors..."
		  << std::endl;
	transferSpec->writeSpectrumForward(EigValForwardFileName,
					   EigVecForwardFileName,
					   fileFormat);
      }
      if (getBackwardEigenvectors) {
	std::cout << "Writing backward eigenvalues and eigenvectors..."
		  << std::endl;
	transferSpec->writeSpectrumBackward(EigValBackwardFileName,
					    EigVecBackwardFileName,
					    fileFormat);
      }
    }
    catch (std::exception &ex) {
      std::cerr << "Error writing spectrum: " << ex.what() << std::endl;
      return EXIT_FAILURE;
    }

    // Free
    delete transferOp;
    delete transferSpec;
  }
  if (initDist)
    gsl_vector_free(initDist);
  if (mask)
    gsl_vector_uint_free(mask);

  // Free configuration
  freeConfig();
  
  return 0;
}
