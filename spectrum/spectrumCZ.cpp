#define DIM 2
#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <libconfig.h++>
#include <ergoPack/transferOperator.hpp>
#include <ergoPack/transferSpectrum.hpp>


/** \file spectrumCZ.cpp
 *  \brief Get spectrum of transfer operators for the Cane-Zebiak model.
 *   
 * Get spectrum of transfer operators for the Cane-Zebiak model.
 * For configuration, the file cfg/transferCZ.cfg is parsed using libconfig C++ library.
 * The transition matrices are then read from matrix files in compressed format.
 * The Eigen problem is then defined and solved using ARPACK++.
 * Finally, the results are written to file.
 */

using namespace libconfig;


// Declarations
/** \brief User defined function to get parameters from a cfg file using libconfig. */
void readConfig(const char *cfgFileName);

/**
 * Structure defining a field.
 */
struct fieldDef {
  /**
   * Column index of the field in the data file.
   */
  int tableIdx;
  /**
   * Short field name used for file names.
   */
  char shortName[256];
  /**
   * Factor to convert the field from adimensional to dimensional units.
   */
  double toDimFactor;
};
const struct fieldDef field_h = {1, "h"};
const struct fieldDef field_T = {2, "T"};

// Paths
const char prefix[] = "zc_1eof_";
const char resDir[] = "../results/";

// Simulation parametrs
const size_t nYears = 500;
const double sampFreq = 5.833333333333333;
const size_t nSeeds = 9;

// Global variables used for configuration 
Config cfg;
std::string indicesName[DIM];
struct fieldDef fieldsDef[DIM];
double mu, eps;
double dt;
gsl_vector_uint *nx;
gsl_vector *nSTDLow, *nSTDHigh;
size_t nLags;
gsl_vector *tauDimRng;
int nev;
configAR config;                //!< Configuration data for the eigen problem
char configFileName[256];       //!< Name of the configuration file
bool stationary;                //!< Whether the problem is stationary or not
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
  // Grid
  size_t N = 1;

  // Lags
  double tauDim;

  // Names and Files
  char obsName[256], srcPostfix[256], gridPostfix[256], cpyBuffer[256],
    postfix[256];
  char forwardTransitionFileName[256], initDistFileName[256];

  char EigValForwardFileName[256], EigVecForwardFileName[256],
    EigValBackwardFileName[256], EigVecBackwardFileName[256];

  // Transition Matrix
  transferOperator *transferOp;
  transferSpectrum *transferSpec;
  
  // Define grid name
  sprintf(srcPostfix, "%smu%04d_eps%04d", prefix, (int) (mu*1000),
	  (int) (eps*1000));
  strcpy(obsName, "");
  strcpy(gridPostfix, "");
  for (size_t d = 0; d < DIM; d++){
    strcpy(cpyBuffer, obsName);
    sprintf(obsName, "%s_%s_%s", cpyBuffer, fieldsDef[d].shortName,
	    indicesName[d].c_str());
    N *= gsl_vector_uint_get(nx, d);
    strcpy(cpyBuffer, gridPostfix);
    sprintf(gridPostfix, "%s_n%dl%dh%d", cpyBuffer, gsl_vector_uint_get(nx, d),
	    (int) gsl_vector_get(nSTDLow, d),
	    (int) gsl_vector_get(nSTDHigh, d));
  }
  strcpy(cpyBuffer, gridPostfix);
  sprintf(gridPostfix, "_%s%s%s", srcPostfix, obsName, cpyBuffer);

  // Scan transition matrices and distributions for different lags
  for (size_t lag = 0; lag < nLags; lag++){
    tauDim = gsl_vector_get(tauDimRng, lag);

    // Get file names
    sprintf(postfix, "%s_tau%03d", gridPostfix, (int) (tauDim * 1000));
    sprintf(forwardTransitionFileName, \
	    "%s/transfer/forwardTransition/forwardTransition%s.coo", resDir, postfix);
    sprintf(initDistFileName, "%s/transfer/initDist/initDist%s.txt",
	    resDir, postfix);
    sprintf(EigValForwardFileName, "%s/spectrum/eigval/eigvalForward_nev%d%s.txt",
	    resDir, nev, postfix);
    sprintf(EigVecForwardFileName, "%s/spectrum/eigvec/eigvecForward_nev%d%s.txt",
	    resDir, nev, postfix);
    sprintf(EigValBackwardFileName, "%s/spectrum/eigval/eigvalBackward_nev%d%s.txt",
	    resDir, nev, postfix);
    sprintf(EigVecBackwardFileName, "%s/spectrum/eigvec/eigvecBackward_nev%d%s.txt",
	    resDir, nev, postfix);

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
  
  return 0;
}


// Definitions
void
readConfig(const char *cfgFileName)
{
  Config cfg;

  // Read the file. If there is an error, report it and exit.
  try {
    std::cout << "Reading config file " << cfgFileName << std::endl;
    cfg.readFile(cfgFileName);
    
    std::cout.precision(6);
    std::cout << "Settings:" << std::endl;
    
    // Get caseDef settings
    std::cout << std::endl << "---caseDef---" << std::endl;
    // Indices and Fields
    const Setting &indicesNameSetting = cfg.lookup("caseDef.indicesName");
    const Setting &fieldsNameSetting = cfg.lookup("caseDef.fieldsName");
    for (size_t d = 0; d < DIM; d++) {
      indicesName[d] = (const char *) indicesNameSetting[d];
      if (!strcmp(fieldsNameSetting[d], "h"))
	fieldsDef[d] = field_h;
      else if (!strcmp(fieldsNameSetting[d], "T"))
	fieldsDef[d] = field_T;
      else {
	std::cerr << "Wrong field name for dimension " << d << std::endl;
	throw std::exception();
      }	
      std::cout << "indicesName[" << d << "]: "
		<< indicesName[d] << std::endl;
      std::cout << "fieldsName[" << d << "]: "
		<< fieldsDef[d].shortName << std::endl;
    }
    mu = cfg.lookup("caseDef.mu");
    eps = cfg.lookup("caseDef.eps");
    std::cout << "mu: " << mu << std::endl
	      << "eps: " << eps << std::endl;

      // Get simulation settings
    dt = cfg.lookup("simulation.dt");
    std::cout << std::endl << "---simulation---" << std::endl
  	      << "dt: " << dt << std::endl;

    // Get grid settings
    std::cout << std::endl << "---grid---" << std::endl;
    const Setting &nxSetting = cfg.lookup("grid.nx");
    const Setting &nSTDLowSetting = cfg.lookup("grid.nSTDLow");
    const Setting &nSTDHighSetting = cfg.lookup("grid.nSTDHigh");
    nx = gsl_vector_uint_alloc(DIM);
    nSTDLow = gsl_vector_alloc(DIM);
    nSTDHigh = gsl_vector_alloc(DIM);
    for (size_t d = 0; d < DIM; d++) {
      gsl_vector_uint_set(nx, d, nxSetting[d]);
      gsl_vector_set(nSTDLow, d, nSTDLowSetting[d]);
      gsl_vector_set(nSTDHigh, d, nSTDHighSetting[d]);
      std::cout << "Grid definition (nSTDLow, nSTDHigh, n):" << std::endl;
      std::cout << "dim " << d+1 << ": ("
		<< gsl_vector_get(nSTDLow, d) << ", "
		<< gsl_vector_get(nSTDHigh, d) << ", "
		<< gsl_vector_uint_get(nx, d) << ")" << std::endl;
    }

    // Get transition settings
    const Setting &tauDimRngSetting = cfg.lookup("transfer.tauDimRng");
    nLags = tauDimRngSetting.getLength();
    tauDimRng = gsl_vector_alloc(nLags);

    std::cout << std::endl << "---transfer---" << std::endl;
    std::cout << "tauDimRng = {";
    for (size_t lag = 0; lag < nLags; lag++) {
      gsl_vector_set(tauDimRng, lag, tauDimRngSetting[lag]);
      std::cout << gsl_vector_get(tauDimRng, lag) << ", ";
    }
    std::cout << "}" << std::endl;

    stationary = cfg.lookup("transfer.stationary");
    std::cout << "Is stationary: " << stationary << std::endl;


    // Get spectrum setting 
    nev = cfg.lookup("spectrum.nev");
    std::cout << std::endl << "---spectrum---" << std::endl;
    // Get eigen problem configuration
    config = defaultCfgAR;
    if (cfg.exists("spectrum.which"))
      {
	strcpy(config.which, (const char *) cfg.lookup("spectrum.which"));
      }
    if (cfg.exists("spectrum.ncv"))
      {
	config.ncv = cfg.lookup("spectrum.ncv");
      }
    if (cfg.exists("spectrum.tol"))
      {
	config.tol = cfg.lookup("spectrum.tol");
      }
    if (cfg.exists("spectrum.maxit"))
	{
	  config.maxit = cfg.lookup("spectrum.maxit");
	}
    if (cfg.exists("spectrum.AutoShift"))
	{
	  config.AutoShift = (bool) cfg.lookup("spectrum.AutoShift");
	}
    std::cout << "nev: " << nev << std::endl;
    std::cout << "which: " << config.which << std::endl;
    std::cout << "ncv: " << config.ncv << std::endl;
    std::cout << "tol: " << config.tol << std::endl;
    std::cout << "maxit: " << config.maxit << std::endl;
    std::cout << "AutoShift: " << config.AutoShift << std::endl;
    std::cout << std::endl;

    if (cfg.exists("spectrum.getForwardEigenvectors"))
      {
	getForwardEigenvectors = cfg.lookup("spectrum.getForwardEigenvectors");
      }
    if (cfg.exists("spectrum.getBackwardEigenvectors"))
      {
	getBackwardEigenvectors = cfg.lookup("spectrum.getBackwardEigenvectors");
      }
    if (cfg.exists("spectrum.makeBiorthonormal"))
      {
	makeBiorthonormal = cfg.lookup("spectrum.makeBiorthonormal");
      }

    
    std::cout << std::endl;

  }
  catch(const SettingNotFoundException &nfex) {
    std::cerr << "Setting " << nfex.getPath() << " not found." << std::endl;
    throw nfex;
  }
  catch(const FileIOException &fioex) {
    std::cerr << "I/O error while reading configuration file." << std::endl;
    throw fioex;
  }
  catch(const ParseException &pex) {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    throw pex;
  }
  catch(const SettingTypeException &stex) {
    std::cerr << "Setting type exception." << std::endl;
    throw stex;
  }
  catch(const std::exception &ex) {
    std::cerr << "Standard exception." << std::endl;
    throw ex;
  }

  return ;
}
