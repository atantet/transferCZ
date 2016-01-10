#define DIM 2
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_uint.h>
#include <gsl/gsl_matrix.h>
#include <Eigen/Sparse>
#include <libconfig.h++>
#include <arpack++/arlnsmat.h>
#include <arpack++/arlsnsym.h>
#include <ATSuite/atio.hpp>
#include <ATSuite/transferOperator.hpp>
#include <ATSuite/atspectrum.hpp>


/** \file spectrumCZ.cpp
 *  \brief Get spectrum of transfer operators for the Cane-Zebiak model.
 *   
 * Get spectrum of transfer operators for the Cane-Zebiak model.
 * For configuration, the file cfg/transferCZ.cfg is parsed using libconfig C++ library.
 * The transition matrices are then read from matrix files in compressed format.
 * The Eigen problem is then defined and solved using ARPACK++.
 * Finally, the results are written to file.
 */

/** \brief Eigen sparse CSR matrix of double type. */
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMatCSR;
/** \brief Eigen sparse CSC matrix of double type. */
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMatCSC;

using namespace libconfig;


// Declarations
/**
 * User defined function to get parameters from config file cfg/transferCZ.cfg using libconfig.
 */
int readConfig(char *);

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
const char indexDir[] = "../observables/";

// Simulation parametrs
const size_t nYears = 500;
const double sampFreq = 5.833333333333333;
const size_t nSeeds = 9;

// Arpack configuration class
configAR cfgAR("LM");

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
int minNumberStates;


// Main program
int main(int argc, char * argv[])
{
  // Read configuration file
  if (readConfig("transferCZ")) {
    std::cerr << "Error reading config file " << "transferCZ" << "."
	      << std::endl;
    return(EXIT_FAILURE);
  }

  
  // Declarations
  // Grid
  size_t N = 1;

  // Lags
  double tauDim;

  // Filtering
  double alpha;

  // Eigen problem
  int nconv;

  // Names and Files
  char obsName[256], srcPostfix[256], gridPostfix[256], cpyBuffer[256],
    postfix[256];
  char forwardTransitionFileName[256], backwardTransitionFileName[256],
    initDistFileName[256], finalDistFileName[256];

  char EigValFileName[256];
  char EigVecFileName[256];
  FILE *forwardTransitionFile, *backwardTransitionFile,
    *initDistFile, *finalDistFile, *EigValFile, *EigVecFile;

  // Transition Matrix
  ARluNonSymMatrix<double, double> *PT, *QT;
  SpMatCSR PCSR, QCSR;
  SpMatCSC PTCSC, QTCSC;
  gsl_vector *initDist, *finalDist;
  
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
  
  // Eigen problem
  double *EigValReal = new double [nev+1];
  double *EigValImag = new double [nev+1];
  double *EigVec = new double [(nev+2)*N];
  
  // Scan transition matrices and distributions for different lags
  initDist = gsl_vector_alloc(N);
  finalDist = gsl_vector_alloc(N);
  for (size_t lag = 0; lag < nLags; lag++){
    tauDim = gsl_vector_get(tauDimRng, lag);
    alpha = minNumberStates * 1. / ((nYears*sampFreq*12 - tauDim)*nSeeds);
    std::cout << "alpha = " << alpha << std::endl;

    // Open source files from which to read the transition matrices and distributions
    sprintf(postfix, "%s_tau%03d", gridPostfix, (int) (tauDim * 1000));
    // Forward transition matrix
    sprintf(forwardTransitionFileName, \
	    "../results/transitionMatrix/forwardTransition%s.csr", postfix);
    if ((forwardTransitionFile = fopen(forwardTransitionFileName, "r")) \
	== NULL){
      fprintf(stderr, "Can't open %s for reading!\n",
	      forwardTransitionFileName);
      return -1;
    }
    // Backward transition matrix
    sprintf(backwardTransitionFileName,
	    "../results/transitionMatrix/backwardTransition%s.csr", postfix);
    if ((backwardTransitionFile = fopen(backwardTransitionFileName, "r"))
	== NULL){
      fprintf(stderr, "Can't open %s for reading!\n",
	      backwardTransitionFileName);
      return -1;
    }
    // Initial distribution
    sprintf(initDistFileName, "../results/transitionMatrix/initDist%s.txt", postfix);
    if ((initDistFile = fopen(initDistFileName, "r")) \
	== NULL){
      fprintf(stderr, "Can't open %s for reading!\n",
	      initDistFileName);
      return -1;
    }
    // Final distribution
    sprintf(finalDistFileName, "../results/transitionMatrix/finalDist%s.txt", postfix);
    if ((finalDistFile = fopen(finalDistFileName, "r")) \
	== NULL){
      fprintf(stderr, "Can't open %s for reading!\n",
	      finalDistFileName);
      return -1;
    }

    // Read transition matrices written in CSR as CSC to take the transpose
    std::cout << std::endl << "Reading transpose forward transition matrix for lag "
	      << tauDim << " in " << forwardTransitionFileName << std::endl;
    // Read initial and final distributions
    gsl_vector_fscanf(initDistFile, initDist);
    fclose(initDistFile);
    gsl_vector_fscanf(finalDistFile, finalDist);
    fclose(finalDistFile);
    // Read forward transition matrix
    PCSR = *Compressed2Eigen(forwardTransitionFile);
    fclose(forwardTransitionFile);
    // Filter and get left stochastic
    filterTransitionMatrix(&PCSR, initDist, finalDist, alpha, 2);
    // Get transpose and convert to arpack
    PTCSC = SpMatCSC(PCSR.transpose());
    PT = Eigen2AR(&PTCSC);

    // Declare eigen problem: real non-symetric (see examples/areig.h)
    std::cout << "Solving eigen problem for the first " << nev
	      << " eigenvalues on a square matrix of size " << PT->nrows()
	      << " with non-zeros " << PT->nzeros()
	      << std::endl;
    // Call
    nconv = getSpectrum(PT, nev, cfgAR, EigValReal, EigValImag, EigVec);
    std::cout << "Found " << nconv << "/" << nev << " eigenvalues."
	      << std::endl;

    // Open destination files and write spectrum
    sprintf(EigValFileName, "../results/spectrum/eigval/eigval_nev%d%s.txt", nev, postfix);
    sprintf(EigVecFileName, "../results/spectrum/eigvec/eigvec_nev%d%s.txt", nev, postfix);
    if ((EigValFile = fopen(EigValFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", EigValFileName);
      return -1;
    }
    if ((EigVecFile = fopen(EigVecFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", EigVecFileName);
      return -1;
    }
    std::cout << "Write eigenvalues to " << EigValFileName << std::endl;
    std::cout << "and eigenvectors to " << EigVecFileName << std::endl;
    writeSpectrum(EigValFile, EigVecFile, EigValReal, EigValImag, EigVec,
		  nev, N);
    fclose(EigValFile);
    fclose(EigVecFile);
    delete PT;


    // Solve the adjoint problem on the transpose (of the transpose)
    // Read backward transition matrix
    std::cout << std::endl << "Reading backward transition matrix for lag "
	      << tauDim << " in " << backwardTransitionFileName << std::endl;
    QCSR = *Compressed2Eigen(backwardTransitionFile);
    fclose(backwardTransitionFile);
    // Filter and get left stochastic
    filterTransitionMatrix(&QCSR, finalDist, initDist, alpha, 2);
    // Get transpose and convert to arpack
    QTCSC = SpMatCSC(QCSR.transpose());
    QT = Eigen2AR(&QTCSC);

    // Declare eigen problem: real non-symetric (see examples/areig.h)
    std::cout << "Solving eigen problem for the first " << nev
	      << " eigenvalues" << std::endl;
    // Call
    nconv = getSpectrum(QT, nev, cfgAR, EigValReal, EigValImag, EigVec);
    std::cout << "Found " << nconv << "/" << nev << " eigenvalues."
	      << std::endl;

    // Open destination files and write spectrum
    sprintf(EigValFileName, "../results/spectrum/eigval/eigvalAdjoint_nev%d%s.txt",
	    nev, postfix);
    sprintf(EigVecFileName, "../results/spectrum/eigvec/eigvecAdjoint_nev%d%s.txt",
	    nev, postfix);
    if ((EigValFile = fopen(EigValFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", EigValFileName);
      return -1;
    }
    if ((EigVecFile = fopen(EigVecFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", EigVecFileName);
      return -1;
    }
    std::cout << "Write adjoint eigenvalues to " << EigValFileName
	      << std::endl;
    std::cout << "and adjoint eigenvectors to " << EigVecFileName
	      << std::endl;
    writeSpectrum(EigValFile, EigVecFile, EigValReal, EigValImag, EigVec,
		  nev, N);
    fclose(EigValFile);
    fclose(EigVecFile);
    delete QT;
  }
  
  // Clean-up
  delete[] EigValReal;
  delete[] EigValImag;
  delete[] EigVec;
  gsl_vector_free(initDist);
  gsl_vector_free(finalDist);
  
  return 0;
}


// Definitions
int readConfig(char *cfgFileNamePrefix)
{
  char cfgFileName[256];
  sprintf(cfgFileName, "cfg/%s.cfg", cfgFileNamePrefix);

  // Read the file. If there is an error, report it and exit.
  try {
    std::cout << "Reading config file " << cfgFileName << std::endl;
    cfg.readFile(cfgFileName);
  }
  catch(const FileIOException &fioex) {
    std::cerr << "I/O error while reading configuration file." << std::endl;
    return(EXIT_FAILURE);
  }
  catch(const ParseException &pex) {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    return(EXIT_FAILURE);
  }

  try {
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
	return(EXIT_FAILURE);
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
    const Setting &tauDimRngSetting = cfg.lookup("transition.tauDimRng");
    nLags = tauDimRngSetting.getLength();
    tauDimRng = gsl_vector_alloc(nLags);

    std::cout << std::endl << "---transition---" << std::endl;
    std::cout << "tauDimRng = {";
    for (size_t lag = 0; lag < nLags; lag++) {
      gsl_vector_set(tauDimRng, lag, tauDimRngSetting[lag]);
      std::cout << gsl_vector_get(tauDimRng, lag) << ", ";
    }
    std::cout << "}" << std::endl;
  }
  catch(const SettingNotFoundException &nfex) {
    std::cerr << "Setting " << nfex.getPath() << " not found." << std::endl;
    return(EXIT_FAILURE);
  }

  // Get spectrum setting 
  nev = cfg.lookup("spectrum.nev");
  minNumberStates = cfg.lookup("spectrum.minNumberStates");
  std::cout << std::endl << "---spectrum---" << std::endl
	    << "nev: " << nev << std::endl
	    << "minNumberStates: " << minNumberStates << std::endl;

  std::cout << std::endl;

  return 0;
}
