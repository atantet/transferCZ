#define DIM 2
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_int.h>
#include <gsl/gsl_matrix.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <libconfig.h++>
#include "arlnsmat.h"
#include "arlsnsym.h"
#include "atio.hpp"
#include "transferOperator.hpp"


typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMatCSC;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMatCSR;

using namespace libconfig;


// Declarations
void writeSpectrum(FILE *, FILE *, double *, double *, double *, int, size_t);
int readConfig(char *);
struct fieldDef {
  int tableIdx;
  char shortName[256];
  double toDimFactor;
};
const char prefix[] = "zc_1eof_";
const char indexDir[] = "../observables/";
const double L = 1.5e7;
const double c0 = 2.;
const double timeDim = L / c0 / (60. * 60. * 24);
const double H = 200.;
const double tau_0 = 1.0922666667e-2;
const double delta_T = 1.;
const struct fieldDef field_h = {1, "h"};
const struct fieldDef field_T = {2, "T"};

const size_t nYears = 500;
const double sampFreq = 5.833333333333333;
const size_t nSeeds = 9;

const std::string& which = "LM"; // Look for eigenvalues of Largest Magnitude
int ncv = 0;
double tol = 0.;
int maxit = 0;
double *resid = NULL;
bool AutoShift = true;

// Configuration 
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


int main(int argc, char * argv[])
{
  // Read configuration file
  if (readConfig("transferCZ")) {
    std::cerr << "Error reading config file " << "transferCZ" << "."
	      << std::endl;
    return(EXIT_FAILURE);
  }

  // Definitions
  // Lags
  double tauDim;

  // Filtering
  double alpha;

  // Names and Files
  size_t N = 1;
  char obsName[256], srcPostfix[256], gridPostfix[256], cpyBuffer[256],
    postfix[256];
  char forwardTransitionFileName[256], backwardTransitionFileName[256],
    initDistFileName[256], finalDistFileName[256];

  char EigValFileName[256];
  char EigVecFileName[256];
  FILE *forwardTransitionFile, *backwardTransitionFile,
    *initDistFile, *finalDistFile, *EigValFile, *EigVecFile;
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
  
  // Transition Matrix
  ARluNonSymMatrix<double, double> *PT, *QT;
  SpMatCSR PCSR, QCSR;
  SpMatCSC PTCSC, QTCSC;
  gsl_vector *initDist, *finalDist;
  
  // Eigen problem
  double *EigValReal = new double [nev+1];
  double *EigValImag = new double [nev+1];
  double *EigVec = new double [(nev+2)*N];
  ARluNonSymStdEig<double> EigProb;
  int nconv;
  
  // Get transition matrices for different lags
  initDist = gsl_vector_alloc(N);
  finalDist = gsl_vector_alloc(N);
  for (size_t lag = 0; lag < nLags; lag++){
    tauDim = gsl_vector_get(tauDimRng, lag);
    alpha = minNumberStates * 1. / ((nYears*sampFreq*12 - tauDim)*nSeeds);
    std::cout << "alpha = " << alpha << std::endl;

    // Open source and destination files
    sprintf(postfix, "%s_tau%03d", gridPostfix, (int) (tauDim * 1000));
    sprintf(forwardTransitionFileName, \
	    "transitionMatrix/forwardTransition%s.csr", postfix);
    if ((forwardTransitionFile = fopen(forwardTransitionFileName, "r")) \
	== NULL){
      fprintf(stderr, "Can't open %s for reading!\n",
	      forwardTransitionFileName);
      return -1;
    }
    sprintf(initDistFileName, "transitionMatrix/initDist%s.txt", postfix);
    if ((initDistFile = fopen(initDistFileName, "r")) \
	== NULL){
      fprintf(stderr, "Can't open %s for reading!\n",
	      initDistFileName);
      return -1;
    }
    sprintf(finalDistFileName, "transitionMatrix/finalDist%s.txt", postfix);
    if ((finalDistFile = fopen(finalDistFileName, "r")) \
	== NULL){
      fprintf(stderr, "Can't open %s for reading!\n",
	      finalDistFileName);
      return -1;
    }

    // Read transition matrix written in CSR as CSC to take the transpose
    std::cout << endl << "Reading transpose forward transition matrix for lag "
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
    EigProb = ARluNonSymStdEig<double>(nev, *PT, which, ncv, tol, maxit, resid,
				       AutoShift);

    // Open destination files and write spectrum
    sprintf(EigValFileName, "spectrum/eigval/eigval_nev%d%s.txt", nev, postfix);
    sprintf(EigVecFileName, "spectrum/eigvec/eigvec_nev%d%s.txt", nev, postfix);
    if ((EigValFile = fopen(EigValFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", EigValFileName);
      return -1;
    }
    if ((EigVecFile = fopen(EigVecFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", EigVecFileName);
      return -1;
    }

    // Find eigenvalues and left eigenvectors
    EigProb.EigenValVectors(EigVec, EigValReal, EigValImag);
    nconv = EigProb.ConvergedEigenvalues();
    std::cout << "Found " << nconv << "/" << nev << " eigenvalues."
	      << std::endl;

    // Write results
    std::cout << "Write eigenvalues to " << EigValFileName << std::endl;
    std::cout << "and eigenvectors to " << EigVecFileName << std::endl;
    writeSpectrum(EigValFile, EigVecFile, EigValReal, EigValImag, EigVec,
		  nev, N);
    fclose(EigValFile);
    fclose(EigVecFile);
    delete PT;


    // Sove the adjoint problem on the transpose (of the transpose)
    // Read transition matrix written in CSR and convert to CSC
    std::cout << std::endl << "Reading backward transition matrix for lag "
	      << tauDim << " in " << backwardTransitionFileName << std::endl;
    sprintf(backwardTransitionFileName,
	    "transitionMatrix/backwardTransition%s.csr", postfix);
    if ((backwardTransitionFile = fopen(backwardTransitionFileName, "r"))
	== NULL){
      fprintf(stderr, "Can't open %s for reading!\n",
	      backwardTransitionFileName);
      return -1;
    }
    // Read backward transition matrix
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
    EigProb = ARluNonSymStdEig<double>(nev, *QT, which, ncv, tol, maxit,
				       resid, AutoShift);

    // Open destination files and write spectrum
    sprintf(EigValFileName, "spectrum/eigval/eigvalAdjoint_nev%d%s.txt",
	    nev, postfix);
    sprintf(EigVecFileName, "spectrum/eigvec/eigvecAdjoint_nev%d%s.txt",
	    nev, postfix);
    if ((EigValFile = fopen(EigValFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", EigValFileName);
      return -1;
    }
    if ((EigVecFile = fopen(EigVecFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", EigVecFileName);
      return -1;
    }
  
    // Find eigenvalues and left eigenvectors
    EigProb.EigenValVectors(EigVec, EigValReal, EigValImag);
    nconv = EigProb.ConvergedEigenvalues();
    std::cout << "Found " << nconv << "/" << nev << " adjoint eigenvalues."
	      << std::endl;

    // Write results
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

// Write complex eigenvalues and eigenvectors
void writeSpectrum(FILE *fEigVal, FILE *fEigVec,
		   double *EigValReal, double *EigValImag,
		   double *EigVec, int nev, size_t N)
{
  size_t vecCount = 0;
  int ev =0;
  // Write real and imaginary parts of each eigenvalue on each line
  // Write on each pair of line the real part of an eigenvector then its imaginary part
  while (ev < nev) {
    // Always write the eigenvalue
    fprintf(fEigVal, "%lf %lf\n", EigValReal[ev], EigValImag[ev]);
    // Always write the real part of the eigenvector ev
    for (size_t i = 0; i < N; i++){
      fprintf(fEigVec, "%lf ", EigVec[vecCount*N+i]);
    }
    fprintf(fEigVec, "\n");
    vecCount++;
    
    // Write its imaginary part or the zero vector
    if (EigValImag[ev] != 0.){
      for (size_t i = 0; i < N; i++)
	fprintf(fEigVec, "%lf ", EigVec[vecCount*N+i]);
      vecCount++;
      // Skip the conjugate
      ev += 2;
    }
    else{
      for (size_t i = 0; i < N; i++)
	fprintf(fEigVec, "%lf ", 0.);
      ev += 1;
    }
    fprintf(fEigVec, "\n");
  }

  return;
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
    std::cout << endl << "---caseDef---" << std::endl;
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
    std::cout << endl << "---simulation---" << std::endl
  	      << "dt: " << dt << std::endl;

    // Get grid settings
    std::cout << endl << "---grid---" << std::endl;
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
  std::cout << endl << "---spectrum---" << std::endl
	    << "nev: " << nev << std::endl
	    << "minNumberStates: " << minNumberStates << std::endl;

  std::cout << std::endl;

  return 0;
}
