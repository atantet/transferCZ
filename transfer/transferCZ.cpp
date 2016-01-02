#define DIM 2
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_uint.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_uint.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <libconfig.h++>
#include "transferOperator.hpp"
#include "atio.hpp"

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMatCSC;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMatCSR;

using namespace libconfig;


// Declarations
int readConfig(char *);
struct fieldDef {
  int tableIdx;
  char shortName[256];
  double toDimFactor;
};
const char prefix[] = "zc_1eof_";
const char indexDir[] = "../observables/";
const size_t dimSeries = 5;
const double L = 1.5e7;
const double c0 = 2.;
const double timeDim = L / c0 / (60. * 60. * 24);
const double H = 200.;
const double tau_0 = 1.0922666667e-2;
const double delta_T = 1.;
const struct fieldDef field_h = {1, "h", H};
const struct fieldDef field_T = {2, "T", delta_T};

// Configuration 
char cfgFileName[256];
Config cfg;
std::string indicesName[DIM];
struct fieldDef fieldsDef[DIM];
double mu, eps;
gsl_vector_uint *seedRng;
size_t nSeeds;
bool readGridMem;
double dt, spinupMonth;
gsl_vector_uint *nx;
gsl_vector *nSTDLow, *nSTDHigh;
size_t nLags;
gsl_vector *tauDimRng;


int main(int argc, char * argv[])
{
  // Read configuration file
  if (readConfig(argv[0])) {
    std::cerr << "Error reading config file " << argv[0] << ".cfg"
	      << std::endl;
    return(EXIT_FAILURE);
  }

  const double sampFreq = 0.35 / dt;
  size_t spinup = (size_t) (spinupMonth * sampFreq);

  // Names and Files
  char obsName[256], srcPostfix[256], srcPostfixSeed[256], postfix[256],
    gridPostfix[256], gridCFG[256], gridPostfixSeed[256], cpyBuffer[256];
  char gridMemFileName[256],
    forwardTransitionFileName[256], backwardTransitionFileName[256],
    initDistFileName[256], finalDistFileName[256];
  char gridFileName[256];
  FILE *gridMemFile, *gridFile,
    *forwardTransitionFile, *backwardTransitionFile,
    *initDistFile, *finalDistFile;
  FILE *indexFile[DIM];
  char indexPath[DIM][256];
  size_t count;
  
  // Vectors and matrices
  size_t tauStep;
  double tauDim;
  vector<gsl_vector_uint *> gridMemSeeds(nSeeds);
  gsl_matrix_uint *gridMemMatrix, *gridMemMatrixSeed;
  SpMatCSR *P, *Q;
  gsl_vector *initDist, *finalDist;
  gsl_matrix *data;
  gsl_vector *xmin;
  gsl_vector *xmax;
  gsl_vector_uint *ntIndex = gsl_vector_uint_alloc(DIM);
  gsl_vector *statesMean = gsl_vector_calloc(DIM);
  gsl_vector *statesSTD = gsl_vector_calloc(DIM);
  std::vector<gsl_matrix *> statesSeeds(nSeeds);
  gsl_vector_uint *ntSeeds = gsl_vector_uint_alloc(nSeeds);
  size_t ntTot = 0;

  // Grid related
  size_t N = 1;
  std::vector<gsl_vector *> *gridBounds;

  // Define grid name
  sprintf(srcPostfix, "%smu%04d_eps%04d", prefix, (int) (mu*1000),
	  (int) (eps*1000));
  sprintf(gridCFG, "");
  for (size_t d = 0; d < DIM; d++) {
    N *= gsl_vector_uint_get(nx, d);
    strcpy(cpyBuffer, gridCFG);
    sprintf(gridCFG, "%s_n%dl%dh%d", cpyBuffer,
	    gsl_vector_uint_get(nx, d),
	    (int) gsl_vector_get(nSTDLow, d),
	    (int) gsl_vector_get(nSTDHigh, d));
  }

  // Get membership vector
  if (! readGridMem) {
    for (size_t seed = 0; seed < nSeeds; seed++) {
      // Read observable
      sprintf(srcPostfixSeed, "%s_seed%d", srcPostfix, (int) seed);
      strcpy(obsName, "");
      for (size_t d = 0; d < DIM; d++){
	strcpy(cpyBuffer, obsName);
	sprintf(obsName, "%s_%s_%s", cpyBuffer, fieldsDef[d].shortName,
		indicesName[d].c_str());
	sprintf(indexPath[d], "%s/%s/%s.txt", indexDir, srcPostfixSeed,
		indicesName[d].c_str());
	// Open observable file d
	if ((indexFile[d] = fopen(indexPath[d], "r")) == NULL){
	  fprintf(stderr, "Can't open %s for reading!\n", indexPath[d]);
	  return(EXIT_FAILURE);
	}
	// Get number of lines
	count = lineCount(indexFile[d]);
	gsl_vector_uint_set(ntIndex, d, count);
	fseek(indexFile[d], 0, SEEK_SET);
      }
      // Get index length and allocate
      gsl_vector_uint_set(ntSeeds, seed, gsl_vector_uint_min(ntIndex) - spinup);
      ntTot += gsl_vector_uint_get(ntSeeds, seed);
      statesSeeds[seed] = gsl_matrix_alloc(gsl_vector_uint_get(ntSeeds, seed),
					   DIM);

      // Read trajectory and get the total mean and standard defination
      for (size_t d = 0; d < DIM; d++){
	data = gsl_matrix_alloc(gsl_vector_uint_get(ntIndex, d), dimSeries);
	// Read trajectory
	std::cout << "Reading trajectory in " << indexPath[d] << std::endl;
	gsl_matrix_fscanf(indexFile[d], data);
	fclose(indexFile[d]);
	for (size_t k = 0; k < gsl_vector_uint_get(ntSeeds, seed); k++) {
	  gsl_matrix_set(statesSeeds[seed], k, d,
			 gsl_matrix_get(data, spinup+k, fieldsDef[d].tableIdx)
			 * fieldsDef[d].toDimFactor);
	  gsl_vector_set(statesMean, d, gsl_vector_get(statesMean, d) \
			 + gsl_matrix_get(statesSeeds[seed], k, d));
	  gsl_vector_set(statesSTD, d, gsl_vector_get(statesSTD, d)
			 + pow(gsl_matrix_get(statesSeeds[seed], k, d), 2));
	}
	gsl_matrix_free(data);
      }
    }
    
    // Define grid bounds
    xmin = gsl_vector_alloc(DIM);
    xmax = gsl_vector_alloc(DIM);
    for (size_t d = 0; d < DIM; d++) {
      gsl_vector_set(statesMean, d, gsl_vector_get(statesMean, d) / ntTot);
      gsl_vector_set(statesSTD, d,
		     sqrt(gsl_vector_get(statesSTD, d) / ntTot
			  - pow(gsl_vector_get(statesMean, d), 2)));
      gsl_vector_set(xmin, d, gsl_vector_get(statesMean, d)
		     - gsl_vector_get(nSTDLow, d)
		     * gsl_vector_get(statesSTD, d));
      gsl_vector_set(xmax, d, gsl_vector_get(statesMean, d)
		     + gsl_vector_get(nSTDHigh, d)
		     * gsl_vector_get(statesSTD, d));
    }

    // Open grid
    sprintf(gridPostfix, "_%s%s%s", srcPostfix, obsName, gridCFG);
    sprintf(gridFileName, "grid/grid%s.txt", gridPostfix);
    if ((gridFile = fopen(gridFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", gridFileName);
      return(EXIT_FAILURE);
    }
    // Define grid
    gridBounds = getGridRect(nx, xmin, xmax);
    // Write Grid
    writeGridRect(gridFile, gridBounds, true);
    // Close
    fclose(gridFile);

    // Get grid membership for each seed
    for (size_t seed = 0; seed < nSeeds; seed++){
      // Open grid membership file
      sprintf(srcPostfixSeed, "%s_seed%d", srcPostfix, (int) seed);
      sprintf(gridPostfixSeed, "_%s%s%s", srcPostfixSeed, obsName, gridCFG);
      sprintf(gridMemFileName, "transitionMatrix/gridMem%s.txt",
	      gridPostfixSeed);
      if ((gridMemFile = fopen(gridMemFileName, "w")) == NULL){
	fprintf(stderr, "Can't open %s for writing!\n", gridMemFileName);
	return(EXIT_FAILURE);
      }

      // Get grid membership
      std::cout << "Getting grid membership vector for seed" << seed
		<< std::endl;
      gridMemSeeds[seed] = getGridMembership(statesSeeds[seed], gridBounds);

      // Write grid membership
      gsl_vector_uint_fprintf(gridMemFile, gridMemSeeds[seed], "%d");
      fclose(gridMemFile);
    }

    // Free
    for (size_t d = 0; d < DIM; d++)
      gsl_vector_free((*gridBounds)[d]);
    delete gridBounds;
    gsl_vector_free(xmin);
    gsl_vector_free(xmax);
    for (size_t seed = 0; seed < nSeeds; seed++){
      gsl_matrix_free(statesSeeds[seed]);
    }
  }
  else {
    for (size_t seed = 0; seed < nSeeds; seed++) {
      gridMemSeeds[seed] = gsl_vector_uint_alloc(gsl_vector_uint_get(ntSeeds,
								     seed));
      // Open grid membership file to read
      std::cout << "Reading grid membership vector..." << std::endl;
      if ((gridMemFile = fopen(gridMemFileName, "r")) == NULL){
	fprintf(stderr, "Can't open %s for writing!\n", gridMemFileName);
	return(EXIT_FAILURE);
      }
      gsl_vector_uint_fscanf(gridMemFile, gridMemSeeds[seed]);
      fclose(gridMemFile);
    }
  }

  // Get transition matrices for ifferent lags
  for (size_t lag = 0; lag < nLags; lag++){
    tauDim = gsl_vector_get(tauDimRng, lag);
    tauStep = (size_t) (tauDim * sampFreq);

    // Open forward and backward files
    sprintf(postfix, "%s_tau%03d", gridPostfix, (int) (tauDim * 1000));
    // Forward transition matrix
    sprintf(forwardTransitionFileName,
	    "transitionMatrix/forwardTransition%s.csr", postfix);
    if ((forwardTransitionFile = fopen(forwardTransitionFileName, "w"))
	== NULL){
      fprintf(stderr, "Can't open %s for writing!\n",
	      forwardTransitionFileName);
      return(EXIT_FAILURE);
    }
    // Backward transition matrix
    sprintf(backwardTransitionFileName,
	    "transitionMatrix/backwardTransition%s.csr", postfix);
    if ((backwardTransitionFile = fopen(backwardTransitionFileName, "w"))
	== NULL){
      fprintf(stderr, "Can't open %s for writing!\n",
	      backwardTransitionFileName); 
      return(EXIT_FAILURE);
    }
    // Initial distribution
    sprintf(initDistFileName, "transitionMatrix/initDist%s.txt", postfix);
    if ((initDistFile = fopen(initDistFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", initDistFileName);
      return(EXIT_FAILURE);
    }
    // Final distribution
    sprintf(finalDistFileName, "transitionMatrix/finalDist%s.txt", postfix);
    if ((finalDistFile = fopen(finalDistFileName, "w")) == NULL){
      fprintf(stderr, "Can't open %s for writing!\n", finalDistFileName);
      return(EXIT_FAILURE);
    }

    // Get full membership matrix
    std::cout << "Getting full membership matrix" << std::endl;
    count = 0;
    gridMemMatrix = gsl_matrix_uint_alloc(ntTot - tauStep * nSeeds, 2);
    for (size_t seed = 0; seed < nSeeds; seed++) {
      gridMemMatrixSeed = getGridMembership(gridMemSeeds[seed], tauStep);
      for (size_t t = 0; t < gsl_vector_uint_get(ntSeeds, seed)-tauStep; t++) {
	gsl_matrix_uint_set(gridMemMatrix, count, 0,
			    gsl_matrix_uint_get(gridMemMatrixSeed, t, 0));
	gsl_matrix_uint_set(gridMemMatrix, count, 1,
			    gsl_matrix_uint_get(gridMemMatrixSeed, t, 1));
	count++;
      }
      gsl_matrix_uint_free(gridMemMatrixSeed);
    }

    // Get transition matrices as CSR 
    P = new SpMatCSR(N, N);
    Q = new SpMatCSR(N, N);
    initDist = gsl_vector_alloc(N);
    finalDist = gsl_vector_alloc(N);
    std::cout << "Getting transition matrix" << std::endl;
    getTransitionMatrix(gridMemMatrix, N, P, Q, initDist, finalDist);
    
    // Write transition matrix as CSR
    std::cout << "Writing forward transition matrix to "
	      << forwardTransitionFileName << std::endl;
    Eigen2Compressed(forwardTransitionFile, P);
    fclose(forwardTransitionFile);
    std::cout << "Writing backward transition matrix to "
	      << backwardTransitionFileName << std::endl;
    Eigen2Compressed(backwardTransitionFile, Q);
    fclose(backwardTransitionFile);
    gsl_vector_fprintf(initDistFile, initDist, "%f");
    fclose(initDistFile);
    gsl_vector_fprintf(finalDistFile, finalDist, "%f");
    fclose(finalDistFile);

    // Free
    gsl_vector_free(initDist);
    gsl_vector_free(finalDist);
    gsl_matrix_uint_free(gridMemMatrix);
    delete P;
    delete Q;
  }

  // Free
  gsl_vector_uint_free(nx);
  gsl_vector_free(nSTDLow);
  gsl_vector_free(nSTDHigh);
  gsl_vector_free(tauDimRng);
  gsl_vector_uint_free(seedRng);
  gsl_vector_uint_free(ntSeeds);
		
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
    const Setting &seedRngSetting = cfg.lookup("caseDef.seedRng");
    nSeeds = seedRngSetting.getLength();
    seedRng = gsl_vector_uint_alloc(nSeeds);
    std::cout << "seedRng = {";
    for (size_t seed = 0; seed < nSeeds; seed++) {
      gsl_vector_uint_set(seedRng, seed, seedRngSetting[seed]);
      std::cout << gsl_vector_uint_get(seedRng, seed) << ", ";
    }
    std::cout << "}" << std::endl;

      // Get simulation settings
    dt = cfg.lookup("simulation.dt");
    spinupMonth = cfg.lookup("simulation.spinupMonth");
    std::cout << endl << "---simulation---" << std::endl
  	      << "dt: " << dt << std::endl
  	      << "spinupMonth: " << spinupMonth << std::endl;

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
    readGridMem = cfg.lookup("grid.readGridMem");
    std::cout << "readGridMem: " << readGridMem << std::endl;


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

  std::cout << std::endl;

  return 0;
}
