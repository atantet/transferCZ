#define DIM 2
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <libconfig.h++>
#include <ATSuite/transferOperator.hpp>

/** \file transferCZ.cpp
 *  \brief Get transition matrices and distributions for the Cane-Zebiak model.
 *   
 * Get transition matrices and distributions from a long time series
 * from the Cane-Zebiak model.
 * For configuration, the file cfg/transferCZ.cfg is parsed using libconfig C++ library.
 * In this case, several time series corresponding to different
 * noise realizations are used.
 * First read the observable and get its mean and standard deviation
 * used to adapt the grid.
 * A rectangular grid is used here.
 * A grid membership vector is calculated for each time series 
 * assigning to each realization a grid box.
 * Then, the membership matrix is calculated for a given lag
 * with a concatenation over the seeds.
 * The forward and backward transition matrices as well as the initial and final distributions
 * are calculated from the membership matrix.
 * Note that, since the transitions are calculated from long time series,
 * these distributions should be equal (to the exception of the first and last realizations).
 * Finally, the results are printed.
 */

using namespace libconfig;


// Declarations
/** \brief User defined function to get parameters from a cfg file using libconfig. */
int readConfig(const char *cfgFileNamePrefix);

/** \briefCount number of lines in file. */
size_t lineCount(FILE *fp);


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

// Paths
const char prefix[] = "zc_1eof_";
const char indexDir[] = "../data/observables/";
const size_t dimSeries = 5;

// Units
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


// Main program
int main(int argc, char * argv[])
{
  // Read configuration file
  if (readConfig("transferCZ")) {
    std::cerr << "Error reading config file " << argv[0] << ".cfg"
	      << std::endl;
    return(EXIT_FAILURE);
  }

  // Declarations
  // Simulation parameters
  const double sampFreq = 0.35 / dt;
  size_t spinup = (size_t) (spinupMonth * sampFreq);

  // Names and Files
  char obsName[256], srcPostfix[256], srcPostfixSeed[256], postfix[256],
    gridPostfix[256], gridCFG[256], gridPostfixSeed[256], cpyBuffer[256];
  char gridMemFileName[256],
    forwardTransitionFileName[256], backwardTransitionFileName[256],
    initDistFileName[256], finalDistFileName[256];
  char gridFileName[256];
  FILE *gridMemFile;
  FILE *indexFile[DIM];
  char indexPath[DIM][256];
  size_t count;
  
  // Vectors and matrices
  size_t tauStep;
  double tauDim;
  std::vector<gsl_vector_uint *> gridMemSeeds(nSeeds);
  gsl_matrix_uint *gridMemMatrix;
  transferOperator *transferOp;
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
  Grid *grid;

  // Define grid name
  sprintf(srcPostfix, "%smu%04d_eps%04d", prefix, (int) (mu*1000),
	  (int) (eps*1000));
  sprintf(gridCFG, "");
  for (size_t d = 0; d < DIM; d++) {
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
			 + gsl_pow_2(gsl_matrix_get(statesSeeds[seed], k, d)));
	}
	gsl_matrix_free(data);
      }
    }

    // Create grid
    // Define grid limits
    xmin = gsl_vector_alloc(DIM);
    xmax = gsl_vector_alloc(DIM);
    for (size_t d = 0; d < DIM; d++) {
      gsl_vector_set(statesMean, d, gsl_vector_get(statesMean, d) / ntTot);
      gsl_vector_set(statesSTD, d,
		     sqrt(gsl_vector_get(statesSTD, d) / ntTot
			  - gsl_pow_2(gsl_vector_get(statesMean, d))));
      gsl_vector_set(xmin, d, gsl_vector_get(statesMean, d)
		     - gsl_vector_get(nSTDLow, d)
		     * gsl_vector_get(statesSTD, d));
      gsl_vector_set(xmax, d, gsl_vector_get(statesMean, d)
		     + gsl_vector_get(nSTDHigh, d)
		     * gsl_vector_get(statesSTD, d));
    }
    // Define grid
    sprintf(gridPostfix, "_%s%s%s", srcPostfix, obsName, gridCFG);
    sprintf(gridFileName, "../results/grid/grid%s.txt", gridPostfix);
    grid = new Grid(nx, xmin, xmax);
    // Print grid
    grid->printGrid(gridFileName, "%.12lf", true);

    // Get grid membership for each seed
    for (size_t seed = 0; seed < nSeeds; seed++){
      // Open grid membership file
      sprintf(srcPostfixSeed, "%s_seed%d", srcPostfix, (int) seed);
      sprintf(gridPostfixSeed, "_%s%s%s", srcPostfixSeed, obsName, gridCFG);
      sprintf(gridMemFileName, "../results/transitionMatrix/gridMem%s.txt",
	      gridPostfixSeed);
      if ((gridMemFile = fopen(gridMemFileName, "w")) == NULL){
	fprintf(stderr, "Can't open %s for writing!\n", gridMemFileName);
	return(EXIT_FAILURE);
      }

      // Get grid membership
      std::cout << "Getting grid membership vector for seed" << seed
		<< std::endl;
      gridMemSeeds[seed] = getGridMemVector(statesSeeds[seed], grid);

      // Write grid membership
      gsl_vector_uint_fprintf(gridMemFile, gridMemSeeds[seed], "%d");
      fclose(gridMemFile);
    }

    // Free
    gsl_vector_free(xmin);
    gsl_vector_free(xmax);
    for (size_t seed = 0; seed < nSeeds; seed++){
      gsl_matrix_free(statesSeeds[seed]);
    }
  }
  else {
    // Create empty grid
    grid = new Grid(nx);
      
    // Read membership vectors
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

  // Get transition matrices for different lags
  for (size_t lag = 0; lag < nLags; lag++){
    tauDim = gsl_vector_get(tauDimRng, lag);
    tauStep = (size_t) (tauDim * sampFreq);

    // Update file names
    sprintf(postfix, "%s_tau%03d", gridPostfix, (int) (tauDim * 1000));
    sprintf(forwardTransitionFileName,
	    "../results/transitionMatrix/forwardTransition%s.coo", postfix);
    sprintf(backwardTransitionFileName,
	    "../results/transitionMatrix/backwardTransition%s.coo", postfix);
    sprintf(initDistFileName, "../results/transitionMatrix/initDist%s.txt", postfix);
    sprintf(finalDistFileName, "../results/transitionMatrix/finalDist%s.txt", postfix);

    // Get full membership matrix
    std::cout << "Getting full membership matrix from the list of membership vecotrs..."
	      << std::endl;
    gridMemMatrix = memVectorList2memMatrix(&gridMemSeeds, tauStep);

    // Get transition matrices as CSR
    std::cout << "Building transfer operator..." << std::endl;
    transferOp = new transferOperator(gridMemMatrix, grid->N);
    
    // Write transition matrix as CSR
    std::cout << "Writing transfer operator..." << std::endl;
    transferOp->printForwardTransition(forwardTransitionFileName, "%.12lf");
    transferOp->printBackwardTransition(backwardTransitionFileName, "%.12lf");
    transferOp->printInitDist(initDistFileName, "%.12lf");
    transferOp->printFinalDist(finalDistFileName, "%.12lf");
	
    // Free
    delete transferOp;
    gsl_matrix_uint_free(gridMemMatrix);
  }

  // Free
  delete grid;
  gsl_vector_uint_free(nx);
  gsl_vector_free(nSTDLow);
  gsl_vector_free(nSTDHigh);
  gsl_vector_free(tauDimRng);
  gsl_vector_uint_free(seedRng);
  gsl_vector_uint_free(ntSeeds);
		
  return 0;
}

// Definitions
int
readConfig(const char *cfgFileNamePrefix)
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
    std::cout << std::endl << "---simulation---" << std::endl
  	      << "dt: " << dt << std::endl
  	      << "spinupMonth: " << spinupMonth << std::endl;

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


/**
 * \brief Count the number of lines in a file.
 *
 * Count the number of lines in a file.
 * \param[in] fp File from which to count lines.
 * \return Number of lines in file.
 */
size_t
lineCount(FILE *fp)
{
  size_t count = 0;
  int ch;

  // Count lines
  do {
    ch = fgetc(fp);
    if (ch == '\n')
      count++;
  } while (ch != EOF);

  return count;
}

