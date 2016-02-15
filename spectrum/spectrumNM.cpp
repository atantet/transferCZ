#define DIM 2
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <libconfig.h++>
#include <ergoPack/transferOperator.hpp>

using namespace libconfig;


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
 * The forward transition matrix as well as the initial distribution
 * are calculated from the membership matrix.
 * Note that, since the transitions are calculated from long time series,
 * the problem must be stationary and there is no need to calculate
 * the backward transition matrix and final distribution.
 * Finally, the results are printed.
 */


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
const char resDir[] = "../results/";
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
char srcPostfix[256], gridCFG[256], obsName[256], gridPostfix[256], gridFileName[256];
std::string indicesName[DIM];
struct fieldDef fieldsDef[DIM];
double mu, eps;
gsl_vector_uint *seedRng;
size_t nSeeds;
bool readGridMem;
size_t N;
double dt, spinupMonth;
gsl_vector_uint *nx;
gsl_vector *nSTDLow, *nSTDHigh;
size_t nLags;
gsl_vector *tauDimRng;
int nev;


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
  char srcPostfixSeed[256], postfix[256], gridPostfixSeed[256];
  char gridMemFileName[256], forwardTransitionFileName[256], initDistFileName[256];
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

  // Read membership vectors
  for (size_t seed = 0; seed < nSeeds; seed++) {
    // Open grid membership file to read
    sprintf(srcPostfixSeed, "%s_seed%d", srcPostfix, (int) seed);
    sprintf(gridPostfixSeed, "_%s%s%s", srcPostfixSeed, obsName, gridCFG);
    sprintf(gridMemFileName, "%s/transitionMatrix/gridMem%s.txt",
	    resDir, gridPostfixSeed);
    std::cout << "Reading grid membership vector for seed " << seed
	      << " at " << gridMemFileName << std::endl;
    if ((gridMemFile = fopen(gridMemFileName, "r")) == NULL){
      fprintf(stderr, "Can't open %s for reading!\n", gridMemFileName);
      return(EXIT_FAILURE);
    }

    // Allocate grid membership vector for seed
    count = lineCount(gridMemFile);
    fseek(gridMemFile, 0, SEEK_SET);
    gridMemSeeds[seed] = gsl_vector_uint_alloc(count);

    gsl_vector_uint_fscanf(gridMemFile, gridMemSeeds[seed]);
    fclose(gridMemFile);
  }

  // Get transition matrices for different lags
  double *EigValReal = new double [nev];
  double *EigValImag = new double [nev];

  size_t count = 0;
  for (size_t ev = 0; ev < (size_t) nev; ev++)
    {
      for (size_t lag = 0; lag < nLags; lag++)
	{
	  tauDim = gsl_vector_get(tauDimRng, lag);
	  tauStep = (size_t) (tauDim * sampFreq);

	  // Update file names
	  sprintf(postfix, "%s_tau%03d", gridPostfix, (int) (tauDim * 1000));
	  sprintf(forwardTransitionFileName,
		  "%s/transitionMatrix/forwardTransition%s.coo", resDir, postfix);
	  sprintf(initDistFileName, "%s/transitionMatrix/initDist%s.txt", resDir, postfix);

	  // Get full membership matrix
	  std::cout << "Getting full membership matrix from the list of membership vecotrs..."
		    << std::endl;
	  gridMemMatrix = memVectorList2memMatrix(&gridMemSeeds, tauStep);

	  // Get transition matrices as CSR
	  std::cout << "Building transfer operator..." << std::endl;
	  transferOp = new transferOperator(gridMemMatrix, N);

	  // Solve eigen value problem with default configuration
	  transferSpec = new transferSpectrum(1, transferOp, config);
	  std::cout << "Solving eigen problem for forward transition matrix..." << std::endl;
	  transferSpec->getSpectrumForward();
	  std::cout << "Found " << transferSpec->EigProbForward.ConvergedEigenvalues()
		    << " eigenvalues." << std::endl;

	  // Get condition number
	  transferSpec->getConditionNumber();

	  if (

	  // Save eigenvalue
	  if ((lag == 0) || (transferSpec->EigValForwardReal[0] > EigValReal[count]))
	    {
	      EigValReal[count] = transferSpec->EigValForwardReal[0];
	      EigValImag[count] = transferSpec->EigValForwardImag[0];
	      
	      if (transferSpec->EigProbForward.ConvergedEigenvalues() == 2)
		{
		  EigValReal[count+1] = transferSpec->EigValForwardReal[0];
		  EigValImag[count+1] = transferSpec->EigValForwardImag[0];
		}
	    }
	  // std::cout << "Solving eigen problem for backward transition matrix..." << std::endl;
	  // transferSpec->getSpectrumBackward();
	  // std::cout << "Found " << transferSpec->EigProbBackward.ConvergedEigenvalues()
	  // 		<< " eigenvalues." << std::endl;
    
	  // Free
	  delete transferOp;
	  gsl_matrix_uint_free(gridMemMatrix);
	}
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
int
readConfig(const char *cfgFileNamePrefix)
{
  Config cfg;
  char cfgFileName[256], cpyBuffer[256];
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
    N = 1;
    nSTDLow = gsl_vector_alloc(DIM);
    nSTDHigh = gsl_vector_alloc(DIM);
    for (size_t d = 0; d < DIM; d++) {
      gsl_vector_uint_set(nx, d, nxSetting[d]);
      N *= gsl_vector_uint_get(nx, d);
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

  // Define grid name
  sprintf(srcPostfix, "%smu%04d_eps%04d", prefix, (int) (mu*1000),
	  (int) (eps*1000));
  sprintf(gridCFG, "");
  strcpy(obsName, "");
  for (size_t d = 0; d < DIM; d++)
    {
      strcpy(cpyBuffer, obsName);
      sprintf(obsName, "%s_%s_%s", cpyBuffer, fieldsDef[d].shortName,
	      indicesName[d].c_str());
      strcpy(cpyBuffer, gridCFG);
      sprintf(gridCFG, "%s_n%dl%dh%d", cpyBuffer,
	      gsl_vector_uint_get(nx, d),
	      (int) gsl_vector_get(nSTDLow, d),
	      (int) gsl_vector_get(nSTDHigh, d));
    }

  sprintf(gridPostfix, "_%s%s%s", srcPostfix, obsName, gridCFG);
  sprintf(gridFileName, "%s/grid/grid%s.txt", resDir, gridPostfix);

  // Get spectrum setting 
  nev = cfg.lookup("spectrum.nev");
  std::cout << std::endl << "---spectrum---" << std::endl
	    << "nev: " << nev << std::endl;

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

