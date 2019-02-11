#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <transferOperator.hpp>
#include "../cfg/readConfig.hpp"

/** \file transferCZ.cpp
 *  \brief Get transition matrices and for the Cane-Zebiak model.
 *   
 * Get transition matrices and distributions from a long time series
 * from the Cane-Zebiak model.
 * For configuration, the file cfg/transferCZ.cfg
 * is parsed using libconfig C++ library.
 * In this case, several time series corresponding to different
 * noise realizations are used.
 * First read the observable and get its mean and standard deviation
 * used to adapt the grid.
 * A rectangular grid is used here.
 * A grid membership vector is calculated for each time series 
 * assigning to each realization a grid box.
 * Then, the membership matrix is calculated for a given lag
 * with a concatenation over the seeds.
 * The transition matrix and the initial distribution
 * are calculated from the membership matrix.
 * Note that, since the transitions are calculated from long time series,
 * these distributions should be equal
 * (to the exception of the first and last realizations).
 * Finally, the results are printed.
 */


/** \briefCount number of lines in file. */
size_t lineCount(FILE *fp);

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
bool makeBiorthonormal;       //!< Whether to make eigenvectors biorthonormal


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
  catch (...)  {
    std::cerr << "Error reading configuration file" << std::endl;
    return(EXIT_FAILURE);
  }

  // Declarations
  // Simulation parameters
  size_t spinup = (size_t) (spinupMonth * sampFreq);

  // Names and Files
  char srcPostfixSeed[256], dstPostfix[256], dstPostfixTau[256],
    gridPostfixSeed[256];
  char gridMemFileName[256], transitionFileName[256],
    maskFileName[256], initDistFileName[256];
  FILE *gridMemStream;
  FILE *indexFile;
  char indexPath[256];
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
  size_t ntIndex;
  size_t d;
  gsl_vector *statesMean = gsl_vector_calloc(DIM);
  gsl_vector *statesSTD = gsl_vector_calloc(DIM);
  std::vector<gsl_matrix *> statesSeeds(nSeeds);
  gsl_vector_uint *ntSeeds = gsl_vector_uint_alloc(nSeeds);
  size_t ntTot = 0;

  // Grid related
  Grid *grid;

  // Get membership vector
  if (! readGridMem) {
    for (size_t seed = 0; seed < nSeeds; seed++) {
      
      // Read observable
      sprintf(srcPostfixSeed, "%s_seed%d", srcPostfix, (int) seed);
      sprintf(indexPath, "%s/%s/indices.txt", indicesDir, srcPostfixSeed);
      // Open observable file d
      if ((indexFile = fopen(indexPath, "r")) == NULL){
	fprintf(stderr, "Can't open %s for reading!\n", indexPath);
	return(EXIT_FAILURE);
      }
      // Get number of lines
      ntIndex = lineCount(indexFile);
      fseek(indexFile, 0, SEEK_SET);

      // Get index length and allocate
      gsl_vector_uint_set(ntSeeds, seed, ntIndex - spinup);
      ntTot += gsl_vector_uint_get(ntSeeds, seed);
      statesSeeds[seed] = gsl_matrix_alloc(gsl_vector_uint_get(ntSeeds, seed),
					   DIM);

      // Read trajectory and get the total mean and standard defination
      data = gsl_matrix_alloc(ntIndex, dimSeries);
      // Read trajectory
      std::cout << "Reading trajectory in " << indexPath << std::endl;
      gsl_matrix_fscanf(indexFile, data);
      fclose(indexFile);
      for (size_t k = 0; k < gsl_vector_uint_get(ntSeeds, seed); k++) {
	// Eastern SST
	d = 0;
	gsl_matrix_set(statesSeeds[seed], k, d,
		       gsl_matrix_get(data, spinup+k, 1) * delta_T);
	gsl_vector_set(statesMean, d, gsl_vector_get(statesMean, d)
		       + gsl_matrix_get(statesSeeds[seed], k, d));
	gsl_vector_set(statesSTD, d, gsl_vector_get(statesSTD, d)
		       + gsl_pow_2(gsl_matrix_get(statesSeeds[seed],
						  k, d)));

	// Western TD
	d = 1;
	gsl_matrix_set(statesSeeds[seed], k, d,
		       gsl_matrix_get(data, spinup+k, 2) * H);
	gsl_vector_set(statesMean, d, gsl_vector_get(statesMean, d)
		       + gsl_matrix_get(statesSeeds[seed], k, d));
	gsl_vector_set(statesSTD, d, gsl_vector_get(statesSTD, d)
		       + gsl_pow_2(gsl_matrix_get(statesSeeds[seed],
						  k, d)));
      }
      gsl_matrix_free(data);
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
    grid = new RegularGrid(nx, xmin, xmax);
    // Print grid
    grid->printGrid(gridFileName, "%.12lf", true);
    sprintf(dstPostfix, "%s", gridPostfix);

    // Get grid membership for each seed
    for (size_t seed = 0; seed < nSeeds; seed++){
      // Open grid membership file
      sprintf(srcPostfixSeed, "%s_seed%d", srcPostfix, (int) seed);
      sprintf(gridPostfixSeed, "_%s%s%s", srcPostfixSeed, obsName, gridCFG);
      sprintf(gridMemFileName, "%s/transfer/gridMem/gridMem%s.%s",
	      resDir, gridPostfixSeed, fileFormat);
      if ((gridMemStream = fopen(gridMemFileName, "w")) == NULL){
	fprintf(stderr, "Can't open %s for writing!\n", gridMemFileName);
	return(EXIT_FAILURE);
      }

      // Get grid membership
      std::cout << "Getting grid membership vector for seed" << seed
		<< std::endl;
      gridMemSeeds[seed] = getGridMemVector(statesSeeds[seed], grid);

      // Write grid membership
      if (strcmp(fileFormat, "bin") == 0)
	gsl_vector_uint_fwrite(gridMemStream, gridMemSeeds[seed]);
      else
	gsl_vector_uint_fprintf(gridMemStream, gridMemSeeds[seed], "%d");
      fclose(gridMemStream);
    }

    // Free
    gsl_vector_free(xmin);
    gsl_vector_free(xmax);
    for (size_t seed = 0; seed < nSeeds; seed++)
      gsl_matrix_free(statesSeeds[seed]);
    delete grid;
  }
  else {
    // Read membership vectors
    for (size_t seed = 0; seed < nSeeds; seed++) {
      // Open grid membership file to read
      sprintf(srcPostfixSeed, "%s_seed%d", srcPostfix, (int) seed);
      sprintf(gridPostfixSeed, "_%s%s%s", srcPostfixSeed, obsName, gridCFG);
      sprintf(gridMemFileName, "%s/transfer/gridMem/gridMem%s.%s",
	      resDir, gridPostfixSeed, fileFormat);
      std::cout << "Reading grid membership vector for seed " << seed
		<< " at " << gridMemFileName << std::endl;
      if ((gridMemStream = fopen(gridMemFileName, "r")) == NULL){
	fprintf(stderr, "Can't open %s for reading!\n", gridMemFileName);
	return(EXIT_FAILURE);
      }

      // Allocate grid membership vector for seed
      count = lineCount(gridMemStream);
      fseek(gridMemStream, 0, SEEK_SET);
      gridMemSeeds[seed] = gsl_vector_uint_alloc(count);

      if (strcmp(fileFormat, "bin") == 0)
	gsl_vector_uint_fread(gridMemStream, gridMemSeeds[seed]);
      else
	gsl_vector_uint_fscanf(gridMemStream, gridMemSeeds[seed]);
      fclose(gridMemStream);
    }
  }

  // Get transition matrices for different lags
  for (size_t lag = 0; lag < (size_t) nLags; lag++) {
    tauDim = gsl_vector_get(tauDimRng, lag);
    tauStep = (size_t) (tauDim * sampFreq);
    sprintf(dstPostfixTau, "%s_tau%03d", dstPostfix,
	    (int) (tauDim * 1000 + 0.1));
    std::cout << "Getting transfer operator for lag "
	      << tauDim << std::endl;

      
    // Get full membership matrix
    std::cout << "Getting full membership matrix from the list of membership vectors..."
	      << std::endl;
    gridMemMatrix = memVectorList2memMatrix(&gridMemSeeds, tauStep);

      
    // Get transition matrices as CSR
    std::cout << "Building transfer operator..." << std::endl;
    transferOp = new transferOperator(gridMemMatrix, N);


    // Wrute results
    // Write transition matrix
    std::cout
      << "Writing transition matrix and initial distribution..."
      << std::endl;
    sprintf(transitionFileName,
	    "%s/transfer/transition/transition%s.crs%s",
	    resDir, dstPostfixTau, fileFormat);
    transferOp->printTransition(transitionFileName, fileFormat, "%.12lf");
      
    // Write mask and initial distribution
    if (lag == 0) {
      sprintf(maskFileName, "%s/transfer/mask/mask%s.%s",
	      resDir, dstPostfix, fileFormat);
      transferOp->printMask(maskFileName, fileFormat, "%.12lf");

      sprintf(initDistFileName, "%s/transfer/initDist/initDist%s.%s",
	      resDir, dstPostfix, fileFormat);
      transferOp->printInitDist(initDistFileName, fileFormat, "%.12lf");
    }
      
	
    // Free
    delete transferOp;
    gsl_matrix_uint_free(gridMemMatrix);
  }
  
  // Free
  freeConfig();
  gsl_vector_uint_free(ntSeeds);
		
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

