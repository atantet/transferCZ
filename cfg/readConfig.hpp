#ifndef READ_CONFIG_HPP
#define READ_CONFIG_HPP
#define DIM 2

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <libconfig.h++>
#include <configAR.hpp>
#include <ergoParam.hpp>

using namespace libconfig;


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

  fieldDef(){}
  fieldDef(const int idx, const char *name, const double factor)
    : tableIdx(idx), toDimFactor(factor) { strcpy(shortName, name); }
  fieldDef(const fieldDef &field)
    : tableIdx(field.tableIdx), toDimFactor(field.toDimFactor) { strcpy(shortName, field.shortName); }
};

// Configuration
// General
extern char indicesDir[256];
extern char resDir[256];
extern char specDir[256];
extern char plotDir[256];
extern char fileFormat[256];
// Units
extern double L;                   //!< Horizontal unit of length
extern double c0;                  //!< Ocean wave celerity
extern double timeDim;             //!< Time dimension in days
extern double sampFreq;            //!< Sampling frequency in months
extern double H;                   //!< Vertical unit of length
extern double tau_0;               //!< Wind-stress unit
extern double delta_T;             //!< Temperature unit
// CaseDef
extern char prefix[256];
extern char simType[256];
extern std::string indicesName[DIM];
extern fieldDef field_h;   //!< Thermocline depth record definition
extern fieldDef field_T;   //!< SST record definition
extern fieldDef fields[DIM];
extern double mu;
extern double eps;
extern gsl_vector_uint *seedRng;
extern size_t nSeeds;
extern char obsName[256];
extern char srcPostfix[256];
// Simulation
extern double dt;
extern double spinupMonth;
extern int dimSeries;
// Grid
extern bool readGridMem;
extern size_t N;
extern gsl_vector_uint *nx;
extern gsl_vector *nSTDLow;
extern gsl_vector *nSTDHigh;
extern char gridPostfix[256];
extern char gridCFG[256];
extern char gridFileName[256];
// Transfer
extern int nLags;
extern gsl_vector *tauDimRng;
// Spectrum
extern int nev;                        //!< Number of eigenvalues to look for
extern configAR config;                //!< Configuration data for the eigen problem
extern bool getForwardEigenvectors;    //!< Whether to get forward eigenvectors
extern bool getBackwardEigenvectors;   //!< Whether to get backward eigenvectors
extern bool makeBiorthonormal;         //!< Whether to make eigenvectors biorthonormal
char methodCfg[] = "LM";


// Definitions
void
readConfig(const char *configFileName)
{
  Config cfg;
  char cpyBuffer[256];
  configAR defaultCfgAR = {methodCfg, 0, 0., 0, NULL, true};
  
  
  // Read the file. If there is an error, report it and exit.
  try {
    std::cout << "Reading config file " << configFileName << std::endl;
    cfg.readFile(configFileName);
    
    std::cout.precision(6);
    std::cout << "Settings:" << std::endl;

    // Get general settings
    std::cout << std::endl << "---general---" << std::endl;
    strcpy(indicesDir, (const char *) cfg.lookup("general.indicesDir"));
    std::cout << "Indices directory: " << indicesDir << std::endl;
    strcpy(resDir, (const char *) cfg.lookup("general.resDir"));
    std::cout << "Results directory: " << resDir << std::endl;
    strcpy(specDir, (const char *) cfg.lookup("general.specDir"));
    std::cout << "Spectrum results directory: " << specDir << std::endl;
    strcpy(plotDir, (const char *) cfg.lookup("general.plotDir"));
    std::cout << "Plots directory: " << plotDir << std::endl;
    strcpy(fileFormat, (const char *) cfg.lookup("general.fileFormat"));
    std::cout << "Files format: " << fileFormat << std::endl;

    
    // Get units
    std::cout << std::endl << "---units---" << std::endl;
    L = cfg.lookup("units.L");
    std::cout << "L = " << L << std::endl;
    c0 = cfg.lookup("units.c0");
    std::cout << "c0 = " << c0 << std::endl;
    H = cfg.lookup("units.H");
    std::cout << "H = " << H << std::endl;
    tau_0 = cfg.lookup("units.tau_0");
    std::cout << "tau_0 = " << tau_0 << std::endl;
    delta_T = cfg.lookup("units.delta_T");
    std::cout << "delta_T = " << delta_T << std::endl;
    timeDim = L / c0 / (60. * 60. * 24); 
    std::cout << "Time dimension (days): " << timeDim << std::endl;
    field_h = fieldDef(1, "h", H);
    field_T = fieldDef(2, "T", delta_T);

    // Get caseDef settings
    std::cout << std::endl << "---caseDef---" << std::endl;
    // Indices and Fields
    strcpy(prefix, (const char *) cfg.lookup("caseDef.prefix"));
    strcpy(simType, (const char *) cfg.lookup("caseDef.simType"));
    const Setting &indicesNameSetting = cfg.lookup("caseDef.indicesName");
    const Setting &fieldsNameSetting = cfg.lookup("caseDef.fieldsName");
    for (size_t d = 0; d < DIM; d++) {
      indicesName[d] = (const char *) indicesNameSetting[d];
      if (!strcmp(fieldsNameSetting[d], "h"))
	fields[d] = field_h;
      else if (!strcmp(fieldsNameSetting[d], "T"))
	fields[d] = field_T;
      else {
	std::cerr << "Wrong field name for dimension " << d << std::endl;
	throw std::exception();
      }	
      std::cout << "indicesName[" << d << "]: "
		<< indicesName[d] << std::endl;
      std::cout << "fieldsName[" << d << "]: "
		<< fields[d].shortName << std::endl;
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
    sampFreq = 1. / timeDim * 30 / dt;
    spinupMonth = cfg.lookup("simulation.spinupMonth");
    dimSeries = cfg.lookup("simulation.dimSeries");
    std::cout << std::endl << "---simulation---" << std::endl
  	      << "dt: " << dt << std::endl
  	      << "spinupMonth: " << spinupMonth << std::endl
  	      << "dimSeries: " << dimSeries << std::endl
	      << "Sampling frequency (1/months): " << sampFreq << std::endl;

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
    std::cout << std::endl << "---transfer---" << std::endl;
    if (cfg.exists("transfer.tauDimRng"))
      {
	const Setting &tauDimRngSetting = cfg.lookup("transfer.tauDimRng");
	nLags = tauDimRngSetting.getLength();
	tauDimRng = gsl_vector_alloc((size_t) nLags);
	for (size_t lag = 0; lag < (size_t) nLags; lag++)
	  gsl_vector_set(tauDimRng, lag, tauDimRngSetting[lag]);
      }
    else if (cfg.exists("transfer.nLags")
	     && cfg.exists("transfer.stepLag")
	     && cfg.exists("transfer.lag0"))
      {
	double stepLag = cfg.lookup("transfer.stepLag");
	double lag0 = cfg.lookup("transfer.lag0");
	nLags = cfg.lookup("transfer.nLags");
	tauDimRng = gsl_vector_alloc((size_t) nLags);
	for (size_t lag = 0; lag < (size_t) nLags; lag++)
	  gsl_vector_set(tauDimRng, lag, lag0 + lag*stepLag);
      }
    else
      {
	std::cerr << "Lags not properly configured." << std::endl;
	throw std::exception();
      }
    std::cout << "tauDimRng = {";
    for (size_t lag = 0; lag < (size_t) nLags; lag++)
      {
	std::cout << gsl_vector_get(tauDimRng, lag) << ", ";
      }
    std::cout << "}" << std::endl;

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

    
    // Define grid name
    sprintf(srcPostfix, "%s%s_mu%04d_eps%04d", prefix, simType,
	    (int) (mu*1000), (int) (eps*1000));
    sprintf(gridCFG, "");
    strcpy(obsName, "");
    for (size_t d = 0; d < DIM; d++)
      {
	strcpy(cpyBuffer, obsName);
	sprintf(obsName, "%s_%s_%s", cpyBuffer, fields[d].shortName,
		indicesName[d].c_str());
	strcpy(cpyBuffer, gridCFG);
	sprintf(gridCFG, "%s_n%dl%dh%d", cpyBuffer,
		gsl_vector_uint_get(nx, d),
		(int) gsl_vector_get(nSTDLow, d),
		(int) gsl_vector_get(nSTDHigh, d));
      }
    sprintf(gridPostfix, "_%s%s%s", srcPostfix, obsName, gridCFG);
    sprintf(gridFileName, "%s/grid/grid%s.txt", resDir, gridPostfix);

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
  
  return;
}


/**
 * Free memory allocated during configuration.
 */
void
freeConfig()
{
  gsl_vector_uint_free(nx);
  gsl_vector_free(nSTDLow);
  gsl_vector_free(nSTDHigh);
  gsl_vector_free(tauDimRng);
  gsl_vector_uint_free(seedRng);

  return;
}

#endif
