#include "PEG.h"

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cfloat>
#include <climits>
#include <vector>
#include <cmath>

// h*c (Planck constant * speed of light), in eV * um.
#define M_HC 1.23984172

/*! \todo

- For the following input, beta2_n is coming out weird (all imaginary, instead of all real) for the outside orders.

./pegSerial --mode constantIncidence --min 100 --max 120 --increment 5 --incidenceAngle 88 --outputFile testOutput.txt --progressFile testProgress.txt --gratingType blazed --gratingPeriod 1 --printDebugOutput --gratingMaterial Au --N 15 --gratingGeometry 2.5,30 --eV

Maybe because at this incidence, there are no outside orders? Just evanescent waves?

*/

// Command-line option variables
////////////////////////////////
enum Mode {InvalidMode, ConstantIncidence, ConstantIncludedAngle, ConstantWavelength};
Mode iMode = InvalidMode;

double iMin = DBL_MAX, 
		iMax = DBL_MAX, 
		iIncrement = DBL_MAX, 
		iIncidenceAngle = DBL_MAX, 
		iIncludedAngle = DBL_MAX, 
		iWavelength = DBL_MAX;
int iToOrder = INT_MAX;

int iN = INT_MAX;

PEGrating::Profile iProfile = PEGrating::InvalidProfile;
double iPeriod = DBL_MAX;
double iGeometry[8] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
std::string iMaterial;

std::string iOutputFile, iProgressFile;

bool iEv = false;
bool iPrintDebugOutput = false;
//////////////////////////////////


// Helper Functions
////////////////////////////////

/// This helper function parses the command-line options and sets the global input variables above.  Returns true on success, and false on failure.
bool parseCommandLineOptions(int argc, char** argv);

/// This helper function writes the header to the output file stream
void writeOutputFileHeader(std::ostream& outputFileStream);

/// This helper function appends a single efficiency result to the output file stream
void writeOutputFileResult(std::ostream& outputFileStream, const PEResult& result);

/// This helper function appends the progress description to the given output stream
void writeOutputFileProgress(std::ostream& outputFileStream, int completedSteps, int totalSteps, bool anySuccesses, bool anyFailures);

////////////////////////////////


/// This main program provides a command-line interface to run a series of sequential grating efficiency calculations. The results are written to an output file, and (optionally) a second file is written to provide information on the status of the calculation.  [This file is only responsible for input processing and output; all numerical details are structured within PEGrating and PESolver.]
/*! 
<b>Command-line options</b>

<b>Required:</b>

Grating specification:

--gratingType <rectangular|blazed|sinusoidal|trapezoidal>
--gratingPeriod <grating period in um>
--gratingGeometry <command-delimited list of geometry parameters, in um and/or degrees>
	TODO
--gratingMaterial <grating substrate material>
	This should be a name corresponding to a refractive index database filename, ex: Au, Ni, C, SiO2, etc.
	
--N <truncation index>
	Specifies the number of positive and negative orders to include in the Fourier expansion. Will also determine the number of orders that are calculated, although if you only want to calculate 3 orders, you will still need a much larger truncation index for accurate results.  In the soft x-ray range, convergence is usually attained with N ~ 15..45.

Operating mode:

--mode <constantIncidence|constantIncludedAngle|constantWavelength>
--min <min>
--max <max>
--increment <increment>

[Required, depending on the \c mode]

--incidenceAngle <incidence angle in degrees>
--includedAngle <deviation angle in degrees> --toOrder <diffraction order for the included angle>
--wavelength <wavelength in um>

	In constant incidence mode, a calculation is performed for wavelengths from --min to --max in steps of --increment, at a fixed incidence angle given by --incidenceAngle.
	In constant included angle mode, the incidence angle is calculated at each wavelength to ensure a constant included angle of --includedAngle between the incident light and the order specified in --toOrder. This is the operating mode for many monochromators. (Inside orders are negative, outside orders are positive.)
	In constant wavelength mode, a calculation is performed for incidence angles from --min to --max in steps of --increment, for a fixed wavelength given by --wavelength.
	
Output:

--outputFile <file name>
	The calculation output will be written to this file.

<b>Optional:</b>

--progressFile <file name>
	If provided, the current status of the calculation will be written in this file; it can be monitored to determine the progress of long calculations.  This provides an interface for other processes to monitor the status of this calculation (for example, a web-based or GUI front-end, etc.).
	
--eV
	If this flag is included, all wavelength inputs (--min, --max, --increment, and --wavelength) will instead be interpreted as photon energies in electron volts (eV).
	
--printDebugOutput
	If this flag is included, each calculation will print intermediate results to standard output.
	
	
	
<b>Output</b>

An example of the output file written to --outputFile is shown below. If the file exists already, it will be overwritten.

=========================
# Input
mode=constantIncidence
incidenceAngle=88
units=eV
min=100
max=300
increment=5
gratingType=blazed
gratingPeriod=1.6
gratingGeometry=3.2,30.0
gratingMaterial=Au
N=5
# Progress
status=succeeded     (inProgress, someFailed, allFailed, succeeded)
completedSteps=41
totalSteps=41
# Output
100[tab]<e-5>,<e-4>,<e-3>,<e-2>,<e-1>,<e0>,<e1>,<e2>,<e3>,<e4>,<e5>
105[tab]<e-5>,<e-4>,<e-3>,<e-2>,<e-1>,<e0>,<e1>,<e2>,<e3>,<e4>,<e5>
110[tab]<e-5>,<e-4>,<e-3>,<e-2>,<e-1>,<e0>,<e1>,<e2>,<e3>,<e4>,<e5>
...
=========================

If a --progressFile is specified, it is written and re-written during the calculation, containing the # Progress section from the main output file:

=========================
# Progress
status=inProgress
completedSteps=3
totalSteps=41
=========================

*/
int main(int argc, char** argv) {
	
	// parse command line options into input variables.
	if(!parseCommandLineOptions(argc, argv)) {
		return -1;
	}
	
	
	// Open the output file:
	std::ofstream outputFile(iOutputFile.c_str(), std::ios::out | std::ios::trunc);
	if(!outputFile.is_open()) {
		std::cerr << "Could not open output file " << iOutputFile;
		return -1;
	}
	
	// Open the progress file, if provided:
	std::ofstream progressFile;
	if(!iProgressFile.empty()) {
		progressFile.open(iProgressFile.c_str(), std::ios::out | std::ios::trunc);
		if(!progressFile.is_open()) {
			std::cerr << "Could not open progress file " << iProgressFile;
			return -1;
		}
	}
	
	// Write the file header:
	writeOutputFileHeader(outputFile);
	// Remember this position in the output file; it is where we will write the progress and output lines
	std::streampos outputFilePosition = outputFile.tellp();
	
	// How many steps do we have?
	int totalSteps = int((iMax - iMin)/iIncrement) + 1;
	
	// Write the initial progress:
	writeOutputFileProgress(outputFile, 0, totalSteps, false, false);
	if(!iProgressFile.empty()) {
		std::ofstream progressFile(iProgressFile.c_str(), std::ios::out | std::ios::trunc);
		writeOutputFileProgress(progressFile, 0, totalSteps, false, false);
	}
	
	// create the grating object.
	PEGrating* grating;
	switch(iProfile) {
	case PEGrating::RectangularProfile:
		grating = new PERectangularGrating(iPeriod, iGeometry[0], iGeometry[1], iMaterial);
		break;
	case PEGrating::BlazedProfile:
		grating = new PEBlazedGrating(iPeriod, iGeometry[0], iGeometry[1], iMaterial);
		break;
	case PEGrating::SinusoidalProfile:
		grating = new PESinusoidalGrating(iPeriod, iGeometry[0], iMaterial);
		break;
	case PEGrating::TrapezoidalProfile:
		grating = new PETrapezoidalGrating(iPeriod, iGeometry[0], iGeometry[1], iGeometry[2], iGeometry[3], iMaterial);
		break;
	default:
		grating = 0;	// this will never happen; input validation assures one of the valid grating types.
		break;
	}
	
	// set math options: truncation index from input.
	PEMathOptions mathOptions(iN);
	
	// output data stored here:
	bool anyFailures = false;
	bool anySuccesses = false;
	std::vector<PEResult> results;
	
	// sequential loop over calculation steps
	for(int i=0; i<totalSteps; ++i) {
		
		double currentValue = iMin + iIncrement*i;
		
		// determine wavelength (um): depends on mode and eV/um setting.
		double wavelength = (iMode == ConstantWavelength) ? iWavelength : currentValue;
		if(iEv)
			wavelength = M_HC / wavelength;	// interpret input wavelength as eV instead, and convert to actual wavelength.  Formula: wavelength = hc / eV.     hc = 1.23984172 eV * um.
		
		// determine incidence angle: depends on mode and possibly wavelength.
		double incidenceAngle;
		switch(iMode) {
		case ConstantIncidence:
			incidenceAngle = iIncidenceAngle;
			break;
		case ConstantIncludedAngle: {
			double ciaRad = iIncludedAngle * M_PI / 180;
			incidenceAngle = (asin(-iToOrder*wavelength/2/iPeriod/cos(ciaRad/2)) + ciaRad/2) * 180 / M_PI;	// formula for constant included angle: satisfies alpha + beta = cia, and grating equation iToOrder*wavelength/d = sin(beta) - sin(alpha).
			break;
			}
		case ConstantWavelength:
			incidenceAngle = currentValue;
			break;
		default:
			incidenceAngle = 0; // never happens: input validation assures valid mode.
			break;
		}
		
		// run calculation
		PEResult result = grating->getEff(incidenceAngle, wavelength, mathOptions);
		if(result.status == PEResult::Success)
			anySuccesses = true;
		else
			anyFailures = true;
		results.push_back(result);
		
		// Print progress and results to output file.
		outputFile.seekp(outputFilePosition);
		writeOutputFileProgress(outputFile, i+1, totalSteps, anySuccesses, anyFailures);
		outputFile << "# Output" << std::endl;
		for(int j=0, cc=results.size(); j<cc; ++j)
			writeOutputFileResult(outputFile, results.at(j));
			
		// Update progress in progressFile, if provided.
		if(!iProgressFile.empty()) {
			std::ofstream progressFile(iProgressFile.c_str(), std::ios::out | std::ios::trunc);
			writeOutputFileProgress(progressFile, i+1, totalSteps, anySuccesses, anyFailures);
		}

	} // end of calculation loop.

	outputFile.close();
	return 0;
}


// This function parses the command-line options and sets the global input variables.  Returns true on success, and false on failure.
bool parseCommandLineOptions(int argc, char** argv) {
	// Parse Command-line options
	try {
		// loop over all input options
		while(1) {
			// define the options available
			static struct option long_options[] = {
				{"mode", required_argument, 0, 1},
				{"min", required_argument, 0, 2},
				{"max", required_argument, 0, 3},
				{"increment", required_argument, 0, 4},
				{"incidenceAngle", required_argument, 0, 5},
				{"includedAngle", required_argument, 0, 6},
				{"toOrder", required_argument, 0, 7},
				{"wavelength", required_argument, 0, 8},
				{"outputFile", required_argument, 0, 9},
				{"progressFile", required_argument, 0, 10},
				{"gratingType", required_argument, 0, 11},
				{"gratingGeometry", required_argument, 0, 12},
				{"gratingPeriod", required_argument, 0, 13},
				{"gratingMaterial", required_argument, 0, 14},
				{"N", required_argument, 0, 15},
				{"eV", no_argument, 0, 16},
				{"printDebugOutput", no_argument, 0, 17},
				{0, 0, 0, 0}
			};
				
			int option_index = 0;
			int c = getopt_long (argc, argv, "", long_options, &option_index);

			/* Detect the end of the options. */
			if (c == -1)
				break;

			switch (c) {
			case 1: // mode
				//if(!optarg) throw "An argument to --mode must be provided: constantIncidence, constantIncludedAngle, or constantWavelength.";
				if(strcmp(optarg, "constantIncidence") == 0) iMode = ConstantIncidence;
				else if(strcmp(optarg, "constantIncludedAngle") == 0) iMode = ConstantIncludedAngle;
				else if(strcmp(optarg, "constantWavelength") == 0) iMode = ConstantWavelength;
				else throw "The argument to --mode must be one of: constantIncidence, constantIncludedAngle, or constantWavelength.";
				break;
				
			case 2: // min
				//if(!optarg) throw "An argument to --min must be provided.";
				iMin = atof(optarg);
				break;
				
			case 3: // max
				//if(!optarg) throw "An argument to --max must be provided.";
				iMax = atof(optarg);
				break;
				
			case 4: // increment
				//if(!optarg) throw "An argument to --increment must be provided.";
				iIncrement = atof(optarg);
				break;
				
			case 5: // incidence angle
				//if(!optarg) throw "An argument to --incidenceAngle must be provided.";
				iIncidenceAngle = atof(optarg);
				break;
				
			case 6: // includedAngle
				//if(!optarg) throw "An argument to --includedAngle must be provided.";
				iIncludedAngle = atof(optarg);
				break;
				
			case 7: // toOrder
				//if(!optarg) throw "An argument to --toOrder must be provided.";
				iToOrder = atol(optarg);
				break;
				
			case 8: // wavelength
				//if(!optarg) throw "An argument to --wavelength must be provided.";
				iWavelength = atof(optarg);
				break;
				
			case 9: // outputFile
				//if(!optarg) throw "An argument to --outputFile must be provided.";
				iOutputFile = optarg;
				break;
				
			case 10: // progressFile
				//if(!optarg) throw "An argument to --progressFile must be provided.";
				iProgressFile = optarg;
				break;
				
			case 11: // gratingType
				//if(!optarg) throw "An argument to --gratingType must be provided: rectangular, blazed, sinusoidal, or trapezoidal.";
				if(strcmp(optarg, "rectangular") == 0) iProfile = PEGrating::RectangularProfile;
				else if(strcmp(optarg, "blazed") == 0) iProfile = PEGrating::BlazedProfile;
				else if(strcmp(optarg, "sinusoidal") == 0) iProfile = PEGrating::SinusoidalProfile;
				else if(strcmp(optarg, "trapezoidal") == 0) iProfile = PEGrating::TrapezoidalProfile;
				else throw "The argument to --gratingType must be one of: rectangular, blazed, sinusoidal, or trapezoidal.";
				break;
			case 12: { // gratingGeometry
				//if(!optarg) throw "An argument to --gratingGeometry must be provided, with a list of geometry parameters.";
				char* saveptr;
				char* token;
				int i=0;
				token = strtok_r(optarg, ",", &saveptr);
				while(token && i<8) {
					iGeometry[i++] = atof(token);
					token = strtok_r(0, ",", &saveptr);
				}
				break;
				}
				
			case 13: // gratingPeriod
				//if(!optarg) throw "An argument to --gratingPeriod must be provided.";
				iPeriod = atof(optarg);
				break;
				
			case 14: // gratingMaterial
				iMaterial = optarg;
				break;
			
			case 15: // N
				iN = atol(optarg);
				break;
				
			case 16: // interpret as eV instead of wavelength
				iEv = true;
				break;
			case 17: // print debug output to stdout
				iPrintDebugOutput = true;
				break;
			}
		} // end of loop over input options.
		
		// Check for missing/invalid input
		if(iMode == InvalidMode) throw "The operating mode --mode must be provided: constantIncidence, constantIncludedAngle, or constantWavelength";
		if(iMin == DBL_MAX) throw "The minimum range --min must be provided.";
		if(iMax == DBL_MAX) throw "The maximum range --max must be provided.";
		if(iIncrement == DBL_MAX) throw "The range increment --increment must be provided.";
		
		if(iMode == ConstantIncidence && iIncidenceAngle == DBL_MAX) throw "In constant incidence mode, the incidence angle --incidenceAngle must be provided.";
		if(iMode == ConstantIncludedAngle && iIncludedAngle == DBL_MAX) throw "In constant included angle mode, the included angle --includedAngle must be provided.";
		if(iMode == ConstantIncludedAngle && iToOrder == INT_MAX) throw "In constant included angle mode, the operating order --toOrder must be provided.";
		if(iMode == ConstantWavelength && iWavelength == DBL_MAX) throw "In constant wavelength mode, the wavelength --wavelength must be provided.";
		
		if(iOutputFile.empty()) throw "The output data file --outputFile must be provided.";
		
		if(iProfile == PEGrating::InvalidProfile) throw "The grating profile --gratingType must be provided.";
		if(iMaterial.empty()) throw "The grating material --gratingMaterial must be provided.";
		// todo: check that material is valid.
		if(iPeriod == DBL_MAX) throw "The grating period --gratingPeriod must be provided.";
		
		if(iProfile == PEGrating::RectangularProfile && (iGeometry[0] == DBL_MAX || iGeometry[1] == DBL_MAX)) throw "The rectangular profile requires two arguments to --gratingGeometry <depth>,<valleyWidth>.";
		if(iProfile == PEGrating::BlazedProfile && (iGeometry[0] == DBL_MAX || iGeometry[1] == DBL_MAX)) throw "The blazed profile requires two arguments to --gratingGeometry <blazeAngleDeg>,<antiBlazeAngleDeg>.";
		if(iProfile == PEGrating::SinusoidalProfile && (iGeometry[0] == DBL_MAX)) throw "The sinusoidal profile requires one arguments to --gratingGeometry <depth>.";
		if(iProfile == PEGrating::TrapezoidalProfile && (iGeometry[0] == DBL_MAX || iGeometry[1] == DBL_MAX || iGeometry[2] == DBL_MAX || iGeometry[3] == DBL_MAX)) throw "The trapezoidal profile requires four arguments to --gratingGeometry <depth>,<valleyWidth>,<blazeAngle>,<antiBlazeAngle>.";
		
		if(iN == INT_MAX) throw "The truncation index --N must be provided.";		
	}
	
	
	catch(const char* errMsg) {
		std::cerr << "Invalid command-line options: " << errMsg << std::endl;
		return false;
	}
	
	return true;
}

/// This helper function writes the header to the output file stream
void writeOutputFileHeader(std::ostream& of) {
	of << "# Input" << std::endl;
	
	switch(iMode) {
		case ConstantIncidence:
		of << "mode=constantIncidence" << std::endl;
		of << "incidenceAngle=" << iIncidenceAngle << std::endl;
		break;
		case ConstantIncludedAngle:
		of << "mode=constantIncludedAngle" << std::endl;
		of << "includedAngle=" << iIncludedAngle << std::endl;
		of << "toOrder=" << iToOrder << std::endl;
		break;
		case ConstantWavelength:
		of << "mode=constantWavelength" << std::endl;
		of << "wavelength/eV=" << iWavelength << std::endl;
		break;
		default:
		of << "mode=invalid" << std::endl;
		break;
	}
	
	if(iEv)
		of << "units=eV" << std::endl;
	else
		of << "units=um" << std::endl;
		
	of << "min=" << iMin << std::endl;
	of << "max=" << iMax << std::endl;
	of << "increment=" << iIncrement << std::endl;
	
	switch(iProfile) {
		case PEGrating::RectangularProfile:
		of << "gratingType=rectangular" << std::endl;
		break;
		case PEGrating::BlazedProfile:
		of << "gratingType=blazed" << std::endl;
		break;
		case PEGrating::SinusoidalProfile:
		of << "gratingType=sinusoidal" << std::endl;
		break;
		case PEGrating::TrapezoidalProfile:
		of << "gratingType=trapezoidal" << std::endl;
		break;
		default:
		of << "gratingType=invalid" << std::endl;
		break;
	}
	
	of << "gratingPeriod=" << iPeriod << std::endl;
	
	of << "gratingGeometry=";
	for(int i=0;i<8;++i) {
		if(iGeometry[i] == DBL_MAX)
			break;
		if(i!=0)
			of << ",";
		of << iGeometry[i];
	}
	of << std::endl;
	
	of << "gratingMaterial=" << iMaterial << std::endl;
	
	of << "N=" << iN << std::endl;
}


// This helper function appends the progress to the given output stream
void writeOutputFileProgress(std::ostream& of, int completedSteps, int totalSteps, bool anySuccesses, bool anyFailures) {
	of << "# Progress" << std::endl;
	// not done yet:
	if(completedSteps < totalSteps) {
		of << "status=inProgress" << std::endl;
	}
	// or all done:
	else {
		if(anySuccesses && anyFailures)
			of << "status=someFailed" << std::endl;
		else if(anyFailures)
			of << "status=allFailed" << std::endl;
		else
			of << "status=succeeded" << std::endl;
	}
	
	of << "completedSteps=" << completedSteps << std::endl;
	of << "totalSteps=" << totalSteps << std::endl;
}

// This helper function appends a single efficiency result to the output file stream
void writeOutputFileResult(std::ostream& of, const PEResult& result) {
	
	// the "x column" depends on the mode. In constant wavelength mode, it should be the incidence angle. In all others, it should be the wavelength, unless iEv is set, in which case we should convert the wavelength to eV.
	if(iMode == ConstantWavelength) {
		of << result.incidenceDeg << "\t";
	}
	else {
		of << (iEv ? M_HC / result.wavelength : result.wavelength) << "\t";
	}
	
	switch(result.status) {
		case PEResult::InvalidGratingFailure:
		of << "Error:InvalidGratingFailure" << std::endl;
		break;
		case PEResult::ConvergenceFailure:
		of << "Error:ConvergenceFailure" << std::endl;
		break;
		case PEResult::InsufficientCoefficientsFailure:
		of << "Error:InsufficientCoefficientsFailure" << std::endl;
		break;
		case PEResult::AlgebraError:
		of << "Error:AlgebraError" << std::endl;
		break;
		case PEResult::OtherFailure:
		of << "Error:OtherFailure" << std::endl;
		break;
		case PEResult::Success:
		for(int i=0, cc=result.insideEff.size(); i<cc; ++i) {
			if(i!=0)
				of << ",";
			of << result.insideEff.at(cc-1-i);	// insideEff is in the order <e0><e-1><e-2>...  But we actually want the opposite order <e-N><e-N+1>...<e0>.  Kinda silly that we split them up like this in PEResult, just to reassemble results in [-N, N] order.
		}
		for(int i=1, cc=result.outsideEff.size(); i<cc; ++i) {	// omit the 0 order, since we already have it from insideEff.
			of << "," << result.outsideEff.at(i);
		}
		of << std::endl;
		break;
	}
}


