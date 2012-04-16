#ifndef PEMAINSUPPORT_H
#define PEMAINSUPPORT_H

#include "PEG.h"
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cfloat>
#include <climits>
#include <vector>
#include <cmath>

/// This file contains common support routines for both the parallel and sequential versions of the PEG command-line app, related to input and output handling.

/// This structure specifies the operating input, and is able to parse it from command-line options \c argc, \c argv.
class PECommandLineOptions {
public:
	enum Mode {InvalidMode, ConstantIncidence, ConstantIncludedAngle, ConstantWavelength};
	
	// Input variables:
	////////////////////////////////
	Mode mode;

	double min, max,increment,incidenceAngle,includedAngle,wavelength;
	int toOrder;

	int N;

	PEGrating::Profile profile;
	double period;
	double geometry[8];
	std::string material;

	std::string outputFile, progressFile;

	bool eV;
	bool printDebugOutput;
	////////////////////////////////
	
	/// Default constructor initializes all input variables to recognizable values. Doubles are set to DBL_MAX, and integers are set to INT_MAX.
	PECommandLineOptions() {
		init();
	}
	
	/// This constructor parses immediately from the command-line input arguments. Check for validity after with isValid().
	PECommandLineOptions(int argc, char** argv) {
		init();
		parseFromCommandLine(argc, argv);
	}
	
	/// Sets options based on command-line input arguments. Returns isValid().
	bool parseFromCommandLine(int argc, char** argv);
	
	/// Returns true if all required input options have been provided correctly. Sets firstErrorMessage() if any problems found.
	bool isValid();
	
	/// If !isValid(), returns a description of the first error found during validation.
	std::string firstErrorMessage() const { return firstErrorMessage_; }
	
protected:
	/// Initializes all input variables to recognizable values. Doubles are set to DBL_MAX, and integers are set to INT_MAX.
	void init();
	/// A description of the first error found during validation.
	std::string firstErrorMessage_;
};


/// This helper function writes the header to the output file stream
void writeOutputFileHeader(std::ostream& outputFileStream, const PECommandLineOptions& io);

/// This helper function appends a single efficiency result to the output file stream
void writeOutputFileResult(std::ostream& outputFileStream, const PEResult& result, const PECommandLineOptions& io);

/// This helper function appends the progress description to the given output stream
void writeOutputFileProgress(std::ostream& outputFileStream, int completedSteps, int totalSteps, bool anySuccesses, bool anyFailures);

#endif
