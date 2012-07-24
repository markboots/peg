/*
Copyright 2012 Mark Boots (mark.boots@usask.ca).

This file is part of the Parallel Efficiency of Gratings project ("PEG").

PEG is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PEG is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PEG.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "PEMainSupport.h"
#include <stdlib.h>
#include <string.h>

void PECommandLineOptions::init() {
	mode = InvalidMode;

	min = max = increment = incidenceAngle = includedAngle = wavelength = DBL_MAX;

	toOrder = N = INT_MAX;

	profile = PEGrating::InvalidProfile;
	period = DBL_MAX;
	geometry[0] = geometry[1] = geometry[2] = geometry[3] = geometry[4] = geometry[5] = geometry[6] = geometry[7] = DBL_MAX;

	eV = false;
	printDebugOutput = false;
	threads = 1;	// by default, just one thread.
	measureTiming = false;
	integrationTolerance = 1e-5;	// default: 1e-5 if not provided.
	coatingThickness = 0;	// default: 0 (no coating) if not provided.
}

// Sets options based on command-line input arguments. Returns isValid().
bool PECommandLineOptions::parseFromCommandLine(int argc, char** argv) {
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
				{"threads", required_argument, 0, 18},
				{"measureTiming", no_argument, 0, 19},
				{"integrationTolerance", required_argument, 0, 20},
				{"coatingMaterial", required_argument, 0, 21},
				{"coatingThickness", required_argument, 0, 22},
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
				if(strcmp(optarg, "constantIncidence") == 0) mode = ConstantIncidence;
				else if(strcmp(optarg, "constantIncludedAngle") == 0) mode = ConstantIncludedAngle;
				else if(strcmp(optarg, "constantWavelength") == 0) mode = ConstantWavelength;
				else throw "The argument to --mode must be one of: constantIncidence, constantIncludedAngle, or constantWavelength.";
				break;
				
			case 2: // min
				//if(!optarg) throw "An argument to --min must be provided.";
				min = atof(optarg);
				break;
				
			case 3: // max
				//if(!optarg) throw "An argument to --max must be provided.";
				max = atof(optarg);
				break;
				
			case 4: // increment
				//if(!optarg) throw "An argument to --increment must be provided.";
				increment = atof(optarg);
				break;
				
			case 5: // incidence angle
				//if(!optarg) throw "An argument to --incidenceAngle must be provided.";
				incidenceAngle = atof(optarg);
				break;
				
			case 6: // includedAngle
				//if(!optarg) throw "An argument to --includedAngle must be provided.";
				includedAngle = atof(optarg);
				break;
				
			case 7: // toOrder
				//if(!optarg) throw "An argument to --toOrder must be provided.";
				toOrder = atol(optarg);
				break;
				
			case 8: // wavelength
				//if(!optarg) throw "An argument to --wavelength must be provided.";
				wavelength = atof(optarg);
				break;
				
			case 9: // outputFile
				//if(!optarg) throw "An argument to --outputFile must be provided.";
				outputFile = optarg;
				break;
				
			case 10: // progressFile
				//if(!optarg) throw "An argument to --progressFile must be provided.";
				progressFile = optarg;
				break;
				
			case 11: // gratingType
				//if(!optarg) throw "An argument to --gratingType must be provided: rectangular, blazed, sinusoidal, or trapezoidal.";
				if(strcmp(optarg, "rectangular") == 0) profile = PEGrating::RectangularProfile;
				else if(strcmp(optarg, "blazed") == 0) profile = PEGrating::BlazedProfile;
				else if(strcmp(optarg, "sinusoidal") == 0) profile = PEGrating::SinusoidalProfile;
				else if(strcmp(optarg, "trapezoidal") == 0) profile = PEGrating::TrapezoidalProfile;
				else throw "The argument to --gratingType must be one of: rectangular, blazed, sinusoidal, or trapezoidal.";
				break;
			case 12: { // gratingGeometry
				//if(!optarg) throw "An argument to --gratingGeometry must be provided, with a list of geometry parameters.";
				char* saveptr;
				char* token;
				int i=0;
				token = strtok_r(optarg, ",", &saveptr);
				while(token && i<8) {
					geometry[i++] = atof(token);
					token = strtok_r(0, ",", &saveptr);
				}
				break;
				}
				
			case 13: // gratingPeriod
				//if(!optarg) throw "An argument to --gratingPeriod must be provided.";
				period = atof(optarg);
				break;
				
			case 14: // gratingMaterial
				material = optarg;
				break;
			
			case 15: // N
				N = atol(optarg);
				break;
				
			case 16: // interpret as eV instead of wavelength
				eV = true;
				break;
			case 17: // print debug output to stdout
				printDebugOutput = true;
				break;
			case 18: // number of threads to use for fine parallelization
				threads = atol(optarg);
				break;
			case 19: // measure timing
				measureTiming = true;
				break;
			case 20: // integrationTolerance
				integrationTolerance = atof(optarg);
				break;
			case 21: // coatingMaterial
				coating = optarg;
				break;
			case 22: // coatingThickness
				coatingThickness = atof(optarg);
				break;
			}
		} // end of loop over input options.
				
	}
	
	catch(const char* errMsg) {
		firstErrorMessage_ = errMsg;
		return false;
	}
	
	return isValid();
}

// Returns true if all required input options have been provided correctly.
bool PECommandLineOptions::isValid() {
	try {
		// Check for missing/invalid input
		if(mode == InvalidMode) throw "The operating mode --mode must be provided: constantIncidence, constantIncludedAngle, or constantWavelength";
		if(min == DBL_MAX) throw "The minimum range --min must be provided.";
		if(max == DBL_MAX) throw "The maximum range --max must be provided.";
		if(increment == DBL_MAX) throw "The range increment --increment must be provided.";
		
		if(mode == ConstantIncidence && incidenceAngle == DBL_MAX) throw "In constant incidence mode, the incidence angle --incidenceAngle must be provided.";
		if(mode == ConstantIncludedAngle && includedAngle == DBL_MAX) throw "In constant included angle mode, the included angle --includedAngle must be provided.";
		if(mode == ConstantIncludedAngle && toOrder == INT_MAX) throw "In constant included angle mode, the operating order --toOrder must be provided.";
		if(mode == ConstantWavelength && wavelength == DBL_MAX) throw "In constant wavelength mode, the wavelength --wavelength must be provided.";
		
		if(outputFile.empty()) throw "The output data file --outputFile must be provided.";
		
		if(profile == PEGrating::InvalidProfile) throw "The grating profile --gratingType must be provided.";
		if(material.empty()) throw "The grating material --gratingMaterial must be provided.";
		// todo: check that material is valid.
		if(period == DBL_MAX) throw "The grating period --gratingPeriod must be provided.";

		if(coatingThickness != 0 && coating.empty()) throw "A coating thickness was specified, but this requires a --coatingMaterial.";
		if(!coating.empty() && coatingThickness == 0) throw "A coating material was specified, but this requires a non-zero --coatingThickness.";

		if(profile == PEGrating::RectangularProfile && (geometry[0] == DBL_MAX || geometry[1] == DBL_MAX)) throw "The rectangular profile requires two arguments to --gratingGeometry <depth>,<valleyWidth>.";
		if(profile == PEGrating::BlazedProfile && (geometry[0] == DBL_MAX || geometry[1] == DBL_MAX)) throw "The blazed profile requires two arguments to --gratingGeometry <blazeAngleDeg>,<antiBlazeAngleDeg>.";
		if(profile == PEGrating::SinusoidalProfile && (geometry[0] == DBL_MAX)) throw "The sinusoidal profile requires one arguments to --gratingGeometry <depth>.";
		if(profile == PEGrating::TrapezoidalProfile && (geometry[0] == DBL_MAX || geometry[1] == DBL_MAX || geometry[2] == DBL_MAX || geometry[3] == DBL_MAX)) throw "The trapezoidal profile requires four arguments to --gratingGeometry <depth>,<valleyWidth>,<blazeAngle>,<antiBlazeAngle>.";
		
		if(N == INT_MAX) throw "The truncation index --N must be provided.";
		if(threads < 1) throw "The number of --threads to use for fine parallelization must be a positive number, at least 1.";
	}
	
	catch(const char* errMsg) {
		firstErrorMessage_ = errMsg;
		return false;
	}
	
	firstErrorMessage_ = "No errors.";
	return true;
}


/// This helper function writes the header to the output file stream
void writeOutputFileHeader(std::ostream& of, const PECommandLineOptions& io) {
	of << "# Input" << std::endl;
	
	switch(io.mode) {
		case PECommandLineOptions::ConstantIncidence:
		of << "mode=constantIncidence" << std::endl;
		of << "incidenceAngle=" << io.incidenceAngle << std::endl;
		break;
		case PECommandLineOptions::ConstantIncludedAngle:
		of << "mode=constantIncludedAngle" << std::endl;
		of << "includedAngle=" << io.includedAngle << std::endl;
		of << "toOrder=" << io.toOrder << std::endl;
		break;
		case PECommandLineOptions::ConstantWavelength:
		of << "mode=constantWavelength" << std::endl;
		of << "wavelength/eV=" << io.wavelength << std::endl;
		break;
		default:
		of << "mode=invalid" << std::endl;
		break;
	}
	
	if(io.eV)
		of << "units=eV" << std::endl;
	else
		of << "units=um" << std::endl;
		
	of << "min=" << io.min << std::endl;
	of << "max=" << io.max << std::endl;
	of << "increment=" << io.increment << std::endl;
	
	switch(io.profile) {
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
	
	of << "gratingPeriod=" << io.period << std::endl;
	
	of << "gratingGeometry=";
	for(int i=0;i<8;++i) {
		if(io.geometry[i] == DBL_MAX)
			break;
		if(i!=0)
			of << ",";
		of << io.geometry[i];
	}
	of << std::endl;
	
	of << "gratingMaterial=" << io.material << std::endl;
	if(io.coatingThickness > 0) {
		of << "coatingMaterial=" << io.coating << std::endl;
		of << "coatingThickness=" << io.coatingThickness << std::endl;
	}
	else {
		of << "coatingMaterial=[none]" << std::endl;
		of << "coatingThickness=0" << std::endl;
	}
	
	of << "N=" << io.N << std::endl;
	of << "integrationTolerance=" << io.integrationTolerance << std::endl;
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
void writeOutputFileResult(std::ostream& of, const PEResult& result, const PECommandLineOptions& io) {
	
	// the "x column" depends on the mode. In constant wavelength mode, it should be the incidence angle. In all others, it should be the wavelength, unless io.eV is set, in which case we should convert the wavelength to eV.
	if(io.mode == PECommandLineOptions::ConstantWavelength) {
		of << result.incidenceDeg << "\t";
	}
	else {
		of << (io.eV ? M_HC / result.wavelength : result.wavelength) << "\t";
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
	case PEResult::AlgebraFailure:
		of << "Error:AlgebraError" << std::endl;
		break;
	case PEResult::OtherFailure:
		of << "Error:OtherFailure" << std::endl;
		break;
	case PEResult::InactiveCalculation:
		of << "Error:InactiveCalculation" << std::endl;
		break;
	case PEResult::MissingRefractiveDataFailure:
		of << "Error:MissingRefractiveIndexData" << std::endl;
		break;
	case PEResult::Success:
		for(int i=0, cc=result.eff.size(); i<cc; ++i) {
			if(i!=0)
				of << ",";
			of << result.eff.at(i);
		}
		of << std::endl;
		break;
	}
}
