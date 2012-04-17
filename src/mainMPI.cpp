#include "PEG.h"
#include "PEMainSupport.h"

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cfloat>
#include <climits>
#include <vector>
#include <algorithm>
#include <cmath>

#include "mpi.h"

/// This main program provides a command-line interface to run a set of parallel grating efficiency calculations. The results are written to an output file, and (optionally) a second file is written to provide information on the status of the calculation.  [This file is only responsible for input processing and output; all numerical details are structured within PEGrating and PESolver.]
/*! 
<b>Command-line options</b>

<b>Required:</b>

Grating specification:

--gratingType <rectangular|blazed|sinusoidal|trapezoidal>
--gratingPeriod <grating period in um>
--gratingGeometry <command-delimited list of geometry parameters, in um and/or degrees>
	Rectangular profile: depth (um),valley width (um)
	Blazed profile: blaze angle (deg),anti-blaze angle (deg)
	Sinusoidal profile: depth (um)
	Trapezoial profile: depth (um),valley width (um),blaze angle (deg),anti-blaze angle (deg)

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
	
	// Initialize MPI
	MPI_Init(&argc, &argv);
	
	// Wait for all processes to be up and running, and then start timing
	MPI_Barrier(MPI_COMM_WORLD);
	double startTime = MPI_Wtime();
	
	int rank, commSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	
	// parse command line options into input variables.
	PECommandLineOptions io;
	// Prevent getopt() from printing error messages except on Process 0.
	if(rank != 0) opterr = 0;
	if(!io.parseFromCommandLine(argc, argv)) {
		if(rank == 0) std::cerr << "Invalid command-line options: " << io.firstErrorMessage() << std::endl;
		MPI_Finalize();
		return -1;
	}
	
	// On Process 0: Open the output file:
	std::ofstream outputFile;
	std::streampos outputFilePosition;
	if(rank == 0) {
		outputFile.open(io.outputFile.c_str(), std::ios::out | std::ios::trunc);
		
		// if(!outputFile.is_open()) {
		// 	std::cerr << "Could not open output file " << io.outputFile;
		// 	return -1;
		// }

		// Check that we can open the progress file, if provided:
		// if(!io.progressFile.empty()) {
		// 	std::ofstream progressFile(io.progressFile.c_str(), std::ios::out | std::ios::trunc);
		// 	if(!progressFile.is_open()) {
		// 		std::cerr << "Could not open progress file " << io.progressFile;
		// 		return -1;
		// 	}
		// }
		
		// Write the file header:
		writeOutputFileHeader(outputFile, io);
		// Remember this position in the output file; it is where we will write the progress and output lines
		outputFilePosition = outputFile.tellp();
	}
	
	// How many steps do we have?
	int totalSteps = int((io.max - io.min)/io.increment) + 1;
	
	// On Process 0: Write the initial progress:
	if(rank == 0) {
		writeOutputFileProgress(outputFile, 0, totalSteps, false, false);
		if(!io.progressFile.empty()) {
			std::ofstream progressFile(io.progressFile.c_str(), std::ios::out | std::ios::trunc);
			writeOutputFileProgress(progressFile, 0, totalSteps, false, false);
		}
	}
	
	// create the grating object.
	PEGrating* grating;
	switch(io.profile) {
	case PEGrating::RectangularProfile:
		grating = new PERectangularGrating(io.period, io.geometry[0], io.geometry[1], io.material);
		break;
	case PEGrating::BlazedProfile:
		grating = new PEBlazedGrating(io.period, io.geometry[0], io.geometry[1], io.material);
		break;
	case PEGrating::SinusoidalProfile:
		grating = new PESinusoidalGrating(io.period, io.geometry[0], io.material);
		break;
	case PEGrating::TrapezoidalProfile:
		grating = new PETrapezoidalGrating(io.period, io.geometry[0], io.geometry[1], io.geometry[2], io.geometry[3], io.material);
		break;
	default:
		grating = 0;	// this will never happen; input validation assures one of the valid grating types.
		break;
	}
	
	// set math options: truncation index from input.
	PEMathOptions mathOptions(io.N);
	
	// On Process 0: output data will be stored here:
	bool anyFailures = false;
	bool anySuccesses = false;
	std::vector<PEResult> results;
	// On process 0: create a buffer for receiving results from other processes
	int resultSize = 2*io.N+1 + 4;
	double* mpiReceiveBuffer = 0;
	if(rank == 0)
		mpiReceiveBuffer = new double[resultSize * commSize];
	// On all processes, create a send buffer for packing up the results
	double* mpiSendBuffer = new double[resultSize];
	
	// Loop over calculation steps.  Loop goes up by commSize each round, since we handle that many steps simultaneously.
	for(int i=0; i<totalSteps; i+=commSize) {
		
		double currentValue = io.min + io.increment*(i+rank);	// creates a cyclic partition. i=0 to P0, i=1 to P1, i=2 to P2... 
		
		// determine wavelength (um): depends on mode and eV/um setting.
		double wavelength = (io.mode == PECommandLineOptions::ConstantWavelength) ? io.wavelength : currentValue;
		if(io.eV)
			wavelength = M_HC / wavelength;	// interpret input wavelength as eV instead, and convert to actual wavelength.  Formula: wavelength = hc / eV.     hc = 1.23984172 eV * um.
		
		// determine incidence angle: depends on mode and possibly wavelength.
		double incidenceAngle;
		switch(io.mode) {
		case PECommandLineOptions::ConstantIncidence:
			incidenceAngle = io.incidenceAngle;
			break;
		case PECommandLineOptions::ConstantIncludedAngle: {
			double ciaRad = io.includedAngle * M_PI / 180;
			incidenceAngle = (asin(-io.toOrder*wavelength/2/io.period/cos(ciaRad/2)) + ciaRad/2) * 180 / M_PI;	// formula for constant included angle: satisfies alpha + beta = cia, and grating equation io.toOrder*wavelength/d = sin(beta) - sin(alpha).
			break;
			}
		case PECommandLineOptions::ConstantWavelength:
			incidenceAngle = currentValue;
			break;
		default:
			incidenceAngle = 0; // never happens: input validation assures valid mode.
			break;
		}
		
		// run the calculation, but only if (i+rank) is still in range.  The last processes will have nothing to do on the last round, if the number of steps does not divided evenly by the number of processes.
		PEResult result = PEResult(PEResult::InactiveCalculation);
		if(i+rank < totalSteps)
			result = grating->getEff(incidenceAngle, wavelength, mathOptions, (io.printDebugOutput && rank == 0), io.threads);	/// Debug output only shown on Process 0?
		
		// Pack up result into the send buffer
		result.toDoubleArray(mpiSendBuffer);
		
		// MPI Gather all results for this round onto Process 0.
		int err = MPI_Gather(mpiSendBuffer, resultSize, MPI_DOUBLE, mpiReceiveBuffer, resultSize, MPI_DOUBLE, 0, MPI_COMM_WORLD); /// \todo Err check

		// On Process 0: collect results and output.
		if(rank == 0) {
			
			for(int j=0; j<commSize; ++j) {
				PEResult result;
				result.fromDoubleArray(mpiReceiveBuffer + j*resultSize);
				
				if(result.status != PEResult::InactiveCalculation) {	// indicates non-calculation for inactive process on last round
					results.push_back(result);
					if(result.status == PEResult::Success)
						anySuccesses = true;
					else
						anyFailures = true;
				}
			}
		
			// Print progress and results to output file.
			outputFile.seekp(outputFilePosition);
			writeOutputFileProgress(outputFile, std::min(i+commSize, totalSteps), totalSteps, anySuccesses, anyFailures);
			outputFile << "# Output" << std::endl;
			for(int j=0, cc=results.size(); j<cc; ++j)	// kinda lame and expensive that we need to do this on each step. Maybe switch to append-only mode, and leave the progress in just progressFile?
				writeOutputFileResult(outputFile, results.at(j), io);
			
			// Update progress in progressFile, if provided.
			if(!io.progressFile.empty()) {
				std::ofstream progressFile(io.progressFile.c_str(), std::ios::out | std::ios::trunc);
				writeOutputFileProgress(progressFile, std::min(i+commSize, totalSteps), totalSteps, anySuccesses, anyFailures);
			}
		}

	} // end of calculation loop.
	
	// Timing: We know we're synchronized here because the last MPI_Gather has ensured that we have everyone's results.
	double runTime = MPI_Wtime() - startTime;
	if(rank == 0)
		std::cout << "Run time (s): " << runTime;

	outputFile.close();
	delete grating;
	delete [] mpiReceiveBuffer;
	delete [] mpiSendBuffer;
	
	// Finalize MPI
	MPI_Finalize();
	return 0;
}
