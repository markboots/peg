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


/// arguments: wavelength, period (um, um)
int main(int argc, char** argv) {
	
	// Initialize MPI
	MPI_Init(&argc, &argv);
	
	// Wait for all processes to be up and running, and then start timing
	MPI_Barrier(MPI_COMM_WORLD);
	double startTime = MPI_Wtime();
	
	int rank, commSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	
	double wl = atof(argv[1]);
	double period = atof(argv[2]);
	double order = -1;
	double N = 15;

	double startingHeight = 0.004;
	double deltaHeight = 0.002;
	double endingHeight = 0.14;

	double maxEff = 0;
	double maxIncidence = -1, maxHeight = -1;

	// manual loop over heights:
	for(int i=0,cc=(endingHeight-startingHeight)/deltaHeight+1; i<cc; ++i) {
//		MPI_Barrier(MPI_COMM_WORLD);

		double height = startingHeight + i*deltaHeight;

		// determine incidence angle based on MPI rank. Assumes 139 processors.
		double incidence = 75 + rank*0.1;

		PERectangularGrating g(period, height, period*0.5, "Pt");
		PEResult r = g.getEff(incidence, wl, 0, PEMathOptions(N));
		// extract 1st-order efficiency
		double eff = -1;
		if(r.status == PEResult::Success)
			eff = r.eff.at(N + order);
		else
			std::cout << "CalculationError: WL: " << wl << " Period: " << period << " Height: " << height << " Incidence: " << incidence << std::endl;

		// on root process, create buffer for eff from everyone
		std::vector<double> allResults;
		if(rank == 0) {
			allResults.resize(commSize);
		}

		int err = MPI_Gather(&eff, 1, MPI_DOUBLE, allResults.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if(rank == 0) {
			// look through all incidences
			for(int j=0; j<commSize; ++j) {
				if(allResults.at(j) > maxEff) {
					maxEff = allResults.at(j);
					maxIncidence = 75 + j*0.1;
					maxHeight = height;
				}
			}
		}
	}


	// that's it.
	if(rank == 0) {
		std::cout << "WL: " << wl << " Period: " << period << " OptimalIncidence: " << maxIncidence << " Height: " << maxHeight << " Eff: " << maxEff << std::endl;
		std::cout << "Search time: " << MPI_Wtime() - startTime << std::endl;
	}
	
	// Finalize MPI
	MPI_Finalize();
	return 0;
}
