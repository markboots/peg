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


/// arguments: wavelength, period (um, um)
int main(int argc, char** argv) {

	if(argc != 3) {
		std::cout << "Usage: blazedIncidenceSearch [wavelength] [period]" << std::endl;
		return 0;
	}

	double wl = atof(argv[1]);
	double period = atof(argv[2]);
	double order = -1;
	double N = 15;

	double startingHeight = 0.8;	// height is blaze here
	double deltaHeight = 0.2;
	double endingHeight = 3;

	double maxEff = 0;
	double maxIncidence = -1, maxHeight = -1;

	// manual loop over heights:
	for(int i=0,cc=(endingHeight-startingHeight)/deltaHeight+1; i<cc; ++i) {

		double height = startingHeight + i*deltaHeight;

		for(int j=0; j<23; ++j) {
			double incidence = 85 + j*0.2;

			PEBlazedGrating g(period, height, 30, "Pt");
			PEResult r = g.getEff(incidence, wl, PEMathOptions(N), false, 4);
			// extract 1st-order efficiency
			double eff = -1;
			if(r.status == PEResult::Success) {
				eff = r.eff.at(N + order);
				std::cout << " Blaze: " << height << " Incidence: " << incidence << " Eff: " << eff << std::endl;
			}
			else
				std::cout << "Calculation error: WL: " << wl << " Blaze: " << period << " Height: " << height << " Incidence: " << incidence << std::endl;

			if(eff > maxEff) {
				maxEff = eff;
				maxIncidence = incidence;
				maxHeight = height;
			}
		}
	}

	std::cout << "DONE: WL: " << wl << " Period: " << period << " OptimalIncidence: " << maxIncidence << " Height: " << maxHeight << " Eff: " << maxEff << std::endl;

	return 0;
}

