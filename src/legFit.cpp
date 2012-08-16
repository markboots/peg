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
#include <limits>

/// Computes sum of squared errors; ignores points where meas == 0.
double computeSSE(double meas1[], double calc1[], double scale1, double meas2[], double calc2[], double scale2, int num) {
	double rv = 0;

	for(int i=0; i<num; ++i) {
		if(meas1[i] != 0) {
			double error = meas1[i] - calc1[i]*scale1;
			rv += error*error;
		}
		if(meas2[i] != 0) {
			double error = meas2[i] - calc2[i]*scale2;
			rv += error*error;
		}
	}

	return rv;
}

// usage: legFit [numThreads]
int main(int argc, char** argv) {

	int numThreads = 1;
	if(argc == 2)
		numThreads = atoi(argv[1]);


	double eV[] = {60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275};

	double meas1[] = {0.281228, 0.310422, 0.331426, 0.340257, 0.350128, 0.353944, 0.365529, 0.369424, 0.365019, 0.360616, 0.350993, 0.342171, 0.330305, 0.316607, 0.300431, 0.281799, 0.263826, 0.247158, 0.228298, 0.20805, 0.188193, 0.169017, 0.148416, 0.12863, 0.111252, 0.0971923, 0.0841875, 0.0731424, 0.0636194, 0.055375, 0.0480748, 0.0413767, 0.0352753, 0.0299969, 0.0255217, 0.0216865, 0.018439, 0.0156958, 0.0129751, 0.0107681, 0.00872383, 0.00701543, 0.00563292, 0.00434082};

	double meas2[] = {0.016586, 0.0289956, 0.043031, 0.05686, 0.0735305, 0.0935087, 0.117535, 0.141632, 0.164234, 0.188717, 0.204973, 0.228301, 0.248337, 0.266696, 0.280911, 0.289931, 0.295556, 0.2986, 0.295044, 0.284276, 0.271328, 0.255574, 0.234172, 0.212218, 0.191616, 0.176201, 0.161513, 0.148914, 0.138309, 0.128554, 0.119271, 0.109894, 0.100966, 0.0929721, 0.0861107, 0.0800226, 0.0746219, 0.0690934, 0.0637878, 0.0585326, 0.0536358, 0.0490325, 0.0443844, 0.0396718};

	const int numEvs = sizeof(eV)/sizeof(double);
	std::cout << "Minimization of SSE over " << numEvs << "energy points" << std::endl;
	double calc1[numEvs];
	double calc2[numEvs];
	double calc1R[numEvs];
	double calc2R[numEvs];

	// Parameters to probe:
	double blazes[] = {1.8, 1.9, 2.0, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.5};
	double antiBlazes[] = { 3, 5, 7.5, 10, 12.5, 30, 60 };
//	double thicknesses[] = { 1, 1.5, 2, 2.5, 3, 3.5, 4 };	// in nm!
	double sigmas[] = { 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, .4, .5 };

	double period = 1000 / 593.02;	// from PSD measurements
	double incidence = 86;

	double bestSSE = std::numeric_limits<double>::infinity();
	double bestBlaze, bestAntiBlaze,/* bestThickness,*/ bestSigma, bestScale1, bestScale2;


	for(int i=0; i<sizeof(blazes)/sizeof(double); ++i) {
		double blaze = blazes[i];

		for(int j=0; j<sizeof(antiBlazes)/sizeof(double); ++j) {
			double antiBlaze = antiBlazes[j];


			// At this point, can calculate all the efficiencies
			PEBlazedGrating g(period, blaze, antiBlaze, "Au");
			for(int e=0; e<numEvs; ++e) {
				PEResult r = g.getEff(incidence, M_HC/eV[e], 0, PEMathOptions(), false, numThreads);
				if(r.status != PEResult::Success) {
					std::cout << "Calculation error:" << blaze << " " << antiBlaze << " " << eV[e] << std::endl;
					calc1[e] = 0;
					calc2[e] = 0;
				}
				else {
					calc1[e] = r.eff.at(14);
					calc2[e] = r.eff.at(13);
				}
			}

			// Loop over the reflectivities:
			for(int l=0; l<sizeof(sigmas)/sizeof(double); ++l) {
				double sigma = sigmas[l];

				for(int e=0; e<numEvs; ++e) {
					double R = PEGrating::roughnessFactor(sigma/1000.0, M_HC/eV[e], "Au", incidence);

					calc1R[e] = calc1[e]*R;
					calc2R[e] = calc2[e]*R;
				}


				// loop over scale factors:
				for(int m=20; m<100; ++m) {
					double scale1 = m/100.0;

					for(int n=10; n<100; ++n) {
						double scale2 = n/100.0;

						// Inner step: calculate SSE difference between calc and meas.
						double sse = computeSSE(meas1, calc1R, scale1, meas2, calc2R, scale2, numEvs);
						std::cout << "IR:\t" << blaze << "\t" << antiBlaze << "\t" << sigma << "\t" << scale1 << "\t" << scale2 << "\t" << sse << std::endl;
						if(sse < bestSSE) {
							bestSSE = sse;
							bestBlaze = blaze;
							bestAntiBlaze = antiBlaze;
//							bestThickness = thickness;
							bestSigma = sigma;
							bestScale1 = scale1;
							bestScale2 = scale2;
						}
					}
				}
			}
		}
	}

	// final results:
	std::cout << std::endl;
	std::cout << "Best:\t" << bestBlaze << "\t" << bestAntiBlaze << "\t" << bestSigma << "\t" << bestScale1 << "\t" << bestScale2 << "\t" << bestSSE << std::endl;

	return 0;
}


