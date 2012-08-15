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

// usage: impFit [numThreads]
int main(int argc, char** argv) {

	int numThreads = 1;
	if(argc == 2)
		numThreads = atoi(argv[1]);

//	for(int i=0; i<11; ++i) {
//		double nm = i*0.5;
//		std::cout << "Roughness for sigma of nm" << nm << std::endl;
//		// test out Sinha reflectivity factor
//		for(int j=0; j<201; ++j) {
//			double ev = 100 + j*5;
//			std::cout << PEGrating::roughnessFactor(nm/1000, M_HC/ev, "Pt", 87) << "\t";
//		}
//		std::cout << "\n" << std::endl;
//	}


	double eV[] = {90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770};

	double meas1[] = {0.0815564, 0.093514, 0.105582, 0.11814, 0.135631, 0.149823, 0.162523, 0.172861, 0.182804, 0.191646, 0.209792, 0.223341, 0.230919, 0.238659, 0.246419, 0.253504, 0.260322, 0.268056, 0.272541, 0.27357, 0.273369, 0.273167, 0.272965, 0.272764, 0.272562, 0.271745, 0.271404, 0.269776, 0.270649, 0.263382, 0.256468, 0.249956, 0.243445, 0.236933, 0.230422, 0.223911, 0.218922, 0.212678, 0.205257, 0.197053, 0.188374, 0.179477, 0.169791, 0.159437, 0.139532, 0.101926, 0.119324, 0.109036, 0.107686, 0.103872, 0.0996485, 0.0962727, 0.0939692, 0.0912157, 0.0894392, 0.0878844, 0.0857625, 0.0834982, 0.0807857, 0.0776163, 0.0741904, 0.070469, 0.0662804, 0.0615201, 0.0562261, 0.0502769, 0.0437055, 0.0364555, 0.0287281};

	double meas2[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0.018, 0.0222968, 0.0262367, 0.029718, 0.0338954, 0.0384051, 0.0431593, 0.0483114, 0.0538261, 0.0590691, 0.0631717, 0.0661967, 0.0692218, 0.0724852, 0.0787447, 0.0855707, 0.0918654, 0.0977083, 0.102966, 0.107498, 0.112175, 0.115077, 0.117292, 0.119508, 0.121723, 0.123939, 0.126154, 0.127448, 0.127343, 0.126495, 0.124794, 0.122432, 0.119618, 0.11586, 0.110805, 0.0967016, 0.068768, 0.0866279, 0.0781171, 0.0766141, 0.0744911, 0.0706644, 0.0670927, 0.0645003, 0.0612718, 0.0590344, 0.0580089, 0.0545655, 0.05217, 0.0493318, 0.0460774, 0.0422342, 0.0385539, 0.0353424, 0.030601, 0.0268897, 0.0233295, 0.0196904, 0.0158567, 0.0111016};

	const int numEvs = sizeof(eV)/sizeof(double);
	std::cout << "Minimization of SSE over " << numEvs << "energy points" << std::endl;
	double calc1[numEvs];
	double calc2[numEvs];
	double calc1R[numEvs];
	double calc2R[numEvs];

	// Parameters to probe:
	double blazes[] = {1.3, 1.4, 1.5, 1.55, 1.6, 1.65, 1.7, 1.8};
	double antiBlazes[] = { 3, 5, 7.5, 10, 12.5, 30 };
	double thicknesses[] = { 1, 1.5, 2, 2.5, 3, 3.5, 4 };	// in nm!
	double sigmas[] = { 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, .4, .5 };

	double period = 1000 / 892.86;	// from PSD measurements
	double incidence = 87;

	double bestSSE = std::numeric_limits<double>::infinity();
	double bestBlaze, bestAntiBlaze, bestThickness, bestSigma, bestScale1, bestScale2;


	for(int i=0; i<sizeof(blazes)/sizeof(double); ++i) {
		double blaze = blazes[i];

		for(int j=0; j<sizeof(antiBlazes)/sizeof(double); ++j) {
			double antiBlaze = antiBlazes[j];

			for(int k=0; k<sizeof(thicknesses)/sizeof(double); ++k) {
				double thickness = thicknesses[k];

				// At this point, can calculate all the efficiencies
				PEBlazedGrating g(period, blaze, antiBlaze, "Ni", "NiO", thickness/1000.0);
				for(int e=0; e<numEvs; ++e) {
					PEResult r = g.getEff(incidence, M_HC/eV[e], 0, PEMathOptions(), false, numThreads);
					if(r.status != PEResult::Success) {
						std::cout << "Calculation error:" << blaze << " " << antiBlaze << " " << thickness << " " << eV[e] << std::endl;
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
						double R = PEGrating::roughnessFactor(sigma/1000.0, M_HC/eV[e], "NiO", incidence);

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
							std::cout << "IR:\t" << blaze << "\t" << antiBlaze << "\t" << thickness << "\t" << sigma << "\t" << scale1 << "\t" << scale2 << "\t" << sse << std::endl;
							if(sse < bestSSE) {
								bestSSE = sse;
								bestBlaze = blaze;
								bestAntiBlaze = antiBlaze;
								bestThickness = thickness;
								bestSigma = sigma;
								bestScale1 = scale1;
								bestScale2 = scale2;
							}
						}
					}
				}
			}
		}
	}

	// final results:
	std::cout << std::endl;
	std::cout << "Best:\t" << bestBlaze << "\t" << bestAntiBlaze << "\t" << bestThickness << "\t" << bestSigma << "\t" << bestSSE << std::endl;

	return 0;
}


