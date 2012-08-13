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

// usage: megFit [numThreads]
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


	double eV[] = {182.314, 191.789, 201.264, 210.739, 220.214, 229.689, 239.164, 248.639, 258.114, 267.589, 277.064, 343.389, 352.864, 362.339, 371.814, 381.289, 390.764, 400.239, 409.714, 419.189, 428.664, 438.139, 447.614, 457.089, 466.564, 476.039, 485.514, 494.989, 504.464, 513.939, 523.414, 532.889, 542.364, 551.839, 561.314, 570.789, 580.264, 589.739, 599.214, 608.689, 618.164, 627.639, 637.114, 646.589, 656.064, 665.539, 675.014, 684.489, 693.964, 703.439, 712.914, 722.389, 731.864, 741.339, 750.814, 760.289, 769.764, 779.239, 788.714, 798.189, 807.664, 817.139, 826.614, 836.089, 845.564, 855.039, 864.514, 873.989, 883.464, 892.939, 902.414, 911.889, 921.364, 930.839};

	double meas1[] = {0.12429, 0.131575, 0.135866, 0.137821, 0.139759, 0.142126, 0.144167, 0.146117, 0.149038, 0.150986, 0.151906, 0.149091, 0.150101, 0.150205, 0.150655, 0.150706, 0.149655, 0.145499, 0.144, 0.145194, 0.141457, 0.141205, 0.141077, 0.139927, 0.138129, 0.136396, 0.133419, 0.130511, 0.127373, 0.123481, 0.118338, 0.0977342, 0.0994115, 0.0979107, 0.093027, 0.0944021, 0.0944473, 0.094026, 0.0941694, 0.0953542, 0.0959479, 0.0974813, 0.098334, 0.0993646, 0.100155, 0.100704, 0.100654, 0.100787, 0.100705, 0.0994293, 0.0999224, 0.0995201, 0.0981428, 0.0966052, 0.0946371, 0.092154, 0.088207, 0.0850923, 0.0807993, 0.0757604, 0.0687788, 0.059506, 0.0451735, 0.0223253, 0.00422916, 0.0250565, 0.0124156, 0.0190871, 0.016859, 0.0179062, 0.0181021, 0.0193052, 0.0205697, 0.023324};

	double meas2[] = {0.0053669, 0.00788878, 0.00938365, 0.0113988, 0.013051, 0.0154457, 0.017439, 0.0195011, 0.0217893, 0.0244275, 0.0263301, 0.0381059, 0.0395776, 0.0412652, 0.0429173, 0.0437528, 0.0445883, 0.0454237, 0.0462592, 0.0470947, 0.0479302, 0.0487656, 0.0492273, 0.0491951, 0.0488673, 0.0481245, 0.0471649, 0.0458994, 0.0443843, 0.0425669, 0.040052, 0.0311961, 0.0296452, 0.0310038, 0.0282867, 0.0282438, 0.0270318, 0.0255272, 0.0247226, 0.0239774, 0.0235368, 0.0233225, 0.023092, 0.0231301, 0.0230617, 0.0230717, 0.0231335, 0.0230556, 0.0232607, 0.0226754, 0.0235178, 0.0234624, 0.0231616, 0.0228089, 0.0224573, 0.0220269, 0.0209773, 0.0202602, 0.0192512, 0.018334, 0.0171155, 0.0151113, 0.012068, 0.00717441, 0.00273204, 0.0133314, 0.00707034, 0.0101854, 0.00904179, 0.00961182, 0.00994689, 0.010846, 0.0111166, 0.011778};

	const int numEvs = sizeof(eV)/sizeof(double);
	std::cout << "Minimization of SSE over " << numEvs << "energy points" << std::endl;
	double calc1[numEvs];
	double calc2[numEvs];
	double calc1R[numEvs];
	double calc2R[numEvs];

	// Parameters to probe:
	double blazes[] = {1.7, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.3};
	double antiBlazes[] = { 3, 5, 7.5, 10, 12.5, 30 };
	double thicknesses[] = { 1, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4 };	// in nm!
	double sigmas[] = { 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, .4, .5 };

	double period = 1000 / 1187.819;	// from PSD measurements
	double incidence = 88;

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
					PEResult r = g.getEff(incidence, M_HC/eV[e], PEMathOptions(), false, numThreads);
					if(r.status != PEResult::Success) {
						std::cout << "Calculation error:" << blaze << " " << antiBlaze << " " << thickness << " " << eV[e] << " Code: " << r.status << std::endl;
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
	std::cout << "Best:\t" << bestBlaze << "\t" << bestAntiBlaze << "\t" << bestThickness << "\t" << bestSigma << "\t" << bestScale1 << "\t" << bestScale2 << "\t" << bestSSE << std::endl;

	return 0;
}


