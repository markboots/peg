#include "PEG.h"
#include "PESolver.h"

#include <gsl/gsl_complex_math.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string.h>

// This packs the current result into a plain double \c array, for easy communication using standard MPI types.  The array must have pre-allocated room for 2N+1 + 4 elements.
void PEResult::toDoubleArray(double* array) const {
	array[0] = double(status);
	array[1] = wavelength;
	array[2] = incidenceDeg;
	if(eff.empty()) {
		array[3] = 0;
	}
	else {
		array[3] = double((eff.size()-1)/2);	// N
		memcpy(array+4, &(eff.at(0)), eff.size()*sizeof(double));
	}
}

// This unpacks the result from a plain double \c array that was filled by toDoubleArray().
void PEResult::fromDoubleArray(const double* array) {
	status = Code(int(array[0]));
	wavelength = array[1];
	incidenceDeg = array[2];
	
	int N = int(array[3]);
	
	eff.resize(2*N+1);
	memcpy(&(eff[0]), array+4, eff.size()*sizeof(double));	
}

PEResult PEGrating::getEff(double incidenceDeg, double wl, const PEMathOptions& mo, bool printDebugOutput, int numThreads, bool measureTiming) const {
	PESolver s(*this, mo, numThreads, measureTiming);
	return s.getEff(incidenceDeg, wl, printDebugOutput);
}

gsl_complex PEGrating::refractiveIndex(double wl, const std::string& material) {
//	return gsl_complex_rect(0.993, 0.00754);
//	return gsl_complex_rect(1.4, -0.);

	// Attempt to open the material database file.
	std::string fileName = std::string(PEG_MATERIALS_DB_PATH) + std::string("/") + material + std::string(".idx");

	std::ifstream matFile;
	matFile.open(fileName.c_str());
	if(matFile.fail())
		return gsl_complex_rect(0,0);

	std::vector<double> wl_nm, delta, beta;
	wl_nm.reserve(3120);	// my standard files have 3111 lines. This avoids vector re-sizing for performance.
	delta.reserve(3120);
	beta.reserve(3120);

	// read whole file intp vectors for fast searching.
	while(!matFile.eof()) {
		double wl_nmIn, deltaIn, betaIn;
		matFile >> wl_nmIn >> deltaIn >> betaIn;

		wl_nm.push_back(wl_nmIn);
		delta.push_back(deltaIn);
		beta.push_back(betaIn);
	}
	matFile.close();

	if(wl_nm.empty())
		return gsl_complex_rect(0,0);	// the database is empty.

	// Need to convert wavelength from um to nm, since our data files are in the format:
	//   wl(nm)	delta	beta
	// where v = 1 - delta + i*beta       [beta assumed to be >= 0]
	wl *= 1000;	// now wl is in nm.


	// Use binary search to find the entry of the target wavelength in O(log(n)).
	////////////////////////

	// 'larger' is an iterator to one entry higher (or possibly equal) to the target wl.
	std::vector<double>::const_iterator larger = std::lower_bound(wl_nm.begin(), wl_nm.end(), wl);

	// Not found? This wl is larger than any we have in the database.
	if(larger == wl_nm.end())
		return gsl_complex_rect(0,0);

	// Did this return the first entry?
	if(larger == wl_nm.begin()) {
		// if equal to wl, then we've found it at the first entry.
		if(*larger == wl)
			return gsl_complex_rect(1-delta[0], beta[0]);
		else
			return gsl_complex_rect(0,0); // otherwise, it's too low a wavelength for us to have in the database.
	}


	// Not the first entry... So we can also grab the preceeding entry (at lowerIndex), and interpolate between them.
	size_t upperIndex = larger - wl_nm.begin();
	size_t lowerIndex = upperIndex - 1;

	double upperWl = wl_nm[upperIndex];
	double lowerWl = wl_nm[lowerIndex];
	double interp = (wl - lowerWl)/(upperWl - lowerWl);

	double upperDelta = delta[upperIndex];
	double lowerDelta = delta[lowerIndex];
	double interpDelta = lowerDelta + interp*(upperDelta-lowerDelta);

	double upperBeta = beta[upperIndex];
	double lowerBeta = beta[lowerIndex];
	double interpBeta = lowerBeta + interp*(upperBeta - lowerBeta);

	return gsl_complex_rect(1-interpDelta, interpBeta);
}


std::ostream& operator<<(std::ostream& os, const PEResult& result) {
	int N = (result.eff.size()-1)/2;
	
	switch(result.status) {
	case PEResult::Success:
		os << "Inside Orders" << std::endl; 
		for(int i=0,cc=N; i<=cc; ++i) {
			os << i << "\t" << result.eff.at(N-i) << std::endl;
		}
		os << "\nOutside Orders" << std::endl; 
		for(int i=0,cc=N; i<=cc; ++i) {
			os << i << "\t" << result.eff.at(N+i) << std::endl;
		}
		break;
	case PEResult::InvalidGratingFailure:
		os << "Error: Invalid Grating" << std::endl;
		break;
	case PEResult::ConvergenceFailure:
		os << "Error: Integration Convergence Failure" << std::endl;
		break;
	case PEResult::InsufficientCoefficientsFailure:
		os << "Error: Insufficient Coefficients" << std::endl;
		break;
	case PEResult::AlgebraFailure:
		os << "Error: Linear Algebra Error" << std::endl;
		break;
	case PEResult::OtherFailure:
		os << "Error: Unknown other Failure" << std::endl;
		break;
	case PEResult::MissingRefractiveDataFailure:
		os << "Error: Missing refractive index data for this material at this wavelength." << std::endl;
		break;
	case PEResult::InactiveCalculation:
		os << "Notice: Inactive Calculation" << std::endl;
		break;
	}
	
	return os;
}

int PEGrating::computeK2StepsAtY(gsl_complex k2_vaccuum, gsl_complex k2_substrate, gsl_complex k2_coating, double *stepsX, double *stepsK2)
{
	return -1;/// \todo
}
