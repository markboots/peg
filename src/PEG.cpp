#include "PEG.h"
#include "PESolver.h"

#include <gsl/gsl_complex_math.h>
#include <math.h>

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

gsl_complex PEGrating::refractiveIndex(double wl) const {
	return gsl_complex_rect(0.993, -0.00754);
}

double PEGrating::height() const {
	switch(profile_) {
	
	case RectangularProfile:
		return geo(0);
	
	case BlazedProfile:
		// geo(0) is blaze, geo(1) is anti-blaze angle, both in degrees
		return period() / (1/tan(geo(0)*M_PI/180) + 1/tan(geo(1)*M_PI/180));
		
	case SinusoidalProfile:
		// geo(0) is the depth, aka height.
		return geo(0);
		
	case TrapezoidalProfile:
		// geo(0) is the depth, aka height.
		return geo(0);
	
	default: return 0;
	}
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
	case PEResult::AlgebraError:
		os << "Error: Linear Algebra Error" << std::endl;
		break;
	case PEResult::OtherFailure:
		os << "Error: Unknown other Failure" << std::endl;
		break;
	case PEResult::InactiveCalculation:
		os << "Notice: Inactive Calculation" << std::endl;
		break;
	}
	
	return os;
}
