#include "PEG.h"
#include "PESolver.h"

#include <gsl/gsl_complex_math.h>
#include <math.h>

PEResult PEGrating::getEff(double incidenceDeg, double wl, const PEMathOptions& mo) const {
	PESolver s(*this, mo);
	return s.getEff(incidenceDeg, wl);
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
	
	default: return 0;
	}
}


std::ostream& operator<<(std::ostream& os, const PEResult& result) {
	switch(result.status) {
	case PEResult::Success:
		os << "Inside Orders" << std::endl; 
		for(int i=0,cc=result.insideEff.size(); i<cc; ++i) {
			os << i << "\t" << result.insideEff.at(i) << std::endl;
		}
		os << "\nOutside Orders" << std::endl; 
		for(int i=0,cc=result.outsideEff.size(); i<cc; ++i) {
			os << i << "\t" << result.outsideEff.at(i) << std::endl;
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
	}
	
	return os;
}
