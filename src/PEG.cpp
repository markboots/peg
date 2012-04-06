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
