#include "PEG.h"
#include "PESolver.h"

#include <gsl/gsl_complex_math.h>

PEResult PEGrating::getEff(double incidenceDeg, double wl, const PEMathOptions& mo) const {
	PESolver s(*this, mo);
	return s.getEff(incidenceDeg, wl);
}

gsl_complex PEGrating::refractiveIndex(double wl) const {
	return gsl_complex_rect(0.993, -0.00754);
}

