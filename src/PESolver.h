#ifndef PESOLVER_H
#define PESOLVER_H

#include "PEG.h"

/// Contains the context (memory structures, etc.) and algorithm for solving the grating efficiency.

class PESolver {
public:
	/// Construct a solver context for the given \c grating and math options \c mo
	PESolver(const PEGrating& grating, const PEMathOptions& mo = PEMathOptions());
	/// Destroy a solver context
	~PESolver();
	
	/// Calculates the efficiency at incidence angle \c incidenceDeg and wavelength \c wl.
	PEResult getEff(double incidenceDeg, double wl) const;
	
protected:
			
};



#endif
