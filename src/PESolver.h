#ifndef PESOLVER_H
#define PESOLVER_H

#include "PEG.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>

/// Contains the context (memory structures, etc.) and algorithm for solving the grating efficiency.

class PESolver {
public:
	/// Construct a solver context for the given \c grating and math options \c mo
	PESolver(const PEGrating& grating, const PEMathOptions& mo = PEMathOptions());
	/// Destroy a solver context
	~PESolver();
	
	/// Calculates the efficiency at incidence angle \c incidenceDeg and wavelength \c wl.  Side effects: sets the refractive index member variable v_1_; modifies the contents of u_, uprime_, alpha_, beta_, etc.
	PEResult getEff(double incidenceDeg, double wl);
	
	
	// Helper functions:
	///////////////////////////////
	
	/// Returns the square root \c w of a complex number \c z, choosing the branch cut so that Im(w) >= 0.  
	/*! Note: This is different than the usual branch cut defined for the "principal square root", which uses the negative real axis so that Re(w) >= 0 [ie: w is in the right complex plane].*/
	static gsl_complex complex_sqrt_upperComplexPlane(gsl_complex z);
	
	/// Calculates the grating fourier expansion for k^2 at a given \c y value and wavelength \c wl, and stores in \c k2.  \c k2 must have space for 4*N_ + 1 coefficients, since we will be computing from n = -2N_ to 2N.   Reads member variables N_, grating refractive index \c v_1_, and grating geometry from \c g_.  Returns PEResult::Success, or PEResult::InvalidGratingFailure if the profile is not supported or \c y is larger than the groove height.
	PEResult::Code computeGratingExpansion(double y, double wl, gsl_complex* k2) const;
	
protected:

	// The number of Fourier coefficients
	int N_;
	// 2*N_ + 1, since this is used a lot
	int twoNp1_;
	
	// The number of integration steps to use
	int niy_;
	
	
	// alpha array (size 2N+1)
	double* alpha_;
	// beta array (size 2N+1)
	gsl_complex* beta2_, * beta1_;
	
	// pre-allocated storage for the grating k^2 fourier coefficients.  Size must be 4N+1, since we need to compute for n from -2N to 2N.  (This only works for the single-threaded version, since this will be written for every trial solution at every y value.)
	gsl_complex* k2_;
	
	// matrix u. Columns are across trial solutions (p), rows are fourier coefficients n
	gsl_matrix_complex* u_;
	// matrix u' = du/dy.  Columns are across trial solutions (p), rows are fourier coefficients n.
	gsl_matrix_complex* uprime_;
	
	// matrix T, used to solve the linear system 
	
	// refractive index of the grating material
	gsl_complex v_1_;
	
	// a reference to the grating we're solving
	const PEGrating& g_;
	
};



#endif
