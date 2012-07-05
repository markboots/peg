#ifndef PESOLVER_H
#define PESOLVER_H

#include "PEG.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_linalg.h>

/// Contains the context (memory structures, etc.) and algorithm for solving the grating efficiency.

class PESolver {
public:
	/// Construct a solver context for the given \c grating and math options \c mo.  \c numThreads specifies how many threads to use for fine parallelization; ideally it should be <= the number of processor cores on your computer / on a single cluster node.
	PESolver(const PEGrating& grating, const PEMathOptions& mo = PEMathOptions(), int numThreads = 1, bool measureTiming = false);
	/// Destroy a solver context
	~PESolver();
	
	/// Calculates the efficiency at incidence angle \c incidenceDeg and wavelength \c wl.  Side effects: sets the refractive index member variable v_1_; modifies the contents of u_, uprime_, alpha_, beta_, etc.
	PEResult getEff(double incidenceDeg, double wl, bool printDebugOutput = false);
	
	
	// Helper functions:
	///////////////////////////////
	
	/// Returns the square root \c w of a complex number \c z, choosing the branch cut so that Im(w) >= 0.  
	/*! Note: This is different than the usual branch cut defined for the "principal square root", which uses the negative real axis so that Re(w) >= 0 [ie: w is in the right complex plane].*/
	static gsl_complex complex_sqrt_upperComplexPlane(gsl_complex z);
	
	/// Calculates the grating fourier expansion for k^2 at a given \c y value and wavelength \c wl, and stores in \c k2.  \c k2 must have space for 4*N_ + 1 coefficients, since we will be computing from n = -2N_ to 2N.   Reads member variables N_, wavelength wl_, grating refractive index \c v_1_, and grating geometry from \c g_.  Returns PEResult::Success, or PEResult::InvalidGratingFailure if the profile is not supported or \c y is larger than the groove height.
	PEResult::Code computeGratingExpansion(double y, gsl_complex* k2) const;
	
	/// integrates the complex vectors \c u and \c uprime from y=0 to y=a, using the differential equation and ______ method.  Calls computeGratingExpansion() at each y value, so reads member variables N_, v_1_, and g_.  Modifies k2_ at each step.  Results are returned in-place.
	PEResult::Code integrateTrialSolutionAlongY(gsl_vector_complex* u, gsl_vector_complex* uprime);
	
	/// The function callback for the integration process.  Must be static so we have an address for it, so \c peSolver will be a pointer to a solver (this).
	static int odeFunctionCB(double y, const double w[], double f[], void* peSolver) { 
		PESolver* s = static_cast<PESolver*>(peSolver);
		return s->odeFunction(y, w, f);
	}
	/// called to compute the values for the integration process.
	int odeFunction(double y, const double w[], double f[]);

	/// Static callback function for computing the jacobian for the ODE integration.
	static int odeJacobianCB(double y, const double w[], double * dfdw, double dfdy[], void * params) {
		PESolver* s = static_cast<PESolver*>(params);
		return s->odeJacobian(y, w, dfdw, dfdy);
	}

	/// Called to compute the jacobian for the ODE integration. \c y is the independent variable, \c dfdu is the jacobian matrix in row-major order, and \c dfdy is the ODE function f'(u, y).
	int odeJacobian(double y, const double w[], double * dfdw, double dfdy[]);


	/// Helper function: returns the condition number of a complex matrix
	double conditionNumber(const gsl_matrix_complex* A) const;
	
protected:

	// Number of threads to use for this calculation
	int numThreads_;
	
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
	
	// pre-allocated storage for the grating k^2 fourier coefficients.  There is one array for each thread to use.  Array size must be 4N+1, since we need to compute for n from -2N to 2N.
	gsl_complex** k2_;
	// Helper function: returns the k^2 array that should be used by a given thread.
	gsl_complex* k2ForCurrentThread();
	
	// matrix u. Columns are across trial solutions (p), rows are fourier coefficients n
	gsl_matrix_complex* u_;
	// matrix u' = du/dy.  Columns are across trial solutions (p), rows are fourier coefficients n.
	gsl_matrix_complex* uprime_;
	// This is a diagonal matrix that will have components 1/(i*beta2_n), required in the calculation of the T matrix.
	gsl_matrix_complex* iBeta2Diag_;
	// This is the T matrix that relates the input and ouput basis expansion of the grating.
	gsl_matrix_complex* T_;
	// This vector describes the incident light in the basis expansion
	gsl_vector_complex* Vincident_;
	// This vector contains the rayleigh coefficients A^(1)_n, which also contain the constants of the linear superposition.  T_ A1_ = Vincident_  is an A x = b linear system, that we need to solve for A1_.
	gsl_vector_complex* A1_;
	
	// A permutation structure (size 2*N + 1) used to solve the linear system.
	gsl_permutation* permutation_;
	
	// This vector contains the rayleigh coefficients B^(2)_n, which are essentially the output of the grating solver.
	gsl_vector_complex* B2_;
	
		
	// wavelength for the current calculation
	double wl_;
	// refractive index of the grating material, at wl_
	gsl_complex v_1_;
	
	// a reference to the grating we're solving
	const PEGrating& g_;
	
	// A flag that indicates that we should measure the time required for all related blocks of operations
	bool measureTiming_;
	// Stores the timing results:
	double timing_[12];	
};



#endif
