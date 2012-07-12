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
	
	







	// Solving implementation functions
	////////////////////////////////////////

	/// Calculates the grating fourier expansion for k^2 at a given \c y value and wavelength \c wl, and stores in \c k2.  \c k2 must have space for 4*N_ + 1 coefficients, since we will be computing from n = -2N_ to 2N.   Reads member variables N_, wavelength wl_, grating refractive index \c v_1_, and grating geometry from \c g_.  Returns PEResult::Success, or PEResult::InvalidGratingFailure if the profile is not supported or \c y is larger than the groove height.
	PEResult::Code computeGratingExpansion(double y, gsl_complex* k2) const;

	/// Initializes an 8N+4 array of double [\c u, \c uprime] to contain the starting integration values of the electric field Fourier components. The first half of the array \c w contains \c u, the second half contains \c uprime, with each entry in {re,im} order.  The u value is set to $\delta_{n,p}$ and the u' value is set to $-i \beta_n^{(M)} \delta_{n,p}$ or $i \beta_n^{(M)} \delta_{n,p}$, depending on whether p > 2N_.  (Note n,p here are using numbering from 0, and that for the delta functions, p is aliased back onto [0, 2N] once it reaches 2N_+1.
	void setIntegrationStartingValues(double* w, int p);

	/// Integrates the electric field Fourier component vectors contained in \c w from y = \c yStart to y = \c yEnd, using the differential equation and ______ method.  Array \c w should contain vector \c u followed by \c uprime, with each entry in {re,im} order. Calls computeGratingExpansion() at each y value, so reads member variables N_, v_1_, and g_.  Modifies k2 (for thread) at each step.  Results are returned in-place.
	PEResult::Code integrateTrialSolutionAlongY(double* w, double yStart, double yEnd);
	/// This is an overloaded function. Integrates the electric field Fourier component vectors \c u and \c uprime from y=0 to y=a, using the differential equation and ______ method.  Calls computeGratingExpansion() at each y value, so reads member variables N_, v_1_, and g_.  Modifies k2 (for thread) at each step.  Results are returned in-place.
	PEResult::Code integrateTrialSolutionAlongY(gsl_vector_complex* u, gsl_vector_complex* uprime, double yStart, double yEnd);

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






	// General Mathematical Helper functions:
	///////////////////////////////
	
	/// Returns the square root \c w of a complex number \c z, choosing the branch cut so that Im(w) >= 0.  
	/*! Note: This is different than the usual branch cut defined for the "principal square root", which uses the negative real axis so that Re(w) >= 0 [ie: w is in the right complex plane].*/
	static gsl_complex complex_sqrt_upperComplexPlane(gsl_complex z);

	/// Helper function: returns the condition number of a complex square matrix. INCOMPLETE!
	static double conditionNumber(const gsl_matrix_complex* A);

	/// Solves the linear system AX = B where \c A, \c X, \c B are all \c n x \c n square matrices (\c X is unknown) using the pre-computed LU decomposition of \c A.  Uses gsl_linalg_LU_solve() to do \c n back-substitutions for each column of \c X.
	static int linalg_LU_complex_solve(const gsl_matrix_complex* LU, const gsl_permutation* P, const gsl_matrix_complex* B, gsl_matrix_complex* X);




protected:

	/// Number of threads to use for this calculation
	int numThreads_;
	
	/// The number of Fourier coefficients
	int N_;
	/// 2*N_ + 1, since this is used a lot
	int twoNp1_;
	
	/// The accuracy to use for numerical integration. Default 1e-8. \todo Get from math options.
	double integrationTolerance_;
	
	
	// alpha array (size 2N+1)
	double* alpha_;
	// beta array (size 2N+1)
	gsl_complex* beta2_, * beta1_;
	
	// pre-allocated storage for the grating k^2 fourier coefficients.  There is one array for each thread to use.  Array size must be 4N+1, since we need to compute for n from -2N to 2N.
	gsl_complex** k2_;
	// Helper function: returns the k^2 array that should be used by a given thread.
	gsl_complex* k2ForCurrentThread();
	
	// This block of storage contains the [u, uprime] electric field Fourier component vectors. They are arranged with each value in {re,im} sequential order, the u vector followed by uprime vector... repeated for each trial solution.  Access the [u, uprime] vector for a given trial solution with wVectorForP().
	double* wVectors_;

	/// Returns the wVector for a trial solution \c p, where \c p is numbered from [0, 2*N_].
	double* wVectorForP(int p) { return wVectors_ + (8*N_+4)*p; }
	/// Returns the complex field component \c u for order \c n and trial solution \c p (both numbered from 0 here), out of wVectors_.
	gsl_complex* u(int n, int p) {
		return (gsl_complex*)(wVectorForP(p) + 2*n);
	}
	/// Returns the complex field component derivative \c u' for order \c n and trial solution \c p (both numbered from 0 here), out of wVectors_.
	gsl_complex* uprime(int n, int p) {
		return (gsl_complex*)(wVectorForP(p) + 4*N_+2 + 2*n);
	}

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
