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

	/// Computes alpha_, beta1_, and betaM_ for all n, based on v_1_ (material refractive index) and wl_.  Also calculates
	void computeAlphaAndBeta(double incidenceDeg);

	/// Calculates how many vertical layers (numLayers_ and M_) are sufficient to keep exponentials from contamination.  Fills y_ with the vertical coordinate at each layer.
	void computeLayers();

	/// Computes the blocks of the T matrix (T11_, T12_, T21_, T22_) for the layer below \c y_[m].  Since \c m = 1 is the top of the substrate, \c m can range from [2, M-1].
	PEResult::Code computeTMatrixBelowLayer(int m, bool printDebugOutput = false);

	/// Calculates the grating fourier expansion for k^2 at a given \c y value and wavelength \c wl, and stores in \c k2.  \c k2 must have space for 4*N_ + 1 coefficients, since we will be computing from n = -2N_ to 2N.   Reads member variables N_, wavelength wl_, grating refractive index \c v_1_, and grating geometry from \c g_.  Returns PEResult::Success, or PEResult::InvalidGratingFailure if the profile is not supported or \c y is larger than the groove height.
	PEResult::Code computeGratingExpansion(double y, gsl_complex* k2) const;

	/// Initializes an 8N+4 array of double [\c u, \c uprime] to contain the starting integration values of the electric field Fourier components. The first half of the array \c w contains \c u, the second half contains \c uprime, with each entry in {re,im} order.  The u value is set to $\delta_{n,p}$ and the u' value is set to $-i \beta_n^{(M)} \delta_{n,p}$ or $i \beta_n^{(M)} \delta_{n,p}$, depending on whether p > 2N_.  (Note n,p here are using numbering from 0, and that for the delta functions, p is aliased back onto [0, 2N] once it reaches 2N_+1.
	/*! If layer \c m = 1, then uses beta1_n instead of betaM_n for the derivative. */
	void setIntegrationStartingValues(double* w, int p, int m);

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


	/// Computes the \c BM_ outgoing reflected Rayleigh coefficients, based on a finished S matrix (S12_ block).
	void computeBMFromSMatrix();




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
	/// 2*N_ + 1, 4*N_+2, and 8*N_+4, since these are used a lot
	int twoNp1_, fourNp2_, eightNp4_;
	
	/// The accuracy to use for numerical integration. Default 1e-8. \todo Get from math options.
	double integrationTolerance_;
	
	
	/// alpha array (size 2N+1)
	double* alpha_;
	/// beta array (size 2N+1).  betaM_ is for the superstrate, beta1_ is for the substrate.
	gsl_complex* betaM_, * beta1_;
	/// B_n^{M} array: outgoing reflected Rayleigh coefficients. (size 2N+1)
	gsl_complex* BM_;

	/// Number of layers to use in the vertical stack to keep the growing exponentials from numerical contamination.
	int numLayers_;
	/// M-2 = numLayers_.
	int M_;
	
	/// pre-allocated storage for the grating k^2 fourier coefficients.  There is one array for each thread to use.  Array size must be 2N+1.
	gsl_complex** k2_;
	/// Helper function: returns the k^2 array that should be used by a given thread.
	gsl_complex* k2ForCurrentThread();
	
	/// This block of storage contains the [u, uprime] electric field Fourier component vectors. They are arranged with each value {re,im}, from order [-N to N], the u vector followed by uprime vector... repeated for each trial solution.  Access the [u, uprime] vector for a given trial solution with wVectorForP().   The size of one w vector is (8*N_+4), and there are (4*N_+2) trial solutions.
	double* wVectors_;

	/// Returns the wVector for a trial solution \c p at index \c j, where \c j is numbered from [0, 4*N_+1].
	double* wVectorForP(int j) { return wVectors_ + eightNp4_*j; }
	/// Returns the electric field Fourier component \c u for order \c n (index \c i) and trial solution \c p (index \c j), out of wVectors_.  i and j are numbered from 0.
	gsl_complex* u(int i, int j) {
		return (gsl_complex*)(wVectorForP(j) + 2*i);
	}
	/// Returns the electric field Fourier component derivative \c u' for order \c n (index \c i) and trial solution \c p (index \c j), out of wVectors_.  i and j are numbered from 0.
	gsl_complex* uprime(int i, int j) {
		return (gsl_complex*)(wVectorForP(j) + fourNp2_ + 2*i);
	}

	/// Blocks of T matrix, used in computation of a single layer.
	gsl_matrix_complex* T11_, *T12_, *T21_, *T22_;
	/// Blocks of S matrix, used in recursive computation of everything up to current layer.
	gsl_matrix_complex* S12_, *S22_;
	/// Inverse of Z-matrix, used in computation of S. (Note: we don't actually compute any inverses; Zinv_ is directly calculated from Zinv^{q+1} = T11^{q+1} + T12^{q+1} S12^{1}, and then we use LU decomp and multiplication to avoid loss of precision in computing Z = Zinv_^{-1}.)
	gsl_matrix_complex* Zinv_;
	/// This is a 2*N_+1 x 2*N_+1 matrix used as a workspace matrix.
	gsl_matrix_complex* Z_, * workMatrix_;

	
	/// A permutation structure (size 2*N + 1) used to solve the linear system.
	gsl_permutation* permutation_;
	
		
	/// wavelength for the current calculation
	double wl_;
	/// refractive index of the grating material, at wl_
	gsl_complex v_1_;

	/// The y-coordinate of the infinitely-thin Rayleigh layer at y_m, with m = [1, M_ - 1].  y_[0] is unused, so that we can take y_m = y_[m].
	double* y_;
	
	/// a reference to the grating we're solving
	const PEGrating& g_;
	
	/// A flag that indicates that we should measure the time required for all related blocks of operations
	bool measureTiming_;
	/// Stores the timing results:
	double timing_[12];	
};



#endif
