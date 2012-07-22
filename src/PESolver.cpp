#include "PESolver.h"

#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <omp.h>
#include <string.h>

// For debug output only:
#include <iostream>

/// Speed of light in um/s.
#define M_c 2.99792458e14

PESolver::PESolver(const PEGrating& grating, const PEMathOptions& mo, int numThreads, bool measureTiming)
	: g_(grating)
{
	numThreads_ = numThreads;
	measureTiming_ = measureTiming;
	time_ = omp_get_wtime();
	
	// set math options.
	N_ = mo.N;
	integrationTolerance_ = mo.integrationTolerance;
	twoNp1_ = 2*N_ + 1;
	fourNp2_ = 4*N_ + 2;
	eightNp4_ = 8*N_ + 4;
	
	// allocate matrices and vectors
	wVectors_ = new double[eightNp4_ * fourNp2_];	// size: (2N_+1) orders * 2(for re,im) * 2(for u,uprime) * 2(2N_+1) trial solutions.
	T11_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	T12_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	T21_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	T22_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	S12_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	S22_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	Zinv_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	Z_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	workMatrix_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	permutation_ = gsl_permutation_alloc(twoNp1_);
	
	alpha_ = new double[twoNp1_];
	betaM_ = new gsl_complex[twoNp1_];
	beta1_ = new gsl_complex[twoNp1_];
	BM_ = new gsl_complex[twoNp1_];
	
	// we need one k2_ array for each thread, since they will be used simultaneously.
	k2_ = new gsl_complex*[numThreads_];
	for(int i=0; i<numThreads_; ++i)
		k2_[i] = new gsl_complex[twoNp1_];

	y_ = 0;
	timing_[0] = omp_get_wtime() - time_;		// time to allocate memory.
}


PESolver::~PESolver() {

	gsl_matrix_complex_free(T11_);
	gsl_matrix_complex_free(T12_);
	gsl_matrix_complex_free(T21_);
	gsl_matrix_complex_free(T22_);
	gsl_matrix_complex_free(S12_);
	gsl_matrix_complex_free(S22_);
	gsl_matrix_complex_free(Zinv_);
	gsl_matrix_complex_free(Z_);
	gsl_matrix_complex_free(workMatrix_);
	gsl_permutation_free(permutation_);
	
	delete [] alpha_;
	delete [] betaM_;
	delete [] beta1_;
	delete [] wVectors_;
	delete [] y_;
	delete [] BM_;

	for(int i=0; i<numThreads_; ++i)
		delete [] k2_[i];
	delete [] k2_;
}

/// \todo Imp.
PEResult PESolver::getEff(double incidenceDeg, double wl, bool printDebugOutput) {

	time_ = omp_get_wtime();
	
	// 1. Setup incidence variables and constants
	/////////////////////////////////////////
	
	// Set this as the current wavelength
	wl_ = wl;
	// get material refractive index at wavelength
	v_1_ = g_.substrateRefractiveIndex(wl_);
	if(GSL_REAL(v_1_) == 0.0 && GSL_IMAG(v_1_) == 0.0) {
		return PEResult::MissingRefractiveDataFailure;
	}
	if(g_.coatingThickness() > 0) {
		v_c_ = g_.coatingRefractiveIndex(wl_);
		if(GSL_REAL(v_c_) == 0.0 && GSL_IMAG(v_c_) == 0.0) {
			return PEResult::MissingRefractiveDataFailure;
		}
	}

	timing_[1] = time_;
	time_ = omp_get_wtime();
	timing_[1] = time_ - timing_[1];	// time to look up refractive index.
	
	// 2. compute all alpha_n and beta1_n, betaM_n.
	///////////////////////////////////
	
	computeAlphaAndBeta(incidenceDeg);	// (This uses a parallel loop).

	// Calculates how many vertical layers we need, and the division into slices at y_.
	computeLayers();

	timing_[2] = time_;
	time_ = omp_get_wtime();
	timing_[2] = time_ - timing_[2];	// time to calculate alpha, beta, and layers.

	if(printDebugOutput) {
		std::cout << "\nWavelength wl (um): " << wl_ << std::endl;
		std::cout << "Refractive index: " << GSL_REAL(v_1_) << ", " << GSL_IMAG(v_1_) << std::endl;
		std::cout << "Grating height a (um): " << g_.totalHeight() << std::endl;
		std::cout << "Grating period (um): " << g_.period() << std::endl;
		std::cout << "Number of layers: " << numLayers_ << std::endl;
		for(int m=1; m<M_; ++m) {
			std::cout << "   y_" << m << " = " << y_[m] << std::endl;
		}

		std::cout << "\nbeta1_n:" << std::endl;
		for(int i=0; i<twoNp1_; ++i) {
			std::cout << i - N_ << ":\t" << GSL_REAL(beta1_[i]) << "\t\t" << GSL_IMAG(beta1_[i]) << std::endl;
		}

		std::cout << "\nbetaM_n:" << std::endl;
		for(int i=0; i<twoNp1_; ++i) {
			std::cout << i - N_ << ":\t" << GSL_REAL(betaM_[i]) << "\t\t" << GSL_IMAG(betaM_[i]) << std::endl;
		}
	}

	// 3. Recursive computation of S-matrix below each layer.
	/////////////////////////////////////////////////////////////

	gsl_complex one = gsl_complex_rect(1,0);
	gsl_complex zero = gsl_complex_rect(0,0);

	// Handle first layer separately, as a special case.
	PEResult::Code status = computeTMatrixBelowLayer(2, printDebugOutput);
	if(status != PEResult::Success)
		return status;

	timing_[3] = time_;
	time_ = omp_get_wtime();
	timing_[3] = time_ - timing_[3];	// timing_[3]: numerical integration.

	// For the first layer, we have Zinv_ = T11_.
	// S12_ = T21_ Zinv_^{-1}
	// S22_ = Zinv_^{-1}
	///////////////
	gsl_matrix_complex_memcpy(Zinv_, T11_);
	int signnum;
	if(gsl_linalg_complex_LU_decomp(Zinv_, permutation_, &signnum) != GSL_SUCCESS)
		return PEResult(PEResult::AlgebraFailure);
	// S22_ = Zinv_^{-1}, ie, Z.
	if(gsl_linalg_complex_LU_invert(Zinv_, permutation_, S22_) != GSL_SUCCESS)
		return PEResult(PEResult::AlgebraFailure);
	// S12_ = T21_ Zinv_^{-1} = T21_ S22_.
	if(gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, T21_, S22_, zero, S12_) != GSL_SUCCESS)
		return PEResult(PEResult::AlgebraFailure);

	timing_[4] = time_;
	time_ = omp_get_wtime();
	timing_[4] = time_ - timing_[4];	// timing_[4]: matrix calcs.


	// First layer done. Handle subsequent layers
	for(int m=3; m<M_; ++m) {

		time_ = omp_get_wtime();

		status = computeTMatrixBelowLayer(m, printDebugOutput);
		if(status != PEResult::Success) return status;

		timing_[3] += omp_get_wtime() - time_;

		time_ = omp_get_wtime();

		// Compute Zinv_ = T11_ + T12_ S12_.
		gsl_matrix_complex_memcpy(Zinv_, T11_);
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, T12_, S12_, one, Zinv_);

		// Invert Zinv_...
		if(gsl_linalg_complex_LU_decomp(Zinv_, permutation_, &signnum) != GSL_SUCCESS) return PEResult(PEResult::AlgebraFailure);
		if(gsl_linalg_complex_LU_invert(Zinv_, permutation_, Z_) != GSL_SUCCESS) return PEResult(PEResult::AlgebraFailure);

		// S12 = (T21 + T22 S12) Z
		gsl_matrix_complex_memcpy(workMatrix_, T21_);
		if(gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, T22_, S12_, one, workMatrix_) != GSL_SUCCESS)
			return PEResult(PEResult::AlgebraFailure);
		if(gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, workMatrix_, Z_, zero, S12_) != GSL_SUCCESS)
			return PEResult(PEResult::AlgebraFailure);
		// S22 = S22 Z
		gsl_matrix_complex_memcpy(workMatrix_, S22_);
		if(gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, workMatrix_, Z_, zero, S22_) != GSL_SUCCESS)
			return PEResult(PEResult::AlgebraFailure);

		timing_[4] += omp_get_wtime() - time_;
	}

	// 4.  Calculate B_n^M from center column of S matrix * exp(...).
	//////////////////////////////////////////////
	time_ = omp_get_wtime();

	computeBMFromSMatrix();

	timing_[5] = time_;
	time_ = omp_get_wtime();
	timing_[5] = time_ - timing_[5];	// timing_[5]: compute BM_ from S-matrix.

	if(printDebugOutput) {
		std::cout << "\nBM_:" << std::endl;
		for(int i=0; i<twoNp1_; ++i) {
			std::cout << i - N_ << ":\t" << GSL_REAL(BM_[i]) << "\t\t" << GSL_IMAG(BM_[i]) << std::endl;
		}
	}

	// 8. Now we have BM_. Compute efficiency and put into result structure.
	////////////////////////////////////////
	
	PEResult result(N_);
	result.wavelength = wl_;
	result.incidenceDeg = incidenceDeg;
	
	double effSum = 0;
	for(int i=0; i<twoNp1_; ++i) {
		// int n = i - N_;
		result.eff[i] =  gsl_complex_abs2(BM_[i])*GSL_REAL(betaM_[i])/GSL_REAL(betaM_[N_]);
		// is this a non-propagating order?  Then the real part of beta2_n will be exactly 0, so the efficiency will come out as 0.
		effSum += result.eff[i];
	}
	
	timing_[6] = time_;
	time_ = omp_get_wtime();
	timing_[6] = time_ - timing_[6];	// timing_[6]: calculate efficiencies from BM_
	
	if(measureTiming_) {
		std::cout << "Timing Profile:" << std::endl;
		std::cout << "   Allocate Memory: " << timing_[0] << std::endl;
		std::cout << "   Look up refractive index: " << timing_[1] << std::endl;
		std::cout << "   Compute alpha, beta values and layers: " << timing_[2] << std::endl;
		std::cout << "   Numerically integrating trial solutions: " << timing_[3] << std::endl;
		std::cout << "   Matrix operations: " << timing_[4] << std::endl;
		std::cout << "   Computing Rayleigh coeffients B_n: " << timing_[5] << std::endl;
		std::cout << "   Compute and package efficiencies: " << timing_[6] << std::endl;
		time_ = timing_[0] + timing_[1] + timing_[2] + timing_[3] + timing_[4] + timing_[5] + timing_[6];
		std::cout << "   Total (solver) time: " << time_ << std::endl << std::endl;
	}

	if(printDebugOutput) {
		std::cout << "Sum of efficiencies (should be <= 1): " << effSum << std::endl;
	}
	
	return result;
}

gsl_complex PESolver::complex_sqrt_upperComplexPlane(gsl_complex z) {

	gsl_complex w = gsl_complex_sqrt(z); // returns w in the right half of complex plane.
	
	if(GSL_IMAG(w) < 0) {
		// add Pi to the phase. This will still have the same square, since e^{i(pi+pi)} = 1, but it will be in the upper complex plane now.
		w = gsl_complex_mul(w, gsl_complex_polar(1, M_PI));
	}
	
	return w;
}

PEResult::Code PESolver::computeGratingExpansion(double y, gsl_complex* k2) const {
	
	// wave number in free space: k_2 = v_2 * w / c.  v_2 = 1 in empty space, so k_2 = 2pi / wl.
	gsl_complex k_M = gsl_complex_rect(2 * M_PI / wl_, 0);
	// wave number in the substrate: k_1 = v_1 * w / c = v_1 * 2pi / wl = v_1 * k_2.
	gsl_complex k_1 = gsl_complex_mul(v_1_, k_M);
	// wave number in the coating layer:
	gsl_complex k_c = gsl_complex_mul(v_c_, k_M);
	
	// square them to get k^2_M, k^2_1 and k^2_c:
	gsl_complex k2_M = gsl_complex_mul(k_M, k_M);
	gsl_complex k2_1 = gsl_complex_mul(k_1, k_1);
	gsl_complex k2_c = gsl_complex_mul(k_c, k_c);

	double stepsX[PEG_MAX_PROFILE_CROSSINGS];
	gsl_complex stepsK2[PEG_MAX_PROFILE_CROSSINGS];

	// Compute multistep function from grating:
	int numSteps = g_.computeK2StepsAtY(y, k2_M, k2_1, k2_c, stepsX, stepsK2);
	if(numSteps < 1)
		return PEResult::InvalidGratingFailure;

	computeGratingExpansion(stepsX, stepsK2, numSteps, k2);
	
	return PEResult::Success;
}

// Now unused:
PEResult::Code PESolver::integrateTrialSolutionAlongY(gsl_vector_complex* u, gsl_vector_complex* uprime, double yStart, double yEnd) {

	// fill starting conditions from u, uprime
	double* w = new double[eightNp4_];

	for(int i=0; i<twoNp1_; ++i) {
		gsl_complex u_n = gsl_vector_complex_get(u, i);
		w[2*i] = GSL_REAL(u_n);
		w[2*i + 1] = GSL_IMAG(u_n);
		
		gsl_complex uprime_n = gsl_vector_complex_get(uprime, i);
		w[2*i + fourNp2_] = GSL_REAL(uprime_n);
		w[2*i + fourNp2_ + 1] = GSL_IMAG(uprime_n);
	}
	
	// integrate it:
	PEResult::Code status = integrateTrialSolutionAlongY(w, yStart, yEnd);
	
	// copy results back into u, uprime
	for(int i=0; i<twoNp1_; ++i) {
		gsl_vector_complex_set(u, i, gsl_complex_rect(w[2*i], w[2*i + 1]));
		gsl_vector_complex_set(uprime, i, gsl_complex_rect(w[2*i + fourNp2_], w[2*i + fourNp2_ + 1]));
	}
	
	delete [] w;
	return status;
}

PEResult::Code PESolver::integrateTrialSolutionAlongY(double *w, double yStart, double yEnd) {
	// define ode solving system, with our function to evaluate dw/dy, the Jacobian, and 8*N_+4 components.
	gsl_odeiv2_system odeSys = {odeFunctionCB, odeJacobianCB, eightNp4_, this};

	// initial starting step in y: choose grating height / 200.
	double hStart = (yEnd - yStart)/200;

	// setup driver
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_standard_new (&odeSys, gsl_odeiv2_step_msadams, hStart, integrationTolerance_, integrationTolerance_, 0.5, 0.5);	// Variable-coefficient linear multistep Adams method in Nordsieck form. Uses explicit Adams-Bashforth (predictor) and implicit Adams-Moulton (corrector) methods in P(EC)^m functional iteration mode.
	//	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_standard_new (&odeSys, gsl_odeiv2_step_rkf45, hStart, integrationTolerance_, integrationTolerance_, 0.5, 0.5); // Explicit embedded Runge-Kutta-Fehlberg (4, 5) method.

	// run it: integrate from y = yStart to y=yEnd.
	double y = yStart;
	int status = gsl_odeiv2_driver_apply (d, &y, yEnd, w);
	//	int status = gsl_odeiv2_driver_apply_fixed_step(d, &y, hStart, 200, w);

	gsl_odeiv2_driver_free(d);

	if (status != GSL_SUCCESS) {
		std::cout << "ODE: Integration failure: Code: " << status << std::endl;
		return PEResult::ConvergenceFailure;
	}

	return PEResult::Success;
}

int PESolver::odeFunction(double y, const double w[], double dwdy[]) {
	
	// w contains the last values of u_n{re, im} and u'_n{re, im}, in that order.
	// need to compute f = dw/dy = u'_n{re, im} followed by u''_n{re, im}
	
	// get k2_n at this y value.
	gsl_complex* localK2 = k2ForCurrentThread();
	if(computeGratingExpansion(y, localK2) != PEResult::Success) {
		std::cout << "ODE: Function Error: Cannot compute grating expansion at y = " << y << std::endl;
		return GSL_EBADFUNC;	// can't calculate here. Invalid profile? y above the profile height?
	}
	
	// fourNp2 divides the top and bottom of the arrays w, dwdy.  Top of w is u; bottom of w is u' = v.   Top of dwdy is u';  bottom of dwdy is u''.
	// total size is (2N+1)x2x2, ie: 8N+4.
	for(int i=0; i<fourNp2_; ++i) {
		// working on computing u'_n [top of dwdy array]. Just copy from u' = [bottom half of w array]
		dwdy[i] = w[i + fourNp2_];
	}

	for(int i=fourNp2_; i<eightNp4_; i+=2) {
		// working on computing u''_n [bottom of f array].  We'll handle real and imaginary components at once.
		int n = (i - fourNp2_)/2 - N_;
		gsl_complex upp_n = gsl_complex_rect(0,0);	// initialize sum to 0.

		// alpha^2_n:
		double alpha2 = alpha_[n + N_];
		alpha2 *= alpha2;

		// loop over m.  Retrieving u_n values from top of array.
		for(int j=0; j<fourNp2_; j+=2) { // also handling real and imaginary components at once.
			int m = j/2 - N_;

			gsl_complex u_m = gsl_complex_rect(w[j], w[j+1]);

			gsl_complex M_nm = gsl_complex_rect(0,0);
			if(n-m >= -N_ && n-m <= N_)
				M_nm = gsl_complex_mul_real(localK2[n-m + N_], -1.0);	// -k^2_{n-m}
			if(n == m)
				M_nm = gsl_complex_add_real(M_nm, alpha2);

			upp_n = gsl_complex_add(upp_n, gsl_complex_mul(M_nm, u_m));
		}

		dwdy[i] = GSL_REAL(upp_n);
		dwdy[i+1] = GSL_IMAG(upp_n);
	}

	// Debugging output:
	//////////////////////
	//	double sum;
	//	for(int i=0; i<fourNp2_; i+=2) {
	//		gsl_complex u_n = gsl_complex_rect(w[i], w[i+1]);
	//		sum += GSL_REAL(gsl_complex_mul(u_n, gsl_complex_conjugate(u_n)));
	//	}
	//	std::cout << "y: " << y << " Power: " << sum << std::endl;
	
	return GSL_SUCCESS;
}


int PESolver::odeJacobian(double y, const double w[], double *dfdw, double dfdy[])
{
	(void)w;

	// y is the indep. variable.
	// u[] contains the current solution u(y): contains u1=u followed by u2=u'.
	// dfdw is the Jacobian; in row-major order.
	// dfdy is the partial time deriv; no idea how to calc.

	// get k2_n at this y value.
	gsl_complex* localK2 = k2ForCurrentThread();
	if(computeGratingExpansion(y, localK2) != PEResult::Success) {
		std::cout << "ODE: Jacobian Error: Cannot compute grating expansion at y = " << y << std::endl;
		return GSL_EBADFUNC;	// can't calculate here. Invalid profile? y above the profile height?
	}

	// fourNp2 divides the top and bottom of the arrays w, f.  Top of w is u; bottom of w is u' = v.   Top of f is u' = v;  bottom of f is v' = u''.
	// total size is (2N+1)x2x2, ie: 8N+4.

	// clear the jacobian
	memset(dfdw, 0, eightNp4_*eightNp4_*sizeof(double));
	memset(dfdy, 0, eightNp4_*sizeof(double));

	// go through top rows of jac [i=0,fourNp2]. Set ident. matrix in upper right-hand block.
	for(int i=0; i<fourNp2_; ++i) {
		dfdw[i*eightNp4_+fourNp2_+i] = 1.0;		// set at index dfdw(i, fourNp2+i).
	}

	// go throgh lower rows of jac [i=fourNp2, eightNp4]. In lower left-hand block, set 4x4 submatrices at once, so go by i+=2.
	for(int i=fourNp2_; i<eightNp4_; i+=2) {
		int n = (i - fourNp2_)/2 - N_;	// n index [row] from -N to N.

		// alpha^2_n:
		double alpha2 = alpha_[n + N_];
		alpha2 *= alpha2;

		// loop over m [col].
		for(int j=0; j<fourNp2_; j+=2) { // also handling real and imaginary components at once.
			int m = j/2 - N_;	// m index [col] from -N to N.

			// get M_{nm}: -(k^2)_{n-m}(y) + alpha^2_n \delta_{nm}
			gsl_complex M_nm = gsl_complex_rect(0,0);
			if(n-m >= -N_ && n-m <= N_)
				M_nm = gsl_complex_mul_real(localK2[n-m + N_], -1.0);	// k2_ ranges from -N to N.
			if(n == m)
				M_nm = gsl_complex_add_real(M_nm, alpha2);

			// set 4x4 matrix here: [M_re, -M_im; M_im, M_re] at (i,j), (i, j+1); (i+1, j), (i+1, j+1)
			dfdw[i*eightNp4_ + j] = GSL_REAL(M_nm);
			dfdw[i*eightNp4_ + j + 1] = -GSL_IMAG(M_nm);
			dfdw[(i+1)*eightNp4_ + j] = GSL_IMAG(M_nm);
			dfdw[(i+1)*eightNp4_ + j + 1] = GSL_REAL(M_nm);
		}
	}

	/// \todo IMPORTANT! Leaving dfdy = 0 for now. This is only true in case of rectangular grating...

	return GSL_SUCCESS;
}

gsl_complex* PESolver::k2ForCurrentThread() {
	/// Return based on OpenMP current thread.
	return k2_[omp_get_thread_num()];
}

double PESolver::conditionNumber(const gsl_matrix_complex *A)
{
	// according to http://en.wikipedia.org/wiki/Condition_number,
	// When using the 2-norm, cond(A) = max(Singular Values) / min(Singular Values)

	// compute Singular Value Decomposition

	// Make copy of A:
	gsl_matrix_complex* A2 = gsl_matrix_complex_alloc(A->size1, A->size2);
	gsl_matrix_complex_memcpy(A2, A);

	// Create other matrices:
	gsl_matrix_complex* Ainv = gsl_matrix_complex_alloc(A->size1, A->size2);
	gsl_permutation* permut = gsl_permutation_alloc(A->size1);

	int signum;
	if(gsl_linalg_complex_LU_decomp(A2, permut, &signum) != GSL_SUCCESS) {
		std::cout << "Warning: LU failed!" << std::endl;
		return 0;
	}

	if(gsl_linalg_complex_LU_invert(A2, permut, Ainv) != GSL_SUCCESS) {
		std::cout << "Warning: Invert failed!" << std::endl;
		return 0;
	}

	// how to get the norm of a complex matrix?

	gsl_matrix_complex_free(A2);
	gsl_matrix_complex_free(Ainv);
	gsl_permutation_free(permut);

	return 1;
}


int PESolver::linalg_LU_complex_solve(const gsl_matrix_complex *LU, const gsl_permutation *P, const gsl_matrix_complex *B, gsl_matrix_complex *X)  {
	int n = LU->size1;
	int status;
	for(int j=0; j<n; ++j) {
		gsl_vector_complex_view x = gsl_matrix_complex_column(X, j);
		gsl_vector_complex_const_view b = gsl_matrix_complex_const_column(B, j);

		status = gsl_linalg_complex_LU_solve(LU, P, &(b.vector), &(x.vector));
		if(status != GSL_SUCCESS)
			return status;
	}
	return GSL_SUCCESS;
}

void PESolver::setIntegrationStartingValues(double *w, int j, int m)
{
	// set all to 0
	memset(w, 0, (eightNp4_)*sizeof(double));

	bool secondRound = false;
	if(j >= twoNp1_) {
		secondRound = true;
		j -= twoNp1_;
	}

	// set u[j] = 1.  Multiplication by 2 is due to {re,im}.
	w[2*j] = 1.0;

	gsl_complex uprime = gsl_complex_mul_imag(m == 1 ? beta1_[j] : betaM_[j], secondRound ? 1 : -1);
	w[fourNp2_ + 2*j] = GSL_REAL(uprime);
	w[fourNp2_ + 2*j+1] = GSL_IMAG(uprime);
}

void PESolver::computeAlphaAndBeta(double incidenceDeg)
{
	// incidence angle in radians
	double theta_2 = incidenceDeg * M_PI / 180;
	// omega: angular frequency. Computed from wavelength in vacuum. w = 2 pi f = 2 pi c / wl
	// double omega = 2 * M_PI * M_c / wl;

	// wave number in free space: k_2 = v_2 * omega / c.  v_2 = 1 in empty space, so k_2 = 2pi / wl.
	double k_2 = 2 * M_PI / wl_;
	// wave number in the grating: k_1 = v_1 * w / c = v_1 * 2pi / wl.
	gsl_complex k_1 = gsl_complex_mul_real(v_1_, k_2);

	// Grating period:
	double d = g_.period();


#pragma omp parallel for num_threads(numThreads_)
	for(int i=0; i<twoNp1_; i++) {
		int n = i - N_;

		double alpha = k_2 * sin(theta_2) + 2 * M_PI * n / d;
		alpha_[i] = alpha;

		// beta2_: rayleigh expansion above grating.
		double k22minusAn2 = k_2*k_2 - alpha*alpha;
		if(k22minusAn2 >= 0)
			betaM_[i] = gsl_complex_rect(sqrt(k22minusAn2), 0);
		else
			betaM_[i] = gsl_complex_rect(0, sqrt(-k22minusAn2));

		// beta1_: rayleigh expansion inside grating
		gsl_complex k12minusAn2 = gsl_complex_sub_real(gsl_complex_mul(k_1, k_1), alpha*alpha);
		beta1_[i] = complex_sqrt_upperComplexPlane(k12minusAn2);

		//		// TEST:
		//		z = gsl_complex_mul(beta1_[i], beta1_[i]);
		//		if(fabs(GSL_REAL(z) - GSL_REAL(k12minusAn2)) > 1e-5 || fabs(GSL_IMAG(z) - GSL_IMAG(k12minusAn2)) > 1e-5) {
		//			std::cout << "Square root error!: " << GSL_REAL(k12minusAn2) << "," << GSL_IMAG(k12minusAn2) << "\t" << GSL_REAL(z) << "," << GSL_IMAG(z) << std::endl;
		//			exit(-1);
		//		}
	}
}

void PESolver::computeLayers()
{
	// How many layers do we need?
	double a = g_.totalHeight();

	double magicNumber = 7;	// should be ln(1e15). However, emperically this is not enough to maintain stability (ex: REIXS LEG).  7 = ln(1e3) seems stable for all tests so far.

	// How many layers to use? In order to keep size of exp(i betaM_{±N}) < 1e15 to avoid losing precision in double values compared with unity-size numbers.
	numLayers_ = std::max( fabs(GSL_IMAG(betaM_[0]))*a/magicNumber, fabs(GSL_IMAG(betaM_[2*N_]))*a/magicNumber );
	if(numLayers_ < 1)
		numLayers_ = 1;	// we need at least one layer, in addition to the substrate.

	M_ = numLayers_+2;

	y_ = new double[M_];	// To be consistent with the text, we number the lowest layer as 1.  y_[0] is therefore unused, so that we can use y_m = y_[m], with lowest m=1, highest m=M-1.

	for(int m=1; m<M_; ++m) {
		y_[m] = double(m-1)/numLayers_*a;
	}
}

void PESolver::computeBMFromSMatrix()
{
	double a = g_.totalHeight();

#pragma omp parallel for num_threads(numThreads_)
	for(int i=0; i<twoNp1_; i++) {
		BM_[i] = gsl_complex_mul(
					gsl_matrix_complex_get(S12_, i, N_),
					gsl_complex_exp(gsl_complex_mul_imag(gsl_complex_add(betaM_[i],
																		 betaM_[N_]),
														 -a)));
	}
}

PEResult::Code PESolver::computeTMatrixBelowLayer(int m, bool printDebugOutput)
{
	bool integrationFailureOccurred = false;

	// We now need 2*(2N+1) trial solutions.  j will be the loop index over p, but ranging from [0,4*N+1].
#pragma omp parallel for num_threads(numThreads_) schedule(dynamic)
	for(int j=0; j<fourNp2_; ++j) {

		// Get a [u,uprime] vector to work with for this trial solution.
		double* w = wVectorForP(j);

		// set its initial values at y_[m-1]
		setIntegrationStartingValues(w, j, m-1);

		////////////////////////////
		if(printDebugOutput && omp_get_thread_num() == 0) {
			std::cout << "Initial value u_{p=" << j-N_ << "}(yStart):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w[2*n] << "," << w[2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << "Initial value u'_{p=" << j-N_ << "}(yStart):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w[fourNp2_ + 2*n] << "," << w[fourNp2_ + 2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << std::endl;
		}
		//////////////////////////

		// Integrate from y_[m-1] to y_[m].
		PEResult::Code status = integrateTrialSolutionAlongY(w, y_[m-1], y_[m]);

		////////////////////////////
		if(printDebugOutput && omp_get_thread_num() == 0) {
			std::cout << "Final value u_{p=" << j-N_ << "}(yEnd):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w[2*n] << "," << w[2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << "Final value u'_{p=" << j-N_ << "}(yEnd):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w[fourNp2_ + 2*n] << "," << w[fourNp2_ + 2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << std::endl;
		}
		//////////////////////////


		if(status != PEResult::Success)
			integrationFailureOccurred = true;
		else {
			// Fill T-matrix at this column. If j >= 2*N+1, we're dealing with the right-side blocks T12_, T22_.
			if(j >= twoNp1_) {
				int jj = j - twoNp1_;

				// loop over orders (n). i is the loop index, ranging from [0, 2*N_]
				for(int i=0; i<twoNp1_; ++i) {
					gsl_complex* u_ij = (gsl_complex*)(w + 2*i);
					gsl_complex* uprime_ij = (gsl_complex*)(w + fourNp2_ + 2*i);
					gsl_complex temp = gsl_complex_div(*uprime_ij, gsl_complex_mul_imag(betaM_[i], 1)); // = u'_ij / (i*betaM_n)

					// T12_ij = 0.5(u_ij - u'_ij / (i*betaM_n) )
					// T22_ij = 0.5(u_ij + u'_ij / (i*betaM_n) )
					gsl_matrix_complex_set(T12_, i, jj, gsl_complex_mul_real(gsl_complex_sub(*u_ij, temp), 0.5));
					gsl_matrix_complex_set(T22_, i, jj, gsl_complex_mul_real(gsl_complex_add(*u_ij, temp), 0.5));
				}
			}
			// Otherwise we're filling T11, T21. (left-side blocks).
			else {
				// loop over orders (n). i is the loop index, ranging from [0, 2*N_]
				for(int i=0; i<twoNp1_; ++i) {
					gsl_complex* u_ij = (gsl_complex*)(w + 2*i);
					gsl_complex* uprime_ij = (gsl_complex*)(w + fourNp2_ + 2*i);
					gsl_complex temp = gsl_complex_div(*uprime_ij, gsl_complex_mul_imag(betaM_[i], 1)); // = u'_ij / (i*betaM_n)

					// T12_ij = 0.5(u_ij - u'_ij / (i*betaM_n) )
					// T22_ij = 0.5(u_ij + u'_ij / (i*betaM_n) )
					gsl_matrix_complex_set(T11_, i, j, gsl_complex_mul_real(gsl_complex_sub(*u_ij, temp), 0.5));
					gsl_matrix_complex_set(T21_, i, j, gsl_complex_mul_real(gsl_complex_add(*u_ij, temp), 0.5));
				}
			}
		}
	}

	if(integrationFailureOccurred)
		return PEResult::ConvergenceFailure;
	else
		return PEResult::Success;
}

// computes the fourier expansion of the multistep function given by values stepsK2 at x-axis locations stepsX, and stores in k2.
void PESolver::computeGratingExpansion(const double *stepsX, const gsl_complex *stepsK2, int numSteps, gsl_complex *k2) const
{
	// period:
	double d = g_.period();
	double K = 2*M_PI/d;

	/// \warning assumes numSteps is in [1, PEG_MAX_PROFILE_CROSSINGS]

	// Optimization for numSteps = 1: f_n = 0 (n!=0).   f_0 = stepsK2[0].
	if(numSteps == 1) {
		for(int i=0; i<twoNp1_; ++i)
			k2[i] = gsl_complex_rect(0,0);
		k2[N_] = stepsK2[0];
		return;
	}

	// compute sigma values at crossings:
	// sigma[p] = stepsK2[p+1] - stepsK2[p] for p<numSteps-1; sigma[numSteps-1]=sigma[0]-sigma[numSteps-1]
	gsl_complex sigma[PEG_MAX_PROFILE_CROSSINGS];
	for(int p=0;p<numSteps-1; ++p)
		sigma[p] = gsl_complex_sub(stepsK2[p+1], stepsK2[p]);
	sigma[numSteps-1] = gsl_complex_sub(stepsK2[0], stepsK2[numSteps-1]);

	// loop over fourier indexes
	for(int i=0; i<twoNp1_; ++i) {
		int n = i - N_;

		if(n == 0) {
			gsl_complex f0 = gsl_complex_mul_real(stepsK2[0], d);
			for(int p=0; p<numSteps; ++p)
				f0 = gsl_complex_sub(f0, gsl_complex_mul_real(sigma[p], stepsX[p]));
			k2[i] = gsl_complex_div_real(f0, d);
		}

		else {
			gsl_complex fn = gsl_complex_rect(0,0);
			for(int p=0; p<numSteps; ++p) {
				double nKx = n*K*stepsX[p];
				fn = gsl_complex_add(fn, gsl_complex_mul(sigma[p], gsl_complex_rect(sin(nKx), cos(nKx))));
			}
			k2[i] = gsl_complex_div_real(fn, -2*M_PI*n);
		}
	}

	// that's it!
	////////////////////////////

// This is potentially an optimization that avoids the inner loops for the common case when numSteps = 2. However, the performance difference isn't enough to be noticed above the random noise.

//	gsl_complex k2_2 = stepsK2[0];
//	gsl_complex k2_1 = stepsK2[1];


//	// sigma1 (s1) and sigma2 (s2) are defined as: s_p = k^2_(p+1) - k^2_p... Except at p = P (max): s_P = k^2_1 - k^2_P.  We only have two crossings, so P = 2.  p = [1,2].
//	// k^2_p is the refractive index value on the left side of x_p crossing.
//	gsl_complex s1 = gsl_complex_sub(k2_1, k2_2);
//	gsl_complex s2 = gsl_complex_sub(k2_2, k2_1);

//	// for shorter calculations, define K as the grating number, 2*pi/d.  We just calculate it once.
//	double K = 2*M_PI/d;

//	// loop over all values of n, from -N to N.  See paper notes for formula for k^2_n [Page 5,6]
//	for(int i=0; i<twoNp1_; ++i) {
//		int n = i - N_;

//		if(n == 0) {
//			gsl_complex kk = gsl_complex_mul_real(k2_2, d);
//			kk = gsl_complex_sub(kk, gsl_complex_mul_real(s1, x1));
//			kk = gsl_complex_sub(kk, gsl_complex_mul_real(s2, x2));

//			kk = gsl_complex_mul_real(kk, 1.0/d);
//			k2[i] = kk;
//		}
//		else {
//			double t1 = n*K*x1, t2 = n*K*x2;
//			gsl_complex kk = gsl_complex_mul(s1, gsl_complex_rect(sin(t1), cos(t1)));
//			kk = gsl_complex_add(kk, gsl_complex_mul(s2, gsl_complex_rect(sin(t2), cos(t2))));

//			kk = gsl_complex_mul_real(kk, -1.0/2/M_PI/n);
//			k2[i] = kk;
//		}
//	}
}

