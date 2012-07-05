#include "PESolver.h"

#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <omp.h>

// For debug output only:
#include <iostream>

/// Speed of light in um/s.
#define M_c 2.99792458e14

PESolver::PESolver(const PEGrating& grating, const PEMathOptions& mo, int numThreads, bool measureTiming)
	: g_(grating)
{
	numThreads_ = numThreads;
	measureTiming_ = measureTiming;
	timing_[0] = omp_get_wtime();
	
	// set math options.
	N_ = mo.N;
	twoNp1_ = 2*N_ + 1;
//	niy_ = mo.niy;
	
	// allocate matrices and vectors
	u_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	uprime_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	iBeta2Diag_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	T_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	gsl_matrix_complex_set_zero(iBeta2Diag_);
	Vincident_ = gsl_vector_complex_alloc(twoNp1_);
	gsl_vector_complex_set_zero(Vincident_);
	A1_ = gsl_vector_complex_alloc(twoNp1_);
	permutation_ = gsl_permutation_alloc(twoNp1_);
	B2_ = gsl_vector_complex_alloc(twoNp1_);
	
	alpha_ = new double[twoNp1_];
	beta2_ = new gsl_complex[twoNp1_];
	beta1_ = new gsl_complex[twoNp1_];
	
	// we need one k2_ array for each thread, since they will be used simultaneously.
	k2_ = new gsl_complex*[numThreads_];
	for(int i=0; i<numThreads_; ++i)
		k2_[i] = new gsl_complex[twoNp1_];
	
	timing_[1] = omp_get_wtime();
	
}


PESolver::~PESolver() {
	gsl_matrix_complex_free(u_);
	gsl_matrix_complex_free(uprime_);
	gsl_matrix_complex_free(iBeta2Diag_);
	gsl_matrix_complex_free(T_);
	gsl_vector_complex_free(Vincident_);
	gsl_vector_complex_free(A1_);
	gsl_permutation_free(permutation_);
	gsl_vector_complex_free(B2_);
	
	delete [] alpha_;
	delete [] beta2_;
	delete [] beta1_;
	for(int i=0; i<numThreads_; ++i)
		delete [] k2_[i];
	delete [] k2_;
}

/// \todo Imp.
PEResult PESolver::getEff(double incidenceDeg, double wl, bool printDebugOutput) {
	
	int errCode;
	
	// 1. Setup incidence variables and constants
	/////////////////////////////////////////
	
	timing_[2] = omp_get_wtime();
	// Set this as the current wavelength
	wl_ = wl;
	// get refractive index at wavelength
	v_1_ = g_.refractiveIndex(wl_);
	
	// incidence angle in radians
	double theta_2 = incidenceDeg * M_PI / 180;
	// omega: angular frequency. Computed from wavelength in vacuum. w = 2 pi f = 2 pi c / wl
	// double omega = 2 * M_PI * M_c / wl;
	
	// wave number in free space: k_2 = v_2 * omega / c.  v_2 = 1 in empty space, so k_2 = 2pi / wl.
	double k_2 = 2 * M_PI / wl;
	// wave number in the grating: k_1 = v_1 * w / c = v_1 * 2pi / wl.
	gsl_complex k_1 = gsl_complex_mul_real(v_1_, k_2);
	// grating period:
	double d = g_.period();
	
	// temporary complex variable for doing calculations
	gsl_complex z;
	// the height of the grating (in um)
	double a = g_.height();
	
	if(printDebugOutput) {
		std::cout << "\nWavelength wl (um):" << wl_ << std::endl;
		std::cout << "Grating height a (um): " << a << std::endl;
		std::cout << "Grating period (um):" << d << std::endl;
	}
	
	// 2. compute all alpha_n and beta1_n, beta2_n.  (This can be a parallel loop)
	///////////////////////////////////
	timing_[3] = omp_get_wtime();
	
#pragma omp parallel for num_threads(numThreads_)
	for(int i=0; i<twoNp1_; i++) {
		int n = i - N_;
		
		double alpha = k_2 * sin(theta_2) + 2 * M_PI * n / d;
		alpha_[i] = alpha;
		
		// beta2_: rayleigh expansion above grating.
		double k22minusAn2 = k_2*k_2 - alpha*alpha;
		if(k22minusAn2 >= 0)
			beta2_[i] = gsl_complex_rect(sqrt(k22minusAn2), 0);
		else
			beta2_[i] = gsl_complex_rect(0, sqrt(-k22minusAn2));
		
		// beta1_: rayleigh expansion inside grating
		gsl_complex k12minusAn2 = gsl_complex_sub_real(gsl_complex_mul(k_1, k_1), alpha*alpha);
		beta1_[i] = complex_sqrt_upperComplexPlane(k12minusAn2);
	}
	
	if(printDebugOutput) {
		std::cout << "\nbeta1_n:" << std::endl;
		for(int i=0; i<twoNp1_; ++i) {
			std::cout << i - N_ << ":\t" << GSL_REAL(beta1_[i]) << "\t\t" << GSL_IMAG(beta1_[i]) << std::endl;
		}
	}
	
	if(printDebugOutput) {
		std::cout << "\nbeta2_n:" << std::endl;
		for(int i=0; i<twoNp1_; ++i) {
			std::cout << i - N_ << ":\t" << GSL_REAL(beta2_[i]) << "\t\t" << GSL_IMAG(beta2_[i]) << std::endl;
		}
	}
	
	// 3. set up trial solutions at y=0, and then integrate their ODEs from y=0 to y=a.
	////////////////////////////////////////
	timing_[4] = omp_get_wtime();
	
	// u_ should be the identify matrix
	gsl_matrix_complex_set_identity(u_);
	// zero uprime_; we will set the components we need to below.
	gsl_matrix_complex_set_zero(uprime_);

	////////////////////////////////
	if(printDebugOutput) {
		std::cout << "Example grating expansion at y = 0.02" << std::endl;
		gsl_complex* localK2 = k2ForCurrentThread();
		computeGratingExpansion(0.01, localK2);

		// wave number in free space: k_2 = v_2 * w / c.  v_2 = 1 in empty space, so k_2 = 2pi / wl.
		gsl_complex k_2 = gsl_complex_rect(2 * M_PI / wl_, 0);
		// wave number in the grating: k_1 = v_1 * w / c = v_1 * 2pi / wl = v_1 * k_2.
		gsl_complex k_1 = gsl_complex_mul(v_1_, k_2);

		// square them to get k^2_2 and k^2_1:
		gsl_complex k2_2 = gsl_complex_mul(k_2, k_2);
		gsl_complex k2_1 = gsl_complex_mul(k_1, k_1);

		std::cout << "k2_2" << GSL_REAL(k2_2) << "," << GSL_IMAG(k2_2) << std::endl;
		std::cout << "k2_1" << GSL_REAL(k2_1) << "," << GSL_IMAG(k2_1) << std::endl;

		for(int i=0; i<twoNp1_; ++i) {
			std::cout << "   " << GSL_REAL(localK2[i]) << "," << GSL_IMAG(localK2[i]);
		}
		std::cout << std::endl;
	}
	/////////////////////////////////
	
	// Loop over all trial solutions p:  [row: n.  col: p]
	// !! This is the time-consuming step, but all of the trial solutions are independent. Every loop iteration is independent, therefore we use an OpenMP parallel-for loop here. We might expect that the integration time will increase or decrease monotonically over p, so we use a cyclic partition for load-balancing. !!
	bool integrationFailureOccurred = false;
#pragma omp parallel for num_threads(numThreads_) schedule(static,1)
	for(int i=0; i<twoNp1_; ++i) {
		// int p = i - N_;
		
		// insert initial values at y=0.  u(0) = 1* \delta_{np}.  (Already set: identity matrix)
		// u'(0) = -i * beta^(1)_n \delta_{np}
		gsl_matrix_complex_set(uprime_, i, i, gsl_complex_mul_imag(beta1_[i], -1));
		
		// create vector views: the column vectors for this trial solution
		gsl_vector_complex_view u = gsl_matrix_complex_column(u_, i);
		gsl_vector_complex_view uprime = gsl_matrix_complex_column(uprime_, i);

		////////////////////////////
		if(printDebugOutput && omp_get_thread_num() == 0) {
			std::cout << "Initial value u_{p=" << i-N_ << "}(0):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << GSL_REAL(gsl_vector_complex_get(&u.vector, n)) << "," << GSL_IMAG(gsl_vector_complex_get(&u.vector, n)) << "    ";
			std::cout << std::endl;
			std::cout << "Initial value u'_{p=" << i-N_ << "}(0):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << GSL_REAL(gsl_vector_complex_get(&uprime.vector, n)) << "," << GSL_IMAG(gsl_vector_complex_get(&uprime.vector, n)) << "    ";
			std::cout << std::endl;
			std::cout << std::endl;
		}
		//////////////////////////
		
		// need to integrate u and uprime along y, using d^2 u/dy^2 = M(y) u
		PEResult::Code err = integrateTrialSolutionAlongY(&u.vector, &uprime.vector);
		if(err != PEResult::Success)
			integrationFailureOccurred = true;	// We used to exit the loop here, but that is not allowed in an OpenMP parallel for.  Set this flag instead.  It is only written, but not read, by different threads, so we should be OK.

		//////////////////////////
		if(printDebugOutput && omp_get_thread_num() == 0) {
			if(err == PEResult::Success) {
				std::cout << "Final value u_{p=" << i-N_ << "}(" << g_.height() << "):" <<std::endl;
				std::cout << "     ";
				for(int n=0; n<twoNp1_; ++n)
					std::cout << GSL_REAL(gsl_vector_complex_get(&u.vector, n)) << "," << GSL_IMAG(gsl_vector_complex_get(&u.vector, n)) << "    ";
				std::cout << std::endl;
				std::cout << "Final value u'_{p=" << i-N_ << "}(" << g_.height() << "):" <<std::endl;
				std::cout << "     ";
				for(int n=0; n<twoNp1_; ++n)
					std::cout << GSL_REAL(gsl_vector_complex_get(&uprime.vector, n)) << "," << GSL_IMAG(gsl_vector_complex_get(&u.vector, n)) << "    ";
				std::cout << std::endl;
				std::cout << std::endl;
			}
			else
				std::cout << "Integration Failure on p=" << i-N_ << std::endl;
		}
		/////////////////////////////
	}
	if(integrationFailureOccurred)
		return PEResult(PEResult::ConvergenceFailure);
	
	// Now we have u(a) and u'(a) in (u, uprime) and in the columns of the actual (u_, uprime_) matrices.
	
	// 4. Need to calculate the T matrix that maps the grating input to grating output (in the basis expansion)
	/////////////////////////////
	timing_[5] = omp_get_wtime();
	
	// Need a diagonal matrix with components 1/(i*beta2_n).
	for(int i=0; i<twoNp1_; ++i) {
		gsl_matrix_complex_set(iBeta2Diag_, i, i, gsl_complex_mul_imag(beta2_[i], 1.0));
	}
	
	// WRONG: The T matrix is = 0.5 * (u_ - iBeta2Diag_ uprime_).  Compute and store back in T_.
	// RIGHT: The actual matrix should be (iBeta2Diag_ u - uprime), where iBeta2Diag_ is diagonals of i*beta2_n. (Not: 1/(i*beta2_n))
	// copy uprime_ into T_.
	gsl_matrix_complex_memcpy(T_, uprime_);
	// [This function computes C = \alpha A B + \beta C, when A is symmetric.  We will compute it for A = iBeta2Diag_, B = u_, \alpha = 1, \beta = -1, and C = uprime_ = T_.]
	errCode = gsl_blas_zsymm(CblasLeft, CblasUpper, gsl_complex_rect(1,0), iBeta2Diag_, u_, gsl_complex_rect(-1.0,0), T_);
	if(errCode) return PEResult(PEResult::AlgebraError);

	//////////////////////////////////////
	if(printDebugOutput) {

		std::cout << "T matrix:" << std::endl;

		for(int i=0; i<twoNp1_; ++i) {
			for(int j=0; j<twoNp1_; ++j) {
				gsl_complex z = gsl_matrix_complex_get(T_, i, j);
				std::cout << GSL_REAL(z) << "+" << GSL_IMAG(z) << "i, ";
			}
			std::cout << ";" << std::endl;
		}

		std::cout << std::endl << std::endl;
	}
	//////////////////////////////////////
	
	// 5. Set up the input (incidence) basis vector Vincident_
	/////////////////////////////
	
	// fill out the incidence matrix Vincident_: when n=0, Vincident_0 = i(beta2_0 + beta1_0)*exp(-i*beta2_0*a). For all other n, Vincident_ = 0 (still).
	z = gsl_complex_mul_imag(gsl_complex_mul(gsl_complex_add(beta2_[N_], beta1_[N_]), gsl_complex_exp(gsl_complex_mul_imag(beta2_[N_], -a))), 1.0);
	if(printDebugOutput) {
		std::cout << "\nIncident vector component Vincident_0: " << GSL_REAL(z) << "," << GSL_IMAG(z) << std::endl;
	}
	gsl_vector_complex_set(Vincident_, N_, z);
	
	// 6. Now we have a linear system:  T_ A1_ = Vincident_.   Solve for A1_, using standard LU decomposition.
	///////////////////////////////
	
	int s;
	errCode = gsl_linalg_complex_LU_decomp(T_, permutation_, &s);
	if(errCode) return PEResult(PEResult::AlgebraError);
	errCode = gsl_linalg_complex_LU_solve(T_, permutation_, Vincident_, A1_);
	if(errCode) return PEResult(PEResult::AlgebraError);
	
	if(printDebugOutput) {
		std::cout << "\nA1_:" << std::endl;
		for(int i=0; i<twoNp1_; ++i) {
			std::cout << i - N_ << ":\t" << GSL_REAL(gsl_vector_complex_get(A1_, i)) << "\t\t" << GSL_IMAG(gsl_vector_complex_get(A1_, i)) << std::endl;
		}
	}
	
	// 7. Solve for the reflected Rayleigh coefficients B2_.
	////////////////////////////////
	timing_[6] = omp_get_wtime();
	
	// Now we have A1_.  Need to get B2_.  From eqn. II.18, \sum_p { A^(1)_p u_{np}(a) } - A2_0 exp(-i beta2_0 a) \delta_{n,0} = B2_n exp(-i beta2_n a)
	// The sum can be computed by the matrix-vector multiplication (u_ A1_).
	errCode = gsl_blas_zgemv(CblasNoTrans, gsl_complex_rect(1,0), u_, A1_, gsl_complex_rect(0,0), B2_);
	if(errCode) return PEResult(PEResult::AlgebraError);
	// now we need to subtract exp(-i beta2_0 a) when n = 0:
	gsl_complex temp = gsl_complex_exp(gsl_complex_mul_imag(beta2_[N_], -a));
	temp = gsl_complex_sub(gsl_vector_complex_get(B2_, N_), temp);
	gsl_vector_complex_set(B2_, N_, temp);
	// and then divide by exp(i beta2_n a) for all n.
	for(int i=0; i<twoNp1_; ++i) {
		temp = gsl_complex_exp(gsl_complex_mul_imag(beta2_[i], a));
		temp = gsl_complex_div(gsl_vector_complex_get(B2_, i), temp);
		gsl_vector_complex_set(B2_, i, temp);
	}
	
	if(printDebugOutput) {
		std::cout << "\nB2_:" << std::endl;
		for(int i=0; i<twoNp1_; ++i) {
			std::cout << i - N_ << ":\t" << GSL_REAL(gsl_vector_complex_get(B2_, i)) << "\t\t" << GSL_IMAG(gsl_vector_complex_get(B2_, i)) << std::endl;
		}
	}
	
	// 8. Now we have B2_. Compute efficiency and put into result structure.
	////////////////////////////////////////
	
	timing_[7] = omp_get_wtime();
	
	PEResult result(N_);
	result.wavelength = wl_;
	result.incidenceDeg = incidenceDeg;
	
	double effSum = 0;
	for(int i=0; i<twoNp1_; ++i) {
		// int n = i - N_;
		result.eff[i] =  gsl_complex_abs2(gsl_vector_complex_get(B2_, i))*GSL_REAL(beta2_[i])/GSL_REAL(beta2_[N_]);
		// is this a non-propagating order?  Then the real part of beta2_n will be exactly 0, so the efficiency will come out as 0.
		
		effSum += result.eff[i];
	}
	
	timing_[8] = omp_get_wtime();
	
	if(measureTiming_) {
		std::cout << "Timing Profile:" << std::endl;
		std::cout << "   Allocate Memory: " << timing_[1] - timing_[0] << std::endl;
		std::cout << "   Setup problem variables: " << timing_[3] - timing_[2] << std::endl;
		std::cout << "   Compute alpha, beta2, and beta1 values: " << timing_[4] - timing_[3] << std::endl;
		std::cout << "   Integrate all trial solutions: " << timing_[5] - timing_[4] << std::endl;
		std::cout << "   Solve linear system for all An: " << timing_[6] - timing_[5] << std::endl;
		std::cout << "   Compute all Bn: " << timing_[7] - timing_[6] << std::endl;
		std::cout << "   Compute and package efficiencies: " << timing_[8] - timing_[7] << std::endl;
		std::cout << "   Total (solver) time: " << timing_[8] - timing_[2] + timing_[1] - timing_[0] << std::endl << std::endl;
	}

	if(printDebugOutput) {
		std::cout << "Sum of efficiencies (should be 1): " << effSum << std::endl;
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
	// For the simple supported profiles, we will always have an x1 value where the line at y enters the grating material, and an x2 value where the line exits the grating material.  Both x1 and x2 > 0, < d.
	double x1, x2;
	// grating period:
	double d = g_.period();

	// Determine x1 and x2, based on the geometry of the profile.
	switch(g_.profile()) {

	case PEGrating::RectangularProfile:
		// For rectangular profiles, geo(0) is the depth, geo(1) is the valley width.
		if(y > g_.geo(0) + 1e-16) {
			std::cout << "Grating Expansion: Error: Height " << y-0.04 << " is above the grating maximum " << g_.geo(0)-0.04 << std::endl;
			return PEResult::InvalidGratingFailure;	// above the grating.
		}
		if(g_.geo(1) > d) {
			std::cout << "Grating Expansion: Error: Valleys " << g_.geo(1) << " are wider than period " << d << std::endl;
			return PEResult::InvalidGratingFailure;	// valleys are wider than the period... Not possible.
		}
		
		// simple, because the groove walls don't depend on y:
		x1 = g_.geo(1);
		x2 = d;
		break;
		
	case PEGrating::BlazedProfile: {
		// For blazed profiles, geo(0) is the blaze angle, and geo(1) is the anti-blaze angle, both in degrees.
		double blaze = g_.geo(0)*M_PI/180.0;
		double antiBlaze = g_.geo(1)*M_PI/180.0;
		
		if(blaze <= 0 || blaze >= M_PI_2 || antiBlaze <= 0 || antiBlaze > M_PI_2)
			return PEResult::InvalidGratingFailure;
		
		x1 = y / tan(blaze);
		x2 = d - y / tan(antiBlaze);
		
		if(x2 < x1 - 1e-16) {
			std::cout << "Grating Expansion: Error: Above the grating: x2 is less than x1 by " << x1-x2 << std::endl;
			return PEResult::InvalidGratingFailure;	// above the grating.
		}
		if(x2 < x1) {
			x2 = x1;
		}
		break;
	}
		
	case PEGrating::SinusoidalProfile: {
		// for sine profile, the surface is described by y = g(x) = depth/2 * ( 1 - cos(2pi x/d) ).
		// geo(0) is the depth.
		
		double depth = g_.geo(0);
		if(y > depth)
			return PEResult::InvalidGratingFailure;	// above the grating.
		
		x1 = depth/2/M_PI * acos(1 - 2*y/depth);	// inverse of grating profile formula above; will always return x in first half of period.
		x2 = d - x1;	// due to symmetry.
		break;
	}

	case PEGrating::TrapezoidalProfile:
	case PEGrating::InvalidProfile:
	case PEGrating::CustomProfile:
	default:
		return PEResult::InvalidGratingFailure;
	}
	
	// wave number in free space: k_2 = v_2 * w / c.  v_2 = 1 in empty space, so k_2 = 2pi / wl.
	gsl_complex k_2 = gsl_complex_rect(2 * M_PI / wl_, 0);
	// wave number in the grating: k_1 = v_1 * w / c = v_1 * 2pi / wl = v_1 * k_2.
	gsl_complex k_1 = gsl_complex_mul(v_1_, k_2);
	
	// square them to get k^2_2 and k^2_1:
	gsl_complex k2_2 = gsl_complex_mul(k_2, k_2);
	gsl_complex k2_1 = gsl_complex_mul(k_1, k_1);
	
	
	// sigma1 (s1) and sigma2 (s2) are defined as: s_p = k^2_(p+1) - k^2_p... Except at p = P (max): s_P = k^2_1 - k^2_P.  We only have two crossings, so P = 2.  p = [1,2].
	// k^2_p is the refractive index value on the left side of x_p crossing.
	gsl_complex s1 = gsl_complex_sub(k2_1, k2_2);
	gsl_complex s2 = gsl_complex_sub(k2_2, k2_1);
	
	// for shorter calculations, define K as the grating number, 2*pi/d.  We just calculate it once.
	double K = 2*M_PI/d;
	
	// loop over all values of n, from -N to N.  See paper notes for formula for k^2_n [Page 5,6]
	for(int i=0; i<twoNp1_; ++i) {
		int n = i - N_;
		
		if(n == 0) {
			gsl_complex kk = gsl_complex_mul_real(k2_2, d);
			kk = gsl_complex_sub(kk, gsl_complex_mul_real(s1, x1));
			kk = gsl_complex_sub(kk, gsl_complex_mul_real(s2, x2));
			
			kk = gsl_complex_mul_real(kk, 1.0/d);
			k2[i] = kk;
		}
		else {
			double t1 = n*K*x1, t2 = n*K*x2;
			gsl_complex kk = gsl_complex_mul(s1, gsl_complex_rect(sin(t1), cos(t1)));
			kk = gsl_complex_add(kk, gsl_complex_mul(s2, gsl_complex_rect(sin(t2), cos(t2))));
			
			kk = gsl_complex_mul_real(kk, -1.0/2/M_PI/n);
			k2[i] = kk;
		}
	}
	
	return PEResult::Success;
}

PEResult::Code PESolver::integrateTrialSolutionAlongY(gsl_vector_complex* u, gsl_vector_complex* uprime) {
	// define ode solving system, with our function to evaluate dw/dy, no Jacobian, and 8*N_+4 components.
	gsl_odeiv2_system odeSys = {odeFunctionCB, odeJacobianCB, 8*N_+4, this};
	
	// initial starting step in y: choose grating height / 200.
	double gratingHeight = g_.height();
	double hStart = gratingHeight / 200.0;
	
	// setup driver
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_standard_new (&odeSys, gsl_odeiv2_step_msadams, hStart, 1e-8, 1e-8, 0.5, 0.5);

	
	// fill starting conditions from u, uprime
	double* w = new double[8*N_+4];
	int fourNp2 = 4*N_+2;
	for(int i=0; i<twoNp1_; ++i) {
		gsl_complex u_n = gsl_vector_complex_get(u, i);
		w[2*i] = GSL_REAL(u_n);
		w[2*i + 1] = GSL_IMAG(u_n);
		
		gsl_complex uprime_n = gsl_vector_complex_get(uprime, i);
		w[2*i + fourNp2] = GSL_REAL(uprime_n);
		w[2*i + fourNp2 + 1] = GSL_IMAG(uprime_n);
	}
	
	// run it: integrate from y = 0 to gratingHeight.
	double y = 0;
	int status = gsl_odeiv2_driver_apply (d, &y, gratingHeight, w);
	if (status != GSL_SUCCESS) {
		std::cout << "ODE: Integration failure: Code: " << status << std::endl;
		return PEResult::ConvergenceFailure;
	}
	
	// copy results back into u, uprime
	for(int i=0; i<twoNp1_; ++i) {
		gsl_vector_complex_set(u, i, gsl_complex_rect(w[2*i], w[2*i + 1]));
		gsl_vector_complex_set(uprime, i, gsl_complex_rect(w[2*i + fourNp2], w[2*i + fourNp2 + 1]));
	}
	
	gsl_odeiv2_driver_free(d);
	delete [] w;
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
	
	int fourNp2 = 4*N_+2;	// fourNp2 divides the top and bottom of the arrays w, f.  Top of w is u; bottom of w is u' = v.   Top of f is u' = v;  bottom of f is v' = u''.
	int eightNp4 = 8*N_+4;
	// size is (2N+1)x2x2, ie: 8N+4.
	for(int i=0; i<fourNp2; ++i) {
		// working on computing u'_n = v. [top of f array]. Just copy from v = u' = [bottom half of w array]
		dwdy[i] = w[i + fourNp2];
	}

	for(int i=fourNp2; i<eightNp4; i+=2) {
		// working on computing v'_n = u''_n [bottom of f array].  We'll handle real and imaginary components at once.
		int n = (i - fourNp2)/2 - N_;
		gsl_complex upp_n = gsl_complex_rect(0,0);	// initialize sum to 0.
		// loop over m.  Retrieving u_n values from top of array.
		for(int j=0; j<fourNp2; j+=2) { // also handling real and imaginary components at once.
			int m = j/2 - N_;

			gsl_complex u_m = gsl_complex_rect(w[j], w[j+1]);

			gsl_complex M_nm = gsl_complex_rect(0,0);
			if(n-m >= -N_ && n-m <= N_)
				M_nm = gsl_complex_mul_real(localK2[n-m + N_], -1.0);	// -k^2_{n-m}
			if(n == m)
				M_nm = gsl_complex_add_real(M_nm, alpha_[n + N_]);

			upp_n = gsl_complex_add(upp_n, gsl_complex_mul(M_nm, u_m));
		}

		dwdy[i] = GSL_REAL(upp_n);
		dwdy[i+1] = GSL_IMAG(upp_n);
	}

	// Debugging output:
	//////////////////////
	//	double sum;
	//	for(int i=0; i<4*N_+2; i+=2) {
	//		gsl_complex u_n = gsl_complex_rect(w[i], w[i+1]);
	//		sum += GSL_REAL(gsl_complex_mul(u_n, gsl_complex_conjugate(u_n)));
	//	}
	//	std::cout << "y: " << y << " Power: " << sum << std::endl;
	
	return GSL_SUCCESS;
}


int PESolver::odeJacobian(double y, const double w[], double *dfdw, double dfdy[])
{
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

	int fourNp2 = 4*N_+2;	// fourNp2 divides the top and bottom of the arrays w, f.  Top of w is u; bottom of w is u' = v.   Top of f is u' = v;  bottom of f is v' = u''.
	int eightNp4 = 8*N_+4;	// size is (2N+1)x2x2, ie: 8N+4.

	// clear the jacobian
	memset(dfdw, 0, eightNp4*eightNp4*sizeof(double));
	memset(dfdy, 0, eightNp4*sizeof(double));

	// go through top rows of jac [i=0,fourNp2]. Set ident. matrix in upper right-hand block.
	for(int i=0; i<fourNp2; ++i) {
		dfdw[i*eightNp4+fourNp2+i] = 1.0;		// set at index dfdw(i, fourNp2+i).
	}

	// go throgh lower rows of jac [i=fourNp2, eightNp4]. In lower left-hand block, set 4x4 submatrices at once, so go by i+=2.
	for(int i=fourNp2; i<eightNp4; i+=2) {
		int n = (i - fourNp2)/2 - N_;	// n index [row] from -N to N.

		// loop over m [col].
		for(int j=0; j<fourNp2; j+=2) { // also handling real and imaginary components at once.
			int m = j/2 - N_;	// m index [col] from -N to N.

			// get M_{nm}: -(k^2)_{n-m}(y) + alpha^2_n \delta_{nm}
			gsl_complex M_nm = gsl_complex_rect(0,0);
			if(n-m >= -N_ && n-m <= N_)
				M_nm = gsl_complex_mul_real(localK2[n-m + N_], -1.0);	// k2_ ranges from -N to N.
			if(n == m)
				M_nm = gsl_complex_add_real(M_nm, alpha_[n + N_]);

			// set 4x4 matrix here: [M_re, -M_im; M_im, M_re] at (i,j), (i, j+1); (i+1, j), (i+1, j+1)
			dfdw[i*eightNp4 + j] = GSL_REAL(M_nm);
			dfdw[i*eightNp4 + j + 1] = -GSL_IMAG(M_nm);
			dfdw[(i+1)*eightNp4 + j] = GSL_IMAG(M_nm);
			dfdw[(i+1)*eightNp4 + j + 1] = GSL_REAL(M_nm);
		}
	}

	/// \todo IMPORTANT! Leaving dfdy = 0 for now. This is only true in case of rectangular grating...

	return GSL_SUCCESS;
}

gsl_complex* PESolver::k2ForCurrentThread() {
	/// Return based on OpenMP current thread.
	return k2_[omp_get_thread_num()];
}

double PESolver::conditionNumber(const gsl_matrix_complex *A) const
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


