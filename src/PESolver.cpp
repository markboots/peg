#include "PESolver.h"

#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>


/// Speed of light in um/s.
#define M_c 2.99792458e14

PESolver::PESolver(const PEGrating& grating, const PEMathOptions& mo)
	: g_(grating)
{
	// set math options.
	N_ = mo.N;
	twoNp1_ = 2*N_ + 1;
	niy_ = mo.niy;
	
	// allocate matrices and vectors
	u_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	uprime_ = gsl_matrix_complex_alloc(twoNp1_, twoNp1_);
	
	alpha_ = new double[twoNp1_];
	beta2_ = new gsl_complex[twoNp1_];
	beta1_ = new gsl_complex[twoNp1_];
	k2_ = new gsl_complex[4*N_+1];
	
}

/// \todo Imp.
PESolver::~PESolver() {
	gsl_matrix_complex_free(u_);
	gsl_matrix_complex_free(uprime_);
	
	delete [] alpha_;
	delete [] beta2_;
	delete [] beta1_;
	delete [] k2_;
}

/// \todo Imp.
PEResult PESolver::getEff(double incidenceDeg, double wl) {
	
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
	
	// compute all alpha_n and beta_n.  (This could be a parallel loop)
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
	
	// set up trial solutions at y=0.
	// u_ should be the identify matrix
	gsl_matrix_complex_set_identity(u_);
	// zero uprime_, before we set the components we need to.
	gsl_matrix_complex_set_zero(uprime_);
	// Loop over all trial solutions p:  [row: n.  col: p]
	for(int i=0; i<twoNp1_; i++) {
		int p = i - N_;
		
		// insert initial values at y=0.  u'(0) = 1* \delta_{np}.  (Already set: identity matrix)
		// u'(0) = -i * beta^(1)_n \delta_{np}
		gsl_matrix_complex_set(uprime_, i, i, gsl_complex_mul_imag(beta1_[i], -1));
		
		// create vector views: the column vectors for this trial solution
		gsl_vector_complex_view u = gsl_matrix_complex_column(u_, i);
		gsl_vector_complex_view uprime = gsl_matrix_complex_column(uprime_, i);
		
		// need to integrate u and uprime along y, using d^2 u/dy^2 = M(y) u
		PEResult::Code err = integrateTrialSolutionAlongY(&u.vector, &uprime.vector);
		if(err != PEResult::Success)
			return PEResult(err);
		
		// Now we have u(a) and u'(a) in (u, uprime) and in the columns of the actual (u_, uprime_) matrices.
	}
	
	// construct Ax = b matrices
	// solve for An
	// construct Ax = b matrices
	// solve for Bn
	
	return PEResult();	
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
		if(y > g_.geo(0))
			return PEResult::InvalidGratingFailure;	// above the grating.
		if(g_.geo(1) > d)
			return PEResult::InvalidGratingFailure;	// valleys are wider than the period... Not possible.
		
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
		
		if(x2 < x1)
			return PEResult::InvalidGratingFailure;	// above the grating.
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
	
	// loop over all values of n, from -2N to 2N.  See paper notes for formula for k^2_n [Page 5,6]
	for(int i=0; i<(4*N_+1); ++i) {
		int n = i - 2*N_;
		
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
	gsl_odeiv2_system odeSys = {odeFunctionCB, 0, 8*N_+4, this};
	
	// initial starting step in y: choose grating height / 200.
	double gratingHeight = g_.height();
	double hStart = gratingHeight / 200.0;
	
	// setup driver
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&odeSys, gsl_odeiv2_step_rkf45, hStart, 1e-8, 0.0);	/// \todo Use standard (y, yp) instead of y error control?  What tolerances?
	
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
		return PEResult::ConvergenceFailure;
	}
	
	// copy results back into u, uprime
	for(int i=0; i<twoNp1_; ++i) {
		gsl_vector_complex_set(u, i, gsl_complex_rect(w[2*i], w[2*i + 1]));
		gsl_vector_complex_set(uprime, i, gsl_complex_rect(w[2*i + fourNp2], w[2*i + fourNp2 + 1]));
	}
	
	delete [] w;	
	return PEResult::Success;
}

int PESolver::odeFunction(double y, const double w[], double f[]) {
	
	// w contains the last values of u_n{re, im} and u'_n{re, im}, in that order.
	// need to compute f = dw/dy.
	
	// get k2_n at this y value.
	if(computeGratingExpansion(y, k2_) != PEResult::Success) {
		return GSL_EBADFUNC;	// can't calculate here. Invalid profile? y above the profile height?
	}
	
	// size is (2N+1)x2x2, ie: 8N+4.
	for(int i=0, eightNp4=8*N_+4; i<eightNp4; ++i) {
		int fourNp2 = 4*N_+2;	// fourNp2 divides the top and bottom of the arrays w, f.  Top of w is u; bottom of w is u' = v.   Top of f is u' = v;  bottom of f is v' = u''.
		
		if(i < fourNp2) {
			// working on computing u'_n = v. [top of f array]. Just copy from v = u' = [bottom half of w array]
			f[i] = w[i + fourNp2];
		}
		
		else {
			// working on computing v'_n = u''_n [bottom of f array].  We'll handle real and imaginary components at once, so only do this if i is even.
			if(i%2 == 0) {
				int n = (i - fourNp2)/2 - N_;
				gsl_complex upp_n = gsl_complex_rect(0,0);	// initialize sum to 0.
				// loop over m.  Retrieving u_n values from top of array.
				for(int j=0; j<fourNp2; ++j) {
					if(j%2 == 0) {	// also handling real and imaginary components at once.
						int m = j/2 - N_;
					
						gsl_complex u_m = gsl_complex_rect(w[j], w[j+1]);
						
						gsl_complex minus_k2 = gsl_complex_mul_real(k2_[n-m + 2*N_], -1.0);	// k2_ ranges from -2N to 2N.
						if(n == m)
							minus_k2 = gsl_complex_add_real(minus_k2, alpha_[n + N_]);
						
						upp_n = gsl_complex_add(upp_n, gsl_complex_mul(minus_k2, u_m));
					}
				}
				
				f[i] = GSL_REAL(upp_n);
				f[i+1] = GSL_IMAG(upp_n);
			}
		}
	}
	
	return GSL_SUCCESS;
}
