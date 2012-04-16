#ifndef PEG_H
#define PEG_H

#define HAVE_INLINE
#include <gsl/gsl_complex.h>

#include <string>
#include <vector>
#include <ostream>

// h*c (Planck constant * speed of light), in eV * um.  Used in conversion from eV to um.
#define M_HC 1.23984172

// Common definitions for the PEG parallel grating efficiency library

/// This type is returned by a single grating efficiency calculation. It contains a status code/error code to indicate the result of the calculation, a vector of inside order efficiencies, and a vector of outside order efficiencies. The first element in the output vectors is the 0 order, and is duplicated over both.
class PEResult {
public:
	enum Code { Success, InvalidGratingFailure, ConvergenceFailure, InsufficientCoefficientsFailure, AlgebraError, OtherFailure };
	
	/// Constructs an empty result with the given \c statusCode
	PEResult(Code statusCode = OtherFailure) { status = statusCode; }
	/// Constructs a successful result, where insideEff and outsideEff have size \c N+1.
	PEResult(int N) : insideEff(N+1), outsideEff(N+1) {
		status = Success;
	}
	
	Code status;
	double wavelength;
	double incidenceDeg;
	std::vector<double> insideEff;
	std::vector<double> outsideEff;
	
	/// Prints the output efficiencies in a table, to standard output.
	friend std::ostream& operator<<(std::ostream& os, const PEResult& result);
};

std::ostream& operator<<(std::ostream& os, const PEResult& result);

/// Represents the numerical options to be used for a single grating calculation.
class PEMathOptions {
public:
	/// Fourier truncation index N. Should be a positive number.
	int N;
	/// Number of integration steps along y
	int niy;
	
	/// Constructor:
	PEMathOptions(int FourierN = 15, int numberOfIntegrationSteps = 401) {
		N = FourierN;
		niy = numberOfIntegrationSteps;
	}
};


/// Represents the parameters for a grating
class PEGrating {
public:
	
	/// Specifies the grating profile, if one of the standard profiles, or CustomProfile
	enum Profile { InvalidProfile, RectangularProfile, BlazedProfile, SinusoidalProfile, TrapezoidalProfile, CustomProfile };
	
	/// Default constructor; does not provide a valid grating. Use the constructors in the subclasses for a valid grating.
	PEGrating() {
		profile_ = InvalidProfile;
		period_ = 1.0;
		material_ = "Au";
	}
	
	virtual ~PEGrating() {}
	
	// Accessor Functions
	////////////////////////
	
	/// Returns the profile of this grating: CustomProfile, RectangularProfile, BlazedProfile, SinusoidalProfile, TrapezoidalProfile
	Profile profile() const { return profile_; }
	/// Returns the grating periodicity, in um
	double period() const { return period_; }
	/// Returns profile-dependent geometry parameters
	double geo(int parameterIndex) const { return geo_[parameterIndex]; }
	/// Returns the height from the bottom of the grooves to the top of the grooves.  Profile-dependent, but the base class implementation handles the default rectangular, blazed, sinusoidal, and trapezoidal.
	virtual double height() const;
	
	/// Returns the grating material code (ex: "Ni", "Au", etc.)
	std::string material() const { return material_; }
	
	
	/// Returns the complex refractive index at a given wavelength \c wl in um. \todo Imp. Current is Pt at 410 eV.
	gsl_complex refractiveIndex(double wl) const;
	
	
	/// Calculates the grating efficiency at a given incidence angle \c incidenceDeg (degrees) and wavelength \c wl (um).
	PEResult getEff(double incidenceDeg, double wl, const PEMathOptions& mo = PEMathOptions(), bool printDebugOutput = false) const;
	
	
protected:
	/// Grating profile type
	Profile profile_;
	/// Grating periodicity, in um
	double period_;
	/// General geometry parameters. Interpretation depends on profile.
	double geo_[8];
	/// Grating material
	std::string material_;
	
};

/// Rectangular grating subclass
class PERectangularGrating : public PEGrating {
public:
	/// Constructs a grating with a rectangular profile. The required geometry parameters are the groove \c height in um, and the \c valleyWidth in um.  The \c valleyWidth is the width of the low part of the groove, and must obviously be less than the period.
	PERectangularGrating(double period = 1.0, double height = 0.05, double valleyWidth = 0.5, const std::string& material = "Au") {
		profile_ = RectangularProfile;
		period_ = period;
		geo_[0] = height;
		geo_[1] = valleyWidth;
		material_ = material;
	}
};

/// Blazed grating subclass
class PEBlazedGrating : public PEGrating {
public:
	/// Constructs a grating with the blazed profile. The required geometry parameters are the blaze angle \c blazeAngleDeg, in deg., and the anti-blaze angle \c antiBlazeAngleDeg.
	PEBlazedGrating(double period = 1.0, double blazeAngleDeg = 2.0, double antiBlazeAngleDeg = 30, const std::string& material = "Au") {
		profile_ = BlazedProfile;
		period_ = period;
		geo_[0] = blazeAngleDeg;
		geo_[1] = antiBlazeAngleDeg;
		material_ = material;
	}
};

/// Sinusoidal grating subclass
class PESinusoidalGrating : public PEGrating {
public:
	/// Constructs a grating with a perfect sinusoidal profile. The only required geometry parameter is the groove \c height, in um.
	PESinusoidalGrating(double period = 1.0, double height = 0.05, const std::string& material = "Au") {
		profile_ = BlazedProfile;
		period_ = period;
		geo_[0] = height;
		material_ = material;
	}
};

/// Trapezoidal grating subclass
class PETrapezoidalGrating : public PEGrating {
public:
	/// Constructs a grating with a trapezoidal profile. The required geometry parameters are the \c height, in um, the \c valleyWidth, in um, the blaze angle \c blazeAngleDeg, in deg., and the anti-blaze angle \c antiBlazeAngleDeg.
	PETrapezoidalGrating(double period = 1.0, double height = 0.05, double valleyWidth = 0.5, double blazeAngleDeg = 30.0, double antiBlazeAngleDeg = 30.0, const std::string& material = "Au") {
		profile_ = BlazedProfile;
		period_ = period;
		geo_[0] = height;
		geo_[1] = valleyWidth;
		geo_[2] = blazeAngleDeg;
		geo_[3] = antiBlazeAngleDeg;
		material_ = material;
	}
};

#endif
