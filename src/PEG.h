#ifndef PEG_H
#define PEG_H

#define HAVE_INLINE
#include <gsl/gsl_complex.h>

#include <string>
#include <vector>

// Common definitions for the PEG parallel grating efficiency library


/// This type is returned by a single grating efficiency calculation. It contains a status code/error code to indicate the result of the calculation, a vector of inside order efficiencies, and a vector of outside order efficiencies. The first element in the output vectors is the 0 order, and is duplicated over both.
class PEResult {
public:
	enum Code { Success, InvalidGratingFailure, ConvergenceFailure, InsufficientCoefficientsFailure };
	
	Code status;
	std::vector<double> insideEff;
	std::vector<double> outsideEff;
};

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
	
	/// Default constructor; does not provide a valid grating. Use the constructors in the subclasses.
	PEGrating() {
		profile_ = InvalidProfile;
		period_ = 1.0;
		material_ = "Au";
	}
	
	// Accessor Functions
	////////////////////////
	
	/// Returns the profile of this grating: CustomProfile, RectangularProfile, BlazedProfile, SinusoidalProfile, TrapezoidalProfile
	Profile profile() const { return profile_; }
	/// Returns the grating periodicity, in um
	double period() const { return period_; }
	/// Returns profile-dependent geometry parameters
	double geo(int parameterIndex) const { return geo_[parameterIndex]; }
	
	/// Returns the grating material code (ex: "Ni", "Au", etc.)
	std::string material() const { return material_; }
	
	
	/// Returns the complex refractive index at a given wavelength \c wl in um. \todo Imp. Current is Pt at 410 eV.
	gsl_complex refractiveIndex(double wl) const;
	
	
	/// Calculates the grating efficiency at a given incidence angle \c incidenceDeg (degrees) and wavelength \c wl (um).
	PEResult getEff(double incidenceDeg, double wl, const PEMathOptions& mo = PEMathOptions()) const;
	
	
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

class PEBlazedGrating : public PEGrating {
public:
	/// Constructs a grating with the blazed profile. The required geometry parameters are \c blazeAngleDeg, the blaze angle in deg., and the anti-blaze angle \c antiBlazeAngle
	PEBlazedGrating(double period = 1.0, double blazeAngleDeg = 2.0, double antiBlazeAngleDeg = 30, const std::string& material = "Au") {
		profile_ = BlazedProfile;
		period_ = period;
		geo_[0] = blazeAngleDeg;
		geo_[1] = antiBlazeAngleDeg;
		material_ = material;
	}
		
};

#endif
