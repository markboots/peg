/*
Copyright 2012 Mark Boots (mark.boots@usask.ca).

This file is part of the Parallel Efficiency of Gratings project ("PEG").

PEG is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PEG is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PEG.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PEG_H
#define PEG_H

#define HAVE_INLINE
#include <gsl/gsl_complex.h>

#include <string>
#include <vector>
#include <ostream>
#include <math.h>

/// h*c (Planck constant * speed of light), in eV * um.  Used in conversion from eV to um.
#define M_HC 1.23984172

/// Path to materials database (folder)
#define PEG_MATERIALS_DB_PATH "materialDatabase"

/// Maximum supported number of interfaces crossed in a horizontal slice through the grating; determines the size of arrays in computeK2StepsAtY().
#define PEG_MAX_PROFILE_CROSSINGS 60

// Common definitions for the PEG parallel grating efficiency library

/// This type is returned by a single grating efficiency calculation. It contains a status code/error code to indicate the result of the calculation, a vector of inside order efficiencies, and a vector of outside order efficiencies. The first element in the output vectors is the 0 order, and is duplicated over both.
class PEResult {
public:
	enum Code { Success, InvalidGratingFailure, ConvergenceFailure, InsufficientCoefficientsFailure, AlgebraFailure, MissingRefractiveDataFailure, OtherFailure, InactiveCalculation };
	
	/// Constructs an empty result with the given \c statusCode
	PEResult(Code statusCode = OtherFailure) { status = statusCode; }
	/// Constructs a successful result, where the eff array has size \c 2*N+1.
	PEResult(int N) : eff(2*N+1) {
		status = Success;
	}
	
	/// Result of the calculation: Success, InvalidGratingFailure, ConvergenceFailure, InsufficientCoefficientsFailure, AlgebraError, or OtherFailure
	Code status;
	/// Wavelength for this calculation
	double wavelength;
	/// Incidence angle for this calculation
	double incidenceDeg;
	/// Array of efficiencies, going from -N order up to N.  Size is 2N+1; the 0-order efficiency can be found at eff[N].
	std::vector<double> eff;
	
	/// This packs the current result into a plain double \c array, for easy communication using standard MPI types.  The array must have pre-allocated room for 2N+1 + 4 elements.
	void toDoubleArray(double* array) const;
	/// This unpacks the result from a plain double \c array that was filled by toDoubleArray().
	void fromDoubleArray(const double* array);
	
	/// Prints the output efficiencies in a table, to standard output.
	friend std::ostream& operator<<(std::ostream& os, const PEResult& result);
};

std::ostream& operator<<(std::ostream& os, const PEResult& result);

/// Represents the numerical options to be used for a single grating calculation.
class PEMathOptions {
public:
	/// Fourier truncation index N. Should be a positive number; Fourier components from [-N, N].
	int N;
	/// Tolerance required (eps) in the numerical integration process at each step
	double integrationTolerance;
	
	/// Constructor
	PEMathOptions(int FourierN = 15, double IntegrationTolerance = 1e-5) {
		N = FourierN;
		integrationTolerance = IntegrationTolerance;
	}
};


/// Represents the parameters and geometry of a grating. Subclassed as required for different profiles.
class PEGrating {
public:
	
	/// Specifies the grating profile, if one of the standard profiles, or CustomProfile
	enum Profile { InvalidProfile, RectangularProfile, BlazedProfile, SinusoidalProfile, TrapezoidalProfile, CustomProfile };
	
	/// Default constructor; does not provide a valid grating. Use the constructors in the subclasses for a valid grating.
	PEGrating() {
		profile_ = InvalidProfile;
		period_ = 1.0;
		substrateMaterial_ = "SiO2";
		coatingMaterial_ = "Au";
		coatingThickness_ = 0.;
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
	/// Returns the height from the substrate surface (y=0) to the highest feature. If the grating has a coating, this should include the coating.
	/*! This is shape-dependent, but the base class implementation handles the default rectangular, blazed, sinusoidal, and trapezoidal profiles. Re-implement for custom profiles.*/
	virtual double totalHeight() const { return profileHeight() + coatingThickness_; }

	/// Returns the substrate material code (ex: "Ni", "Au", etc.)
	std::string substrateMaterial() const { return substrateMaterial_; }
	/// Returns the coating material code (ex: "Ni", "Au", etc.)
	std::string coatingMaterial() const { return coatingMaterial_; }
	/// Returns the coating thickness, or 0 if there is o coating.
	double coatingThickness() const { return coatingThickness_; }
	
	/// Returns the complex refractive index of the substrate at a given wavelength \c wl in um.  Returns gsl_complex_rect(0,0) if the substrate material's database was not found.
	gsl_complex substrateRefractiveIndex(double wl) const { return refractiveIndex(wl, substrateMaterial_); }
	/// Returns the complex refractive index of the coating at a given wavelength \c wl in um.  Returns gsl_complex_rect(0,0) if the coating material's database was not found.
	gsl_complex coatingRefractiveIndex(double wl) const { return refractiveIndex(wl, coatingMaterial_); }

	/// Looks up the complex refractive index of \c material at a given wavelength \c wl in um.  Returns gsl_complex_rect(0,0) if the  material's database was not found.
	static gsl_complex refractiveIndex(double wl, const std::string& material);

	
	
	/// Calculates the grating efficiency at a given incidence angle \c incidenceDeg (degrees) and wavelength \c wl (um). \c numThreads is the number of threads to use for fine parallelization; ideally it should be <= the number of processor cores on your computer / on a single cluster node.
	PEResult getEff(double incidenceDeg, double wl, const PEMathOptions& mo = PEMathOptions(), bool printDebugOutput = false, int numThreads = 1, bool measureTiming = false) const;


	// Detailed geometry
	/////////////////////////

	/// Specific to each grating shape, this computes the values of the multistep function of the impedence k^2 at a height \c y, where \c y can be from [0, totalHeight()].  (This is the step function created by slicing horizontally through the structure.) Going from x = 0 to x = period(), every intersection with a new medium should result in an x-coordinate stored in \c stepsX, and the k^2 value at (x - eps, y) stored in stepsK2. This function should write values into \c stepsX, \c stepsY, and return the number of steps.  Return -1 to indicate an error (invalid geometry).
	/*!
	- For simple profiles with no coating, there are usually 2 steps: entering the material, and leaving the material.
	- Coated gratings may have 4 steps at certain y values, if the are <i>interpenetrating</i>, i.e., the coating is thinner than the uncoated profileHeight().
	- Profiles with multiple local maxima could have more steps at some \c y values, from entering and exiting the material multiple times.
	- Within a homogeneous layer, there is only "one" step, with \c stepsX[0] = 0, \c stepsK2[0] = <the layer impedance>.

	This function is called at each integration step, so we avoid dynamic memory for performance here. Note that \c stepsX and \c stepsK2 only have storage for a maximum of PEG_MAX_PROFILE_CROSSINGS.

	The base class implementation handles simple profile shapes (those with a single local maximum), with or without an interpenetrating or thick coating. For this to work, the subclass must implement xIntersection1() and xIntersection2().
	*/
	virtual int computeK2StepsAtY(double y, gsl_complex k2_vaccuum, gsl_complex k2_substrate, gsl_complex k2_coating, double* stepsX, gsl_complex* stepsK2) const;


	// Computational Geometry. The following geometry functions describe the basic, bare profile, assuming there is no coating.
	////////////////

	/// Returns the height of the bare bump, ignoring any coating that might exist.
	/*! This is shape-dependent, but the base class implementation handles the default rectangular, blazed, sinusoidal, and trapezoidal profiles. Re-implement for custom profiles.*/
	virtual double profileHeight() const { return 0.; }
	/// Returns the x-coordinate of the first intersection with the bump [i.e., entering the material], at height \c y (assuming no coating).
	/*! This is used by computeK2StepsAtY() for simple bump shapes with a single maximum; if re-implementing computeK2StepsAtY(), you can omit this.

	\c y will range from 0 to profileHeight().

	Base class returns negative number to indicate geometry failure; must re-implement.
*/
	virtual double xIntersection1(double y) const { (void)y; return -1; }
	/// Returns the x-coordinate of the second intersection with the bump [i.e., leaving the material], at height \c y (assuming no coating).
	/*! This is used by computeK2StepsAtY() for simple bump shapes with a single maximum; if re-implementing computeK2StepsAtY(), you can omit this.

	\c y will range from 0 to profileHeight().

	Base class returns negative number to indicate geometry failure; must re-implement.
*/
	virtual double xIntersection2(double y) const { (void)y; return -1; }

	////////////////////////////
	
	
protected:
	/// Grating profile type
	Profile profile_;
	/// Grating periodicity, in um
	double period_;
	/// General geometry parameters. Interpretation depends on profile.
	double geo_[8];
	/// Substrate material
	std::string substrateMaterial_;
	/// Coating material
	std::string coatingMaterial_;
	/// The thickness of the coating (um), or 0 for no coating.
	double coatingThickness_;

	// Helper functions:
	/////////////////////////////

	/// Implements computeK2StepsAtY() for no coating.
	int computeK2StepsAtY_noCoating(double y, gsl_complex k2_vaccuum, gsl_complex k2_substrate, gsl_complex k2_coating, double* stepsX, gsl_complex* stepsK2) const;
	/// Implements computeK2StepsAtY() for an interpenetrating coating.
	int computeK2StepsAtY_interpenetratingCoating(double y, gsl_complex k2_vaccuum, gsl_complex k2_substrate, gsl_complex k2_coating, double* stepsX, gsl_complex* stepsK2) const;
	/// Implements computeK2StepsAtY() for a thick (non-interpenetrating) coating.
	int computeK2StepsAtY_thickCoating(double y, gsl_complex k2_vaccuum, gsl_complex k2_substrate, gsl_complex k2_coating, double* stepsX, gsl_complex* stepsK2) const;
};

/// Rectangular grating subclass
class PERectangularGrating : public PEGrating {
public:
	/// Constructs a grating with a rectangular profile. The required geometry parameters are the groove \c height in um, and the \c valleyWidth in um.  The \c valleyWidth is the width of the low part of the groove, and must obviously be less than the period.
	PERectangularGrating(double period = 1.0, double height = 0.05, double valleyWidth = 0.5, const std::string& material = "Au", const std::string& coating = "Au", double coatingThickness = 0) {
		profile_ = RectangularProfile;
		period_ = period;
		geo_[0] = height;
		geo_[1] = valleyWidth;
		substrateMaterial_ = material;
		coatingMaterial_ = coating;
		coatingThickness_ = coatingThickness;
	}

	/// geo_[0] is the height, directly.
	virtual double profileHeight() const { return geo(0); }

	/// Returns the x-coordinate of the first intersection with the surface at \c y. Simple, because the grating doesn't change with height.
	virtual double xIntersection1(double y) const { (void)y; return geo(1); }
	/// Returns the x-coordinate of the second intersection with the surface at \c y. Simple, because the grating doesn't change with height.
	virtual double xIntersection2(double y) const { (void)y; return period(); }
};

/// Blazed grating subclass
class PEBlazedGrating : public PEGrating {
public:
	/// Constructs a grating with the blazed profile. The required geometry parameters are the blaze angle \c blazeAngleDeg, in deg., and the anti-blaze angle \c antiBlazeAngleDeg.
	PEBlazedGrating(double period = 1.0, double blazeAngleDeg = 2.0, double antiBlazeAngleDeg = 30, const std::string& material = "Au", const std::string& coating = "Au", double coatingThickness = 0) {
		profile_ = BlazedProfile;
		period_ = period;
		geo_[0] = blazeAngleDeg;
		geo_[1] = antiBlazeAngleDeg;
		substrateMaterial_ = material;
		coatingMaterial_ = coating;
		coatingThickness_ = coatingThickness;
	}

	/// geo(0) is blaze, geo(1) is anti-blaze angle, both in degrees.
	virtual double profileHeight() const { return period() / (1/tan(geo(0)*M_PI/180) + 1/tan(geo(1)*M_PI/180)); }

	/// Returns the x-coordinate of the first intersection with the surface at \c y.
	virtual double xIntersection1(double y) const { return y / tan(geo(0)*M_PI/180.0); }
	/// Returns the x-coordinate of the second intersection with the surface at \c y.
	virtual double xIntersection2(double y) const { return period() - y / tan(geo(1)*M_PI/180.0); }
};

/// Sinusoidal grating subclass
class PESinusoidalGrating : public PEGrating {
public:
	/// Constructs a grating with a perfect sinusoidal profile. The only required geometry parameter is the groove \c height, in um.
	PESinusoidalGrating(double period = 1.0, double height = 0.05, const std::string& material = "Au", const std::string& coating = "Au", double coatingThickness = 0) {
		profile_ = BlazedProfile;
		period_ = period;
		geo_[0] = height;
		substrateMaterial_ = material;
		coatingMaterial_ = coating;
		coatingThickness_ = coatingThickness;
	}

	/// geo(0) is the depth, aka height.
	virtual double profileHeight() const { return geo(0); }


	/// Returns the x-coordinate of the first intersection with the surface at \c y.
	/*! For sine profile, the surface is described by y = g(x) = depth/2 * ( 1 - cos(2pi x/d) ).
  // geo(0) is the depth.

	This returns inverse of grating profile formula above; will always return x in first half of period.*/
	virtual double xIntersection1(double y) const { return geo(0)/2/M_PI * acos(1 - 2*y/geo(0)); }
	/// Returns the x-coordinate of the second intersection with the surface at \c y. Due to symmetry, x2 = period() - x1().
	virtual double xIntersection2(double y) const { return period() - xIntersection1(y); }
};

/// Trapezoidal grating subclass
class PETrapezoidalGrating : public PEGrating {
public:
	/// Constructs a grating with a trapezoidal profile. The required geometry parameters are the \c height, in um, the \c valleyWidth, in um, the blaze angle \c blazeAngleDeg, in deg., and the anti-blaze angle \c antiBlazeAngleDeg.
	PETrapezoidalGrating(double period = 1.0, double height = 0.05, double valleyWidth = 0.5, double blazeAngleDeg = 30.0, double antiBlazeAngleDeg = 30.0, const std::string& material = "Au", const std::string& coating = "Au", double coatingThickness = 0) {
		profile_ = BlazedProfile;
		period_ = period;
		geo_[0] = height;
		geo_[1] = valleyWidth;
		geo_[2] = blazeAngleDeg;
		geo_[3] = antiBlazeAngleDeg;
		substrateMaterial_ = material;
		coatingMaterial_ = coating;
		coatingThickness_ = coatingThickness;
	}

	/// geo(0) is the depth, aka height.
	virtual double profileHeight() const { return geo(0); }

	/// Returns the x-coordinate of the first intersection with the surface at \c y. \todo Not yet implemented!
	virtual double xIntersection1(double y) const { (void)y; return -1; }
	/// Returns the x-coordinate of the second intersection with the surface at \c y. \todo Not yet implemented!
	virtual double xIntersection2(double y) const { (void)y; return -1; }
};

#endif
