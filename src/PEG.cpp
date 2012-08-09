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

#include "PEG.h"
#include "PESolver.h"

#include <gsl/gsl_complex_math.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string.h>

// This packs the current result into a plain double \c array, for easy communication using standard MPI types.  The array must have pre-allocated room for 2N+1 + 4 elements.
void PEResult::toDoubleArray(double* array) const {
	array[0] = double(status);
	array[1] = wavelength;
	array[2] = incidenceDeg;
	if(eff.empty()) {
		array[3] = 0;
	}
	else {
		array[3] = double((eff.size()-1)/2);	// N
		memcpy(array+4, &(eff.at(0)), eff.size()*sizeof(double));
	}
}

// This unpacks the result from a plain double \c array that was filled by toDoubleArray().
void PEResult::fromDoubleArray(const double* array) {
	status = Code(int(array[0]));
	wavelength = array[1];
	incidenceDeg = array[2];
	
	int N = int(array[3]);
	
	eff.resize(2*N+1);
	memcpy(&(eff[0]), array+4, eff.size()*sizeof(double));	
}

PEResult PEGrating::getEff(double incidenceDeg, double wl, const PEMathOptions& mo, bool printDebugOutput, int numThreads, bool measureTiming) const {
	PESolver s(*this, mo, numThreads, measureTiming);
	return s.getEff(incidenceDeg, wl, printDebugOutput);
}

gsl_complex PEGrating::refractiveIndex(double wl, const std::string& material) {
	// Override for testing:
//	return gsl_complex_rect(0.993, 0.00754);	// Pt at 410 eV
//	return gsl_complex_rect(1.4, -0.);	// plain glass in visible range.
//	return gsl_complex_rect(0.993, 0);

	// Attempt to open the material database file.
	std::string fileName = std::string(PEG_MATERIALS_DB_PATH) + std::string("/") + material + std::string(".idx");

	std::ifstream matFile;
	matFile.open(fileName.c_str());
	if(matFile.fail())
		return gsl_complex_rect(0,0);

	std::vector<double> wl_nm, delta, beta;
	wl_nm.reserve(3120);	// my standard files have 3111 lines. This avoids vector re-sizing for performance.
	delta.reserve(3120);
	beta.reserve(3120);

	// read whole file intp vectors for fast searching.
	while(!matFile.eof()) {
		double wl_nmIn, deltaIn, betaIn;
		matFile >> wl_nmIn >> deltaIn >> betaIn;

		wl_nm.push_back(wl_nmIn);
		delta.push_back(deltaIn);
		beta.push_back(betaIn);
	}
	matFile.close();

	if(wl_nm.empty())
		return gsl_complex_rect(0,0);	// the database is empty.

	// Need to convert wavelength from um to nm, since our data files are in the format:
	//   wl(nm)	delta	beta
	// where v = 1 - delta + i*beta       [beta assumed to be >= 0]
	wl *= 1000;	// now wl is in nm.


	// Use binary search to find the entry of the target wavelength in O(log(n)).
	////////////////////////

	// 'larger' is an iterator to one entry higher (or possibly equal) to the target wl.
	std::vector<double>::const_iterator larger = std::lower_bound(wl_nm.begin(), wl_nm.end(), wl);

	// Not found? This wl is larger than any we have in the database.
	if(larger == wl_nm.end())
		return gsl_complex_rect(0,0);

	// Did this return the first entry?
	if(larger == wl_nm.begin()) {
		// if equal to wl, then we've found it at the first entry.
		if(*larger == wl)
			return gsl_complex_rect(1-delta[0], beta[0]);
		else
			return gsl_complex_rect(0,0); // otherwise, it's too low a wavelength for us to have in the database.
	}


	// Not the first entry... So we can also grab the preceeding entry (at lowerIndex), and interpolate between them.
	size_t upperIndex = larger - wl_nm.begin();
	size_t lowerIndex = upperIndex - 1;

	double upperWl = wl_nm[upperIndex];
	double lowerWl = wl_nm[lowerIndex];
	double interp = (wl - lowerWl)/(upperWl - lowerWl);

	double upperDelta = delta[upperIndex];
	double lowerDelta = delta[lowerIndex];
	double interpDelta = lowerDelta + interp*(upperDelta-lowerDelta);

	double upperBeta = beta[upperIndex];
	double lowerBeta = beta[lowerIndex];
	double interpBeta = lowerBeta + interp*(upperBeta - lowerBeta);

	return gsl_complex_rect(1-interpDelta, 0.9*interpBeta);
}


std::ostream& operator<<(std::ostream& os, const PEResult& result) {
	int N = (result.eff.size()-1)/2;
	
	switch(result.status) {
	case PEResult::Success:
		os << "Inside Orders" << std::endl; 
		for(int i=0,cc=N; i<=cc; ++i) {
			os << i << "\t" << result.eff.at(N-i) << std::endl;
		}
		os << "\nOutside Orders" << std::endl; 
		for(int i=0,cc=N; i<=cc; ++i) {
			os << i << "\t" << result.eff.at(N+i) << std::endl;
		}
		break;
	case PEResult::InvalidGratingFailure:
		os << "Error: Invalid Grating" << std::endl;
		break;
	case PEResult::ConvergenceFailure:
		os << "Error: Integration Convergence Failure" << std::endl;
		break;
	case PEResult::InsufficientCoefficientsFailure:
		os << "Error: Insufficient Coefficients" << std::endl;
		break;
	case PEResult::AlgebraFailure:
		os << "Error: Linear Algebra Error" << std::endl;
		break;
	case PEResult::OtherFailure:
		os << "Error: Unknown other Failure" << std::endl;
		break;
	case PEResult::MissingRefractiveDataFailure:
		os << "Error: Missing refractive index data for this material at this wavelength." << std::endl;
		break;
	case PEResult::InactiveCalculation:
		os << "Notice: Inactive Calculation" << std::endl;
		break;
	}
	
	return os;
}

int PEGrating::computeK2StepsAtY(double y, gsl_complex k2_vaccuum, gsl_complex k2_substrate, gsl_complex k2_coating, double *stepsX, gsl_complex *stepsK2) const
{
	if(coatingThickness_ <= 0)
		return computeK2StepsAtY_noCoating(y, k2_vaccuum, k2_substrate, k2_coating, stepsX, stepsK2);
	else if(coatingThickness_ < profileHeight())
		return computeK2StepsAtY_interpenetratingCoating(y, k2_vaccuum, k2_substrate, k2_coating, stepsX, stepsK2);
	else
		return computeK2StepsAtY_thickCoating(y, k2_vaccuum, k2_substrate, k2_coating, stepsX, stepsK2);
}

int PEGrating::computeK2StepsAtY_noCoating(double y, gsl_complex k2_vaccuum, gsl_complex k2_substrate, gsl_complex k2_coating, double *stepsX, gsl_complex *stepsK2) const {

	(void)k2_coating;	// unused.

	double x1 = xIntersection1(y);
	double x2 = xIntersection2(y);
	double d = period();

	// Validity checks:
	////////////////////////////
	if(x1 < 0 || x2 < 0) {	// explicitly signals invalid geometry.
		std::cout << "Grating Expansion: Error: Invalid Geometry." << std::endl;
		return -1;	// invalid geometry.
	}
	if(x1 > d || x2 > d) {
		std::cout << "Grating Expansion: Error: Intersections are wider than period " << d << std::endl;
		return -1;	// invalid geometry.
	}

	// Sometimes x2 might drop below x1 if they're very close, and the integration routine goes slightly past the top of the grating. Let this slide.
	if(x2 < x1 - 1e-10) {
		std::cout << "Grating Expansion: Error: x2 is less than x1 by " << x1-x2 << std::endl;
		return -1;	// above the grating.
	}
	if(x2 < x1)
		x2 = x1;
	//////////////////////////////

	stepsX[0] = x1;
	stepsK2[0] = k2_vaccuum;

	stepsX[1] = x2;
	stepsK2[1] = k2_substrate;

	return 2;
}

int PEGrating::computeK2StepsAtY_interpenetratingCoating(double y, gsl_complex k2_vaccuum, gsl_complex k2_substrate, gsl_complex k2_coating, double *stepsX, gsl_complex *stepsK2) const {

	double d = period();

	// three regions based on y: less than the coating thickness, from coating thickness to profile height, and above.

	if(y <= coatingThickness_) {	// two crossings.  k^2 = k2_coating at first, k2_substrate at second.

		stepsX[0] = xIntersection1(y);
		stepsK2[0] = k2_coating;

		stepsX[1] = xIntersection2(y);
		stepsK2[1] = k2_substrate;

		// validity check:
		if(stepsX[0] < 0 || stepsX[1] < 0 || stepsX[0] > d || stepsX[1] > d || stepsX[1] < stepsX[0] - 1e-10)
			return -1;
		if(stepsX[1] < stepsX[0]) stepsX[1] = stepsX[0];

		return 2;
	}
	else if(y <= profileHeight()) {	// four crossings: vacuum, coating, substrate, coating.
		stepsX[0] = xIntersection1(y-coatingThickness_);
		stepsK2[0] = k2_vaccuum;

		stepsX[1] = xIntersection1(y);
		stepsK2[1] = k2_coating;

		stepsX[2] = xIntersection2(y);
		stepsK2[2] = k2_substrate;

		stepsX[3] = xIntersection2(y-coatingThickness_);
		stepsK2[3] = k2_coating;

		for(int i=0; i<4; ++i)
			if(stepsX[i] < 0 || stepsX[i] > d)
				return -1;

		return 4;
	}
	else {	// two crossings: vacuum, coating.
		stepsX[0] = xIntersection1(y-coatingThickness_);
		stepsK2[0] = k2_vaccuum;

		stepsX[1] = xIntersection2(y-coatingThickness_);
		stepsK2[1] = k2_coating;

		if(stepsX[0] < 0 || stepsX[0] > d || stepsX[1] < 0 || stepsX[1] > d || stepsX[1] < stepsX[0] - 1e-10)
			return -1;
		if(stepsX[1] < stepsX[0]) stepsX[1] = stepsX[0];

		return 2;
	}
}

int PEGrating::computeK2StepsAtY_thickCoating(double y, gsl_complex k2_vaccuum, gsl_complex k2_substrate, gsl_complex k2_coating, double *stepsX, gsl_complex *stepsK2) const {

	double d = period();

	// three regions based on y: crossing the substrate bump (y < profileHeight()), inside the homogeneous coating layer (y < coatingThickness_), and crossing the coating bump.

	if(y <= profileHeight()) {	// two crossings: coating, substrate.
		stepsX[0] = xIntersection1(y);
		stepsK2[0] = k2_coating;

		stepsX[1] = xIntersection2(y);
		stepsK2[1] = k2_substrate;

		// validity check:
		if(stepsX[0] < 0 || stepsX[1] < 0 || stepsX[0] > d || stepsX[1] > d || stepsX[1] < stepsX[0] - 1e-10) {
//			std::cerr << "err 1:" << stepsX[0] << " " << stepsX[1];
			return -1;
		}
		if(stepsX[1] < stepsX[0]) stepsX[1] = stepsX[0];

		return 2;
	}
	else if(y <= coatingThickness_) {	// homogeneous layer in coating.
		stepsX[0] = 0;
		stepsK2[0] = k2_coating;
		return 1;
	}
	else {	// two crossings: vaccuum, coating.
		stepsX[0] = xIntersection1(y-coatingThickness_);
		stepsK2[0] = k2_vaccuum;

		stepsX[1] = xIntersection2(y-coatingThickness_);
		stepsK2[1] = k2_coating;

		// validity check:
		if(stepsX[0] < 0 || stepsX[1] < 0 || stepsX[0] > d || stepsX[1] > d || stepsX[1] < stepsX[0] - 1e-10) {
//			std::cerr << "err 1:" << stepsX[0] << " " << stepsX[1];
			return -1;
		}
		if(stepsX[1] < stepsX[0]) stepsX[1] = stepsX[0];

		return 2;
	}
}

int PECustomProfileGrating::computeK2StepsAtY(double y, gsl_complex k2_vaccuum, gsl_complex k2_substrate, gsl_complex k2_coating, double *stepsX, gsl_complex *stepsK2) const
{
	// coatings are not supported.
	(void)k2_coating;
	if(coatingThickness_ != 0)
		return -1;

	if(!isValid())
		return -1;

	int numCrossings = 0;

	// move along points until y_i > y if searching for an entering point, or y_i <= y if searching for an exit point.
	for(int i=0,cc=y_.size(); i<cc; ++i) {

		if(numCrossings%2 == 0) {	// looking for an entry point.
			if(y_[i] > y) {	// found it
				double slope = (x_[i] - x_[i-1])/(y_[i]-y_[i-1]);
				double x = x_[i-1] + slope*(y - y_[i-1]);

				stepsK2[numCrossings] = k2_vaccuum;
				stepsX[numCrossings++] = x;
			}
		}
		else {// looking for an exit point.
			if(y_[i] <= y) {	// found it
				double slope = (x_[i] - x_[i-1])/(y_[i]-y_[i-1]);
				double x = x_[i-1] + slope*(y - y_[i-1]);

				stepsK2[numCrossings] = k2_substrate;
				stepsX[numCrossings++] = x;
			}
		}
	}

	if(numCrossings%2 != 0)
		return -1;	// there must be an even number of crossings: none, in-out, in-out-in-out, etc.

	// No crossing: we are above the grating. Homogeneous vacuum.
	if(numCrossings == 0) {
		numCrossings = 1;
		stepsK2[0] = k2_vaccuum;
		stepsX[0] = 0;
	}


	return numCrossings;
}

double PEGrating::roughnessFactor(double sigma, double wl, const std::string &material, double incidence)
{
	gsl_complex v = refractiveIndex(wl, material);

	// convert to radians:
	double sinTheta = sin(incidence * M_PI / 180.0);

	// calculate complex factor

	gsl_complex c = gsl_complex_sqrt(gsl_complex_sub(gsl_complex_mul(v,v), gsl_complex_rect(sinTheta*sinTheta,0)));

	return exp(-pow(4*M_PI*sigma/wl, 2)*sinTheta*GSL_REAL(c));
}
