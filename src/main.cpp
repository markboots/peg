#include "PEG.h"

/// This function tests a single calculation of blazed grating efficiency at 410 eV, 2.5deg. incidence, on a default BlazedGrating.
int main(int argc, const char** argv) {
	
	PEBlazedGrating blazedGrating;
	
	PEResult r = blazedGrating.getEff(2.5, 1.240/410);
	
	return 0;
}
