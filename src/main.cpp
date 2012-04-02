#include "PEG.h"

/// This function tests a single calculation of blazed grating efficiency.
int main(int argc, const char** argv) {
	
	PEBlazedGrating blazedGrating;
	
	PEResult r = blazedGrating.getEff();
	
	return 0;
}
