#include "PEG.h"

/// This function tests a serial calculation of blazed grating efficiency.
int main(int argc, const char** argv) {
	
	PEBlazedGrating blazedGrating();
	
	PEResult r = blazedGrating.getEff();
	
	return 0;
}