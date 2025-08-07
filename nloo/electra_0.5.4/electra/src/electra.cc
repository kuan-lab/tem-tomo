/*
	electra.cc
	Show a list of the commands available
	Author: Giovanni Cardone
	Created: 20050329 	Modified: 20070719
*/

#include <stdio.h>

int		main(int argc, char* argv[])
{
   
	printf("\n---------------------------------------------------------\n");
	printf("       Electra package: list of available commands\n");
	printf("---------------------------------------------------------\n\n");
	printf("* et-fsceo\n  --------\n");
	printf("  Estimate resolution by even/odd Fourier Shell Correlation (FSCe/o).\n\n");
	printf("* et-fsceo_setup\n  --------------\n");
	printf("  Generate the tomograms needed for resolution analysis by FSCe/o.\n\n");
	printf("* et-nloo\n  -------\n");
	printf("  Estimate resolution by Noise-compensated Leave-One-Out (NLOO3D).\n\n");
	printf("* et-nloo_setup\n  -------------\n");
	printf("  Generate the projections needed for resolution analysis by NLOO3D.\n\n");
	printf("* et-tiltnloo\n  -----------\n");
	printf("  Evaluate dependence of resolution on the tilt angle by NLOO2D.\n\n");
	printf("* et-findbeads\n  -----------\n");
	printf("  Locate gold markers in tomograms and remove from corresponding tilt series.\n\n");
	printf("   (Warning: high memory requirements)\n\n");

	return(0);
}
