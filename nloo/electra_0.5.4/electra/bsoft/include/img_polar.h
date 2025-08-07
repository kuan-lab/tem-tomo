/*
	img_polar.h
	Library routines for polar and spherical transformations and calculations
	Author: Bernard Heymann
	Created: 19990904	Modified: 20030721
*/

#include "rwimg.h"

// Function prototypes
int 		img_cartesian_to_spherical(Bimage* p,
					int nannuli, int ntheta, int nphi);
int		 	img_radial_shells(Bimage* p);
Bimage* 	img_polar_power_spectrum(Bimage* p, int nannuli, int npsi, 
				int transform_flag, int norm_flag);
Bimage*		img_polar2D(Bimage* p, int nannuli, int nseg);
