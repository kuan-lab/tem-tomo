/*
	img_complex.h
	Functions to handle complex data sets
	Author: Bernard Heymann
	Created: 19990424  	    Modified: 20041113
*/

#include "rwimg.h"

// Function prototypes
int			img_to_complex_short(Bimage* p);
int			img_to_complex_int(Bimage* p);
int			img_to_complex_float(Bimage* p);
int			img_to_polar(Bimage* p);
int 		img_complex_product(Bimage* p, Bimage* p2);
int 		img_complex_apply_mask(Bimage* p, Bimage* pmask);
Bimage*		img_pack_two_in_complex(Bimage* p1, Bimage* p2);
Bimage* 	img_unpack_combined_transform(Bimage* p);
int 		img_combined_complex_product(Bimage* p, float hires, float lores);
int 		img_combined_masked_complex_product(Bimage* p, float hires, 
				float lores, Bimage* pmask);
int 		img_complex2real(Bimage* p);
int 		img_complex2amplitudes(Bimage* p);
int 		img_complex2intensities(Bimage* p);
int 		img_complex2phases(Bimage* p);
int 		img_complex2polar(Bimage* p);
int 		img_polar2complex(Bimage* p);
double		img_merge_amplitudes_and_phases(Bimage* pamp, Bimage* p);

