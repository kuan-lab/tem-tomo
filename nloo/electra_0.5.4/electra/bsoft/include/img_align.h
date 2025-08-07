/*
	img_align.h
	Header file for functions to align images.
	Author: Bernard Heymann
	Created: 20000505	Modified: 20041129
*/

#include "rwimg.h"
#include "matrix.h"

// Function prototypes
float*		img_align(Bimage* p1, Bimage* p2, VectorInt3 tile_size, 
				float res_lo, float res_hi, int filter_flag, int refine_flag);
Bimage*		img_sum_two(Bimage* p1, Bimage* p2);


