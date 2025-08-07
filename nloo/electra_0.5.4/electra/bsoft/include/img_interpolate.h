/*
	interpolate.h
	Library routines for interpolating images
	Author: Bernard Heymann
	Created: 19990904	Modified: 20040204
*/

#include "rwimg.h"
#include "matrix.h"

// Function prototypes
int 		img_transform(Bimage* p, VectorInt3 newsize, Vector3 scale, 
				Vector3 origin, Vector3 translate, Matrix3 mat, float fill);
Bimage*		img_rotate_to_view(Bimage* p, View view);
int			img_rotate_using_euler_angles(Bimage* p, Euler euler);
int			img_shift(Bimage* p, Vector3* shift);
int			img_unique_shift_global_rotate(Bimage* p, Vector3* origin, Vector3* shift, float angle);
int 		img_integer_interpolation(Bimage* p, int integer_factor);
int 		img_bin(Bimage* p, VectorInt3 bin);
int 		img_median_bin(Bimage* p, int binning);
