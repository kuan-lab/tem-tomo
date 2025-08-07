/*
	img_symmetry.h
	Header file for symmetry functions
	Author: Bernard Heymann
	Created: 19990509 	Modified: 20030421
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "rwimg.h"
#include "matrix.h"

// Function prototypes in symmetry.c and rwsymop.c
int 		img_apply_point_group(Bimage* p, char* symmetry_string, View ref_view);
int 		img_add_with_rotation(Bimage* p, Bimage* psum, Vector3 axis, float angle);
int			img_write_symmetry_views(Bimage* p, char* symmetry_string, View ref_view,
				char* filename, DataType datatype, float avg, float std);
float* 		img_check_space_group(Bimage* p);
float* 		img_apply_space_group(Bimage* p);
int 		img_apply_friedel(Bimage* p);
float 		img_check_handedness(Bimage* p, float threshold);
