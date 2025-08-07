/*
	img_util.h
	Header file for image utilities
	Author: Bernard Heymann
	Created: 20000430 	Modified: 20040622
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include "rwimg.h"
#include "matrix.h"

// Function prototypes
float		img_get_recip_resolution(Bimage* p);
int 	 	img_average(Bimage* p);
int 		img_reslice(Bimage* p, char* order);
int 		img_select(Bimage* p, int n);
int 		img_delete(Bimage* p, char* delete_list);
Bimage* 	img_copy_selected(Bimage* p, int n);
int 		img_translate(Bimage* p, VectorInt3 shift, int wrap, float fill);
Bimage*		img_rotate_by_pi(Bimage* p);
int 		img_adjust_FOM(Bimage* p, float fommin, float fommax, float fomcut);
int 		img_fill_upper_triangle(Bimage* p);
int 		img_transpose(Bimage* p);
int 		img_zero_origin(Bimage* p);
