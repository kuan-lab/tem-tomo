/*
	img_combine.h
	Header file for combining two images in various ways
	Author: Bernard Heymann
	Created: 20000430 	Modified: 20041225
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "rwimg.h"

// Function prototypes
int 		img_check_if_same_size(Bimage* p1, Bimage* p2);
int 		img_compatibility(Bimage* p1, Bimage* p2);
int 		add_images(Bimage* p1, Bimage* p2, float scale, float shift);
Bimage* 	img_add_weighed(int nfiles, char* file_list, float* weight,
					float newavg, float newstd);
int 		img_place(Bimage* p, Bimage* psmall, Vector3 origin, float radius, float scale, float shift);
int 		multiply_images(Bimage* p1, Bimage* p2, float scale, float shift);
int 		largest_img(Bimage* p1, Bimage* p2, float scale, float shift);
int 		smallest_img(Bimage* p1, Bimage* p2, float scale, float shift);
float*		img_linear_fit(Bimage* p1, Bimage* p2, Bimage* pmask, float max_exclude);
float		img_R_factor(Bimage* p1, Bimage* p2);
double		img_one_correlation_coefficient(Bimage* p1, Bimage* p2);
float*		img_correlation_coefficient(Bimage* p1, Bimage* p2);
float		img_correlation_coefficient_radial(Bimage* pA, Bimage* pB, float radius_min, float radius_max);
double		correlate_images(Bimage* p1, Bimage* p2, float bas1, float bas2);
double		correlate_images_within_radii(Bimage* p1, Bimage* p2, 
				float bas1, float bas2, float minr, float maxr);
int 		histomatch_images(Bimage* p1, Bimage* p2, int bins);
int 		img_mask_with_image(Bimage* p, Bimage* m, float threshold, float fill);
Bimage* 	img_morph_blend(Bimage* p1, Bimage* p2, int number);
