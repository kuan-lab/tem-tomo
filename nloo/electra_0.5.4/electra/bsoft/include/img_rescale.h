/*
	img_rescale.h
	Header file for image utilities
	Author: Bernard Heymann
	Created: 20000430 	Modified: 20031202
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
int 		img_invert_data(Bimage* p);
int 		img_rescale(Bimage* p, float scale, float shift);
int 		img_rescale_to_min_max(Bimage* p, float min, float max);
int 		img_rescale_to_avg_std(Bimage* p, float avg, float std);
int 		img_square(Bimage* p);
int 		img_logarithm(Bimage* p);
int 		img_truncate(Bimage* p, float min, float max, float setmin, float setmax);
int 		img_truncate_to_min_max(Bimage* p, float min, float max);
int 		img_truncate_to_avg(Bimage* p, float min, float max);
int 		img_limit_levels(Bimage* p, int nlevels);
