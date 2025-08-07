/*
	img_random.h
	Header file for functions for creating random images
	Author: Bernard Heymann
	Created: 19990703 	Modified: 20030915
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
Bimage* 	img_random_uniform(int nimg, VectorInt3 size, float min, float max);
Bimage* 	img_random_gaussian(int nimg, VectorInt3 size, float avg, float std);
Bimage* 	img_random_poisson(int nimg, VectorInt3 size, float avg);
Bimage* 	img_random_logistical(int nimg, VectorInt3 size, float avg, float std);
int 		img_add_gaussian_noise(Bimage* p, float snr);
int 		img_add_poisson_noise(Bimage* p);
