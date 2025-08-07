/*
	img_random.h
	Header file for functions for creating random images
	Author: Bernard Heymann
	Created: 19990703 	Modified: 20030812
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

// Function prototypes
float* 		random_array_uniform(long n, float min, float max);
float 		random_gaussian(float avg, float std);
float* 		random_array_gaussian(long n, float avg, float std);
float		random_poisson(float avg);
float*		random_array_poisson(int n, float avg);
float 		random_logistical(float avg, float std);
float* 		random_array_logistical(long n, float avg, float std);
