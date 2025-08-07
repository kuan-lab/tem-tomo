/*
	moving_average.h
	Header file moving average calculations
	Author: Bernard Heymann
	Created: 20000430 	Modified: 20030414
*/

#include "utilities.h"

// Function prototypes 
float*		moving_average(int number, float* x, int window);
complex_float*	moving_average_complex(int number, complex_float* x, int window);
polar*		moving_average_polar(int number, polar* x, int window);

