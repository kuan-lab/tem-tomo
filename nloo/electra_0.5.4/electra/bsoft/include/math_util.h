/* 
	math_util.h 
	Header file for general utilities 
	Author: Bernard Heymann 
	Created: 19990722 	Modified: 20040507
*/ 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#ifdef SUNOS
#include <ieeefp.h>
#endif

// Function prototypes 
double			bfloor(double value, int places);
double			bround(double value, int places);
float			gamma_ln(float x);
double			factorial(int n);
unsigned long*  prime_factors(unsigned long number, int* nfactors);
unsigned long	smallest_prime(unsigned long number);

