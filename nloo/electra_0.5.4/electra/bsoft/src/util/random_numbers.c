/*
	random_numbers.c
	Functions for generating random sequences
	Author: Bernard Heymann
	Created: 19990703 	Modified: 20031007
*/

#include "math_util.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

/************************************************************************
@Function: random_array_uniform
@Description:
	Generates a series with a uniform random distribution.
@Algorithm:
	An array of floating point numbers is generated distributed uniformly 
	in the range of the given minimum and maximum:
		value = random_value*(max - min) + min
	where random_value is between 0 and 1.
	The average and standard deviation are:
		average = (max + min)/2
		standard deviation = 0.5*sqrt(1/3)*(max - min).
@Arguments:
	long n				number of values.
	float min 			minimum value.
	float max 			maximum value.
@Returns:
	float*				array with uniform random numbers.
**************************************************************************/
float* 		random_array_uniform(long n, float min, float max)
{
	float* 			r = (float *) balloc(n*sizeof(float));
	 
	long 			i;
	double			range = (max - min)/get_rand_max();
	
	if ( verbose & VERB_FULL ) {
		printf("Generating a random series with a uniform distribution:\n");
		printf("Size:                           %ld\n", n);
		printf("Minimum & maximum:              %g %g\n\n", min, max);
	}
		
//	srandom(time(NULL));
	srandom(getpid());
	
	for ( i=0; i<n; i++ )
		r[i] = random()*range + min;
	
	return(r);
}

/************************************************************************
@Function: random_gaussian
@Description:
	Generates a value with a gaussian random distribution.
@Algorithm:
	A floating point number is generated with a gaussian 
	distribution with a given average and standard deviation:
		value = average + std_dev*sqrt(-2*log(random_value))*
						cos(2*PI*random_value);
	where random_value is between 0 and 1.
@Arguments:
	float avg 			average.
	float std 			standard deviation.
@Returns:
	float				the random value.
**************************************************************************/
float 		random_gaussian(float avg, float std)
{
	float			r, invrandmax = 1.0/get_rand_max();
	
	if ( verbose & VERB_FULL ) {
		printf("Generating a random value with a gaussian distribution:\n");
		printf("Avg and std:                    %g %g\n\n", avg, std);
	}
		
	r = avg + std*sqrt(-2*log(random()*invrandmax))*
				cos(2*PI*random()*invrandmax);
	
	return(r);
}

/************************************************************************
@Function: random_array_gaussian
@Description:
	Generates a series with a gaussian random distribution of values.
@Algorithm:
	An array of floating point numbers is generated with a gaussian 
	distribution with a given average and standard deviation:
		value = average + std_dev*sqrt(-2*log(random_value))*
						cos(2*PI*random_value);
	where random_value is between 0 and 1.
@Arguments:
	long n				number of values.
	float avg 			average.
	float std 			standard deviation.
@Returns:
	float*				the array of values.
**************************************************************************/
float* 		random_array_gaussian(long n, float avg, float std)
{
	float* 			r = (float *) balloc(n*sizeof(float));
	 
	long 			i;
	float			invrandmax = 1.0/get_rand_max();
	
	if ( verbose & VERB_FULL ) {
		printf("Generating a random series with a gaussian distribution:\n");
		printf("Size:                           %ld\n", n);
		printf("Avg and std:                    %g %g\n\n", avg, std);
	}
		
//	srandom(time(NULL));
	srandom(getpid());
	
	for ( i=0; i<n; i++ )
		r[i] = avg + std*sqrt(-2*log(random()*invrandmax))*
				cos(2*PI*random()*invrandmax);
	
	return(r);
}

/************************************************************************
@Function: random_poisson
@Description:
	Generates a value deviating from the average based on a poisson distribution.
@Algorithm:
	Algorithm taken from Numerical Recipes in C.
	The poisson distribution is given for j = 0,1,... by:
		        avg^j * exp(-avg)
		P(j) = -----------------
		               j!
	Note that only positive integer values are defined for j and sum(P(j)) = 1.
	A value is generated with a poisson distribution with a given average.
	If the average <= 0, the return value is zero.
@Arguments:
	float avg 			average.
@Returns:
	float				value.
**************************************************************************/
float		random_poisson(float avg)
{
	float 			r = 0;
	
	if ( avg <= 0 ) return(r);
	 
	double			irm = 1.0/get_rand_max();
	double			sq, lnavg, g, t, y;

	if ( avg < 12.0 ) {
		g = exp(-avg);
		for ( r = -1, t = 1; t > g; r++ ) t *= random()*irm;
	} else {
		sq = sqrt(2.0*avg);
		lnavg = log(avg);
		g = avg*lnavg - gamma_ln(avg + 1);
		do {
			for ( r = -1; r < 0; ) {
				y = tan(PI*random()*irm);
				r = sq*y + avg;
			}
			r = floor(r);
			t = 0.9*(1 + y*y)*exp(r*lnavg - gamma_ln(r + 1) - g);
		} while ( random()*irm > t );
	}
	
	return(r);
}

/************************************************************************
@Function: random_array_poisson
@Description:
	Generates a series with a poisson random distribution of values.
@Algorithm:
	Algorithm taken from Numerical Recipes in C.
	The poisson distribution is given for j = 0,1,... by:
		        avg^j * exp(-avg)
		P(j) = -----------------
		               j!
	Note that only positive integer values are defined for j and sum(P(j)) = 1.
	An array of floating point numbers is generated with a poisson 
	distribution with a given average. The standard deviation is:
		std = sqrt(avg)
	If the average <= 0, the return array contains only zeroes.
@Arguments:
	long n				number of values.
	float avg 			average.
@Returns:
	float*				the array of values.
**************************************************************************/
float*		random_array_poisson(int n, float avg)
{
	float* 			r = (float *) balloc(n*sizeof(float));
	
	if ( avg <= 0 ) return(r);
	 
	long 			i;
	double			irm = 1.0/get_rand_max();
	double			sq, lnavg, g, t, y;

//	srandom(time(NULL));
	srandom(getpid());
	
	if ( avg < 12.0 ) {
		g = exp(-avg);
		for ( i=0; i<n; i++ )
			for ( r[i] = -1, t = 1; t > g; r[i]++ ) t *= random()*irm;
	} else {
		sq = sqrt(2.0*avg);
		lnavg = log(avg);
		g = avg*lnavg - gamma_ln(avg + 1);
		for ( i=0; i<n; i++ ) {
			do {
				for ( r[i] = -1; r[i] < 0; ) {
					y = tan(PI*random()*irm);
					r[i] = sq*y + avg;
				}
				r[i] = floor(r[i]);
				t = 0.9*(1 + y*y)*exp(r[i]*lnavg - gamma_ln(r[i] + 1) - g);
			} while ( random()*irm > t );
		}
	}
	
	return(r);
}

/************************************************************************
@Function: random_logistical
@Description:
	Generates a value with a logistical random distribution.
@Algorithm:
	A floating point number is generated with a logistical 
	differential distribution with a given average and standard deviation:
		value = average + (std_dev/golden)*ln(1/random_value - 1)
	where random_value is between 0 and 1 and:
		golden  = (sqrt(5) + 1)/2
@Arguments:
	float avg 			average.
	float std 			standard deviation.
@Returns:
	float				the random value.
**************************************************************************/
float 		random_logistical(float avg, float std)
{
	double			r, rand_max = get_rand_max();
	double			golden = 0.5*(sqrt(5.0)+1);
	double			prefac = std/golden;
	
	if ( verbose & VERB_FULL ) {
		printf("Generating a value with a logistical differential distribution:\n");
		printf("Avg, std and range:             %g %g\n", avg, std);
	}
	
	r = avg + prefac*log(rand_max/(random()+1) - 1);
	
	return(r);
}

/************************************************************************
@Function: random_array_logistical
@Description:
	Generates an array with a logistical random distribution.
@Algorithm:
	An array of floating point numbers is generated with a logistical 
	differential distribution with a given average and standard deviation:
		value = average + (std_dev/golden)*ln(1/random_value - 1)
	where random_value is between 0 and 1 and:
		golden  = (sqrt(5) + 1)/2
@Arguments:
	long n				number of values.
	float avg 			average.
	float std 			standard deviation.
@Returns:
	float*				array of values.
**************************************************************************/
float* 		random_array_logistical(long n, float avg, float std)
{
	float* 			r = (float *) balloc(n*sizeof(float));
	 
	long 			i;
	double			rand_max = get_rand_max();
	double			golden = 0.5*(sqrt(5.0)+1);
	double			prefac = std/golden;
	
	if ( verbose & VERB_FULL ) {
		printf("Generating a random series with a logistical differential distribution:\n");
		printf("Size:                           %ld\n", n);
		printf("Avg, std and range:             %g %g\n", avg, std);
	}
	
//	srandom(time(NULL));
	srandom(getpid());
	
	for ( i=0; i<n; i++ )
		r[i] = avg + prefac*log(rand_max/(random()+1) - 1);
	
	return(r);
}


