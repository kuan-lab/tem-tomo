/*
	simplex.h
	Nelder and Mead downhill simplex method for generalized parameter fitting
	Adapted from Numerical Recipes, 2nd edition, Press et al. 1992
	Author: Bernard Heymann
	Created: 20000426	Modified: 20041011
	
	The function "funk" is user-defined and references the "Bsimplex" structure.
	It returns the "R" value and is called as:
		float	R = (funk)(simplex_struct);
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/************************************************************************
@Object: struct Bsimplex
@Description:
	Structure used in the downhill simplex method.
@Features:
	Nelder and Mead downhill simplex method for generalized parameter fitting
	Adapted from Numerical Recipes, 2nd edition, Press et al. 1992
	The structure is set up to accommodate any number of variables, parameters,
	constants and points.
	The structure is very flexible in the sense that only some fields
	are absolutely required and with a fixed meaning for the simplex method.
	The required fields are:
		nparam, param, lolimit, hilimit
	The other fields may be recast and used as desired in the user function.
	Intended sizes:
		param		nparam.
		lolimit 	nparam.
		hilimit 	nparam.
		constant	nconstant.
		x			npoint*nvar.
		fx			npoint.
	x or fx can be recast as a different pointer, as long as it is handled
	by the user before calling kill_simplex.
*************************************************************************/
#ifndef _Bsimplex_
struct Bsimplex {
	int		nvar;		// Number of variables
	int		nparam;		// Number of parameters
	int		nconstant;	// Number of constants
	int		npoint; 	// Number of function values
	float*	param;		// Parameter values
	float*	lolimit;	// Lower limits on parameter values
	float*	hilimit;	// Upper limits on parameter values
	float*	constant;	// Constant values
	float*	x;			// Independent variables: npoint*nvar array
	float*	fx;			// Function values
} ;
#define _Bsimplex_
#endif

// Function prototypes
double		simplex(Bsimplex* simp, int maxcycles, float tolerance, double (funk)(Bsimplex *));
Bsimplex*	init_simplex(int nvar, int nparam, int nconstant, int npoint);
int 		kill_simplex(Bsimplex* simp);
