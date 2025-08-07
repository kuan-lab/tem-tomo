/*
	intersect.h
	Header file for handling
	intersection between convex polygons
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004
*/

#ifndef _INTERSECT_H__
#define _INTERSECT_H__

#include "compgeo.h"
#include "ptinpoly.h"

typedef enum { Pin, Qin, Unknown } tInFlag;

/*---------------------------------------------------------------------
Function prototypes.
---------------------------------------------------------------------*/
//Poly* poly_old_intersect(Poly *P, Poly *Q);
int 		poly_intersect(Poly *P, Poly *Q, Poly *Pol);
void    SubVec( tPointi a, tPointi b, tPointi c );
double  Dot( tPointi a, tPointi b );
char    SegSegInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p, tPointd q );
char    ParallelInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p, tPointd q );
bool    Between( tPointi a, tPointi b, tPointi c );
void    Assigndi( tPointd p, tPointi a );
tInFlag InOut( tPointd p, tInFlag inflag, int aHB, int bHA );
int     Advance( int a, int *aa, int n, bool inside, tPointi v, Poly* P );

#endif  // #ifndef _INTERSECT_H__

