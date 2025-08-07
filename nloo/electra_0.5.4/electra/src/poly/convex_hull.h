/*
	convex_hull.h
	Header file for convex hull algorithm
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004
*/

#ifndef _CONVEXHULL_H__
#define _CONVEXHULL_H__

#include <stdlib.h>
#include <stdio.h>
#include "poly_def.h"
#include "poly_util.h"


/*----------Function Prototypes-------------*/
int poly_create_convex_hull(Poly *P);
int ccw(coord *P, int i, int j, int k);
int cmpl(const void *a, const void *b);
int cmph(const void *a, const void *b);
int make_chain(coord* V, int n, int (*cmp)(const void*, const void*));

#endif  // #ifndef _CONVEXHULL_H__
