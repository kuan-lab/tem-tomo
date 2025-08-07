/*
	convex_hull.cc
	Functions for handling polygons
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004

	All the functions are adapted from existing code.
*/
/*
 * Ken Clarkson wrote this.  Copyright (c) 1996 by AT&T..
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */
/*
 * 	Algorihm somewhat like Procedure 8.2 in Edelsbrunner's "Algorithms in Combinatorial
 *	Geometry", and very close to:
 *	    A.M. Andrew, "Another Efficient Algorithm for Convex Hulls in Two Dimensions",
 *		Info. Proc. Letters 9, 216-219 (1979)
 *	(See also http://geometryalgorithms.com/Archive/algorithm_0110/algorithm_0110.htm)
*/

#include "convex_hull.h"

/************************************************************************
@Function: poly_create_convex_hull
@Description:
	Determines the convex hull associated to given points
@Algorithm:
	a) sort the points from left to right
	b) starting from the leftmost point (A), find next point (B) such that
	   all the set lies above the line (AB).
	   Some "inner" points are discarded, as they are violating the
	   convexity property.
	   Then, start over with point B until you've reach the rightmost point.
	   You then have the bottom "wrapping" of the convex hull.
	c) Repeat the two passes, only inverting the sorting order,
	   to find the upper boundary. 

	The original array storing the points is overwritten.
	IMPORTANT: the size of the original array must exceed by one the size
			   of the input points.
@Arguments:
	Poly* P		input: coordinates of points to be contained in the convex hull
						output: convex hull vertices
@Returns:
	int					number of output vertices.
**************************************************************************/
int poly_create_convex_hull(Poly *P)
{
	int nin = P->n;
	if (!nin) return 0;
	
	if ( verbose & VERB_FULL ) {
		printf("starting points:\n");
		poly_print(P);
	}
	
	/* lower bounds */
	int nout = make_chain(P->vertex, nin, cmpl);
//	P[nin] = P[0];
	cp_coord(P->vertex[0],P->vertex[nin]);

	/* upper bounds */
	nout += make_chain(P->vertex+nout, nin-nout+1, cmph);
	P->n = nout;

	if ( verbose & VERB_FULL ) {
		printf("convex polygon created:\n");
		poly_print(P);
	}
	
	poly_init_bounding_box(P);

	return nout;
}

int ccw(coord *P, int i, int j, int k) {
	int	a = P[i][X] - P[j][X],
		b = P[i][Y] - P[j][Y],
		c = P[k][X] - P[j][X],
		d = P[k][Y] - P[j][Y];
	return a*d - b*c <= 0;	   /* true if points i, j, k counterclockwise */
}

int cmpl(const void *a, const void *b) {
	double v; 
	v = (*(coord*)a)[X] - (*(coord*)b)[X];
	if (v>0) return 1;
	if (v<0) return -1;
	v = (*(coord*)b)[Y] - (*(coord*)a)[Y];
	if (v>0) return 1;
	if (v<0) return -1;
	return 0;
}

int cmph(const void *a, const void *b) {return cmpl(b,a);}


int make_chain(coord* V, int n, int (*cmp)(const void*, const void*)) {
	int i, j, s = 1;
	coord t;

	qsort(V, n, sizeof(coord), cmp);
	for (i=2; i<n; i++) {
		for (j=s; j>=1 && ccw(V, i, j, j-1); j--){}
		s = j+1;
//		t = V[s]; V[s] = V[i]; V[i] = t;
		cp_coord(V[s],t);
		cp_coord(V[i],V[s]);
		cp_coord(t,V[i]);
	}
	return s;
}
