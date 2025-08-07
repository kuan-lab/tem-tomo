/*
	poly_def.h
	Header file for coordinates definition
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004
*/

#ifndef _POLYDEF_H__
#define _POLYDEF_H__

#define X   0
#define Y   1

#define DIM 2               /* Dimension of points */
//typedef int coord;
typedef int coord[DIM];

/*----------Polygon(s) Structure-------------*/
/************************************************************************
@Object: struct poly
@Description:
	Convex polygon parameter structure.
@Features:
	This contains all the information pertinent to a convex polygon.
*************************************************************************/
struct Poly {
	coord*		vertex; 	    // Vertices coordinates
	int			n;				// Number of vertices
	int			nt;				// Number of maximum number of allowed vertex elements
	int			bbxmin;			// Bounding box - minimum x
	int			bbymin;			// Bounding box - minimum y
	int			bbxmax;			// Bounding box - maximum x
	int			bbymax;			// Bounding box - maximum y
} ;

#endif  // #ifndef _POLYDEF_H__
