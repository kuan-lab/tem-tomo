/*
	distance.h
	Header file for distance measures
	between point andpolygons
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004
*/

#ifndef _DISTANCE_H__
#define _DISTANCE_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "poly_def.h"


/*----------Function Prototypes-------------*/
float poly_distance_to_polygon( Poly* P, coord pnt);
float poly_distance_to_segment(int x0, int y0, int x1, int y1, int x, int y);

#endif  // #ifndef _DISTANCE_H__
