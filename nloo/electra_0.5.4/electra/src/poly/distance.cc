/*
	distance.cc
	Functions for measuring distance of a point
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004
*/

#include "distance.h"

/************************************************************************
@Function: poly_distance_to_polygon
@Description:
	Evaluates the distance of a point from a polygon.
@Algorithm:
	The distance between the test point and each line segment composing
	the polygon is calculated.
	Then the lowest value is returned.
@Arguments:
	Poly P				reference polygon.
	coord pnt			test point.
@Returns:
	float				distance measure.
**************************************************************************/
float poly_distance_to_polygon( Poly* P, coord pnt)
{
	int n = P->n;
	float dd = poly_distance_to_segment(P->vertex[n-1][X],P->vertex[n-1][Y],
                                        P->vertex[0][X],P->vertex[0][Y],
										pnt[X],pnt[Y]);
	for (int i=0; i<n-1; ++i) {
		float nd = poly_distance_to_segment(P->vertex[i][X], P->vertex[i][Y],
                                             P->vertex[i+1][X], P->vertex[i+1][Y],
                                             pnt[X],pnt[Y]);
		if (nd<dd)
		dd = nd;
	}

	return dd;
}
 
/************************************************************************
@Function: poly_distance_to_segment
@Description:
	Evaluates the distance of a point from a segment.
@Algorithm:
	The distance between the test point and each line segment composing
	the polygon is calculated.
	Then the lowest value is returned.
@Arguments:
	Poly P				reference polygon.
	coord pnt			test point.
@Returns:
	float				distance measure.
**************************************************************************/
float poly_distance_to_segment(int x0, int y0,
                                int x1, int y1,
                                int x, int y)
{
//	printf("%d,%d %d,%d\n",x0,y0,x1,y1);
	// squared distance between endpoints :
	long ddh = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0);
 
	// squared distance to endpoints :
	long dd0 = (x-x0)*(x-x0) + (y-y0)*(y-y0);
	long dd1 = (x-x1)*(x-x1) + (y-y1)*(y-y1);
 
	// if closest to the start point :
	if (dd1 > ddh + dd0)
		return sqrt((float) dd0);
 
	// if closest to the end point :
	if (dd0 > ddh + dd1)
		return sqrt((float) dd1);
 
	// squared perpendicular distance to line :
	float a = y0-y1;
	float b = x1-x0;
	float c = x0*y1-x1*y0;
	float ddn = (a*x + b*y + c)*(a*x + b*y + c)/(a*a + b*b);
//	printf("%f,%f,%f\n",a,b,a*a + b*b);
	return sqrt(ddn);
}
