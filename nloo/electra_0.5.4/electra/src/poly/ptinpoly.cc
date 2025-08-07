/*
	inpoly.cc
	Functions for determining
	if a point is inside a convex polygon
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004

	The routine is adapted from Haines code.
*/

#include "ptinpoly.h"

/************************************************************************
@Function: poly_point_inside
@Description:
	Test of inclusion of a point in a convex polygon.
@Algorithm:
    Counts the crossing made by a ray from the test point.
	Shoot a test ray along +X axis.  The strategy, from MacMartin, is to
	compare vertex Y values to the testing point's Y and quickly discard
	edges which are entirely to one side of the test ray.
	Haines, Eric, "Point in Polygon Strategies,"
	 Graphics Gems IV, ed. Paul Heckbert, Academic Press, p. 24-46, 1994
@Arguments:
	Poly pgon			2D convex polygon
	int* point			test point
@Returns:
	int					1 if inside, 0 if outside
**************************************************************************/
int poly_point_inside( Poly* pgon, int* point )
{
	register int	j, yflag0, yflag1, inside_flag, xflag0 ;
//	register double	ty, tx, *vtx0, *vtx1 ;
	register int	ty, tx;
	register int    *vtx0, *vtx1 ;
	register int	line_flag ;

    tx = point[X];
    ty = point[Y];

    vtx0 = pgon->vertex[pgon->n-1] ;
    /* get test bit for above/below X axis */
    yflag0 = ( vtx0[Y] >= ty ) ;
    vtx1 = pgon->vertex[0] ;

    inside_flag = 0 ;
    line_flag = 0 ;

    for ( j = pgon->n+1 ; --j ; ) {

		yflag1 = ( vtx1[Y] >= ty ) ;
		/* check if endpoints straddle (are on opposite sides) of X axis
		* (i.e. the Y's differ); if so, +X ray could intersect this edge.
		*/
		if ( yflag0 != yflag1 ) {
			xflag0 = ( vtx0[X] >= tx ) ;
			/* check if endpoints are on same side of the Y axis (i.e. X's
			* are the same); if so, it's easy to test if edge hits or misses.
			*/
			if ( xflag0 == ( vtx1[X] >= tx ) ) {
			/* if edge's X values both right of the point, must hit */
			if ( xflag0 ) inside_flag = !inside_flag ;
			} else {
			/* compute intersection of pgon segment with +X ray, note
			* if >= point's X; if so, the ray hits it.
			*/
				if ( (vtx1[X] - (vtx1[Y]-ty)*
					(double)(vtx0[X]-vtx1[X])
					/(double)(vtx0[Y]-vtx1[Y])) >= tx ) {
					inside_flag = !inside_flag ;
				}
			}
			/* if this is second edge hit, then done testing */
			if ( line_flag ) break;

			/* note that one edge has been hit by the ray's line */
			line_flag = 1;
		}

		/* move to next pair of vertices, retaining info as possible */
		yflag0 = yflag1 ;
		vtx0 = vtx1 ;
		vtx1 += 2 ;
    }
    return( inside_flag ) ;
}
