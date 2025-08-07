/*
	poly_util.cc
	Basic functions for manipulating polygons
	Author: Giovanni Cardone
	Created: 2004 	Modified: 20050131
*/

#include "poly_util.h"

/************************************************************************
@Function: cp_coord
@Description:
	copy coordinate values from an input point to an output point
@Algorithm:
@Arguments:
	coord csrc		input coordinates
	coord cdest		output coordinates
@Returns:
	void				
**************************************************************************/
void cp_coord(int* csrc, int* cdest) {
	cdest[X] = csrc[X];
	cdest[Y] = csrc[Y];
}

/************************************************************************
@Function: double2coord
@Description:
	cast coordinates from double to coord (int)
@Algorithm:
@Arguments:
	double* dc		coordinates to transform
@Returns:
	coord*			transformed coordinates
**************************************************************************/
coord* double2coord(double* dc)
{
	coord* ic = (coord *) balloc(sizeof(coord));
	*ic[X] = (int) dc[X];
	*ic[Y] = (int) dc[Y];
	return (ic);
}

int  poly_add_coord(coord p, Poly* Pol)
{
	if(Pol->n >= Pol->nt) return(-1);
	
	cp_coord( p, Pol->vertex[Pol->n]);
	Pol->n++;
	poly_update_bounding_box(p, Pol);
	
	return (0);
}

int  poly_add_point(int px, int py, Poly* Pol)
{
	coord pp;
	pp[X] = px;
	pp[Y] = py;
	
	return poly_add_coord(pp, Pol);
}

int  poly_init_bounding_box(Poly* P)
{
	int xmin = BBHUGE;
	int xmax = -BBHUGE;
	int ymin = BBHUGE;
	int ymax = -BBHUGE;

	for ( int i=0; i<P->n; i++ ) {
		int px = P->vertex[i][X];
		int py = P->vertex[i][Y];
		if(px > xmax) xmax = px;
		if(px < xmin) xmin = px;
		if(py > ymax) ymax = py;
		if(py < ymin) ymin = py;
	}

	P->bbxmin = xmin;
	P->bbxmax = xmax;
	P->bbymin = ymin;
	P->bbymax = ymax;

	return (0);
}

int poly_update_bounding_box(coord p, Poly* P)
{
	int xmin = P->bbxmin;
	int xmax = P->bbxmax;
	int ymin = P->bbymin;
	int ymax = P->bbymax;

	int px = p[X];
	int py = p[Y];
	if(px > xmax) xmax = px;
	if(px < xmin) xmin = px;
	if(py > ymax) ymax = py;
	if(py < ymin) ymin = py;
	
	P->bbxmin = xmin;
	P->bbxmax = xmax;
	P->bbymin = ymin;
	P->bbymax = ymax;

	return (0);
}

void   poly_print( Poly *P )
{
   int   i;

   printf("Polygon vertices:\n");
   printf("  i\tx   y\n");
   for( i = 0; i < P->n; i++ )
      printf("%3d\t%4d %4d\n", i, P->vertex[i][X], P->vertex[i][Y]);
}

/*
   Returns true iff c is strictly to the left of the directed
   line through a to b.
*/
bool    Left( coord a, coord b, coord c )
{
        return  AreaSign( a, b, c ) > 0;
}

bool    LeftOn( coord a, coord b, coord c )
{
        return  AreaSign( a, b, c ) >= 0;
}

bool    Collinear( coord a, coord b, coord c )
{
        return  AreaSign( a, b, c ) == 0;
}

int	AreaSign( coord a, coord b, coord c )
{
    double area2;

    area2 = ( b[0] - a[0] ) * (double)( c[1] - a[1] ) -
            ( c[0] - a[0] ) * (double)( b[1] - a[1] );

    /* The area should be an integer. */
    if      ( area2 >  0.5 ) return  1;
    else if ( area2 < -0.5 ) return -1;
    else                     return  0;
}

Poly* poly_init(int n)
{
	char*		ptr = NULL; 
	
	if ( ( ptr = (char *) malloc(sizeof(Poly)*sizeof(char)) ) == NULL ) { 
		fprintf(stderr,"\nError: Memory allocation of %ld bytes failed!\n",  sizeof(Poly)); 
		return(NULL); 
	}

	// Allocate memory for the image parameter structure
	Poly* 	p = (Poly *) ptr;

	// Set parameter defaults
	if ( (p->vertex = (coord *) malloc(n*sizeof(coord)) ) == NULL ) { 
		fprintf(stderr,"\nError: Memory allocation of %ld bytes failed!\n", n*sizeof(coord)); 
		return(NULL); 
	}
	p->n = 0;
	p->nt = n;
	
	memset(p->vertex, 0, n*sizeof(coord)); 

	poly_init_bounding_box(p);

	if ( verbose & VERB_DEBUG )
		printf("DEBUG poly_init: polygon setup done\n");
	
	return(p);


}

int poly_destroy(Poly* p)
{

	if ( p == NULL ) return(0);

	if ( p->vertex ) free((void *)p->vertex);
	p->vertex = NULL;
	
	free((void *)p);
	p = NULL; 
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG poly_destroy completed\n");
	
	return(0);
}

int poly_set_to_zero(Poly* p)
{

	if ( p == NULL ) return(-1);

	p->n = 0;
	poly_init_bounding_box(p);
	
	return(0);
}

Poly* poly_copy(Poly* pi)
{
	if (pi==NULL) return(NULL);
	
	Poly* po = poly_init(pi->nt);
	for ( int i =0; i<pi->n; i++)
		poly_add_coord(pi->vertex[i], po);
		
	return(po);
}

int poly_bin(Poly* p, int bin)
{
	if (p==NULL) return(-1);
	if (bin==1) return(0);
	
	coord pc;	
	int n = p->n;
	p->n = 0;
	poly_init_bounding_box(p);
	for ( int i =0; i<n; i++) {
		pc[X] = (p->vertex[i][X] + bin - 1) / bin;
		pc[Y] = (p->vertex[i][Y] + bin - 1) / bin;
		poly_add_coord(pc, p);
	}

//	poly_create_convex_hull(p);	

	return(0);
}
	
float poly_area(Poly* p)
{

	if (p->n == 0) return(0.);
	
	float parea = 0.;

	if (p->n>2) {

		// 2 A(P) = sum_{i=0}^{n-1} ( x_i  (y_{i+1} - y_{i-1}) )

		// i=0
		parea = p->vertex[0][X] * ( p->vertex[1][Y] - p->vertex[p->n-1][Y] );

		// i=1:n-2
         for ( int i = 1; i<p->n-2; i++)
			// parea = parea + Pol(x,k)*(Pol(y,k+1)-Pol(y,k-1))
			parea += p->vertex[i][X] * ( p->vertex[i+1][Y] - p->vertex[i-1][Y] );

		//i=n-1
		parea += p->vertex[p->n-1][X] * ( p->vertex[0][Y] - p->vertex[p->n-2][Y]);
                 
		parea = 0.5 * fabs(parea);

	} else if (p->n == 1) {
		// degenerate polygon: one vertex
        parea = 1.0;

	} else if (p->n == 2) {
		// degenerate polygon: a segment
		// parea = SQRT((Pol(x,1)-Pol(x,2))**2+(Pol(y,1)-Pol(y,2))**2)
		parea = sqrt((float)((p->vertex[0][X]-p->vertex[1][X])*(p->vertex[0][X]-p->vertex[1][X]) +
					(p->vertex[0][Y]-p->vertex[1][Y])*(p->vertex[0][Y]-p->vertex[1][Y])));
	}

	return(parea);
}
