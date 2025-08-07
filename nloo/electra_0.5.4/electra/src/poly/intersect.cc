/*
	intersect.cc
	Functions for determining
	the intersection between convex polygons
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004

	All the functions are adapted from existing code.
*/
/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 7.

Written by Joseph O'Rourke.
Last modified: December 1997
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
redistributed in its entirety provided that this copyright notice is
not removed.
--------------------------------------------------------------------
*/

#include "intersect.h"

/************************************************************************
@Function: poly_intersect
@Description:
	Determines the intersection between two convex polygons
	as a new convex polygon.
@Algorithm:
    Walks around the boundaries of the two polygons in concert, intersecting
    edges.
@Arguments:
	Poly* P				first polygon to compare
	Poly* Q				second polygon to compare
@Returns:
	Poly*				polygon intersection
**************************************************************************/
int poly_intersect(Poly *P, Poly *Q, Poly *Pol)
{
	int     a, b;           /* indices on P and Q (resp.) */
	int     a1, b1;         /* a-1, b-1 (resp.) */
	tPointi A;           /* directed edges on P and Q (resp.) */
	tPointi B;           /* directed edges on P and Q (resp.) */
	int     cross;          /* sign of z-component of A x B */
	int     bHA, aHB;       /* b in H(A); a in H(b). */
	tPointi Origin = {0,0}; /* (0,0) */
	tPointd p;              /* double point of intersection */
	tPointi pc;             /* point of intersection transformed to coord*/
	tPointd q;              /* second point of intersection */
	tInFlag inflag;         /* {Pin, Qin, Unknown}: which inside */
	int     aa, ba;         /* # advances on a & b indices (after 1st inter.) */
	bool    FirstPoint;     /* Is this the first point? (used to initialize).*/
//	tPointd p0;             /* The first point. */
	int     code;           /* SegSegInt return code. */ 

//	Poly*   Pol = poly_init(P->n+Q->n);
	if ( Pol->n > 0 )   poly_set_to_zero(Pol);
	
   /* Initialize variables. */
   a = 0; b = 0; aa = 0; ba = 0;
   inflag = Unknown; FirstPoint = true;

   do {
	  if (verbose & VERB_FULL) {
		printf("Before Advances:a=%d, b=%d; aa=%d, ba=%d; inflag=%d\n", a, b, aa, ba, inflag);
		poly_print(Pol);
	  }
      /* Computations of key variables. */
      a1 = (a + P->n - 1) % P->n;
      b1 = (b + Q->n - 1) % Q->n;

      SubVec( P->vertex[a], P->vertex[a1], A );
      SubVec( Q->vertex[b], Q->vertex[b1], B );

      cross = AreaSign( Origin, A, B );
      aHB   = AreaSign( Q->vertex[b1], Q->vertex[b], P->vertex[a] );
      bHA   = AreaSign( P->vertex[a1], P->vertex[a], Q->vertex[b] );

	  if (verbose & VERB_FULL)
		printf("cross=%d, aHB=%d, bHA=%d\n", cross, aHB, bHA );

      /* If A & B intersect, update inflag. */
      code = SegSegInt( P->vertex[a1], P->vertex[a], Q->vertex[b1], Q->vertex[b], p, q );
	  if (verbose & VERB_FULL)
		printf("SegSegInt: code = %c\n", code );
 
      if ( code == '1' || code == 'v' ) {
         if ( inflag == Unknown && FirstPoint ) {
            aa = ba = 0;
            FirstPoint = false;
//            p0[X] = p[X]; p0[Y] = p[Y];
//            printf("%8.2lf %8.2lf moveto\n", p0[X], p0[Y] );
         }
		 pc[X] = (int) p[X];
		 pc[Y] = (int) p[Y];
		 poly_add_coord( pc, Pol );
         inflag = InOut( p, inflag, aHB, bHA );
		 if (verbose & VERB_FULL) {
			printf("new point added: %d\t%d %d\n", Pol->n, Pol->vertex[Pol->n-1][X], Pol->vertex[Pol->n-1][Y]);
			printf("InOut sets inflag=%d\n", inflag);
		 }
	  }

      /*-----Advance rules-----*/

      /* Special case: A & B overlap and oppositely oriented. */
      if ( ( code == 'e' ) && (Dot( A, B ) < 0) ) {
			pc[X] = (int) q[X];
			pc[Y] = (int) q[Y];
			poly_add_coord( pc, Pol );
			if (verbose & VERB_FULL) {
				printf("new point added: %d\t%d %d\n", Pol->n, Pol->vertex[Pol->n-1][X], Pol->vertex[Pol->n-1][Y]);
				printf("Completed: A & B overlap and oppositely oriented\n");
			}
//			exit(EXIT_SUCCESS);
//			return(Pol);
			return(0);

	  }

      /* Special case: A & B parallel and separated. */
      if ( (cross == 0) && ( aHB < 0) && ( bHA < 0 ) ) {
			Pol->n = 0;
			if (verbose & VERB_FULL)
				printf("P and Q are disjoint.\n");
//			exit(EXIT_SUCCESS);
//			return(Pol);
			return(0);
	  }
      /* Special case: A & B collinear. */
      else if ( (cross == 0) && ( aHB == 0) && ( bHA == 0 ) ) {
            /* Advance but do not output point. */
            if ( inflag == Pin )
               b = Advance( b, &ba, Q->n, inflag == Qin, Q->vertex[b], Pol );
            else
               a = Advance( a, &aa, P->n, inflag == Pin, P->vertex[a], Pol );
      }
      /* Generic cases. */
      else if ( cross >= 0 ) {
         if ( bHA > 0)
            a = Advance( a, &aa, P->n, inflag == Pin, P->vertex[a], Pol );
         else
            b = Advance( b, &ba, Q->n, inflag == Qin, Q->vertex[b], Pol );
      }
      else /* if ( cross < 0 ) */{
         if ( aHB > 0)
            b = Advance( b, &ba, Q->n, inflag == Qin, Q->vertex[b], Pol );
         else
            a = Advance( a, &aa, P->n, inflag == Pin, P->vertex[a], Pol );
      }
      if (verbose & VERB_FULL)
		printf("After advances:a=%d, b=%d; aa=%d, ba=%d; inflag=%d\n", a, b, aa, ba, inflag);

   /* Quit when both adv. indices have cycled, or one has cycled twice. */
   } while ( ((aa < P->n) || (ba < Q->n)) && (aa < 2*P->n) && (ba < 2*Q->n) );

//   if ( !FirstPoint ) /* If at least one point output, close up. */
//            printf("%8.2lf %8.2lf lineto\n", p0[X], p0[Y] );

   /* Deal with special cases */
   if ( Pol->n > 0 ) {
   /* verify if first and last point coincidents */
		if ( Pol->vertex[0][X] == Pol->vertex[Pol->n-1][X] &&
			Pol->vertex[0][Y] == Pol->vertex[Pol->n-1][Y] ) {
			Pol->n -=1;
			if (verbose & VERB_FULL)
				printf("FULL: Eliminated last point because redundant.\n");
		}
	} else {
	/* test which polygon is inside the other */
		pc[X] = Q->vertex[0][X];
		pc[Y] = Q->vertex[0][Y];
		if (poly_point_inside(P, pc)) {
			for (int i=0; i<Q->n; i++)
//				poly_add_point(Q->vertex[i][X],Q->vertex[i][Y],Pol);
				poly_add_coord(Q->vertex[i],Pol);
			if (verbose & VERB_FULL)
				printf("FULL: Intersection coincident with second polygon.\n");
		} else {
			pc[X] = P->vertex[0][X];
			pc[Y] = P->vertex[0][Y];
			if(poly_point_inside(Q, pc)){
				for (int i=0; i<P->n; i++)
//					poly_add_point(P->vertex[i][X],P->vertex[i][Y],Pol);
					poly_add_coord(P->vertex[i],Pol);
				if (verbose & VERB_FULL)
					printf("FULL: Intersection coincident with second polygon.\n");
			}
		}
	}		
/*   if ( inflag == Unknown)
      if (verbose & VERB_FULL)
		printf("The boundaries of P and Q do not cross.\n");
*/
//	return(Pol);
	return(0);
}

/*---------------------------------------------------------------------
Prints out the double point of intersection, and toggles in/out flag.
---------------------------------------------------------------------*/
tInFlag InOut( tPointd p, tInFlag inflag, int aHB, int bHA )
{
//   printf("%8.2lf %8.2lf lineto\n", p[X], p[Y] );

   /* Update inflag. */
   if      ( aHB > 0)
      return Pin;
   else if ( bHA > 0)
      return Qin;
   else    /* Keep status quo. */
      return inflag;
}
/*---------------------------------------------------------------------
   Advances and prints out an inside vertex if appropriate.
---------------------------------------------------------------------*/
int     Advance( int a, int *aa, int n, bool inside, tPointi v, Poly* P )
{
   if ( inside ) {
	  poly_add_coord( v, P);
	  if (verbose & VERB_FULL)
		printf("new point added: %d\t%d %d\n", P->n, P->vertex[P->n-1][X], P->vertex[P->n-1][Y]);
   }
   (*aa)++;
   return  (a+1) % n;
}

/*---------------------------------------------------------------------
a - b ==> c.
---------------------------------------------------------------------*/
void    SubVec( tPointi a, tPointi b, tPointi c )
{
   int i;

   for( i = 0; i < DIM; i++ )
      c[i] = a[i] - b[i];
}

/*---------------------------------------------------------------------
SegSegInt: Finds the point of intersection p between two closed
segments ab and cd.  Returns p and a char with the following meaning:
   'e': The segments collinearly overlap, sharing a point.
   'v': An endpoint (vertex) of one segment is on the other segment,
        but 'e' doesn't hold.
   '1': The segments intersect properly (i.e., they share a point and
        neither 'v' nor 'e' holds).
   '0': The segments do not intersect (i.e., they share no points).
Note that two collinear segments that share just one point, an endpoint
of each, returns 'e' rather than 'v' as one might expect.
---------------------------------------------------------------------*/
char	SegSegInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p, tPointd q )
{
   double  s, t;       /* The two parameters of the parametric eqns. */
   double num, denom;  /* Numerator and denoninator of equations. */
   char code = '?';    /* Return char characterizing intersection. */

	if ( verbose & VERB_FULL )
		printf("SegSegInt: a,b,c,d: (%d,%d), (%d,%d), (%d,%d), (%d,%d)\n",
			a[X],a[Y], b[X],b[Y], c[X],c[Y], d[X],d[Y]);

   denom = a[X] * (double)( d[Y] - c[Y] ) +
           b[X] * (double)( c[Y] - d[Y] ) +
           d[X] * (double)( b[Y] - a[Y] ) +
           c[X] * (double)( a[Y] - b[Y] );

   /* If denom is zero, then segments are parallel: handle separately. */
   if (denom == 0.0)
      return  ParallelInt(a, b, c, d, p, q);

   num =    a[X] * (double)( d[Y] - c[Y] ) +
            c[X] * (double)( a[Y] - d[Y] ) +
            d[X] * (double)( c[Y] - a[Y] );
   if ( (num == 0.0) || (num == denom) ) code = 'v';
   s = num / denom;
   /*printf("num=%lf, denom=%lf, s=%lf\n", num, denom, s);*/

   num = -( a[X] * (double)( c[Y] - b[Y] ) +
            b[X] * (double)( a[Y] - c[Y] ) +
            c[X] * (double)( b[Y] - a[Y] ) );
   if ( (num == 0.0) || (num == denom) ) code = 'v';
   t = num / denom;
   /*printf("num=%lf, denom=%lf, t=%lf\n", num, denom, t);*/

   if      ( (0.0 < s) && (s < 1.0) &&
             (0.0 < t) && (t < 1.0) )
     code = '1';
   else if ( (0.0 > s) || (s > 1.0) ||
             (0.0 > t) || (t > 1.0) )
     code = '0';

   p[X] = a[X] + s * ( b[X] - a[X] );
   p[Y] = a[Y] + s * ( b[Y] - a[Y] );

   return code;
}

char   ParallelInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p, tPointd q )
{
/*   
   printf("ParallelInt: a,b,c,d: (%d,%d), (%d,%d), (%d,%d), (%d,%d)\n",
	a[X],a[Y], b[X],b[Y], c[X],c[Y], d[X],d[Y]);
*/

   if ( !Collinear( a, b, c) )
      return '0';

   if ( Between( a, b, c ) && Between( a, b, d ) ) {
      Assigndi( p, c );
      Assigndi( q, d );
      return 'e';
   }
   if ( Between( c, d, a ) && Between( c, d, b ) ) {
      Assigndi( p, a );
      Assigndi( q, b );
      return 'e';
   }
   if ( Between( a, b, c ) && Between( c, d, b ) ) {
      Assigndi( p, c );
      Assigndi( q, b );
      return 'e';
   }
   if ( Between( a, b, c ) && Between( c, d, a ) ) {
      Assigndi( p, c );
      Assigndi( q, a );
      return 'e';
   }
   if ( Between( a, b, d ) && Between( c, d, b ) ) {
      Assigndi( p, d );
      Assigndi( q, b );
      return 'e';
   }
   if ( Between( a, b, d ) && Between( c, d, a ) ) {
      Assigndi( p, d );
      Assigndi( q, a );
      return 'e';
   }
   return '0';
}

void	Assigndi( tPointd p, tPointi a )
{
   int i;
   for ( i = 0; i < DIM; i++ )
      p[i] = (double) a[i];
}
/*---------------------------------------------------------------------
Returns TRUE iff point c lies on the closed segement ab.
Assumes it is already known that abc are collinear.
---------------------------------------------------------------------*/
bool    Between( tPointi a, tPointi b, tPointi c )
{

   /* If ab not vertical, check betweenness on x; else on y. */
   if ( a[X] != b[X] )
      return ((a[X] <= c[X]) && (c[X] <= b[X])) ||
             ((a[X] >= c[X]) && (c[X] >= b[X]));
   else
      return ((a[Y] <= c[Y]) && (c[Y] <= b[Y])) ||
             ((a[Y] >= c[Y]) && (c[Y] >= b[Y]));
}

/*---------------------------------------------------------------------
Returns the dot product of the two input vectors.
---------------------------------------------------------------------*/
double  Dot( tPointi a, tPointi b )
{
    int i;
    double sum = 0.0;

    for( i = 0; i < DIM; i++ )
       sum += a[i] * b[i];

    return  sum;
}

