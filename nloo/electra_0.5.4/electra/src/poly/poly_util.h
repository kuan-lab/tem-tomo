/*
	poly_util.h
	Header file for coordinates standard
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004
*/

#ifndef _POLYUTIL_H__
#define _POLYUTIL_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bsoft.h"
#include "poly_def.h"

#define BBHUGE 32768

/*----------Function Prototypes-------------*/
void cp_coord(int* csrc, int* cdest);
coord* double2coord(double* dc); 
Poly* poly_init(int n);
int poly_destroy(Poly* p);
int poly_set_to_zero(Poly* p);
void poly_print( Poly *P );
int poly_add_coord(coord p, Poly* Pol);
int poly_add_point(int px, int py, Poly* Pol);
int poly_init_bounding_box(Poly* P);
int poly_update_bounding_box(coord p, Poly* P);
int	AreaSign( coord a, coord b, coord c );
bool LeftOn( coord a, coord b, coord c );
bool Left( coord a, coord b, coord c );
bool Collinear( coord a, coord b, coord c );
Poly* poly_copy(Poly* pi);
float poly_area(Poly* p);

#endif  // #ifndef _POLYUTIL_H__
