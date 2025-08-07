/*
	symmetry.h
	Header file for general symmetry functions
	Author: Bernard Heymann
	Created: 20010420 	Modified: 20041118
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include "matrix.h"

#ifndef _Bsymmetry_
/************************************************************************
@Object: struct UnitCell
@Description:
	Unit cell parameters.
@Features:
	The convention is as follows:
		a lies on the x-axis, 
		b lies in the x,y plane. 
*************************************************************************/
struct UnitCell {
	float   a, b, c;		// Major axis lengths
	float   alf, bet, gam;	// Angles between major axes
} ;

/************************************************************************
@Object: struct Bsymop
@Description:
	General symmetry operator.
@Features:
	Both point groups and space groups are covered.
	A symmetry operator describes the complete series of rotations around 
		an axis, with the number of operations given by the order field
		and the rotation angle given by the angle field. The shift field
		is used in gliding operations, such as needed for helical symmetry.
*************************************************************************/
struct Bsymop {
	Vector3	axis;		// Symmetry axis
	int 	order;		// Number of times the operation is applied, 1=identity
	float	angle;		// Rotation angle, derived from order when 0
	float	shift;		// Shift up symmetry axis for helix
} ;

/************************************************************************
@Object: struct Bsymmetry
@Description:
	General symmetry structure.
@Features:
	Both point groups and space groups are covered.
*************************************************************************/
struct Bsymmetry {
	int 	point;		// Point group (< 1 if a crystal)
	int 	space;		// Space group (< 1 if not a crystal)
	char	label[80];	// Label
	int 	n;			// Number of symmetry operators
	Bsymop* op; 		// Symmetry operators
} ;
#define _Bsymmetry_
#endif

// Function prototypes
Bsymmetry* 	init_point_group_symmetry(char* symmetry_string);
int 	 	kill_symmetry(Bsymmetry* sym);
int			symmetry_clean_pointgroup(char* symmetry_string);
View*		symmetry_get_all_views(Bsymmetry* sym, View asu_view);
View* 		asymmetric_unit_views(Bsymmetry* sym, float theta_step, float phi_step, int full);
int 		change_views_to_asymmetric_unit(Bsymmetry* sym, View* view);
int 		find_asymmetric_unit_view(Bsymmetry* sym, View* theview);
float*  	calc_skew_matrix(UnitCell unit_cell);
float*		sym_matrices_from_text_list(int nsym, char* symop, int line_len);
int			sym_show_matrices(Bsymmetry* sym);

