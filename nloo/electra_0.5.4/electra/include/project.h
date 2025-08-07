/*
	tilt.h
	Header file to handle projection functions
	Author: Giovanni Cardone
	Created: 20050119 	Modified: 20050608
*/

#ifndef _PROJECTION_H__
#define _PROJECTION_H__

#include "bsoft.h"
#include "shape.h"
#include "proj.h"
#include "ioroutines.h"

// Global Variables

// Function prototypes
Proj* 	prjct_generate_one(Bimage* p, int nv, View* view, Shape3d* volut, int taper,
							int bground_type, int exp_tag, Buff* bf);
int		prjct_rotate_and_project(Proj* proj, Bimage* p, int ip, Matrix3 rotat_mat,
							Vector3 translat_vec, Shape3d* volut, int exp_tag, Buff* bf);
Vector3 	prjct_locate_point(Vector3 p, Vector3 origin, View view);

#endif  // #ifndef _PROJECTION_H__
