/*
	img_edit.h
	Header file for functions for editing image contents
	Author: Bernard Heymann
	Created: 19980520 	Modified: 20040927
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "rwimg.h"
#include "matrix.h"

// Function prototypes
int			img_delete_within_box(Bimage* p, Vector3 origin, VectorInt3 box, 
				float width, float fill);
int			img_delete_within_radius(Bimage* p, Vector3 origin, float radius, 
				float width, float fill);
int			img_delete_within_cylinder(Bimage* p, Vector3 origin, float radius, 
				float height, float width, float fill);
int			img_delete_range(Bimage* p, float min, float max, float fill);
int			img_edge_oval(Bimage* p, VectorInt3 size, Vector3 origin, float width, float fill);
int			img_edge_cylinder(Bimage* p, VectorInt3 size, Vector3 origin, float width, float fill);
int			img_edge_rect(Bimage* p, VectorInt3 size, Vector3 origin, float width, float fill);
int 		img_gradient_correction(Bimage* p);
float		img_shell_average(Bimage* p, Vector3 origin, float radius, float width);
int  		img_extract_tetrahedron(Bimage* p, Vector3* vec, float fill);
