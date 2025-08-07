/*
	img_project.h
	Header file for projection functions
	Author: Bernard Heymann
	Created: 20010420 	Modified: 20041108
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "rwimg.h"
#include "matrix.h"

// Function prototypes
int 		img_project(Bimage* p);
Bimage*		img_rotate_project(Bimage* p, Matrix3 mat, Vector3 translate,
					float radial_cutoff);
Bimage* 	img_project_views(Bimage* p, View* view);
int			img_project_write_views(Bimage* p, View* view, char* proj_file);
int 		img_compare_projections(Bimage* p);

