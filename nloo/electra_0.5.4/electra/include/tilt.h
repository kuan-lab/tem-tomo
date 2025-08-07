/*
	tilt.h
	Header file to handle tilt functions
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004
*/

#ifndef _TILT_H__
#define _TILT_H__

#include "bsoft.h"

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <ctype.h>
//#include <math.h>
//#include <time.h>
//#include <unistd.h>

// Constants
const int   MAXNANGLES=500;  /* max number of tilt angles allowed */ 
	
// Function prototypes
float*		tlt_load_tilt_angles(char* filename, int* ntlt_angs);
float*		tlt_create_series_from_range(float amin, float amax, float astep, int* ntlt_angs);
int			tlt_save_angles_to_file(float* tlt_angs, int ntlt_angs, char* filename);
View*		tlt_views_from_angles(int nviews, View tlt0_view, float* tlt_angs);

#endif  // #ifndef _TILT_H__
