/*
	views.h
	Header file for general view functions
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

// Function prototypes
View* 		tilt_views(float ang_min, float ang_max, float ang_step);
View*		random_views(int nviews);


