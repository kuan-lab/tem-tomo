/*
	img_background.h
	Header file for manipulating image backgrounds
	Author: Bernard Heymann
	Created: 20000430 	Modified: 20040112
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include "rwimg.h"

// Function prototypes
int			img_calculate_background(Bimage* p);
int			img_correct_background(Bimage* p);
int			img_subtract_background(Bimage* p);
int			img_set_background(Bimage* p, float background);

