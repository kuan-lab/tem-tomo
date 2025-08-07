/*
	img_resize.h
	Header file for image resizing
	Author: Bernard Heymann
	Created: 20000430 	Modified: 20030927
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
int 		img_resize(Bimage* p, VectorInt3 newsize, VectorInt3 translate, int wrap, float fill);
int 		img_shrink_wrap(Bimage* p, VectorInt3 newsize, VectorInt3 translate);
int 		img_resize_to_next_power2(Bimage* p, float fill);
int 		img_pad(Bimage* p, int size, float fill);
