/*
	img_arithmetic.h
	Header file for doing simple image arithmetic
	Author: Bernard Heymann
	Created: 20040727 	Modified: 20040915
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
int			img_add_constant(Bimage* p, float constant);
int			img_multiply_with_constant(Bimage* p, float constant);

