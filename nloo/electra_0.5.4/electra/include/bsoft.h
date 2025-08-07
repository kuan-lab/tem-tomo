/*
	bsoft.h
	Header file for linking bsoft library
	Author: Giovanni Cardone
	Created: 20040323   Modified: 20040323
*/
	
#ifndef _BSOFT_H__
#define _BSOFT_H__


#include "img_combine.h"
#include "img_complex.h"
#include "img_datatypes.h"
#include "img_extract.h"
#include "img_background.h"
#include "img_util.h"
#include "img_interpolate.h"
#include "img_fourier.h"
#include "img_pick.h"
#include "img_project.h"
#include "img_random.h"
#include "img_rescale.h"
#include "img_resize.h"
#include "img_resolution.h"
#include "matrix.h"
#include "matrix_linear.h"
#include "options.h"
#include "random_numbers.h" 
#include "spline.h"
#include "string_util.h"
#include "symmetry.h"
#include "timer.h"
#include "utilities.h"
#include "rwimg.h"
#include "rwMRC.h"
#include "rwCCP4.h"

#ifndef M_PI
#define M_PI    3.14159265358979323846264338327950288
#endif  // #ifndef M_PI

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

#endif  // #ifndef _BSOFT_H__
