/*
	shape3d.h
	Header file for normalizing functions
	Author: Giovanni Cardone
	Created: 20050113 	Modified: 
*/

#ifndef _NORM_H__
#define _NORM_H__

#include <math.h>
#include "bsoft.h"
#include "shape.h"


// Function prototypes
int norm3d_rescale_to_avg_std(Bimage* p, float avg, float std, float* background, Shape3d* mask);
int	norm3d_rescale_linregression(Bimage* p, Bimage* pref, Shape3d* mask);
int norm2d_rescale_to_avg_std(Bimage* p, float avg, float std, float* background, Shape2d* mask);
int	norm2d_rescale_linregression(Bimage* p, Bimage* pref, Shape2d* mask);


#endif  // #ifndef _NORM_H__


