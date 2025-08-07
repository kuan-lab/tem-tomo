/*
	wedge.h
	Header file for missing wedge handling
	Author: Giovanni Cardone
	Created: 20050116 	Modified:
*/

#ifndef _WEDGE_H__
#define _WEDGE_H__

#include "bsoft.h"
#include "imod.h"


/************************************************************************
@Object: struct Wedge
@Description:
	Wedge parameters structure.
	Note: routine to handle xyrot not yet implemented.
@Features:
	Encapsulates all the angles needed to define a missing wedge:
	minimum and maximum tilt angle,
	rotation around the x axis,
	rotation of the tilt axis (y) on the xy plane.
*************************************************************************/
struct Wedge {
	float			tilt_min; 	// minimum tilt angle
	float			tilt_max; 	// maximum tilt angle
	float			xtilt;		// x-axis tilt
	float			xyrot;		// rotation on xy plane
};

// Function prototypes
Wedge* 	wedge_init(float tmin, float tmax, float xt, float xyrot);
int 	wedge_kill(Wedge* w);
int wedge_point_inside(Wedge* w, float x, float y, float z);
Wedge* wedge_init_from_imod(Imod_tilt * im_tilt);
int wedge_mask_fourier_vol(Bimage* b, Wedge* w);

#endif  // #ifndef _WEDGE_H__

