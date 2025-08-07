/*
	proj.h
	Header file for projection structure handling
	Author: Giovanni Cardone
	Created: 20041215 	Modified:
*/

#ifndef _PROJ_H__
#define _PROJ_H__

#include "bsoft.h"
#include "tilt.h"
#include "shape.h"
#include "norm.h"


/************************************************************************
@Object: struct Proj
@Description:
	Projection image parameter structure.
@Features:
	Encapsulate array of Bimage, one for each projections, and adds some
	parameters needed for tomography
	All floating point coordinates are in angstroms and converted 
	using the voxel units parameters.
*************************************************************************/
struct Proj {
	Bimage**		img;				// images list (pointers array)
	char			filename[256];	// File name
	unsigned int		dataflag;	// Flag to force reading of the data
	unsigned long	x, y;		// Dimensions, xy
	unsigned long	nmax;		// Number of maximum images allowed
	unsigned long	nimg;		// Number of actual images allocated
	DataType			datatype;	// Data type
	float			resolution; 	// Resolution limit of data
								// - used for low-pass filtering
	float			ox, oy, oz; 	// Origin (position of tilt axis origin
								// in voxel coordinates)
	float			vx, vy, vz; 	// View orientation of z-axis,
								// i.e. axis perpendicular to
								// untilted projection plane
	float			angle;		// Rotation around view vector in radians
	float			ux, uy;		// Voxel units (angstrom/pixel edge)
	float			avg;			// overall projections average
	float			std;			// overall projections standard deviation
	float			min;			// overall projections standard deviation
	float			max;			// overall projections standard deviation
	float *			tilt;		// tilt angles
	float			xtilt;		// x-axis tilt
	Shape2d **		footprint;	// footprint of selected area (pointers array)
} ;

// Function prototypes
Proj* 	proj_init(int np);
int 		proj_kill(Proj* p);
int 	proj_init_footprint(Proj* p, Vector3 origin, Shape3d* s, int taper);
int 	proj_copy_footprint(Proj* pout, Proj* pin);
int  	proj_cpy_img_sett(Proj* p, Bimage* b);
int  	proj_puttilt(Proj* p, float tilt_angle, unsigned int idx);
int  	proj_putview(Proj* p, View v, unsigned int idx);
View  	proj_getview(Proj* p, unsigned int idx);
int  	proj_putnprojs(Proj* p, int np);
int 		proj_display(Proj* p);
int  	proj_putimage(Proj* p, Bimage* b, unsigned int idx, int pin);
int  	proj_deleteimage(Proj* p, unsigned int idx);
int  	proj_write_angles(char* tfile, Proj* p, int ntlt);
int 		proj_stats(Proj* p, int bkg_type);
int 		proj_mask(Proj* p);
int	 	proj_to_datatype(Proj* p, DataType dt);
int 		proj_rescale_to_avg_std(Proj* p, float avg, float std, float* background);
int		proj_rescale_linregression(Proj* p, Proj* pref);

#endif  // #ifndef _PROJ_H__

