/*
	shape.h
	Header file for 2D/3D shape definitions
	Author: Giovanni Cardone
	Created: 20050112 	Modified: 
*/

#ifndef _SHAPE_H__
#define _SHAPE_H__

#include "bsoft.h"
#include "poly.h"
//#include "img_average.h"

/*----------Shape Structure-------------*/
/************************************************************************
@Object: struct shape3d
@Description:
	General volume shape structure.
@Features:
	This contains all the information pertinent to a general volume shape.
*************************************************************************/
struct Shape3d {
	int			type;			// 0=entire volume 1=box 2=sphere
	Vector3		start;
	Vector3		size;
	int          taper;
	float        tapsize;
	Shape3d*     inner;
} ;

/************************************************************************
@Object: struct shape2d
@Description:
	General area shape structure.
@Features:
	This contains all the information pertinent to a general area shape.
*************************************************************************/
struct Shape2d {
	int			type;			// 0=entire area 1=polygon 2=circle
	Poly*		pol;
	Vector3		center;
	float		radius;
	int          taper;
	float        tapsize;
} ;

// constants
const int MAX_NPVERT = 12;
const float SHAPE_TAP_WINDOW_PERC = 0.15;
const int VALL = 0, VBOX = 1, VSPHERE = 2, VMAXSPHERE = 3, VNONE = 0; 
const int VOUTSIDE = 0, VFRAME = 1, VINSIDE = 2; 
const int SALL = 0, SPOLY = 1, SCIRCLE = 2; 
const int SOUTSIDE = 0, SFRAME = 1, SINSIDE = 2, SNONE = -1; 

// Function prototypes
Shape3d* shape3d_init( Shape3d* s, int type, float ox, float oy, float oz,
					  float sx, float sy, float sz, int inner); 
Shape3d* shape3d_maxinit( Shape3d* s, Bimage* b, int type);
int 		shape3d_kill(Shape3d* s);
int 		shape3d_mask_image( Bimage* p, Shape3d* mask, int taper);
int 		shape3d_taper_image( Bimage* p, Shape3d* mask);
int 		shape3d_point_inside(Shape3d* mask, int ix, int iy, int iz);
float 	shape3d_point_distance(Shape3d* mask, int ix, int iy, int iz);
int		shape3d_image_stats(Bimage* p, Shape3d* mask, int bground_type);
int 		shape3d_taper_init(Shape3d* mask, int taper);
float	shape3d_fill_ratio(Bimage* p, Shape3d* mask);
Shape3d* shape3d_copy( Shape3d* s);
int		shape3d_bin( Shape3d* s, int bin);
int		shape3d_check(Bimage* p, Shape3d* mask);
Shape2d* shape_proj3dto2d(Shape3d* mask, int px, int py, Vector3 origin, View view);
int 		shape2d_taper_init(Shape2d* mask, int taper);
Shape2d* shape2d_init_to_box(int ox, int oy, int sizex, int sizey);
int 		shape2d_kill(Shape2d* s);
int 		shape2d_point_inside(Shape2d* mask, int ix, int iy);
int		shape2d_point_inframe(Shape2d* mask, int ix, int iy);
float 	shape2d_point_distance(Shape2d* mask, int ix, int iy);
int		shape2d_image_stats(Bimage* p, Shape2d* mask, int bground_type);
int		shape2d_mask_image(Bimage* p, Shape2d* mask, int taper);
Shape2d* shape2d_copy( Shape2d* s);
int		shape2d_bin( Shape2d* s, int bin);
float		shape2d_fill_ratio(Bimage* p, Shape2d* mask);

#endif  // #ifndef _SHAPE_H__


