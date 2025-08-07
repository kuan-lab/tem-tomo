/*
	wedge.cc
	Functions for handling missing wedge
	Author: Giovanni Cardone
	Created: 20050116 	Modified: 20050602
*/

#include "wedge.h"

/************************************************************************
@Function: wedge_init
@Description:
	initialize a wedge structure.
@Algorithm:
@Arguments:
	float    tmin	minimum tilt angle
	float    tmax	maximum tilt angle
	float    xt		x-tilt angle
	float    xyrot	tilt rotation on xy plane
@Returns:
	Wedge*         wedge initialized
**************************************************************************/
Wedge* wedge_init( float tmin, float tmax, float xt, float xyrot)
{

	if (tmin == tmax) return(NULL);
	
	Wedge*  w = (Wedge *) balloc(sizeof(Wedge));

	if (tmin > tmax) {
		float ftmp = tmax;
		tmax = tmin;
		tmin = ftmp;
	}
	w->tilt_min = tmin / 180. * M_PI;
	w->tilt_max = tmax / 180. * M_PI;
	
	w->xtilt = xt / 180. * M_PI;
	
	w->xyrot = 0.;
	
	if (xyrot!=0.) {
		fprintf(stderr,"Warning(wedge_init): in-plane rotation not yet implemented:");
		fprintf(stderr," xyrot set to zero!\n");
	}
		
	return(w);
}

/************************************************************************
@Function: wedge_kill
@Description:
	General wedge structure destruction.
@Algorithm:
	This function deallocates all memory associated to the structure
@Arguments:
	Wedge*	w			wedge structure
@Returns:
	int					error code (<0 means failure)
**************************************************************************/
int 	wedge_kill(Wedge* w)
{
	if ( w == NULL ) return(0);
	
	bfree(w, sizeof(Wedge));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG wedge_kill memory = %ld\n", memory);
	
	return(0);
}

/************************************************************************
@Function: wedge_point_inside
@Description:
	verify if a voxel is inside a wedge.
	Voxel coordinates are given in the reciprocal space
@Algorithm:
@Arguments:
	Wedge*  w		   reference wedge.
	float   x, y, z     point coords (1/A)
@Returns:
	int					>0 if test verified
**************************************************************************/
int wedge_point_inside(Wedge* w, float x, float y, float z) {

	int is_inside = 0;
	if(x==0. && y==0. && z==0.) return(is_inside);
	
	float tangle = 0.;
	if (x==0.) {
		if (z==0.) tangle = 0.;
		else tangle = 0.5*M_PI;
	}
	else		tangle = atan((z*cos(w->xtilt)+y*sin(w->xtilt))/x);

	is_inside = tangle > w->tilt_max || tangle < w->tilt_min;

	return(is_inside);	
}

/************************************************************************
@Function: wedge_init_from_imod
@Description:
	initialize a wedge structure from a imod parameter file.
@Algorithm:
@Arguments:
	Imod_tilt * im_tilt	imod parameter file
@Returns:
	Wedge*         wedge initialized
**************************************************************************/
Wedge* wedge_init_from_imod(Imod_tilt * im_tilt)
{

	Wedge*  w = (Wedge *) balloc(sizeof(Wedge));

	int nt = 0;
	float* tangles =	tlt_load_tilt_angles(imd_tilt_gettiltfile(im_tilt), &nt);

	w->tilt_min = 0.;
	w->tilt_max = 0.;
	for ( int i = 0; i<nt; i++) {
		if (tangles[i] > w->tilt_max)
			w->tilt_max = tangles[i];
		else if (tangles[i] < w->tilt_min)
			w->tilt_min = tangles[i];
	}
	
	w->xtilt = imd_tilt_getxtilt(im_tilt);
	
	w->xyrot = 0.;

	bfree(tangles, nt*sizeof(float));
		
	return(w);
}

/************************************************************************
@Function: wedge_mask_fourier_vol
@Description:
	mask volume in the missing wedge region.
	Volume needs to be already in the fourier domain
@Algorithm:
@Arguments:
	Bimage*  b	volume to mask
	Wedge*   w  missing wedge
@Returns:
	int          >0 if successful
**************************************************************************/
int wedge_mask_fourier_vol(Bimage* b, Wedge* w)
{

	if ( b->transform != Standard ) {
		fprintf(stderr, "Error: map from %s must be given in the fourier domain!\n", b->filename); 
		return(-1); 
	}

	float		freq_scale[3], rx, ry, rz;
	
	// The frequency scaling is linked to the different dimensions of the
	// data set (important when the x, y and z dimensions are different)
	// The radius scaling is set to a value consistent with one pixel width
	// in reciprocal space for the longest dimension.
	freq_scale[0] = 1.0/(b->x*b->ux);
	freq_scale[1] = 1.0/(b->y*b->uy);
	freq_scale[2] = 1.0/(b->z*b->uz);
	
	complex_float*	data = (complex_float *) b->data;
	
	int n, z, zz, x, xx, y, yy;
	unsigned long i;
		
	for ( n=0; n< (int) b->n; n++ ) { 
		for ( z=0; z< (int) b->z; z++ ) { 
			zz = z; 
			if ( z > (int) (b->z - 1)/2 ) zz -= b->z;
			rz = zz*freq_scale[2];
			for ( y=0; y< (int) b->y; y++ ) { 
				yy = y; 
				if ( y > (int) (b->y - 1)/2 ) yy -= b->y;
				ry = yy*freq_scale[1];
				for ( x=0; x< (int) b->x; x++ ) {
					i = ((n*b->z + z)*b->y + y)*b->x + x;
					xx = x; 
					if ( x > (int) (b->x - 1)/2 ) xx -= b->x;
					rx = xx*freq_scale[0];
					if(wedge_point_inside(w, rx, ry, rz)) {
						data[i].re = 0.;
						data[i].im = 0.;
					}
				}
			} 
		}
	}

	return(0);
}
