/*
	shape.cc
	Functions for handling masking by shapes
	Author: Giovanni Cardone
	Created: 20050112 	Modified: 20050628
*/

#include "shape.h"

/************************************************************************
@Function: shape3d_init
@Description:
	initialize a volume structure.
@Algorithm:
	Struct is allocated if needed. According to the value of inner, external
	or internal surface of the volume is initialized. 
@Arguments:
	Shape3d* s		input shape (NULL if to be allocated)
	int      type    volume type
	float    ox		size values (6)
	float    oy           their meaning depends on the value of the type argument
	float    oz
	float    sx
	float    sy
	float    sz
	int      inner	flag for selecting external (0) or internal (1) surface
@Returns:
	Shape3d*         shape initialized
**************************************************************************/
Shape3d* shape3d_init( Shape3d* s, int type, float ox, float oy, float oz,
					  float sx, float sy, float sz, int inner)
{
	Shape3d*  sh = NULL;
	Shape3d*  sp = NULL;
	
	if ( type != VBOX && type != VSPHERE ) {
		fprintf(stderr,"Error(shape3d_init): volume type not allowed!\n");
		return(NULL);
	}

	if (s==NULL) {
		if (inner==1) {
			fprintf(stderr,"Error(shape3d_init): can not initialize inner shape of volume not allocated!\n");
			return(NULL);
		}
		sh = (Shape3d *) balloc(sizeof(Shape3d));	
		sh->inner = NULL;
	} else
		sh = s;

	if (inner==1) {
		if(sh->inner==NULL) {
			sh->inner = (Shape3d *) balloc(sizeof(Shape3d));	
			(sh->inner)->inner =NULL;
		}
		sp = sh->inner;
	} else
		sp = sh;

	sp->type = type;
	switch ( type ) {
		case VBOX:
			sp->start = vector3_from_3_values(ox, oy, oz);
			sp->size = vector3_from_3_values(sx, sy, sz);
		break;
		case VSPHERE:
			sp->start = vector3_from_3_values(ox, oy, oz);
			sp->size = vector3_from_3_values(sx, 0., 0.);
		break;
	}

	sp->taper = 0;
	sp->tapsize = 0.;
	
	return(sh);
}

/************************************************************************
@Function: shape3d_maxinit
@Description:
	initialize a volume structure to maximum extension.
@Algorithm:
	Maxima dimensions needed to determine the size of the volume
	are extracted from the reference image. 
@Arguments:
	Shape3d* s		input shape (NULL if to be allocated)
	Bimage*  b       reference image
	int      type    volume type
@Returns:
	Shape3d*         shape initialized
**************************************************************************/
Shape3d* shape3d_maxinit( Shape3d* s, Bimage* b, int type)
{

	Shape3d*  sh = NULL;
	
	if ( type != VALL && type != VMAXSPHERE ) {
		fprintf(stderr,"Error(shape3d_maxinit): volume type not allowed!\n");
		return(NULL);
	}
		
	if (s==NULL) {
		sh = (Shape3d *) balloc(sizeof(Shape3d));	
		sh->inner = NULL;
	} else
		sh = s;

	switch ( type ) {
		case VALL:
			sh->start = vector3_from_3_values(0, 0, 0);
			sh->size = vector3_from_3_values(b->x, b->y, b->z);
		break;	
		case VMAXSPHERE:
			sh->start = vector3_from_3_values(0.5*(b->x-1.0), 0.5*(b->y-1.0), 0.5*(b->z-1.0));
			unsigned long min_diameter = b->x;
			if (b->y < min_diameter )
				min_diameter = b->y;
			if (b->z < min_diameter )
				min_diameter = b->z;
			sh->size = vector3_from_3_values(0.5*min_diameter-1., 0., 0.);
		break;
	}

	if (type==VALL) sh->type = VALL;
	else if (type==VMAXSPHERE) sh->type = VSPHERE;

	sh->taper = 0;
	sh->tapsize = 0.;
	
	return(sh);
}

/************************************************************************
@Function: shape3d_kill
@Description:
	General shape3d structure destruction.
@Algorithm:
	This function deallocates all memory associated to the structure
@Arguments:
	Shape3d*				projection structure
@Returns:
	int					error code (<0 means failure)
**************************************************************************/
int 	shape3d_kill(Shape3d* s)
{
	if ( s == NULL ) return(0);
	
	if (s->inner != NULL)
		bfree(s->inner, sizeof(Shape3d));
			
	bfree(s, sizeof(Shape3d));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG shape3d_kill memory = %ld\n", memory);
	
	return(0);
}

/************************************************************************
@Function: shape3d_mask_image
@Description:
	mask an image according to given shape.
@Algorithm:
	The background value is applied to the volume outside the given shape.
	If requested, a tapering is applied.
@Arguments:
	Bimage*  p		input image - modified.
	Shape3d* mask	mask to apply.
	int      taper   taper flag
@Returns:
	int				<0 if not successfull.
**************************************************************************/
int		shape3d_mask_image(Bimage* p, Shape3d* mask, int taper)
{
	if ( !p->data ) return(-1);
	if ( mask->type==VALL && taper==0 ) return(1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    		fdata = (float *) p->data;
    complex_short* 	csdata = (complex_short *) p->data;
    complex_int* 	cidata = (complex_int *) p->data;
    complex_float*  	cfdata = (complex_float *) p->data;
    polar*  			pdata = (polar *) p->data;
    
	float	background, shp_dist;
    unsigned int  n, x, y, z;
    unsigned long i;
	float   shp_tapsize = mask->tapsize;
	
	for ( n = 0; n < p->n; n++) {
	
		background = p->image[n].background;
		
// mask with background value, and taper if needed
		for ( z=0; z<p->z; z++ ) {
			for ( y=0; y<p->y; y++ ) {
				for ( x=0; x<p->x; x++ ) {
					if ( !shape3d_point_inside(mask, x, y, z) ) {
						i = ((n*p->z+z)*p->y+y)*p->x+x;
						switch ( p->datatype ) {
							case UChar:
								udata[i] = (unsigned char) background;
								break;
							case SChar:
								cdata[i] = (signed char) background;
								break;
							case UShort:
								usdata[i] = (unsigned short) background;
								break;
							case Short:
								sdata[i] = (short) background;
								break;
							case Int:
								idata[i] = (int) background;
								break;
							case Float:
								fdata[i] = background;
								break;
							case ComplexShort:
								csdata[i].re = (short) (sqrt(0.5)*background);
								csdata[i].im = (short) (sqrt(0.5)*background);
								break;
							case ComplexInt:
								cidata[i].re = (int) (sqrt(0.5)*background);
								cidata[i].im = (int) (sqrt(0.5)*background);
								break;
							case ComplexFloat:
								cfdata[i].re = (sqrt(0.5)*background);
								cfdata[i].im = (sqrt(0.5)*background);
								break;
							case Polar:
								pdata[i].amp = background;
								break;
							default: break;
						}
					} else if (taper == 1) {
						shp_dist = shape3d_point_distance(mask,x,y,z);
						if (shp_dist < shp_tapsize) {
							i = ((n*p->z+z)*p->y+y)*p->x+x;
							switch ( p->datatype ) {
							case UChar:
								udata[i] = (unsigned char) ((udata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case SChar:
								cdata[i] = (signed char) ((cdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case UShort:
								usdata[i] = (unsigned short) ((usdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case Short:
								sdata[i] = (short) ((sdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case Int:
								idata[i] = (int) ((idata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case Float:
								fdata[i] = ((fdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case ComplexShort:
								csdata[i].re = (short) (sqrt(0.5)*((csdata[i].re - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								csdata[i].im = (short) (sqrt(0.5)*((csdata[i].im - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								break;
							case ComplexInt:
								cidata[i].re = (int) (sqrt(0.5)*((cidata[i].re - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								cidata[i].im = (int) (sqrt(0.5)*((cidata[i].im - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								break;
							case ComplexFloat:
								cfdata[i].re = (sqrt(0.5)*((cfdata[i].re - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								cfdata[i].im = (sqrt(0.5)*((cfdata[i].im - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								break;
							case Polar:
								pdata[i].amp = ((pdata[i].amp - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							default: break;
							}
						}
					}
				}
			}
		}

		if ( verbose & VERB_FULL )
			printf("\n Image # %d masked with background value: %f \n", n, background);
	}
	
	return(1);
}

/************************************************************************
@Function: shape3d_taper_image
@Description:
	taper an image around the given shape.
@Algorithm:
@Arguments:
	Bimage*  p			input image - modified.
	Shape3d* mask		reference mask.
@Returns:
	int					<0 if not successfull.
**************************************************************************/
int		shape3d_taper_image(Bimage* p, Shape3d* mask)
{
	if ( !p->data || mask == NULL) return(-1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    		fdata = (float *) p->data;
    complex_short* 	csdata = (complex_short *) p->data;
    complex_int* 	cidata = (complex_int *) p->data;
    complex_float*  	cfdata = (complex_float *) p->data;
    polar*  			pdata = (polar *) p->data;
    
	float		background, shp_tapsize, shp_dist;
    unsigned int   n, x, y, z;
    unsigned long  i;

	shp_tapsize = mask->tapsize;

	for ( n = 0; n < p->n; n++) {
	
		background = p->image[n].background;
		
		for ( z=0; z<p->z; z++ ) {
			for ( y=0; y<p->y; y++ ) {
				for ( x=0; x<p->x; x++ ) {
					if ( shape3d_point_inside(mask, x, y, z)) {
						shp_dist = shape3d_point_distance(mask,x,y,z);
						if ( shp_dist < shp_tapsize ) {
							i = ((n*p->z+z)*p->y+y)*p->x+x;
							switch ( p->datatype ) {
							case UChar:
								udata[i] = (unsigned char) ((udata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case SChar:
								cdata[i] = (signed char) ((cdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case UShort:
								usdata[i] = (unsigned short) ((usdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case Short:
								sdata[i] = (short) ((sdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case Int:
								idata[i] = (int) ((idata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case Float:
								fdata[i] = ((fdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case ComplexShort:
								csdata[i].re = (short) (sqrt(0.5)*((csdata[i].re - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								csdata[i].im = (short) (sqrt(0.5)*((csdata[i].im - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								break;
							case ComplexInt:
								cidata[i].re = (int) (sqrt(0.5)*((cidata[i].re - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								cidata[i].im = (int) (sqrt(0.5)*((cidata[i].im - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								break;
							case ComplexFloat:
								cfdata[i].re = (sqrt(0.5)*((cfdata[i].re - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								cfdata[i].im = (sqrt(0.5)*((cfdata[i].im - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								break;
							case Polar:
								pdata[i].amp = ((pdata[i].amp - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							default: break;
							}
						}
					}
				}
			}
		}

		if ( verbose & VERB_FULL )
			printf("\nTapering applied\n");
	}
	
	return(1);
}

/************************************************************************
@Function: shape3d_point_inside
@Description:
	verify if a voxel is inside a mask.
@Algorithm:
@Arguments:
	Shape3d*  mask		  reference mask.
	int       ix, iy, iz  point coords
@Returns:
	int					>0 if test verified
**************************************************************************/
int shape3d_point_inside(Shape3d* mask, int ix, int iy, int iz) {

	if (mask==NULL) return(1);
	
	int is_inside = 0;
	float x = ix;
	float y = iy;
	float z = iz;

	if ( mask->type == VALL ) is_inside = 1;
	else {
		switch ( mask->type ) {
		case VBOX:
			if ( x>=mask->start.x && x<=(mask->start.x+mask->size.x-1.) &&
				 y>=mask->start.y && y<=(mask->start.y+mask->size.y-1.) &&
				 z>=mask->start.z && z<=(mask->start.z+mask->size.z-1.) )
					is_inside = 1;
		break;
		case VSPHERE:
			if ( (x-mask->start.x)*(x-mask->start.x) +
				 (y-mask->start.y)*(y-mask->start.y) +
				 (z-mask->start.z)*(z-mask->start.z) <=
				  mask->size.x*mask->size.x) 
					is_inside = 1;
		break;
		}
	}

	if ( is_inside == 1 && mask->inner != NULL && (mask->inner)->type != 0 ) {
		Shape3d * m = mask->inner;
		switch ( m->type ) {
		case VBOX:
			if ( x>=m->start.x && x<=(m->start.x+m->size.x-1.) &&
				  y>=m->start.y && y<=(m->start.y+m->size.y-1.) &&
				   z>=m->start.z && z<=(m->start.z+m->size.z-1.) ) 
					is_inside = 0;
			break;
		case VSPHERE:
			if ( (x-m->start.x)*(x-m->start.x) +
				  (y-m->start.y)*(y-m->start.y) +
				   (z-m->start.z)*(z-m->start.z) <= m->size.x*m->size.x)
					is_inside = 0;
			break;
		}
	}

	return(is_inside);	
}

/************************************************************************
@Function: shape3d_point_inframe
@Description:
	verify if a voxel is in a frame around the mask surface(s).
@Algorithm:
@Arguments:
	Shape3d*  mask		  reference mask.
	int       ix, iy, iz  point coords
@Returns:
	int					>0 if test verified
**************************************************************************/
int shape3d_point_inframe(Shape3d* mask, int ix, int iy, int iz) {

	return(shape3d_point_distance(mask,ix,iy,iz) < (int) mask->tapsize);

}

/************************************************************************
@Function: shape3d_point_distance
@Description:
	calculate the distance of a point from the border of the ROI
@Algorithm:
	it evaluate the minimum distance from each plane/sphere.
	Warning: it works only if point inside the ROI
@Arguments:
	Shape   mask		reference mask.
	int     ix, iy, iz  point coords
@Returns:
	float				distance value
**************************************************************************/
float shape3d_point_distance(Shape3d* mask, int ix, int iy, int iz) {

	float distance=-1., idistance=-1.;
	Vector3 vdist;
	float x = ix;
	float y = iy;
	float z = iz;
	
	switch ( mask->type ) {
	case VALL: case VBOX:
		distance = fabs(mask->start.x+mask->size.x-1. - x);
		if ( fabs(mask->start.x - x) < distance )
			distance = fabs(mask->start.x - x);
		if ( fabs(mask->start.y - y) < distance )
			distance = fabs(mask->start.y - y);
		if ( fabs(mask->start.y+mask->size.y-1. - y) < distance )
			distance = fabs(mask->start.y+mask->size.y-1. - y);
		if ( fabs(mask->start.z - z) < distance )
			distance = fabs(mask->start.z - z);
		if ( fabs(mask->start.z+mask->size.z-1. - z) < distance )
			distance = fabs(mask->start.z+mask->size.z-1. - z);
		break;
	case VSPHERE:
		distance = fabs(mask->size.x - sqrt(
			(x-mask->start.x)*(x-mask->start.x) +
			(y-mask->start.y)*(y-mask->start.y) +
			(z-mask->start.z)*(z-mask->start.z)));
		break;
	}

	if ( mask->inner != NULL ) {
		Shape3d * m = mask->inner;
		switch ( m->type ) {
		case VBOX:
			vdist.x = fabs(m->start.x+m->size.x-1. - x);
			if ( fabs(m->start.x - x) < vdist.x ) {
				vdist.x = fabs(m->start.x - x);
			}
			if ((y>=m->start.y && y<=(m->start.y+m->size.y-1.)) &&
			    (z>=m->start.z || z<=(m->start.z+m->size.z-1.)))
				vdist.y = vdist.z = 0.;
			else {
				vdist.y = fabs(m->start.y+m->size.y-1. - y);
				if ( fabs(m->start.y - y) < vdist.y ) {
					vdist.y = fabs(m->start.y - y);
				}
				if ((x>=m->start.x && x<=(m->start.x+m->size.x-1.)) &&
				    (z>=m->start.z || z<=(m->start.z+m->size.z-1.)))
					vdist.x = vdist.z = 0.;
				else {
					vdist.z = fabs(m->start.z+m->size.z-1. - z);
					if ( fabs(m->start.z - z) < vdist.z ) {
						vdist.z = fabs(m->start.z - z);
					}
					if ((x>=m->start.x && x<=(m->start.x+m->size.x-1.)) &&
					    (y>=m->start.y || y<=(m->start.y+m->size.y-1.)))
						vdist.x = vdist.y = 0.;
				}
			}
			idistance = sqrt(vdist.x*vdist.x+vdist.y*vdist.y+vdist.z*vdist.z);
		break;
		case VSPHERE:
			idistance = fabs(m->size.x - sqrt(
				(x-m->start.x)*(x-m->start.x) +
				(y-m->start.y)*(y-m->start.y) +
				(z-m->start.z)*(z-m->start.z)));
		break;
		}
		if (idistance < distance) distance = idistance; 
	}
	
	return(distance);	
}

/************************************************************************
@Function: shape3d_image_stats
@Description:
	evaluate statistics on image
@Algorithm:
	The background value, average and standard deviation
	are calculates in the region defined by the mask.
	Average and standard deviation refer to all the images,
	if a multi-image format provided
@Arguments:
	Bimage*  p				input image - stats modified.
	Shape3d* mask			image mask.
	int 		bground_type		background evaluation method
@Returns:
	int					<0 if not successfull.
**************************************************************************/
int		shape3d_image_stats(Bimage* p, Shape3d* mask, int bground_type)
{
	if ( !p->data ) return(-1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* usdata = (unsigned short *) p->data;
    short* 	    	sdata = (short *) p->data;
    int* 	    	idata = (int *) p->data;
    float*  	    fdata = (float *) p->data;
    complex_short* 	csdata = (complex_short *) p->data;
    complex_int* 	cidata = (complex_int *) p->data;
    complex_float*  cfdata = (complex_float *) p->data;
    polar*  		pdata = (polar *) p->data;
    
	double		sum, sum2, background;
	float		amp;
    long     	i, number, bgnumber;
	unsigned int x, y, z;
	int			isinside;


	sum = sum2 = 0.;
	number = 0;
	
	for ( unsigned int n = 0; n < p->n; n++) {
	
		bgnumber = 0;
		background = 0.;
		for ( z=0; z<p->z; z++ ) {
			for ( y=0; y<p->y; y++ ) {
				for ( x=0; x<p->x; x++ ) {
					i = ((n*p->z+z)*p->y+y)*p->x+x;
					isinside = shape3d_point_inside(mask, x, y, z);
					// average and standard deviation
					if ( isinside ) {
						switch ( p->datatype ) {
							case UChar:
								sum += udata[i];
								sum2 += udata[i]*1.0*udata[i];
								break;
							case SChar:
								sum += cdata[i];
								sum2 += cdata[i]*1.0*cdata[i];
								break;
							case UShort:
								sum += usdata[i];
								sum2 += usdata[i]*1.0*usdata[i];
								break;
							case Short:
								sum += sdata[i];
								sum2 += sdata[i]*1.0*sdata[i];
								break;
							case Int:
								sum += idata[i];
								sum2 += idata[i]*1.0*idata[i];
								break;
							case Float:
								sum += fdata[i];
								sum2 += fdata[i]*fdata[i];
								break;
							case ComplexShort:
								amp = sqrt(1.0*csdata[i].re*csdata[i].re +
    	    								csdata[i].im*csdata[i].im);
								sum += amp;
								sum2 += amp*amp;
								break;
							case ComplexInt:
								amp = sqrt(1.0*cidata[i].re*cidata[i].re +
    	    								cidata[i].im*cidata[i].im);
								sum += amp;
								sum2 += amp*amp;
								break;
							case ComplexFloat:
								amp = sqrt(cfdata[i].re*cfdata[i].re +
    	    								cfdata[i].im*cfdata[i].im);
								sum += amp;
								sum2 += amp*amp;
								break;
							case Polar:
								sum += pdata[i].amp;
								sum2 += pdata[i].amp*pdata[i].amp;
								break;
							default: break;
						}
						number++;
					}

					// background
					if ( (bground_type==VOUTSIDE && !isinside) || 
					   (bground_type==VINSIDE && isinside) ||
					   (isinside && bground_type==VFRAME &&
					   shape3d_point_inframe(mask, x, y, z))) {
						switch ( p->datatype ) {
							case UChar:
								background += udata[i];
								break;
							case SChar:
								background += cdata[i];
								break;
							case UShort:
								background += usdata[i];
								break;
							case Short:
								background += sdata[i];
								break;
							case Int:
								background += idata[i];
								break;
							case Float:
								background += fdata[i];
								break;
							case ComplexShort:
								background += sqrt(1.0*csdata[i].re*csdata[i].re +
    	    								csdata[i].im*csdata[i].im);
								break;
							case ComplexInt:
								background += sqrt(1.0*cidata[i].re*cidata[i].re +
    	    								cidata[i].im*cidata[i].im);
								break;
							case ComplexFloat:
								background += sqrt(cfdata[i].re*cfdata[i].re +
    	    								cfdata[i].im*cfdata[i].im);
								break;
							case Polar:
								background += pdata[i].amp;
								break;
							default: break;
						}
						bgnumber++;
					}
				}
			}
		}
		if ( bgnumber ) background = background/bgnumber;
		p->image[n].background = background;
	}

	if (number) {
		p->avg = sum/number;
		p->std = sum2/number - p->avg*p->avg;
		if ( p->std > 0 ) p->std = sqrt(p->std);
		else p->std = 0.;
	}

	return(0);
}

/************************************************************************
@Function: shape3d_taper_init
@Description:
	Set tapering and evaluate the size of tapering.
	The size of the tapering window is always calculated.
	The shape is set to be tapered, depending on the value of 'taper'
@Algorithm:
	The window is dimensioned as a percentile (value fixed into the code) of
	the bounding box the mask.
	The dimension with minimum size is selected.
@Arguments:
	Shape*   mask       volume under test
	int      taper      tapering flag
@Returns:
	int					<0 if not successfull.
**************************************************************************/
int shape3d_taper_init(Shape3d* mask, int taper) {

	Shape3d*    sh;
	if (mask==NULL) return(-1);

	mask->taper = taper;
	if (mask->inner!=NULL) (mask->inner)->taper = taper;

	float tapsize = 0.;
	float in_max = 0.;
	Vector3 in = {0.,0.,0.};

	if (mask->inner != NULL) {
		if (verbose & VERB_DEBUG)
			printf("DEBUG(shape3d_taper_init): dimensioning inner volume\n"); 
		sh = mask->inner;
		switch (sh->type) {
			case VBOX:
				in.x = sh->size.x;
				in.y = sh->size.y;
				in.z = sh->size.z;
			break;
			case VSPHERE:
				in.x = 2.0*sh->size.x;
				in.y = 2.0*sh->size.x;
				in.z = 2.0*sh->size.x;
			break;
		}
		if (in.x >= in.y && in.x >= in.z)
			in_max = in.x;
		else if (in.y >= in.x && in.y >= in.z)
			in_max = in.y;
		else if (in.z >= in.x && in.z >= in.y)
			in_max = in.z;
	}

	sh = mask;
	switch (sh->type) {
		case VALL: case VBOX:
			if (verbose & VERB_DEBUG)
				printf("DEBUG(shape3d_taper_init): box dimensioning\n"); 
			if( (sh->size.x-in.x) > 1 && (sh->size.x-in.x) <= (sh->size.y-in.y) &&
			    (sh->size.x-in.x) <= (sh->size.z-in.z))
				tapsize = (0.5*(sh->size.x-in.x)*SHAPE_TAP_WINDOW_PERC);
			else if ( (sh->size.y-in.y) > 1 &&
			          (sh->size.y-in.y) <= (sh->size.x-in.x) &&
		    		      (sh->size.y-in.y) <= (sh->size.z-in.z))
				tapsize = (0.5*(sh->size.y-in.y)*SHAPE_TAP_WINDOW_PERC);
			else if ( (sh->size.z-in.z) > 1 &&
					 (sh->size.z-in.z) <= (sh->size.x-in.x) &&
					 (sh->size.z-in.z) <= (sh->size.y-in.y))
				tapsize = (0.5*(sh->size.z-in.z)*SHAPE_TAP_WINDOW_PERC);
		break;
		case VSPHERE:
			if (verbose & VERB_DEBUG)
				printf("DEBUG(shape3d_taper_init): spher dimensioning\n"); 
			tapsize = (0.5*(2.0*sh->size.x - in_max)*SHAPE_TAP_WINDOW_PERC);
		break;
	}

	mask->tapsize = tapsize;
	if(mask->inner!=NULL) (mask->inner)->tapsize = tapsize;

	return(0);
}

/************************************************************************
@Function: shape3d_fill_ratio
@Description:
	return the fill ratio of the shape with respect to the given volume.
@Algorithm:
@Arguments:
	Bimage*  p		input image - modified.
	Shape3d* mask	reference mask.
@Returns:
	float			fill ratio.
**************************************************************************/
float		shape3d_fill_ratio(Bimage* p, Shape3d* mask)
{
	
	if (p==NULL || mask == NULL) return(0.);
	
	float pvol = 1.;
	if (p->x>1) pvol*= p->x;
	if (p->y>1) pvol*= p->y;
	if (p->z>1) pvol*= p->z;
	
	float svol = 1.;
	switch (mask->type) {
		case VALL:
			svol = pvol;
		break;	
		case VBOX:
			if ( mask->size.x>1 ) svol*= mask->size.x;
			if ( mask->size.y>1 ) svol*= mask->size.y;
			if ( mask->size.z>1 ) svol*= mask->size.z;
		break;
		case VSPHERE:
			// area of cube containing the sphere
			svol = 8.0*mask->size.x*mask->size.x*mask->size.x;
		break;
	}
	
	return(svol/pvol);
	
}

/************************************************************************
@Function: shape3d_copy
@Description:
	copy a volume structure to a new variable.
@Algorithm:
@Arguments:
	Shape3d* s		input shape
@Returns:
	Shape3d*         copy of input shape
**************************************************************************/
Shape3d* shape3d_copy( Shape3d* s) {

	Shape3d*		sh = (Shape3d *) balloc(sizeof(Shape3d));	
	
	sh->inner = NULL;
	
	Shape3d* ps = s;
	Shape3d* psh = sh;
	
	for ( int i = 0; i<2; i++) {
		if (psh==NULL) continue;
		psh->type = ps->type;
		psh->start.x = ps->start.x;
		psh->start.y = ps->start.y;
		psh->start.z = ps->start.z;
		psh->size.x = ps->size.x;
		psh->size.y = ps->size.y;
		psh->size.z = ps->size.z;
		psh->taper = ps->taper;
		psh->tapsize = ps->tapsize;
		if(ps->inner!=NULL) sh->inner = (Shape3d *) balloc(sizeof(Shape3d));
		psh = sh->inner;
		ps = ps->inner;
	}

	return(sh);
}

/************************************************************************
@Function: shape3d_bin
@Description:
	Bins by an integer size.
@Algorithm:
@Arguments:
	Shape3d*   s		input shape
	int		  bin	bin factor.
@Returns:
	int 			0.
**************************************************************************/
int shape3d_bin( Shape3d* s, int bin) {

	if (s==NULL) return(-1);
	
	Shape3d* sh = s;
	
	sh->start.x = (sh->start.x + bin - 1) / bin;
	sh->start.y = (sh->start.y + bin - 1) / bin;
	sh->start.z = (sh->start.z + bin - 1) / bin;
	sh->size.x = (sh->size.x + bin - 1) / bin;;
	sh->size.y = (sh->size.y + bin - 1) / bin;
	sh->size.z = (sh->size.z + bin - 1) / bin;

	if(s->inner!=NULL) {
		sh = s->inner;
		sh->start.x = (sh->start.x + bin - 1) / bin;
		sh->start.y = (sh->start.y + bin - 1) / bin;
		sh->start.z = (sh->start.z + bin - 1) / bin;
		sh->size.x = (sh->size.x + bin - 1) / bin;;
		sh->size.y = (sh->size.y + bin - 1) / bin;
		sh->size.z = (sh->size.z + bin - 1) / bin;
	}

	shape3d_taper_init(sh, sh->taper);

	return(0);
}

/************************************************************************
@Function: shape3d_check
@Description:
	Verify if shape properly contained into the given image volume.
@Algorithm:
@Arguments:
	Bimage*  p		input image.
	Shape3d* mask	reference mask.
@Returns:
	int					error code (<0 means failure)
**************************************************************************/
int		shape3d_check(Bimage* p, Shape3d* mask)
{
	
	if (p==NULL || mask == NULL) return(-1);
	
	int ic = 1;

	switch ( mask->type ) {
		case VALL: case VBOX:
			ic = ic && (mask->start.x >= 0);
			ic = ic && (mask->start.y >= 0);
			ic = ic && (mask->start.z >= 0);
			ic = ic && (mask->start.x+mask->size.x-1. < p->x);
			ic = ic && (mask->start.y+mask->size.y-1. < p->y);
			ic = ic && (mask->start.z+mask->size.z-1. < p->z);
			if (ic!=1)
				fprintf(stderr,"Warning: selected box out of image size!\n");
		break;
		case VSPHERE:
			ic = ic && (mask->start.x-mask->size.x+1. >= 0);
			ic = ic && (mask->start.y-mask->size.x+1. >= 0);
			ic = ic && (mask->start.z-mask->size.x+1. >= 0);
			ic = ic && (mask->start.x+mask->size.x-1. <p->x);
			ic = ic && (mask->start.y+mask->size.x-1. <p->y);
			ic = ic && (mask->start.z+mask->size.x-1. <p->z);
			if (ic!=1)
				fprintf(stderr,"Warning: selected sphere out of image size!\n");
		break;
	}
	
	return(ic==1);
}

/************************************************************************
@Function: shape_proj3dto2d
@Description:
	Project a 3d shape along a given view
@Algorithm:
@Arguments:
	Shape3d*		mask		input 3D shape
	int			px, py	maximum 2D box size allowed
						(values read only if the mask is a box) 
	Vector3 		origin	origin coordinates
	View			view		view vector
@Returns:
	Shape2d*    			projected area
**************************************************************************/
Shape2d* shape_proj3dto2d(Shape3d* mask, int px, int py, Vector3 origin, View view)
{
	
	Vector3      new_origin = { 0.5*(px-1), 0.5*(py-1), 0.};
	int 			x, y, z, xo, yo;
	double 		oldx, oldy, oldz;
	double		newx, newy;
	VectorInt3  box_origin, box_size;

	Poly * pol_acq = NULL;
	Poly * pol_prj = NULL;
	Shape2d* 	s = (Shape2d*) balloc(sizeof(Shape2d));
	Matrix3		rotat_mat;

	if ( verbose & VERB_FULL ) {
		printf("Initialing footprint for projection\n");
		printf("with view (%f,%f,%f,%f)\n", view.x, view.y, view.z, view.a*180.0/M_PI);
	}

	rotat_mat = matrix3_from_view(view);
	// The rotation matrix is transposed because the orientation parameters
	// give the orientation as rotated from the reference view
	rotat_mat = matrix3_transpose(rotat_mat);

	switch ( mask->type ) {
		case VALL: case VBOX:
			s->type = SPOLY;
			box_origin = vectorint3_from_vector3(mask->start);
			box_size = vectorint3_from_vector3(mask->size);

			s->pol = poly_init(MAX_NPVERT);

			// initialize polygons describing projection and acquisition area
			pol_acq = poly_init(MAX_NPVERT);

			// loop over vertices of 2D acquisition  area
			for ( y=0; y<py; y+=py-1 )
				for ( x=0; x<px; x+=px-1 )
					poly_add_point(x, y, pol_acq);
							
			if ( verbose & VERB_FULL ) {
				printf("Loaded points for acquisition polygon\n");
				poly_print(pol_acq);
			}
			poly_create_convex_hull(pol_acq);

			pol_prj = poly_init(MAX_NPVERT);
	
			// loop over vertices of 3D volume
			for ( z=box_origin.z; z<box_origin.z+box_size.z; z+=box_size.z-1 ) {
				oldz = z - origin.z;
				for ( y=box_origin.y; y<box_origin.y+box_size.y; y+=box_size.y-1 ) {
					oldy = y - origin.y;
					for ( x=box_origin.x; x<box_origin.x+box_size.x; x+=box_size.x-1 ) {
						oldx = x - origin.x;
						newy = oldx*rotat_mat.r10 + oldy*rotat_mat.r11 +
								oldz*rotat_mat.r12 + new_origin.y;
						yo = (int) newy;
						newx = oldx*rotat_mat.r00 + oldy*rotat_mat.r01 +
								oldz*rotat_mat.r02 + new_origin.x;
						xo = (int) newx;
						poly_add_point(xo, yo, pol_prj);
					}
				}
			}
			if ( verbose & VERB_FULL ) {
				printf("Loaded points for projection polygon\n");
				poly_print(pol_prj);
			}
			poly_create_convex_hull(pol_prj);

			poly_intersect(pol_acq, pol_prj, s->pol);

			if ( verbose & VERB_FULL ) {
				printf("Footprint with %d vertices initialized\n", (s->pol)->n);
				poly_print(s->pol);
			}

			bfree(pol_prj,sizeof(Poly));
			bfree(pol_acq,sizeof(Poly));
		break;
		case VSPHERE:
			s->type = SCIRCLE;
			oldz = mask->start.z - origin.z;
			oldy = mask->start.y - origin.y;
			oldx = mask->start.x - origin.x;
			s->center.y = oldx*rotat_mat.r10 + oldy*rotat_mat.r11 + 
					oldz*rotat_mat.r12 + new_origin.y;
			s->center.x = oldx*rotat_mat.r00 + oldy*rotat_mat.r01 + 
					oldz*rotat_mat.r02 + new_origin.x;
			s->center.z = 0.;
			s->radius = mask->size.x;
		break;
	}

	return(s);
}

/************************************************************************
@Function: shape2d_taper_init
@Description:
	Set tapering and evaluate the size of tapering.
	The size of the tapering window is always calculated.
	The shape is set to be tapered, depending on the value of 'taper'
@Algorithm:
	The window is dimensioned as a percentile (value fixed into the code) of
	the bounding box of the polygon or as a percentile of the circle radius.
	The dimension with minimum size is selected.
@Arguments:
	Shape2d* mask       area under test
	int      taper      tapering flag
@Returns:
	int					<0 if not successfull.
**************************************************************************/
int shape2d_taper_init(Shape2d* mask, int taper) {

	if (mask==NULL) return(-1);

	float tapsize = 0.;
	int xrange, yrange;
	Shape2d* sh = mask;
	
	sh->taper = taper;

	switch (sh->type) {
		case SALL: case SPOLY:
			if (verbose & VERB_DEBUG)
				printf("DEBUG(shape2d_taper_init): box dimensioning\n"); 
			xrange = ((sh->pol)->bbxmax-(sh->pol)->bbxmin);
			yrange = ((sh->pol)->bbymax-(sh->pol)->bbymin);
	
			if( xrange < yrange)
				tapsize = (0.5*xrange*SHAPE_TAP_WINDOW_PERC);
			else
				tapsize = (0.5*yrange*SHAPE_TAP_WINDOW_PERC);
		break;
		case SCIRCLE:
			if (verbose & VERB_DEBUG)
				printf("DEBUG(shape2d_taper_init): circle dimensioning\n"); 
			tapsize =  (sh->radius*SHAPE_TAP_WINDOW_PERC);
		break;
	}

	sh->tapsize = tapsize;

	return(0);
}

/************************************************************************
@Function: shape2d_init_to_box
@Description:
	Initialize a 2D to all the area covered by a given map
	3rd dimension of the map is ignored
@Algorithm:
@Arguments:
	int 		ox, oy       lower coordinates of the box
	int 		sizex, sizey box size
@Returns:
	Shape2d*    			projected area
**************************************************************************/
Shape2d* shape2d_init_to_box(int ox, int oy, int sizex, int sizey)
{
	Shape2d* 	s = (Shape2d*) balloc(sizeof(Shape2d));

	s->type = SALL;
	s->pol = poly_init(MAX_NPVERT);

	// loop over vertices of 2D box
	for ( int y=oy; y<oy+sizey; y+=sizey-1 )
		for ( int x=ox; x<ox+sizex; x+=sizex-1 )
			poly_add_point(x, y, s->pol);
					
	if ( verbose & VERB_FULL ) {
		printf("Loaded points for defining a box polygon\n");
		poly_print(s->pol);
	}
	poly_create_convex_hull(s->pol);

	return(s);
}				

/************************************************************************
@Function: shape2d_kill
@Description:
	General shape2d structure destruction.
@Algorithm:
	This function deallocates all memory associated to the structure
@Arguments:
	Shape2d*				projection structure
@Returns:
	int					error code (<0 means failure)
**************************************************************************/
int 	shape2d_kill(Shape2d* s)
{
	if ( s == NULL ) return(0);
	
	if (s->pol != NULL)
		poly_destroy(s->pol);
			
	bfree(s, sizeof(Shape2d));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG shape2d_kill memory = %ld\n", memory);
	
	return(0);
}

/************************************************************************
@Function: shape2d_point_inside
@Description:
	verify if a voxel is inside a mask.
@Algorithm:
@Arguments:
	Shape2d*  mask		 reference mask.
	int       ix, iy, iz  point coords
@Returns:
	int					>0 if test verified
**************************************************************************/
int shape2d_point_inside(Shape2d* mask, int ix, int iy)
{

	if (mask==NULL) return(1);

	int is_inside = 0;

	if ( mask->type == SALL ) is_inside = 1;
	else {
		switch ( mask->type ) {
		case SPOLY:
		{
			coord pnt;
			pnt[X] = ix;
			pnt[Y] = iy;
			is_inside = poly_point_inside(mask->pol, pnt);
		}
		break;
		case SCIRCLE:
		{
			float dd = mask->radius - sqrt((ix-mask->center.x)*(ix-mask->center.x) +
								 (iy-mask->center.y)*(iy-mask->center.y));
			if (dd>=0.) is_inside = 1;
		}
		break;
		}
	}

	return(is_inside);	
}

/************************************************************************
@Function: shape2d_point_inframe
@Description:
	verify if a voxel is in a frame around the mask surface(s).
@Algorithm:
@Arguments:
	Shape2d*  mask		  reference mask.
	int       ix, iy		  point coords
@Returns:
	int					>0 if test verified
**************************************************************************/
int shape2d_point_inframe(Shape2d* mask, int ix, int iy)
{
	return(shape2d_point_distance(mask,ix,iy) < (int) mask->tapsize);
}

/************************************************************************
@Function: shape2d_point_distance
@Description:
	calculate the distance of a point from the border of the ROI
@Algorithm:
	it evaluate the minimum distance from each line/circle.
	Warning: it works only if point inside the ROI
@Arguments:
	Shape2d*   mask		reference mask.
	int     ix, iy		point coords
@Returns:
	float				distance value
**************************************************************************/
float shape2d_point_distance(Shape2d* mask, int ix, int iy)
{

	float distance = 0.;
	
	switch ( mask->type ) {
		case SALL: case SPOLY:
		{
			coord pnt;
			pnt[X] = ix;
			pnt[Y] = iy;
			distance = poly_distance_to_polygon(mask->pol, pnt);
		}
		break;
		case SCIRCLE:
			distance = fabs(mask->radius - sqrt((ix-mask->center.x)*(ix-mask->center.x) +
					 (iy-mask->center.y)*(iy-mask->center.y)));
		break;
	}

	
	return(distance);	
}

/************************************************************************
@Function: shape2d_image_stats
@Description:
	evaluate statistics on image
@Algorithm:
	The background value, average and standard deviation
	are calculates in the region defined by the mask.
	If input images are 3D, the 2d mask is applied to all its depth slices
	Average and standard deviation refer to all the images,
	if a multi-image format provided
@Arguments:
	Bimage*  p				input image - stats modified.
	Shape2d* mask			image mask.
	int 		bground_type		background evaluation method
@Returns:
	int					<0 if not successfull.
**************************************************************************/
int		shape2d_image_stats(Bimage* p, Shape2d* mask, int bground_type)
{
	if ( !p->data ) return(-1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short*	usdata = (unsigned short *) p->data;
    short* 		    	sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    		fdata = (float *) p->data;
    complex_short* 	csdata = (complex_short *) p->data;
    complex_int* 	cidata = (complex_int *) p->data;
    complex_float*	cfdata = (complex_float *) p->data;
    polar*  			pdata = (polar *) p->data;
    
	double		sum, sum2, background;
	float		amp;
    long     	i, number, bgnumber;
	unsigned int x, y, z, isinside;


	sum = sum2 = 0.;
	number = 0;
	
	for ( unsigned int n = 0; n < p->n; n++) {
	
		bgnumber = 0;
		background = 0.;
		for ( z=0; z<p->z; z++ ) {
			for ( y=0; y<p->y; y++ ) {
				for ( x=0; x<p->x; x++ ) {
					i = ((n*p->z+z)*p->y+y)*p->x+x;
					isinside = shape2d_point_inside(mask, x, y);
					// average and standard deviation
					if ( isinside ) {
						switch ( p->datatype ) {
							case UChar:
								sum += udata[i];
								sum2 += udata[i]*1.0*udata[i];
								break;
							case SChar:
								sum += cdata[i];
								sum2 += cdata[i]*1.0*cdata[i];
								break;
							case UShort:
								sum += usdata[i];
								sum2 += usdata[i]*1.0*usdata[i];
								break;
							case Short:
								sum += sdata[i];
								sum2 += sdata[i]*1.0*sdata[i];
								break;
							case Int:
								sum += idata[i];
								sum2 += idata[i]*1.0*idata[i];
								break;
							case Float:
								sum += fdata[i];
								sum2 += fdata[i]*fdata[i];
								break;
							case ComplexShort:
								amp = sqrt(1.0*csdata[i].re*csdata[i].re +
    	    								csdata[i].im*csdata[i].im);
								sum += amp;
								sum2 += amp*amp;
								break;
							case ComplexInt:
								amp = sqrt(1.0*cidata[i].re*cidata[i].re +
    	    								cidata[i].im*cidata[i].im);
								sum += amp;
								sum2 += amp*amp;
								break;
							case ComplexFloat:
								amp = sqrt(cfdata[i].re*cfdata[i].re +
    	    								cfdata[i].im*cfdata[i].im);
								sum += amp;
								sum2 += amp*amp;
								break;
							case Polar:
								sum += pdata[i].amp;
								sum2 += pdata[i].amp*pdata[i].amp;
								break;
							default: break;
						}
						number++;
					}

					// background
					if ( (bground_type==SOUTSIDE && !isinside) || 
					   (bground_type==SINSIDE && isinside) ||
					   (isinside && bground_type==SFRAME &&
					   shape2d_point_inframe(mask, x, y))) {
						switch ( p->datatype ) {
							case UChar:
								background += udata[i];
								break;
							case SChar:
								background += cdata[i];
								break;
							case UShort:
								background += usdata[i];
								break;
							case Short:
								background += sdata[i];
								break;
							case Int:
								background += idata[i];
								break;
							case Float:
								background += fdata[i];
								break;
							case ComplexShort:
								background += sqrt(1.0*csdata[i].re*csdata[i].re +
    	    								csdata[i].im*csdata[i].im);
								break;
							case ComplexInt:
								background += sqrt(1.0*cidata[i].re*cidata[i].re +
    	    								cidata[i].im*cidata[i].im);
								break;
							case ComplexFloat:
								background += sqrt(cfdata[i].re*cfdata[i].re +
    	    								cfdata[i].im*cfdata[i].im);
								break;
							case Polar:
								background += pdata[i].amp;
								break;
							default: break;
						}
						bgnumber++;
					}
				}
			}
		}
		if ( bgnumber ) background = background/bgnumber;
		p->image[n].background = background;
	}

	if (number) {
		p->avg = sum/number;
		p->std = sum2/number - p->avg*p->avg;
		if ( p->std > 0 ) p->std = sqrt(p->std);
		else p->std = 0.;
	}

	return(0);
}

/************************************************************************
@Function: shape2d_mask_image
@Description:
	mask an image according to given shape.
@Algorithm:
	The background value is applied to the area outside the given shape.
	If requested, a tapering is applied.
	If input images are 3D, the 2d mask is applied to all its depth slices
@Arguments:
	Bimage*  p		input image - modified.
	Shape2d* mask	mask to apply (tapering flag included).
	int      taper   taper flag
@Returns:
	int				<0 if not successfull.
**************************************************************************/
int		shape2d_mask_image(Bimage* p, Shape2d* mask, int taper)
{
	if ( !p->data ) return(-1);
	if ( mask->type==SALL && mask->taper==0 ) return(1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    		fdata = (float *) p->data;
    complex_short* 	csdata = (complex_short *) p->data;
    complex_int* 	cidata = (complex_int *) p->data;
    complex_float*  	cfdata = (complex_float *) p->data;
    polar*  			pdata = (polar *) p->data;
    
	float		background, shp_dist;
    unsigned int   n, x, y, z;
    unsigned long  i;
	float   		shp_tapsize = mask->tapsize;
//	int		taper = mask->taper;
	
	for ( n = 0; n < p->n; n++) {
	
		background = p->image[n].background;
		
// mask with background value, and taper if needed
		for ( z=0; z<p->z; z++ ) {
			for ( y=0; y<p->y; y++ ) {
				for ( x=0; x<p->x; x++ ) {
					if ( !shape2d_point_inside(mask, x, y) ) {
						i = ((n*p->z+z)*p->y+y)*p->x+x;
						switch ( p->datatype ) {
							case UChar:
								udata[i] = (unsigned char) background;
								break;
							case SChar:
								cdata[i] = (signed char) background;
								break;
							case UShort:
								usdata[i] = (unsigned short) background;
								break;
							case Short:
								sdata[i] = (short) background;
								break;
							case Int:
								idata[i] = (int) background;
								break;
							case Float:
								fdata[i] = background;
								break;
							case ComplexShort:
								csdata[i].re = (short) (sqrt(0.5)*background);
								csdata[i].im = (short) (sqrt(0.5)*background);
								break;
							case ComplexInt:
								cidata[i].re = (int) (sqrt(0.5)*background);
								cidata[i].im = (int) (sqrt(0.5)*background);
								break;
							case ComplexFloat:
								cfdata[i].re = (sqrt(0.5)*background);
								cfdata[i].im = (sqrt(0.5)*background);
								break;
							case Polar:
								pdata[i].amp = background;
								break;
							default: break;
						}
					} else if (taper == 1) {
						shp_dist = shape2d_point_distance(mask,x,y);
						if (shp_dist < shp_tapsize) {
							i = ((n*p->z+z)*p->y+y)*p->x+x;
							switch ( p->datatype ) {
							case UChar:
								udata[i] = (unsigned char) ((udata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case SChar:
								cdata[i] = (signed char) ((cdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case UShort:
								usdata[i] = (unsigned short) ((usdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case Short:
								sdata[i] = (short) ((sdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case Int:
								idata[i] = (int) ((idata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case Float:
								fdata[i] = ((fdata[i] - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							case ComplexShort:
								csdata[i].re = (short) (sqrt(0.5)*((csdata[i].re - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								csdata[i].im = (short) (sqrt(0.5)*((csdata[i].im - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								break;
							case ComplexInt:
								cidata[i].re = (int) (sqrt(0.5)*((cidata[i].re - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								cidata[i].im = (int) (sqrt(0.5)*((cidata[i].im - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								break;
							case ComplexFloat:
								cfdata[i].re = (sqrt(0.5)*((cfdata[i].re - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								cfdata[i].im = (sqrt(0.5)*((cfdata[i].im - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background));
								break;
							case Polar:
								pdata[i].amp = ((pdata[i].amp - background)*(.5 - .5 * cos(2*M_PI*(shp_dist)/(2*shp_tapsize)))+background);
								break;
							default: break;
							}
						}
					}
				}
			}
		}

		if ( verbose & VERB_FULL )
			printf("\n Image # %d masked with background value: %f \n", n, background);
	}
	
	return(1);
}

/************************************************************************
@Function: shape2d_copy
@Description:
	copy a volume structure to a new variable.
@Algorithm:
@Arguments:
	Shape2d* s		input shape
@Returns:
	Shape2d*         copy of input shape
**************************************************************************/
Shape2d* shape2d_copy( Shape2d* s) {

	Shape2d*		sh = (Shape2d *) balloc(sizeof(Shape2d));	
	
	sh->type = s->type;
	sh->pol = poly_copy(s->pol);
	sh->center.x = s->center.x;
	sh->center.y = s->center.y;
	sh->center.z = s->center.z;
	sh->radius = s->radius;
	sh->taper = s->taper;
	sh->tapsize = s->tapsize;

	return(sh);
}

/************************************************************************
@Function: shape2d_bin
@Description:
	Bins by an integer size.
@Algorithm:
@Arguments:
	Shape2d*   s		input shape
	int		  bin	bin factor
@Returns:
	int 			0.
**************************************************************************/
int shape2d_bin( Shape2d* s, int bin) {

	if (s==NULL) return(-1);
	
	Shape2d* sh = s;
	

	sh->center.x = (sh->center.x + bin - 1) / bin;
	sh->center.y = (sh->center.y + bin - 1) / bin;
	sh->radius = (sh->radius + bin - 1) / bin;

	shape2d_taper_init(sh, sh->taper);

	return(0);
}

/************************************************************************
@Function: shape2d_fill_ratio
@Description:
	return the fill ratio of the shape with respect to the given projecion.
@Algorithm:
@Arguments:
	Bimage*  p		input image - modified.
	Shape2d* mask	reference mask.
@Returns:
	float			fill ratio.
**************************************************************************/
float		shape2d_fill_ratio(Bimage* p, Shape2d* mask)
{
	
	if (p==NULL || mask == NULL) return(0.);
	
	float parea = 1.;
	if (p->x>1) parea*= p->x;
	if (p->y>1) parea*= p->y;
	
	float sarea = 1.;
	switch (mask->type) {
		case SALL:
			sarea = parea;
		break;	
		case SPOLY:
			sarea = poly_area(mask->pol);
		break;
		case SCIRCLE:
			// area of square containing the circle
			sarea = 4.0*mask->radius*mask->radius;
		break;
	}
	
	return(sarea/parea);
	
}

