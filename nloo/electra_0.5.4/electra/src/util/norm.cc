/*
	norm.cc
	Functions for normalizing images
	Author: Giovanni Cardone
	Created: 20050113 	Modified:
*/

#include "norm.h"


/************************************************************************
@Function: norm3d_rescale_to_avg_std
@Description:
	Rescales the image data to a given average and standard deviation.
@Algorithm:
	The new data is calculated as:
		                            new_std_dev
		new_datum = (datum - avg) * ----------- + new_avg
		                              std_dev
	The new data replaces the old data.
	The rescaling is applied only inside the region of intereste,
	while outside of it the given background value is applied
	Image statistics are not recalculated.
@Arguments:
	Bimage* p		image.
	float avg		new average.
	float std		new standard deviation.
	float* background	background value (optional).
	Shape3d* mask	reference mask
@Returns:
	int 			0.
**************************************************************************/
int norm3d_rescale_to_avg_std(Bimage* p, float avg, float std, float* background, Shape3d* mask)
{
	if ( !p->dataflag ) return(-1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    		fdata = (float *) p->data;


	if ( std < 0 ) {
		fprintf(stderr, "Warning: Cannot use a negative standard deviation to scale to! (%g)\n", std);
		return(-1);
	}
	
	if ( p->std == 0 ) {
		fprintf(stderr, "Error in previous rescaling: Standard deviation is zero!\n");
		return(-1);
	}
	
    float			scale = std/p->std;
	float			shift = avg - p->avg*scale;
    
	if ( verbose & VERB_LABEL )
	    printf("Rescaling inside mask to average and stdev: %g %g\n", avg, std);
	
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	if ( p->datatype >= ComplexShort ) datasize *= 2;
	
    double			tdata, min, max;
	p->scale = scale;
    
	if ( verbose & VERB_PROCESS )
    	printf("Scale and shift:    %g %g\n\n", scale, shift);
	
	// Adjust the background values
	if ( background != NULL) {
		for ( i=0; i<p->n; i++ )
			p->image[i].background = *background;
	} else {
		for ( i=0; i<p->n; i++ )
			p->image[i].background = p->image[i].background*scale+shift;
	}
	
	// GC: implementation not efficient: rescaling applied to all points,
	// then a mask is applied
    switch ( p->datatype ) {
    	case UChar:
    	    for ( i=0; i<datasize; i++ ) {
    	    		tdata = udata[i]*scale + shift;
			if ( tdata < 0 ) udata[i] = 0;
			else if ( tdata > 255 ) udata[i] = 255;
    	    		else udata[i] = (unsigned char)tdata;
		}
    	    break;
    	case SChar:
    	    for ( i=0; i<datasize; i++ ) {
    	    		tdata = cdata[i]*scale + shift;
			if ( tdata < -128 ) cdata[i] = -128;
			else if ( tdata > 127 ) cdata[i] = 127;
    	    		else cdata[i] = (signed char)tdata;
		}
    	    break;
    	case UShort:
		max = USHRT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
    	    		tdata = usdata[i]*scale + shift;
			if ( tdata < 0 ) usdata[i] = 0;
			else if ( tdata > max ) usdata[i] = (unsigned short)max;
    	    		else usdata[i] = (unsigned short)tdata;
		}
    	    break;
    	case ComplexShort:
    	case Short:
		min = SHRT_MIN;
		max = SHRT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
    	    		tdata = sdata[i]*scale + shift;
			if ( tdata < min ) sdata[i] = (short)min;
			else if ( tdata > max ) sdata[i] = (short)max;
    	    		else sdata[i] = (short)tdata;
		}
    	    break;
    	case ComplexInt:
    	case Int:
		min = INT_MIN;
		max = INT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
    	    		tdata = idata[i]*scale + shift;
			if ( tdata < min ) idata[i] = (int)min;
			else if ( tdata > max ) idata[i] = (int)max;
    	    		else idata[i] = (int)tdata;
		}
    	    break;
    	case ComplexFloat:
    	case Float:
    	    for ( i=0; i<datasize; i++ )
    	    		fdata[i] = fdata[i]*scale + shift;
    	    break;
    	case Polar:
    	    for ( i=0; i<datasize; i+=2 ) {
    	    		fdata[i] = fdata[i]*scale + shift;
			if ( fdata[i] < 0 ) fdata[i] = 0;
		}
    	    break;
   		default: break;
    }
    
	// apply background outer of the mask
    if (mask != NULL )
		shape3d_mask_image(p, mask, 0);
    
	return(0);
}

/************************************************************************
@Function: norm3d_rescale_linregression
@Description:
	normalize an image by least squares apporach (linear regression).
@Algorithm:
	minimize the least squares error inside the mask domain between
	the given image and a reference.
	Find the two values a and b in order to apply the transformation
	xnew = a*x+b
	Out of the mask te background value from the reference is assigned.
@Arguments:
	Bimage* p		input image - modified.
	Bimage* pref		reference image.
	Shape3d* mask	reference mask
@Returns:
	int				<0 if not successfull.
**************************************************************************/
int		norm3d_rescale_linregression(Bimage* p, Bimage* pref, Shape3d* mask)
{
	if ( !p->data ) return(-1);
	if ( !pref->data ) return(-1);

	if (p->n!=1 || pref->n!=1) {
		fprintf(stderr,"Warning(norm3d_rescale_linregression):");
		fprintf(stderr," multiple-images format not allowed!\n");
		fprintf(stderr,"  normalization not performed\n");
		return(-1);
	}
    
	Bimage*		bp = NULL;
	Bimage*		bpr = NULL;
	
	VectorInt3	bin = {4,4,4};
	float		a=0., b=0.;
	double		sx, sy, sx2, sxy;
    unsigned int  n, x, y, z, ns;
    unsigned long  i;

	float*  		bfdata = NULL;
	float*  		bfrefdata = NULL;
/*
	double fr = (double) shape3d_fill_ratio(p,mask);
	fr = pow(fr, 1.0/3.0);
	bin.x = (int) ( bin.x/fr + 0.5 );
	bin.y = (int) ( bin.y/fr + 0.5 );
	bin.z = (int) ( bin.z/fr + 0.5 );
*/	
	Shape3d* mask_bin = shape3d_copy(mask);
	shape3d_bin(mask_bin, bin.x);

	bp = copy_img(p);
	img_bin(bp,bin);
	img_to_float(bp);

	bpr = copy_img(pref);
	img_bin(bpr,bin);
	img_to_float(bpr);

	bfdata = (float *) bp->data;
	bfrefdata = (float *) bpr->data;

	for ( n = 0; n < p->n; n++) {
	
// evaluate terms to calculate a and b
		sx = sx2 = sy = sxy = 0.;
		ns = 0;

		for ( z=0; z<bp->z; z++ ) {
			for ( y=0; y<bp->y; y++ ) {
				for ( x=0; x<bp->x; x++ ) {
					if ( shape3d_point_inside(mask_bin, x, y, z) ) {
						i = ((n*p->z+z)*bp->y+y)*bp->x+x;
						sx += bfdata[i];
						sx2 += (bfdata[i]*bfdata[i]);
						sy += bfrefdata[i];
						sxy += (bfdata[i]*bfrefdata[i]);
						ns++;
					}
				}
			}
		}

		a = (ns*sxy-sx*sy)/(ns*sx2-sx*sx);
		b = (sy-a*sx)/ns;

		if ( verbose & VERB_FULL )
			printf("\n # %d Least squares normalization - a & b values: %f , %f\n", n, a, b);

	}

	float avg = p->avg * a + b;
	float std = p->std * fabs(a);

	norm3d_rescale_to_avg_std(p, avg, std, (float*) NULL, mask);	

	shape3d_kill(mask_bin);
	kill_img(bp);
	kill_img(bpr);
	
	return(1);
}

/************************************************************************
@Function: norm2d_rescale_to_avg_std
@Description:
	Rescales the image data to a given average and standard deviation.
@Algorithm:
	The new data is calculated as:
		                            new_std_dev
		new_datum = (datum - avg) * ----------- + new_avg
		                              std_dev
	The new data replaces the old data.
	The rescaling is applied only inside the region of intereste,
	while outside of it the given background value is applied
	Image statistics are not recalculated.
@Arguments:
	Bimage* p		image.
	float avg		new average.
	float std		new standard deviation.
	float* background	background value (optional).
	Shape2d* mask	reference mask
@Returns:
	int 			0.
**************************************************************************/
int norm2d_rescale_to_avg_std(Bimage* p, float avg, float std, float* background, Shape2d* mask)
{
	if ( !p->dataflag ) return(-1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    		fdata = (float *) p->data;


	if ( std < 0 ) {
		fprintf(stderr, "Warning: Cannot use a negative standard deviation to scale to! (%g)\n", std);
		return(-1);
	}
	
	if ( p->std == 0 ) {
		fprintf(stderr, "Error in previous rescaling: Standard deviation is zero!\n");
		return(-1);
	}
	
    float			scale = std/p->std;
	float			shift = avg - p->avg*scale;
    
	if ( verbose & VERB_LABEL )
	    printf("Rescaling inside mask to average and stdev: %g %g\n", avg, std);
	
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	if ( p->datatype >= ComplexShort ) datasize *= 2;
	
    double			tdata, min, max;
	p->scale = scale;
    
	if ( verbose & VERB_PROCESS )
    	printf("Scale and shift:    %g %g\n\n", scale, shift);
	
	// Adjust the background values
	if ( background != NULL) {
		for ( i=0; i<p->n; i++ )
			p->image[i].background = *background;
	} else {
		for ( i=0; i<p->n; i++ )
			p->image[i].background = p->image[i].background*scale+shift;
	}
	
	// GC: implementation not efficient: rescaling applied to all points,
	// then a mask is applied
    switch ( p->datatype ) {
    	case UChar:
    	    for ( i=0; i<datasize; i++ ) {
    	    		tdata = udata[i]*scale + shift;
			if ( tdata < 0 ) udata[i] = 0;
			else if ( tdata > 255 ) udata[i] = 255;
    	    		else udata[i] = (unsigned char)tdata;
		}
    	    break;
    	case SChar:
    	    for ( i=0; i<datasize; i++ ) {
    	    		tdata = cdata[i]*scale + shift;
			if ( tdata < -128 ) cdata[i] = -128;
			else if ( tdata > 127 ) cdata[i] = 127;
    	    		else cdata[i] = (signed char)tdata;
		}
    	    break;
    	case UShort:
		max = USHRT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
    	    		tdata = usdata[i]*scale + shift;
			if ( tdata < 0 ) usdata[i] = 0;
			else if ( tdata > max ) usdata[i] = (unsigned short)max;
    	    		else usdata[i] = (unsigned short)tdata;
		}
    	    break;
    	case ComplexShort:
    	case Short:
		min = SHRT_MIN;
		max = SHRT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
    	    		tdata = sdata[i]*scale + shift;
			if ( tdata < min ) sdata[i] = (short)min;
			else if ( tdata > max ) sdata[i] = (short)max;
    	    		else sdata[i] = (short)tdata;
		}
    	    break;
    	case ComplexInt:
    	case Int:
		min = INT_MIN;
		max = INT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
    	    		tdata = idata[i]*scale + shift;
			if ( tdata < min ) idata[i] = (int)min;
			else if ( tdata > max ) idata[i] = (int)max;
    	    		else idata[i] = (int)tdata;
		}
    	    break;
    	case ComplexFloat:
    	case Float:
    	    for ( i=0; i<datasize; i++ )
    	    		fdata[i] = fdata[i]*scale + shift;
    	    break;
    	case Polar:
    	    for ( i=0; i<datasize; i+=2 ) {
    	    		fdata[i] = fdata[i]*scale + shift;
			if ( fdata[i] < 0 ) fdata[i] = 0;
		}
    	    break;
   		default: break;
    }
    
	// apply background outer of the mask
    if (mask != NULL )
		shape2d_mask_image(p, mask, 0);
    
	return(0);
}

/************************************************************************
@Function: norm2d_rescale_linregression
@Description:
	normalize an image by least squares apporach (linear regression).
@Algorithm:
	minimize the least squares error inside the mask domain between
	the given image and a reference.
	Find the two values a and b in order to apply the transformation
	xnew = a*x+b
	Out of the mask the background value from the reference is assigned.
@Arguments:
	Bimage* p		input image - modified.
	Bimage* pref		reference image.
	Shape2d* mask	reference mask
@Returns:
	int				<0 if not successfull.
**************************************************************************/
int		norm2d_rescale_linregression(Bimage* p, Bimage* pref, Shape2d* mask)
{
	if ( !p->data ) return(-1);
	if ( !pref->data ) return(-1);

	if (p->n!=1 || pref->n!=1) {
		fprintf(stderr,"Warning(norm2d_rescale_linregression):");
		fprintf(stderr," multiple-images format not allowed!\n");
		fprintf(stderr,"  normalization not preformed\n");
		return(-1);
	}
    
	Bimage*		bp = NULL;
	Bimage*		bpr = NULL;
	
	VectorInt3	bin = {4,4,4};
	float		a=0., b=0.;
	double		sx, sy, sx2, sxy;
    unsigned int n, x, y, z, ns;
    unsigned long i;

	float*  		bfdata = NULL;
	float*  		bfrefdata = NULL;
/*
	double fr = (double) shape2d_fill_ratio(p,mask);
	fr = sqrt(fr);
	bin.x = (int) ( bin.x/fr + 0.5 );
	bin.y = (int) ( bin.y/fr + 0.5 );
	bin.z = (int) ( bin.z/fr + 0.5 );
*/	
	Shape2d* mask_bin = shape2d_copy(mask);
	shape2d_bin(mask_bin, bin.x);

	bp = copy_img(p);
	img_bin(bp,bin);
	img_to_float(bp);

	bpr = copy_img(pref);
	img_bin(bpr,bin);
	img_to_float(bpr);

	bfdata = (float *) bp->data;
	bfrefdata = (float *) bpr->data;

	for ( n = 0; n < p->n; n++) {
	
// evaluate terms to calculate a and b
		sx = sx2 = sy = sxy = 0.;
		ns = 0;

		for ( z=0; z<bp->z; z++ ) {
			for ( y=0; y<bp->y; y++ ) {
				for ( x=0; x<bp->x; x++ ) {
					if ( shape2d_point_inside(mask_bin, x, y) ) {
						i = ((n*p->z+z)*bp->y+y)*bp->x+x;
						sx += bfdata[i];
						sx2 += (bfdata[i]*bfdata[i]);
						sy += bfrefdata[i];
						sxy += (bfdata[i]*bfrefdata[i]);
						ns++;
					}
				}
			}
		}

		a = (ns*sxy-sx*sy)/(ns*sx2-sx*sx);
		b = (sy-a*sx)/ns;

		if ( verbose & VERB_FULL )
			printf("\n # %d Least squares normalization - a & b values: %f , %f\n", n, a, b);

	}

	float avg = p->avg * a + b;
	float std = p->std * fabs(a);

	norm2d_rescale_to_avg_std(p, avg, std, (float*) NULL, mask);	

	shape2d_kill(mask_bin);
	kill_img(bp);
	kill_img(bpr);
	
//	img_stats(p);
	
	return(1);
}

