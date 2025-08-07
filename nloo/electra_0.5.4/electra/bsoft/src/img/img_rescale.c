/*
	img_rescale.c
	Library routines used for images
	Author: Bernard Heymann
	Created: 19990321 	Modified: 20040415
*/

#include "rwimg.h"
#include "img_rescale.h"
#include "img_datatypes.h"
#include "img_complex.h"
#include "matrix.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated

/************************************************************************
@Function: img_invert_data
@Description:
	Inverts the image density data.
@Algorithm:
	The voxel values are negated and re-adjusted to fit into the data type
	range for the UChar and UShort types.
	The new data replaces the old data.
	Image statistics are recalculated.
@Arguments:
	Bimage* p		image.
@Returns:
	int 			0.
**************************************************************************/
int 		img_invert_data(Bimage* p)
{
	if ( !p->dataflag ) return(-1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* usdata = (unsigned short *) p->data;
    short* 	    	sdata = (short *) p->data;
    int* 	    	idata = (int *) p->data;
    float*  	    fdata = (float *) p->data;
    
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	int 			short_max = 256*256;
	
	// Adjust the background values
	for ( i=0; i<p->n; i++ )
		p->image[i].background = -p->image[i].background;
	
 	if ( verbose & VERB_PROCESS ) printf("Inverting:\n");
 	else if ( verbose & VERB_LABEL ) printf("Inverting\n\n");
	
    switch ( p->datatype ) {
    	case UChar:
			if ( p->colormodel == Bit ) {
				for ( i=0; i<datasize/8; i++ ) udata[i] = ~udata[i];
    	    	for ( i=0; i<p->n; i++ ) p->image[i].background = 1 + p->image[i].background;
			} else {
				for ( i=0; i<datasize; i++ ) udata[i] = 255 - udata[i];
    	    	for ( i=0; i<p->n; i++ ) p->image[i].background = 255 + p->image[i].background;
			}
    	    break;
    	case SChar:
    	    for ( i=0; i<datasize; i++ ) cdata[i] = -cdata[i];
    	    break;
    	case UShort:
    	    for ( i=0; i<datasize; i++ ) usdata[i] = short_max - usdata[i];
    	    for ( i=0; i<p->n; i++ ) p->image[i].background = short_max + p->image[i].background;
    	    break;
    	case ComplexShort:
			datasize *= 2;
    	case Short:
    	    for ( i=0; i<datasize; i++ ) sdata[i] = -sdata[i];
    	    break;
    	case ComplexInt:
			datasize *= 2;
    	case Int:
    	    for ( i=0; i<datasize; i++ ) idata[i] = -idata[i];
    	    break;
    	case ComplexFloat:
    	case Polar:
			datasize *= 2;
    	case Float:
    	    for ( i=0; i<datasize; i++ ) fdata[i] = -fdata[i];
    	    break;
    	default: break;
    }
    
	img_stats(p);

 	if ( verbose & VERB_PROCESS )
		printf("New min, max, avg and std:      %g %g %g %g\n\n",
			p->min, p->max, p->avg, p->std);

	return(0);
}

/************************************************************************
@Function: img_rescale
@Description:
	Rescales the image data to a given minimum and maximum.
@Algorithm:
	The new data is calculated as:
		new_datum = datum*scale + shift
	The new data replaces the old data.
	Image statistics are recalculated.
@Arguments:
	Bimage* p		image.
	float scale		multiplier.
	float shift		addition.
@Returns:
	int 			0.
**************************************************************************/
int 		img_rescale(Bimage* p, float scale, float shift)
{
	if ( !p->dataflag ) return(-1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* usdata = (unsigned short *) p->data;
    short* 	    	sdata = (short *) p->data;
    int* 	    	idata = (int *) p->data;
    float*  	    fdata = (float *) p->data;
    
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	if ( p->datatype >= ComplexShort ) datasize *= 2;
	
    double			tdata, min, max;
	p->scale = scale;
    
	if ( verbose & VERB_PROCESS )
    	printf("Scale and shift:                %g %g\n\n", scale, shift);
	
	// Adjust the background values
	for ( i=0; i<p->n; i++ )
		p->image[i].background = p->image[i].background*scale + shift;
	
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
    
	img_stats(p);
    
	return(0);
}

/************************************************************************
@Function: img_rescale_to_min_max
@Description:
	Rescales the image data to a given minimum and maximum.
@Algorithm:
	The new data is calculated as:
		                            new_max - new_min
		new_datum = (datum - min) * ----------------- + new_min
		                                max - min
	The new data replaces the old data.
	Image statistics are recalculated.
@Arguments:
	Bimage* p		image.
	float min		new minimum.
	float max		new maximum.
@Returns:
	int 			0.
**************************************************************************/
int 		img_rescale_to_min_max(Bimage* p, float min, float max)
{
	if ( !p->dataflag ) return(-1);
	
    float			scale = (max - min)/(p->max - p->min);
	float			shift = min - p->min*scale;
    
	if ( verbose & VERB_LABEL )
	    printf("Rescaling to:                   %g %g\n", min, max);
	
	return(img_rescale(p, scale, shift));
}

/************************************************************************
@Function: img_rescale_to_avg_std
@Description:
	Rescales the image data to a given average and standard deviation.
@Algorithm:
	The new data is calculated as:
		                            new_std_dev
		new_datum = (datum - avg) * ----------- + new_avg
		                              std_dev
	The new data replaces the old data.
	Image statistics are recalculated.
@Arguments:
	Bimage* p		image.
	float avg		new average.
	float std		new standard deviation.
@Returns:
	int 			0.
**************************************************************************/
int 		img_rescale_to_avg_std(Bimage* p, float avg, float std)
{
	if ( !p->dataflag ) return(-1);
	
	if ( std < 0 ) {
		fprintf(stderr, "Warning: Cannot use a negative standard deviation to scale to! (%g)\n", std);
		return(-1);
	}
	
	if ( p->std < 1e-30 ) if ( img_stats(p) < 0 ) {
		fprintf(stderr, "Error in img_stats: Statistics not recalculated!\n");
		return(-1);
	}
	
	if ( p->std == 0 ) {
		fprintf(stderr, "Error in rescaling: Standard deviation is zero!\n");
		return(-1);
	}
	
    float			scale = std/p->std;
	float			shift = avg - p->avg*scale;
    
	if ( verbose & VERB_LABEL )
	    printf("Rescaling to average and stdev: %g %g\n", avg, std);
	
	return(img_rescale(p, scale, shift));
}

/************************************************************************
@Function: img_square
@Description:
	Calculates the square of the image.
@Algorithm:
	The new data is calculated as:
		new_datum = datum*datum
	A new image is generated.
	Image statistics are recalculated.
@Arguments:
	Bimage* p		image.
@Returns:
	Bimage* 		squared image.
**************************************************************************/
int 		img_square(Bimage* p)
{
	if ( !p->dataflag ) return(-1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* usdata = (unsigned short *) p->data;
    short* 	    	sdata = (short *) p->data;
    int* 	    	idata = (int *) p->data;
    float*  	    fdata = (float *) p->data;
    
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	if ( p->datatype >= ComplexShort ) datasize *= 2;
	
    float  	    	tdata, min, max;
    
	if ( verbose & VERB_LABEL ) {
	    printf("Squaring the image\n\n");
	}
	
	// Adjust the background values
	for ( i=0; i<p->n; i++ )
		p->image[i].background *= p->image[i].background;
	
    switch ( p->datatype ) {
    	case UChar:
    	    for ( i=0; i<datasize; i++ ) {
    	    	tdata = udata[i]*udata[i];
    	    	udata[i] = (unsigned char)tdata;
				if ( tdata < 0 ) udata[i] = 0;
				if ( tdata > 255 ) udata[i] = 255;
			}
    	    break;
    	case SChar:
    	    for ( i=0; i<datasize; i++ ) {
    	    	tdata = cdata[i]*cdata[i];
    	    	cdata[i] = (signed char)tdata;
				if ( tdata < -128 ) cdata[i] = -128;
				if ( tdata > 127 ) cdata[i] = 127;
			}
    	    break;
    	case UShort:
			max = USHRT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
    	    	tdata = usdata[i]*usdata[i];
    	    	usdata[i] = (unsigned short)tdata;
				if ( tdata < 0 ) usdata[i] = 0;
				if ( tdata > max ) usdata[i] = (unsigned short)max;
			}
    	    break;
    	case ComplexShort:
    	case Short:
			min = SHRT_MIN;
			max = SHRT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
    	    	tdata = sdata[i]*sdata[i];
    	    	sdata[i] = (short)tdata;
				if ( tdata < min ) sdata[i] = (short)min;
				if ( tdata > max ) sdata[i] = (short)max;
			}
    	    break;
    	case ComplexInt:
    	case Int:
			min = INT_MIN;
			max = INT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
    	    	tdata = idata[i]*idata[i];
    	    	idata[i] = (int)tdata;
				if ( tdata < min ) idata[i] = (int)min;
				if ( tdata > max ) idata[i] = (int)max;
			}
    	    break;
    	case ComplexFloat:
    	case Float:
    	    for ( i=0; i<datasize; i++ )
    	    	fdata[i] = fdata[i]*fdata[i];
    	    break;
     	case Polar:
    	    for ( i=0; i<datasize; i+=2 )
    	    	fdata[i] = fdata[i]*fdata[i];
    	    break;
   		default: break;
    }
    
	img_stats(p);
    
	return(0);
}

/************************************************************************
@Function: img_logarithm
@Description:
	Calculates the logarithm of the image data.
@Algorithm:
	The image is first converted to floating point, complex floating point
	for complex data, or left polar for polar data types. The logarithm is 
	calculated to ensure that the whole range of values are included and
	that the average is at zero:
		               data - min + fudge
		new_data = log ------------------
		               avg - min + fudge
	where:
		fudge = (avg - min) / 1000
	The new data replaces the old data and the image is converted back to
	the input data type.
	Image statistics are recalculated.
@Arguments:
	Bimage* p		image.
@Returns:
	int 			0.
**************************************************************************/
int 		img_logarithm(Bimage* p)
{
	if ( !p->dataflag ) return(-1);
	
	DataType		olddatatype = p->datatype;
	
	if ( p->datatype < Float ) img_to_float(p);
	else if ( p->datatype >= ComplexShort && p->datatype <= ComplexInt )
		img_to_complex_float(p);
	
	unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	float			fudge = (p->avg - p->min)/1000;
    float			norm = 1.0/(p->avg - p->min + fudge);
    float*  	    fdata = (float *) p->data;
	
	if ( verbose & VERB_LABEL )
	    printf("Calculating the logarithm of the image\n\n");
		
    switch ( p->datatype ) {
    	case ComplexFloat:
			datasize *= 2;
    	case Float:
    	    for ( i=0; i<datasize; i++ )
				fdata[i] = log((fdata[i] - p->min + fudge)*norm);
    	    break;
    	case Polar:
			datasize *= 2;
    	    for ( i=0; i<datasize; i+=2 )
				fdata[i] = log((fdata[i] - p->min + fudge)*norm);
    	    break;
    	default: break;
    }
	
	img_stats(p);

	img_to_datatype(p, olddatatype);
	
	return(0);
}

/************************************************************************
@Function: img_truncate
@Description:
	Truncates image data to a given minimum and maximum.
@Algorithm:
	All values smaller than the new minimum are set to the new minimum
	and all values larger than the new maximum are set to the new maximum.
	In cases where the given minimum and maximum are outside the data
	type ranges, they are reset to the data type range limits.
	The new data replaces the old data.
	Image statistics are recalculated.
@Arguments:
	Bimage* p		image.
	float min		minimum.
	float max		maximum.
	float setmin	value to set voxels smaller than minimum.
	float setmax	value to set voxels larger than maximum.
@Returns:
	int 			0.
**************************************************************************/
int 		img_truncate(Bimage* p, float min, float max, float setmin, float setmax)
{
	if ( !p->dataflag ) return(-1);
	
    unsigned char 	usetmin = (unsigned char) (setmin + 0.5);
    signed char 	csetmin = (signed char) (setmin + 0.5);
    unsigned short  ussetmin = (unsigned short) (setmin + 0.5);
    short 	    	ssetmin = (short) (setmin + 0.5);
    int 	    	isetmin = (int) (setmin + 0.5);
	
    unsigned char 	usetmax = (unsigned char) (setmax + 0.5);
    signed char 	csetmax = (signed char) (setmax + 0.5);
    unsigned short  ussetmax = (unsigned short) (setmax + 0.5);
    short 	    	ssetmax = (short) (setmax + 0.5);
    int 	    	isetmax = (int) (setmax + 0.5);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* usdata = (unsigned short *) p->data;
    short* 	    	sdata = (short *) p->data;
    int* 	    	idata = (int *) p->data;
    float*  	    fdata = (float *) p->data;
    polar*  	    pdata = (polar *) p->data;
    
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
    
	// Adjust the background values
	for ( i=0; i<p->n; i++ ) {
		if ( p->image[i].background < min ) p->image[i].background = setmin;
		if ( p->image[i].background > max ) p->image[i].background = setmax;
	}
	
    switch ( p->datatype ) {
    	case UChar:
			if ( min < 0 ) min = 0;
			if ( max > 255 ) max = 255;
			if ( setmin < 0 ) usetmin = 0;
			if ( setmax > 255 ) usetmax = 255;
    	    for ( i=0; i<datasize; i++ ) {
				if ( udata[i] < min ) udata[i] = usetmin;
				if ( udata[i] > max ) udata[i] = usetmax;
			}
    	    break;
    	case SChar:
			if ( min < -128 ) min = -128;
			if ( max > 127 ) max = 127;
			if ( setmin < -128 ) csetmin = -128;
			if ( setmax > 127 ) csetmax = 127;
    	    for ( i=0; i<datasize; i++ ) {
				if ( cdata[i] < min ) cdata[i] = csetmin;
				if ( cdata[i] > max ) cdata[i] = csetmax;
			}
    	    break;
    	case UShort:
			if ( min < 0 ) min = 0;
			if ( max > USHRT_MAX ) max = USHRT_MAX;
			if ( setmin < 0 ) ussetmin = 0;
			if ( setmax > USHRT_MAX ) ussetmax = USHRT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
				if ( usdata[i] < min ) usdata[i] = ussetmin;
				if ( usdata[i] > max ) usdata[i] = ussetmax;
			}
    	    break;
    	case ComplexShort:
			datasize *= 2;
    	case Short:
			if ( min < SHRT_MIN ) min = SHRT_MIN;
			if ( max > SHRT_MAX ) max = SHRT_MAX;
			if ( setmin < SHRT_MIN ) ssetmin = SHRT_MIN;
			if ( setmax > SHRT_MAX ) ssetmax = SHRT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
				if ( sdata[i] < min ) sdata[i] = ssetmin;
				if ( sdata[i] > max ) sdata[i] = ssetmax;
			}
    	    break;
    	case Int:
			if ( min < INT_MIN ) min = INT_MIN;
			if ( max > INT_MAX ) max = INT_MAX;
			if ( setmin < INT_MIN ) isetmin = INT_MIN;
			if ( setmax > INT_MAX ) isetmax = INT_MAX;
    	    for ( i=0; i<datasize; i++ ) {
				if ( idata[i] < min ) idata[i] = isetmin;
				if ( idata[i] > max ) idata[i] = isetmax;
			}
    	    break;
    	case ComplexFloat:
			datasize *= 2;
    	case Float:
    	    for ( i=0; i<datasize; i++ ) {
				if ( fdata[i] < min ) fdata[i] = setmin;
				if ( fdata[i] > max ) fdata[i] = setmax;
			}
    	    break;
		case Polar:
    	    for ( i=0; i<datasize; i++ ) {
				if ( pdata[i].amp < min ) pdata[i].amp = setmin;
				if ( pdata[i].amp > max ) pdata[i].amp = setmax;
			}
    	    break;
    	default: break;
    }
	
	if ( verbose & VERB_PROCESS ) {
	    printf("Truncating to:                  %g %g\n", min, max);
	    printf("Min and max replacement values: %g %g\n\n", setmin, setmax);
	}
	
	img_stats(p);

	return(0);
}

/************************************************************************
@Function: img_truncate_to_min_max
@Description:
	Truncates image data to a given minimum and maximum.
@Algorithm:
	All values smaller than the new minimum are set to the new minimum
	and all values larger than the new maximum are set to the new maximum.
	In cases where the given minimum and maximum are outside the data
	type ranges, they are reset to the data type range limits.
	The new data replaces the old data.
	Image statistics are recalculated.
@Arguments:
	Bimage* p		image.
	float min		minimum.
	float max		maximum.
@Returns:
	int 			0.
**************************************************************************/
int 		img_truncate_to_min_max(Bimage* p, float min, float max)
{
	return(img_truncate(p, min, max, min, max));
}

/************************************************************************
@Function: img_truncate_to_avg
@Description:
	Sets voxels in image data exceeding a given minimum and maximum to
	the average.
@Algorithm:
	All values smaller than the new minimum or larger than the new 
	maximum are set to the image average.
	The new data replaces the old data.
	Image statistics are recalculated.
@Arguments:
	Bimage* p		image.
	float min		minimum.
	float max		maximum.
@Returns:
	int 			0.
**************************************************************************/
int 		img_truncate_to_avg(Bimage* p, float min, float max)
{
	return(img_truncate(p, min, max, p->avg, p->avg));
}

/************************************************************************
@Function: img_limit_levels
@Description:
	Converts a full gray scale image to a limited level image.
@Algorithm:
	The dynamic range of the image is decreased to a given number of gray
	scale levels.
	The new data replaces the old data.
	Image statistics are recalculated.
@Arguments:
	Bimage* p		image.
	int nlevels		number of levels.
@Returns:
	int 			0.
**************************************************************************/
int 		img_limit_levels(Bimage* p, int nlevels)
{
	if ( !p->dataflag ) return(-1);
	
	if ( nlevels > 255 ) {
		fprintf(stderr, "Error: There should be less than 255 gray levels (%d)\n", nlevels);
		return(-1);
	}
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* usdata = (unsigned short *) p->data;
    short* 	    	sdata = (short *) p->data;
    int* 	    	idata = (int *) p->data;
    float*  	    fdata = (float *) p->data;
    polar*  	    pdata = (polar *) p->data;
    
    unsigned long   i, datasize = (unsigned long) p->c*p->x*p->y*p->z*p->n;
    float			unit = 255.0/(nlevels - 1);
    float			scale = 0.999*nlevels/(p->max - p->min);
	
    char* 			newdata = (char *) balloc(datasize*sizeof(char));
    
	if ( verbose & VERB_LABEL )
	    printf("Restricting to %d levels.\n\n",nlevels);
	
    switch ( p->datatype ) {
    	case UChar:
    	    for ( i=0; i<datasize; i++ )
    	    	newdata[i] = (char) (unit*floor((udata[i] - p->min)*scale));
    	    break;
    	case SChar:
    	    for ( i=0; i<datasize; i++ )
    	    	newdata[i] = (char) (unit*floor((cdata[i] - p->min)*scale));
    	    break;
    	case UShort:
    	    for ( i=0; i<datasize; i++ )
    	    	newdata[i] = (char) (unit*floor((usdata[i] - p->min)*scale));
    	    break;
    	case ComplexShort:
			datasize *= 2;
    	case Short:
    	    for ( i=0; i<datasize; i++ )
    	    	newdata[i] = (char) (unit*floor((sdata[i] - p->min)*scale));
    	    break;
    	case ComplexInt:
			datasize *= 2;
    	case Int:
    	    for ( i=0; i<datasize; i++ )
    	    	newdata[i] = (char) (unit*floor((idata[i] - p->min)*scale));
    	    break;
    	case ComplexFloat:
			datasize *= 2;
    	case Float:
    	    for ( i=0; i<datasize; i++ )
    	    	newdata[i] = (char) (unit*floor((fdata[i] - p->min)*scale));
    	    break;
    	case Polar:
    	    for ( i=0; i<datasize; i++ )
    	    	newdata[i] = (char) (unit*floor((pdata[i].amp - p->min)*scale));
    	    break;
    	default: break;
    }
	
	p->data = (char *) newdata;
	
	bfree(udata, datasize*gettypesize(p->datatype));
	
	p->datatype = UChar;
	img_stats(p);
	
	return(0);
}
