/*
	img_complex.c
	Routines to convert complex data sets
	Author: Bernard Heymann
	Created: 19990424  	    Modified: 20040714
*/
	
#include "rwimg.h"
#include "img_complex.h"
#include "img_combine.h"
#include "img_datatypes.h"
#include "img_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

/************************************************************************
@Function: img_to_complex_short
@Description:
	Converts any data type to a complex short.
@Algorithm:
	No statistics are calculated.
@Arguments:
	Bimage* p			image (replaced by new image).
@Returns:
	int					error code.
**************************************************************************/
int			img_to_complex_short(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
    if ( p->datatype == ComplexShort ) return(0);
	
	if ( p->colormodel == RGB ) img_RGB2gray(p);
    
    unsigned long   	i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	unsigned long		elementsize = gettypesize(p->datatype);
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    complex_int*  		cidata = (complex_int *) p->data;
    complex_float*  	cfdata = (complex_float *) p->data;
    polar*		  		pdata = (polar *) p->data;
    complex_short*  	newdata = (complex_short *) balloc(datasize*sizeof(complex_short));
	if ( !newdata )
		return(error_show("img_to_complex_short", __FILE__, __LINE__));
    
    switch ( p->datatype ) {
    	case UChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned char -> complex_short\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = udata[i];
    	    break;
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> complex_short\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = cdata[i];
    	    break;
    	case UShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned short -> complex_short\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = usdata[i];
    	    break;
    	case Short:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: short -> complex_short\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = sdata[i];
    	    break;
    	case Int:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: int -> complex_short\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = idata[i];
    	    break;
    	case ComplexInt:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_int -> complex_short\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				newdata[i].re = cidata[i].re;
				newdata[i].im = cidata[i].im;
			}
    	    break;
    	case Float:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: float -> complex_short\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = (short) fdata[i];
    	    break;
    	case ComplexFloat:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_float -> complex_short\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				newdata[i].re = (short) cfdata[i].re;
				newdata[i].im = (short) cfdata[i].im;
			}
    	    break;
    	case Polar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: polar -> complex_short\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				newdata[i].re = (short) (pdata[i].amp*cos(pdata[i].phi));
				newdata[i].im = (short) (pdata[i].amp*sin(pdata[i].phi));
			}
    	    break;
    	default: break;
    }
    
    p->data = (char *) newdata;
	p->datatype = ComplexShort;
    
    bfree(fdata, datasize*elementsize);    
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_to_complex_int
@Description:
	Converts any data type to a complex integer.
@Algorithm:
	No statistics are calculated.
@Arguments:
	Bimage* p			image (replaced by new image).
@Returns:
	int					error code.
**************************************************************************/
int			img_to_complex_int(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
    if ( p->datatype == ComplexInt ) return(0);
	
	if ( p->colormodel == RGB ) img_RGB2gray(p);
    
    unsigned long		i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	unsigned long		elementsize = gettypesize(p->datatype);
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    complex_short*  	csdata = (complex_short *) p->data;
    complex_float*  	cfdata = (complex_float *) p->data;
    polar*		  		pdata = (polar *) p->data;
    complex_int*  		newdata = (complex_int *) balloc(datasize*sizeof(complex_int));
	if ( !newdata )
		return(error_show("img_to_complex_int", __FILE__, __LINE__));
    
    switch ( p->datatype ) {
    	case UChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned char -> complex_int\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = udata[i];
    	    break;
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> complex_int\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = cdata[i];
    	    break;
    	case UShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned short -> complex_int\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = usdata[i];
    	    break;
    	case Short:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: short -> complex_int\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = sdata[i];
    	    break;
    	case ComplexShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_short -> complex_int\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				newdata[i].re = csdata[i].re;
				newdata[i].im = csdata[i].im;
			}
    	    break;
    	case Int:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: int -> complex_int\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = idata[i];
    	    break;
    	case Float:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: float -> complex_int\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = (int) fdata[i];
    	    break;
    	case ComplexFloat:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_float -> complex_int\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				newdata[i].re = (int) cfdata[i].re;
				newdata[i].im = (int) cfdata[i].im;
			}
    	    break;
    	case Polar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: polar -> complex_int\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				newdata[i].re = (int) (pdata[i].amp*cos(pdata[i].phi));
				newdata[i].im = (int) (pdata[i].amp*sin(pdata[i].phi));
			}
    	    break;
    	default: break;
    }
    
    p->data = (char *) newdata;
	p->datatype = ComplexInt;
    
    bfree(fdata, datasize*elementsize);    
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_to_complex_float
@Description:
	Converts any data type to a complex floating point.
@Algorithm:
	No statistics are calculated.
@Arguments:
	Bimage* p			image (replaced by new image).
@Returns:
	int					error code.
**************************************************************************/
int			img_to_complex_float(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
    if ( p->datatype == ComplexFloat ) return(0);
	
	if ( p->colormodel == RGB ) img_RGB2gray(p);
    
    unsigned long		i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	unsigned long		elementsize = gettypesize(p->datatype);
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    complex_short*  	csdata = (complex_short *) p->data;
    complex_int*  		cidata = (complex_int *) p->data;
    polar*		  		pdata = (polar *) p->data;
    complex_float*  	newdata = (complex_float *) balloc(datasize*sizeof(complex_float));
	if ( !newdata )
		return(error_show("img_to_complex_float", __FILE__, __LINE__));
    
    switch ( p->datatype ) {
    	case UChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned char -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = udata[i];
    	    break;
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = cdata[i];
    	    break;
    	case UShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned short -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = usdata[i];
    	    break;
    	case Short:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: short -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = sdata[i];
    	    break;
    	case ComplexShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_short -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				newdata[i].re = csdata[i].re;
				newdata[i].im = csdata[i].im;
			}
    	    break;
    	case Int:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: int -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = idata[i];
    	    break;
    	case ComplexInt:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_int -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				newdata[i].re = cidata[i].re;
				newdata[i].im = cidata[i].im;
			}
    	    break;
    	case Float:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: float -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = fdata[i];
    	    break;
    	case Polar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: polar -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				newdata[i].re = pdata[i].amp*cos(pdata[i].phi);
				newdata[i].im = pdata[i].amp*sin(pdata[i].phi);
			}
    	    break;
    	default: break;
    }
    
    p->data = (char *) newdata;
	p->datatype = ComplexFloat;
    
    bfree(fdata, datasize*elementsize);    
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_to_polar
@Description:
	Converts any data type to polar.
@Algorithm:
	No statistics are calculated.
@Arguments:
	Bimage* p			image (replaced by new image).
@Returns:
	int					error code.
**************************************************************************/
int			img_to_polar(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
    if ( p->datatype == Polar ) return(0);
	
	if ( p->colormodel == RGB ) img_RGB2gray(p);
    
    unsigned long		i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	double				amp;
	unsigned long		elementsize = gettypesize(p->datatype);
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    complex_short*  	csdata = (complex_short *) p->data;
    complex_int*  		cidata = (complex_int *) p->data;
    complex_float*  	cfdata = (complex_float *) p->data;
    polar*  			newdata = (polar *) balloc(datasize*sizeof(polar));
	if ( !newdata )
		return(error_show("img_to_polar", __FILE__, __LINE__));
    
    switch ( p->datatype ) {
    	case UChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned char -> polar\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].amp = udata[i];
    	    break;
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> polar\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].amp = cdata[i];
    	    break;
    	case UShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned short -> polar\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].amp = usdata[i];
    	    break;
    	case Short:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: short -> polar\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].amp = sdata[i];
    	    break;
    	case ComplexShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_short -> polar\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				amp = csdata[i].re*csdata[i].re + csdata[i].im*csdata[i].im;
				newdata[i].amp = sqrt(amp);
				newdata[i].phi = atan2(1.0*csdata[i].im, 1.0*csdata[i].re);
			}
    	    break;
    	case Int:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: int -> polar\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].amp = idata[i];
    	    break;
    	case ComplexInt:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_int -> polar\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				amp = cidata[i].re*cidata[i].re + cidata[i].im*cidata[i].im;
				newdata[i].amp = sqrt(amp);
				newdata[i].phi = atan2(1.0*cidata[i].im, 1.0*cidata[i].re);
			}
    	    break;
    	case Float:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: float -> polar\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].amp = fdata[i];
    	    break;
    	case ComplexFloat:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_float -> polar\n\n");
    	    for ( i=0; i<datasize; i++ ) {
				newdata[i].amp = sqrt(cfdata[i].re*cfdata[i].re
						+ cfdata[i].im*cfdata[i].im);
				newdata[i].phi = atan2(cfdata[i].im, cfdata[i].re);
			}
    	    break;
    	default: break;
    }
    
    p->data = (char *) newdata;
	p->datatype = Polar;
    
    bfree(fdata, datasize*elementsize);    
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_complex_product
@Description:
	Calculates the complex product of two complex images.
@Algorithm:
	Complex product:
		(a + ib)*(c + id) = (a*c - b*d) + i(b*c + a*d).
	Requirement: The two images must be the same size.
@Arguments:
	Bimage* p			first complex image (replaced by new image).
	Bimage* p2			second complex image.
@Returns:
	int					error code.
**************************************************************************/
int 		img_complex_product(Bimage* p, Bimage* p2)
{
	if ( !p ) return(-1);
	
	img_to_complex_float(p);
	img_to_complex_float(p2);
		
	if ( verbose & VERB_PROCESS )
		printf("Calculating the complex product\n\n");
	
	unsigned long	i, datasize = (unsigned long) p->x*p->y*p->z*p->n*p->c;
	complex_float	temp;
	complex_float*	data = (complex_float *) p->data;
	complex_float*	data2 = (complex_float *) p2->data;
	
	for ( i=0; i<datasize; i++ ) {
		temp = data[i];
		data[i].re = temp.re*data2[i].re - temp.im*data2[i].im;
		data[i].im = data2[i].re*temp.im + temp.re*data2[i].im;
	}
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_complex_apply_mask
@Description:
	Applies a mask to a complex image.
@Algorithm:
	The mask must consist of 0 and 1.
	Requirement: The two images must be the same size.
@Arguments:
	Bimage* p			complex image (replaced by new image).
	Bimage* pmask		mask (converted to floating point).
@Returns:
	int					error code.
**************************************************************************/
int 		img_complex_apply_mask(Bimage* p, Bimage* pmask)
{
	if ( !p ) return(-1);
	if ( !pmask ) return(-1);
	
	img_to_float(pmask);
	
	if ( verbose & VERB_PROCESS )
		printf("Masking a complex image\n");
	
	unsigned long	i, datasize = (unsigned long) p->x*p->y*p->z*p->n*p->c;
	complex_float*	data = (complex_float *) p->data;
	float*			mask = (float *) pmask->data;
	
	for ( i=0; i<datasize; i++ ) {
		data[i].re *= mask[i];
		data[i].im *= mask[i];
	}
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_pack_two_in_complex
@Description:
	Packs two real images into one complex image.
@Algorithm:
	Two real space images are packed into the real and imaginary parts
	of a new data block with conversion from the original non-complex type.
	The header values of the first image are adopted for the new image.
	No statistics are calculated.
@Arguments:
	Bimage* p1		first image.
	Bimage* p2		second image.
@Returns:
	Bimage*			the new complex image.
**************************************************************************/
Bimage*		img_pack_two_in_complex(Bimage* p1, Bimage* p2)
{
	if ( p1->dataflag < 1 ) {
		fprintf(stderr, "Error: No data for %s!\n", p1->filename);
		exit(-1);
	}
	
	if ( p2->dataflag < 1 ) {
		fprintf(stderr, "Error: No data for %s!\n", p2->filename);
		exit(-1);
	}
	
    if ( p1->datatype >= ComplexShort ) {
		fprintf(stderr, "Error: Data for %s may not be complex!\n", p1->filename);
		exit(-1);
	}
	
    if ( p2->datatype >= ComplexShort ) {
		fprintf(stderr, "Error: Data for %s may not be complex!\n", p2->filename);
		exit(-1);
	}
	
	img_check_if_same_size(p1, p2);
		
	if ( p1->colormodel == RGB ) img_RGB2gray(p1);
	if ( p2->colormodel == RGB ) img_RGB2gray(p2);
    
    unsigned long   	i, datasize = p1->x*p1->y*p1->z*p1->n;
    unsigned char* 		udata = (unsigned char *) p1->data;
    signed char* 		cdata = (signed char *) p1->data;
    unsigned short* 	usdata = (unsigned short *) p1->data;
    short* 	    		sdata = (short *) p1->data;
    int* 	    		idata = (int *) p1->data;
    float*  	    	fdata = (float *) p1->data;
    complex_float*  	newdata = (complex_float *) balloc(datasize*sizeof(complex_float));
	if ( !newdata ) {
		error_show("img_pack_two_in_complex", __FILE__, __LINE__);
		return(NULL);
	}
    
	if ( verbose & VERB_LABEL )
    	printf("Packing two real images into a complex image\n");
	
    switch ( p1->datatype ) {
    	case UChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned char -> complex_float\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = udata[i];
    	    break;
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> complex_float\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = cdata[i];
    	    break;
    	case Short:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: short -> complex_float\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = sdata[i];
    	    break;
    	case UShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned short -> complex_float\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = usdata[i];
    	    break;
    	case Int:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: int -> complex_float\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = idata[i];
    	    break;
    	case Float:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: float -> complex_float\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].re = fdata[i];
    	    break;
    	default: break;
    }
    
    udata = (unsigned char *) p2->data;
    cdata = (signed char *) p2->data;
    usdata = (unsigned short *) p2->data;
    sdata = (short *) p2->data;
    idata = (int *) p2->data;
    fdata = (float *) p2->data;
    
    switch ( p2->datatype ) {
    	case UChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned char -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].im = udata[i];
    	    break;
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].im = cdata[i];
    	    break;
    	case Short:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: short -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].im = sdata[i];
    	    break;
    	case UShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned short -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].im = usdata[i];
    	    break;
    	case Int:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: int -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].im = idata[i];
    	    break;
    	case Float:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: float -> complex_float\n\n");
    	    for ( i=0; i<datasize; i++ ) newdata[i].im = fdata[i];
    	    break;
    	default: break;
    }
	
	Bimage* 	p = copy_img_header(p1, -1);
	
    p->data = (char *) newdata;
	p->dataflag = 1;
	p->datatype = ComplexFloat;
	p->transform = NoTransform;
	
//	img_stats(p);

	p->std = p1->std*p2->std;
	
	return(p);
}

/************************************************************************
@Function: img_unpack_combined_transform
@Description:
	Unpacks a complex transform obtained from two real images.
@Algorithm:
	The complex image must be a Fourier transform obtained from two real 
	space images which were packed into the real and imaginary parts
	of the complex image before Fourier transformation.
	The input image is used to hold the transform of the first image.
	A new image is created to hold the transform of the second image.
	Both these images are complex floating point.
	Note: Images with even dimensions cannot be unpacked exactly because
		the inverse is not present when x=nx/2 or y=ny/2 or z=nz/2.
@Arguments:
	Bimage* p		complex Fourier transform from two real space images.
@Returns:
	Bimage*			transform of first image.
**************************************************************************/
Bimage* 	img_unpack_combined_transform(Bimage* p)
{
	int				skip = 0, xmid, ymid, zmid;
	unsigned long   i, j, x, y, z, n;
	long			ix, iy, iz;
	complex_float	d1, d2;
	
	Bimage* 		p2 = copy_img(p);
	
	complex_float*	data = (complex_float *) p->data;
	complex_float*	data2 = (complex_float *) p2->data;
	
	if ( verbose & VERB_PROCESS )
		printf("Unpacking a combined Fourier transform\n\n");
	
	for ( n=0; n<p->n; n++ ) { 
		for ( z=0; z<p->z; z++ ) { 
			iz = -z;
			if ( iz < 0 ) iz += p->z;
			zmid = ( z > 0 && z == iz )? 1: 0;
			for ( y=0; y<p->y; y++ ) { 
				iy = -y;
				if ( iy < 0 ) iy += p->y;
				ymid = ( y > 0 && y == iy )? 1: 0;
				for ( x=0; x<p->x/2+1; x++ ) { 
					ix = -x;
					if ( ix < 0 ) ix += p->x; 
					i = ((n*p->z + z)*p->y + y)*p->x + x;
					j = ((n*p->z + iz)*p->y + iy)*p->x + ix;
					xmid = ( x > 0 && x == ix )? 1: 0;
					if ( xmid+ymid+zmid ) {
						data[i].re = data[i].im = data[j].re = data[j].im = 0;
						data2[i].re = data2[i].im = data2[j].re = data2[j].im = 0;
					} else {
						if ( x == 0 ) {
							if ( y > (p->y - 1)/2 ) skip = 1;
							if ( y == 0 && z > (p->z - 1)/2 ) skip = 1;
						}
						if ( skip ) {
							skip = 0;
						} else {
							d1.re = 0.5*(data[i].re + data[j].re);
							d1.im = 0.5*(data[i].im - data[j].im);
							d2.re = 0.5*(data[i].im + data[j].im);
							d2.im = -0.5*(data[i].re - data[j].re);
							data[i] = data[j] = d1;
							data[j].im = -d1.im;
							data2[i] = data2[j] = d2;
							data2[j].im = -d2.im;
						}
					}
				}
			}
		}
	}
	
	return(p2);
}


/************************************************************************
@Function: img_combined_complex_product
@Description:
	Calculates the complex conjugate product of a complex image resulting
	from combining and Fourier transforming two real space images.
@Algorithm:
	Requirement: Fourier transform of two images packed into one complex
		data block with the function img_pack_two_in_complex and then
		transformed with the function img_fft.
	The Friedel relationships in transforms from real space images are
	exploited to transform two images simultaneously and then extract
	the individual transforms from the complex data set.
	This function extracts the individual transforms and calculates the 
	complex conjugate product used in cross-correlation.
	The result is scaled by the total power of the two transforms within
	the resolution limits, yielding the correlation coefficient when
	the product is backtransformed into the cross-correlation map.
@Arguments:
	Bimage* p			image (replaced by new image).
	float hires			high resolution limit.
	float lores			low resolution limit.
@Returns:
	int					error code.
**************************************************************************/
int 		img_combined_complex_product(Bimage* p, float hires, float lores)
{
	if ( !p ) return(-1);
	
	int 			skip = 0, xeven = (int) (1 - p->x%2);
	unsigned long 	i, j, n, x, y, z;
	long 			xx, yy, zz, ix, iy, iz;
	float			sx2, sy2, sz2, s2;
	
	if ( lores <= 0 ) lores = 1e30;
	if ( hires <= 0 ) hires = 2*p->ux;
	if ( lores < hires ) swap_floats(&lores, &hires);
	float			s2hi = 1/(hires*hires);
	float			s2lo = 1/(lores*lores);
	
	float			xscale = 1.0/(p->x*p->ux);
	float			yscale = 1.0/(p->y*p->uy);
	float			zscale = 1.0/(p->z*p->uz);
	
	unsigned long   datasize = (unsigned long) p->x*p->y*p->z*p->n;

	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		printf("Calculating the complex conjugate product\n");
	
	if ( verbose & VERB_PROCESS )
		printf("Resolution range:               %g - %g A\n", 1/sqrt(s2hi), 1/sqrt(s2lo));
		
	complex_float	temp1, temp2;
	complex_float*	data = (complex_float *) p->data;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_combined_complex_product: F0 = %g %g\n", data[0].re, data[0].im);
	
	double			value, sum = 0;
	
    for ( n=0; n<p->n; n++ ) {
		for ( z=0; z<p->z; z++ ) {
			zz = z;
			if ( z > (p->z - 1)/2 ) zz -= p->z;
			iz = -z;
			if ( iz < 0 ) iz += p->z;
			sz2 = zz*zscale;
			sz2 *= sz2;
			for ( y=0; y<p->y; y++ ) {
				yy = y;
				if ( y > (p->y - 1)/2 ) yy -= p->y;
				iy = -y;
				if ( iy < 0 ) iy += p->y;
				sy2 = yy*yscale;
				sy2 *= sy2;
				for ( x=0; x<p->x/2+1; x++ ) {
					skip = 0;
					if ( x == 0 || ( xeven && x == p->x/2 ) ) {
						if ( yy < 0 ) skip = 1;
						if ( y == 0 && zz < 0 ) skip = 1;
					}
					if ( !skip ) {
						xx = x;
						if ( x > (p->x - 1)/2 ) xx -= p->x;
						ix = -x;
						if ( ix < 0 ) ix += p->x; 
						sx2 = xx*xscale;
						sx2 *= sx2;
						i = ((n*p->z + z)*p->y + y)*p->x + x;
						j = ((n*p->z + iz)*p->y + iy)*p->x + ix;
						s2 = sx2 + sy2 + sz2;
						if ( s2 <= s2hi && s2 >= s2lo ) {
							temp1.re = 0.5*(data[i].re + data[j].re);
							temp1.im = 0.5*(data[i].im - data[j].im);
							temp2.re = 0.5*(data[i].im + data[j].im);
							temp2.im = -0.5*(data[i].re - data[j].re);
							data[i].re = temp1.re*temp2.re + temp1.im*temp2.im;
							data[i].im = temp2.re*temp1.im - temp1.re*temp2.im;
							data[j].re = data[i].re;
							data[j].im = -data[i].im;
							value = (double)data[i].re*data[i].re + (double)data[i].im*data[i].im;
//							if ( !finite(value) ) {
//								printf("Value %ld not finite: %g %g\n", i, data[i].re, data[i].im);
//								exit(-1);
//							}
							if ( value > 0 ) sum += sqrt(value);
						} else {
							data[i].re = 0;
							data[i].im = 0;
							data[j].re = 0;
							data[j].im = 0;
						}
					}
				}
			}
		}
	}

//	printf("sum = %g\n", sum);
	
	double			scale = 0.5/sum;
	
	if ( verbose & VERB_PROCESS )
		printf("Scale:                          %g\n\n", scale);
	
	for ( i=0; i<datasize; i++ ) {
		data[i].re *= scale;
		data[i].im *= scale;
	}
	
	return(0);
}

/************************************************************************
@Function: img_combined_masked_complex_product
@Description:
	Calculates the complex conjugate product within a mask of a complex 
	image resulting from combining and Fourier transforming two real space images.
@Algorithm:
	Requirement: Fourier transform of two images packed into one complex
		data block with the function img_pack_two_in_complex and then
		transformed with the function img_fft.
	The Friedel relationships in transforms from real space images are
	exploited to transform two images simultaneously and then extract
	the individual transforms from the complex data set.
	This function extracts the individual transforms and calculates the 
	complex conjugate product used in cross-correlation.
	The result is scaled by the total power of the two transforms within
	the resolution limits, yielding the correlation coefficient when
	the product is backtransformed into the cross-correlation map.
@Arguments:
	Bimage* p			image (replaced by new image).
	float hires			high resolution limit.
	float lores			low resolution limit.
	Bimage* pmask		reciprocal space mask (converted to floating point).
@Returns:
	int					error code.
**************************************************************************/
int 		img_combined_masked_complex_product(Bimage* p, float hires, 
				float lores, Bimage* pmask)
{
	if ( !p ) return(-1);
	
	if ( pmask ) img_complex_apply_mask(p, pmask);
	
	img_combined_complex_product(p, hires, lores);
	
	return(0);
}

/************************************************************************
@Function: img_complex2real
@Description:
	Gets the real part of a complex floating point image.
@Algorithm:
	The new image contains the real part of the complex image.
	Statistics are calculated for the imaginary part of the complex image
	to assess whether it is sufficiently small - large values may
	indicate improper treatment of the image or transform.
	No statistics are calculated.
@Arguments:
	Bimage* p			image (replaced by new image).
@Returns:
	int					error code.
**************************************************************************/
int 		img_complex2real(Bimage* p)
{
	if ( p->datatype < ComplexShort ) return(0);
	
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
    complex_short*  csdata = (complex_short *) p->data;
    complex_int*  	cidata = (complex_int *) p->data;
    complex_float*  cfdata = (complex_float *) p->data;
    polar*			pdata = (polar *) p->data;
    float*  	    newdata = (float *) balloc(datasize*sizeof(float));
	if ( !newdata )
		return(error_show("img_complex2real", __FILE__, __LINE__));

    float   	    min = 1e30;
    float   	    max = -1e30;
    double  	    imaginary, sum = 0,ssum = 0;
	
    switch ( p->datatype ) {
    	case ComplexShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_short -> float\n\n");
			for ( i=0; i<datasize; i++ ) {
				newdata[i] = csdata[i].re;
				if ( min > csdata[i].im ) min = csdata[i].im;
				if ( max < csdata[i].im ) max = csdata[i].im;
				sum += csdata[i].im;
				ssum += csdata[i].im*csdata[i].im;
 			}
    	    break;
    	case ComplexInt:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_int -> float\n\n");
			for ( i=0; i<datasize; i++ ) {
				newdata[i] = cidata[i].re;
				if ( min > cidata[i].im ) min = cidata[i].im;
				if ( max < cidata[i].im ) max = cidata[i].im;
				sum += cidata[i].im;
				ssum += cidata[i].im*cidata[i].im;
 			}
    	    break;
    	case ComplexFloat:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: complex_float -> float\n\n");
			for ( i=0; i<datasize; i++ ) {
				newdata[i] = cfdata[i].re;
				if ( min > cfdata[i].im ) min = cfdata[i].im;
				if ( max < cfdata[i].im ) max = cfdata[i].im;
				sum += cfdata[i].im;
				ssum += cfdata[i].im*cfdata[i].im;
 			}
    	    break;
    	case Polar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: polar -> float\n\n");
			for ( i=0; i<datasize; i++ ) {
				newdata[i] = pdata[i].amp*cos(pdata[i].phi);
				imaginary = pdata[i].amp*sin(pdata[i].phi);
				if ( min > imaginary ) min = imaginary;
				if ( max < imaginary ) max = imaginary;
				sum += imaginary;
				ssum += imaginary*imaginary;
 			}
    	    break;
    	default: break;
    }
    
    if ( verbose & VERB_PROCESS ) {
    	printf("Imaginary min, max, mean, rms:  %g %g ", min, max);
    	if ( ssum/datasize > (sum/datasize)*(sum/datasize) )
			printf("%g %g\n\n", sum/datasize, sqrt(ssum/datasize - (sum/datasize)*(sum/datasize)));
		else
			printf("%g 0\n\n", sum/datasize);
	}

    bfree(csdata, datasize*gettypesize(p->datatype));
	
    p->data = (char *) newdata;
	p->datatype = Float;
	p->transform = NoTransform;
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_complex2amplitudes
@Description:
	Converts a complex image to an amplitude image.
@Algorithm:
	The new image contains the amplitudes of the complex image.
	No statistics are calculated.
@Arguments:
	Bimage* p			image (replaced by new image).
@Returns:
	int					error code.
**************************************************************************/
int 		img_complex2amplitudes(Bimage* p)
{
	if ( p->datatype < ComplexShort ) return(0);
	
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
    complex_short*  csdata = (complex_short *) p->data;
    complex_int*  	cidata = (complex_int *) p->data;
    complex_float*  cfdata = (complex_float *) p->data;
    polar*			pdata = (polar *) p->data;
    float*  	    newdata = (float *) balloc(datasize*sizeof(float));
	if ( !newdata )
		return(error_show("img_complex2amplitudes", __FILE__, __LINE__));

    switch ( p->datatype ) {
    	case ComplexShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to amplitudes: complex_short -> float\n\n");
			for ( i=0; i<datasize; i++ )
    			newdata[i] = sqrt(1.0*csdata[i].re*csdata[i].re + 
						csdata[i].im*csdata[i].im);
    	    break;
    	case ComplexInt:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to amplitudes: complex_int -> float\n\n");
			for ( i=0; i<datasize; i++ )
    			newdata[i] = sqrt(1.0*cidata[i].re*cidata[i].re + 
						cidata[i].im*cidata[i].im);
    	    break;
    	case ComplexFloat:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to amplitudes: complex_float -> float\n\n");
			for ( i=0; i<datasize; i++ )
    			newdata[i] = sqrt(cfdata[i].re*cfdata[i].re + 
						cfdata[i].im*cfdata[i].im);
    	    break;
    	case Polar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to amplitudes: polar -> float\n\n");
			for ( i=0; i<datasize; i++ )
    			newdata[i] = pdata[i].amp;
    	    break;
    	default: break;
    }
    
    p->data = (char *) newdata;
	
    bfree(csdata, datasize*gettypesize(p->datatype));
	
	p->datatype = Float;
	p->transform = NoTransform;
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_complex2intensities
@Description:
	Converts a complex image to an intensity image.
@Algorithm:
	The new image contains the intensities of the complex image.
	No statistics are calculated.
@Arguments:
	Bimage* p			image (replaced by new image).
@Returns:
	int					error code.
**************************************************************************/
int 		img_complex2intensities(Bimage* p)
{
	if ( p->datatype < ComplexShort ) return(0);
	
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
    complex_short*  csdata = (complex_short *) p->data;
    complex_int*  	cidata = (complex_int *) p->data;
    complex_float*  cfdata = (complex_float *) p->data;
    polar*			pdata = (polar *) p->data;
    float*  	    newdata = (float *) balloc(datasize*sizeof(float));
	if ( !newdata )
		return(error_show("img_complex2intensities", __FILE__, __LINE__));

    switch ( p->datatype ) {
    	case ComplexShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to intensities: complex_short -> float\n\n");
			for ( i=0; i<datasize; i++ )
    			newdata[i] = csdata[i].re*csdata[i].re + 
						csdata[i].im*csdata[i].im;
    	    break;
    	case ComplexInt:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to intensities: complex_int -> float\n\n");
			for ( i=0; i<datasize; i++ )
    			newdata[i] = cidata[i].re*cidata[i].re + 
						cidata[i].im*cidata[i].im;
    	    break;
    	case ComplexFloat:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to intensities: complex_float -> float\n\n");
			for ( i=0; i<datasize; i++ )
    			newdata[i] = cfdata[i].re*cfdata[i].re + 
						cfdata[i].im*cfdata[i].im;
    	    break;
    	case Polar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to intensities: polar -> float\n\n");
			for ( i=0; i<datasize; i++ )
    			newdata[i] = pdata[i].amp*pdata[i].amp;
    	    break;
    	default: break;
    }
    
    bfree(csdata, datasize*gettypesize(p->datatype));
	
    p->data = (char *) newdata;
	p->datatype = Float;
	p->transform = NoTransform;
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_complex2phases
@Description:
	Converts a complex image to a phase image.
@Algorithm:
	The new image contains the phases of the complex image.
	No statistics are calculated.
@Arguments:
	Bimage* p			image (replaced by new image).
@Returns:
	int					error code.
**************************************************************************/
int 		img_complex2phases(Bimage* p)
{
	if ( p->datatype < ComplexShort ) return(0);
	
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
    complex_short*  csdata = (complex_short *) p->data;
    complex_int*  	cidata = (complex_int *) p->data;
    complex_float*  cfdata = (complex_float *) p->data;
    polar*			pdata = (polar *) p->data;
    float*  	    newdata = (float *) balloc(datasize*sizeof(float));
	if ( !newdata )
		return(error_show("img_complex2phases", __FILE__, __LINE__));

    switch ( p->datatype ) {
    	case ComplexShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to phases: complex_short -> float\n\n");
			for ( i=0; i<datasize; i++ )
				newdata[i] = atan2(1.0*csdata[i].im,1.0*csdata[i].re);
    	    break;
    	case ComplexInt:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to phases: complex_int -> float\n\n");
			for ( i=0; i<datasize; i++ )
				newdata[i] = atan2(1.0*cidata[i].im,1.0*cidata[i].re);
    	    break;
    	case ComplexFloat:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to phases: complex_float -> float\n\n");
			for ( i=0; i<datasize; i++ )
				newdata[i] = atan2(cfdata[i].im,cfdata[i].re);
    	    break;
    	case Polar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting to phases: polar -> float\n\n");
			for ( i=0; i<datasize; i++ )
				newdata[i] = pdata[i].phi;
    	    break;
    	default: break;
    }
    
	for ( i=0; i<datasize; i++ ) {
		if ( newdata[i] <= -PI ) newdata[i] += PI;
		if ( newdata[i] >   PI ) newdata[i] -= PI;
	}
	
    bfree(csdata, datasize*gettypesize(p->datatype));
	
    p->data = (char *) newdata;
	p->datatype = Float;
	p->transform = NoTransform;
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_complex2polar
@Description:
	Converts a complex image to a polar image in place.
@Algorithm:
	The new image contains the polar representation of the complex image:
		(re,im) ===> (amp,phi)
	No statistics are calculated.
@Arguments:
	Bimage* p			image (replaced by new image).
@Returns:
	int					error code.
**************************************************************************/
int 		img_complex2polar(Bimage* p)
{
	if ( p->datatype != ComplexFloat ) return(-1);
	
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	polar			temp;
    polar*  	    data = (polar *) p->data;

    if ( verbose & VERB_LABEL ) printf("Converting: complex -> polar\n\n");
    
	// The complex float data type is here represented as polar to
	// make the VMS compiler happy
    for ( i=0; i<datasize; i++ ) {
    	temp.amp = sqrt(data[i].amp*data[i].amp + data[i].phi*data[i].phi);
		temp.phi = atan2(data[i].phi,data[i].amp);
    	data[i] = temp;
		if ( data[i].phi <= -PI ) data[i].phi += PI;
		if ( data[i].phi >   PI ) data[i].phi -= PI;
	}
    
	p->datatype = Polar;
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_polar2complex
@Description:
	Converts a polar image to a complex image in place.
@Algorithm:
	The new image contains the complex representation of the polar image:
		(amp.phi) ===> (re,im)
	No statistics are calculated.
@Arguments:
	Bimage* p			image (replaced by new image).
@Returns:
	int					error code.
**************************************************************************/
int 		img_polar2complex(Bimage* p)
{
	if ( p->datatype != Polar ) return(-1);
	
    unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	complex_float	temp;
    complex_float*  data = (complex_float *) p->data;

    if ( verbose & VERB_LABEL ) printf("Converting: polar -> complex\n\n");
    
	// The polar data type is here represented as complex float to
	// make the VMS compiler happy
    for ( i=0; i<datasize; i++ ) {
		temp.re = data[i].re*cos(data[i].im);
		temp.im = data[i].re*sin(data[i].im);
		data[i] = temp;
	}
    
	p->datatype = ComplexFloat;
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_merge_amplitudes_and_phases
@Description:
	Merges the amplitudes from one map with the phases of another.
@Algorithm:
	The amplitude image can be a floating point image or a complex image.
	The phase image must be complex and its amplitudes are replaced by
	the values from the amplitude image.
	No statistics are calculated.
@Arguments:
	Bimage* pamp		amplitude image (simple or complex).
	Bimage* p			image (replaced by new image).
@Returns:
	double				RMSD of amplitudes.
**************************************************************************/
double		img_merge_amplitudes_and_phases(Bimage* pamp, Bimage* p)
{
	if ( pamp->datatype >= ComplexShort ) {
		if ( pamp->datatype != Polar ) img_to_polar(pamp);
	} else if ( pamp->datatype != Float ) {
		img_to_float(pamp);
	}
	
	img_complex2polar(p);
	
	if ( verbose & VERB_PROCESS )
		printf("Merging amplitudes and phases\n\n");
	
	unsigned long	i, datasize = p->x*p->y*p->z;
	double			d, R = 0;
	polar*			data = (polar *) p->data;
	float*			amp = (float *) pamp->data;
	polar*			ampp = (polar *) pamp->data;
	
	for ( i=0; i<datasize; i++ ) {
		if ( pamp->datatype == Float ) d = data[i].amp - amp[i];
		else d = data[i].amp - ampp[i].amp;
		R += d*d;
		if ( pamp->datatype == Float ) data[i].amp = amp[i];
		else data[i].amp = ampp[i].amp;
	}
	
	R = sqrt(R/datasize);
		
	return(R);
}
