/*
	img_op.cc
	Functions for manipulation of images
	Author: Giovanni Cardone
	Created: 20050314 	Modified:
*/

#include "img_op.h"


/************************************************************************
@Function: img_op_ln
@Description:
	Calculate the natural logarithm of an image.
@Algorithm:
@Arguments:
	Bimage* p		image.
	float   lbase	base value to add before the logarithm.
@Returns:
	int 			0.
**************************************************************************/
int img_op_ln(Bimage* p, float lbase)
{
	if ( !p->dataflag ) return(-1);
	
	img_to_float(p);
	
    float*	fdata = (float *) p->data;

	if ( verbose & VERB_PROCESS )
	    printf("Calculating logarithm of image, after adding %f\n", lbase);
	
	unsigned long   i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	
	for ( i=0; i<datasize; i++ ) {
    		if (fdata[i]+lbase <= 0.) {
    				printf("(img_op_ln): negative values found. Exiting...");
    				return(-1);
    		}
    		fdata[i] = log(fdata[i] + lbase);
	}
     
    img_stats(p);
      
	return(0);
}
