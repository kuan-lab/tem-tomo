/*
	img_combine.c
    Functions to combine two images in various ways
	Author: Bernard Heymann
	Created: 19990219 	Modified: 20041225
*/

#include "rwimg.h"
#include "img_combine.h"
#include "img_rescale.h"
#include "img_util.h"
#include "img_datatypes.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;		// Total memory allocated 

/************************************************************************
@Function: img_check_if_same_size
@Description:
	Check if two images are the same size.
@Algorithm:
	The function returns the answer and leaves the calling function to 
	deal with the result.
@Arguments:
	Bimage* p1			first image.
	Bimage* p2			second image.
@Returns:
	int 				0 if yes, -1 if no.
**************************************************************************/
int 		img_check_if_same_size(Bimage* p1, Bimage* p2)
{
	if ( p1->x != p2->x || p1->y != p2->y || p1->z != p2->z ||
			p1->c != p2->c || p1->n != p2->n ) {
		fprintf(stderr, "Image dimensions not the same:\n");
		fprintf(stderr, "Image 1 (c, x, y, z, n):        %ld %ld %ld %ld %ld\n",
				p1->c, p1->x, p1->y, p1->z, p1->n);
		fprintf(stderr, "Image 2 (c, x, y, z, n):        %ld %ld %ld %ld %ld\n",
				p2->c, p2->x, p2->y, p2->z, p2->n);
		return(-1);
	}
    
	return(0);
}

/************************************************************************
@Function: img_compatibility
@Description:
	Check if two images have the same number of channels and data type.
@Algorithm:
	The function returns the answer and leaves the calling function to 
	deal with the result.
@Arguments:
	Bimage* p1			first image.
	Bimage* p2			second image.
@Returns:
	int 				0 if yes, <0 if no.
**************************************************************************/
int 		img_compatibility(Bimage* p1, Bimage* p2)
{
	int 	non = 0;

	if ( p1->c != p2->c ) {
		printf("%s has an incompatible number of channels:  %ld\n", 
				p2->filename, p2->c);
		non += -1;
	}
	
	if ( p1->datatype != p2->datatype ) {
		printf("%s has an incompatible data type:   %d\n", 
				p2->filename, p2->datatype);
		non += -1;
	}
	
	return(non);
}

/************************************************************************
@Function: add_images
@Description:
	Adds two images together.
@Algorithm:
	The second image is added to the first:
		image1 = image1 + image2*scale + shift
	Both images are converted to floating point.
@Arguments:
	Bimage* p1			first image (modified).
	Bimage* p2			second image.
	float scale 		density scale to apply to second image
	float shift 		density shift to apply to second image.
@Returns:
	int 				0, <0 if error.
**************************************************************************/
int 		add_images(Bimage* p1, Bimage* p2, float scale, float shift)
{
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("add_images", __FILE__, __LINE__);
		return(-1);
	}
	
	img_to_float(p1);
	img_to_float(p2);
	
    float*  	    data1 = (float *) p1->data;
    float*  	    data2 = (float *) p2->data;
    
    unsigned long   i, datasize = p1->n*p1->x*p1->y*p1->z*p1->c;
	
	if ( p1->datatype == ComplexFloat ) datasize *= 2;
    
    if ( verbose & VERB_LABEL )
		printf("Adding/subtracting two images\n\n");

    for ( i=0; i<datasize; i++ )
    	data1[i] += data2[i]*scale + shift;

    for ( i=0; i<p1->n; i++ )
		p1->image[i].background += p2->image[i].background*scale + shift;
    
	return(0);
}

/************************************************************************
@Function: img_add_weighed
@Description:
	Adds multiple images together with given weights.
@Algorithm:
	Images are read from a number files and added to each other, using
	the given weights to determine each contribution.
	The images are rescaled to a new average and standard deviation before 
	weighted addition. If the given standard deviation is zero or less,
	this step is omitted.
	The weighed average and standard deviation is calculated, the latter
	written into the FOM block.
	All images are converted to floating point.
@Arguments:
	int nfiles			number of files to add together.
	char* file_list 	list of file names.
	float* weight		list of weights.
	float newavg		new average for rescaling.
	float newstd		new standard deviation for rescaling.
@Returns:
	Bimage* 			resultant image (floating point).
**************************************************************************/
Bimage* 	img_add_weighed(int nfiles, char* file_list, float* weight,
					float newavg, float newstd)
{
	unsigned long   i;
	int             n;
	char*			filename = file_list;
	Bimage*			p = NULL;
	Bimage*			psum = init_img();

	for ( n=0; n<nfiles; n++ ) {
		p = read_img(filename, 0, -1);
		if ( p != NULL ) {
			if ( psum->c < p->c ) psum->c = p->c;
			if ( psum->x < p->x ) psum->x = p->x;
			if ( psum->y < p->y ) psum->y = p->y;
			if ( psum->z < p->z ) psum->z = p->z;
			if ( psum->n < p->n ) psum->n = p->n;
			if ( n == 0 ) {
				if ( p->ux > 0 && p->ux != 1 ) psum->ux = p->ux;
				if ( p->uy > 0 && p->uy != 1 ) psum->uy = p->uy;
				if ( p->uz > 0 && p->uz != 1 ) psum->uz = p->uz;
			}
			kill_img(p);
		}
		filename += strlen(filename) + 1;
	}
	
	if ( verbose & VERB_PROCESS ) {
		printf("\nAdding %d images together\n", nfiles);
		printf("New image size:                 %ld %ld %ld\n", psum->x, psum->y, psum->z);
	} else if ( verbose & VERB_LABEL )
		printf("\nAdding %d images together\n\n", nfiles);
	
	psum->dataflag = 1;
	psum->datatype = Float;
	psum->data = (char *) balloc(psum->x*psum->y*psum->z*psum->n*sizeof(float));
	psum->image = (Bsub_image *) balloc(psum->n*sizeof(Bsub_image));
	float*			sumdata = (float *) psum->data;
	float*			data;
	
	double			weightsum = 0;
	unsigned long   datasize = psum->x*psum->y*psum->z*psum->n;
	float*			fom = img_init_fom(psum);
	
	filename = file_list;
	for ( n=0; n<nfiles; n++ ) {
		p = read_img(filename, 1, -1);
		if ( p != NULL ) {
			if ( verbose & VERB_LABEL )
				printf("Adding image %d with weight %g\n", n, weight[n]);
			weightsum += weight[n];
			if ( n == 0 )
				memcpy(psum->image, p->image, p->n*sizeof(Bsub_image));
			img_to_float(p);
			if ( newstd > 0 ) img_rescale_to_avg_std(p, newavg, newstd);
			data = (float *) p->data;
			for ( i=0; i<datasize; i++ ) {
				sumdata[i] += weight[n]*data[i];
				fom[i] += weight[n]*data[i]*data[i];
			}
			kill_img(p);
		}
		filename += strlen(filename) + 1;
	}
	
	for ( i=0; i<datasize; i++ ) {
		sumdata[i] /= weightsum;
		fom[i] = fom[i]/weightsum - sumdata[i]*sumdata[i];
		if ( fom[i] > 0 ) fom[i] = sqrt(fom[i]);
		else fom[i] = 0;
	}
	
	return(psum); 
}

/************************************************************************
@Function: img_place
@Description:
	Places a small image into a large image.
@Algorithm:
	The small image is placed with its origin at given origin in the large image.
	The second image is scaled and shifted before placing into the first:
		image1 = image1 + image2*scale + shift
	Both images are converted to floating point.
@Arguments:
	Bimage* p			large image (destination, modified).
	Bimage* psmall		small image (source).
	Vector3 origin		location in large image of small image origin.
	float radius		radial mask to transfer small image.
	float scale 		density scale to apply to second image
	float shift 		density shift to apply to second image.
@Returns:
	int 				0, <0 if error.
**************************************************************************/
int 		img_place(Bimage* p, Bimage* psmall, Vector3 origin, float radius, float scale, float shift)
{
	if ( p->n > 1 ) {
		error_show("Only single gray scale images can be used for placement!\n", __FILE__, __LINE__);
		return(-1);
	}
	
	img_to_float(p);
	img_to_float(psmall);

    unsigned long   i, si, x, y, z, sz, sy, sx;
	double			sx2, sy2, sz2, r2 = radius*radius;
	
	VectorInt3		lo = {(int) (origin.x - psmall->image->ox), 
							(int) (origin.y - psmall->image->oy), 
							(int) (origin.z - psmall->image->oz)};
	lo = vectorint3_scalar_max(lo, 0);
	
	VectorInt3		hi = {(int) (lo.x+psmall->x), (int) (lo.y+psmall->y), (int) (lo.z+psmall->z)};
	if ( hi.x > p->x ) hi.x = (int) (p->x);
	if ( hi.y > p->y ) hi.y = (int) (p->y);
	if ( hi.z > p->z ) hi.z = (int) (p->z);
	
    float*  	    data = (float *) p->data;
    float*  	    small = (float *) psmall->data;
    
    if ( verbose & VERB_LABEL )
		printf("Placing a small image into a larger image\n\n");
		
	for ( z=lo.z; z<hi.z; z++ ) {
		sz = (unsigned long) (psmall->image->oz + z - origin.z);
		sz2 = sz - psmall->image->oz;
		sz2 *= sz2;
		if ( sz >= p->z ) printf("sz = %ld\n", sz);
		for ( y=lo.y; y<hi.y; y++ ) {
			sy = (unsigned long) (psmall->image->oy + y - origin.y);
			sy2 = sy - psmall->image->oy;
			sy2 *= sy2;
			if ( sy >= p->y ) printf("sy = %ld\n", sy);
			for ( x=lo.x; x<hi.x; x++ ) {
				sx = (unsigned long) (psmall->image->ox + x - origin.x);
				sx2 = sx - psmall->image->ox;
				sx2 *= sx2;
				if ( sx >= p->x ) printf("sx = %ld\n", sx);
				if ( sx2 + sy2 + sz2 <= r2 ) {
					i = (z*p->y + y)*p->x + x;
					si = (sz*psmall->y + sy)*psmall->x + sx;
					if ( small[si] < data[i] )
						data[i] = small[si];
				}
			}
		}
	}

	return(0); 
}

/************************************************************************
@Function: multiply_images
@Description:
	Multiplies two images.
@Algorithm:
	The second image is multiplied with the first:
		image1 = image1 * (image2*scale + shift)
	Both images are converted to floating point.
@Arguments:
	Bimage* p1			first image (modified).
	Bimage* p2			second image.
	float scale 		density scale to apply to second image
	float shift 		density shift to apply to second image.
@Returns:
	int 				0, <0 if error.
**************************************************************************/
int 		multiply_images(Bimage* p1, Bimage* p2, float scale, float shift)
{
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("multiply_images", __FILE__, __LINE__);
		return(-1);
	}
	
	img_to_float(p1);
	img_to_float(p2);
	
    float*  	    data1 = (float *) p1->data;
    float*  	    data2 = (float *) p2->data;
    
    unsigned long   i, datasize = (unsigned long) p1->n*p1->x*p1->y*p1->z*p1->c;
	float			avprod = p1->avg*(p2->avg*scale+shift);
    
    if ( verbose & VERB_LABEL )
	    printf("Multiplying two images\n\n");

    for ( i=0; i<datasize; i++ )
    	data1[i] *= data2[i]*scale + shift;
	
	img_stats(p1);
    
	if ( verbose )
		printf("Correspondence:                 %g\n\n", p1->avg/avprod);
	
	return(0);
}

/************************************************************************
@Function: largest_img
@Description:
	Selects the largest of each pixel from two images.
@Algorithm:
	Both images are converted to floating point.
@Arguments:
	Bimage* p1			first image.
	Bimage* p2			second image.
	float scale 		density scale to apply to second image
	float shift 		density shift to apply to second image.
@Returns:
	int 				0, <0 if error.
**************************************************************************/
int 		largest_img(Bimage* p1, Bimage* p2, float scale, float shift)
{
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("largest_img", __FILE__, __LINE__);
		return(-1);
	}
	
	img_to_float(p1);
	img_to_float(p2);
	
    float*  	    data1 = (float *) p1->data;
    float*  	    data2 = (float *) p2->data;
    
    unsigned long   i, datasize = (unsigned long) p1->n*p1->x*p1->y*p1->z*p1->c;
    
    if ( verbose & VERB_LABEL )
	    printf("Largest of two images\n\n");

    for ( i=0; i<datasize; i++ )
    	if ( data1[i] < data2[i]*scale + shift )
			data1[i] = data2[i]*scale + shift;
    
	return(0);
}

/************************************************************************
@Function: smallest_img
@Description:
	Selects the smallest of each pixel from two images.
@Algorithm:
	Both images are converted to floating point.
@Arguments:
	Bimage* p1			first image.
	Bimage* p2			second image.
	float scale 		density scale to apply to second image
	float shift 		density shift to apply to second image.
@Returns:
	int 				0, <0 if error.
**************************************************************************/
int 		smallest_img(Bimage* p1, Bimage* p2, float scale, float shift)
{
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("smallest_img", __FILE__, __LINE__);
		return(-1);
	}
	
	img_to_float(p1);
	img_to_float(p2);
	
    float*  	    data1 = (float *) p1->data;
    float*  	    data2 = (float *) p2->data;
    
    unsigned long   i, datasize = (unsigned long) p1->n*p1->x*p1->y*p1->z*p1->c;
    
    if ( verbose & VERB_LABEL )
	    printf("Smallest of two images\n\n");

    for ( i=0; i<datasize; i++ )
    	if ( data1[i] > data2[i]*scale + shift )
			data1[i] = data2[i]*scale + shift;
    
	return(0);
}

/************************************************************************
@Function: img_linear_fit
@Description:
	Linear least squares fit of two images.
@Algorithm:
	The data blocks from two images are fit by a simple linear least squares
	regression algorithm with exclusion of a percentage of outliers:
		image2 = intercept + slope * image1
	The first image is modified to return the difference:
		image1_new = intercept + slope * image1 - image2
	The two data blocks must have the same size and are converted to
	floating point.
	Note: A linear fit is not symmetric with respect to the two input
	data sets - the order of the input images determine the output.
	Both images are converted to floating point.
@Arguments:
	Bimage* p1			first image (modified).
	Bimage* p2			second image.
	Bimage* pmask		mask to limit calculation to a certain region.
	float max_exclude	maximum percentage of outlying points to exclude.
@Returns:
	float* 				array of R factors, NULL if not run.
**************************************************************************/
float*		img_linear_fit(Bimage* p1, Bimage* p2, Bimage* pmask, float max_exclude)
{
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("img_linear_fit", __FILE__, __LINE__);
		return(NULL);
	}
	
	img_to_float(p1);
	img_to_float(p2);
	
	if ( max_exclude >= 1 ) max_exclude /= 100;	// Assume it is a percentage
	if ( max_exclude > 0.5 ) max_exclude = 0.5; // At least half of the data must be fit
	
    unsigned long   i, j, n, nm, nfit, best_excl = 0;
    unsigned long   imgsize = p1->x*p1->y*p1->z*p1->c;
//    unsigned long   datasize = p1->n*imgsize;
	unsigned long   num = imgsize;
    int     	    bin;
    double	  		sx, sx2, sy, sxy, sd, dy, denominator, a, b, R, best_a = 0, best_b = 0;
	float			excl_voxels, diff, maxdiff, diff_cutoff;
    
    float*  	    data1 = (float *) p1->data;
    float*  	    data2 = (float *) p2->data;
	int*			hist = (int *) balloc(1000*sizeof(int));
	int*			inc_mask = (int *) balloc(imgsize*sizeof(int));
	for ( i=0; i<imgsize; i++) inc_mask[i] = 1;
	float*			bestR = (float *) balloc(p1->n*sizeof(float));
    
	unsigned char*  mask = NULL;
	if ( pmask ) {
		img_to_byte(pmask);
		mask = (unsigned char *) pmask->data;
	}
    
    if ( verbose & VERB_LABEL ) {
	    printf("Linear least squares fit:       image2 = intercept + slope x image1\n");
	    printf("Maximum voxels to exclude:      %g %%\n", 100.0*max_exclude);
	}

	if ( verbose & VERB_FULL )
		printf(" %%Set\t Voxels  \t%%Excluded\tIntercept\tSlope\t\tR\n");
	else if ( verbose & VERB_PROCESS )
		printf("Intercept\t  Slope\t  R-factor\t%%Excluded\n");
	
	for ( n=0; n<p1->n; n++) {
		best_excl = 0;
		best_a = best_b = 0;
		bestR[n] = 1e30;
		if ( mask ) {
			nm = n;
			if ( pmask->n < p1->n ) nm = 0;
			for ( num=j=0, i=nm*imgsize; j<imgsize; j++,i++ ) {
				if ( mask[i] ) {
					inc_mask[i] = 1;
					num++;
				} else inc_mask[j] = -1;
			}
		}
		for ( excl_voxels = 0; excl_voxels < max_exclude+0.01; excl_voxels += 0.01 ) {
			nfit = 0;
			a = b = sx = sx2 = sy = sxy = sd = dy = maxdiff = 0;
			for ( j=0, i=n*imgsize; j<imgsize; j++,i++ ) if ( inc_mask[j] > 0 ) {
				sx  += data1[i];
				sx2 += data1[i]*data1[i];
				sy  += data2[i];
				sxy += data1[i]*data2[i];
				nfit++;
			}
			denominator = nfit*sx2 - sx*sx;
			if ( denominator ) {
				a = (sx2*sy   - sx*sxy)/denominator;
				b = (nfit*sxy - sx*sy) /denominator;
			}
			sy /= nfit;
    
			for ( j=0, i=n*imgsize; j<imgsize; j++,i++ ) if ( inc_mask[j] >= 0 ) {
				diff = a + b*data1[i] - data2[i];
				if ( inc_mask[j] ) {
					dy += diff*diff;
					sd += (data2[i]-sy)*(data2[i]-sy);
				}
				diff = fabs(diff);
				if ( maxdiff < diff ) maxdiff = diff;
			}
			R = 1e30;
			if ( dy < 1e-30 ) R = 0;
			else if ( sd ) R = sqrt(dy/sd);
			if ( bestR[n] > R ) {
				bestR[n] = R;
				best_a = a;
				best_b = b;
				best_excl = num - nfit;
			}
		
			if ( verbose & VERB_FULL )
				printf("%6.2f\t%12ld\t%6.2f\t\t%9g\t%9g\t%g\n", 100*excl_voxels, 
						nfit, 100.0-nfit*100.0/num, a, b, R);
	
			if ( verbose & VERB_DEBUG )
				printf("DEBUG img_linear_fit: maxdiff = %g\n", maxdiff);
			if ( maxdiff < 1e-30 ) maxdiff = 1;
			memset(hist, 0, 1000*sizeof(int));
			for ( j=0, i=n*imgsize; j<imgsize; j++,i++ ) if ( inc_mask[j] >= 0 ) {
				diff = fabs(a + b*data1[i] - data2[i]);
				bin = (int) (1000*diff/maxdiff);
				if ( bin > 999 ) bin = 999;
				hist[bin]++;
			}
		
			i = 999;
			j = hist[i];
			while ( i > 1 && j < num*excl_voxels ) {
				i--;
				j += hist[i];
			}
			diff_cutoff = i*maxdiff/1000.0;
//			printf("excl_voxels=%7.3f i=%d j=%d maxdiff=%g diff_cutoff=%g\n", excl_voxels, i, j, maxdiff, diff_cutoff);
		
			for ( j=0, i=n*imgsize; j<imgsize; j++,i++ ) if ( inc_mask[j] >= 0 ) {
				diff = fabs(a + b*data1[i] - data2[i]);
				inc_mask[j] = 1;
				if ( diff > diff_cutoff ) inc_mask[j] = 0;
			}
		}
	
		for ( j=0, i=n*imgsize; j<imgsize; j++,i++ ) {
			if ( inc_mask[j] >= 0 )
				data1[i] = best_a + best_b*data1[i] - data2[i];
			else
				data1[i] = 0;
		}
		if ( verbose & VERB_PROCESS )
			printf("%g\t%g\t%g\t%ld (%5.2f %%)\n", best_a, best_b, bestR[n], best_excl, best_excl*100.0/num);
	}
	
	bfree(hist, 1000*sizeof(int));
	bfree(inc_mask, imgsize*sizeof(int));
	
//	if ( verbose & VERB_PROCESS ) {
//		printf("Best linear fit:                image2 = %g + %g x image1\n", best_a, best_b);
//		printf("R factor:                       %g\n", bestR);
//		printf("Voxels excluded:                %ld (%5.2f %%)\n", 
//				best_excl, best_excl*100.0/num);
//	}
	if ( verbose >= VERB_PROCESS ) printf("\n");
    
	return(bestR);
}

/************************************************************************
@Function: img_R_factor
@Description:
	Calculates an R factor between two images.
@Algorithm:
	The difference between two images is calculated and normalized as:
		                      sum(image1 - image2)^2
		R = sqrt(-------------------------------------------------)
		         sqrt(sum(image1 - avg1)^2 * sum(image2 - avg2)^2)
	Both images are converted to floating point.
@Arguments:
	Bimage* p1			first image.
	Bimage* p2			second image.
@Returns:
	float 				R factor, -1 if not run.
**************************************************************************/
float		img_R_factor(Bimage* p1, Bimage* p2)
{
	if ( p1->c > 1 || p2->c > 1 ) {
		fprintf(stderr, "Error: Only single channel images are supported for R factor calculation!\n");
		return(-1.0);
	}
	
	if ( p1->n > 1 || p2->n > 1 ) {
		fprintf(stderr, "Error: Only single images are supported for R factor calculation!\n");
		return(-1.0);
	}
	
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("img_R_factor", __FILE__, __LINE__);
		return(-1);
	}
	
	img_to_float(p1);
	img_to_float(p2);
	
    unsigned long   i, datasize = (unsigned long) p1->x*p1->y*p1->z;
	float			R = 0;
	double			sx = 0, sy = 0, sx2 = 0, sy2 = 0, sxy = 0;
	float*			data1 = (float *) p1->data;
	float*			data2 = (float *) p2->data;
	
	for ( i=0; i<datasize; i++ ) {
		sx += data1[i];
		sy += data2[i];
		sx2 += data1[i]*data1[i];
		sy2 += data2[i]*data2[i];
		sxy += data1[i]*data2[i];
	}
	
	R = sqrt((sx2 - 2*sxy + sy2)/sqrt((sx2 - sx*sx/datasize)*(sy2 - sy*sy/datasize)));
	
	return(R);
}

/************************************************************************
@Function: img_one_correlation_coefficient
@Description:
	Calculates correlation coefficient between two images.
@Algorithm:
	The correlation between two images is calculated and normalized as:
		           sum((image1 - avg1)*(image2 - avg2))
		CC = -------------------------------------------------
		     sqrt(sum(image1 - avg1)^2 * sum(image2 - avg2)^2)
	.
	Both images are converted to floating point.
	Only the first image is used.
@Arguments:
	Bimage* p1			first image.
	Bimage* p2			second image.
@Returns:
	double 				correlation coefficient, -1 if not run.
**************************************************************************/
double		img_one_correlation_coefficient(Bimage* p1, Bimage* p2)
{
	if ( p1->c > 1 || p2->c > 1 ) {
		fprintf(stderr, "Error: Only single channel images are supported for correlation!\n");
		return(-1.0);
	}
	
	if ( p1->n > 1 || p2->n > 1 ) {
		fprintf(stderr, "Error: Only single images are supported for correlation!\n");
		return(-1.0);
	}
	
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("img_one_correlation_coefficient", __FILE__, __LINE__);
		return(-1);
	}
	
	img_to_float(p1);
	img_to_float(p2);
	
    unsigned long   i, datasize = p1->x*p1->y*p1->z;
	double			cc = 0;
	double			sx = 0, sy = 0, sx2 = 0, sy2 = 0, sxy = 0;
	float*			data1 = (float *) p1->data;
	float*			data2 = (float *) p2->data;
	
	for ( i=0; i<datasize; i++ ) {
		sx += data1[i];
		sy += data2[i];
		sx2 += data1[i]*data1[i];
		sy2 += data2[i]*data2[i];
		sxy += data1[i]*data2[i];
	}
	
	double			varx = sx2 - sx*sx/datasize;
	double			vary = sy2 - sy*sy/datasize;
	
	if ( varx > 0 && vary > 0 )
		cc = (sxy - sx*sy/datasize)/sqrt(varx*vary);
	
	return(cc);
}

/************************************************************************
@Function: img_correlation_coefficient
@Description:
	Calculates correlation coefficient between two images.
@Algorithm:
	The correlation between two images is calculated and normalized as:
		           sum((image1 - avg1)*(image2 - avg2))
		CC = -------------------------------------------------
		     sqrt(sum(image1 - avg1)^2 * sum(image2 - avg2)^2)
	.
	Both images are converted to floating point.
@Arguments:
	Bimage* p1			first image.
	Bimage* p2			second image.
@Returns:
	float* 				array of correlation coefficients, NULL if error.
**************************************************************************/
float*		img_correlation_coefficient(Bimage* p1, Bimage* p2)
{
	if ( p1->c > 1 || p2->c > 1 ) {
		fprintf(stderr, "Error: Only single channel images are supported for correlation!\n");
		return(NULL);
	}
	
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("img_correlation_coefficient", __FILE__, __LINE__);
		return(NULL);
	}
	
	img_to_float(p1);
	img_to_float(p2);
	
    unsigned long   i, n, datasize = (unsigned long) p1->x*p1->y*p1->z;
	double			sx = 0, sy = 0, sx2 = 0, sy2 = 0, sxy = 0, v1, v2;
	float			*data1, *data2;
	float*			cc = (float *) balloc(p1->n*sizeof(float));
	
	for ( n=0; n<p1->n; n++ ) {
		data1 = (float *) (p1->data + n*datasize*sizeof(float));
		data2 = (float *) (p2->data + n*datasize*sizeof(float));
		sx = sy = sx2 = sy2 = sxy = 0;
		for ( i=0; i<datasize; i++ ) {
			sx += data1[i];
			sy += data2[i];
			sx2 += data1[i]*data1[i];
			sy2 += data2[i]*data2[i];
			sxy += data1[i]*data2[i];
		}
		v1 = sx2 - sx*sx/datasize;
		v2 = sy2 - sy*sy/datasize;
		cc[n] = (sxy - sx*sy/datasize)/sqrt(v1*v2);
//		printf("%g %g %g\n", v1, v2, cc[n]);
	}
	
	return(cc);
}

/************************************************************************
@Function: img_correlation_coefficient_radial
@Author: David Belnap
@Description:
	Calculates the correlation coefficient between two images and 
	applies a radial filter.
@Algorithm:
	Calculates correlation coefficient between two images using the 
	equation:

                      Sum(Ai*Bi) - [Sum(Ai)*Sum(Bi)]/n
 CC = -----------------------------------------------------------------
      SQRT[Sum(Ai^2) - (Sum(Ai))^2/n] * SQRT[Sum(Bi^2) - (Sum(Bi))^2/n]

	where Ai and Bi are pixel values of images A and B, respectively,
	and n is the total number of pixels sampled.
	   A radial filter is applied.  The center of the image is needed and
	is read from the pA Bimage structure.  The two images must have the 
	same center (apply a translation shift if necessary).  The program 
	tests that each pixel is within the radii range.  If the pixel is
	within the range, then it is included in the CC calculation.
@Arguments:
	Bimage* pA           first image.
	Bimage* pB           second image.
	float   radius_min   minimum radius for CC calculation, in pixels
	float   radius_max   maximum radius for CC calculation, in pixels
@Returns:
	float				correlation coefficient.
**************************************************************************/
float		img_correlation_coefficient_radial(Bimage* pA, Bimage* pB, float radius_min, float radius_max)
{
	img_to_float(pA);   // converts to floating point if it isn't already
	img_to_float(pB);
	unsigned int     i, x, y, z;
	long    n = 0;
	float   cc = 0;
	float   xdiff, ydiff, zdiff, xsq, ysq, zsq;
	float   radminsq, radmaxsq, radsq;
	float   xcen, ycen, zcen;
	double  sumA = 0, sumB = 0, sumAsq = 0, sumBsq = 0, sumAB = 0;
	float*  dataA = (float *) pA->data;
	float*  dataB = (float *) pB->data;

	xcen = pA->image->ox;
	ycen = pA->image->oy;
	zcen = pA->image->oz;

	if ( pA->c > 1 || pB->c > 1 ) {
		fprintf(stderr, "Error: Only single channel images are supported for correlation!\n");
		return(-1.0);
	}
	if ( pA->n > 1 || pB->n > 1 ) {
		fprintf(stderr, "Error: Only single images are supported for correlation!\n");
		return(-1.0);
	}
	if ( (xcen != pB->image->ox) || (ycen != pB->image->oy) || (zcen != pB->image->oz) )  {
		printf("Centers of the two images are not the same.\n");
		printf("center of image 1:  %f %f %f\n",xcen,ycen,zcen);
		printf("center of image 2:  %f %f %f\n",pB->image->ox,pB->image->oy,pB->image->oz);
		printf("Align images or set centers before using function ""img_correlation_coefficient_radial"".\n");
		return(-1.0);
	}

	radminsq = radius_min*radius_min;
	radmaxsq = radius_max*radius_max;

	for (z=0; z < pA->z; z++)  {
		zdiff = z-zcen;
		zsq   = zdiff*zdiff;
		for (y=0; y < pA->y; y++)  {
			ydiff = y-ycen;
			ysq   = ydiff*ydiff;
			for (x=0; x < pA->x; x++)  {
				xdiff = x-xcen;
				xsq   = xdiff*xdiff;
				radsq = xsq + ysq + zsq;
				if ( (radsq >= radminsq) && (radsq <= radmaxsq) )  {
					i = (unsigned int) (x + (pA->x)*y + (pA->x)*(pA->y)*z);
					sumA   += dataA[i];
					sumB   += dataB[i];
					sumAB  += dataA[i] * dataB[i];
					sumAsq += dataA[i]*dataA[i];
					sumBsq += dataB[i]*dataB[i];
					n      += 1;
				}
			}
		}
	}

	cc = (sumAB - sumA*sumB/n) / (sqrt(sumAsq - sumA*sumA/n) * sqrt(sumBsq - sumB*sumB/n));
	
	return(cc);
}

/************************************************************************
@Function: correlate_images
@Description:
	Calculates a correlation image and a correlation coefficient between 
	two images.
@Algorithm:
	The correlation between two images is calculated and normalized as:
		         sum((image1 - basis1)*(image2 - basis2))
		CC = -------------------------------------------------
		     sqrt(sum(image1 - avg1)^2 * sum(image2 - avg2)^2)
	The basis values can be anything, but is typically either the
	average or the background in an image.
	Both images are converted to floating point.
@Arguments:
	Bimage* p1			first image (modified).
	Bimage* p2			second image.
	float bas1			basis value for first image.
	float bas2			basis value for second image.
@Returns:
	double 				correlation coefficient, -1 if error.
**************************************************************************/
double		correlate_images(Bimage* p1, Bimage* p2, float bas1, float bas2)
{
	if ( p1->c > 1 || p2->c > 1 ) {
		fprintf(stderr, "Error: Only single channel images are supported for correlation!\n");
		return(-1.0);
	}
	
	if ( p1->n > 1 || p2->n > 1 ) {
		fprintf(stderr, "Error: Only single images are supported for correlation!\n");
		return(-1.0);
	}
	
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("correlate_images", __FILE__, __LINE__);
		return(-1);
	}
	
	if ( verbose & VERB_PROCESS )
		printf("Basis values for correlation:   %g %g\n", bas1, bas2);

	img_to_float(p1);
	img_to_float(p2);
	
	if ( p1->std == 0 ) img_stats(p1);
	if ( p2->std == 0 ) img_stats(p2);
	
    float*  	    data1 = (float *) p1->data;
    float*  	    data2 = (float *) p2->data;
    
    unsigned long   i, datasize = (unsigned long) p1->x*p1->y*p1->z;
    double  		cc = 0;
    
    for ( i=0; i<datasize; i++ ) {
		data1[i] = (data1[i] - bas1)*(data2[i] - bas2)/(p1->std*p2->std);
		cc += data1[i];
	}
	cc /= datasize;
	
	img_stats(p1);
    
    if ( verbose & VERB_RESULT )
		printf("Correlation:                    %g (%ld)\n\n", cc, datasize);
	
	return(cc);
}

/************************************************************************
@Function: correlate_images_within_radii
@Description:
	Calculates a correlation image and a correlation coefficient between 
	two images, only within the given radii.
@Algorithm:
	The correlation between two images is calculated and normalized as:
		         sum((image1 - basis1)*(image2 - basis2))
		CC = -------------------------------------------------
		     sqrt(sum(image1 - avg1)^2 * sum(image2 - avg2)^2)
	The basis values can be anything, but is typically either the
	average or the background in an image.
	Both images are converted to floating point.
@Arguments:
	Bimage* p1			first image (modified).
	Bimage* p2			second image.
	float bas1			basis value for first image.
	float bas2			basis value for second image.
	float minr			minimum radius (pixel units).
	float maxr			maximum radius (pixel units).
@Returns:
	double 				correlation coefficient, -1 if error.
**************************************************************************/
double		correlate_images_within_radii(Bimage* p1, Bimage* p2, 
				float bas1, float bas2, float minr, float maxr)
{
	if ( p1->c > 1 || p2->c > 1 ) {
		fprintf(stderr, "Error: Only single channel images are supported for correlation!\n");
		return(-1.0);
	}
	
	if ( p1->n > 1 || p2->n > 1 ) {
		fprintf(stderr, "Error: Only single images are supported for correlation!\n");
		return(-1.0);
	}
	
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("correlate_images_within_radii", __FILE__, __LINE__);
		return(-1);
	}
	
	img_to_float(p1);
	img_to_float(p2);
	
	if ( p1->std == 0 ) img_stats_within_radii(p1, minr, maxr);
	if ( p2->std == 0 ) img_stats_within_radii(p2, minr, maxr);
	
    float*  	    data1 = (float *) p1->data;
    float*  	    data2 = (float *) p2->data;
    
    unsigned long   i, x, y, z, n;
	float			dx, dy, dz, d2, rmax2 = maxr*maxr, rmin2 = minr*minr;
    double  		cc = 0, w = 0, norm = 1.0/(p1->std*p2->std);
    
	for ( n=0; n<p1->n; n++ ) {
		for ( z=0; z<p1->z; z++ ) {
			dz = z - p1->image[n].oz;
			dz *= dz;
			for ( y=0; y<p1->y; y++ ) {
				dy = y - p1->image[n].oy;
				dy *= dy;
				for ( x=0; x<p1->x; x++ ) {
					dx = x - p1->image[n].ox;
					dx *= dx;
					d2 = dx + dy + dz;
					i = ((n*p1->z + z)*p1->y + y)*p1->x + x;
					if ( d2 <= rmax2 && d2 >= rmin2 ) {
						data1[i] = (data1[i] - bas1)*(data2[i] - bas2)*norm;
						cc += data1[i];
						w += 1;
					} else {
						data1[i] = 0;
					}
				}
			}
		}
	}
	
	cc /= w;
	
	img_stats(p1);
    
    if ( verbose & VERB_RESULT )
		printf("Correlation:                    %g (%g)\n\n", cc, w);
	
	return(cc);
}

/************************************************************************
@Function: histomatch_images
@Description:
	Fits two images by matching the histogram of the second to the first.
@Algorithm:
	Both images are converted to floating point.
@Arguments:
	Bimage* p1			first image.
	Bimage* p2			second image.
	float bins			number of bins in the histograms.
@Returns:
	int 				0, <0 if error.
**************************************************************************/
int 		histomatch_images(Bimage* p1, Bimage* p2, int bins)
{
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("histomatch_images", __FILE__, __LINE__);
		return(-1);
	}
	
	img_to_float(p1);
	img_to_float(p2);
	
    float*  	    data1 = (float *) p1->data;
    float*  	    data2 = (float *) p2->data;
    unsigned long   i, datasize = (unsigned long) p1->n*p1->x*p1->y*p1->z;
	int 			h1, h2;
	float 			denom, fraction;
	int*			hist1 = (int *) balloc(bins*sizeof(int));
	int*			hist2 = (int *) balloc(bins*sizeof(int));
	int*			map = (int *) balloc(bins*sizeof(int));
	float*			dens = (float *) balloc(bins*sizeof(float));
    
    if ( verbose & VERB_LABEL )
	    printf("Matching histograms in %d bins:\n\n", bins);
	
    for ( i=0; i<datasize; i++ ) {
		h1 = (int) (bins*(data1[i] - p1->min)/(p1->max - p1->min));
		h2 = (int) (bins*(data2[i] - p2->min)/(p2->max - p2->min));
		if ( h1 >= bins ) h1 = bins - 1;
		if ( h2 >= bins ) h2 = bins - 1;
		hist1[h1]++;
		hist2[h2]++;
	}
	
	if ( verbose & VERB_DEBUG )
		printf("Histograms done\n");

	for ( i=bins-1; i>0; i-- ) {
		hist1[i-1] += hist1[i];
		hist2[i-1] += hist2[i];
	}
	
	if ( verbose & VERB_DEBUG )
		printf("Integration done\n");	
	
	for ( h2=0; h2<bins; h2++ ) {
		h1 = 0;
		while ( ( h1 < bins ) && ( hist1[h1] > hist2[h2] ) ) h1++;
		if ( h1 > 0 ) h1--;
		if ( h2 < bins - 1 )
			denom = hist1[h1+1] - hist1[h1] - hist2[h2+1] + hist2[h2];
		else
			denom = hist1[h1] - hist1[h1-1] - hist2[h2] + hist2[h2-1];
		if ( denom != 0 )
			fraction = (hist2[h2] - hist1[h1])/denom;
		else
			fraction = 0.5;
		map[h2] = h1;
		dens[h2] = (h2 + fraction)*(p1->max - p1->min)/bins + p1->min;
	}
	
	if ( verbose & VERB_DEBUG ) {
		printf("Mapping calculated\n");
		printf("Bin\tHist1\tHist2\tMapping\tDensity\n");
		for ( i=0; i<bins; i++ )
			printf("%4ld\t%7d\t%7d\t%4d\t%g\n", i, hist1[i], hist2[i], map[i], dens[i]);
	}	
	
    for ( i=0; i<datasize; i++ ) {
		h2 = (int) (bins*(data2[i] - p2->min)/(p2->max - p2->min));
		if ( h2 >= bins ) h2 = bins - 1;
		data1[i] = dens[map[h2]];
	}
	
	img_stats(p1);
	
	if ( verbose & VERB_DEBUG )
		printf("Histogram matching done\n");
	
	bfree(hist1, bins*sizeof(int));
	bfree(hist2, bins*sizeof(int));
	bfree(map, bins*sizeof(int));
	bfree(dens, bins*sizeof(float));
	
	return(0);
}

/************************************************************************
@Function: img_mask_with_image
@Description:
	Masks the first image with the second at a threshold level of the second.
@Algorithm:
	Both images are converted to floating point.
@Arguments:
	Bimage* p			image to mask.
	Bimage* m			masking image.
	float threshold 	threshold on second image.
	float fill			fill value for masking.
@Returns:
	int 				0, <0 if error.
**************************************************************************/
int 		img_mask_with_image(Bimage* p, Bimage* m, float threshold, float fill)
{
	if ( img_check_if_same_size(p, m) < 0 ) {
		error_show("img_mask_with_image", __FILE__, __LINE__);
		return(-1);
	}
	
	img_to_float(p);
	img_to_float(m);
	
    unsigned long   i, datasize = (unsigned long) p->n*p->x*p->y*p->z*p->c;
    float*			data = (float *) p->data;
    float*			mdata = (float *) m->data;
    
    if ( ( threshold < m->min ) || ( threshold > m->max ) ) 
    	threshold = (m->max + m->min)/2;
	
    if ( verbose & VERB_LABEL ) {
	    printf("Mask density at level:          %g\n",threshold);
    	printf("Fill value:                     %g\n\n",fill);
	}
	
    for ( i=0; i<datasize; i++ )
    	if ( mdata[i] < threshold ) data[i] = fill;

    return(0);
}

/************************************************************************
@Function: img_morph_blend
@Description:
	Blends the two images, creating a new set of sub-images.
@Algorithm:
	A number of images are created by blending the two input images in
	different ratios. The input images become the first and last sub-images of 
	the new image structure, with the intermediate images changing over from
	the first to the last:
			new_data = (1-fraction)*data1 + fraction*data2
	where fraction = index/(number - 1)
	At least 3 new sub-images are packed into the new image. 
	All of the header information in the first image is copied into the new image.
	Both images are converted to floating point.
@Arguments:
	Bimage* p1			first image.
	Bimage* p2			second image.
	int number			number of images in the series.
@Returns:
	Bimage* 			new image structure, NULL if error.
**************************************************************************/
Bimage* 	img_morph_blend(Bimage* p1, Bimage* p2, int number)
{
	if ( img_check_if_same_size(p1, p2) < 0 ) {
		error_show("img_morph_blend", __FILE__, __LINE__);
		return(NULL);
	}
	
	img_to_float(p1);
	img_to_float(p2);
	
    unsigned long   datasize = (unsigned long) p1->x*p1->y*p1->z*p1->c;
	
    if ( number < 3 ) number = 3;
    if ( sizeof(float)*datasize*number > 1e9 ) {
		number = (int) (1e9/(sizeof(float)*datasize));
		printf("The requested number of sub-images to create is too large!\n");
		printf("Changed to:                     %d\n", number);
	}
	
    if ( verbose & VERB_LABEL )
	    printf("Blending two images to create a series of %d sub-images\n", number);
	
    int     	i, j;
	float		fraction1, fraction2;
    float*  	data1 = (float *) p1->data;
    float*  	data2 = (float *) p2->data;
    
	// Copy the first image into the new image structure with a new number of sub-images
	Bimage* 	pnew = copy_img_header(p1, number);
	float*		newdata = (float *) balloc(number*datasize*sizeof(float));
	pnew->data = (char *) newdata;
	pnew->dataflag = 1;
	img_check_param(pnew);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_morph_blend: datasize = %ld\n", datasize);
	
    for ( i=0; i<number; i++ ) {
		fraction2 = i*1.0/(number - 1.0);
		fraction1 = 1 - fraction2;
		if ( verbose & VERB_PROCESS )
			printf("Calculating image %d: dnew = %g x d1 + %g x d2\n", 
					i, fraction1, fraction2);
		for ( j=0; j<datasize; j++ )
    		newdata[i*datasize+j] = fraction1*data1[j] + fraction2*data2[j];
		pnew->image[i].background = fraction1*p1->image[0].background + 
				fraction2*p2->image[0].background;
	}
	if ( verbose & VERB_PROCESS ) printf("\n");
	
	img_stats(pnew);

    return(pnew);
}
