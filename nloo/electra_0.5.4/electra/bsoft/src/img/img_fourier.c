/*
	img_fourier.c
	General FFT for n-dimensional data
		Implementing the FFTW library
	Author: Bernard Heymann
	Created: 19980805  	    Modified: 20050128
*/

#include <fftw.h>

#include "img_fourier.h"
#include "img_util.h"
#include "img_datatypes.h"
#include "matrix.h"
#include "matrix_linear.h"
#include "math_util.h"
#include "utilities.h"

#include <sys/time.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

/************************************************************************
@Function: img_fft
@Description:
	Fast Fourier transforms an image.
@Algorithm:
	FFTW library (www.fftw.org).
	A multi-image 1D, 2D and 3D data set is transformed forward or backward 
	and rescaled by 1/sqrt(N).  The forward transformation has a negative 
	signed exponent in the kernel and the backward transform a positive
	sign. The transformation is done in place and the resultant data are 
	returned within the original image structure.  For forward transforms, 
	the resultant data type is ComplexFloat, while for backward transforms, 
	it is Float.
@Arguments:
	fftw_direction dir	direction of transformation (either FFTW_FORWARD 
						or FFTW_BACKWARD)
	Bimage* p			image (replaced by transformed image).
@Returns:
	int 				error code.
**************************************************************************/
int 		img_fft(fftw_direction dir, Bimage* p)
{
	if ( img_fft_complex(dir, p, 1) )
		return(error_show("img_fft", __FILE__, __LINE__));
		
	if ( dir == FFTW_BACKWARD )
	    img_complex2real(p);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_fft: Transformation finished!!\n\n");
	
	return(0);
}

/************************************************************************
@Function: img_fft_complex
@Description:
	Fast Fourier transforms an image.
@Algorithm:
	FFTW library (www.fftw.org).
	A multi-image 1D, 2D and 3D data set is transformed forward or backward 
	and optionally rescaled by 1/sqrt(N) or N.  The forward transformation 
	has a negative signed exponent in the kernel and the backward transform 
	a positive sign. The transformation is done in place and the resultant 
	data are returned within the original image structure.
	For both directions the resultant image is complex.
@Arguments:
	fftw_direction dir	direction of transformation (either FFTW_FORWARD 
						or FFTW_BACKWARD)
	Bimage* p			image (replaced by transformed image).
	int norm_flag		normalization: 0=none, 1=sqrtN, 2=N.
@Returns:
	int 				error code.
**************************************************************************/
int 		img_fft_complex(fftw_direction dir, Bimage* p, int norm_flag)
{
	if ( !p->data ) {
		fprintf(stderr, "Error: Cannot Fourier transform - the data block is empty!\n");
		return(error_show("img_fft_complex", __FILE__, __LINE__));
	}
	
	if ( img_to_complex_float(p) )
		return(error_show("img_fft_complex", __FILE__, __LINE__));
	
    // Do FFT
    int				j, n[3], rank = 3;
    n[2] = (int) p->x;
    n[1] = (int) p->y;
    n[0] = (int) p->z;
    if ( p->z < 2 ) {
		rank = 2;
		n[1] = (int) p->x;
		n[0] = (int) p->y;
	}
    if ( p->y*p->z < 2 ) {
		rank = 1;
		n[0] = (int) p->x;
	}
    
    // Organize dimensions
    unsigned long   i, datasize = 1;
    for ( j = 0; j < rank; j++) datasize *= n[j];

    int     	    in_place = 1;
    int     	    flags = FFTW_ESTIMATE;
    if (in_place) flags |= FFTW_IN_PLACE;

    if (rank > MAX_RANK) return(-1);
	
	if ( sizeof(fftw_complex) != gettypesize(p->datatype) ) {
		printf("float size = %ld  fftw_real size = %ld\n", sizeof(float), sizeof(fftw_real));
		fprintf(stderr, "Error: img_fft_complex: The FFTW complex number size = %ld (should be %ld)!\n", 
			sizeof(fftw_complex), gettypesize(p->datatype));
		bexit(-1);
	}

	if ( verbose & VERB_PROCESS ) {
    	if ( dir == FFTW_FORWARD ) printf("Doing a forward FFT:\n");
    	else printf("Doing a backward FFT:\n");
//		if ( in_place ) printf("In place\n");
#ifdef FFTW_VERSION
		printf("Using FFTW:                     Version %s\n", FFTW_VERSION);
#endif
    	printf("Complex size:                   %ld\n", sizeof(fftw_complex));
    	printf("Rank and size:                  %d %ldx%ldx%ld\n", rank, p->x, p->y, p->z);
		printf("\n");
	} else if ( verbose & VERB_LABEL )
		printf("FFTW\n\n");
	
    // Specify FFTW arrays
    fftw_complex*		in = (fftw_complex *) p->data;
    fftw_complex*		out = NULL;
    if ( in_place ) out = in;
    else out = (fftw_complex *) fftw_malloc(datasize*sizeof(fftw_complex));

    if ( !in || !out ) {
	    printf("Error: Out of memory!!!\n");
	    exit(-1);
    }

    if ( rank < 2 ) {
    	fftw_plan     plan = fftw_create_plan( (int) p->x, dir, flags);
		for ( i=0; i<p->n; i++ ) {
			in = (fftw_complex *) (p->data+i*datasize*sizeof(complex_float));
	    	if (in_place) fftw(plan, 1, in, 1, 1, 0, 0, 0);
    		else fftw(plan, 1, in, 1, 1, out, 1, 1);
    	}
    	fftw_destroy_plan(plan);
    } else {
    	fftwnd_plan   plan_nd = fftwnd_create_plan(rank, n, dir, flags);
		for ( i=0; i<p->n; i++ ) {
			in = (fftw_complex *) (p->data+i*datasize*sizeof(complex_float));
    		if (in_place) fftwnd(plan_nd, 1, in, 1, 1, 0, 0, 0);
    		else fftwnd(plan_nd, 1, in, 1, 1, out, 1, 1);
    	}
    	fftwnd_destroy_plan(plan_nd);
    }
	
	// Scale data
	double		scale = 1.0/datasize;
	if ( norm_flag == 1 ) scale = sqrt(scale);
	if ( norm_flag ) {
		for ( i=0; i<datasize*p->n; i++ ) {
			out[i].re *= scale;
			out[i].im *= scale;
		}
		if ( verbose & VERB_DEBUG )
			printf("DEBUG img_fft_complex: scale = %g\n", scale);
	}
	
	p->data = (char *) out;
	
	if ( !in_place ) fftw_free(in);
    
    fftw_check_memory_leaks();
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_fft_complex: FFT done!\n\n");
	
	p->datatype = ComplexFloat;
	p->transform = Standard;
	if ( dir == FFTW_BACKWARD ) p->transform = NoTransform;
	
	return(0);
}

/************************************************************************
@Function: img_fft_times
@Description:
	Tests the execution times for a series of image sizes.
@Algorithm:
	FFTW library (www.fftw.org).
	Blank complex floating point images are created and transformed.
	Only the call to the complex FFT function is timed.
@Arguments:
	int ndim			number of dimensions (1,2,3).
	int minsize			minimum image size.
	int maxsize			maximum image size.
@Returns:
	int 				0.
**************************************************************************/
int 		img_fft_times(int ndim, int minsize, int maxsize)
{
	int				i, size, nprime;
	unsigned long*  prime = NULL;
    Bimage*			p = NULL;
	timeval			tvs, tvf;
	
	if ( ndim < 1 ) ndim = 1;
	if ( ndim > 3 ) ndim = 3;
	if ( minsize < 1 ) minsize = 1;
	if ( maxsize < minsize ) maxsize = minsize;
	
	printf("\nTiming the execution of fast Fourier transforms:\n");
	printf("Number of dimensions:           %d\n", ndim);
	printf("Size range:                     %d - %d\n\n", minsize, maxsize);
	printf("Size\tTime\tPrime factors\n");
	for ( size = minsize; size <= maxsize; size++ ) {
		switch ( ndim ) {
			case 1:
				p = init_img_with_parameters(ComplexFloat, 1, size, 1, 1, 1);
				break;
			case 2:
				p = init_img_with_parameters(ComplexFloat, 1, size, size, 1, 1);
				break;
			case 3:
				p = init_img_with_parameters(ComplexFloat, 1, size, size, size, 1);
				break;
		}
		prime = prime_factors(size, &nprime);
		gettimeofday(&tvs, NULL);
		if ( img_fft_complex(FFTW_FORWARD, p, 0) )
			return(error_show("img_fft_times", __FILE__, __LINE__));
		gettimeofday(&tvf, NULL);
		printf("%d\t%lg", size, (double) (tvf.tv_sec - tvs.tv_sec + 1e-6L*(tvf.tv_usec - tvs.tv_usec)));
		for ( i=0; i<nprime; i++ ) printf("\t%ld", prime[i]);
		printf("\n");
		bfree(prime, nprime*sizeof(unsigned long));
		kill_img(p);
	}
	printf("\n");

	return(0);
}

/************************************************************************
@Function: img_auto_correlate
@Description:
	Calculates an autocorrelation map by Fast Fourier transformation.
@Algorithm:
	FFTW library (www.fftw.org).
	A multi-image 1D, 2D and 3D data set is transformed forward, the 
	transform multiplied by its complex conjugate, followed by backward 
	transformation and rescaling by 1/(N*N). Data beyond the resolution 
	set in the image structure are zeroed. Therefore the correct setting 
	of units and resolution in the image are required. Defaults for the 
	units are usually 1 Angstrom/voxel and a zero resolution would
	include the whole image (i.e., no resolution limitation).
	The resultant image is Float.
@Arguments:
	Bimage* p			image (replaced by autocorrelation image).
@Returns:
	int 				0.
**************************************************************************/
int 		img_auto_correlate(Bimage* p)
{
	unsigned long   n, j, x, y, z;
	long			xx, yy, zz;
	float			sx2, sy2, sz2, res2;
	if ( p->resolution > 0 ) 
		res2 = 1/(p->resolution*p->resolution);
	else {
		res2 = 0;
		if ( p->x > 1 ) res2 += 1.0/(4.0*p->ux*p->ux);
		if ( p->y > 1 ) res2 += 1.0/(4.0*p->uy*p->uy);
		if ( p->z > 1 ) res2 += 1.0/(4.0*p->uz*p->uz);
	}

	if ( verbose & VERB_LABEL )
		printf("Auto-correlation at a resolution of %g A\n\n", 1/sqrt(res2));
	
	img_fft(FFTW_FORWARD, p);
	
	complex_float	temp;
	complex_float*	data = (complex_float *) p->data;
	
    for ( n=0; n<p->n; n++ ) {
		for ( z=0; z<p->z; z++ ) {
			zz = z;
			if ( z > (p->z - 1)/2 ) zz -= p->z;
			sz2 = zz/(p->z*p->uz);
			sz2 *= sz2;
			for ( y=0; y<p->y; y++ ) {
				yy = y;
				if ( y > (p->y - 1)/2 ) yy -= p->y;
				sy2 = yy/(p->y*p->uy);
				sy2 *= sy2;
				for ( x=0; x<p->x; x++ ) {
					xx = x;
					if ( x > (p->x - 1)/2 ) xx -= p->x;
					sx2 = xx/(p->x*p->ux);
					sx2 *= sx2;
					j = ((n*p->z + z)*p->y + y)*p->x + x;
					if ( sx2 + sy2 + sz2 <= res2 ) {
						temp = data[j];
    					data[j].re = temp.re*temp.re + temp.im*temp.im;
					} else {
						data[j].re = 0;
					}
    				data[j].im = 0;
				}
			}
		}
	}

	img_fft(FFTW_BACKWARD, p);
		
	return(0);
}

/************************************************************************
@Function: img_cross_correlate
@Description:
	Calculates an masked cross-correlation map by Fast Fourier transformation.
@Algorithm:
	FFTW library (www.fftw.org).
	Two equally sized multi-image 1D, 2D and 3D real space data sets are 
	packed into a complex data set and transformed forward. The transform
	is unpacked and masked with mask image before the first transform is
	multiplied with the complex conjugate of the second transform. This is 
	then back-transformed to obtain the cross-correlation map in real space.
	The mask must be composed of 0 and 1, and is converted to floating point.
	The mask can be omitted (NULL).
	The resultant cross-correlation image data type is floating point.
@Arguments:
	Bimage* p1			first image.
	Bimage* p2			second image.
	float hires			high resolution limit.
	float lores			low resolution limit.
	Bimage* pmask		binary mask (only 0 and 1), NULL if not desired.
@Returns:
	Bimage* 			cross-correlation image.
**************************************************************************/
Bimage* 	img_cross_correlate(Bimage* p1, Bimage* p2, float hires, 
				float lores, Bimage* pmask)
{
	Bimage* 	pc = img_pack_two_in_complex(p1, p2);
	if ( !pc ) return(NULL);
	
	unsigned long		i;
	for ( i=0; i<pc->n; i++ ) {
		pc->image[i].ox = p2->image[i].ox - p1->image[i].ox;
		pc->image[i].oy = p2->image[i].oy - p1->image[i].oy;
		pc->image[i].oz = p2->image[i].oz - p1->image[i].oz;
	}
	
	if ( hires < pc->ux ) hires = 2*pc->ux;
	if ( lores > pc->ux*pc->x ) lores = pc->ux*pc->x;
	
	if ( verbose & VERB_PROCESS ) {
		if ( pmask )
			printf("Cross-correlation with a mask:  %s\n", pmask->filename);
		else
			printf("Cross-correlation:\n");
		printf("Resolution range:               %g - %g A\n\n", hires, lores);
	}
	
	img_fft_complex(FFTW_FORWARD, pc, 0);
	
	img_combined_masked_complex_product(pc, hires, lores, pmask);
	
	img_fft_complex(FFTW_BACKWARD, pc, 0);
	
	img_complex2real(pc);
	
	return(pc);
}

/************************************************************************
@Function: img_find_shift
@Description:
	Calculates a cross-correlation map to find the shift for the pair of images.
@Algorithm:
	FFTW library (www.fftw.org).
	Two equally sized multi-image 1D, 2D and 3D data sets are transformed 
	forward, the first transform multiplied by the complex conjugate of
	the second transform, followed by backward transformation and 
	rescaling by 1/(N*N). Data beyond the resolution set in the first 
	image structure are not used. Therefore the correct setting 
	of units and resolution in the image are required. Defaults for the 
	units are usually 1 Angstrom/voxel and a zero resolution would
	include the whole image (i.e., no resolution limitation).
	A shift vector for each pair of images is calculated to
	determine the cross-correlation peak to sub-pixel resolution.
	Note: The first image is the reference and the shift returned is to
		transform the second to fit the first.
@Arguments:
	Bimage* pref		reference image.
	Bimage* p2			image to be shifted.
	Bimage* pmask		binary mask (only 0 and 1).
	float hires			high resolution limit.
	float lores			low resolution limit.
	float radius		search radius (if < 1, default 1e30).
	int refine_flag 	set to refine shift to subpixel resolution.
@Returns:
	Vector3* 			a set of three-value shift vectors, each for a
						pair of images.
**************************************************************************/
Vector3* 	img_find_shift(Bimage* pref, Bimage* p, Bimage* pmask, float hires, 
				float lores, float radius, int refine_flag)
{
	if ( verbose & VERB_PROCESS ) {
		if ( refine_flag )
			printf("Finding shift by cross-correlation and polynomial fitting\n\n");
		else
			printf("Finding shift by cross-correlation\n\n");
	}
	
	Bimage* 		pc = img_cross_correlate(pref, p, hires, lores, pmask);
	if ( !pc ) return(NULL);
	
	if ( verbose & VERB_DEBUG )
		write_img("cc.map", pc);
	
	VectorInt3		ivec;
	Vector3			origin = {0,0,0};
	Vector3*		shift = img_find_peak(pc, origin, radius);
	unsigned long	i, n;
	float*			cc = (float *) balloc(pc->n*sizeof(float));
	float*			data = (float *) pc->data;

	for ( n=0; n<pc->n; n++ ) {
		ivec = vectorint3_from_vector3(shift[n]);
		i = (ivec.z*pc->y + ivec.y)*pc->x + ivec.x;
		cc[n] = data[i];
		if ( shift[n].x > pc->x/2 ) shift[n].x -= pc->x;
		if ( shift[n].y > pc->y/2 ) shift[n].y -= pc->y;
		if ( shift[n].z > pc->z/2 ) shift[n].z -= pc->z;
	}
	
	if ( refine_flag ) img_refine_peak(pc, shift);
	
	if ( verbose & VERB_PROCESS ) {
		printf("Image\t   dx\t   dy\t   dz\t  CC\n");
		for ( n=0; n<pc->n; n++ )
			printf("%ld\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n", n+1, shift[n].x, shift[n].y, shift[n].z, cc[n]);
		printf("\n");
	}
	
	bfree(cc, pc->n*sizeof(float));
	
	kill_img(pc);
	
	return(shift);
}

/************************************************************************
@Function: img_find_peak
@Description:
	Finds the peak in an image to the nearest voxel.
@Algorithm:
	An image is searched for the global maximum (typically used to find
	the shift vector in a cross-correlation map).
	The peak vectors returned are in actual pixel coordinates (no wrapping).
@Arguments:
	Bimage* p			image (not altered).
	Vector3 origin		search origin.
	float radius		search radius (if < 1, default 1e30).
@Returns:
	Vector3* 			three-value peak vector for every sub-image.
**************************************************************************/
Vector3*	img_find_peak(Bimage* p, Vector3 origin, float radius)
{
	if ( !p->data ) return(NULL);
	
	if ( radius < 1 ) radius = 1e30;
	
	unsigned long   i, n, x, y, z;
	long			xh = p->x/2, yh = p->y/2, zh = p->z/2;
	float 			x2, y2, z2, r2 = radius*radius;
	double			value = 0, max;
	Vector3*		peak = (Vector3 *) balloc(p->n*sizeof(Vector3));
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* usdata = (unsigned short *) p->data;
    short* 	    	sdata = (short *) p->data;
    int* 	    	idata = (int *) p->data;
    float*  	    fdata = (float *) p->data;
    complex_short* 	csdata = (complex_short *) p->data;
    complex_int* 	cidata = (complex_int *) p->data;
    complex_float*  cfdata = (complex_float *) p->data;
    polar*  	    pdata = (polar *) p->data;
    
	if ( verbose & VERB_FULL ) {
	    printf("Finding the brightest voxel\n");
		printf("n\tx\ty\tz\n");
	}
	
	for ( n=0; n<p->n; n++ ) {
		max = -1e37;
		for ( z=0; z<p->z; z++ ) {
			z2 = origin.z - z;
			if ( z2 < -zh ) z2 += p->z;
			if ( z2 >  zh ) z2 -= p->z;
			z2 *= z2;
			for ( y=0; y<p->y; y++ ) {
				y2 = origin.y - y;
				if ( y2 < -yh ) y2 += p->y;
				if ( y2 >  yh ) y2 -= p->y;
				y2 *= y2;
				for ( x=0; x<p->x; x++ ) {
					x2 = origin.x - x;
					if ( x2 < -xh ) x2 += p->x;
					if ( x2 >  xh ) x2 -= p->x;
					x2 *= x2;
					if ( x2 + y2 + z2 <= r2 ) {
						i = ((n*p->z + z)*p->y + y)*p->x + x;
    					switch ( p->datatype ) {
    						case UChar: value = udata[i]; break;
    						case SChar: value = cdata[i]; break;
    						case UShort: value = usdata[i]; break;
    						case Short: value = sdata[i]; break;
    						case Int: value = idata[i]; break;
	    					case Float: value = fdata[i]; break;
    						case ComplexShort:
								value = csdata[i].re*csdata[i].re + csdata[i].im*csdata[i].im;
								break;
    						case ComplexInt:
								value = cidata[i].re*cidata[i].re + cidata[i].im*cidata[i].im;
								break;
	    					case ComplexFloat:
								value = cfdata[i].re*cfdata[i].re + cfdata[i].im*cfdata[i].im;
								break;
    						case Polar: value = pdata[i].amp; break;
    						default: break;
						}
						if ( max < value ) {
							max = value;
							peak[n].x = x;
							peak[n].y = y;
							peak[n].z = z;
//							printf("x,y,z = %d %d %d, max = %g\n", x, y, z, max/voxel_size);
						}
					}
				}
			}
		}
		if ( verbose & VERB_FULL )
			printf("%ld\t%.7g\t%.7g\t%.7g\n", n+1, peak[n].x, peak[n].y, peak[n].z);
    }
	if ( verbose & VERB_FULL ) printf("\n");
	
	return(peak);
}

/************************************************************************
@Function: img_refine_peak
@Description:
	Refines the position of a peak to sub-voxel resolution.
@Algorithm:
	The sub-voxel resolution peak in the vicinity of a voxel is defined 
	by fitting a 2D/3D second order function around the voxel.
	(typically used to find the shift vector in a cross-correlation map).
@Arguments:
	Bimage* pc			cross-correlation image (not altered).
	Vector3* shift		three-value shift vector for every sub-image.
@Returns:
	int 				0.
**************************************************************************/
int			img_refine_peak(Bimage* pc, Vector3* shift)
{
	if ( !pc->data ) return(-1);
	
	img_to_float(pc);
	
	unsigned long   i, j, n, kernel_size = 5, nterm = 7;
	long			x, y, z, xlo, xhi, ylo, yhi, zlo, zhi, ix, iy, iz;
	VectorInt3		int_shift;
	float*			data = (float *) pc->data;
	double*			a = (double *) balloc(nterm*nterm*sizeof(double));
	double*			b = (double *) balloc(nterm*sizeof(double));
	double*			vec = (double *) balloc(nterm*sizeof(double));

	long			xh = pc->x/2;
	long			yh = pc->y/2;
	long			zh = pc->z/2;
	
	if ( verbose & VERB_FULL ) {
	    printf("Finding the shift in each image\n");
		printf("n\tx\ty\tz\n");
	}
	
	for ( n=0; n<pc->n; n++ ) {
		int_shift.x = (int) shift[n].x;
		int_shift.y = (int) shift[n].y;
		int_shift.z = (int) shift[n].z;
		xlo = (long) (shift[n].x - kernel_size);
		xhi = (long) (shift[n].x + kernel_size);
		if ( labs(xlo) > pc->x ) xlo = -xh;
		if ( labs(xhi) > pc->x ) xhi =  xh;
		ylo = (long) (shift[n].y - kernel_size);
		yhi = (long) (shift[n].y + kernel_size);
		if ( labs(ylo) > pc->y ) ylo = -yh;
		if ( labs(yhi) > pc->y ) yhi =  yh;
		zlo = (long) (shift[n].z - kernel_size);
		zhi = (long) (shift[n].z + kernel_size);
		if ( labs(zlo) > pc->z ) zlo = -zh;
		if ( labs(zhi) > pc->z ) zhi =  zh;
		if ( verbose & VERB_DEBUG )
			printf("DEBUG img_refine_peak: %ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", 
					xlo, xhi, ylo, yhi, zlo, zhi);
		memset(a, 0, nterm*nterm*sizeof(double));
		memset(b, 0, nterm*sizeof(double));
		vec[0] = 1;
		for ( z=zlo; z<=zhi; z++ ) {
			iz = z;
			if ( iz < 0 ) iz += pc->z;
			if ( iz >= (long)pc->z ) iz -= pc->z;
			vec[3] = z;
			vec[6] = z*z;
			for ( y=ylo; y<=yhi; y++ ) {
				iy = y;
				if ( iy < 0 ) iy += pc->y;
				if ( iy >= (long)pc->y ) iy -= pc->y;
				vec[2] = y;
				vec[5] = y*y;
				for ( x=xlo; x<=xhi; x++ ) {
					ix = x;
					if ( ix < 0 ) ix += pc->x;
					if ( ix >= (long)pc->x ) ix -= pc->x;
					i = ((n*pc->z + iz)*pc->y + iy)*pc->x + ix;
//					printf("%ld\t%ld\t%ld\t%ld\n", ix, iy, iz, i);
					vec[1] = x;
					vec[4] = x*x;
					for ( j=0; j<nterm; j++ ) b[j] += vec[j]*data[i];
					for ( i=0; i<nterm; i++ )
						for ( j=0; j<=i; j++ ) a[nterm*i+j] += vec[i]*vec[j];
				}
			}
		}
		for ( i=0; i<nterm-1; i++ )
			for ( j=i+1; j<nterm; j++ ) a[nterm*i+j] = a[nterm*j+i];
		for ( i=0; i<nterm; i++ )
			if ( fabs(a[nterm*i+i]) < 1e-37 ) a[nterm*i+i] = 1;
		if ( verbose & VERB_DEBUG ) {
			printf("DEBUG img_refine_peak: Input matrix:\n");
			show_matrix((int) nterm, (int) nterm, a);
			printf("DEBUG img_refine_peak: Input vector:");
			for ( i=0; i<nterm; i++ ) printf("\t%g", b[i]);
			printf("\n");
		}
		matrix_LU_decomposition((int) nterm, a, b);
		if ( verbose & VERB_DEBUG )
			printf("DEBUG img_refine_peak: Coefficients: %.7g\t%.7g\t%.7g\t%.7g\t%.7g\t%.7g\t%.7g\n",
				b[0], b[1], b[2], b[3], b[4], b[5], b[6]);
		if ( b[4] ) shift[n].x = -0.5*b[1]/b[4];
		if ( b[5] ) shift[n].y = -0.5*b[2]/b[5];
		if ( b[6] ) shift[n].z = -0.5*b[3]/b[6];
		if ( fabs(int_shift.x - shift[n].x) > 2 ) shift[n].x = int_shift.x;
		if ( fabs(int_shift.y - shift[n].y) > 2 ) shift[n].y = int_shift.y;
		if ( fabs(int_shift.z - shift[n].z) > 2 ) shift[n].z = int_shift.z;
		while ( shift[n].x < -xh ) shift[n].x += pc->x;
		while ( shift[n].y < -yh ) shift[n].y += pc->y;
		while ( shift[n].z < -zh ) shift[n].z += pc->z;
		while ( shift[n].x >  xh ) shift[n].x -= pc->x;
		while ( shift[n].y >  yh ) shift[n].y -= pc->y;
		while ( shift[n].z >  zh ) shift[n].z -= pc->z;
		if ( verbose & VERB_FULL )
			printf("%ld\t%.7g\t%.7g\t%.7g\n", n+1, shift[n].x, shift[n].y, shift[n].z);
	}
	if ( verbose & VERB_FULL ) printf("\n");

	bfree(a, nterm*nterm*sizeof(double));
	bfree(b, nterm*sizeof(double));
	bfree(vec, nterm*sizeof(double));

	return(0);
}

/************************************************************************
@Function: img_change_transform_size
@Description:
	Resizes a "standard" transform.
@Algorithm:
	A standard transform is resized by inserting or removing zero
	rows or columns in the middle of the data set.
@Arguments:
	Bimage* p			complex or polar image.
	VectorInt3 size		new image size.
@Returns:
	int 				0.
**************************************************************************/
int 		img_change_transform_size(Bimage* p, VectorInt3 size)
{
    if ( size.x*size.y*size.z < 1 ) return(-1);
    if ( size.x == (long)p->x && size.y == (long)p->y && size.z == (long)p->z ) return(0);
	if ( p->datatype != ComplexFloat && p->datatype != Polar ) return(-1);
    
    unsigned long   i, j, n, x, y, z;
    long     	    xold, yold, zold, zerox, zeroy, zeroz;
	VectorInt3		hold, h, d;
    unsigned long   datasize = (unsigned long)p->n*size.x*size.y*size.z;
    complex_float*  data = (complex_float *) p->data;
    complex_float*  newdata = (complex_float *) balloc(datasize*sizeof(complex_float));
	float*			fom = p->fom;
	float*			newfom = NULL;
	if ( p->fomflag ) newfom = (float *) balloc(datasize*sizeof(float));
	
	hold.x = (int) ((p->x - 1)/2);
	hold.y = (int) ((p->y - 1)/2);
	hold.z = (int) ((p->z - 1)/2);
	h.x = (size.x - 1)/2;
	h.y = (size.y - 1)/2;
	h.z = (size.z - 1)/2;
	d.x = (int) (size.x - p->x);
	d.y = (int) (size.y - p->y);
	d.z = (int) (size.z - p->z);
	
	if ( verbose & VERB_LABEL )
	    printf("Changing to size:               %d %d %d voxels\n\n",
				size.x, size.y, size.z);

    for ( n=0; n<p->n; n++ ) {
	    for ( z=0; z<size.z; z++ ) {
    		zold = z;
			if ( z > h.z ) zold -= d.z;
			zeroz = 0;
			if ( z > hold.z && z <= d.z + hold.z ) zeroz += 1;
    		for ( y=0; y<size.y; y++ ) {
			 	yold = y;
				if ( y > h.y ) yold -= d.y;
				zeroy = zeroz;
				if ( y > hold.y && y <= d.y + hold.y ) zeroy += 1;
				for ( x=0; x<size.x; x++ ) {
					xold = x;
					if ( x > h.x ) xold -= d.x;
					zerox = zeroy;
					if ( x > hold.x && x <= d.x + hold.x ) zerox += 1;
					if ( zerox < 1 ) {
						i = ((n*p->z + zold)*p->y + yold)*p->x + xold;
						j = ((n*p->z + z)*size.y + y)*size.x + x;
						newdata[j].re = data[i].re;
						newdata[j].im = data[i].im;
						if ( p->fomflag ) newfom[j] = fom[i];
					}
				}
	    	}
		}
    }
    
    p->data = (char *) newdata;
	if ( p->fomflag ) p->fom = newfom;
    
    bfree(data, p->n*p->x*p->y*p->z*sizeof(complex_float));
    if ( fom ) bfree(fom, p->n*p->x*p->y*p->z*sizeof(float));
	
	p->ux = p->x*p->ux/size.x;
	p->uy = p->y*p->uy/size.y;
	p->uz = p->z*p->uz/size.z;
	p->x = p->px = size.x;
	p->y = p->py = size.y;
	p->z = p->pz = size.z;
	
	return(0);
}

/************************************************************************
@Function: img_phase_difference
@Description:
	Calculates the cosine of the phase difference between two images.
@Algorithm:
	Both images are Fourier transformed and the cosine of the phase
	difference calculated.
@Arguments:
	Bimage* p1			real space image.
	Bimage* p2			real space image.
	float resolution_hi upper resolution limit.
	float resolution_lo lower resolution limit.
@Returns:
	Bimage* 			phase difference image.
**************************************************************************/
Bimage* 	img_phase_difference(Bimage* p1, Bimage* p2, float resolution_hi, float resolution_lo)
{
	if ( resolution_hi <= 0 ) resolution_hi = 0.1;
	if ( resolution_lo < resolution_hi + 1 ) resolution_lo = resolution_hi + 1;
	
	Bimage* 		ppd = img_pack_two_in_complex(p1, p2);
	
	img_fft(FFTW_FORWARD, ppd);
	
	if ( verbose & VERB_PROCESS )
		printf("Cosine phase difference in the resolution range %g - %g\n", resolution_hi, resolution_lo);
	
	unsigned long   i, j, x, y, z, n;
	long			ix, iy, iz;
	unsigned long   datasize = (unsigned long)ppd->x*ppd->y*ppd->z*ppd->n;
	Vector3			shift;
	double			res_hi2 = 1.0/(resolution_hi*resolution_hi);
	double			res_lo2 = 1.0/(resolution_lo*resolution_lo);
	double			dy, dz, sx2, sy2, sz2, s2;
	Vector3			scale = {1.0/(ppd->x*ppd->ux),1.0/(ppd->y*ppd->uy),1.0/(ppd->z*ppd->uz)};
	complex_float	data1, data2;
	complex_float*	data = (complex_float *) ppd->data;
	
	float*			dphi = (float *) balloc(datasize*sizeof(float));
	
	for ( n=0; n<ppd->n; n++ ) {
		ppd->image[n].ox = p1->image[n].ox - p2->image[n].ox;
		ppd->image[n].oy = p1->image[n].oy - p2->image[n].oy;
		ppd->image[n].oz = p1->image[n].oz - p2->image[n].oz;
		shift.x = ppd->image[n].ox/ppd->x;
		shift.y = ppd->image[n].oy/ppd->y;
		shift.z = ppd->image[n].oz/ppd->z;
		if ( verbose & VERB_FULL )
			printf("%ld\t%7.4f\t%7.4f\t%7.4f\n", n+1, shift.x, shift.y, shift.z);
		for ( z=0; z<ppd->z; z++ ) { 
			iz = -z;
			if ( iz < 0 ) iz += ppd->z;
			dz = z;
			if ( z > (ppd->z - 1)/2 ) dz -= ppd->z;
			sz2 = dz*scale.z;
			sz2 *= sz2;
			for ( y=0; y<ppd->y; y++ ) { 
				iy = -y;
				if ( iy < 0 ) iy += ppd->y;
				dy = y;
				if ( y > (ppd->y - 1)/2 ) dy -= ppd->y;
				sy2 = dy*scale.y;
				sy2 *= sy2;
				for ( x=0; x<(ppd->x+1)/2; x++ ) { 
					ix = -x;
					if ( ix < 0 ) ix += ppd->x; 
					sx2 = x*scale.x;
					sx2 *= sx2;
					s2 = sx2 + sy2 + sz2;
					if ( s2 >= res_lo2 && s2 <= res_hi2 ) {
						i = ((n*ppd->z + z)*ppd->y + y)*ppd->x + x;
						j = ((n*ppd->z + iz)*ppd->y + iy)*ppd->x + ix;
						data1.re = 0.5*(data[i].re + data[j].re);
						data1.im = -0.5*(data[i].im - data[j].im);
						data2.re = 0.5*(data[i].im + data[j].im);
						data2.im = 0.5*(data[i].re - data[j].re);
						dphi[i] = atan2(data1.im,data1.re) - atan2(data2.im,data2.re);
						dphi[i] -= 2*PI*(x*shift.x + y*shift.y + z*shift.z);
						dphi[j] = dphi[i] = cos(dphi[i]);
					}
				}
			}
		}
	}
	
	ppd->data = (char *) dphi;
	ppd->datatype = Float;
	bfree(data, datasize*sizeof(complex_float));
	
	return(ppd);
}

/************************************************************************
@Function: img_average_phase_difference
@Description:
	Calculates the average of the cosine of the phase difference between 
	two images within given resolution shells.
@Algorithm:
	Both images are Fourier transformed and the cosine of the phase
	difference calculated and averaged.
@Arguments:
	Bimage* p1			real space image.
	Bimage* p2			real space image.
	float resolution_hi upper resolution limit.
	float resolution_lo lower resolution limit.
	int weighting		weighting type: 0, none; 1, p1-amp; 2, p2-amp; 3, both amps
@Returns:
	float	 			average of the cosine of the phase difference.
**************************************************************************/
float 		img_average_phase_difference(Bimage* p1, Bimage* p2, 
				float resolution_hi, float resolution_lo, int weighting)
{
	Bimage* 		ppd = img_pack_two_in_complex(p1, p2);
	
	img_fft(FFTW_FORWARD, ppd);
	
	if ( resolution_hi <= 0 ) resolution_hi = 0.1;
	if ( resolution_lo < resolution_hi + 1 ) resolution_lo = resolution_hi + 1;
	if ( weighting < 0 ) weighting = 0;
	if ( weighting > 3 ) weighting = 3;
	
	unsigned long   i, j, x, y, z, n;
	long			ix, iy, iz;
	double			res_hi2 = 1.0/(resolution_hi*resolution_hi);
	double			res_lo2 = 1.0/(resolution_lo*resolution_lo);
	double			dy, dz, sx2, sy2, sz2, szy2, szyx2, dphi, amp;
	Vector3			scale = {1.0/(ppd->x*ppd->ux),1.0/(ppd->y*ppd->uy),1.0/(ppd->z*ppd->uz)};
	Vector3			shift;
	complex_float	data1, data2;
	complex_float*	data = (complex_float *) ppd->data;
	
	double			sum_amp = 0;
	double			avg_pd = 0;
	
	for ( n=0; n<ppd->n; n++ ) {
		ppd->image[n].ox = p1->image[n].ox - p2->image[n].ox;
		ppd->image[n].oy = p1->image[n].oy - p2->image[n].oy;
		ppd->image[n].oz = p1->image[n].oz - p2->image[n].oz;
		shift.x = ppd->image[n].ox/ppd->x;
		shift.y = ppd->image[n].oy/ppd->y;
		shift.z = ppd->image[n].oz/ppd->z;
//		printf("Shift: %g %g %g\n", shift.x, shift.y, shift.z);
		for ( z=0; z<ppd->z; z++ ) {
			iz = -z;
			if ( iz < 0 ) iz += ppd->z;
			dz = z;
			if ( z > (ppd->z - 1)/2 ) dz -= ppd->z;
			sz2 = dz*scale.z;
			sz2 *= sz2;
			if ( sz2 <= res_hi2 ) for ( y=0; y<ppd->y; y++ ) { 
				iy = -y;
				if ( iy < 0 ) iy += ppd->y;
				dy = y;
				if ( y > (ppd->y - 1)/2 ) dy -= ppd->y;
				sy2 = dy*scale.y;
				sy2 *= sy2;
				szy2 = sz2 + sy2;
				if ( szy2 <= res_hi2 ) for ( x=0; x<(ppd->x+1)/2; x++ ) {
					ix = -x;
					if ( ix < 0 ) ix += ppd->x; 
					sx2 = x*scale.x;
					sx2 *= sx2;
					szyx2 = szy2 + sx2;
					if ( szyx2 <= res_hi2 && szyx2 >= res_lo2 ) {
						i = ((n*ppd->z + z)*ppd->y + y)*ppd->x + x;
						j = ((n*ppd->z + iz)*ppd->y + iy)*ppd->x + ix;
						data1.re = 0.5*(data[i].re + data[j].re);
						data1.im = -0.5*(data[i].im - data[j].im);
						data2.re = 0.5*(data[i].im + data[j].im);
						data2.im = 0.5*(data[i].re - data[j].re);
						switch ( weighting ) {
							case 1:
								amp = sqrt(data1.re*data1.re + data1.im*data1.im);
								break;
							case 2:
								amp = sqrt(data2.re*data2.re + data2.im*data2.im);
								break;
							case 3:
								amp = sqrt(data1.re*data1.re + data1.im*data1.im) +
										sqrt(data2.re*data2.re + data2.im*data2.im);
								break;
							default:						
								amp = 1;
						}
						dphi = atan2(data1.im,data1.re) - atan2(data2.im,data2.re);
						dphi -= 2*PI*(x*shift.x + y*shift.y + z*shift.z);
						while ( dphi <= -PI ) dphi += 2*PI; 
						while ( dphi >   PI ) dphi -= 2*PI;
						sum_amp += amp;
//						avg_pd += amp*dphi*dphi;
						avg_pd += amp*fabs(dphi);
					}
				}
			}
		}
	}
	
//	if ( sum_amp ) avg_pd = sqrt(avg_pd/sum_amp);
	if ( sum_amp ) avg_pd /= sum_amp;
	else avg_pd = PI;
	
	kill_img(ppd);
	
	return(avg_pd);
}

/************************************************************************
@Function: img_flip_phases
@Description:
	Flips the phases of an image based on a phase difference map. 
@Algorithm:
	.
@Arguments:
	Bimage* p			real space or reciprocal space image.
	Bimage* p2			reciprocal space phase difference map.
@Returns:
	int					0.
**************************************************************************/
int			img_flip_phases(Bimage* p, Bimage* pd)
{
	img_fft(FFTW_FORWARD, p);
	
	unsigned long   i;
	unsigned long   datasize = p->x*p->y*p->z*p->n;
	complex_float*	data = (complex_float *) p->data;
	float*			pdiff = (float *) pd->data;
	
	for ( i=0; i<datasize; i++ ) {
		if ( pdiff[i] < 0 ) {
			data[i].re = -data[i].re;
			data[i].im = -data[i].im;
		}
	}

	img_fft(FFTW_BACKWARD, p);
	
	return(0);
}

