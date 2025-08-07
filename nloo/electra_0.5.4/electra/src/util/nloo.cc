/* 
	nloo.cc
	Functions to analyze maps by NLOO method
	Author: Giovanni Cardone
	Created: 20060613	Modified:
*/ 

#include "nloo.h"


/************************************************************************
@Function: nloo_corr_estimate
@Description:
	Calculates correlation coefficient by Noise-compensated Leave-One_out (NLOO) method.
@Algorithm:
	The three images must be of the same dimensions.
	Two of them  are packed into a combined complex data set (the first
	image in the real part and the second image in the imaginary part).
	This is Fourier transformed and the individual transforms calculated
	based on the Hermitian or Friedel symmetry of transforms of real space
	data sets.

	Formula for NLOO2D 
	--------------------------------------------------------------------------------
	Cardone et al. (2005) J. Struct. Biology XXX, XXX-XXX 

	              sum(|Fmiss|*|Finput|*cos(phi_miss - phi_ref)) 
	              --------------------------------------------- 
	                sqrt( sum(|Fmiss|^2) * sum(|Finput|^2) ) 
	NLOO2D = -------------------------------------------------------- 
	              sum(|Ffull|*|Finput|*cos(phi_full - phi_ref)) 
	              --------------------------------------------- 
	                sqrt( sum(|Ffull|^2) * sum(|Finput|^2) ) 

	With respect to a NLOO resolution estimate, the sum is extendend to all
	frequency values within the given range.

@Arguments:
	Bimage* pref		reference image (input projeciton).
	Bimage* pmiss	missing image   (reprojection from partial reconstruction).
	Bimage* pfull	full image   (reprojection from entire reconstruction).
	float   lores	lower resolution limit
	float   hires	higher resolution limit
@Returns:
	float* 			correlation coefficient. Failure results in all 0.
**************************************************************************/
float* 	nloo_corr_estimate(Bimage* pref, Bimage* pmiss, Bimage* pfull, float lores, float hires)
{ 
	// Resolution estimates
    float*	corr_coeff = (float *) balloc(pref->n*sizeof(float));
	
	if ( pref->transform != NoTransform ) {
		fprintf(stderr, "Error(nloo_corr_estimate): File %s must be a real space map!\n", pref->filename); 
		return(corr_coeff); 
	}
	if ( pmiss->transform != NoTransform ) {
		fprintf(stderr, "Error(nloo_corr_estimate): File %s must be a real space map!\n", pmiss->filename); 
		return(corr_coeff); 
	}
	if ( pfull->transform != NoTransform ) {
		fprintf(stderr, "Error(nloo_corr_estimate): File %s must be a real space map!\n", pfull->filename); 
		return(corr_coeff); 
	}
	
	Bimage* pfmiss = copy_img(pmiss);
	img_fft(FFTW_FORWARD, pfmiss);

	Bimage* pffull = copy_img(pfull);
	img_fft(FFTW_FORWARD, pffull);
	
	Bimage* pfref = copy_img(pref);
	img_fft(FFTW_FORWARD, pfref);

	if ( pfmiss == NULL  || pfref == NULL || pffull == NULL) return(corr_coeff);

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG(nloo_corr_estimate) correlation: FFT done\n");
	}
	
	if ( lores <= 0 ) lores = 1e30;
	if ( hires <= 0 ) hires = 2*pref->ux;
	if ( lores < hires ) swap_floats(&lores, &hires);
	float	s2hi = 1/(hires*hires);
	float	s2lo = 1/(lores*lores);


//	float		deltar = bin_size/2.0;
	unsigned long i;
	unsigned int n;
	int 			x, y, z, xx, yy, zz, ix, iy, iz;
	float		Amiss, Afull, Aref, A2miss, A2full, A2ref, dphi_refmiss, dphi_reffull;
	float		freq_scale[3], rx, ry, rz;
	float		rx2, ry2, rz2, s2;
	
	// The frequency scaling is linked to the different dimensions of the
	// data set (important when the x, y and z dimensions are different)
	// The radius scaling is set to a value consistent with one pixel width
	// in reciprocal space for the longest dimension.
	freq_scale[0] = 1.0/(pref->x*pref->ux);  
	freq_scale[1] = 1.0/(pref->y*pref->uy);  
	if (pref->z >1)
		freq_scale[2] = 1.0/(pref->z*pref->uz);
	else
		freq_scale[2] = 0.;	
	
//	res_axis_dim(pref, &maxrad, &rad_scale, (float*) NULL);
	
	double		Frefmiss = 0.;
	double		Freffull = 0.; 
	double		F2full = 0.; 
	double		F2miss = 0.; 
//	float*		Frefmiss = (float *) balloc(maxrad*sizeof(float)); 
//	float*		Freffull = (float *) balloc(maxrad*sizeof(float)); 
//	float*		F2full = (float *) balloc(maxrad*sizeof(float)); 
//	float*		F2miss = (float *) balloc(maxrad*sizeof(float)); 
//	float*		NLOO = (float *) balloc(maxrad*sizeof(float)); 
//	float*		s = (float *) balloc(maxrad*sizeof(float)); 
//	float*		res = (float *) balloc(maxrad*sizeof(float)); 
	
	complex_float	refproj, missproj, fullproj;
	complex_float*	missdata = (complex_float *) pfmiss->data;
	complex_float*	fulldata = (complex_float *) pffull->data;
	complex_float*	refdata = (complex_float *) pfref->data;
	
	 
//	if ( verbose & VERB_DEBUG ) {
//		printf("DEBUG(nloo_corr_estimate) resolution: rad_scale = %g  maxrad = %d\n\n", 
//				rad_scale, maxrad);
//	}
/*	
	FILE*		fpt = NULL;
	if ( txtfile )
		if ( strlen(txtfile) > 0 )
			fpt = fopen(txtfile, "a");
*/
	for ( n=0; n<pref->n; n++ ) { 
		Frefmiss = 0.; 
		Freffull = 0.; 
		F2full = 0.; 
		F2miss = 0.; 
		
		for ( z=0; z < (int) pref->z; z++ ) { 
			zz = z; 
			if ( z > (int) (pref->z - 1)/2 ) zz -= pref->z;
			rz = zz*freq_scale[2];
			rz2 = rz*rz;
			iz = -z;
			if ( iz < 0 ) iz += pref->z;
			for ( y=0; y< (int) pref->y; y++ ) { 
				yy = y; 
				if ( y > (int) (pref->y - 1)/2 ) yy -= pref->y;
				ry = yy*freq_scale[1];
				ry2 = ry*ry;
				iy = -y;
				if ( iy < 0 ) iy += pref->y;
				for ( x=0; x< (int) pref->x; x++ ) { 
					i = ((n*pref->z + z)*pref->y + y)*pref->x + x;
					xx = x; 
					if ( x > (int) (pref->x - 1)/2 ) xx -= pref->x;
					rx = xx*freq_scale[0];
					rx2 = rx*rx;
					ix = -x;
					if ( ix < 0 ) ix += pref->x; 
					
					s2 = rx2 + ry2 + rz2;

					if ( s2 <= s2hi && s2 >= s2lo ) {
						missproj.re = missdata[i].re;
						missproj.im = missdata[i].im;
						fullproj.re = fulldata[i].re;
						fullproj.im = fulldata[i].im;
						refproj.re = refdata[i].re;
						refproj.im = refdata[i].im;
						A2miss = missproj.re*missproj.re + missproj.im*missproj.im; 
						A2full = fullproj.re*fullproj.re + fullproj.im*fullproj.im; 
						A2ref = refproj.re*refproj.re + refproj.im*refproj.im; 
						Amiss = sqrt(A2miss); 
						Afull = sqrt(A2full); 
						Aref = sqrt(A2ref); 
						dphi_refmiss = atan2(refproj.im,refproj.re) - atan2(missproj.im,missproj.re);
						dphi_reffull = atan2(refproj.im,refproj.re) - atan2(fullproj.im,fullproj.re);
						if ( !finite(dphi_refmiss) ) {
							fprintf(stderr, "%d %d %d: dphi not finite: ", x, y, z );
							fprintf(stderr, "reference = %g + %gi  missing projection = %g + %gi\n",
									refproj.re, refproj.im, missproj.re, missproj.im);
							dphi_refmiss = 0;
						}
						if ( !finite(dphi_reffull) ) {
							fprintf(stderr, "%d %d %d: dphi not finite: ", x, y, z );
							fprintf(stderr, "reference = %g + %gi  full projection = %g + %gi\n",
									refproj.re, refproj.im, fullproj.re, fullproj.im);
							dphi_reffull = 0;
						}
						while ( dphi_refmiss <= -M_PI ) dphi_refmiss += 2*M_PI; 
						while ( dphi_refmiss >   M_PI ) dphi_refmiss -= 2*M_PI;
						while ( dphi_reffull <= -M_PI ) dphi_reffull += 2*M_PI; 
						while ( dphi_reffull >   M_PI ) dphi_reffull -= 2*M_PI;
						Frefmiss += Aref*Amiss*cos(dphi_refmiss); 
						Freffull += Aref*Afull*cos(dphi_reffull); 
						F2miss += A2miss; 
						F2full += A2full; 
					}
				} 
			} 
		}
	 
/*		if ( verbose & VERB_PROCESS ) { 
			if ( pref->z == 1 )
				printf("\nNLOO2D Resolution measure for particle %d:\n", n); 
			else
				printf("\nNLOO2D Resolution measure:\n"); 
		} 
		if ( verbose & VERB_PROCESS ) { 
			if ( pref->z == 1 )
				printf("Radius\ts(1/A)\tRes(A)\tNLOO2D\n"); 
			else
				printf("Radius\ts(1/A)\tRes(A)\tNLOO2D\n"); 
		} 
*/	 
		corr_coeff[n] = 1.;
//		printf("%7.5f\t%7.5f\n",Frefmiss , Freffull);
		if (fabs(Frefmiss) < SMALLFLOAT || fabs(Freffull) < SMALLFLOAT)
			corr_coeff[n] = 0.;
		else
			corr_coeff[n] = (Frefmiss*sqrt(F2full))/(Freffull*sqrt(F2miss));
//		if ( corr_coeff[n] > 1. ) corr_coeff[n] = 1.;
//		if ( corr_coeff[n] < 0. ) corr_coeff[n] = 0.;
		if ( verbose & VERB_PROCESS ) {
				printf("%d\t%7.5f\n",n, corr_coeff[n]); 
		}

	}
	
	kill_img(pffull);
	kill_img(pfmiss);
	kill_img(pfref);

	return(corr_coeff); 
}
