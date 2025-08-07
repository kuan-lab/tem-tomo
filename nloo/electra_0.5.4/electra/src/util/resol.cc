/*
	resol.cc
	Functions to estimate resolution
	Author: Giovanni Cardone
	Created: 2004	Modified: 2004
*/

#include "resol.h"

/************************************************************************
@Function: res_fsc_estimate
@Description:
	Calculates resolution by Fourier Shell Correlation (FSC).
	Alternatively, it can evaluate the even/odd FSC (FSCeo).
@Algorithm:
	The two images (3d maps) must be of the same dimensions.
	They are packed into a combined complex data set (the first
	image in the real part and the second image in the imaginary part).
	This is Fourier transformed and the individual transforms calculated
	based on the Hermitian or Friedel symmetry of transforms of real space
	data sets.

	An optional input mask can be used to exclude parts of the transforms.
	The mask must consist of 0's and 1's, where 1 indicates inclusion.

	Formula for FSC
	-----------------------------------
	Saxton & Baumeister (1982) J. Microscopy 127, 127-138
	de la Fraga et al. (1995) Ultramicroscopy 60, 385-391
	           sum(|F1|*|F2|*cos(phi1 - phi2))
	FSC = ---------------------------------
	          sqrt( sum(|F1|^2) * sum(|F2|^2) )

	Formula for even/odd FSC
	-----------------------------------
	Cardone et al. (2005) J. Struct. Biology
	                   2 * FSC
	FSCeo = ---------------------------------
	                   FSC + 1

@Arguments:
	Bimage* p		image.
	Bimage* pmod	    (model) image.
	int bin_size	    annulus width of resolution region
	Bimage* pmask    reciprocal space mask (0 & 1, indicating inclusion of structure factors).
	char* txtfile	file name for text file output with all values.
	float* res_thr	resolution threshold array (many thresholds allowed)
	Wedge* mss_wdg   missing wedge
	int    fsceo 	evaluate FSCeo instead of FSC
	Bimage* vres     volumetric fourier correlation in the Fourier domain.
@Returns:
	float* 			resolution estimate for each image.
					Failure results in all zeroes.
**************************************************************************/
float* 	res_fsc_estimate(Bimage* p, Bimage* pmod, int bin_size, Bimage* pmask,
		  char* txtfile, float * res_thr, Wedge* mss_wdg, int fsceo, Bimage* vres)
{
	// Resolution estimates
	float*	res_est = (float *) balloc(p->n*RES_NTHRESHOLDS*sizeof(float));

	if ( p->transform != NoTransform ) {
		fprintf(stderr, "Error: File %s must be a real space map!\n", p->filename);
		return(res_est);
	}

	if ( pmod->transform != NoTransform ) {
		fprintf(stderr, "Error: File %s must be a real space map!\n", pmod->filename);
		return(res_est);
	}

	int nthr = res_thr_dim(res_thr);

	if ( pmask ) img_to_byte(pmask);
	if ( vres ) img_to_float(vres);

	// Pack the two images into one complex block
	Bimage* 	pc = img_pack_two_in_complex(p, pmod);
	if ( pc == NULL ) return(res_est);

	img_fft_complex(FFTW_FORWARD, pc,2);

	img_stats(pc);

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG fsc resolution: FFT done\n");
		printf("DEBUG fsc resolution: Average = %g, StDev = %g\n\n", pc->avg, pc->std);
	}

	float		deltar = 0.5*bin_size;
	unsigned long i;
	unsigned int	 n;
	int          x, y, z;
	long         j, nn, xx, yy, zz, ix, iy, iz;
	unsigned int	iradius, iradius2, maxrad, use_this;
	int 			ir, ir_low, ir_hi;
	float		radius, I1, I2, dphi;
	float		freq_scale[3], rx, ry, rz, rad_scale, ires, ires_1;
	float		invsigma2 = -0.5/(pc->std*pc->std);

	// The frequency scaling is linked to the different dimensions of the
	// data set (important when the x, y and z dimensions are different)
	// The radius scaling is set to a value consistent with one pixel width
	// in reciprocal space for the longest dimension.
	freq_scale[0] = 1.0/(pc->x*pc->ux);
	freq_scale[1] = 1.0/(pc->y*pc->uy);
	freq_scale[2] = 1.0/(pc->z*pc->uz);

//	maxrad is the main limit for the length of all the arrays,
//	and is derived from the set resolution
	res_axis_dim(pc, &maxrad, &rad_scale, (float*) NULL);

	if ( vres ) {
		vres->ux = freq_scale[0];
		vres->uy = freq_scale[1];
		vres->uz = freq_scale[2];
	}

	float*		num = (float *) balloc(maxrad*sizeof(float));
	float*		F1 = (float *) balloc(maxrad*sizeof(float));
	float*		F2 = (float *) balloc(maxrad*sizeof(float));
	float*		F1F2 = (float *) balloc(maxrad*sizeof(float));
	float*		s = (float *) balloc(maxrad*sizeof(float));
	float*		res = (float *) balloc(maxrad*sizeof(float));
	float*		FSC = (float *) balloc(maxrad*sizeof(float));

	complex_float	image, model;
	complex_float*	data = (complex_float *) pc->data;

	char*			mask = NULL;
	if ( pmask ) mask = (char *) pmask->data;

	float*			vdata = NULL;
	if ( vres ) vdata = (float *) vres->data;

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG fsc resolution: invsigma2 = %g\n\n", invsigma2 );
		printf("DEBUG fsc resolution: rad_scale = %g  maxrad = %d\n\n",
				rad_scale, maxrad);
	}

	FILE*		fpt = NULL;
	if ( txtfile )
		if ( strlen(txtfile) > 0 )
			fpt = fopen(txtfile, "w");

	char label[10];
	if (pc->z==1) {
		if (fsceo) strcpy(label, "FRCeo");
		else strcpy(label, "FRC");
	} else {
		if (fsceo) strcpy(label, "FSCeo");
		else strcpy(label, "FSC");
	}

	for ( n=0; n< pc->n; n++ ) {
		nn =n*RES_NTHRESHOLDS;
		if ( n>0 ) {
			memset(num, '\0', maxrad*sizeof(float));
			memset(F1, '\0', maxrad*sizeof(float));
			memset(F2, '\0', maxrad*sizeof(float));
			memset(F1F2, '\0', maxrad*sizeof(float));
			memset(FSC, '\0', maxrad*sizeof(float));
		}
		for ( z=0; z< (int) pc->z; z++ ) {
			zz = z;
			if ( z > (int) (pc->z - 1)/2 ) zz -= pc->z;
			rz = zz*freq_scale[2];
			iz = -z;
			if ( iz < 0 ) iz += pc->z;
			for ( y=0; y< (int) pc->y; y++ ) {
				yy = y;
				if ( y > (int) (pc->y - 1)/2 ) yy -= pc->y;
				ry = yy*freq_scale[1];
				iy = -y;
				if ( iy < 0 ) iy += pc->y;
				for ( x=0; x< (int) pc->x; x++ ) {
					use_this = 1;
					i = ((n*pc->z + z)*pc->y + y)*pc->x + x;
					if ( pmask ) use_this = mask[i];
					xx = x;
					if ( x > (int) (pc->x - 1)/2 ) xx -= pc->x;
					rx = xx*freq_scale[0];
					ix = -x;
					if ( ix < 0 ) ix += pc->x;
					if(mss_wdg)
						use_this = use_this && !wedge_point_inside(mss_wdg, rx, ry, rz);
					radius = rad_scale*sqrt(rx*rx + ry*ry + rz*rz);
					iradius = (int) radius;
					iradius2 = iradius + 1;
					if ( use_this && iradius2 < maxrad ) {
						ir_low = (int) (radius - deltar + 1.);
						if ( ir_low < 0 ) ir_low = 0;
						ir_hi = (int) (radius + deltar);
						if ( ir_hi >= (int) maxrad ) ir_hi = maxrad-1;
						// the following two commands ensure that resolution measure at zero component
						// is evaluated only with zero components of data, which furthermore do
						// not contribute to the other components
						if ( ir_low <= 0 && i>n*(int)pc->z) ir_low = 1;
						if ( i == n*pc->z ) ir_hi = 0;
						j = ((n*pc->z + iz)*pc->y + iy)*pc->x + ix;
						image.re = 0.5*(data[i].re + data[j].re);
						image.im = 0.5*(data[i].im - data[j].im);
						model.re = 0.5*(data[i].im + data[j].im);
						model.im = -0.5*(data[i].re - data[j].re);
						I1 = image.re*image.re + image.im*image.im;
						I2 = model.re*model.re + model.im*model.im;
						dphi = atan2(image.im,image.re) - atan2(model.im,model.re);
						if ( !finite(dphi) ) {
							fprintf(stderr, "%d %d %d: dphi not finite: ", x, y, z );
							fprintf(stderr, "image = %g + %gi  model = %g + %gi\n",
									image.re, image.im, model.re, model.im);
							dphi = 0;
						}
						while ( dphi <= -M_PI ) dphi += 2*M_PI;
						while ( dphi >   M_PI ) dphi -= 2*M_PI;
						if ( !finite(sqrt(I1)) )
							fprintf(stderr, "%d %d %d: I1 not finite\n", x, y, z);
						if ( !finite(sqrt(I2)) )
							fprintf(stderr, "%d %d %d: I2 not finite\n", x, y, z);
						for ( ir = ir_low; ir <= ir_hi; ir++) {
							num[ir] += 1;
							F1[ir] += I1;
							F2[ir] += I2;
							F1F2[ir] += sqrt(I1)*sqrt(I2)*cos(dphi);
						}
						if (vres ) {
							vdata[i] = 1.;
							if ( sqrt(I1)*sqrt(I2)*cos(dphi) < SMALLFLOAT )
								vdata[i] = 0.;
							else if ( I1*I2 > SMALLFLOAT )
								vdata[i] = sqrt(I1)*sqrt(I2)*cos(dphi)/sqrt(I1*I2);
							if ( vdata[i] > 1. ) vdata[i] = 1.;
							if (fsceo) vdata[i] = 2 * vdata[i] / (vdata[i] + 1);
						}
					} else if (vres) vdata[i]=0.;
				}
			}
		}

		if ( verbose & VERB_PROCESS ) {
			if ( pc->z == 1 )
				printf("\nResolution measures for particle %d:\n", n);
			else
				printf("\nResolution measures:\n");
		}
		if ( verbose & VERB_PROCESS )
			printf("Radius\ts(1/A)\tRes(A)\tSigma2\t%s\n",label);

		// frequency zero evaluated separately
		FSC[0] = 1.;
		s[0] = 0.;
		res[0] = 10.*rad_scale;
		if ( num[0] ) num[0] = 2.0/sqrt(2*num[0]);
		ires_1 = s[0];
		if ( verbose & VERB_PROCESS ) {
			printf("0\t%7.5f\t%9.3f\t%8.4f\t%8.4f\n",
				s[0], res[0], num[0], FSC[0]);
		}
		for ( i=1; i< maxrad; i++ ) {
			s[i] = i/rad_scale;
			ires = s[i];
			res[i] = 1./ires;
			if ( num[i] ) num[i] = 2.0/sqrt(2*num[i]);
			FSC[i] = 1.;
			if ( F1F2[i] < SMALLFLOAT ) FSC[i] = 0.;
			else if ( F1[i]*F2[i] > SMALLFLOAT ) FSC[i] = F1F2[i]/sqrt(F1[i]*F2[i]);
			if ( FSC[i] > 1. ) FSC[i] = 1.;
			if (fsceo) FSC[i] = 2 * FSC[i] / (FSC[i] + 1);
			for (int ii = 0; ii<nthr; ii++)
//				if ( res_est[nn+ii]==0.0 && FSC[i-1] >= res_thr[ii] && FSC[i] < res_thr[ii] )
				if ( FSC[i-1] >= res_thr[ii] && FSC[i] < res_thr[ii] )
					res_est[nn+ii] = ires + (ires_1 - ires)*(res_thr[ii] - FSC[i])/(FSC[i-1] - FSC[i]);
			if ( verbose & VERB_PROCESS ) {
				printf("%ld\t%7.5f\t%9.3f\t%8.4f\t%8.4f\n",
					i, s[i], res[i], num[i], FSC[i]);
			}
			ires_1 = ires;
		}

		if ( fpt ) {
			fprintf(fpt,"Radius\ts(1/A)\tRes(A)\t%s\n",label);
			for ( i=0; i<maxrad; i++ ) {
				fprintf(fpt,"%ld\t%7.5f\t%9.3f\t%8.4f\n",
						i, s[i], res[i], FSC[i]);
			}
		}


		for (int ii = 0; ii<nthr; ii++)
			if ( res_est[nn+ii] > SMALLFLOAT ) res_est[nn+ii] = 1./res_est[nn+ii];

		if ( verbose & VERB_RESULT ) {
			printf("Resolution estimates(A): # %d", n);
			for (int ii = 0; ii<nthr; ii++)
				printf("\t%s(%3.1f)=%8.4f", label, res_thr[ii], res_est[nn+ii]);
			printf("\n");
		}

	}

	if ( fpt ) fclose(fpt);

	bfree(num, maxrad*sizeof(float));
	bfree(F1, maxrad*sizeof(float));
	bfree(F2, maxrad*sizeof(float));
	bfree(F1F2, maxrad*sizeof(float));
	bfree(res, maxrad*sizeof(float));
	bfree(s, maxrad*sizeof(float));
	bfree(FSC, maxrad*sizeof(float));

	kill_img(pc);

	return(res_est);
}


/************************************************************************
@Function: res_nloo3d_init
@Description:
	Initialize integrated resolution arrays for nloo3d estimate.
@Algorithm:
@Arguments:
	Bimage* p		reference image for resolution dimensioning.
	float*  res_thr resolution thresholds
@Returns:
	Ires_nloo3d*		integrated resolution structure.
**************************************************************************/
Ires_nloo3d* res_nloo3d_init(Bimage* p, float* res_thr)
{

	float rad_scale;
	unsigned int maxrad;

	res_axis_dim(p, &maxrad, &rad_scale, (float*) NULL);

	Ires_nloo3d* int_res = (Ires_nloo3d *) balloc(sizeof(Ires_nloo3d));
	int_res->n = maxrad;
	int_res->np = 0;
	int_res->nthr = 0;
	int_res->rad_scale = rad_scale;
	int_res->nloo = (double *) balloc(maxrad*sizeof(double));
	int_res->nloo_missinput = (double *) balloc(maxrad*sizeof(double));
	int_res->nloo_miss = (double *) balloc(maxrad*sizeof(double));
	int_res->nloo_fullinput = (double *) balloc(maxrad*sizeof(double));
	int_res->nloo_full = (double *) balloc(maxrad*sizeof(double));

	int_res->nthr = res_thr_dim(res_thr);

	for (int i = 0; i<int_res->nthr; i++)
		 int_res->thresh[i]=res_thr[i];

	return (int_res);
}

/************************************************************************
@Function: res_nloo3d_estimate
@Description:
	Evaluate the resolution according to the 3D Noise-compensated Leave-One-Out
	(NLOO3D) method and writes the results to file.
@Algorithm:
	The measure is the result of the integration of the 2D measures over all
	the projections.
@Arguments:
	Ires_nloo3d* ir			integrated resolution structure.
	char* curve_file	curve	resolution function file name.
@Returns:
	int			0.
**************************************************************************/
int res_nloo3d_estimate(Ires_nloo3d* ir, char* curve_file)
{

	int i;
	int n = ir->n;

	FILE*		fc = NULL;
	FILE*		fr = NULL;

	float*		s = (float *) balloc(n*sizeof(float));
	float*		res = (float *) balloc(n*sizeof(float));

	if ( curve_file ) {
		if ( strlen(curve_file) > 0 ) {
			fc = fopen(curve_file, "w");
			fr = fopen(insert_in_filename(curve_file,"est",'_'), "w");
		}
	}
	if (fc) fprintf(fc,"Radius\ts(1/A)\tRes(A)\tNLOO3D\n");

	if (fr) {
		fprintf(fr,"Resolution estimates (A)\n");
		for ( int ii = 0; ii<ir->nthr; ii++)
			fprintf(fr,"NLOO3D(%3.1f)\t",ir->thresh[ii]);
		fprintf(fr,"\n");
	}

	//correct anomalous values
	if ( ir->np < 1 ) ir->np = 1;
	for ( i=0; i < n; i++ ) {
		if (fabs(ir->nloo_fullinput[i]) < SMALLFLOAT || fabs(ir->nloo_miss[i]) < SMALLFLOAT)
			ir->nloo[i] = 0.;
		else
			ir->nloo[i] = (ir->nloo_missinput[i]*sqrt(ir->nloo_full[i]))/
						(ir->nloo_fullinput[i]*sqrt(ir->nloo_miss[i]));
		if ( ir->nloo[i] < SMALLFLOAT ) ir->nloo[i] = 0.;
		if ( ir->nloo[i] > 1.0 ) ir->nloo[i] = 1.;
	}

	if ( verbose & VERB_PROCESS ) {
		printf("\nNLOO3D Resolution function:\n");
		printf("Radius\ts(1/A)\tRes(A)\tNLOO3D\n");
	}

	// initialize resolution axis arrays
	res[0] = 10.0*ir->rad_scale;
	for ( i=0; i<n; i++ ) {
		s[i] = i/ir->rad_scale;
		if(i>0) res[i] = ir->rad_scale/i;
		if ( verbose & VERB_PROCESS )
			printf("%d\t%7.5f\t%9.3f\t%8.4f\n", i, s[i], res[i], ir->nloo[i]);
	}

	// write curve to file
	if (fc) {
		for ( i=0; i<n; i++ )
			fprintf(fc,"%d\t%7.5f\t%9.3f\t%8.4f\n", i, s[i], res[i], ir->nloo[i]);
	}

	// compute resolution
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG(res_nloo3d_estimate): integrated resolution estimates running\n");
		printf("Radius\ts(1/A)\tRes(A)");
		for (int ii= 0; ii<ir->nthr; ii++)
			printf("\tNLOO3D(%3.1f)",ir->thresh[ii]);
		printf("\n");
	}

	for ( i=1; i<n; i++ ) {
		if ( verbose & VERB_DEBUG )
			printf("%d\t%7.5f\t%9.3f", i, s[i], res[i]);
		for (int ii = 0; ii<ir->nthr; ii++) {
			if ( ir->nloo[i-1] >= ir->thresh[ii] && ir->nloo[i] < ir->thresh[ii] )
				ir->estimate[ii] = s[i] + (s[i-1] - s[i]) *
					(ir->thresh[ii] - ir->nloo[i])/(ir->nloo[i-1] - ir->nloo[i]);
			if ( verbose & VERB_DEBUG )
				printf("\t%8.4f", 1./ir->estimate[ii]);
		}
		if ( verbose & VERB_DEBUG )
			printf("\n");
	}

	//convert the estimate into Angstroms
	for ( i = 0; i<RES_NTHRESHOLDS; i++)
		if ( ir->estimate[i] > 0.0 ) ir->estimate[i] = 1./ir->estimate[i];

	if (fr)	{
		for ( int ii = 0; ii<ir->nthr; ii++)
			fprintf(fr,"%8.2f\t",ir->estimate[ii]);
		fprintf(fr,"\n");
	}

	if ( verbose & VERB_RESULT ) {
		printf("Resolution estimates(A):");
		for (int ii = 0; ii<ir->nthr; ii++)
			printf("\tNLOO3D(%3.1f)=%8.4f", ir->thresh[ii], ir->estimate[ii]);
		printf("\n");
	}

	if ( fc ) fclose(fc);
	if ( fr ) fclose(fr);

	return (0);
}

/************************************************************************
@Function: res_nloo3d_kill
@Description:
	Free all memory allocated from Ires_nloo3d structure.
@Algorithm:
@Arguments:
	Ires_nloo3d* ir			integrated resolution structure.
@Returns:
	int			0.
**************************************************************************/
int	res_nloo3d_kill(Ires_nloo3d* ir)
{

	if ( ir == NULL) return(0);

	int n = ir->n;

	bfree(ir->nloo,n*sizeof(double));
	bfree(ir->nloo_missinput,n*sizeof(double));
	bfree(ir->nloo_miss,n*sizeof(double));
	bfree(ir->nloo_fullinput,n*sizeof(double));
	bfree(ir->nloo_full,n*sizeof(double));

	bfree(ir,sizeof(Ires_nloo3d));

	return (0);
}

/************************************************************************
@Function: res_nloo2d_estimate
@Description:
	Calculates resolution by 2D Noise-compensated Leave-One_out (NLOO2D).
@Algorithm:
	The three images (2d maps) must be of the same dimensions.
	Two of them  are packed into a combined complex data set (the first
	image in the real part and the second image in the imaginary part).
	This is Fourier transformed and the individual transforms calculated
	based on the Hermitian or Friedel symmetry of transforms of real space
	data sets.

	An optional input mask can be used to exclude parts of the transforms.
	The mask must consist of 0's and 1's, where 1 indicates inclusion.

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

@Arguments:
	int     pn		projection number
	float   tilt_ang tilt angle (degrees)
	Bimage* pref		reference image (input projeciton).
	Bimage* pmiss	missing image   (reprojection from partial reconstruction).
	Bimage* pfull	full image   (reprojection from entire reconstruction).
	int bin_size		annulus width of resolution region
	Bimage* pmask   	reciprocal space mask (0 & 1, indicating inclusion of structure factors).
	char* txtfile	file name for text file output with all values (written in append mode).
	float* res_thr	resolution threshold array (many thresholds allowed)
	Ires_nloo3d* int_res  integrated resolution (in order to evaluate NLOO3D)
@Returns:
	float* 			resolution estimate for each image.
					Failure results in all zeroes.
**************************************************************************/
float* 	res_nloo2d_estimate(int pn, float tilt_ang, Bimage* pref, Bimage* pmiss,
    Bimage* pfull, int bin_size, Bimage* pmask, char* txtfile, float * res_thr, Ires_nloo3d* int_res)
{
	// Resolution estimates
	float*	res_est = (float *) balloc(pref->n*RES_NTHRESHOLDS*sizeof(float));

	if ( pref->transform != NoTransform ) {
		fprintf(stderr, "Error(res_2dnloo_estimate): File %s must be a real space map!\n", pref->filename);
		return(res_est);
	}
	if ( pmiss->transform != NoTransform ) {
		fprintf(stderr, "Error(res_2dnloo_estimate): File %s must be a real space map!\n", pmiss->filename);
		return(res_est);
	}
	if ( pfull->transform != NoTransform ) {
		fprintf(stderr, "Error(res_2dnloo_estimate): File %s must be a real space map!\n", pfull->filename);
		return(res_est);
	}

	int nthr = res_thr_dim(res_thr);

	if ( pmask ) img_to_byte(pmask);

	Bimage* pfmiss = copy_img(pmiss);
	img_fft(FFTW_FORWARD, pfmiss);

	Bimage* pffull = copy_img(pfull);
	img_fft(FFTW_FORWARD, pffull);

	Bimage* pfref = copy_img(pref);
	img_fft(FFTW_FORWARD, pfref);

	if ( pfmiss == NULL  || pfref == NULL || pffull == NULL) return(res_est);

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG(res_2dnloo_estimate) resolution: FFT done\n");
	}

	float		deltar = bin_size/2.0;
	unsigned long i;
	unsigned int n;
	int 			nn, x, y, z, xx, yy, zz;
	unsigned int	iradius, iradius2, maxrad, use_this=1;
	int 			ir, ir_low, ir_hi;
	float		radius, Amiss, Afull, Aref, A2miss, A2full, A2ref, dphi_refmiss, dphi_reffull;
	float		freq_scale[3], rx, ry, rz, rad_scale, ires, ires0, ires_1;

	// The frequency scaling is linked to the different dimensions of the
	// data set (important when the x, y and z dimensions are different)
	// The radius scaling is set to a value consistent with one pixel width
	// in reciprocal space for the longest dimension.
	freq_scale[0] = 1.0/(pref->x*pref->ux);
	freq_scale[1] = 1.0/(pref->y*pref->uy);
	freq_scale[2] = 0.;

	res_axis_dim(pref, &maxrad, &rad_scale, (float*) NULL);

	float*		Frefmiss = (float *) balloc(maxrad*sizeof(float));
	float*		Freffull = (float *) balloc(maxrad*sizeof(float));
	float*		F2full = (float *) balloc(maxrad*sizeof(float));
	float*		F2miss = (float *) balloc(maxrad*sizeof(float));
	float*		NLOO = (float *) balloc(maxrad*sizeof(float));
	float*		s = (float *) balloc(maxrad*sizeof(float));
	float*		res = (float *) balloc(maxrad*sizeof(float));

	complex_float	refproj, missproj, fullproj;
	complex_float*	missdata = (complex_float *) pfmiss->data;
	complex_float*	fulldata = (complex_float *) pffull->data;
	complex_float*	refdata = (complex_float *) pfref->data;

	char*			mask = NULL;
	if ( pmask ) mask = (char *) pmask->data;

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG(res_2dnloo_estimate) resolution: rad_scale = %g  maxrad = %d\n\n",
				rad_scale, maxrad);
	}

	FILE*		fpt = NULL;
	if ( txtfile )
		if ( strlen(txtfile) > 0 )
			fpt = fopen(txtfile, "a");

	for ( n=0; n<pref->n; n++ ) {
		nn =n*RES_NTHRESHOLDS;
		if ( n>0 ) {
			memset(Frefmiss, '\0', maxrad*sizeof(float));
			memset(Freffull, '\0', maxrad*sizeof(float));
			memset(F2full, '\0', maxrad*sizeof(float));
			memset(F2miss, '\0', maxrad*sizeof(float));
			memset(NLOO, '\0', maxrad*sizeof(float));
		}

		for ( z=0; z < (int) pref->z; z++ ) {
			zz = z;
			if ( z > (int) (pref->z - 1)/2 ) zz -= pref->z;
			rz = zz*freq_scale[2];
			for ( y=0; y< (int) pref->y; y++ ) {
				yy = y;
				if ( y > (int) (pref->y - 1)/2 ) yy -= pref->y;
				ry = yy*freq_scale[1];
				for ( x=0; x< (int) pref->x; x++ ) {
					i = ((n*pref->z + z)*pref->y + y)*pref->x + x;
					if ( pmask ) use_this = mask[i];
					xx = x;
					if ( x > (int) (pref->x - 1)/2 ) xx -= pref->x;
					rx = xx*freq_scale[0];
					radius = rad_scale*sqrt(rx*rx + ry*ry + rz*rz);
					iradius = (int) radius;
					iradius2 = iradius + 1;
					if ( use_this && iradius2 < maxrad ) {
						ir_low = (int) (radius - deltar + 1.);
						if ( ir_low < 0 ) ir_low = 0;
						ir_hi = (int) (radius + deltar);
						if ( ir_hi >= (int) maxrad ) ir_hi = maxrad-1;
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
						for ( ir = ir_low; ir <= ir_hi; ir++) {
							Frefmiss[ir] += Aref*Amiss*cos(dphi_refmiss);
							Freffull[ir] += Aref*Afull*cos(dphi_reffull);
							F2miss[ir] += A2miss;
							F2full[ir] += A2full;
						}
					}
				}
			}
		}

		if ( verbose & VERB_PROCESS ) {
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

		if ( int_res ) int_res->np++;

		if (fpt && int_res) fprintf(fpt,"Projection # %d - Tilt angle: %5.2f\n",pn, tilt_ang);

		ires_1 = ires0 = 0.1/rad_scale;
		res[0] = 1./ires0;
		for ( i=0; i<maxrad; i++ ) {
			if ( int_res ) {
			// update integrated resolution
				int_res->nloo_missinput[i] += Frefmiss[i];
				int_res->nloo_miss[i] += F2miss[i];
				int_res->nloo_fullinput[i] += Freffull[i];
				int_res->nloo_full[i] += F2full[i];
			}
			s[i] = i/rad_scale;
			ires = s[i];
			if ( i > 0 ) res[i] = 1./ires;
			NLOO[i] = 1.;
			if (fabs(Frefmiss[i]) < SMALLFLOAT || fabs(Freffull[i]) < SMALLFLOAT)
				NLOO[i] = 0.;
			else
				NLOO[i] = (Frefmiss[i]*sqrt(F2full[i]))/(Freffull[i]*sqrt(F2miss[i]));
			if ( NLOO[i] > 1. ) NLOO[i] = 1.;
			if ( NLOO[i] < 0. ) NLOO[i] = 0.;
			if ( i > 0 ) {
				for (int ii = 0; ii<nthr; ii++)
	//				if ( res_est[nn+ii]==0.0 && NLOO[i-1] >= res_thr[ii] && NLOO[i] < res_thr[ii] )
					if ( NLOO[i-1] >= res_thr[ii] && NLOO[i] < res_thr[ii] )
						res_est[nn+ii] = ires + (ires_1 - ires)*(res_thr[ii] - NLOO[i])/(NLOO[i-1] - NLOO[i]);
			}
			if ( verbose & VERB_PROCESS ) {
					printf("%ld\t%7.5f\t%9.3f\t%8.4f\n",i, s[i], res[i], NLOO[i]);
			}
			ires_1 = ires;
		}

		if ( fpt ) {
			fprintf(fpt,"Radius\ts(1/A)\tRes(A)\tNLOO2D\n");
			for ( i=0; i<maxrad; i++ ) {
				fprintf(fpt,"%ld\t%7.5f\t%9.3f\t%8.4f\n",
						i, s[i], res[i], NLOO[i]);
			}
		}

		for (int ii = 0; ii<nthr; ii++)
			if ( res_est[nn+ii] > SMALLFLOAT ) res_est[nn+ii] = 1./res_est[nn+ii];

		if ( verbose & VERB_RESULT ) {
			printf("Resolution estimates(A): # %d", n);
			for (int ii = 0; ii<nthr; ii++)
				printf("\tNLOO2D(%3.1f)=%8.4f", res_thr[ii], res_est[nn+ii]);
			printf("\n");
		}
	}


	if ( fpt ) fclose(fpt);

	bfree(Freffull, maxrad*sizeof(float));
	bfree(Frefmiss, maxrad*sizeof(float));
	bfree(F2miss, maxrad*sizeof(float));
	bfree(F2full, maxrad*sizeof(float));
	bfree(NLOO, maxrad*sizeof(float));
	bfree(s, maxrad*sizeof(float));
	bfree(res, maxrad*sizeof(float));

	kill_img(pffull);
	kill_img(pfmiss);
	kill_img(pfref);

	return(res_est);
}

/************************************************************************
@Function: res_rd_nloo2d_estimate
@Description:
	Calculates resolution by 2D Noise-compensated Leave-One_out (NLOO2D).
@Algorithm:
	The three images (2d maps) must be of the same dimensions.
	Two of them  are packed into a combined complex data set (the first
	image in the real part and the second image in the imaginary part).
	This is Fourier transformed and the individual transforms calculated
	based on the Hermitian or Friedel symmetry of transforms of real space
	data sets.

	An optional input mask can be used to exclude parts of the transforms.
	The mask must consist of 0's and 1's, where 1 indicates inclusion.

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

@Arguments:
	Bimage* pref		reference image (input projeciton).
	Bimage* pmiss	missing image   (reprojection from partial reconstruction).
	Bimage* pfull	full image   (reprojection from entire reconstruction).
	int bin_size		annulus width of resolution region
	Bimage* pmask   	reciprocal space mask (0 & 1, indicating inclusion of structure factors).
	char* txtfile	file name for text file output with all values (written in append mode).
	float* res_thr	resolution threshold array (many thresholds allowed)
	Ires_nloo3d* int_res  integrated resolution (in order to evaluate NLOO3D)
	Vnloo*	vnl		volumetric resolution
@Returns:
	float* 			resolution estimate for each image.
					Failure results in all zeroes.
**************************************************************************/
float* 	res_rd_nloo2d_estimate(Bimage* pref, Bimage* pmiss, Bimage* pfull, int bin_size, Bimage* pmask, char* txtfile, float * res_thr, Ires_nloo3d* int_res, Vnloo* vnl)
{
	// Resolution estimates
	float*	res_est = (float *) balloc(pref->n*RES_NTHRESHOLDS*sizeof(float));

	if ( pref->transform != NoTransform ) {
		fprintf(stderr, "Error(res_2dnloo_estimate): File %s must be a real space map!\n", pref->filename);
		return(res_est);
	}
	if ( pmiss->transform != NoTransform ) {
		fprintf(stderr, "Error(res_2dnloo_estimate): File %s must be a real space map!\n", pmiss->filename);
		return(res_est);
	}
	if ( pfull->transform != NoTransform ) {
		fprintf(stderr, "Error(res_2dnloo_estimate): File %s must be a real space map!\n", pfull->filename);
		return(res_est);
	}

	int nthr = res_thr_dim(res_thr);

	if ( pmask ) img_to_byte(pmask);

	Bimage* pfmiss = copy_img(pmiss);
	img_fft(FFTW_FORWARD, pfmiss);

	Bimage* pffull = copy_img(pfull);
	img_fft(FFTW_FORWARD, pffull);

	Bimage* pfref = copy_img(pref);
	img_fft(FFTW_FORWARD, pfref);

	if ( pfmiss == NULL  || pfref == NULL || pffull == NULL) return(res_est);

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG(res_2dnloo_estimate) resolution: FFT done\n");
	}

	float		deltar = bin_size/2.0;
	unsigned long i;
	unsigned int n;
	int 			nn, x, y, z, xx, yy, zz, ix, iy, iz;
	unsigned int	iradius, iradius2, maxrad, use_this=1;
	int 			ir, ir_low, ir_hi;
	float		radius, Amiss, Afull, Aref, A2miss, A2full, A2ref, dphi_refmiss, dphi_reffull;
	float		Armdp, Arfdp;
	float		freq_scale[3], rx, ry, rz, rad_scale, ires, ires0, ires_1;

	// The frequency scaling is linked to the different dimensions of the
	// data set (important when the x, y and z dimensions are different)
	// The radius scaling is set to a value consistent with one pixel width
	// in reciprocal space for the longest dimension.
	freq_scale[0] = 1.0/(pref->x*pref->ux);
	freq_scale[1] = 1.0/(pref->y*pref->uy);
	freq_scale[2] = 0.;

	res_axis_dim(pref, &maxrad, &rad_scale, (float*) NULL);

	float*		Frefmiss = (float *) balloc(maxrad*sizeof(float));
	float*		Freffull = (float *) balloc(maxrad*sizeof(float));
	float*		F2full = (float *) balloc(maxrad*sizeof(float));
	float*		F2miss = (float *) balloc(maxrad*sizeof(float));
	float*		NLOO = (float *) balloc(maxrad*sizeof(float));
	float*		s = (float *) balloc(maxrad*sizeof(float));
	float*		res = (float *) balloc(maxrad*sizeof(float));

	complex_float	refproj, missproj, fullproj;
	complex_float*	missdata = (complex_float *) pfmiss->data;
	complex_float*	fulldata = (complex_float *) pffull->data;
	complex_float*	refdata = (complex_float *) pfref->data;

	char*			mask = NULL;
	if ( pmask ) mask = (char *) pmask->data;

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG(res_2dnloo_estimate) resolution: rad_scale = %g  maxrad = %d\n\n",
				rad_scale, maxrad);
	}

	FILE*		fpt = NULL;
	if ( txtfile )
		if ( strlen(txtfile) > 0 )
			fpt = fopen(txtfile, "a");

	for ( n=0; n<pref->n; n++ ) {
		nn =n*RES_NTHRESHOLDS;
		if ( n>0 ) {
			memset(Frefmiss, '\0', maxrad*sizeof(float));
			memset(Freffull, '\0', maxrad*sizeof(float));
			memset(F2full, '\0', maxrad*sizeof(float));
			memset(F2miss, '\0', maxrad*sizeof(float));
			memset(NLOO, '\0', maxrad*sizeof(float));
		}

		for ( z=0; z < (int) pref->z; z++ ) {
			zz = z;
			if ( z > (int) (pref->z - 1)/2 ) zz -= pref->z;
			rz = zz*freq_scale[2];
			iz = -z;
			if ( iz < 0 ) iz += pref->z;
			if (vnl) vnloo_reset_tmparr(vnl);
			for ( y=0; y< (int) pref->y; y++ ) {
				yy = y;
				if ( y > (int) (pref->y - 1)/2 ) yy -= pref->y;
				ry = yy*freq_scale[1];
				iy = -y;
				if ( iy < 0 ) iy += pref->y;
				for ( x=0; x< (int) pref->x; x++ ) {
					i = ((n*pref->z + z)*pref->y + y)*pref->x + x;
					if ( pmask ) use_this = mask[i];
					xx = x;
					if ( x > (int) (pref->x - 1)/2 ) xx -= pref->x;
					rx = xx*freq_scale[0];
					ix = -x;
					if ( ix < 0 ) ix += pref->x;
					radius = rad_scale*sqrt(rx*rx + ry*ry + rz*rz);
					iradius = (int) radius;
					iradius2 = iradius + 1;
					if ( use_this && iradius2 < maxrad ) {
						ir_low = (int) (radius - deltar + 1.);
						if ( ir_low < 0 ) ir_low = 0;
						ir_hi = (int) (radius + deltar);
						if ( ir_hi >= (int) maxrad ) ir_hi = maxrad-1;
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
						Armdp = Aref*Amiss*cos(dphi_refmiss);
						Arfdp = Aref*Afull*cos(dphi_reffull);
						if (vnl) vnloo_put(vnl, x, y, Armdp, Arfdp, A2miss, A2full);
						for ( ir = ir_low; ir <= ir_hi; ir++) {
							Frefmiss[ir] += Armdp;
							Freffull[ir] += Arfdp;
							F2miss[ir] += A2miss;
							F2full[ir] += A2full;
						}
					}
				}
			}
			if (vnl) vnloo_add_slice(vnl, pref);
		}

		if ( verbose & VERB_PROCESS ) {
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

		if ( int_res ) int_res->np++;

		if (fpt && int_res) fprintf(fpt,"Projection # %d\n",int_res->np);

		ires_1 = ires0 = 0.1/rad_scale;
		res[0] = 1./ires0;
		for ( i=0; i<maxrad; i++ ) {
			if ( int_res ) {
			// update integrated resolution
				int_res->nloo_missinput[i] += Frefmiss[i];
				int_res->nloo_miss[i] += F2miss[i];
				int_res->nloo_fullinput[i] += Freffull[i];
				int_res->nloo_full[i] += F2full[i];
			}
			s[i] = i/rad_scale;
			ires = s[i];
			if ( i > 0 ) res[i] = 1./ires;
			NLOO[i] = 1.;
			if (fabs(Frefmiss[i]) < SMALLFLOAT || fabs(Freffull[i]) < SMALLFLOAT)
				NLOO[i] = 0.;
			else
				NLOO[i] = (Frefmiss[i]*sqrt(F2full[i]))/(Freffull[i]*sqrt(F2miss[i]));
			if ( NLOO[i] > 1. ) NLOO[i] = 1.;
			if ( NLOO[i] < 0. ) NLOO[i] = 0.;
			if ( i > 0 ) {
				for (int ii = 0; ii<nthr; ii++)
	//				if ( res_est[nn+ii]==0.0 && NLOO[i-1] >= res_thr[ii] && NLOO[i] < res_thr[ii] )
					if ( NLOO[i-1] >= res_thr[ii] && NLOO[i] < res_thr[ii] )
						res_est[nn+ii] = ires + (ires_1 - ires)*(res_thr[ii] - NLOO[i])/(NLOO[i-1] - NLOO[i]);
			}
			if ( verbose & VERB_PROCESS ) {
					printf("%ld\t%7.5f\t%9.3f\t%8.4f\n",i, s[i], res[i], NLOO[i]);
			}
			ires_1 = ires;
		}

		if ( fpt ) {
			fprintf(fpt,"Radius\ts(1/A)\tRes(A)\tNLOO2D\n");
			for ( i=0; i<maxrad; i++ ) {
				fprintf(fpt,"%ld\t%7.5f\t%9.3f\t%8.4f\n",
						i, s[i], res[i], NLOO[i]);
			}
		}

		for (int ii = 0; ii<nthr; ii++)
			if ( res_est[nn+ii] > SMALLFLOAT ) res_est[nn+ii] = 1./res_est[nn+ii];

		if ( verbose & VERB_RESULT ) {
			printf("Resolution estimates(A): # %d", n);
			for (int ii = 0; ii<nthr; ii++)
				printf("\tNLOO2D(%3.1f)=%8.4f", res_thr[ii], res_est[nn+ii]);
			printf("\n");
		}
	}


	if ( fpt ) fclose(fpt);

	bfree(Freffull, maxrad*sizeof(float));
	bfree(Frefmiss, maxrad*sizeof(float));
	bfree(F2miss, maxrad*sizeof(float));
	bfree(F2full, maxrad*sizeof(float));
	bfree(NLOO, maxrad*sizeof(float));
	bfree(s, maxrad*sizeof(float));
	bfree(res, maxrad*sizeof(float));

	kill_img(pffull);
	kill_img(pfmiss);
	kill_img(pfref);

	return(res_est);
}

/************************************************************************
@Function: res_axis_dim
@Description:
	Dimension the resolution axis in the reciprocal space.
@Algorithm:
	Evaluate the number of point in the axis, which is determined
	from the size and the sampling of the given image, and from
	the resolution limit 	defined in the image parameters.
	It also evaluate and returns the scaling and the max resolution
@Arguments:
	Bimage* p			reference image for resolution dimensioning.
	int*	    rad_size		resolution axis size
	float*	rad_sampl	radial frequency sampling
	float*	max_res		maximum resolution values
@Returns:
	int				<0 if error.
**************************************************************************/
int res_axis_dim(Bimage* p, unsigned int* rad_size, float* rad_sampl, float* max_res)
{

	float rad_scale, resolution;
	unsigned int maxrad;

	rad_scale = p->x*p->ux;
	resolution = p->resolution;
	if ( rad_scale < p->y*p->uy ) rad_scale = p->y*p->uy;
	if ( rad_scale < p->z*p->uz ) rad_scale = p->z*p->uz;
	if ( p->x > 1 && resolution < 2*p->ux ) resolution = 2*p->ux;
	if ( p->y > 1 && resolution < 2*p->uy ) resolution = 2*p->uy;
	if ( p->z > 1 && resolution < 2*p->uz ) resolution = 2*p->uz;
	if ( resolution > 0.5*p->x*p->ux ) resolution = 0.5*p->x*p->ux;

//	maxrad is the main limit for the length of all the arrays,
//	and is derived from the set resolution
	maxrad = (int) (2 + rad_scale / resolution);

	if (rad_size) *rad_size = maxrad;
	if (rad_sampl) *rad_sampl = rad_scale;
	if (max_res) *max_res = resolution;

	return(0);

}

/************************************************************************
@Function: res_thr_dim
@Description:
	Evaluate number of effective thresholds.
@Algorithm:
	The first sequence of non-zero values in the given array
	determines the number of effective thresholds
@Arguments:
	float*	thr		threshold array
@Returns:
	int				number of effective thresholds
**************************************************************************/
int	res_thr_dim(float* thr)
{
	int nthr = 0;

	for (int i = 0; i<RES_NTHRESHOLDS; i++) {
		if (thr[i] > 0.) nthr++;
		else break;
	}

	return(nthr);
}
