/*
	resol.h
	Header file for resolution routines
	Author: Giovanni Cardone
	Created: 20050112 	Modified: 
*/

#ifndef _RESOL_H__
#define _RESOL_H__

#include "bsoft.h"
#include "wedge.h"
#include "vnloo.h"

// Constants
#define SMALLFLOAT  1e-37

const int   RES_NTHRESHOLDS=3;  // number of resolution thresholds 

// Structures
/************************************************************************
@Object: struct Ires_nloo3d
@Description:
	structure for integrated resolution measure.
@Features:
	This contains all the auxiliary arrays needed to evaluate the integrated resolution.
*************************************************************************/
struct Ires_nloo3d {
	int			n;				// arrays size
	int			np;				// number of projections
	int			nthr;			// number of thresholds
	float		rad_scale;      // resolution sampling
	float		thresh[RES_NTHRESHOLDS];	//resolution thresholds
	float		estimate[RES_NTHRESHOLDS];	//resolution estimations
	double* 		nloo;
	double* 		nloo_missinput;
	double* 		nloo_miss; 
	double* 		nloo_input; 
	double* 		nloo_fullinput;
	double* 		nloo_full;
} ;

// Function prototypes
float* 	res_fsc_estimate(Bimage* p, Bimage* pmod, int bin_size, Bimage* pmask,
		  char* txtfile, float* res_crit, Wedge* mss_wdg, int fsceo, Bimage* vres);
Ires_nloo3d* res_nloo3d_init(Bimage* p, float* res_crit);
int 		res_nloo3d_estimate(Ires_nloo3d* ir, char* curve_file);
int 		res_nloo3d_kill(Ires_nloo3d* ir);
float* 	res_nloo2d_estimate(int pn, float tilt_ang, Bimage* pref, Bimage* pmiss, Bimage* pfull, int bin_size, Bimage* pmask, char* txtfile, float * res_thr, Ires_nloo3d* int_res);
float* 	res_rd_nloo2d_estimate(Bimage* pref, Bimage* pmiss, Bimage* pfull, int bin_size, Bimage* pmask, char* txtfile, float * res_thr, Ires_nloo3d* int_res, Vnloo* vnl);
int 		res_axis_dim(Bimage* p, unsigned int* rad_size, float* rad_sampl, float* max_res);
int		res_thr_dim(float* thr);

#endif  // #ifndef _RESOL_H__
