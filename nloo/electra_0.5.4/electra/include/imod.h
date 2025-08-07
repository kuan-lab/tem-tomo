/*
	imod.h
	Header file for imod interface
	Author: Giovanni Cardone
	Created: 2004 	Modified: 20050628
*/

#ifndef _IMOD_H__
#define _IMOD_H__

#include "bsoft.h"
#include "ext_proc.h"
#include "string_tools.h"
#include "ioroutines.h"
#include <string.h>

// Constants
const int   FLEN = 256;
const int   XCLDN = 361;

//const char*  IPARSNAME="tilt.com";
const int   IPARALLEL=1;
const int   IMODE=1;
const float ILOGBASE=0.0;
const float IRADIAL_MAX=0.35;
const float IRADIAL_FALL=0.05;
const float ISCALE_FLEVL=1.39;
const float ISCALE_SCALE=1000.0;
const int   ISUB_NX=0;
const int   ISUB_NY=0;
const float IXAXISTILT=0.0;


// Structures
/*----------IMOD tilt parameters(s) Structure-------------*/
/************************************************************************
@Object: struct Imod_tilt
@Description:
	parameter structure for tilt program in imod.
@Features:
	This contains all the parameters needed to launch tilt application.
*************************************************************************/
struct Imod_tilt {
	char		parsname[FLEN];	// name of the tilt parametes file
	char		in_proj[FLEN];	// input projections file name
	char		out_vol[FLEN];	// output volume reconstruction
	unsigned int	nx;			// full size along x
	unsigned int	ny;			// full size along y
	int		setlog;			// apply log to input projected density
	float	logbase;		    // value to add before taking the logarithm of input density
	int		parallel;		// flag for producing sections parallel to the
							// plane of the zero tilt projection
	float	radial_max;		// maximum frequency where to apply radial weighting function (cycles/pixel)
	float	radial_fall;	    // std deviation of gaussian fall-off (cycles/pixel)
	int      mode;           // data mode of the output file (0=byte; 1=short; 2=float)
	int		setfbpinterp;    // order of interpolation flag
	int		fbpinterp;		// fast back-projection interpolation order
	int		setscale;		// scale tomogram flag
	float	scale_flevl;	    // bias to apply to density in the reconstructed image
	float	scale_scale;	    // scale factor to apply to density in the reconstructed image
	int		setshift;		// offset flag
	int		shift_x;		    // offset of final tomogram along x
	int		shift_z;		    // offset of final tomogram along z
	int		sub_nx;			// subset start along x
	int		sub_ny;			// subset start along y
	unsigned int	thickness;	// thickness (pixels)
	char		tiltfile[FLEN];	// file containing tilt angles
	int		exclude_list[XCLDN];   // list of excluded tilts from reconstruction
	int		nxl;			// number of actual excluded projections
	float	xaxistilt;		// common tilt around x (degrees)
	int		use_gpu;		// use gpu for reconstruction
} ;


// Function prototypes
Imod_tilt* imd_tilt_read(char* parsname);
int imd_tilt_parfile_exist(Imod_tilt * im_tilt);
int imd_tilt_checkvol(Imod_tilt * im_tilt);
char* imd_tilt_gettiltfile(Imod_tilt * im_tilt);
char* imd_tilt_getprojfile(Imod_tilt * im_tilt);
char* imd_tilt_getvolfile(Imod_tilt * im_tilt);
int imd_tilt_putvolfile(Imod_tilt * im_tilt, char* filename);
float imd_tilt_getxtilt(Imod_tilt * im_tilt);
float imd_tilt_getresolutionlimit(Imod_tilt * im_tilt);
int imd_tilt_getlog(Imod_tilt * im_tilt);
float imd_tilt_getlogbase(Imod_tilt * im_tilt);
int imd_tilt_readpars(Imod_tilt * im_tilt);
int imd_tilt_setdefault(Imod_tilt * im_tilt);
int imd_tilt_setparallel(Imod_tilt * im_tilt);
int imd_tilt_isparallel(Imod_tilt * im_tilt);
int imd_tilt_putpars(Imod_tilt * im_tilt, char * in_file, char * out_file, int nx, int ny, int nz, int exp_tag, float xtilt, char* tilt_file);
int imd_tilt_cpy_and_modify(Imod_tilt * im_tilt_dst, Imod_tilt* im_tilt_src, char * parsname,
	char * in_file, char * out_file, int nx, int ny, int nz, float xtilt, char* tilt_file);
int imd_tilt_cpy(Imod_tilt * im_tilt_dst, Imod_tilt* im_tilt_src);
int imd_tilt_writepars(Imod_tilt * im_tilt);
int imd_tilt_execute(Imod_tilt * im_tilt);
int imd_tilt_angle_missing(Imod_tilt * im_tilt, int a);
int	imd_tilt_add2xcldlist( Imod_tilt* im_tilt, int a);
int	imd_tilt_rmxcldlist( Imod_tilt* im_tilt);
int	imd_tilt_rmallfiles_butvol( Imod_tilt* im_tilt);
int	imd_tilt_rmalloutfiles( Imod_tilt* im_tilt);
int	imd_tilt_rmparfile( Imod_tilt* im_tilt);
int	imd_tilt_copyxcldlist(int* xcld_list, int* nxl, Imod_tilt* imd_tpars);
int	imd_tilt_makelist( int it, int np, int tstep, Imod_tilt* imd_pars, int* lp);
int imd_tilt_putdatatype(Imod_tilt * im_tilt, DataType dtype);

int imd_clip_flip(char* fmode, char* in_file, char* out_file);

#endif  // #ifndef _IMOD_H__
