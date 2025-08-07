/*
	vnloo.h
	Header file for volumetric resolution by nloo
	Author: Giovanni Cardone
	Created: 20050411 	Modified: 
*/

#ifndef _VNLOO_H__
#define _VNLOO_H__

#include "bsoft.h"
#include "ioroutines.h"
#include "proj.h"

// Constants

// Structures
/************************************************************************
@Object: struct Vnloo
@Description:
	structure for volumetric resolution measure.
@Features:
	This contains all the auxiliary arrays needed to evaluate the integrated resolution.
*************************************************************************/
struct Vnloo {
	char*        fprefix;			// file prefix
	char* 		fvres;			// volumetric resolution file
	char* 		fmissinput;		// partial reproj vs input reproj file
	char* 		fmiss; 			// partial reproj file
	char* 		ffullinput;		// full reproj vs input reproj file
	char* 		ffull;			// full reproj file
	unsigned long nx, ny, nz;		// volume size
	unsigned int  bin;			// volume binning
	float		zb;				// further binning factor along z 
	float* 		vres;
	float* 		missinput;
	float* 		miss; 
	float* 		fullinput;
	float* 		full;
	VectorInt3*  vcrd;			// volume coords correspondent to slice coords
} ;

// Function prototypes
Vnloo* vnloo_init(char* pfx_file, Bimage* b, int vbin);
int vnloo_kill( Vnloo* vnl);
int vnloo_reset_tmparr( Vnloo* vnl);
int vnloo_put(Vnloo* vnl, int x, int y, float Armdp, float Arfdp, float A2m, float A2f);
int vnloo_add_slice(Vnloo* vnl, Bimage* b);
int vnloo_compute(Vnloo* vnl);
#endif  // #ifndef _VNLOO_H__
