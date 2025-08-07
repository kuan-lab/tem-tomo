/*
	resol.h
	Header file for resolution routines
	Author: Giovanni Cardone
	Created: 20060613 	Modified: 
*/

#ifndef _NLOO_H__
#define _NLOO_H__

#include "bsoft.h"
//#include "resol.h"

// Function prototypes
float* 	nloo_corr_estimate(Bimage* pref, Bimage* pmiss, Bimage* pfull, float lores, float hires);

#endif  // #ifndef _NLOO_H__
