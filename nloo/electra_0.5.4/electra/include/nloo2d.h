/*
	nloo2d.h
	Header file for nloo2d evaluations
	Author: Giovanni Cardone
	Created: 20050317 	Modified: 
*/

#ifndef _NLOO2D_H__
#define _NLOO2D_H__

#include "bsoft.h"

// Constants

// Structures
/************************************************************************
@Object: struct N2d_curve
@Description:
	structure for nloo2d curves.
@Features:
	This contains all the auxiliary arrays needed for defining nloo2d
	resolution curves.
*************************************************************************/
struct N2d_curve {
	int			nc;				// number of curves
	int			ns;				// number of radial frequencies
    float*       tlt_angs;		// tilt angles array (size: nc)
	float*		axis;            // radial frequency array (size: ns)
	float*		res;             // resolution array
	float*		data;            // correlation curves (size: ns x nc)
};

// Function prototypes
N2d_curve* n2d_init(int nc, int ns);
N2d_curve* n2d_read_file(char* fname);
//n2d_load
//n2d_average
int n2d_kill(N2d_curve* n2d);
//n2d_write_file

#endif  // #ifndef _NLOO2D_H__
