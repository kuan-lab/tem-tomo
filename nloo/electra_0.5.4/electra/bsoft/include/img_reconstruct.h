/*
	img_reconstruct.h
	Header file for 2D and 3D reconstruction
	Authors: Bernard Heymann
	Created: 20010403	Modified: 20041105
*/

#include "rwimg.h"
#include "mg_processing.h"

// constants
const int INT_ZPWNN = 1, INT_WNN = 2, INT_NN = 3; 

// Function prototypes: Note: rweigh = threshold
Bimage* 	myimg_fourier_reconstruction(Bproject* project, Bsymmetry* sym,
				int select, float hi_res, float scale, VectorInt3 size, int twoD_flag, int interp_type);
Bimage* 	myimg_fourier_reconstruction_fact(Bproject* project, Bsymmetry* sym,
				int select, float hi_res, float scale, VectorInt3 size, int twoD_flag); 
//float		img_pack_2D_in_recip_space(Bimage* p, Bimage* pmap, int* num,
//				float* weight, Bsymmetry* sym, float hi_res, float scale);
float		img_pack_2D_in_recip_space(Bimage* p, Bimage* pmap, int* num,
				float* weight, float* weight2, Bsymmetry* sym, float hi_res, float scale, int interp_type);

