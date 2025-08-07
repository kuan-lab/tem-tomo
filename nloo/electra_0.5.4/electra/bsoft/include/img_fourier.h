/*
	img_transform.h
	General FFT for n-dimensional data
		Implementing the FFTW library
	Author: Bernard Heymann
	Created: 19980805  	    Modified: 20041208
*/

#include <fftw.h>

#include "rwimg.h"
#include "img_complex.h"
#include "matrix.h"
#include "utilities.h"

#define MAX_RANK 	3	// Only up to 3D transforms

/* Function prototypes */
int 		img_fft(fftw_direction dir, Bimage* p);
int 		img_fft_complex(fftw_direction dir, Bimage* p, int norm_flag);
int 		img_fft_times(int ndim, int minsize, int maxsize);
int 		img_auto_correlate(Bimage* p);
Bimage* 	img_cross_correlate(Bimage* p1, Bimage* p2, float hires, 
				float lores, Bimage* pmask);
Vector3* 	img_find_shift(Bimage* p1, Bimage* p2, Bimage* pmask, float hires, 
				float lores, float radius, int refine_flag);
Vector3*	img_find_peak(Bimage* p, Vector3 origin, float radius);
int			img_refine_peak(Bimage* pc, Vector3* shift);
int 		img_change_transform_size(Bimage* p, VectorInt3 size);
Bimage* 	img_phase_difference(Bimage* p1, Bimage* p2, float resolution_hi, float resolution_lo);
float 		img_average_phase_difference(Bimage* p1, Bimage* p2, 
				float resolution_hi, float resolution_lo, int weighting);
int			img_flip_phases(Bimage* p, Bimage* pd);
