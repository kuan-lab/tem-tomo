/*
	img_mask.h
	Header file for binary mask operation tools
	Author: Samuel Payne and Bernard Heymann
	Created: 20010710 	Modified: 20040308 (BH)
*/

#include <stdio.h>
#include <stdlib.h>
#include "rwimg.h"

//Function prototypes
int			img_mask_invert(Bimage* p);
int			img_mask_combine(Bimage* p1, Bimage* p2, int operation);
Bimage* 	img_thresholded_mask(Bimage* p, float threshold);
Bimage*		img_min_max_mask(Bimage* p, float min, float max);
int			img_rectangular_mask(Bimage* p, float length, float width, 
				float rect_angle, int wrap, int operation);
int			img_shell_mask(Bimage* p, Vector3 origin, float rad_min, 
				float rad_max, int wrap, int operation);
int			img_missing_wedge_mask(Bimage* p, Vector3 origin, 
				float x_ang_neg, float x_ang_pos, float y_ang_neg, float y_ang_pos, 
				int wrap, int operation);
int			img_missing_cone_mask(Bimage* p, Vector3 origin, float mis_ang, 
				int wrap, int operation);
int 		img_apply_binary_mask(Bimage* p, Bimage* pmask, float fill);
int 		img_weigh_with_mask(Bimage* p, Bimage* pmask, float fill);
int 		img_mask_difference(Bimage* pmask1, Bimage* pmask2);
int 		img_mask_dilate(Bimage* pmask, int times);
int 		img_mask_erode(Bimage* pmask, int times);
int 		img_mask_write(Bimage* pmask, char* filename);
