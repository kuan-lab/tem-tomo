/*
	img_find.h 
	Header file for functions to find a template in an image
	Author: Bernard Heymann
	Created: 20021027 	Modified: 20041030
*/

#include "rwimg.h"
#include "matrix.h"

// Function prototypes
float		img_find_density(Bimage* p, Bimage* ptemp, View* view,
				float alpha, float alpha_step, float hires, float lores,
				Bimage* pmask, View* currview, Vector3* currshift);
View*		refine_views(View theview, float alpha_step);

