/*
	img_extract.h
	Header file for functions to extract parts of images
	Author: Bernard Heymann
	Created: 20000430 	Modified: 20040406
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "matrix.h"
#include "rwimg.h"

// Function prototypes
Bimage*		img_extract(Bimage* p, int img_number, Vector3 origin, VectorInt3 size);
Bimage*		img_extract_tiles(Bimage* p, int img_number, VectorInt3 origin,
				VectorInt3 size, VectorInt3 tile_size);
Bimage*		img_extract_tiles_at_coords(Bimage* p, int img_number,
				int ntiles, VectorInt3* coords, VectorInt3 tile_size);
Bimage*		img_extract_for_show(Bimage* p, int img_num, int z, float scale);
Bimage*		img_extract_for_show2(Bimage* p, int img_num, int z, float scale);
Bimage*		img_extract_particles(Bimage* p, float radius, float bad_radius,
				int ncoord, int nbad, float* coords);
Bimage*		img_extract_filament(Bimage* p, int img_num, double width,
				int ncoord, Vector3* coords);
Bimage*		img_extract_panels(Bimage* p, int ncoord, float* coords);

