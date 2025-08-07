/*
	img_pick.h
	Header file for functions to pick particles in images
	Author: Bernard Heymann
	Created: 20030428 	Modified: 20041212
*/

#include "rwimg.h"

// Function prototypes
VectorInt3*		img_find_particles(Bimage* p, float inner_radius, float outer_radius, 
					double* threshold, unsigned long* npart);
VectorInt3*		img_pick_particles(Bimage* p, float particle_radius, float gauss_edge,
					double* threshold, unsigned long* npart);
VectorInt3*		img_find_particles_in_ccmap(Bimage* pcc, float particle_radius, 
					double* threshold, unsigned long* ncoord);

