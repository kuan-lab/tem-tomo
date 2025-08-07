/*
	img_pick.c
	Library routines to pick particles in images
	Author: Bernard Heymann
	Created: 20030428 	Modified: 20041212
*/

#include "img_pick.h"
//#include "img_feature.h"
#include "img_fourier.h"
#include "img_rescale.h"
#include "img_edit.h"
//#include "img_filter.h"
#include "img_datatypes.h"
#include "img_util.h"
#include "rwimg.h"
//#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated

struct LVectorInt3 {
	LVectorInt3*	next;
	int				x, y, z;
} ;

/************************************************************************
@Function: img_find_particles
@Description:
	Finds particles using a Hough transform with a simple ring.
@Algorithm:
	A ring defined by inner and outer radii is used to find circular
	shapes within an image.
@Arguments:
	Bimage* p				image.
	float inner_radius		inner radius of ring.
	float outer_radius		outer radius of ring.
	double* threshold		pointer to threshold to accept peaks.
	unsigned long* npart	pointer to number of praticles returned.
@Returns:
	VectorInt3* 			an array of coordinates.
**************************************************************************/
/*VectorInt3*		img_find_particles(Bimage* p, float inner_radius, float outer_radius, 
					double* threshold, unsigned long* npart)
{
	if ( outer_radius < inner_radius ) outer_radius = inner_radius;
	
	unsigned long   i, x, y, z, n;
	long			in, ix, iy, iz;
	double			dx2, dy2, dz2, r2;
	double			inrad2 = inner_radius*inner_radius;
	double			outrad2 = outer_radius*outer_radius;
	
	long			intrad = (long) (outer_radius*2 + 1.5);
	VectorLong3		size = {intrad, intrad, intrad};
	if ( size.x > p->x ) size.x = p->x;
	if ( size.y > p->y ) size.y = p->y;
	if ( size.z > p->z ) size.z = p->z;
	VectorLong3		origin = {size.x/2, size.y/2, size.z/2};
	
	img_filter_mg_extremes(p);
	
	if ( verbose & VERB_PROCESS ) {
		printf("Finding particles:\n");
		printf("Radii:                          %g - %g\n", inner_radius, outer_radius);
	}

	Bimage* 		ptf = init_img_with_parameters(Float, 1, p->x, p->y, p->z, p->n);
	
//	unsigned long   datasize = p->n*p->x*p->y*p->z*p->c;
	float*			data = (float *) p->data;
	float*			tf = (float *) ptf->data;
	LVectorInt3*	v = NULL;
	LVectorInt3*	fv = NULL;
	
	for ( i=z=0; z<size.z; z++ ) {
		iz = (long)z - origin.z;
		dz2 = (double)iz*iz;
		for ( y=0; y<size.y; y++ ) {
			iy = (long)y - origin.y;
			dy2 = (double)iy*iy;
			for ( x=0; x<size.x; x++, i++ ) {
				ix = (long)x - origin.x;
				dx2 = (double)ix*ix;
				r2 = dx2 + dy2 + dz2;
				if ( r2 >= inrad2 && r2 <= outrad2 ) {
					v = (LVectorInt3 *) add_item((char **) &v, sizeof(LVectorInt3));
					if ( !fv ) fv = v;
					v->x = ix;
					v->y = iy;
					v->z = iz;
				}
			}
		}
	}
	
//	if ( verbose & VERB_DEBUG ) {
//		printf("DEBUG img_find_particles: radius = %g  edge = %g\n", particle_radius, gauss_edge);
//		write_img("particle.map", p);
//	}

    for ( i=n=0; n<p->n; n++ ) {
		in = n*p->z;
	    for ( z=0; z<p->z; z++ ) {
    		for ( y=0; y<p->y; y++ ) {
				for ( x=0; x<p->x; x++,i++ ) {
					for ( v=fv; v; v=v->next ) {
						iz = v->z + (long)z;
						if ( iz >= 0 && iz < p->z ) {
							iz = (in + iz)*p->y;
							iy = v->y + (long)y;
							if ( iy >= 0 && iy < p->y ) {
								iy = (iz + iy)*p->x;
								ix = v->x + (long)x;
								if ( ix >= 0 && ix < p->x ) {
									ix = iy + ix;
									tf[i] -= data[ix];
								}
							}
						}
					}
				}
			}
		}
	}
	
	kill_list((char *) fv, sizeof(LVectorInt3));
	
	img_stats(ptf);
	
	if ( verbose & VERB_FULL ) {
		printf("Writing transform map to tf.map\n\n");
		write_img("tf.map", ptf);
	}
	
	if ( *threshold < 0 ) *threshold = ptf->avg + 2 * ptf->std;
	
	VectorInt3*		coord = img_find_particles_in_ccmap(ptf, outer_radius, threshold, npart);
	
	kill_img(ptf);
	
	if ( verbose & VERB_LABEL )
		printf("Number of particles found = %ld\n", *npart);
	
	return(coord);
}
*/
/************************************************************************
@Function: img_find_particles
@Description:
	Finds particles using a Hough transform with a simple ring.
@Algorithm:
	A ring defined by inner and outer radii is used to find circular
	shapes within an image.
@Arguments:
	Bimage* p				image.
	float inner_radius		inner radius of ring.
	float outer_radius		outer radius of ring.
	double* threshold		pointer to threshold to accept peaks.
	unsigned long* npart	pointer to number of praticles returned.
@Returns:
	VectorInt3* 			an array of coordinates.
**************************************************************************/
/*VectorInt3*		img_find_particles2(Bimage* p, float inner_radius, float outer_radius, 
					double* threshold, unsigned long* npart)
{
	if ( outer_radius < inner_radius ) outer_radius = inner_radius;
	
	unsigned long   i, x, y, z, n, xx, yy, zz;
	unsigned long	in, ix, iy, iz, jx, jy, jz;
	double			dx2, dy2, dz2, r2;
	double			inrad2 = inner_radius*inner_radius;
	double			outrad2 = outer_radius*outer_radius;
	
	long			intrad = (long) (outer_radius*2 + 1.5);
	VectorLong3		size = {intrad, intrad, intrad};
	if ( size.x > p->x ) size.x = p->x;
	if ( size.y > p->y ) size.y = p->y;
	if ( size.z > p->z ) size.z = p->z;
	VectorLong3		origin = {size.x/2, size.y/2, size.z/2};
	VectorLong3		lo, hi, start;
	
	img_filter_mg_extremes(p);

	if ( verbose & VERB_PROCESS ) {
		printf("Finding particles:\n");
		printf("Radii:                          %g - %g\n", inner_radius, outer_radius);
	}

	Bimage* 		ptf = init_img_with_parameters(Float, 1, p->x, p->y, p->z, p->n);
	
//	unsigned long   datasize = p->n*p->x*p->y*p->z*p->c;
	float*			data = (float *) p->data;
	float*			tf = (float *) ptf->data;
	char*			ring = (char *) balloc(size.x*size.y*size.z*sizeof(char));
	
	for ( i=z=0; z<size.z; z++ ) {
		dz2 = (double)z - origin.z;
		dz2 *= dz2;
		for ( y=0; y<size.y; y++ ) {
			dy2 = (double)y - origin.y;
			dy2 *= dy2;
			for ( x=0; x<size.x; x++, i++ ) {
				dx2 = (double)x - origin.x;
				dx2 *= dx2;
				r2 = dx2 + dy2 + dz2;
				if ( r2 >= inrad2 && r2 <= outrad2 )
					ring[i] = 1;
			}
		}
	}
	
//	if ( verbose & VERB_DEBUG ) {
//		printf("DEBUG img_find_particles: radius = %g  edge = %g\n", particle_radius, gauss_edge);
//		write_img("particle.map", p);
//	}

    for ( i=n=0; n<p->n; n++ ) {
	    for ( z=0; z<p->z; z++ ) {
			lo.z = 0;
			hi.z = size.z;
			start.z = (long)z - origin.z;
			if ( start.z < 0 ) lo.z = -start.z;
			if ( start.z + size.z > p->z ) hi.z = p->z - start.z;
			in = n*p->z + start.z;
    		for ( y=0; y<p->y; y++ ) {
				lo.y = 0;
				hi.y = size.y;
				start.y = (long)y - origin.y;
				if ( start.y < 0 ) lo.y = -start.y;
				if ( start.y + size.y > p->y ) hi.y = p->y - start.y;
    	    	for ( x=0; x<p->x; x++,i++ ) {
					lo.x = 0;
					hi.x = size.x;
					start.x = (long)x - origin.x;
					if ( start.x < 0 ) lo.x = -start.x;
					if ( start.x + size.x > p->x ) hi.x = p->x - start.x;
					for ( zz=lo.z; zz<hi.z; zz++ ) {
						iz = (in + zz)*p->y + start.y;
						jz = zz*size.y;
						for ( yy=lo.y; yy<hi.y; yy++ ) {
							iy = (iz + yy)*p->x + start.x;
							jy = (jz + yy)*size.x;
							for ( xx=lo.x; xx<hi.x; xx++ ) {
								ix = iy + xx;
								jx = jy + xx;
								if ( ring[jx] ) {
									tf[i] -= data[ix];
								}
							}
						}
					}
				}
			}
		}
	}
	
	bfree(ring, size.x*size.y*size.z*sizeof(char));
	
	img_stats(ptf);
	
//	if ( verbose & VERB_FULL ) {
		printf("Writing transform map to tf.map\n\n");
		write_img("tf.map", ptf);
//	}
	
	if ( *threshold < 0 ) *threshold = ptf->avg + 2 * ptf->std;
	
	VectorInt3*		coord = img_find_particles_in_ccmap(ptf, outer_radius, threshold, npart);
	
	kill_img(ptf);
	
	if ( verbose & VERB_LABEL )
		printf("Number of particles found = %ld\n", *npart);
	
	return(coord);
}
*/
/************************************************************************
@Function: img_pick_particles
@Description:
	Picks particles using cross-correlation and segmentation.
@Algorithm:
	A flat gaussian template with the given radius and gaussian edge is
	generated and cross-correlated with the input image including
	bandpass filtering to target the size of the template particle.
@Arguments:
	Bimage* p				image.
	float particle_radius	radius of particle.
	float gauss_edge		width of edge smoothing.
	double* threshold		pointer to cross-correlation threshold to accept peaks.
	unsigned long* npart	pointer to number of praticles returned.
@Returns:
	VectorInt3* 			an array of coordinates.
**************************************************************************/
/*VectorInt3*		img_pick_particles(Bimage* p, float particle_radius, float gauss_edge,
					double* threshold, unsigned long* npart)
{
	int				n;
	
	Bimage* 		ppart = init_img_with_parameters(Float, 1, p->x, p->y, p->z, p->n);
	
	Vector3 		origin = {p->x/2,p->y/2,p->z/2};

	img_delete_within_radius(ppart, origin, particle_radius+gauss_edge, gauss_edge, 1);
	
	img_delete_within_radius(ppart, origin, particle_radius, gauss_edge, -1);
	
	for ( n=0; n<ppart->n; n++ ) {
		ppart->image[n].ox = origin.x;
		ppart->image[n].oy = origin.y;
		ppart->image[n].oz = origin.z;
	}

	ppart->ux = p->ux;
	ppart->uy = p->uy;
	ppart->uz = p->uz;
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG img_pick_particles: radius = %g  edge = %g\n", particle_radius, gauss_edge);
		write_img("particle.map", p);
	}

	img_filter_mg_extremes(p);
	
	Bimage*			pcc = img_cross_correlate(p, ppart, 4*p->ux, 3*particle_radius*p->ux, NULL);
	
	VectorInt3		shift;
	shift.x = (int) ppart->image[0].ox;
	shift.y = (int) ppart->image[0].oy;
	shift.z = (int) ppart->image[0].oz;
	
	img_translate(pcc, shift, 1, 0);
	
	img_stats(pcc);
	
	if ( verbose & VERB_FULL ) {
		printf("Writing CC map to part.map\n");
		write_img("part.map", ppart);
		printf("Writing CC map to cc.map\n\n");
		write_img("cc.map", pcc);
	}
	
	kill_img(ppart);
	
	if ( *threshold < 0 ) *threshold = pcc->avg + 2 * pcc->std;
	
	VectorInt3*		coord = img_find_particles_in_ccmap(pcc, particle_radius, threshold, npart);
	
	kill_img(pcc);
	
	if ( verbose & VERB_LABEL )
		printf("Number of particles found = %ld\n", *npart);
	
	return(coord);
}
*/
/************************************************************************
@Function: img_find_particles_in_ccmap
@Description:
	Finds the peaks in a cross-correlation map corresponding to particles.
@Algorithm:
	An image of the size of an input sub-image is generated with a
	particle (black or negative density) at its center.
	Each sub-image is cross-correlated with the particle image.
@Arguments:
	Bimage* pcc				cross-correlation map.
	float particle_radius	low resolution limit.
	double* threshold		pointer to cross-correlation threshold to accept peaks.
	unsigned long* ncoord	number of coordinates returned.
@Returns:
	VectorInt3* 			array of coordinates.
**************************************************************************/
VectorInt3*		img_find_particles_in_ccmap(Bimage* pcc, float particle_radius, 
					double* threshold, unsigned long* ncoord)
{
	unsigned long	rad = (unsigned long) particle_radius;
	VectorInt3		tile = {(int) pcc->x,(int) pcc->y,(int) pcc->z};
	VectorInt3		start = {(int) rad,(int) rad,(int) rad};
	VectorInt3		end = {(int) (pcc->x-rad),(int) (pcc->y-rad),(int) (pcc->z-rad)};
	tile = vectorint3_scalar_min(tile, (int) rad);
	if ( start.x >= pcc->x ) start.x = 0;
	if ( start.y >= pcc->y ) start.y = 0;
	if ( start.z >= pcc->z ) start.z = 0;
	end = vectorint3_scalar_max(end, 1);
	
	if ( verbose & VERB_LABEL )
		printf("Finding peaks\n\n");
	
	unsigned long	i, j, k, l, nc, ncp, x, y, z, xx, yy, zz, n;
	unsigned long	nt = ((pcc->x-1)/tile.x+1)*((pcc->y-1)/tile.y+1)*((pcc->z-1)/tile.z+1);
	float			max;
	VectorInt3*		ct = (VectorInt3 *) balloc(nt*sizeof(VectorInt3));
	
	float*			data = (float *) pcc->data;
	
	for ( n=0, nc=0; n<pcc->n; n++ ) {
		for ( z=start.z; z<end.z; z+=rad ) {
			for ( y=start.y; y<end.y; y+=rad ) {
				for ( x=start.x; x<end.x; x+=rad ) {
					max = -1e30;
					for ( zz=z; zz<z+tile.z && zz<pcc->z; zz++ ) {
						for ( yy=y; yy<y+tile.y && yy<pcc->y; yy++ ) {
							for ( xx=x; xx<x+tile.x && xx<pcc->x; xx++ ) {
								i = ((n*pcc->z + zz)*pcc->y + yy)*pcc->x + xx;
								if ( data[i] > *threshold && max < data[i] ) {
									max = data[i];
									ct[nc].x = (int) xx;
									ct[nc].y = (int) yy;
									ct[nc].z = (int) zz;
								}
							}
						}
					}
					if ( max > *threshold ) nc++;
				}
			}
		}
	}

	for ( k=0, ncp=0; k<nc; k++ ) {
		for ( l=0; l<nc; l++ ) if ( k!=l ) {
			if ( vectorint3_length(vectorint3_subtract(ct[k], ct[l])) < rad ) {
				i = (ct[k].z*pcc->y + ct[k].y)*pcc->x + ct[k].x;
				j = (ct[l].z*pcc->y + ct[l].y)*pcc->x + ct[l].x;
				if ( data[i] > data[j] ) ct[l].x = -1000;
				else ct[k].x = -1000;
			}
		}
		if ( ct[k].x > 0 ) ncp++;
	}
	
	VectorInt3*		coord = NULL;
	if ( ncp > 0 ) {
		coord = (VectorInt3 *) balloc(ncp*sizeof(VectorInt3));
		for ( i=0, j=0; i<nc; i++ ) if ( ct[i].x > 0 ) coord[j++] = ct[i];
	}
	
	bfree(ct, nt*sizeof(VectorInt3));
	
	*ncoord = ncp;
	
	return(coord);
}

