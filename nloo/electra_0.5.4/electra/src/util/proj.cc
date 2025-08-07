/*
	proj.cc
	Functions for handling projection structure
	Author: Giovanni Cardone
	Created: 20040606 	Modified: 20050628
*/

#include "proj.h"

/************************************************************************
@Function: proj_init
@Description:
	General projection structure initialization.
@Algorithm:
	This function allocates memory and sets up a number of defaults.
@Arguments:
	int					number of projections
						(0=do not initialize arrays dependent on this number)
@Returns:
	Proj*				new projection structure, NULL if initialization failed.
**************************************************************************/
Proj* 	proj_init(int np)
{
	// Allocate memory for the image parameter structure
	Proj* 	p = (Proj *) balloc(sizeof(Proj));
	if ( !p ) return(p);

	// Set parameter defaults
	p->img = NULL;
	p->dataflag = 0;
	p->x = p->y = 1;
	p->nimg = 0;
	p->nmax = np;
	p->datatype = UChar;
	p->ox = p->oy = p->oz = 0.;
	p->vx = p->vy = 0.;
	p->vz = 1.;
	p->angle = 0.;
	p->ux = p->uy = 1.;
	p->avg = 0.;
	p->std = 0.;
	p->tilt = NULL;
	p->xtilt = 0.;
	p->footprint = NULL;

	if (np > 0) {
		p->tilt = (float *) balloc(np*sizeof(float));
		p->img = (Bimage **) balloc(np*sizeof(Bimage*));
		p->footprint = (Shape2d **) balloc(np*sizeof(Shape2d*));
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG proj_init: Projection setup done\n");
	
	return(p);
}

/************************************************************************
@Function: proj_kill
@Description:
	General projection structure destruction.
@Algorithm:
	This function deallocates all memory associated to the structure
@Arguments:
	Proj*				projection structure
@Returns:
	int					error code (<0 means failure)
**************************************************************************/
int 	proj_kill(Proj* p)
{
	if ( p == NULL ) return(0);
	
	if ( p->tilt ) bfree(p->tilt, p->nmax*sizeof(float));
	
	if (p->img) {
		for (unsigned int i=0; i<p->nmax; i++)
			kill_img(p->img[i]);
		bfree(p->img,p->nmax*sizeof(Bimage*));
	}
		
	if (p->footprint) {
		for (unsigned int i=0; i<p->nmax; i++)
			shape2d_kill(p->footprint[i]);
		bfree(p->footprint,p->nmax*sizeof(Shape2d*));
	}

	bfree(p, sizeof(Proj));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG proj_kill memory = %ld\n", memory);
	
	return(0);
}

/************************************************************************
@Function: proj_display
@Description:
	Show general projection structure information.
@Algorithm:
@Arguments:
	Proj*				projection structure
@Returns:
	int					error code (<0 means failure)
**************************************************************************/
int 	proj_display(Proj* p)
{
	if ( p == NULL ) {
		printf("display: no projection allocated!\n");
		return(0);
	}
	
	if ( p->tilt ) {
		printf("display: tilt angles:\n");
		for ( unsigned int i=0; i<p->nmax; i++)
			printf(" %d: %f\n",i,p->tilt[i]);
	} else {
		printf("display: tilt angles not allocated.\n");
	}
	
	if ( p->img ) {
		printf("display: projections:\n");
		for (unsigned int i=0; i<p->nmax; i++)
			if (p->img[i])
				printf(" %d: allocated\n",i);
			else
				printf(" %d: not allocated\n",i);
	} else {
		printf("display: projections array not allocated.\n");
	}

	return(0);
}

/************************************************************************
@Function: proj_init_footprint
@Description:
	Calculates the footprint of the projections as the reprojection
	from the selected sub-volume shape
@Algorithm:
	Determines the area covered by the projection as a convex polygon or a circle,
	depending on the volume shape.
	The acquisition area is associated to another polygon.
	For each projection a different footprint is calculated.
@Arguments:
	Proj*		p		projection structure
	Bimage*		b		full volume to project
	Shape3d*		s		selected sub-volume
	int			taper	tapering flag
@Returns:
	int					<0 if error
**************************************************************************/
int 	proj_init_footprint(Proj* p, Vector3 origin, Shape3d* s, int taper)
{

	if (p==NULL || p->img==NULL) return(-1);

	View		view;
	for (unsigned int i=0; i<p->nmax; i++) {
		if (p->img[i]) {
			view = proj_getview(p, i);
			if (s==NULL)
				// footprint initialized to all the projection area	
				p->footprint[i] = shape2d_init_to_box(0, 0, (unsigned int) (p->img[i])->x, (unsigned int) (p->img[i])->y);
			else
				p->footprint[i] = shape_proj3dto2d(s, (unsigned int) (p->img[i])->x, (unsigned int) (p->img[i])->y, origin, view);
			shape2d_taper_init(p->footprint[i], taper);	
		}
	}
	
	return (0);
}

/************************************************************************
@Function: proj_copy_footprint
@Description:
	Copy the footprint from a projection structure
@Algorithm:
	All the footprints are copied. The number of rojections in the
	two structure must be the same.
@Arguments:
	Proj*		pout		modified projection structure (footprint destination)
	Proj*		pin		projection structure (footprint origin)
@Returns:
	int					<0 if error
**************************************************************************/
int 	proj_copy_footprint(Proj* pout, Proj* pin)
{

	if (pout==NULL || pin==NULL || pin->nmax!=pout->nmax) return(-1);

	for (unsigned int i=0; i<pin->nmax; i++) {
		if (pin->footprint[i])
			pout->footprint[i] = shape2d_copy(pin->footprint[i]);					
	}
	
	return(0);
}

/************************************************************************
@Function: proj_cpy_img_sett
@Description:
	Copy the genral setting from an image
@Algorithm:
@Arguments:
	Proj*				projection structure
	Bimage*				reference image 
@Returns:
	int					<0 if error
**************************************************************************/
int  	proj_cpy_img_sett(Proj* p, Bimage* b)
{
	
	if (p==NULL || b==NULL) return(-1);
	
	p->x = b->x;
	p->y = b->y;
	p->datatype = b->datatype;
	p->ox = 0.5*( b->x - 1.0 );
	p->oy = 0.5*( b->y - 1.0 );
	p->oz = 0.;
	p->ux = b->ux;
	p->uy = b->uy;
	p->avg = b->avg;
	p->std = b->std;
	p->min = b->min;
	p->max = b->max;
		
	return(0);
}

/************************************************************************
@Function: proj_puttilt
@Description:
	assign the tilt angle to the projection with the given index
@Algorithm:
@Arguments:
	Proj*				projection structure
	float				tilt angle
	int					projection index
@Returns:
	int					<0 if error
**************************************************************************/
int  	proj_puttilt(Proj* p, float tilt_angle, unsigned int idx)
{
	
	if (p==NULL || p->nmax<idx || p->tilt==NULL) return(-1);
	
	p->tilt[idx] = tilt_angle;
	
	return(0);
}

/************************************************************************
@Function: proj_getview
@Description:
	get the view vector assigned to the projection with the given index
@Algorithm:
@Arguments:
	Proj*				projection structure
	int					projection index
@Returns:
	View					view vector
**************************************************************************/
View  	proj_getview(Proj* p, unsigned int idx)
{
	
	View vw = view_from_4_values(0.,0.,1.,0.);
	
	if (p==NULL || p->nmax<idx || p->img[idx]==NULL || (p->img[idx])->image==NULL)
		return(vw);
	
	vw = view_from_4_values(	(p->img[idx])->image[0].vx,
								(p->img[idx])->image[0].vy,
								(p->img[idx])->image[0].vz,
								(p->img[idx])->image[0].angle);
	
	return(vw);
}

/************************************************************************
@Function: proj_putview
@Description:
	assign the view to the projection with the given index
@Algorithm:
@Arguments:
	Proj*				projection structure
	View					view vector
	int					projection index
@Returns:
	int					<0 if error
**************************************************************************/
int  	proj_putview(Proj* p, View v, unsigned int idx)
{
	
	if (p==NULL || p->nmax<idx || p->img[idx]==NULL) return(-1);
	
	(p->img[idx])->image[0].vx = v.x;
	(p->img[idx])->image[0].vy = v.y;
	(p->img[idx])->image[0].vz = v.z;
	(p->img[idx])->image[0].angle = v.a;
	
	return(0);
}

/************************************************************************
@Function: proj_putimage
@Description:
	assign the image data to the projection with the given index
@Algorithm:
	the projection sequence number is the index in the original projection
	series: if <0 it is assigned equal to the projection index
@Arguments:
	Proj*				projection structure
	Bimage*              image to assign
	int					projection index
	int					projection index number 
@Returns:
	int					<0 if error
**************************************************************************/
int  	proj_putimage(Proj* p, Bimage* b, unsigned int idx, int pin)
{
	
	if (p==NULL || p->nmax<idx || p->img==NULL) return(-1);
	
	if (p->img[idx] !=NULL) {
		if (verbose && VERB_PROCESS)
			printf("Warning(proj_putimage): relinking image pointer without freeing memory before!\n");
	} else
		p->nimg++;
	
	p->img[idx] = b;
	
	(p->img[idx])->i = pin<0?idx:pin;
	
	return(0);
}

/************************************************************************
@Function: proj_deleteimage
@Description:
	deallocate the image data from the projection with the given index
@Algorithm:
@Arguments:
	Proj*				projection structure
	int					projection index
@Returns:
	int					<0 if error
**************************************************************************/
int  	proj_deleteimage(Proj* p, unsigned int idx)
{
	
	if (p==NULL || p->nmax<idx || p->img==NULL) return(-1);
	
	if (p->img[idx] == NULL) {
		if (verbose && VERB_PROCESS)
			printf("Warning(proj_deleteimage): data not present!\n");
	} else {
		p->nimg--;
		kill_img(p->img[idx]);
		p->img[idx] = NULL;
	}
		
	if (p->footprint && p->footprint[idx])
		shape2d_kill(p->footprint[idx]);
	
	return(0);
}

/************************************************************************
@Function: proj_write_angles
@Description:
	write to file the tilt angles
@Algorithm:
	if the the number of tilt angles given as argument is zero,
	then determine the number from the allocated images
@Arguments:
	char*				output file
	Proj*				projection structure
	int					number of angles
@Returns:
	int					<0 if error
**************************************************************************/
int  	proj_write_angles(char* tfile, Proj* p, int ntlt)
{
	
	if (p==NULL || p->tilt==NULL) return(-1);
	
	int nt = ntlt;
	
	if (nt==0) 
		nt = (int) p->nimg;
//		while (nt<p->nmax && p->img[nt]) nt++;
	
	tlt_save_angles_to_file(p->tilt, nt, tfile);
	
	return(0);
}

/************************************************************************
@Function: proj_stats
@Description:
	Evaluate the statistics for the projection.
@Algorithm:
	Evaluate the statistics for each single projection in the given structure.
	Region on which evaluating the background is selected by means of bkg_type value.
	If a footprint is allocated, average and standard deviation are evaluated by
	only considering the points inside the footprint.
@Note:
	Calculation of average, standard deviation, min and max values
	for all the series of projections, without taking into account 
	any footprint, is not implemented.
@Arguments:
	Proj*	p			projection structure
	int		bkg_type		background evaluation method
@Returns:
	int					<0 if error
**************************************************************************/
int 	proj_stats(Proj* p, int bkg_type)
{
	if (p==NULL || p->img==NULL) return(-1);
	

	for (unsigned int i=0; i<p->nmax; i++) {
		if (p->img[i])
			shape2d_image_stats(p->img[i], p->footprint[i], bkg_type);
	}
	
	return(0);
}

/************************************************************************
@Function: proj_mask
@Description:
	mask a set of projections.
@Algorithm:
	The background value is applied to the volume outside the footprint.
	Tapering, if related flag is set into the footprint configuration, is applied.
	If related footprint is not allocated, the projection is not masked	
@Arguments:
	Proj*  p			projection.
@Returns:
	int				<0 if not successfull.
**************************************************************************/
int 	proj_mask(Proj* p)
{
	if (p==NULL || p->img==NULL) return(-1);
	
	for (unsigned int i=0; i<p->nmax; i++) {
		if (p->img[i])
			shape2d_mask_image(p->img[i], p->footprint[i],(p->footprint[i])->taper);
	}
	
	return(0);
}

/************************************************************************
@Function: proj_to_datatype
@Description:
	Convert projections to datatype.
@Algorithm:
@Arguments:
	Proj*		p		projection structure
	DataType		dt		new data type
@Returns:
	int					<0 if error
**************************************************************************/
int 	proj_to_datatype(Proj* p, DataType dt)
{
	if (p==NULL || p->img==NULL) return(-1);
	
	p->datatype = dt;
	
	for (unsigned int i=0; i<p->nmax; i++) {
		if (p->img[i]) {
			img_stats(p->img[i]);
			img_to_datatype(p->img[i], dt);
		}
	}
	
	return(0);
}

/************************************************************************
@Function: proj_rescale_to_avg_std
@Description:
	Rescales the projections to a given average and standard deviation.
@Algorithm:
	Rescaling is applied only inside the footprint.
	If given, a background value is applied outside the footprint.
	The proutine norm2d_rescale_to_avg_std is called for each projection.
@Note:
	Image statistics are not recalculated.
@Arguments:
	Proj* p			projection.
	float avg		new average.
	float std		new standard deviation.
	float* background	background value (optional).
@Returns:
	int 			0.
**************************************************************************/
int proj_rescale_to_avg_std(Proj* p, float avg, float std, float* background)
{
	if (p==NULL || p->img==NULL) return(-1);
	
	for (unsigned int i=0; i<p->nmax; i++) {
		if (p->img[i])
			norm2d_rescale_to_avg_std(p->img[i], avg, std, background, p->footprint[i]);
	}

	return(0);
}

/************************************************************************
@Function: proj_rescale_linregression
@Description:
	normalize a projection by least squares apporach (linear regression).
@Algorithm:
	For each projection the routine proj_rescale_linregression is called
@Arguments:
	Proj* p		input projection - modified.
	Proj* pref	reference projection.
@Returns:
	int				<0 if not successfull.
**************************************************************************/
int		proj_rescale_linregression(Proj* p, Proj* pref)
{
	if (p==NULL || p->img==NULL || pref==NULL || pref->img==NULL) return(-1);
	
	for (unsigned int i=0; i<p->nmax; i++) {
		if (p->img[i] && pref->img[i])
			norm2d_rescale_linregression(p->img[i], pref->img[i], p->footprint[i]);
	}

	return(0);
}
