/* 
	vnloo.cc
	Functions to handle volumetric resolution by nloo
	Author: Giovanni Cardone
	Created: 20050411	Modified:
*/ 

#include "vnloo.h"
int vround (double value);

/************************************************************************
@Function: vnloo_init
@Description:
	initialize a Vnloo structure and generate related files.
@Arguments:
	char*	pfx_file		prefix file name
	Bimage*	b 			reconstruction volume file
	int		vbin			vnloo volume binning
@Returns:
	Vnloo*  Vnloo struct initialized
**************************************************************************/
Vnloo* vnloo_init(char* pfx_file, Bimage* b, int vbin)
{
	Vnloo*  vnl = (Vnloo *) balloc(sizeof(Vnloo));

	vnl->fprefix = filename_change_type(pfx_file,"mrc");
	vnl->fvres = insert_in_filename(vnl->fprefix,"vnloo",'_');
	vnl->fmissinput = insert_in_filename(vnl->fprefix,"vnl_mi",'_');
	vnl->fmiss = insert_in_filename(vnl->fprefix,"vnl_m",'_');
	vnl->ffullinput = insert_in_filename(vnl->fprefix,"vnl_fi",'_');
	vnl->ffull = insert_in_filename(vnl->fprefix,"vnl_f",'_');
	vnl->bin = (unsigned int) vbin;
	vnl->zb = (float) b->x / (float) b->z;
	vnl->nx = (b->x - 1) / vbin + 1;
	vnl->ny = (b->y - 1) / vbin + 1;
	vnl->nz = (b->z - 1) / vbin + 1;
	vnl->vres = (float *) balloc(vnl->nx*vnl->ny*sizeof(float));
	vnl->missinput = (float *) balloc(vnl->nx*vnl->ny*sizeof(float));
	vnl->miss = (float *) balloc(vnl->nx*vnl->ny*sizeof(float));
	vnl->fullinput = (float *) balloc(vnl->nx*vnl->ny*sizeof(float));
	vnl->full = (float *) balloc(vnl->nx*vnl->ny*sizeof(float));
	vnl->vcrd = (VectorInt3 *) balloc(vnl->nx*vnl->ny*sizeof(VectorInt3));
		
	Proj* p = proj_init(1);
	proj_cpy_img_sett(p, b);
	p->x = vnl->nx;
	p->y = vnl->ny;
	proj_to_datatype(p, Float);
	io_pwrite(vnl->fvres,p,IO_ZEROS,(unsigned int)vnl->nz);
	io_pwrite(vnl->fmissinput,p,IO_ZEROS,(unsigned int)vnl->nz);
	io_pwrite(vnl->fmiss,p,IO_ZEROS,(unsigned int)vnl->nz);
	io_pwrite(vnl->ffullinput,p,IO_ZEROS,(unsigned int)vnl->nz);
	io_pwrite(vnl->ffull,p,IO_ZEROS,(unsigned int)vnl->nz);

	proj_kill(p);

	return(vnl);
}

/************************************************************************
@Function: vnloo_kill
@Description:
	General Vnloo structure destruction and removal of related files.
@Algorithm:
	This function deallocates all memory associated to the structure
	and remove the auxiliary files containing temporary values
@Arguments:
	Vnloo*  vnl		Vnloo structure
@Returns:
	int		error code (<0 means failure)
**************************************************************************/
int vnloo_kill( Vnloo* vnl)
{
	
	if ( vnl == NULL ) return(0);
	
	remove(vnl->fmissinput);
	remove(vnl->fmiss);
	remove(vnl->ffullinput);
	remove(vnl->ffull);

	bfree(vnl->vres,vnl->nx*vnl->ny*sizeof(float));
	bfree(vnl->missinput, vnl->nx*vnl->ny*sizeof(float));
	bfree(vnl->miss, vnl->nx*vnl->ny*sizeof(float));
	bfree(vnl->fullinput, vnl->nx*vnl->ny*sizeof(float));
	bfree(vnl->full, vnl->nx*vnl->ny*sizeof(float));
	bfree(vnl->vcrd, vnl->nx*vnl->ny*sizeof(VectorInt3));

	bfree_string(vnl->fprefix);
	bfree_string(vnl->fvres);
	bfree_string(vnl->fmissinput);
	bfree_string(vnl->fmiss);
	bfree_string(vnl->ffullinput);
	bfree_string(vnl->ffull);

	bfree(vnl, sizeof(Vnloo));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG vnloo_kill memory = %ld\n", memory);

	return(0);
}

/************************************************************************
@Function: vnloo_reset_tmparr
@Description:
	set arrays to zero.
@Algorithm:
@Arguments:
	Vnloo*  vnl		Vnloo structure
@Returns:
	int		error code (<0 means failure)
**************************************************************************/
int vnloo_reset_tmparr( Vnloo* vnl)
{
	if ( vnl == NULL ) return(-1);

	memset(vnl->missinput, '\0', vnl->nx*vnl->ny*sizeof(float)); 
	memset(vnl->miss, '\0', vnl->nx*vnl->ny*sizeof(float)); 
	memset(vnl->fullinput, '\0', vnl->nx*vnl->ny*sizeof(float)); 
	memset(vnl->full, '\0', vnl->nx*vnl->ny*sizeof(float)); 
/*	for ( int i = 0; i < (int) (vnl->nx*vnl->ny); i++) {
		vnl->miss[i] = 1.;
		vnl->missinput[i] = 1.;
		vnl->full[i] = 1.;
		vnl->fullinput[i] = 1.;
	}		
*/
	return(0);
}

/************************************************************************
@Function: vnloo_put
@Description:
	store values into temp arrays.
@Algorithm:
@Arguments:
	Vnloo*  vnl		Vnloo structure
	int		x, y 	array coordinates
	float	Armdp	miss-input cross product
	float	Afmdp	full-input cross product
	float	A2m		miss squared amplitude
	float	A2f		full squared amplitude
@Returns:
	int		error code (<0 means failure)
**************************************************************************/
int vnloo_put(Vnloo* vnl, int x, int y, float Armdp, float Arfdp, float A2m, float A2f)
{
	if ( vnl == NULL ) return(-1);

	unsigned long i = vnl->nx*(y/vnl->bin) + x/vnl->bin;
	vnl->missinput[i] += Armdp; 
	vnl->miss[i] += A2m; 
	vnl->fullinput[i] += Arfdp;
	vnl->full[i] += A2f;

	return(0);
}

/************************************************************************
@Function: vnloo_add_slice
@Description:
	add values from slice to intermediate volumes.
@Algorithm:
	an array containing the coords corresponding to the position of the points in the
	slice with respect to the volume is calculated, then the volumes are loaded from
	the files, slice by slice, and the calculated values are added.  
@Arguments:
	Vnloo*  vnl		Vnloo structure
	Bimage* b		image reference
@Returns:
	int		error code (<0 means failure)
**************************************************************************/
int vnloo_add_slice(Vnloo* vnl, Bimage* b)
{
	if ( vnl == NULL || b == NULL) return(-1);

	Proj* p = NULL;
	Bimage *bf = NULL, *bfi = NULL, *bm =NULL, *bmi = NULL;
	float *fdata = NULL, *fidata = NULL, *mdata =NULL, *midata = NULL;
	int		i, ix, iy, iz;
	unsigned long j;
	int 		x, y, z;
	Vector3	m, iv = {0.,0.,0.};
	View		v = view_from_4_values((double) b->image[0].vx, (double) b->image[0].vy, (double) b->image[0].vz, (double) b->image[0].angle);
	Matrix3	mat = matrix3_from_view(v);
	int		hx = (int) ((vnl->nx - 1)/2), hy = (int) ((vnl->ny - 1)/2);	
	int 		hz = (int) ((vnl->nz - 1)/2);
	int		sx = (int) (vnl->nx - 1 - hx);
	int		sy = (int) (vnl->ny - 1 - hy);
	int		sz = (int) (vnl->nz - 1 - hz);
	
	// determine transformation of coords from the slice to the volume
	for ( y=i=0; y< (int) vnl->ny; y++ ) {
		iv.y = y;
		if ( y > hy ) iv.y -= vnl->ny;
		for ( x=0; x< (int) vnl->nx; x++, i++ ) {
			iv.x = x;
			if ( x > hx ) iv.x -= vnl->nx;
			m = vector3_matrix3_multiply(mat, iv);

			iz = vround(m.z/vnl->zb);      // Nearest neighbour
//			if ( iz < 0 ) iz += vnl->nz;
			iz += sz;

			iy = vround(m.y);
//			if ( iy < 0 ) iy += vnl->ny;
			iy += sy;
			if (iy<0 || iy>=(int)vnl->ny) iz = (int) (2*vnl->nz); // set a false value in the z
			ix = vround(m.x);
//			if ( ix < 0 ) ix += vnl->nx;
			ix += sx;
			if (ix<0 || ix>=(int)vnl->nx) iz = (int) (2*vnl->nz); // set a false value in the z
			vnl->vcrd[i] = vectorint3_from_3_values(ix, iy, iz);
		}
	}

	p = proj_init(1);
	proj_cpy_img_sett(p, b);
	p->x = vnl->nx;
	p->y = vnl->ny;
	proj_to_datatype(p, Float);

	int z_min = (int) vnl->nz;
	int z_max = 0;
	for (i=0; i< (int) (vnl->ny*vnl->nx)-1; i++) {
		if(vnl->vcrd[i].z < z_min && vnl->vcrd[i].z >=0 ) z_min=vnl->vcrd[i].z;
		if(vnl->vcrd[i].z > z_max && vnl->vcrd[i].z <(int)vnl->nz) z_max=vnl->vcrd[i].z;
	}		 
	// load z-slices, assign values and write back to file
	for ( z = z_min; z <= z_max; z++) {
		bf = io_read_image(vnl->ffull, 1, z);
		bfi = io_read_image(vnl->ffullinput, 1, z);
		bm = io_read_image(vnl->fmiss, 1, z);
		bmi = io_read_image(vnl->fmissinput, 1, z);
		fdata = (float*) bf->data;
		fidata = (float*) bfi->data;
		mdata = (float*) bm->data;
		midata = (float*) bmi->data;
		for ( y=i=0; y < (int) vnl->ny; y++ ) {
			for ( x=0; x < (int) vnl->nx; x++, i++ ) {
				if ( vnl->vcrd[i].z == (int) z ) {
					ix = vnl->vcrd[i].x;
					iy = vnl->vcrd[i].y;
					j = iy*vnl->nx+ix;
					fdata[j] += vnl->full[i];
					fidata[j] += vnl->fullinput[i];
					mdata[j] += vnl->miss[i];
					midata[j] += vnl->missinput[i];
				}
			}
		}
		proj_putimage(p, bf, 0, z);
		io_pwrite(vnl->ffull, p, z, 1);
		proj_deleteimage(p, 0);

		proj_putimage(p, bfi, 0, z);
		io_pwrite(vnl->ffullinput, p, z, 1);
		proj_deleteimage(p, 0);

		proj_putimage(p, bm, 0, z);
		io_pwrite(vnl->fmiss, p, z, 1);
		proj_deleteimage(p, 0);

		proj_putimage(p, bmi, 0, z);
		io_pwrite(vnl->fmissinput, p, z, 1);
		proj_deleteimage(p, 0);
	}

	proj_kill(p);

	return(0);
}

/************************************************************************
@Function: vnloo_compute
@Description:
	evaluate vloumetric resolution from auxiliary volumes.
@Algorithm:
	The auxiliary volumes are loaded from 	the files, slice by slice,
	the volumetric nloo is evaluated and then stored into the final file.  
@Arguments:
	Vnloo*  vnl		Vnloo structure
@Returns:
	int		error code (<0 means failure)
**************************************************************************/
int vnloo_compute(Vnloo* vnl)
{
	if ( vnl == NULL ) return(-1);

	Bimage *bv = NULL, *bf = NULL, *bfi = NULL, *bm =NULL, *bmi = NULL;
	float *resdata = NULL, *fdata = NULL, *fidata = NULL, *mdata =NULL, *midata = NULL;

	bv = io_read_image(vnl->ffull,0,-1);
	Proj* p = proj_init(1);
	proj_cpy_img_sett(p, bv);
	kill_img(bv);


	// load z-slices, evaluate resolution and write result to file
	for ( int z = 0; z < (int) vnl->nz; z++) {
		bv = io_read_image(vnl->fvres, 1, z);
		bf = io_read_image(vnl->ffull, 1, z);
		bfi = io_read_image(vnl->ffullinput, 1, z);
		bm = io_read_image(vnl->fmiss, 1, z);
		bmi = io_read_image(vnl->fmissinput, 1, z);
		fdata = (float*) bf->data;
		fidata = (float*) bfi->data;
		mdata = (float*) bm->data;
		midata = (float*) bmi->data;
		resdata = (float*) bv->data;
		for ( int i=0; i< (int) (vnl->nx*vnl->ny); i++ ) {
			resdata[i] = 0.;
			if (fabs(midata[i]) < SMALLFLOAT || fabs(fidata[i]) < SMALLFLOAT) resdata[i] = 0.;
			else if (sqrt(mdata[i]) > SMALLFLOAT && fabs(fdata[i]) > SMALLFLOAT)
				resdata[i] = midata[i]/(sqrt(mdata[i])) / (fidata[i]/(sqrt(fdata[i])));
			if ( resdata[i] > 1. ) resdata[i] = 1.;
			if ( resdata[i] < -1. ) resdata[i] = -1.;
		}	
		proj_putimage(p, bv, 0, z);
		io_pwrite(vnl->fvres, p, z, 1);
		proj_deleteimage(p, 0);
		kill_img(bf);
		kill_img(bfi);
		kill_img(bm);
		kill_img(bmi);
	}

	io_update_stats(vnl->fvres);
	
	return(0);
}


//**** local auxiliary functions ***
int vround (double value)
{
	if (value < 0)
	//or as an alternate use:
	//return (int) -(floor(-value + 0.5));
	 	return (int) ceil ( value - 0.5);
	else
		return (int) floor( value + 0.5);
}
