/*
	imd_tilt.cc
	Functions to launch reconstruction program by imod - tilt
	Author: Giovanni Cardone
	Created: 2004	Modified: 20050628
*/

#include "imod.h"

// global variables
const char* ITILT_DEFAULT = "tilt.com";

/************************************************************************
@Function: imd_tilt_read
@Description:
	Init from file the parameters for tilt program in imod.
@Algorithm:
	If filename not given, look for tilt.com in the current directory.
	If one of the files exists read the parameters
@Arguments:
	char* parsname		parameters file name.
@Returns:
	Imod_tilt*			imod-tilt parameters
*************************************************************************/
Imod_tilt* imd_tilt_read(char* parsname)
{
	Imod_tilt* im_tilt = (Imod_tilt *) balloc(sizeof(Imod_tilt));

	imd_tilt_setdefault(im_tilt);

	if ( parsname == NULL) {
		strncpy(im_tilt->parsname, ITILT_DEFAULT, strlen(ITILT_DEFAULT));
	} else {
		strncpy(im_tilt->parsname, parsname, strlen(parsname));
	}

	if ( ! imd_tilt_parfile_exist(im_tilt) ) {
		fprintf(stderr,"Error: file %s does not exist\n", im_tilt->parsname);
		return (NULL);
	}

	if (verbose & VERB_PROCESS)
		printf("Opening file %s for reading parameters\n", im_tilt->parsname);

	if (imd_tilt_readpars(im_tilt) < 0) {
		fprintf(stderr,"Error while reading imod paramters from file %s\n", im_tilt->parsname);
		return (NULL);
	}

	return(im_tilt);
}

/************************************************************************
@Function: imd_tilt_parfile_exist
@Description:
	Check if given tilt paramter file already present in the directory.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	int				>0 if file exists
*************************************************************************/
int imd_tilt_parfile_exist(Imod_tilt * im_tilt)
{
	FILE * fd;
	if ( ( fd = fopen(im_tilt->parsname, "r" ) ) == NULL) {
		return(0);
	} else {
		fclose(fd);
		return(1);
	}

}

/************************************************************************
@Function: imd_tilt_checkvol
@Description:
	Check consistency between parameters and volume size.
	Update PARALLEL/PERPENDICULAR flag according to the result.
	Exit with error if sizes not consistent.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	int				<0 if error
*************************************************************************/
int imd_tilt_checkvol(Imod_tilt * im_tilt)
{
	if ( im_tilt == NULL) return (-1);

	char* out_vol = catenate2strings(imd_tilt_getvolfile(im_tilt),":mrc");

	Bimage* bvol = io_read_image(out_vol, 0, -1);

	if ( im_tilt->nx==bvol->x && im_tilt->ny==bvol->y && im_tilt->thickness==bvol->z) {
		// parallel mode
		if( verbose && VERB_PROCESS && im_tilt->parallel!=IPARALLEL)
			printf("Wrong parameter in %s: changing from PERPENDICULAR to PARALLEL mode\n",im_tilt->parsname);
		im_tilt->parallel = IPARALLEL;
	} else if ( im_tilt->nx==bvol->x && im_tilt->ny==bvol->z && im_tilt->thickness==bvol->y) {
		// perpendicular mode
		if( verbose && VERB_PROCESS && im_tilt->parallel==IPARALLEL)
			printf("Wrong parameter in %s: changing from PARALLEL to PERPENDICUALR mode\n",im_tilt->parsname);
		im_tilt->parallel = !IPARALLEL;
	} else {
		printf("ERROR(imd_tilt_checkvol): parameters in %s and volume size of %s not consistent!",
			im_tilt->parsname,im_tilt->out_vol);
			exit(1);
	}

	bfree_string(out_vol);

	kill_img(bvol);

	return(0);
}

/************************************************************************
@Function: imd_tilt_getlog
@Description:
	Get the logarithm setting.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	int				logarithm flag (rad)
*************************************************************************/
int imd_tilt_getlog(Imod_tilt * im_tilt)
{
	return(im_tilt->setlog);

}

/************************************************************************
@Function: imd_tilt_getlogbase
@Description:
	Get the base value to add before the log.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	float			base value
*************************************************************************/
float imd_tilt_getlogbase(Imod_tilt * im_tilt)
{
	return(im_tilt->logbase);
}

/************************************************************************
@Function: imd_tilt_getxtilt
@Description:
	Get the tilt around x axis.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	float			x tilt angle (rad)
*************************************************************************/
float imd_tilt_getxtilt(Imod_tilt * im_tilt)
{
	return(im_tilt->xaxistilt*M_PI/180.0);
}

/************************************************************************
@Function: imd_tilt_getresolutionlimit
@Description:
	Get the resolution limit used for the reconstruction.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	float			resolution limit (Nyquist fraction, between 0 and 0.5)
*************************************************************************/
float imd_tilt_getresolutionlimit(Imod_tilt * im_tilt)
{
	return(im_tilt->radial_max+0.5*im_tilt->radial_fall);
}

/************************************************************************
@Function: imd_tilt_gettiltfile
@Description:
	Get the name of the tilt file.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	char*			filename
*************************************************************************/
char* imd_tilt_gettiltfile(Imod_tilt * im_tilt)
{
	char *  input_string = im_tilt->tiltfile;

	char *	string = (char *) balloc((strlen(input_string)+1)*sizeof(char));

	memcpy(string, input_string, strlen(input_string));

	return(string);

}

/************************************************************************
@Function: imd_tilt_getprojfile
@Description:
	Get the name of the input projections file.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	char*			filename
*************************************************************************/
char* imd_tilt_getprojfile(Imod_tilt * im_tilt)
{
	char *  input_string = im_tilt->in_proj;

	char *	string = (char *) balloc((strlen(input_string)+1)*sizeof(char));

	memcpy(string, input_string, strlen(input_string));

	return(string);

}

/************************************************************************
@Function: imd_tilt_getvolfile
@Description:
	Get the name of the output volume file.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	char*			filename
*************************************************************************/
char* imd_tilt_getvolfile(Imod_tilt * im_tilt)
{
	char *  input_string = im_tilt->out_vol;

	char *	string = (char *) balloc((strlen(input_string)+1)*sizeof(char));

	memcpy(string, input_string, strlen(input_string));

	return(string);

}

/************************************************************************
@Function: imd_tilt_putvolfile
@Description:
	Assign the name of the output volume file.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
	char*			filename
@Returns:
	int 				<0 if error
*************************************************************************/
int imd_tilt_putvolfile(Imod_tilt * im_tilt, char* filename)
{

	strcpy(im_tilt->out_vol, filename);

	return(0);

}

/************************************************************************
@Function: imd_tilt_readpars
@Description:
	Read parameters for tilt code from file.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	int				<0 if data not read
*************************************************************************/
int imd_tilt_readpars(Imod_tilt * im_tilt)
{
	const int NCHARS = 250;
	char  fbuf[NCHARS];
	char s1[NCHARS];
	char s2[NCHARS];
	FILE * fd;
	if ( ( fd = fopen(im_tilt->parsname, "r" ) ) == NULL) {
		error_show(im_tilt->parsname, __FILE__, __LINE__);
		return(-1);
	}

	im_tilt->setlog = 0;

//	while ( strcmp("DONE",s1) ) {

//		fgets(fbuf, NCHARS, fd);
	while ( fgets(fbuf, NCHARS, fd), !feof(fd) ) {

		if ( fbuf[0] !='#') {
			sscanf(fbuf,"%s\n",s1);
			if ( !strcmp(s1,"$tilt") ) {
// read input projections file name
				fgets(fbuf, NCHARS, fd);
				sscanf(fbuf,"%s\n",s2);
				if ( !strcmp(s2,"InputProjections") )
					sscanf(fbuf,"%s %s\n",s2,im_tilt->in_proj);
				else
					sscanf(fbuf,"%s\n",im_tilt->in_proj);
// read output tomogram file name
				fgets(fbuf, NCHARS, fd);
				sscanf(fbuf,"%s\n",s2);
				if ( !strcmp(s2,"OutputFile") )
					sscanf(fbuf,"%s %s\n",s2,im_tilt->out_vol);
				else
					sscanf(fbuf,"%s\n",im_tilt->out_vol);
//				fscanf(fd,"%s %s\n",s2,im_tilt->in_proj);
//				fscanf(fd,"%s %s\n",s2,im_tilt->out_vol);
			}
			else if ( !strcmp(s1,"FULLIMAGE") )
				sscanf(fbuf,"%s %d %d\n", s1, &im_tilt->nx, &im_tilt->ny);
			else if ( !strcmp(s1,"RADIAL") )
				sscanf(fbuf,"%s %f %f\n", s1, &im_tilt->radial_max, &im_tilt->radial_fall);
			else if ( !strcmp(s1,"SCALE") ) {
				im_tilt->setscale = 1;
				sscanf(fbuf,"%s %f %f\n", s1, &im_tilt->scale_flevl, &im_tilt->scale_scale);
			} else if ( !strcmp(s1,"SHIFT") ) {
				im_tilt->setshift = 1;
				sscanf(fbuf,"%s %d %d\n", s1, &im_tilt->shift_x, &im_tilt->shift_z);
			} else if ( !strcmp(s1,"MODE") )
				sscanf(fbuf,"%s %d\n", s1, &im_tilt->mode);
			else if ( !strcmp(s1,"THICKNESS") )
				sscanf(fbuf,"%s %d\n", s1, &im_tilt->thickness);
			else if ( !strcmp(s1,"TILTFILE") )
				sscanf(fbuf,"%s %s\n", s1, im_tilt->tiltfile);
			else if ( !strcmp(s1,"FBPINTERP") ) {
				im_tilt->setfbpinterp = 1;
				sscanf(fbuf,"%s %d\n", s1, &im_tilt->fbpinterp);
			} else if ( !strcmp(s1,"PARALLEL") )
				im_tilt->parallel = IPARALLEL;
			else if ( !strcmp(s1,"PERPENDICULAR") )
				im_tilt->parallel = !IPARALLEL;
			else if ( !strcmp(s1,"LOG") ) {
				im_tilt->setlog = 1;
				sscanf(fbuf,"%s %f\n", s1, &im_tilt->logbase);
			} else if ( !strcmp(s1,"XAXISTILT") )
				sscanf(fbuf,"%s %f\n", s1, &im_tilt->xaxistilt);
			else if ( !strcmp(s1,"EXCLUDELIST") ) {
				unsigned long i = strlen("EXCLUDELIST");
				im_tilt->nxl = str_get_list(im_tilt->exclude_list, XCLDN, &fbuf[i]);
				if ( verbose & VERB_DEBUG ) {
					printf("Tilt images excluded from reconstruction (%d):\n",im_tilt->nxl);
					for ( int it = 0; it < im_tilt->nxl; it++)
						printf("%d %d\n",it+1,im_tilt->exclude_list[it]);
				}
			}
			else if ( !strcmp(s1,"UseGPU") ) {
				sscanf(fbuf,"%s %d\n", s1, &im_tilt->use_gpu);
			}
		}
	}

	fclose(fd);

	return(0);
}

/************************************************************************
@Function: imd_tilt_setdefault
@Description:
	Set parameters for tilt code to default values.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	int				<0 if data not set
*************************************************************************/
int imd_tilt_setdefault(Imod_tilt * im_tilt)
{
	if ( im_tilt == NULL) return (-1);

	strcpy(im_tilt->parsname,ITILT_DEFAULT);
	im_tilt->nx = 0;
	im_tilt->ny = 0;
	im_tilt->setlog = 0;
	im_tilt->logbase = ILOGBASE;
	im_tilt->parallel = IPARALLEL;
	im_tilt->radial_max = IRADIAL_MAX;
	im_tilt->radial_fall = IRADIAL_FALL;
	im_tilt->mode = IMODE;
	im_tilt->setfbpinterp = 0;
	im_tilt->fbpinterp = 1;
	im_tilt->setshift = 0;
	im_tilt->shift_x = 0;
	im_tilt->shift_z = 0;
	im_tilt->setscale = 0;
	im_tilt->scale_flevl = ISCALE_FLEVL;
	im_tilt->scale_scale = ISCALE_SCALE;
	im_tilt->sub_nx = ISUB_NX;
	im_tilt->sub_ny = ISUB_NY;
	im_tilt->nxl = 0;
	im_tilt->xaxistilt = IXAXISTILT;
	im_tilt->use_gpu = -1;

	return(0);
}

/************************************************************************
@Function: imd_tilt_setparallel
@Description:
	Set reconstruction to PARALLEL mode
	(produces sections parallel to the plane of the zero tilt projection)
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	int				<0 if data not set
*************************************************************************/
int imd_tilt_setparallel(Imod_tilt * im_tilt)
{
	if ( im_tilt == NULL) return (-1);

	im_tilt->parallel = IPARALLEL;

	return(0);
}

/************************************************************************
@Function: imd_tilt_isparallel
@Description:
	Check if reconstruction set to PARALLEL mode
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	int				1 if parallel flag set
*************************************************************************/
int imd_tilt_isparallel(Imod_tilt * im_tilt)
{
	if ( im_tilt == NULL) return (-1);

	return(im_tilt->parallel == IPARALLEL);
}

/************************************************************************
@Function: imd_tilt_putpars
@Description:
	Assign given values to tilt code parameters.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
	char * in_file  input projections file
	char * out_file output reconstruction file
	int nx			x size
	int ny			y size
	int nz			thickness
	int exp_tag			exp relationship between intensity and density
	float xtilt		tilt angle around x
	char* tilt_file tilt angles file name
@Returns:
	int				<0 if data not set
*************************************************************************/
int imd_tilt_putpars(Imod_tilt * im_tilt, char * in_file, char * out_file, int nx, int ny, int nz, int exp_tag, float xtilt, char* tilt_file)
{
	if (im_tilt == NULL) return(-1);

	if ( in_file == NULL) {
		printf("ERROR(im_tilt_putpars): file with input projections needed!");
		return(-1);
	}
	else
		strcpy(im_tilt->in_proj,in_file);
	if ( out_file == NULL) {
		printf("ERROR(im_tilt_putpars): output reconstruction file needed!");
		return(-1);
	}
	else
		strcpy(im_tilt->out_vol,out_file);
	im_tilt->nx = nx;
	im_tilt->ny = ny;
	im_tilt->setlog = exp_tag;
	if ( nz>0 ) im_tilt->thickness = nz;
	else if ( nx > ny ) im_tilt->thickness = nx;
	else im_tilt->thickness = ny;
	im_tilt->xaxistilt = xtilt;
	if ( tilt_file == NULL) {
		printf("ERROR(im_tilt_putpars): file with tilt angles needed!");
		return(-1);
	}
	else
		strcpy(im_tilt->tiltfile,tilt_file);


	return(0);
}

/************************************************************************
@Function: imd_tilt_cpy_and_modify
@Description:
	Copy a tilt parameter file, modifying its values, when provided.
@Algorithm:
@Arguments:
	Imod_tilt* im_tilt_dst	destinary mod-tilt parameters
	Imod_tilt* im_tilt_src	source mod-tilt parameters
	char * parsname			input projections file
	char * in_file			input projections file
	char * out_file			output reconstruction file
	int nx					x size
	int ny					y size
	int nz					thickness
	float xtilt				tilt angle around x
	char* tilt_file			tilt angles file name
@Returns:
	int				<0 if copy unsuccessfull
*************************************************************************/
int imd_tilt_cpy_and_modify(Imod_tilt * im_tilt_dst, Imod_tilt* im_tilt_src, char * parsname,
	char * in_file, char * out_file, int nx, int ny, int nz, float xtilt, char* tilt_file)
{
	if (im_tilt_dst == NULL || im_tilt_src == NULL) return(-1);

	imd_tilt_cpy(im_tilt_dst, im_tilt_src);

	if ( parsname != NULL)
		strcpy(im_tilt_dst->parsname,parsname);

	if ( in_file != NULL) strcpy(im_tilt_dst->in_proj,in_file);

	if ( out_file != NULL) strcpy(im_tilt_dst->out_vol,out_file);

	if ( nx > 0) im_tilt_dst->nx = nx;

	if ( ny > 0) im_tilt_dst->ny = ny;

	if ( nz > 0) im_tilt_dst->thickness = nz;

	if ( xtilt<90. && xtilt>-90. ) im_tilt_dst->xaxistilt = xtilt;

	if ( tilt_file != NULL)	strcpy(im_tilt_dst->tiltfile,tilt_file);

	return(0);
}

/************************************************************************
@Function: imd_tilt_cpy
@Description:
	Copy a tilt parameter file.
@Algorithm:
@Arguments:
	Imod_tilt* im_tilt_dst	destinary mod-tilt parameters
	Imod_tilt* im_tilt_src	source mod-tilt parameters
@Returns:
	int				<0 if copy unsuccessfull
*************************************************************************/
int imd_tilt_cpy(Imod_tilt * im_tilt_dst, Imod_tilt* im_tilt_src)
{

	strcpy(im_tilt_dst->parsname,im_tilt_src->parsname);
	strcpy(im_tilt_dst->in_proj,im_tilt_src->in_proj);
	strcpy(im_tilt_dst->out_vol,im_tilt_src->out_vol);
	im_tilt_dst->nx = im_tilt_src->nx;
	im_tilt_dst->ny = im_tilt_src->ny;
	im_tilt_dst->thickness = im_tilt_src->thickness;
	im_tilt_dst->setlog = im_tilt_src->setlog;
	im_tilt_dst->logbase = im_tilt_src->logbase;
	im_tilt_dst->parallel = im_tilt_src->parallel;
	im_tilt_dst->radial_max = im_tilt_src->radial_max;
	im_tilt_dst->radial_fall = im_tilt_src->radial_fall;
	im_tilt_dst->mode = im_tilt_src->mode;
	im_tilt_dst->setfbpinterp = im_tilt_src->setfbpinterp;
	im_tilt_dst->fbpinterp = im_tilt_src->fbpinterp;
	im_tilt_dst->setscale = im_tilt_src->setscale;
	im_tilt_dst->scale_flevl = im_tilt_src->scale_flevl;
	im_tilt_dst->scale_scale = im_tilt_src->scale_scale;
	im_tilt_dst->setshift = im_tilt_src->setshift;
	im_tilt_dst->shift_x = im_tilt_src->shift_x;
	im_tilt_dst->shift_z = im_tilt_src->shift_z;
	im_tilt_dst->sub_nx = im_tilt_src->sub_nx;
	im_tilt_dst->sub_ny = im_tilt_src->sub_ny;
	im_tilt_dst->nxl = im_tilt_src->nxl;
	im_tilt_dst->use_gpu = im_tilt_src->use_gpu;
	for (int i = 0; i<im_tilt_src->nxl; i++ )
		im_tilt_dst->exclude_list[i] = im_tilt_src->exclude_list[i];
	im_tilt_dst->xaxistilt = im_tilt_src->xaxistilt;
	strcpy(im_tilt_dst->tiltfile,im_tilt_src->tiltfile);


	return(0);
}

/************************************************************************
@Function: imd_tilt_writepars
@Description:
	Write parameters for tilt code to file.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	int				<0 if data not written
*************************************************************************/
int imd_tilt_writepars(Imod_tilt * im_tilt)
{
	FILE * fd;
	if ( ( fd = fopen(im_tilt->parsname, "w" ) ) == NULL) return(-1);

	fprintf(fd,"# tilt command file automatically generated by the Electra package\n");
	fprintf(fd,"#     NIAMS/NIH Laboratory of Structural Biology Research\n");
	fprintf(fd,"#\n");
	fprintf(fd,"$tilt\n");
	fprintf(fd,"%s\n",im_tilt->in_proj);
	fprintf(fd,"%s\n",im_tilt->out_vol);
	fprintf(fd,"FULLIMAGE %d %d\n",im_tilt->nx,im_tilt->ny);
	if (im_tilt->setlog)
		fprintf(fd,"LOG %4.2f\n",im_tilt->logbase);
	fprintf(fd,"MODE %d\n",im_tilt->mode);
	if (im_tilt->parallel) fprintf(fd,"PARALLEL\n");
	else fprintf(fd,"PERPENDICULAR\n");
	fprintf(fd,"RADIAL %7.3f %7.3f\n",im_tilt->radial_max,im_tilt->radial_fall);
	if (im_tilt->setscale)
		fprintf(fd,"SCALE %7.3f %7.3f\n",im_tilt->scale_flevl,im_tilt->scale_scale);
	if (im_tilt->setshift)
		fprintf(fd,"SHIFT %d %d\n",im_tilt->shift_x,im_tilt->shift_z);
	if (im_tilt->setfbpinterp)
		fprintf(fd,"FBPINTERP %d\n",im_tilt->fbpinterp);
	fprintf(fd,"SUBSETSTART %d %d\n",im_tilt->sub_nx,im_tilt->sub_ny);
	fprintf(fd,"THICKNESS %d\n",im_tilt->thickness);
	fprintf(fd,"XAXISTILT %7.3f\n",im_tilt->xaxistilt);
	fprintf(fd,"TILTFILE %s\n",im_tilt->tiltfile);
	if ( im_tilt->nxl > 0 ) {
		fprintf(fd,"EXCLUDE");
		for ( int i = 0; i < im_tilt->nxl; i++ )
			fprintf(fd," %d", im_tilt->exclude_list[i]);
		fprintf(fd,"\n");
	}
	if (im_tilt->use_gpu >= 0)
		fprintf(fd,"UseGPU %d\n",im_tilt->use_gpu);
	fprintf(fd,"DONE\n");

	fclose(fd);

	return(0);
}

/************************************************************************
@Function: imd_tilt_execute
@Description:
	Execute tilt program by IMOD.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
@Returns:
	int				<0 if execution failed
*************************************************************************/
int imd_tilt_execute(Imod_tilt * im_tilt) {

	if ( im_tilt == NULL) return (-1);

	Bimage* b;
	int isok = 0;

	while (!isok) {
//  char * exec_string = catenate2strings(copystring("submfg "),im_tilt->parsname);
//	if ( system(exec_string) <0 ) return(-1);
		if ( ext_execute("submfg",im_tilt->parsname) <0 ) return(-1);

		b = io_read_image(im_tilt->out_vol,0,-1);

		if ( b->x == 0 )
			imd_tilt_rmalloutfiles(im_tilt);
		else isok = 1;

		kill_img(b);
	}

	return 0;
}

/************************************************************************
@Function: imd_tilt_angle_missing
@Description:
	Check if given angle omitted from reconstruction.
@Algorithm:
@Arguments:
	Imod_tilt*	im_tilt		imod-tilt parameters
	int			a			angle index
@Returns:
	int			1 if angle is omitted, 0 otherwise
*************************************************************************/
int imd_tilt_angle_missing(Imod_tilt * im_tilt, int a) {

	for ( int i = 0; i < im_tilt->nxl; i++ )
		if (im_tilt->exclude_list[i] == a ) return 1;

	return 0;
}

/************************************************************************
@Function: imd_tilt_add2xcldlist
@Description:
	Add an angle to exclude list
@Algorithm:
@Arguments:
	Imod_tilt*	im_tilt		imod-tilt parameters
	int			a			angle index
@Returns:
	int			<0 if operation failed
*************************************************************************/
int	imd_tilt_add2xcldlist( Imod_tilt* im_tilt, int a) {

	if ( im_tilt == NULL) return (-1);
	if (im_tilt->nxl >= XCLDN) return (-1);

	im_tilt->exclude_list[im_tilt->nxl] = a;
	im_tilt->nxl ++;

	return 0;
}

/************************************************************************
@Function: imd_tilt_rmxcldlist
@Description:
	Remove all the entries in the exclude list
@Algorithm:
@Arguments:
	Imod_tilt*	im_tilt		imod-tilt parameters
@Returns:
	int			<0 if operation failed
*************************************************************************/
int	imd_tilt_rmxcldlist( Imod_tilt* im_tilt) {

	if ( im_tilt == NULL) return (-1);

	im_tilt->nxl = 0;

	return 0;
}

/************************************************************************
@Function: imd_tilt_rmallfiles_butvol
@Description:
	delete all the input and generated files, except for the reconstruction
@Algorithm:
@Arguments:
	Imod_tilt*	im_tilt		imod-tilt parameters
@Returns:
	int			<0 if operation failed
*************************************************************************/
int	imd_tilt_rmallfiles_butvol( Imod_tilt* im_tilt) {

	int ierr = 0;

	ierr += remove(im_tilt->parsname);
	ierr += remove(filename_change_type(im_tilt->parsname,"log"));

	ierr += remove(im_tilt->tiltfile);
	ierr += remove(im_tilt->in_proj);

	if (ierr<0) fprintf(stderr,"Error while deleting files!");

	return(0);
}

/************************************************************************
@Function: imd_tilt_rmalloutfiles
@Description:
	Remove all the files possibly generated by the parameter file.
@Algorithm:
@Arguments:
	Imod_tilt*	im_tilt		imod-tilt parameters
@Returns:
	int			<0 if operation failed
*************************************************************************/
int	imd_tilt_rmalloutfiles( Imod_tilt* im_tilt) {

	int ierr = 0;

	ierr += remove(filename_change_type(im_tilt->parsname,"log"));

	ierr += remove(im_tilt->out_vol);

	if (ierr<0) fprintf(stderr,"Error while deleting files!");

	return(0);
}

/************************************************************************
@Function: imd_tilt_rmparfile
@Description:
	Delete the parameter file.
@Algorithm:
@Arguments:
	Imod_tilt*	im_tilt		imod-tilt parameters
@Returns:
	int			<0 if operation failed
*************************************************************************/
int	imd_tilt_rmparfile( Imod_tilt* im_tilt) {

	int ierr = 0;

	ierr += remove(im_tilt->parsname);

	if (ierr<0) fprintf(stderr,"Error while deleting files!");

	return(0);
}

/************************************************************************
@Function: imd_tilt_copyxcldlist
@Description:
	Copy the list of excluded projections from the imod parameter file.
@Algorithm:
@Arguments:
	int*			xcld_list	exclude list (output)
	int*			nxl			number of exluded projections (output)
	Imod_tilt*	im_tilt		imod-tilt parameters
@Returns:
	int			<0 if operation failed
*************************************************************************/
int	imd_tilt_copyxcldlist(int* xcld_list, int* nxl, Imod_tilt* im_tilt)
{

	*nxl = im_tilt->nxl;

	for (int i = 0; i < im_tilt->nxl; i++)
		xcld_list[i] = im_tilt->exclude_list[i];

	return(0);
}

/************************************************************************
@Function: imd_tilt_makelist
@Description:
	Initialize list of projections to be included in the analysis
	(i.e. to be excluded from the reconstruction)
@Algorithm:
@Arguments:
	int			it			iteration number
	int			np			total number of projections
	int			tstep		separation between projections to be concurrently evaluated
	Imod_tilt*	imd_pars		imod-tilt parameters
	int*			lp			list of projections (output)
@Returns:
	int						number of projections included
*************************************************************************/
int	imd_tilt_makelist( int it, int np, int tstep, Imod_tilt* imd_pars, int* lp) {

	int	nl = 0;
	int	tmod = it % tstep;
	for ( int tt=0; tt<np; tt++) {
		if ((tt % tstep == tmod) && !imd_tilt_angle_missing(imd_pars,tt+1))
			lp[nl++] = tt;
	}

	return(nl);
}

/************************************************************************
@Function: imd_tilt_putdatatype
@Description:
	Change the data format of the output volume file.
@Algorithm:
@Arguments:
	Imod_tilt*		imod-tilt parameters
	DataType			new data type
@Returns:
	int 				<0 if error
*************************************************************************/
int imd_tilt_putdatatype(Imod_tilt * im_tilt, DataType dtype)
{

	if ( dtype == UChar ) im_tilt->mode = 0;
	else if ( dtype == Short ) im_tilt->mode = 1;
	else if ( dtype == Float ) im_tilt->mode = 2;
	else {
		fprintf(stderr,"ERROR(imd_tilt_putdatatype): data format not supported!\n");
		exit(-1);
	}

	return(0);

}


