/* 
	ioroutines.cc
	Routines for reading/writing images from/to files 
	Author: Giovanni Cardone
	Created: 2004	Modified: 20050628
*/ 

#include "ioroutines.h"

/************************************************************************
@Function: io_read_image
@Description:
	General driver function to read multiple image formats
@Algorithm:
	This is the only image reading function that should be called
	from programs.
	If image selection not requested, it just calls the standard routine read_img.
	If requested, verifies if multi-images supported by the
	additional routines included in this package.
	At the moment the additional formats supported are:
	MRC
	The numbering goes from 0 to n-1.
@Arguments:
	char* filename		file name (plus any tags for the RAW format).
	int readdataflag	flag to activate reading of image data.
	int select_img		image selection in multi-image file (-1 = all images).
@Returns:
	Bimage*				the image structure, NULL if reading failed.
**************************************************************************/
Bimage* io_read_image(char* filename, int readdataflag, int select_img)
{

	char*		ext = NULL;	
	Bimage*  	p = NULL;

	if (select_img < 0) {
		p = read_img(filename, readdataflag, select_img);
		if(readdataflag) img_stats(p);
		return(p);
	}

	if ( filename && strlen(filename) ) ext = extension(filename);

	if ( !ext || strstr(ext,"raw") || strstr(ext,"mif") || strstr(ext,"tif") ||
		strstr(ext,"pif") || strstr(ext,"sf") || strstr(ext,"spi") ||
		strstr(ext,"img") || strstr(ext,"hed") || strstr(ext,"di") ) {
			p = read_img(filename, readdataflag, select_img);
			img_stats(p);
			return(p);
		}

	if ( ! strstr(ext,"mrc") ) {
		fprintf(stderr,"File format not supported by this specific program! Please convert into another format\n");
		error_show(filename, __FILE__, __LINE__);
		exit(-1);
	}
		
	if ( verbose & VERB_DEBUG )
		printf("DEBUG read_img: Filename = %s\n", filename);
	
	int 		err = 0;
	if ( filename && strlen(filename) ) ext = extension(filename);
	
	if ( !filename || !strlen(filename) ) {
		fprintf(stderr, "Error: No input image file given!\n");
		return(NULL);
	}
	
	p = init_img();
	p->dataflag = 0;
	if ( readdataflag ) p->dataflag = 1;
	filename_clean_copy(filename, p->filename, 256);

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG read_img: Transferred filename = %s\n", p->filename);
		printf("DEBUG read_img: Extension = %s\n", ext);
	}
	
	if ( strstr(ext,"mrc") )
		err = io_readMRC(p, select_img);
	else {
		fprintf(stderr, "Error: File format with extension \"%s\" not supported!\n", ext);
		err = -1;
	}
	
	bfree_string(ext);
		
	if ( err < 0 ) {
		error_show(filename, __FILE__, __LINE__);
		kill_img(p);
		return(NULL);
	}
	
	img_stats(p);
	img_check_param(p);
		
	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) ) {
		if ( select_img < 0 )
			printf("Reading file:                   %s\n", p->filename);
		else
			printf("Reading image %4d from file:   %s\n", select_img, p->filename);
	}
	
	if ( verbose & VERB_PROCESS ) {
		img_info(p);
		if ( p->dataflag < 1 ) printf("Data not read!\n\n");
	}
	
	if ( verbose & VERB_STATS )
		sub_img_info(p);
	
	return(p);
}

/************************************************************************
@Function: io_read_image_buff
@Description:
	Specific driver function to read data in buffer mode
@Algorithm:
	This function is to be used when handling large volumes. The data
	are buffered, meaning that only part of the data are read and stored.
	The Bimage can have been already opened, in this case just the content
	of the data array is update, according to the specifications of the
	Buff struct.
	If the file has been already opened once, the filename is no more needed.
	The size of the buffer and the portion of data stored are set in the
	Buff structure associated.
	Only one single image is assumed to be present in the file.
	At the moment the formats supported are:
	MRC
@Arguments:
	Bimage*				image structure (NULL if file never opened)
	char* filename		file name.
	int select_img		image selection: if -1, all images, if >= 0, one image.
						(not supported)
	Buff* bf				buffer info.
@Returns:
	Bimage*				the image structure, NULL if reading failed.
**************************************************************************/
Bimage* io_read_image_buff(Bimage* p, char* filename, int select_img, Buff* bf)
{
	char*		ext = NULL;	

	char* fname = NULL;
	
	if (filename) fname = copystring(filename);
	else if (p && p->filename) fname = catenate2strings(copystring(p->filename),":mrc");
	else {
		fprintf(stderr,"ERROR io_read_image_buff: file name not given, and not present in the image info!\n");
		exit(-1);
	}
		
	if ( fname && strlen(fname) ) ext = extension(fname);
	if ( ! strstr(ext,"mrc") ) {
		fprintf(stderr,"File format not supported by this specific program! Please convert into mrc format\n");
		error_show(fname, __FILE__, __LINE__);
		exit(-1);
	}
		
	if ( verbose & VERB_DEBUG )
		printf("DEBUG io_read_image_buff: Filename = %s\n", fname);
	
	int 		err = 0;

	
	if (!p) {
		p = init_img();
		p->dataflag = 1;
		filename_clean_copy(fname, p->filename, 256);
		if ( verbose & VERB_DEBUG ) {
			printf("DEBUG io_read_image_buff: Transferred filename = %s\n", p->filename);
			printf("DEBUG io_read_image_buff: Extension = %s\n", ext);
		}
	} else if ((p && !bf ) || (p && bf && bf->last==0)) {
		fprintf(stderr,"ERROR io_read_image_buff: image and buffer settings not correct!\n");
		error_show(fname, __FILE__, __LINE__);
		exit(-1);
	}
		
	
	if ( strstr(ext,"mrc") )
		err = io_readMRC_buff(p, bf);
	else {
		fprintf(stderr, "Error: File format with extension \"%s\" not supported!\n", ext);
		err = -1;
	}
	
	bfree_string(ext);
	bfree_string(fname);
			
	if ( err < 0 ) {
		error_show(fname, __FILE__, __LINE__);
		kill_img(p);
		return(NULL);
	}
	
/*	GCA: stats and check can not be done on buffered images
 	img_stats(p);
	img_check_param(p);
*/
	
	return(p);
}

/************************************************************************
@Function: io_readMRC_header
@Description:
	Reading the header of a MRC map image file format.
@Algorithm:
	A 3D image format used in electron microscopy.
	Header size:				1024 bytes followed by the symmetry 
		operator table which is composed of 80 character lines, each 
		line for a symmetry operator.
	File format extensions:  	.mrc
	The identifier is a 4-byte machine stamp (same as for CCP4 maps):
			1	Big-endian IEEE 		17 17 00 00
             2	VAX 					34 65 00 00
			3	Cray					-
             4	Little-endian IEEE	68 65 00 00
             5	Convex				85 17 00 00
			6	Fijitsu VP			-
				(Note: not always implemented - so unreliable)
	Byte order determination:	Data type and third dimension values
								must be less than 256*256.
	Data types: 				0 = byte, 1 = short, 2 = float,
							3 = complex short, 4 = complex float.
	Transform type: 			Centered hermitian
								The x-dimension contains the x-size
								of the full transform
@Arguments:
	FILE*	fimg			file pointer.
	Bimage*	p			the image structure.
	int 		select_img	image selection index.
	int*		swap			bbyte swap flag.
	int*     vax			vax machine flag.
@Returns:
	FILE*				r code (<0 means failure).
**************************************************************************/
 int 	io_readMRC_header(FILE *fimg, Bimage* p, int select_img, int *oswap, int *ovax)
 {
 	
	MRChead*	header = (MRChead *) balloc(sizeof(MRChead));
	if ( fread( header, MRCSIZE, 1, fimg ) < 1 ) return(-2);
	
    // Determine byte order and swap bytes if from little-endian machine
    char*   	b = (char *) header;
    int     	i, swap = 0, vax = 0;
    if ( ( abs( header->mode ) > SWAPTRIG ) || ( abs(header->nz) > SWAPTRIG ) ) {
    	if ( verbose & VERB_PROCESS )
	    	fprintf(stderr, "Warning: Swapping header byte order for 4-byte types\n");
    	swap = 1;
		int 	extent = MRCSIZE - 800; // exclude labels from swapping
    	for ( i=0; i<extent; i+=4 ) swapbytes(b+i, 4);
    }
    
	if ( verbose & VERB_DEBUG )
		printf("DEBUG rwMRC: Min, max: %g %g\n", header->amin, header->amax);

    // Convert VAX floating point types if necessary
    if ( header->amin > header->amax ) {
    	if ( verbose & VERB_PROCESS )
	    	fprintf(stderr, "Warning: Converting VAX floating point\n");
    	vax = 1;
		for ( i=40; i<64; i+=4 )		// Unit cell parameters
    	    vax2ieee(b+i, swap);
		for ( i=76; i<88; i+=4 )		// Min, max, mean
    	    vax2ieee(b+i, swap);
		for ( i=96; i<220; i+=4 )		// Extra and origin
    	    vax2ieee(b+i, swap);
    }
    
	if ( verbose & VERB_DEBUG )
		printf("DEBUG rwMRC: Any VAX conversion done\n");

	// Map the parameters
	p->x = header->nx;
	p->y = header->ny;
	if ( select_img<0 ) {
		p->z = header->nz;
	} else {
		p->z = 1;
		if ( select_img > header->nz - 1) select_img = header->nz - 1;
	}	
	p->n = 1;

	p->c = 1;
	p->colormodel = Gray;
	switch ( header->mode ) {
		case 0: p->datatype = UChar; break;
		case 1: p->datatype = Short; break;
		case 2: p->datatype = Float; break;
		case 3: p->datatype = ComplexShort; break;
		case 4: p->datatype = ComplexFloat; break;
		default: p->datatype = UChar; break;
	}
	p->offset = MRCSIZE + header->nsymbt;
	if ( header->mode > 2 && header->mode < 5) {
	fprintf(stderr,"MRC: mode = %d; p->x= %lu\n",header->mode,p->x);
		p->transform = CentHerm;
		fseeko(fimg, 0, SEEK_END);
		if ( ftello(fimg) > (off_t) p->offset + 0.8*p->x*p->y*p->z*gettypesize(p->datatype) )
			p->x = 2*(p->x - 1);
	fprintf(stderr,"MRC: mode = %d; p->x= %lu\n",header->mode,p->x);
		if ( header->mx%2 == 1 ) p->x += 1;	// Quick fix for odd x-size maps
	fprintf(stderr,"MRC: mode = %d; p->x= %lu\n",header->mode,p->x);
	}
	p->min = header->amin;
	p->max = header->amax;
	p->avg = header->amean;
	p->std = header->arms;
	p->ua = header->a;
	p->ub = header->b;
	p->uc = header->c;
	p->alf = header->alpha;
	p->bet = header->beta;
	p->gam = header->gamma;
	if ( header->mx ) p->ux = header->a/header->mx;
	if ( header->my ) p->uy = header->b/header->my;
	if ( header->mz ) p->uz = header->c/header->mz;
	p->spacegroup = header->ispg;
	memcpy(p->label, header->labels, 800);
	i = 800;
	while ( --i > 0 && ( p->label[i] == ' ' || p->label[i] == 0 ) )
		p->label[i] = 0;
	
	// Allocating the single sub-image and setting its origin
	p->image = (Bsub_image *) balloc(p->n*sizeof(Bsub_image));
	for ( i = 0; i < (int) p->n; i++ ) {
		p->image[i].ox = -header->nxStart;
		p->image[i].oy = -header->nyStart;
		p->image[i].oz = -header->nzStart;
		p->image[i].ox = header->xOrigin;
		p->image[i].oy = header->yOrigin;
		p->image[i].oz = header->zOrigin;
	}
	bfree(header, sizeof(MRChead));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG rwMRC: Header setup done\n");

	*oswap = swap;
	*ovax = vax;
	return(0);
}

/************************************************************************
@Function: io_readMRC
@Description:
	Reading a MRC map image file format.
@Algorithm:
	A 3D image format used in electron microscopy.
	Header size:				1024 bytes followed by the symmetry 
		operator table which is composed of 80 character lines, each 
		line for a symmetry operator.
	File format extensions:  	.mrc
	The identifier is a 4-byte machine stamp (same as for CCP4 maps):
				1	Big-endian IEEE 	17 17 00 00
             	2	VAX 				34 65 00 00
				3	Cray				-
             	4	Little-endian IEEE	68 65 00 00
             	5	Convex				85 17 00 00
				6	Fijitsu VP			-
				(Note: not always implemented - so unreliable)
	Byte order determination:	Data type and third dimension values
								must be less than 256*256.
	Data types: 				0 = byte, 1 = short, 2 = float,
								3 = complex short, 4 = complex float.
	Transform type: 			Centered hermitian
								The x-dimension contains the x-size
								of the full transform
@Arguments:
	Bimage*	p			the image structure.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 	io_readMRC(Bimage* p, int select_img)
{

    int     	swap = 0, vax = 0;

	FILE        *fimg;
    if ( ( fimg = fopen(p->filename, "r") ) == NULL ) return(-1);
 
    io_readMRC_header(fimg, p, select_img, &swap, &vax);

	p->data = io_read_data( fimg, p, select_img, swap, vax, 0 );
	
	fclose(fimg);
	
	return(0);
}

/************************************************************************
@Function: io_readMRC_buff
@Description:
	Reading a MRC map image file format in buffer mode.
@Algorithm:
	A 3D image format used in electron microscopy.
	Header size:				1024 bytes followed by the symmetry 
		operator table which is composed of 80 character lines, each 
		line for a symmetry operator.
	File format extensions:  	.mrc
	The identifier is a 4-byte machine stamp (same as for CCP4 maps):
				1	Big-endian IEEE 	17 17 00 00
             	2	VAX 				34 65 00 00
				3	Cray				-
             	4	Little-endian IEEE	68 65 00 00
             	5	Convex				85 17 00 00
				6	Fijitsu VP			-
				(Note: not always implemented - so unreliable)
	Byte order determination:	Data type and third dimension values
								must be less than 256*256.
	Data types: 				0 = byte, 1 = short, 2 = float,
								3 = complex short, 4 = complex float.
	Transform type: 			Centered hermitian
								The x-dimension contains the x-size
								of the full transform
@Arguments:
	Bimage*	p			the image structure.
	Buff*	bf			the buffer structure.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 	io_readMRC_buff(Bimage* p, Buff* bf)
{

 	Vector3 vol_origin;
    int     	swap = 0, vax = 0;

	FILE        *fimg;
    if ( ( fimg = fopen(p->filename, "r") ) == NULL ) return(-1);
 
 	if (bf && bf->last!=0)
 		vol_origin = vector3_from_3_values(p->image[0].ox,p->image[0].oy,p->image[0].oz);;

    io_readMRC_header(fimg, p, -1, &swap, &vax);

 	if (bf && bf->last!=0) {
		p->image[0].ox = vol_origin.x;
		p->image[0].oy = vol_origin.y;
		p->image[0].oz = vol_origin.z;
 	}

	p->data = buff_img_read_data(fimg, p, -1, swap, vax, 0, bf);
	
	fclose(fimg);
	
	return(0);
}

/************************************************************************
@Function: io_pwrite
@Description:
	General function to write projection data to file
@Algorithm:
	Depending on the value of the writeflag,
	a different operation is performed:
	 -3 = IO_HEADER = write/modify only the header
	 -2 = IO_ZEROS  = write the header and zeroed images
	 -1 = IO_ALL    = write header and images from the given projection
	  i >=0         = insert all the given projection images into a existent file,
	  				 from the position i indicated by writeflag 
	In IO_HEADER mode, if the file allready exists, only its header is modified.
	When the header-only mode IO_HEADER or the zero-values image writing mode
	IO_ZEROS are selected, then the number of slices is read from nzp
@Note:
	Only mrc format is supported
@Arguments:
	char* filename	   file name.
	Proj* p             projection structure
	int 	  writeflag     flag selecting what to write.
	int	  nzp		   number of projections
					   (it is read only for IO_HEADER and IO_ZEROS cases)
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int io_pwrite(char* filename, Proj* p, int writeflag, int nzp)
{
	if ( !p ) return(-1);
	
	if ( !filename || !strlen(filename) ) {
		fprintf(stderr, "Error: No output image file name given!\n");
		return(-1);
	}

	int 	err = 0;
	off_t ipos = writeflag;
	
	filename_clean_copy(filename, p->filename, 256);
	
	char*	ext = extension(filename);
	if ( !ext ) {
		fprintf(stderr, "Error(io_pwrite): No extension found in the file name \"%s\"!\n\n", filename);
		error_show(filename, __FILE__, __LINE__);
		return(-2);
	}
	
	if ( !strstr(ext,"mrc") ) {
		fprintf(stderr, "Warning(io_pwrite): Only MRC file format supported for writing!\n");
	}

	if ( verbose & VERB_DEBUG )
		printf("DEBUG io_pwrite: Image setup done for file %s %s (extension %s)\n", 
				filename, p->filename, ext);

// determine number of images to write
	int inp = 0;
	int pz = 0;
	switch ( writeflag ) {
		case IO_HEADER:
			inp = nzp;
			pz  = 0;
		break;
		case IO_ZEROS:
			inp = nzp;
			pz  = nzp;
		break;
		case IO_ALL:
			inp = (int) p->nimg;
			pz  = (int) p->nimg;
		break;
		default:
			inp = (int) p->nimg;
			pz  = (int) p->nimg;
		break;
	}

	MRChead*	header = NULL;
	
	if ( writeflag < 0 ) {
		header = (MRChead *) balloc(sizeof(MRChead));
	
		// Map the parameters
		strncpy(header->map, "MAP ", 4);
		set_CCP4_machine_stamp(header->machst);
		header->nx = (int) p->x;
		header->ny = (int) p->y;
		header->nz = inp;
		switch ( p->datatype ) {
			case UChar: header->mode = 0; break;
			case Short: header->mode = 1; break;
			case Float: header->mode = 2; break;
			case ComplexShort: header->mode = 3; break;
			case ComplexFloat: header->mode = 4; break;
			default: header->mode = 0; break;
		}
		header->nxStart = (int) (-p->ox - 0.5);
		header->nyStart = (int) (-p->oy - 0.5);
		header->nzStart = (int) (-p->oz - 0.5);
		header->mx = (int) p->x;
		header->my = (int) p->y;
		header->mz = inp;
		header->mapc = 1;
		header->mapr = 2;
		header->maps = 3;
		header->amin = p->min;
		header->amax = p->max;
		header->amean = p->avg;
		header->arms = p->std;
		header->a = p->ux*header->mx;
		header->b = p->uy*header->my;
		header->c = p->ux*header->mz;
//		header->c = 1.*header->mz;
		header->xOrigin = p->ox;
		header->yOrigin = p->oy;
		header->zOrigin = p->oz;
	
		header->alpha = 90.0;
		header->beta = 90.0;
		header->gamma = 90.0;
		header->ispg = 0;
	
		int			nsym = 0;

		header->nsymbt = nsym*80;

		header->nlabl = 10;
		if (p->img[0])	memcpy(header->labels, p->img[0]->label, 800);
	
	}  // end of header initialization

	long			   datatypesize = gettypesize(p->datatype);
	unsigned long   projsize = (unsigned long) p->x*p->y*pz*datatypesize;
	unsigned long   imgsize = (unsigned long) p->x*p->y*datatypesize;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG io_pwrite: Typesize = %ld,  Projection datasize = %ld\n", 
				datatypesize, projsize);

	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) )
		printf("Writing %d images to file:   %s\n", pz, p->filename);
	
    FILE      *fimg;
    char      *data_pnt = NULL;
    
 	if (writeflag < 0 && writeflag != IO_HEADER) {
		if ( ( fimg = fopen(p->filename, "w") ) == NULL ) {
			fprintf(stderr,"ERROR(io_pwrite): can not open file %s in 'w' mode\n",p->filename);
			return(-1);
		}
 	} else {
    		if ( ( fimg = fopen(p->filename, "r+") ) == NULL ) {
			fprintf(stderr,"ERROR(io_pwrite): can not open file %s in 'r+' mode\n",p->filename);
			return(-1);
		}
		rewind(fimg);
 	}
 	
	if (writeflag<0)
		fwrite(header, MRCSIZE, 1, fimg);
	else {
		fseeko(fimg, (off_t) MRCSIZE+ipos*imgsize, SEEK_SET);
	}
	
	if (writeflag==IO_ZEROS) {
		data_pnt = (char *) balloc(imgsize);		
	}

	if ( writeflag != IO_HEADER) {
	
		for (int i = 0, ii=0; i<inp; i++) {

			if ( writeflag != IO_ZEROS) {

				while ( p->img[ii] == NULL) ii++;

				img_RGB2gray(p->img[ii]);
	
			    switch ( p->datatype ) {
    					case SChar: img_to_byte((p->img[ii])); break;
    					case UShort: img_to_short((p->img[ii])); break;
    					case Int: img_to_float((p->img[ii])); break;
    					case ComplexInt: img_to_complex_float((p->img[ii])); break;
    					case Polar: img_polar2complex((p->img[ii])); break;
    					default: break;
    				}
				data_pnt = (p->img[ii])->data;
				ii++;
			}
			fwrite( data_pnt, imgsize, 1, fimg );
		}
	}
	fclose(fimg);
	
	if (writeflag == IO_ZEROS)
		bfree(data_pnt, imgsize);
		
	if (writeflag < 0) bfree(header, sizeof(MRChead));
	
	bfree_string(ext);
	
	if ( err < 0 ) {
		error_show(filename, __FILE__, __LINE__);
		return(err);
	}
	
//	if ( verbose & VERB_PROCESS )
//		proj_info(p);
	
 	// open file as image, eval stats, then rewrite only header
 	// GC: to be eliminated after a proj_stats function 
 	// calculating the overall average, standard deviation, minimum and
 	// maximum value be defined
 	if (writeflag == IO_ALL)
		io_update_stats(filename);
 	
	return(err);
}

/************************************************************************
@Function: io_update_stats
@Description:
	Recalculate stats of an image from a file
@Algorithm:
	Open the file, eval stats, then rewrite
@Arguments:
	char* filename	   file name.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int io_update_stats(char* filename)
{
 	int ierr;
 	
 	Bimage* b = io_read_image(filename, 1, -1);
 	ierr = img_stats(b);

	if (verbose & VERB_FULL) printf("Stats: %s mean %f max %f\n",filename,b->avg,b->max);
	Proj*	p = proj_init(1);
	proj_cpy_img_sett(p, b);
	proj_putimage(p, b, 0, 0);

	int np = (int) ((b->n==1 && b->z>1)?b->z:b->n); 
	
	ierr = io_pwrite(filename, p, IO_HEADER, np);

	proj_kill(p);
	
	return(ierr);
	
}
/************************************************************************
@Function: io_update_projfile_stats
@Description:
	Recalculate stats of an set of projections from a file
@Algorithm:
	Open the file, eval stats, then rewrite its header 
    Fom statistics are not calculated
@Arguments:
	char* filename	   file name.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int io_update_projfile_stats(char* filename)
{
 	int ierr;
 	
    unsigned char* 	udata = (unsigned char *) NULL;
    signed char* 	cdata = (signed char *) NULL;
    unsigned short* 	usdata = (unsigned short *) NULL;
    short* 	    		sdata = (short *) NULL;
    int* 	    		idata = (int *) NULL;
    float*  	    		fdata = (float *) NULL;
    complex_short* 	csdata = (complex_short *) NULL;
    complex_int* 	cidata = (complex_int *) NULL;
    complex_float*  	cfdata = (complex_float *) NULL;
    polar*  	    		pdata = (polar *) NULL;

 	Bimage* b = io_read_image(filename, 0, -1);
	int np = (b->n==1 && b->z>1)?b->z:b->n; 

	Bimage* p;
	char					bit;
    	unsigned long		i, j, n, datasize;
    	unsigned long		fomsize = (unsigned long) b->x*b->y*np;
	double				amp, tsum, tssum, avg, std;
	double*				min = (double *) balloc(np*sizeof(double));
	double*				max = (double *) balloc(np*sizeof(double));
	double*				sum = (double *) balloc(np*sizeof(double));
	double*				ssum = (double *) balloc(np*sizeof(double));

	for ( n=0; n<np; n++ ) {
		min[n] = 1e37;
		max[n] = -1e37;
	}


	for ( int n = 0; n < np; n++) {

//load a single projection
		p = io_read_image(filename, 1, n);
		datasize = (unsigned long) p->x*p->y*p->c;
//calculate the stats
    		switch ( p->datatype ) {
    		case UChar:
			if ( p->colormodel == Bit ) {
				for ( j=0; j<datasize; i++, j++ ) {
//					bit = 0x80 & ( udata[i/8] << i%8 );
					bit = 0;
					if ( 0x80 & ( udata[i/8] << i%8 ) ) bit = 1;
					sum[n] += bit;
					ssum[n] += bit*bit;
				}
				min[n] = 0;
				max[n] = 1;
				p->min = 0;
				p->max = 1;
			} else {
				for ( j=0; j<datasize; i++, j++ ) {
					if ( udata[i] < min[n] ) min[n] = udata[i];
					if ( udata[i] > max[n] ) max[n] = udata[i];
					sum[n] += udata[i];
					ssum[n] += (double)udata[i]*udata[i];
				}
    	    		}
    	    		break;
    		case SChar:
			for ( j=0; j<datasize; i++, j++ ) {
				if ( cdata[i] < min[n] ) min[n] = cdata[i];
				if ( cdata[i] > max[n] ) max[n] = cdata[i];
				sum[n] += cdata[i];
				ssum[n] += (double)cdata[i]*cdata[i];
			}
    	    		break;
    		case UShort:
			for ( j=0; j<datasize; i++, j++ ) {
				if ( usdata[i] < min[n] ) min[n] = usdata[i];
				if ( usdata[i] > max[n] ) max[n] = usdata[i];
				sum[n] += usdata[i];
				ssum[n] += (double)usdata[i]*usdata[i];
			}
	    	    break;
    		case Short:
			for ( j=0; j<datasize; i++, j++ ) {
				if ( sdata[i] < min[n] ) min[n] = sdata[i];
				if ( sdata[i] > max[n] ) max[n] = sdata[i];
				sum[n] += sdata[i];
				ssum[n] += (double)sdata[i]*sdata[i];
			}
	    	    break;
		case Int:
			for ( j=0; j<datasize; i++, j++ ) {
				if ( idata[i] < min[n] ) min[n] = idata[i];
				if ( idata[i] > max[n] ) max[n] = idata[i];
				sum[n] += idata[i];
				ssum[n] += (double)idata[i]*idata[i];
			}
			break;
	    	case Float:
			for ( j=0; j<datasize; i++, j++ ) {
				if ( fdata[i] < min[n] ) min[n] = fdata[i];
				if ( fdata[i] > max[n] ) max[n] = fdata[i];
				sum[n] += fdata[i];
				ssum[n] += (double)fdata[i]*fdata[i];
			}
	    	    break;
    		case ComplexShort:
			for ( j=0; j<datasize; i++, j++ ) {
				amp = (double)csdata[i].re*csdata[i].re + (double)csdata[i].im*csdata[i].im;
				amp = sqrt(amp);
				if ( amp < min[n] ) min[n] = amp;
				if ( amp > max[n] ) max[n] = amp;
				sum[n] += amp;
				ssum[n] += amp*amp;
			}
	    	    break;
		case ComplexInt:
			for ( j=0; j<datasize; i++, j++ ) {
				amp = (double)cidata[i].re*cidata[i].re + (double)cidata[i].im*cidata[i].im;
				amp = sqrt(amp);
				if ( amp < min[n] ) min[n] = amp;
				if ( amp > max[n] ) max[n] = amp;
				sum[n] += amp;
				ssum[n] += amp*amp;
			}
	    	    break;
    		case ComplexFloat:
			for ( j=0; j<datasize; i++, j++ ) {
				amp = sqrt(cfdata[i].re*cfdata[i].re + cfdata[i].im*cfdata[i].im);
				if ( amp < min[n] ) min[n] = amp;
				if ( amp > max[n] ) max[n] = amp;
				sum[n] += amp;
				ssum[n] += amp*amp;
			}
	    	    break;
	    	case Polar:
			for ( j=0; j<datasize; i++, j++ ) {
				if ( pdata[i].amp < min[n] ) min[n] = pdata[i].amp;
				if ( pdata[i].amp > max[n] ) max[n] = pdata[i].amp;
				sum[n] += pdata[i].amp;
				ssum[n] += (double)pdata[i].amp*pdata[i].amp;
    	    		}
    	    		break;
    		default:
				fprintf(stderr, "Warning: Data type %d not supported!", p->datatype);
				p->min = p->max = 0;
				break;
	    }
	    
		kill_img(p);
	}
	b->min = min[0];
	b->max = max[0];
	tsum = tssum = 0;
	for ( n=0; n<np; n++ ) {
		if ( b->min > min[n] ) b->min = min[n];
		if ( b->max < max[n] ) b->max = max[n];
		tsum += sum[n];
		tssum += ssum[n];
	}

    b->avg = tsum/(np*datasize);
    b->std = tssum/(np*datasize) - b->avg*b->avg;
    if ( b->std > 0 ) b->std = sqrt(b->std);
	else b->std = 0;
	b->smin = b->min;
	b->smax = b->max;

	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_stats: tsum = %g tssum = %g\n", tsum, tssum);
	
	if ( verbose & VERB_STATS ) {
	    printf("Data size:                      %ld x %ld x %ld x %ld x %ld = %ld\n",
				b->x, b->y, b->z, b->c, b->n, b->n*datasize);
	    printf("Min, max, avg, std:             %g %g %g %g\n", 
				b->min, b->max, b->avg, b->std);
		printf("\nImage\tMin\tMax\tAvg\tStd\n");
		for ( n=0; n<np; n++ ) {
			avg = sum[n]/datasize;
			std = ssum[n]/datasize - avg*avg;
			if ( std > 0 ) std = sqrt(std);
			else std = 0;
			printf("%ld\t%g\t%g\t%g\t%g\n", n+1, min[n], max[n], avg, std);
		}
		printf("\n");
	}
	
	bfree(min, np*sizeof(double));
	bfree(max, np*sizeof(double));
	bfree(sum, np*sizeof(double));
	bfree(ssum, np*sizeof(double));

// save stats into header of file
	Proj*	pp = proj_init(1);
	proj_cpy_img_sett(pp, b);
	proj_putimage(pp, b, 0, 0);

	
	ierr = io_pwrite(filename, pp, IO_HEADER, np);

	proj_kill(pp);
	kill_img(b);
	
	return(ierr);
}
/************************************************************************
@Function: io_get_index
@Description:
	Calculate data index associated to voxel.
@Algorithm:
	If the buffer is initialized, the index is compared to the range
	of those loaded into buffer. If out of the range, the buffer is
	refreshed.
@Arguments:
	int x			voxel coordinate.
	int y 			voxel coordinate.
	int z 			voxel coordinate.
	Bimage* p		image structure.
	Buff* bf			buffer structure.
@Returns:
	long long	data index.
**************************************************************************/
long long io_get_index(int x, int y, int z, Bimage* p, Buff* bf)
{
	long long i;
	
	i = (z*p->y + y)*p->x + x;
	if(bf) {
		if(i>bf->last || i<bf->first) {
			bf->first = i;
			bf->last = i+1;
			if ( verbose & VERB_DEBUG ) {
				printf("DEBUG io_get_index: index data missing: %lld (%lld). Refreshing buffer\n", bf->first, bf->first*gettypesize(p->datatype));
			}

			p = io_read_image_buff(p, (char*) NULL, -1, bf);

			i = 0;				
		} else
			i -= bf->first;
	}
	
	return(i);	
}

/************************************************************************
@Function: io_read_data
@Description:
	Read image data in a generalized style.
@Algorithm:
	The whole file or a single image from a file may be read.
	The data is read in the largest blocks possible for efficiency.
	Any interspersed padding and page sizes not matching the data size
	contribute to inefficiency in reading.
	Swapping:	swap = 1:	swap bytes the size of the data type
				swap > 1:	swap these number of bytes regardless of the data type
@Arguments:
	FILE* fimg			file descriptor: file opened in calling function.
	Bimage* p			image structure: header parameters filled by calling function.
	int select_img		image selection: if -1, all images, if >= 0, one image.
	int swap			flag activates byte swapping.
	int vax 			indicate vax style floating point - activates conversion.
	int pad 			any interspersed control or separation bytes.
@Returns:
	char*				data in one big contiguous block, NULL if reading failed.
**************************************************************************/
char*		io_read_data(FILE* fimg, Bimage* p, int select_img, int swap, int vax, int pad)
{
	if ( !p ) return(NULL);
	
	if ( p->dataflag < 1 ) return(NULL);
	
	if ( p->px < 1 ) p->px = p->x;
	if ( p->py < 1 ) p->py = p->y;
	if ( p->pz < 1 ) p->pz = p->z;
	
	// If only half of a transform is stored, it need to be handled
	unsigned long		xstore = p->x;
	unsigned long		xpage = p->px;
	if ( p->transform == Hermitian || p->transform == CentHerm ) {
		xstore = p->x/2 + 1;
		if ( p->px > xstore ) xpage = xstore;
	}
	
	// Reset select to get the correct offset
	if ( select_img < 0 ) select_img = 0;
	
	unsigned long		i, iz, iy, j, jz, jy, n, x, y, z, px, py, pz, npages;
	unsigned long		datatypesize = gettypesize(p->datatype);
	unsigned long		valuesize = datatypesize*p->c;
	long long			datasize = p->n*xstore*p->y*p->z*valuesize;
	long long			pagesize = xpage*p->py*p->pz*valuesize;
	char*				data = (char *) balloc(datasize*sizeof(char));
	char*				page = NULL;

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG io_read_data: Pointer = %p  Data type size = %ld\n", p, datatypesize);
		printf("DEBUG io_read_data: Data size: %ld %ld %ld %ld %ld (%lld)\n", p->x, p->y, p->z, p->c, p->n, datasize);
		printf("DEBUG io_read_data: Page size: %ld %ld %ld (%lld)\n", p->px, p->py, p->pz, pagesize);
		printf("DEBUG io_read_data: Swap = %d  Vax = %d  Pad = %d  Offset = %ld\n", swap, vax, pad, p->offset);
	}
	
	if ( !pad && p->x == p->px && p->y == p->py && p->z == p->pz ) { // 3D block
		if ( verbose & VERB_DEBUG )
			printf("DEBUG io_read_data: Reading 3D blocks\n");
		fseeko( fimg, (off_t) p->offset + select_img*(pagesize + pad), SEEK_SET );
		if ( pad ) page = (char *) balloc(pad*sizeof(char));
		for ( n=0; n<p->n; n++ ) {
			fread( data + n*pagesize, pagesize, 1, fimg );
			if ( pad ) fread( page, pad, 1, fimg);
		}
		if ( page ) bfree(page, pad*sizeof(char));
	} else if ( p->x == p->px && p->y == p->py && p->pz == 1 ) {	// 2D page (MFF)
		if ( verbose & VERB_DEBUG )
			printf("DEBUG io_read_data: Reading 2D pages\n");
		fseeko( fimg, (off_t) p->offset + select_img*p->z*(pagesize + pad), SEEK_SET );
		if ( pad ) page = (char *) balloc(pad*sizeof(char));
		for ( n=0; n<p->n; n++ ) {
			for ( pz=0; pz<p->z; pz+=p->pz ) {
				i = (n*p->z + pz)*pagesize;
				fread( data+i, pagesize, 1, fimg );
				if ( pad ) fread( page, pad, 1, fimg);
			}
		}
		if ( page ) bfree(page, pad*sizeof(char));
	} else {					// More complicated paging (DSN6 and BRIX)
		if ( verbose & VERB_DEBUG )
			printf("DEBUG io_read_data: Reading single values in pages of size %lld\n", 
					pagesize);
		npages = 1;
		npages *= (p->z - 1)/p->pz + 1;
		npages *= (p->y - 1)/p->py + 1;
		npages *= (xstore - 1)/xpage + 1;
		fseeko( fimg, (off_t) p->offset + select_img*npages*(pagesize + pad), SEEK_SET );
		page = (char *) balloc(pagesize*sizeof(char));
		for ( n=0; n<p->n; n++ ) {
			for ( pz=0; pz<p->z; pz+=p->pz ) {
				for ( py=0; py<p->y; py+=p->py ) {
					for ( px=0; px<xstore; px+=xpage ) {//printf("%d %d %d\n", px, py, pz);
						fread( page, pagesize, 1, fimg );
						if ( swap > 1 )
							 for ( i=0; i<pagesize; i+=swap ) swapbytes(page+i, swap);
						for ( z=0; z<p->pz && z<p->z-pz; z++ ) {
							iz = (n*p->z + z + pz)*p->y;
							jz = z*p->py;
							for ( y=0; y<p->py && y<p->y-py; y++ ) {
								iy = (iz + y + py)*xstore;
								jy = (jz + y)*xpage;
								for ( x=0; x<xstore && x<xstore-px; x++ ) {
									i = (iy + x + px)*valuesize;
									j = (jy + x)*valuesize;
									memcpy(data+i, page+j, valuesize);
								}
							}
						}
						if ( pad ) fread( page, pad, 1, fimg);
					}
				}
			}
		}
		if ( page ) bfree(page, pagesize*sizeof(char));
	}
			
	if ( verbose & VERB_DEBUG )
		printf("DEBUG io_read_data: Converting image data\n");

	// Convert vax format and swap bytes if required
    if ( vax && swap < 2 && ( p->datatype == Float ) )
		for ( i=0; i<datasize; i+=4 )
    	    vax2ieee(data+i, 1-swap);
    else if ( swap == 1 )
		for ( i=0; i<datasize; i+=datatypesize ) 
			swapbytes(data+i, (unsigned int) datatypesize);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG io_read_data: Finished reading and converting data\n");

	return(data);
}

/************************************************************************
@Function: io_write_image
@Description:
	Driver function to write images to file
@Algorithm:
	This is the only image writing function that should be called
	from programs.
@Arguments:
	char* filename		file name.
	Bimage*				the image structure.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 		io_write_image(char* filename, Bimage* p)
{
	if ( !p ) return(-1);
	
	if ( !filename || !strlen(filename) ) {
		fprintf(stderr, "Error: No output image file name given!\n");
		return(-1);
	}
	
	int 		err = 0;

	filename_clean_copy(filename, p->filename, 256);
	
	img_stats(p);
	
	img_check_param(p);
	
    p->time = time(NULL);
	
	if ( strlen(p->label) < 2 )
		strcpy(p->label, "Written by Electra");
	
	char*		ext = extension(filename);
	if ( !ext ) {
		fprintf(stderr, "Error: No extension found in the file name \"%s\"!\n\n", filename);
		error_show(filename, __FILE__, __LINE__);
		return(-2);
	}
	
	if ( p->colormodel == RGB ) {
		fprintf(stderr, "Warning: RGB output format not accepted!\n");
//		error_show(filename, __FILE__, __LINE__);
//		return(-3);
	}

	if ( verbose & VERB_DEBUG )
		printf("DEBUG io_write_image: Image setup done for file %s %s (extension %s)\n", 
				filename, p->filename, ext);

//	if ( strstr(ext,"mrc") )
	err = writeMRC(p);
//	else {
//		printf("File format with extension \"%s\" not supported!\n\n", ext);
//		err = -1;
//	}
	
	bfree_string(ext);
	
	if ( err < 0 ) {
		error_show(filename, __FILE__, __LINE__);
		return(err);
	}
	
	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) )
		printf("Writing file:                   %s\n", p->filename);
	
	if ( verbose & VERB_PROCESS )
		img_info(p);
	
	if ( verbose & VERB_STATS )
		sub_img_info(p);
	
	return(err);
}
