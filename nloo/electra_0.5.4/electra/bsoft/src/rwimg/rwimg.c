/*
	rwimg.c
	Library for 2D and 3D image I/O
	Author: Bernard Heymann
	Created: 19990321 	Modified: 20050203
*/

#include "rwimg.h"
#include "symmetry.h"
#include "img_complex.h"
#include "img_util.h"
#include "matrix.h"
#include "utilities.h"

/*#include "rwRAW.h"
#include "rwASCII.h"
#include "rwBIORAD.h"
#include "rwBRIX.h"
*/
#include "rwCCP4.h"
/*#include "rwDI.h"
#include "rwDM.h"
#include "rwDSN6.h"
#include "rwEM.h"
#include "rwGOODFORD.h"
#include "rwGRD.h"
#include "rwHKL.h"
#include "rwIMAGIC.h"
#include "rwIP.h"
#include "rwJPEG.h"
#include "rwMFF.h"
#include "rwMIFF.h"
*/
#include "rwMRC.h"
/*#include "rwPIC.h"
#include "rwPIF.h"
#include "rwPNG.h"
#include "rwPostScript.h"
#include "rwSPIDER.h"
#include "rwSUPRIM.h"
#include "rwTIFF.h"
#include "rwXPLOR.h"
#include "rwimg_star.h"
*/

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

// Internal function prototypes

/************************************************************************
@Function: init_img
@Description:
	General image structure initialization.
@Algorithm:
	Whenever a new image is read or created, this function should be
	called first. It allocates memory and sets up a number of defaults.
@Arguments:
	.
@Returns:
	Bimage*				new image structure, NULL if initialization failed.
**************************************************************************/
Bimage* 	init_img()
{
	// Allocate memory for the image parameter structure
	Bimage* 	p = (Bimage *) balloc(sizeof(Bimage));
	if ( !p ) return(p);

	// Set parameter defaults
	p->time = time(NULL);
	p->x = p->y = p->z = p->c = p->n = 1;
	p->i = p->px = p->py = p->pz = 0;
	p->ux = p->uy = p->uz = 1;
	p->alf = p->bet = p->gam = PI/2;
	p->datatype = UChar;
	p->transform = NoTransform;
	p->colormodel = Gray;
	p->spacegroup = 0;
	p->dataflag = 0;
	p->label = (char *) balloc(IMGLABELSIZE*sizeof(char));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG init_img: Image setup done\n");
	
	return(p);
}

/************************************************************************
@Function: init_img_header
@Description:
	General image header structure initialization.
@Algorithm:
	Whenever a new image is read or created, this function should be
	called first. It allocates memory and sets up a number of defaults.
	The sub-image structures are also allocated, but not the data block.
@Arguments:
	DataType datatype	the new data type.
	unsigned long c		number of channels.
	unsigned long x		x-dimension.
	unsigned long y		y-dimension.
	unsigned long z		z-dimension.
	unsigned long n		number of images.		
@Returns:
	Bimage*				the new image structure, NULL if initialization failed.
**************************************************************************/
Bimage* 	init_img_header(DataType datatype, unsigned long c, unsigned long x, 
				unsigned long y, unsigned long z, unsigned long n)
{
	if ( n < 1 ) return(NULL);
	
	// Call the more primitive initialization
	Bimage* 	p = init_img();
	if ( !p ) return(p);

	// Set new parameters
	p->datatype = datatype;
	p->x = p->px = x;
	p->y = p->py = y;
	p->z = p->pz = z;
	p->c = c;
	p->n = n;
	p->image = (Bsub_image *) balloc(n*sizeof(Bsub_image));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG init_img_header: Image created\n");
	
	return(p);
}

/************************************************************************
@Function: init_img_with_parameters
@Description:
	General image structure initialization with data memory allocation.
@Algorithm:
	Whenever a new image is read or created, this function should be
	called first. It allocates memory and sets up a number of defaults.
	The sub-image structures and data block are also allocated.
@Arguments:
	DataType datatype	the new data type.
	unsigned long c		number of channels.
	unsigned long x		x-dimension.
	unsigned long y		y-dimension.
	unsigned long z		z-dimension.
	unsigned long n		number of images.		
@Returns:
	Bimage*				the new image structure, NULL if initialization failed.
**************************************************************************/
Bimage* 	init_img_with_parameters(DataType datatype, unsigned long c, unsigned long x, 
				unsigned long y, unsigned long z, unsigned long n)
{
	if ( n < 1 ) return(NULL);
	
	// Call the more primitive initialization
	Bimage* 	p = init_img_header(datatype, c, x, y, z, n);
	if ( !p ) return(p);
	unsigned long   datasize = c*x*y*z*n*gettypesize(datatype);

	// Allocate data memory
	p->dataflag = 1;
	p->data = (char *) balloc(datasize*sizeof(char));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG init_img_with_parameters: Image created with data allocation\n");
	
	return(p);
}

/************************************************************************
@Function: img_init_fom
@Description:
	Initializes the FOM block.
@Algorithm:
	If a FOM block exists, it is deleted.
@Arguments:
	Bimage* p			image.
@Returns:
	float*				the new FOM block.
**************************************************************************/
float*		img_init_fom(Bimage* p)
{
	unsigned long   datasize = (unsigned long)p->x*p->y*p->z*p->n*sizeof(float);
	
	if ( p->fomflag && p->fom ) bfree(p->fom, datasize);
	
	p->fomflag = 1;
	p->fom = (float *) balloc(datasize);
	
	return(p->fom);
}

/************************************************************************
@Function: read_img
@Description:
	General driver function to read multiple image formats
@Algorithm:
	This is the only image reading function that should be called
	from programs.
	A Bimage structure is initialized with default values.
	The file format is deduced from the file name extension.
	Every file format has its own funtion to read the file.
	The selection argument is used to read only one image from a multi-image
	file if it is greater than -1. This selection must be handled within 
	each format to ensure the correct allocation of the sub-image structure.
	If the requested selection is equal or larger than the number of
	images, the selection is set to the last image.
@Arguments:
	char* filename		file name (plus any tags for the RAW format).
	int readdataflag	flag to activate reading of image data.
	int select_img		image selection in multi-image file (-1 = all images).
@Returns:
	Bimage*				the image structure, NULL if reading failed.
**************************************************************************/
Bimage* 	read_img(char* filename, int readdataflag, int select_img)
{
	if ( verbose & VERB_DEBUG )
		printf("DEBUG read_img: Filename = %s\n", filename);
	
	if ( !filename || !strlen(filename) ) {
		fprintf(stderr, "Error: No input image file given!\n");
		return(NULL);
	}
	
	int 		err = 0;
	char*		ext = NULL;
	char		clean_filename[256];
	if ( filename && strlen(filename) ) ext = extension(filename);
	filename_clean_copy(filename, clean_filename, 255);
	
	if ( access(clean_filename, F_OK) ) {
		error_show(clean_filename, __FILE__, __LINE__);
		bfree_string(ext);
		return(NULL);
	}

	Bimage* 	p = init_img();
	if ( strstr(filename, "#") )
		strncpy(p->filename, filename, 255);
	else
		strncpy(p->filename, clean_filename, 255);
	if ( readdataflag ) p->dataflag = 1;
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG read_img: Transferred filename = %s\n", p->filename);
		printf("DEBUG read_img: Extension = %s\n", ext);
	}
	
/*	if ( strstr(filename, "#") || !ext )
		err = readRAW(p, select_img);
	else if ( strstr(ext,"raw") )
		err = readRAW(p, select_img);
	else if ( strstr(ext,"star") || strstr(ext,"ent") )
		err = read_img_star(p);
	else if ( strstr(ext,"asc") || strstr(ext, "txt") )
		err = readASCII(p);
	else if ( strstr(ext,"pic") )
		err = readBIORAD(p);
	else if ( strstr(ext,"brx" ) )
		err = readBRIX(p);
	else if ( strstr(ext,"ccp") || strstr(ext,"map") )
		err = readCCP4(p);
	else if ( strstr(ext,"di") )
		err = readDI(p, select_img);
	else if ( strstr(ext,"dm") )
		err = readDM(p);
	else if ( strstr(ext,"omap" ) || strstr(ext,"dsn6") || strstr(ext,"dn6") )
		err = readDSN6(p);
	else if ( strstr(ext,"em") )
		err = readEM(p);
	else if ( strstr(ext,"pot") )
		err = readGOODFORD(p);
	else if ( strstr(ext,"grd") )
		err = readGRD(p);
	else if ( strstr(ext,"hkl") )
		err = readHKL(p);
	else if ( strstr(ext,"img") || strstr(ext,"hed") )
		err = readIMAGIC(p, select_img);
	else if ( strstr(ext,"ip") )
		err = readIP(p);
	else if ( strstr(ext,"jpg") || strstr(ext,"jpeg") )
		err = readJPEG(p);
	else if ( strstr(ext,"mif") )
		err = readMIFF(p, select_img);
	else if ( strstr(ext,"mff") )
		err = readMFF(p);
	else*/ if ( strstr(ext,"mrc") )
		err = readMRC(p);
/*	else if ( strstr(ext,"pif") || strstr(ext,"sf") )
		err = readPIF(p, select_img);
	else if ( strstr(ext,"bp") || strstr(ext,"bq") )
		err = readPIC(p);
	else if ( strstr(ext,"png") )
		err = readPNG(p);
	else if ( strstr(ext,"spi") )
		err = readSPIDER(p, select_img);
	else if ( strstr(ext,"spm") || strstr(ext,"sup") )
		err = readSUPRIM(p);
	else if ( strstr(ext,"tif") )
		err = readTIFF(p, select_img);
	else if ( strstr(ext,"xpl") || strstr(ext,"rfl") )
		err = readXPLOR(p);
*/	else {
		fprintf(stderr, "Error: File format with extension \"%s\" not supported!\n", ext);
		err = -1;
	}
	
	bfree_string(ext);
		
	if ( err < 0 ) {
		error_show(filename, __FILE__, __LINE__);
		kill_img(p);
		return(NULL);
	}
	
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
@Function: img_read_data
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
char*		img_read_data(FILE* fimg, Bimage* p, int select_img, int swap, int vax, int pad)
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
	unsigned long		datasize = (unsigned long) p->n*xstore*p->y*p->z*valuesize;
	unsigned long		pagesize = (unsigned long) xpage*p->py*p->pz*valuesize;
	char*				data = (char *) balloc(datasize*sizeof(char));
	char*				page = NULL;

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG img_read_data: Pointer = %p  Data type size = %ld\n", p, datatypesize);
		printf("DEBUG img_read_data: Data size: %ld %ld %ld %ld %ld (%ld)\n", p->x, p->y, p->z, p->c, p->n, datasize);
		printf("DEBUG img_read_data: Page size: %ld %ld %ld (%ld)\n", p->px, p->py, p->pz, pagesize);
		printf("DEBUG img_read_data: Swap = %d  Vax = %d  Pad = %d  Offset = %ld\n", swap, vax, pad, p->offset);
	}
	
	if ( !pad && p->x == p->px && p->y == p->py && p->z == p->pz ) { // 3D block
		if ( verbose & VERB_DEBUG )
			printf("DEBUG img_read_data: Reading 3D blocks\n");
		fseek( fimg, p->offset + select_img*(pagesize + pad), SEEK_SET );
		if ( pad ) page = (char *) balloc(pad*sizeof(char));
		for ( n=0; n<p->n; n++ ) {
			fread( data + n*pagesize, pagesize, 1, fimg );
			if ( pad ) fread( page, pad, 1, fimg);
		}
		if ( page ) bfree(page, pad*sizeof(char));
	} else if ( p->x == p->px && p->y == p->py && p->pz == 1 ) {	// 2D page (MFF)
		if ( verbose & VERB_DEBUG )
			printf("DEBUG img_read_data: Reading 2D pages\n");
		fseek( fimg, p->offset + select_img*p->z*(pagesize + pad), SEEK_SET );
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
			printf("DEBUG img_read_data: Reading single values in pages of size %ld\n", 
					pagesize);
		npages = 1;
		npages *= (p->z - 1)/p->pz + 1;
		npages *= (p->y - 1)/p->py + 1;
		npages *= (xstore - 1)/xpage + 1;
		fseek( fimg, p->offset + select_img*npages*(pagesize + pad), SEEK_SET );
		page = (char *) balloc(pagesize*sizeof(char));
		for ( n=0; n<p->n; n++ ) {
			for ( pz=0; pz<p->z; pz+=p->pz ) {
				for ( py=0; py<p->y; py+=p->py ) {
					for ( px=0; px<xstore; px+=xpage ) {//printf("%d %d %d\n", px, py, pz);
						fread( page, pagesize, 1, fimg );
						if ( swap > 1 )
							 for ( i=0; i<pagesize; i+=swap ) swapbytes(page+i, (int) swap);
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
		printf("DEBUG img_read_data: Converting image data\n");

	// Convert vax format and swap bytes if required
    if ( vax && swap < 2 && ( p->datatype == Float ) )
		for ( i=0; i<datasize; i+=4 )
    	    vax2ieee(data+i, 1-swap);
    else if ( swap == 1 )
		for ( i=0; i<datasize; i+=datatypesize ) 
			swapbytes(data+i, (int) datatypesize);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_read_data: Finished reading and converting data\n");

	return(data);
}

/************************************************************************
@Function: write_img
@Description:
	General driver function to write multiple image formats
@Algorithm:
	This is the only image writing function that should be called
	from programs.
	The file format is deduced from the file name extension.
	Every file format has its own funtion to write the file.
@Arguments:
	char* filename		file name (plus any tags for the RAW format).
	Bimage*				the image structure.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 		write_img(char* filename, Bimage* p)
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
		strcpy(p->label, "Written by Bsoft");
	
	char*		ext = extension(filename);
	if ( !ext ) {
		fprintf(stderr, "Error: No extension found in the file name \"%s\"!\n\n", filename);
		error_show(filename, __FILE__, __LINE__);
		return(-2);
	}
	
	if ( p->colormodel == RGB ) {
		if ( !strstr(ext,"tif") && !strstr(ext,"jp") && !strstr(ext,"mif") && !strstr(ext,"ps") ) {
			fprintf(stderr, "Warning: Only TIFF, JPEG, MIFF or postscript files are acceptable as RGB output formats!\n");
//			error_show(filename, __FILE__, __LINE__);
//			return(-3);
		}
	}

	if ( verbose & VERB_DEBUG )
		printf("DEBUG write_img: Image setup done for file %s %s (extension %s)\n", 
				filename, p->filename, ext);

/*	if ( strstr(ext,"raw") )
		err = writeRAW(p);
	else if ( strstr(ext,"star") )
		err = write_img_star(p);
	else if ( strstr(ext,"asc") || strstr(ext, "txt") )
		err = writeASCII(p);
	else if ( strstr(ext,"pic") )
		err = writeBIORAD(p);
	else if ( strstr(ext,"brx" ) )
		err = writeBRIX(p);
	else if ( ( strstr(ext,"ccp") || strstr(ext,"map") ) )
		err = writeCCP4(p);
	else if ( strstr(ext,"omap" ) || strstr(ext,"dsn6") || strstr(ext,"dn6") )
		err = writeDSN6(p);
	else if ( strstr(ext,"em") )
		err = writeEM(p);
	else if ( strstr(ext,"dm") )
		err = writeDM(p);
	else if ( strstr(ext,"pot") )
		err = writeGOODFORD(p);
	else if ( strstr(ext,"grd") )
		err = writeGRD(p);
	else if ( strstr(ext,"hkl") )
		err = writeHKL(p);
	else if ( strstr(ext,"img") || strstr(ext,"hed") )
		err = writeIMAGIC(p);
	else if ( strstr(ext,"ip") )
		err = writeIP(p);
	else if ( strstr(ext,"jpg") || strstr(ext,"jpeg") )
		err = writeJPEG(p);
	else if ( strstr(ext,"mff") )
		err = writeMFF(p);
	else if ( strstr(ext,"mif") )
		err = writeMIFF(p);
	else*/ if ( strstr(ext,"mrc") )
		err = writeMRC(p);
/*	else if ( strstr(ext,"pif") || strstr(ext,"sf") )
		err = writePIF(p);
	else if ( strstr(ext,"bp") || strstr(ext,"bq") )
		err = writePIC(p);
	else if ( strstr(ext,"png") )
		err = writePNG(p);
	else if ( strstr(ext,"ps") )
		err = writePostScriptImage(p);
	else if ( strstr(ext,"spi") )
		err = writeSPIDER(p);
	else if ( strstr(ext,"spm") || strstr(ext,"sup") )
		err = writeSUPRIM(p);
	else if ( strstr(ext,"tif") )
		err = writeTIFF(p, p->z-1);
	else if ( strstr(ext,"xpl") || strstr(ext,"rfl") )
		err = writeXPLOR(p);
*/	else {
		printf("File format with extension \"%s\" not supported!\n\n", ext);
		err = -1;
	}
	
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

/************************************************************************
@Function: write_fom_as_image
@Description:
	Converts a FOM block to an image for output
@Algorithm:
	A new image is generated with the FOM data block.
@Arguments:
	char* filename		file name.
	Bimage*				the image structure with FOM defined.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int			write_fom_as_image(char* filename, Bimage* p)
{
	if ( !p->fomflag || !p->fom ) {
		fprintf(stderr, "Error: No FOM block found in %s!\n", p->filename);
		return(-1);
	}
	
	Bimage* 	pfom = copy_img_header(p, -1);
	pfom->c = 1;
	pfom->datatype = Float;
	pfom->dataflag = 1;
	pfom->data = (char *) p->fom;
	pfom->fomflag = 0;
	pfom->fom = NULL;
	
	write_img(filename, pfom);
	
	pfom->dataflag = 0;
	pfom->data = NULL;
	kill_img(pfom);
	
	return(0);
}

/************************************************************************
@Function: write_data_block_as_image
@Description:
	Converts a data block to an image for output.
@Algorithm:
	A new image is generated with the data block.
	The data block is not deallocated.
@Arguments:
	char* filename		file name.
	char* data			data block
	DataType datatype	the new data type.
	unsigned long c		number of channels.
	unsigned long x		x-dimension.
	unsigned long y		y-dimension.
	unsigned long z		z-dimension.
	unsigned long n		number of images.		
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int			write_data_block_as_image(char* filename, char* data, DataType datatype,
				unsigned long c, unsigned long x, unsigned long y, unsigned long z, unsigned long n)
{
	if ( !data ) {
		fprintf(stderr, "Error: No data block found!\n");
		return(-1);
	}
	
	Bimage* 	p = init_img_header(datatype, c, x, y, z, n);
	p->dataflag = 1;
	p->data = data;
	
	write_img(filename, p);
	
	p->dataflag = 0;
	p->data = NULL;
	kill_img(p);
	
	return(0);
}

/************************************************************************
@Function: img_write_data
@Description:
	Writes a datablock with optional swapping.
@Algorithm:
	If the swap flag is set, a temporary data block is created, swapped,
	written and freed. This may cause memory problems on small computers.
@Arguments:
	FILE* fimg				file descriptor.
	char* data				data pointer.
	unsigned long datasize	data elements.
	unsigned long typesize	data type size (used for swapping).
	int swap				swap flag.
@Returns:
	int						error code (<0 means failure).
**************************************************************************/
int			img_write_data(FILE* fimg, char* data, unsigned long datasize, 
				unsigned long typesize, int swap)
{
	char			*temp, *aptr;
	unsigned long	i;
	
	if ( swap ) {
		temp = (char *) balloc(datasize*typesize);
		memcpy(temp, data, datasize*typesize);
		for ( aptr=temp, i=0; i<datasize*typesize; i+=typesize, aptr+=typesize )
			swapbytes(aptr, (int) typesize);
		fwrite(temp, datasize, typesize, fimg);
		bfree(temp, datasize*typesize);
	} else {
		fwrite(data, datasize, typesize, fimg);
	}
	
	return(0);
}

/************************************************************************
@Function: copy_img_header
@Description:
	Copy an image structure into a new one.
@Algorithm:
	All information from an old image structure is copied into a new one.
	This takes care of the internal arrays in the structure.
	The data pointer is left to point to the original data.
@Arguments:
	Bimage*	p			the original image structure.
	int new_nimg		new number of images (if < 1, keep the old number)
@Returns:
	Bimage*				the new image structure, NULL if copy failed.
**************************************************************************/
Bimage* 	copy_img_header(Bimage* p, int new_nimg)
{
	if ( !p ) return(NULL);
	
	unsigned int	i;
	
	// Allocate memory for the image parameter structure
	Bimage*			pnew = (Bimage *) balloc(sizeof(Bimage));
	if ( !pnew ) return(pnew);
	memcpy(pnew, p, sizeof(Bimage));

	// Clear data and FOM flags and pointers
	pnew->dataflag = 0;
	pnew->data = NULL;
	pnew->fomflag = 0;
	pnew->fom = NULL;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG copy_img_header: New image allocated and flags cleared\n");
	
	// Create and copy new sub-images
	if ( new_nimg > 0 ) pnew->n = new_nimg;
	else pnew->n = p->n;
	if ( verbose & VERB_DEBUG )
		printf("DEBUG copy_img_header: New number of images = %ld\n", pnew->n);
	pnew->image = (Bsub_image *) balloc(pnew->n*sizeof(Bsub_image));
	if ( pnew->n == p->n )
		for ( i=0; i<pnew->n; i++ ) pnew->image[i] = p->image[i];
	else
		for ( i=0; i<pnew->n; i++ ) pnew->image[i] = p->image[0];
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG copy_img_header: %ld sub-images allocated and set\n", pnew->n);
	
	// Allocate and copy internal arrays
	pnew->label = (char *) balloc(IMGLABELSIZE*sizeof(char));
	memcpy(pnew->label, p->label, IMGLABELSIZE);
	
	if ( p->colormap ) {
		pnew->colormap = (char *) balloc(p->colors*3*sizeof(char));
		memcpy(pnew->colormap, p->colormap, p->colors*3*sizeof(char));
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG copy_img_header: Image copy with %ld sub-images done\n", pnew->n);
	
	return(pnew);
}

/************************************************************************
@Function: copy_img
@Description:
	Copies the header information and data of an image into a new image structure.
@Algorithm:

@Arguments:
	Bimage* p			image to be copied.
@Returns:
 	Bimage* 			copied image, NULL if copy failed.
**************************************************************************/
Bimage* 	copy_img (Bimage* p)
{
	Bimage* 	pcopy = copy_img_header(p, -1);
	if ( !pcopy ) return(pcopy);

	unsigned long   datasize = p->x*p->y*p->z*p->n*p->c*gettypesize(p->datatype);
	unsigned long   fomsize = p->x*p->y*p->z*p->n*sizeof(float);

	pcopy->dataflag = 1;
	pcopy->data = (char *) balloc(datasize*sizeof(char));
	memcpy(pcopy->data, p->data, datasize);
	
	if ( p->fomflag && p->fom ) {
		img_init_fom(pcopy);
		memcpy(pcopy->fom, p->fom, fomsize);
	}

	return(pcopy);
}

/************************************************************************
@Function: move_img
@Description:
	Moves an image from one pointer to another.
@Algorithm:
	The object is to use a predefined pointer to move an image into.
	The original image pointer is deallocated.
@Arguments:
	Bimage* pfrom		original image pointer (deallocated).
	Bimage* pto			new image pointer.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 		move_img (Bimage* pfrom, Bimage* pto)
{
	if ( !pfrom || !pto ) return(-1);
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG move_img: pfrom = %p  pto = %p\n", pfrom, pto);
		printf("DEBUG move_img: pfrom->data = %p  pto->data = %p\n", 
				pfrom->data, pto->data);
	}

	kill_all_but_main_img(pto);
	
	memcpy(pto, pfrom, sizeof(Bimage));
	
	bfree(pfrom, sizeof(Bimage));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG move_img: pfrom->data = %p  pto->data = %p\n", 
				pfrom->data, pto->data);

	return(0);
}

/************************************************************************
@Function: img_slices_to_images
@Description:
	Changes the slices in a 3D image into a set of 2D images.
@Algorithm:
	.
@Arguments:
	Bimage* p			image pointer.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 		img_slices_to_images(Bimage* p)
{
	if ( !p ) return(-1);
	
	if ( p->n > 1 ) {
		fprintf(stderr, "Error: %s is already a multi-image file!\n", p->filename);
		exit(-1);
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_slices_to_images: %ld slices to images\n", p->z);
	
	unsigned int	n;
	Bsub_image		subimg = p->image[0];

	if ( p->image ) bfree(p->image, sizeof(Bsub_image));
	
	p->n = p->z;
	p->z = 1;
	
	p->image = (Bsub_image *) balloc(p->n*sizeof(Bsub_image));
	
	for ( n=0; n<p->n; n++ ) p->image[n] = subimg;
	
	return(0);
}

/************************************************************************
@Function: img_images_to_slices
@Description:
	Changes the 2D images to slices in a 3D image.
@Algorithm:
	.
@Arguments:
	Bimage* p			image pointer.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 		img_images_to_slices(Bimage* p)
{
	if ( !p ) return(-1);
	
	if ( p->z > 1 ) {
		fprintf(stderr, "Error: %s is already a 3D image file!\n", p->filename);
		exit(-1);
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_images_to_slices: %ld images to slices\n", p->n);
	
	Bsub_image		subimg = p->image[0];

	if ( p->image ) bfree(p->image, p->n*sizeof(Bsub_image));
	
	p->z = p->n;
	p->n = 1;
	
	p->image = (Bsub_image *) balloc(p->n*sizeof(Bsub_image));
	
	p->image[0] = subimg;
	
	return(0);
}

/************************************************************************
@Function: kill_img
@Description:
	General image structure destruction.
@Algorithm:
	All allocated blocks within a structure is freed, as well as the structure.
	All image structures down the linked list are also deallocated.
@Arguments:
	Bimage*				the image structure to be destroyed.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 		kill_img(Bimage* p)
{
	Bimage*		p2;
	
	while ( p ) {
		p2 = p->next;
		kill_all_but_main_img(p);
		bfree(p, sizeof(Bimage));
		p = p2;
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG kill_img: memory = %ld\n", memory);
	
	return(0);
}

/************************************************************************
@Function: kill_all_but_main_img
@Description:
	Deallocate all pointers in an image structure.
@Algorithm:
	All allocated blocks within a structure is freed, but not the structure.
@Arguments:
	Bimage*				the image structure to be deallocated.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 		kill_all_but_main_img(Bimage* p)
{
	if ( p == NULL ) return(0);
	
	if ( p->image ) bfree(p->image, p->n*sizeof(Bsub_image));
	
	if ( p->colormap ) bfree(p->colormap, p->colors*3*sizeof(char));
	
	unsigned long 	xstore = p->x;
	if ( p->transform == Hermitian || p->transform == CentHerm )
		xstore = p->x/2 + 1;
	unsigned long	datasize = (unsigned long) xstore*p->y*p->z*p->n*p->c*gettypesize(p->datatype);
	if ( p->colormodel == Bit ) datasize /= 8;
	if ( p->dataflag && p->data ) bfree(p->data, datasize);
	p->dataflag = 0;
	
	if ( p->label ) bfree(p->label, IMGLABELSIZE);
	
	if ( p->fomflag && p->fom )
		bfree(p->fom, p->x*p->y*p->z*p->n*sizeof(float));
	p->fomflag = 0;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG kill_all_but_main_img: memory = %ld\n", memory);
	
	return(0);
}

/************************************************************************
@Function: img_clear_data
@Description:
	Sets the data block to zero.
@Algorithm:
	.
@Arguments:
	Bimage*				the image structure to be cleared.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 		img_clear_data(Bimage* p)
{
	if ( p == NULL ) return(0);
	
	unsigned long	datasize = p->x*p->y*p->z*p->n*p->c*gettypesize(p->datatype);
	
	if ( p->dataflag ) memset(p->data, 0, datasize);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_clear_data: datasize = %ld\n", datasize);
	
	return(0);
}

/************************************************************************
@Function: img_check_param
@Description:
	Checks an image structure for badly set parameters.
@Algorithm:
	Check parameters for conformance and set defaults if necessary.
@Arguments:
	Bimage*				the image structure.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 		img_check_param(Bimage* p)
{
	if ( !p ) return(-1);
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG img_check_param: Checking parameters\n");
		printf("DEBUG img_check_param: dataflag = %d, p->data = %p\n", 
				p->dataflag, p->data);
		printf("DEBUG img_check_param: fomflag = %d, p->fom = %p\n", 
				p->fomflag, p->fom);
	}
	
	if ( p->dataflag ) if ( !p->data ) p->dataflag = 0;
	if ( p->fomflag ) if ( !p->fom ) p->fomflag = 0;
	
	unsigned int	i;
	
	img_convert_fourier(p, Standard);

	if ( p->c < 1 ) p->c = 1;
	if ( p->px < 1 ) p->px = p->x;
	if ( p->py < 1 ) p->py = p->y;
	if ( p->pz < 1 ) p->pz = p->z;
	if ( p->px > 2*p->x ) p->px = p->x;
	if ( p->py > 2*p->y ) p->py = p->y;
	if ( p->pz > 2*p->z ) p->pz = p->z;
	if ( p->n < 1 ) p->n = 1;
	if ( p->ux <= 1e-6 ) p->ux = 1;
	if ( p->uy <= 1e-6 ) p->uy = 1;
	if ( p->uz <= 1e-6 ) p->uz = 1;
	if ( p->ux > 1e6 ) p->ux = 1;
	if ( p->uy > 1e6 ) p->uy = 1;
	if ( p->uz > 1e6 ) p->uz = 1;
	if ( fabs(p->ux - p->uy) > 100 ) p->ux = p->uy = p->uz = 1;
	float 		maxres = 0;
	if ( p->x > 1 ) maxres += 1/(p->ux*p->ux);
	if ( p->y > 1 ) maxres += 1/(p->uy*p->uy);
	if ( p->z > 1 ) maxres += 1/(p->uz*p->uz);
	maxres = 2/sqrt(maxres);
	if ( p->resolution < maxres ) p->resolution = maxres;	// Check if the resolution is set
	if ( p->colormodel < Bit ) p->colormodel = Gray;
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG img_check_param: Data type: %d\n", p->datatype);
		printf("DEBUG img_check_param: Checking size parameters\n");
		printf("DEBUG img_check_param: Images & channels: %ld %ld\n", p->n, p->c);
		printf("DEBUG img_check_param: Size: %ld %ld %ld\n", p->x, p->y, p->z);
		printf("DEBUG img_check_param: Page: %ld %ld %ld\n", p->px, p->py, p->pz);
		printf("DEBUG img_check_param: Units: %g %g %g\n", p->ux, p->uy, p->uz);
		printf("DEBUG img_check_param: Resolution: %g\n", p->resolution);
		printf("DEBUG img_check_param: Checking statistics:\n");
		printf("DEBUG img_check_param: min = %g  max = %g  avg = %g  std = %g\n",
				p->min, p->max, p->avg, p->std);
	}
	
	if ( !finite(p->min) || !finite(p->max) || !finite(p->avg) || !finite(p->std) )
		img_stats(p);
	if ( (p->min>=p->max) || (p->std<=0) || (p->avg<p->min) || (p->avg>p->max) )
		img_stats(p);
	if ( p->smin == 0 && p->smax == 0 ) {
		p->smin = p->min;
		p->smax = p->max;
/*		p->smin = p->avg - 5*p->std;
		p->smax = p->avg + 5*p->std;
		if ( p->smin < p->min ) p->smin = p->min;
		if ( p->smax > p->max ) p->smax = p->max; */
	}
	if ( p->scale == 0 ) p->scale = 1;
	
	if ( p->image == NULL ) {
		if ( verbose & VERB_DEBUG )
			printf("DEBUG img_check_param: Allocating memory for %ld sub-images\n", p->n);
		p->image = (Bsub_image *) balloc(p->n*sizeof(Bsub_image));
		for ( i=0; i<p->n; i++ ) p->image[i].background = p->avg;
	} else
		if ( verbose & VERB_DEBUG )
			printf("DEBUG img_check_param: Memory for %ld sub-images already allocated\n", p->n);
	
	for ( i=0; i<p->n; i++ ) {
		if ( !finite(p->image[i].background) ) p->image[i].background = 0;
		else if ( fabs((double)p->image[i].background) < 1e-30 ) p->image[i].background = 0;
		vector_normalize(3, &(p->image[i].vx));
	}
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG img_check_param: Min, max, mean, stdev: %g %g %g %g\n", 
				p->min, p->max, p->avg, p->std);
		printf("DEBUG img_check_param: Background:");
		for ( i=0; i<p->n; i++ ) printf(" %g", p->image[i].background);
		printf("\n");
		printf("DEBUG img_check_param: Checking unit cell parameters:\n");
		printf("DEBUG img_check_param: %g %g %g %g %g %g\n", p->ua, p->ub, p->uc, p->alf, p->bet, p->gam);
	}
	
	if ( p->ua < 0.001 ) p->ua = p->ux*p->x;
	if ( p->ub < 0.001 ) p->ub = p->uy*p->y;
	if ( p->uc < 0.001 ) p->uc = p->uz*p->z;
	if ( p->ua > 100*p->x*p->ux ) p->ua = p->ux*p->x;
	if ( p->ub > 100*p->y*p->uy ) p->ub = p->uy*p->y;
	if ( p->uc > 100*p->z*p->uz ) p->uc = p->uz*p->z;
	if ( p->alf > PI ) p->alf *= PI/180.0;
	if ( p->bet > PI ) p->bet *= PI/180.0;
	if ( p->gam > PI ) p->gam *= PI/180.0;
	if ( p->alf < 0.001 || p->alf > PI ) p->alf = PI/2.0;
	if ( p->bet < 0.001 || p->bet > PI ) p->bet = PI/2.0;
	if ( p->gam < 0.001 || p->gam > PI ) p->gam = PI/2.0;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_check_param: Unit cell: %g %g %g %g %g %g\n", p->ua, p->ub, p->uc, p->alf, p->bet, p->gam);

	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_check_param: reading symmetry operator table\n");

	int			j = (int) (strlen(p->label) - 1);
	
	if ( j > 0 ) {
		if ( verbose & VERB_DEBUG )
			printf("DEBUG img_check_param: label length = %d last char = %c\n", j, p->label[j]);

		 while ( j > 0 && isspace(p->label[j]) ) {
			p->label[j] = 0;
			j--;
		}
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_check_param: If a transform, converting to standard format\n");
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_check_param: Finished checking parameters\n");
	
	return(0);
}

/************************************************************************
@Function: img_set_origin
@Description:
	Sets the origins (offsets from the first voxel) in an image header.
@Algorithm:
	. 
@Arguments:
	Bimage* p			image.
	Vector3 origin	 	3-value origin vector.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int  		img_set_origin(Bimage* p, Vector3 origin)
{	
	if ( !p ) return(-1);
	
	unsigned int	i;
	
	for ( i=0; i<p->n; i++ ) {
		p->image[i].ox = origin.x;
		p->image[i].oy = origin.y;
		p->image[i].oz = origin.z;
	}
	
    return(0);
}

/************************************************************************
@Function: img_set_sampling
@Description:
	Sets the sampling values (voxel size) in an image header.
@Algorithm:
	. 
@Arguments:
	Bimage* p			image.
	Vector3 sampling 	3-value sampling vector.
	int uc_flag 		indicating resetting of unit cell parameters.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int  		img_set_sampling(Bimage* p, Vector3 sampling, int uc_flag)
{	
	if ( !p ) return(-1);
	
	if ( sampling.x ) p->ux = sampling.x;
	if ( sampling.y ) p->uy = sampling.y;
	if ( sampling.z ) p->uz = sampling.z;
	if ( uc_flag ) {
		p->ua = p->x*p->ux;
		p->ub = p->y*p->uy;
		p->uc = p->z*p->uz;
		p->alf = p->bet = p->gam = PI/2;
	}
	
    return(0);
}

/************************************************************************
@Function: img_set_isotropic_sampling
@Description:
	Sets the sampling values (voxel size) in an image header.
@Algorithm:
	. 
@Arguments:
	Bimage* p			image.
	float sampling		isotropic sampling.
	int uc_flag 		indicating resetting of unit cell parameters.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int  		img_set_isotropic_sampling(Bimage* p, float sampling, int uc_flag)
{	
	if ( !p ) return(-1);
	if ( sampling <= 0 ) return(-1);
	
	p->ux = p->uy = p->uz = sampling;
	if ( uc_flag ) {
		p->ua = p->x*p->ux;
		p->ub = p->y*p->uy;
		p->uc = p->z*p->uz;
		p->alf = p->bet = p->gam = PI/2;
	}
	
    return(0);
}

/************************************************************************
@Function: img_max_radius
@Description:
	Returns the radius of the largest sphere or circle that will fit 
	inside a three- or two-dimensional image, respectively, or the 
	midpoint of a one-dimensional image
@Algorithm:
	Assumes that the x, y, and z dimensions are the number of data points
	in those respective directions, and that the length (in pixels)
	of the respective dimensions is x, y, or z minus one.  The minimum 
	length is two times the radius of the largest sphere that will fit 
	inside the image.
@Arguments:
	Bimage*  p			image.
@Returns:
	float				radius of largest sphere that will fit in an image.
**************************************************************************/
float		img_max_radius(Bimage* p)
{
	unsigned int    min_dimension;    // minimum non-zero length of an image
	float  max_radius;

	min_dimension = (unsigned int) p->x;
	if ( (p->y > 1) && (p->y < min_dimension) )  min_dimension = (unsigned int) p->y;
	if ( (p->z > 1) && (p->z < min_dimension) )  min_dimension = (unsigned int) p->z;

	max_radius = (float) ((min_dimension - 1)/2);

	return(max_radius);
}


/************************************************************************
@Function: img_stats
@Description:
	Calculates the statistics for an image.
@Algorithm:
	. 
@Arguments:
	Bimage* p			image.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 		img_stats(Bimage *p)
{
	if ( !p ) {
		fprintf(stderr, "Error: No image in memory!\n");
		return(-1);
	}
	
	if ( p->dataflag < 1 || !p->data ) {
//		fprintf(stderr, "Error: No data for image %s in memory!\n", p->filename);
		return(-1);
	}
	
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    complex_short* 	    csdata = (complex_short *) p->data;
    complex_int* 	    cidata = (complex_int *) p->data;
    complex_float*  	cfdata = (complex_float *) p->data;
    polar*  	    	pdata = (polar *) p->data;
    
	char				bit;
    unsigned long		i, j, n, datasize = (unsigned long) p->x*p->y*p->z*p->c;
    unsigned long		fomsize = (unsigned long) p->x*p->y*p->z*p->n;
	double				amp, tsum, tssum, avg, std;
	double*				min = (double *) balloc(p->n*sizeof(double));
	double*				max = (double *) balloc(p->n*sizeof(double));
	double*				sum = (double *) balloc(p->n*sizeof(double));
	double*				ssum = (double *) balloc(p->n*sizeof(double));
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG img_stats: Data size:  %ld x %ld x %ld x %ld x %ld = %ld\n",
				p->x, p->y, p->z, p->c, p->n, datasize);
		printf("DEBUG img_stats: datatype = %d\n", p->datatype);
		printf("DEBUG img_stats: min = %g,   max = %g\n", p->min, p->max);
		printf("DEBUG img_stats: ave = %g,   std = %g\n", p->avg, p->std);
		printf("DEBUG img_stats: fomflag = %d,   fom = %p\n", p->fomflag, p->fom);
	}

	if ( p->fomflag && p->fom )
		for ( i=0, p->fommax=p->fom[0]; i<fomsize; i++ )
			if ( p->fommax < p->fom[i] ) p->fommax = p->fom[i];

	for ( n=0; n<p->n; n++ ) {
		min[n] = 1e37;
		max[n] = -1e37;
	}

    switch ( p->datatype ) {
    	case UChar:
			if ( p->colormodel == Bit ) {
				for ( n=i=0; n<p->n; n++ ) {
					for ( j=0; j<datasize; i++, j++ ) {
//						bit = 0x80 & ( udata[i/8] << i%8 );
						bit = 0;
						if ( 0x80 & ( udata[i/8] << i%8 ) ) bit = 1;
						sum[n] += bit;
						ssum[n] += bit*bit;
					}
					min[n] = 0;
					max[n] = 1;
				}
				p->min = 0;
				p->max = 1;
			} else {
				for ( n=i=0; n<p->n; n++ ) {
					for ( j=0; j<datasize; i++, j++ ) {
						if ( udata[i] < min[n] ) min[n] = udata[i];
						if ( udata[i] > max[n] ) max[n] = udata[i];
						sum[n] += udata[i];
						ssum[n] += (double)udata[i]*udata[i];
					}
				}
    	    }
    	    break;
    	case SChar:
			for ( n=i=0; n<p->n; n++ ) {
				for ( j=0; j<datasize; i++, j++ ) {
					if ( cdata[i] < min[n] ) min[n] = cdata[i];
					if ( cdata[i] > max[n] ) max[n] = cdata[i];
					sum[n] += cdata[i];
					ssum[n] += (double)cdata[i]*cdata[i];
				}
    	    }
    	    break;
    	case UShort:
			for ( n=i=0; n<p->n; n++ ) {
				for ( j=0; j<datasize; i++, j++ ) {
					if ( usdata[i] < min[n] ) min[n] = usdata[i];
					if ( usdata[i] > max[n] ) max[n] = usdata[i];
					sum[n] += usdata[i];
					ssum[n] += (double)usdata[i]*usdata[i];
				}
			}
    	    break;
    	case Short:
			for ( n=i=0; n<p->n; n++ ) {
				for ( j=0; j<datasize; i++, j++ ) {
					if ( sdata[i] < min[n] ) min[n] = sdata[i];
					if ( sdata[i] > max[n] ) max[n] = sdata[i];
					sum[n] += sdata[i];
					ssum[n] += (double)sdata[i]*sdata[i];
				}
    	    }
    	    break;
		case Int:
			for ( n=i=0; n<p->n; n++ ) {
				for ( j=0; j<datasize; i++, j++ ) {
					if ( idata[i] < min[n] ) min[n] = idata[i];
					if ( idata[i] > max[n] ) max[n] = idata[i];
					sum[n] += idata[i];
					ssum[n] += (double)idata[i]*idata[i];
				}
    	    }
			break;
    	case Float:
			for ( n=i=0; n<p->n; n++ ) {
				for ( j=0; j<datasize; i++, j++ ) {
					if ( fdata[i] < min[n] ) min[n] = fdata[i];
					if ( fdata[i] > max[n] ) max[n] = fdata[i];
					sum[n] += fdata[i];
					ssum[n] += (double)fdata[i]*fdata[i];
				}
    	    }
    	    break;
    	case ComplexShort:
			for ( n=i=0; n<p->n; n++ ) {
				for ( j=0; j<datasize; i++, j++ ) {
					amp = (double)csdata[i].re*csdata[i].re + (double)csdata[i].im*csdata[i].im;
					amp = sqrt(amp);
					if ( amp < min[n] ) min[n] = amp;
					if ( amp > max[n] ) max[n] = amp;
					sum[n] += amp;
					ssum[n] += amp*amp;
				}
    	    }
    	    break;
		case ComplexInt:
			for ( n=i=0; n<p->n; n++ ) {
				for ( j=0; j<datasize; i++, j++ ) {
					amp = (double)cidata[i].re*cidata[i].re + (double)cidata[i].im*cidata[i].im;
					amp = sqrt(amp);
					if ( amp < min[n] ) min[n] = amp;
					if ( amp > max[n] ) max[n] = amp;
					sum[n] += amp;
					ssum[n] += amp*amp;
				}
    	    }
    	    break;
    	case ComplexFloat:
			for ( n=i=0; n<p->n; n++ ) {
				for ( j=0; j<datasize; i++, j++ ) {
					amp = sqrt(cfdata[i].re*cfdata[i].re + cfdata[i].im*cfdata[i].im);
					if ( amp < min[n] ) min[n] = amp;
					if ( amp > max[n] ) max[n] = amp;
					sum[n] += amp;
					ssum[n] += amp*amp;
				}
    	    }
    	    break;
    	case Polar:
			for ( n=i=0; n<p->n; n++ ) {
				for ( j=0; j<datasize; i++, j++ ) {
					if ( pdata[i].amp < min[n] ) min[n] = pdata[i].amp;
					if ( pdata[i].amp > max[n] ) max[n] = pdata[i].amp;
					sum[n] += pdata[i].amp;
					ssum[n] += (double)pdata[i].amp*pdata[i].amp;
				}
    	    }
    	    break;
    	default:
			fprintf(stderr, "Warning: Data type %d not supported!", p->datatype);
			p->min = p->max = 0;
			break;
    }
	
	p->min = min[0];
	p->max = max[0];
	tsum = tssum = 0;
	for ( n=0; n<p->n; n++ ) {
		if ( p->min > min[n] ) p->min = min[n];
		if ( p->max < max[n] ) p->max = max[n];
		tsum += sum[n];
		tssum += ssum[n];
	}

    p->avg = tsum/(p->n*datasize);
    p->std = tssum/(p->n*datasize) - p->avg*p->avg;
    if ( p->std > 0 ) p->std = sqrt(p->std);
	else p->std = 0;
	p->smin = p->min;
	p->smax = p->max;

	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_stats: tsum = %g tssum = %g\n", tsum, tssum);
	
	if ( verbose & VERB_STATS ) {
	    printf("Data size:                      %ld x %ld x %ld x %ld x %ld = %ld\n",
				p->x, p->y, p->z, p->c, p->n, p->n*datasize);
	    printf("Min, max, avg, std:             %g %g %g %g\n", 
				p->min, p->max, p->avg, p->std);
	    if ( p->fomflag ) printf("FOM maximum:                    %g\n", p->fommax);
		printf("\nImage\tMin\tMax\tAvg\tStd\n");
		for ( n=0; n<p->n; n++ ) {
			avg = sum[n]/datasize;
			std = ssum[n]/datasize - avg*avg;
			if ( std > 0 ) std = sqrt(std);
			else std = 0;
			printf("%ld\t%g\t%g\t%g\t%g\n", n+1, min[n], max[n], avg, std);
		}
		printf("\n");
	}
	
	bfree(min, p->n*sizeof(double));
	bfree(max, p->n*sizeof(double));
	bfree(sum, p->n*sizeof(double));
	bfree(ssum, p->n*sizeof(double));
	
	return(0);
}

/************************************************************************
@Function: img_stats_within_radii
@Description:
	Calculates the statistics for an image with given radii.
@Algorithm:
	If a voxel lies within the specified radii, it is included in
	the statistical calculations.
@Arguments:
	Bimage* p			image.
	float rad_min    	minimum radius (pixel units).
	float rad_max    	maximum radius (pixel units).
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 		img_stats_within_radii(Bimage *p, float rad_min, float rad_max)
{
	if ( !p ) {
		fprintf(stderr, "Error: No image in memory!\n");
		return(-1);
	}
	
	if ( p->dataflag < 1 || !p->data ) {
//		fprintf(stderr, "Error: No data for image %s in memory!\n", p->filename);
		return(-1);
	}
	
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    polar*  	    	pdata = (polar *) p->data;
    
    unsigned long		i, x, y, z, n;
	unsigned long		datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	float				dx, dy, dz, d2, rmax2 = rad_max*rad_max, rmin2 = rad_min*rad_min;
	double				value, sum = 0, ssum = 0, w = 0;
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG img_stats: Data size:  %ld x %ld x %ld x %ld x %ld = %ld\n",
				p->x, p->y, p->z, p->c, p->n, datasize);
		printf("DEBUG img_stats: datatype = %d\n", p->datatype);
		printf("DEBUG img_stats: min = %g,   max = %g\n", p->min, p->max);
		printf("DEBUG img_stats: ave = %g,   std = %g\n", p->avg, p->std);
		printf("DEBUG img_stats: fomflag = %d,   fom = %p\n", p->fomflag, p->fom);
	}

	if ( p->fomflag )
		for ( i=0; i<datasize; i++ )
			if ( p->fommax < p->fom[i] ) p->fommax = p->fom[i];

	if ( p->datatype >= ComplexShort && p->datatype <= ComplexFloat ) datasize *= 2;
	
    p->min = 1e37;
    p->max = -1e37;
	for ( n=0; n<p->n; n++ ) {
		for ( z=0; z<p->z; z++ ) {
			dz = z - p->image[n].oz;
			dz *= dz;
			if ( dz <= rmax2 ) for ( y=0; y<p->y; y++ ) {
				dy = y - p->image[n].oy;
				dy *= dy;
				if ( dy <= rmax2 ) for ( x=0; x<p->x; x++ ) {
					dx = x - p->image[n].ox;
					dx *= dx;
					d2 = dx + dy + dz;
					if ( d2 <= rmax2 && d2 >= rmin2 ) {
						value = 0;
						i = ((n*p->z + z)*p->y + y)*p->x + x;
						switch ( p->datatype ) {
							case UChar:
								if ( p->colormodel == Bit ) {
									if ( 0x80 & ( udata[i/8] << i%8 ) ) value = 1;
								} else {
									value = udata[i];
								}
								break;
							case SChar:
								value = cdata[i];
								break;
							case Short:
							case ComplexShort:
								value = sdata[i];
								break;
							case UShort:
								value = usdata[i];
								break;
							case Int:
							case ComplexInt:
								value = idata[i];
								break;
							case Float:
							case ComplexFloat:
								value = fdata[i];
								break;
							case Polar:
								value = pdata[i].amp;
								break;
							default:
								value = 0;
								break;
						}
						sum += value;
						ssum += value*value;
						w += 1;
						if ( p->min > value ) p->min = value;
						if ( p->max < value ) p->max = value;
					}
				}
			}
		}
    }

    p->avg = sum/w;
    p->std = ssum/w - p->avg*p->avg;
    if ( p->std > 0 ) p->std = sqrt(p->std);
	else p->std = 0;
	p->smin = p->min;
	p->smax = p->max;

	if ( verbose & VERB_STATS ) {
	    printf("Data size:                      %ld x %ld x %ld x %ld x %ld = %lg\n",
				p->x, p->y, p->z, p->c, p->n, w);
	    printf("Min, max, avg, std:             %g %g %g %g\n", 
				p->min, p->max, p->avg, p->std);
	}
	
	return(0);
}

/************************************************************************
@Function: img_convert_fourier
@Description:
	Converts Fourier transform types.
@Algorithm:
	Fourier transform classification:
	0=NoTransform:	No transform: Just a complex data set
	1=Standard:		Standard transform with origin = (0,0,0)
	2=Centered:		Centered transform with origin = (nx/2,ny/2,nz/2)
						(Imagic)
	3=Hermitian:	Hermitian transform with origin = (0,0,0) and size (nx/2+1,ny,nz)
						(EM)
	4=CentHerm:		Centered hermitian transform with origin = (0,ny/2,nz/2)
					and size (nx/2+1,ny,nz)
						(MRC)
	Assumption: 	The correct dimensions for a standard transform is stored
					in the x, y, and z fields.
@Arguments:
	Bimage* p				the image.
	Transform newtransform	new transform type.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int			img_convert_fourier(Bimage* p, Transform newtransform)
{
	if ( !p ) return(-1);
	
	if ( p->dataflag < 1 ) return(-1);
	
	if ( p->datatype < ComplexShort ) {
		p->transform = NoTransform;
		return(0);
	}
	
	if ( p->transform == NoTransform )
		p->transform = Standard;
	
	if ( newtransform == p->transform )
		return(0);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_convert_fourier: ");
	
	if ( verbose & VERB_PROCESS ) {
		printf("Converting Fourier transform:   ");
		switch ( p->transform ) {
			case Standard: printf("Standard -> "); break;
			case Centered: printf("Centered -> "); break;
			case Hermitian: printf("Hermitian -> "); break;
			case CentHerm: printf("Centered hermitian -> "); break;
			default: break;
		}
		switch ( newtransform ) {
			case Standard: printf("Standard"); break;
			case Centered: printf("Centered"); break;
			case Hermitian: printf("Hermitian"); break;
			case CentHerm: printf("Centered hermitian"); break;
			default: break;
		}
		printf("\n\n");
	}
		
	unsigned int	hermx = (unsigned int) (p->x/2 + 1);
	unsigned int	friedel = 0;
	unsigned int	oldx = (unsigned int) (p->x), newx = (unsigned int) (p->x);
	if ( p->transform == Hermitian || p->transform == CentHerm ) {
		oldx = hermx;
		if ( newtransform == Standard || newtransform == Centered )
			friedel = 1;		// Flag => fill in friedel pairs
	}
	if ( newtransform == Hermitian || newtransform == CentHerm )
		newx = hermx;
	
	unsigned int	xo1, yo1, zo1, xo2, yo2, zo2;	// Origins
	switch ( p->transform ) {
		case Centered:
			xo1 = (unsigned int) (p->x/2);
			yo1 = (unsigned int) (p->y/2);
			zo1 = (unsigned int) (p->z/2);
			break;
		case CentHerm:
			xo1 = 0;
			yo1 = (unsigned int) (p->y/2);
			zo1 = (unsigned int) (p->z/2);
			break;
		default:
			xo1 = yo1 = zo1 = 0;
			break;
	}
	
	switch ( newtransform ) {
		case Centered:
			xo2 = (unsigned int) (p->x/2);
			yo2 = (unsigned int) (p->y/2);
			zo2 = (unsigned int) (p->z/2);
			break;
		case CentHerm:
			xo2 = 0;
			yo2 = (unsigned int) (p->y/2);
			zo2 = (unsigned int) (p->z/2);
			break;
		default:
			xo2 = yo2 = zo2 = 0;
			break;
	}
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG img_convert_fourier: hermx = %d, oldx = %d, newx = %d\n", 
				hermx, oldx, newx);
		printf("DEBUG img_convert_fourier: input origins = %d %d %d\n",
				xo1, yo1, zo1);
		printf("DEBUG img_convert_fourier: output origins = %d %d %d\n",
				xo2, yo2, zo2);
	}
	
    complex_short*	sdata = (complex_short *) p->data;
    complex_int*	idata = (complex_int *) p->data;
    complex_float*	fdata = (complex_float *) p->data;
    polar*			pdata = (polar *) p->data;
	
	unsigned long 	x, y, z, i, j, n;
	long 			xf, yf, zf, xi, yi, zi;
    unsigned long   datasize = (unsigned long) p->n*newx*p->y*p->z*p->c;
	int 			phisgn = 1;
	if ( p->datatype == ComplexShort ) datasize *= sizeof(complex_short);
	else datasize *= sizeof(complex_float);
    complex_float*	newfdata = (complex_float *) balloc(datasize*sizeof(char));
    complex_int*	newidata = (complex_int *) newfdata;
    complex_short*	newsdata = (complex_short *) newfdata;
    polar*			newpdata = (polar *) newfdata;
    
	int 			z_even = 0, y_even = 0;
	if ( 2*(p->z/2) == p->z ) z_even = 1;
	if ( 2*(p->y/2) == p->y ) y_even = 1;
	
	for ( n=0; n<p->n; n++ ) {
		for ( z=0; z<p->z; z++ ) {
			zi = z + zo1 - zo2;
			if ( zi < 0 ) zi += p->z;
			if ( zi >= (int)p->z ) zi -= p->z;
			if ( z_even )
				zf = (zi)? p->z-zi: zi;
			else
				zf = 2*zo1 - zi;
			if ( zf < 0 ) zf += p->z;
			for ( y=0; y<p->y; y++ ) {
    	    	yi = y + yo1 - yo2;
	    		if ( yi < 0 ) yi += p->y;
	    		if ( yi >= (int)p->y ) yi -= p->y;
				if ( y_even )
					yf = (yi)? p->y-yi: yi;
				else
					yf = 2*yo1 - yi;
				if ( yf < 0 ) yf += p->y;
				for ( x=0; x<newx; x++ ) { //printf("%d %d %d\n", x, y, z);
					xi = x + xo1 - xo2;
	    			if ( xi < 0 ) xi += p->x;
	    			if ( xi >= (int)p->x ) xi -= p->x;
					xf = (xi)? p->x-xi: xi;
					if ( !friedel || ( xi < (int)hermx ) ) {
						i = ((n*p->z + zi)*p->y + yi)*oldx + xi;
						phisgn = 1;
					} else {
						i = ((n*p->z + zf)*p->y + yf)*oldx + xf;
						phisgn = -1;
					}
					j = ((n*p->z + z)*p->y + y)*newx + x;
					if ( p->datatype == ComplexShort ) {
						newsdata[j] = sdata[i];
						newsdata[j].im *= phisgn;
					} else if ( p->datatype == ComplexInt ) {
						newidata[j] = idata[i];
						newidata[j].im *= phisgn;
					} else if ( p->datatype == ComplexFloat ) {
						newfdata[j] = fdata[i];
						newfdata[j].im *= phisgn;
					} else {
						newpdata[j] = pdata[i];
						newpdata[j].phi *= phisgn;
					}
				}
			}
		}
	}
    
	p->data = (char *) newfdata;
    datasize = (long) p->n*oldx*p->y*p->z*p->c;
	if ( p->datatype == ComplexShort )
		bfree(fdata, datasize*sizeof(complex_short));
	else
		bfree(fdata, datasize*sizeof(complex_float));
	
	p->transform = newtransform;
	
	return(0);
}

/************************************************************************
@Function: img_info
@Description:
	Prints out header information for an image.
@Algorithm:
	.
@Arguments:
	Bimage* p			the image.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int		img_info(Bimage* p)
{
	if ( !p ) return(-1);
	
	printf("Data offset:                    %ld\n", p->offset);
	printf("Number of images:               %ld\n", p->n);
	printf("Dimensions:                     %ld %ld %ld voxels\n", p->x, p->y, p->z);
	printf("Page dimensions:                %ld %ld %ld voxels\n", p->px, p->py, p->pz);
	printf("Channels:                       %ld\n", p->c);
	char*		datatype_string = get_string_from_datatype(p->datatype);
	printf("Data type:                      %s (size = %ld)\n", 
			datatype_string, gettypesize(p->datatype));
	bfree_string(datatype_string);
	if ( p->transform ) {
		printf("Transform:                      ");
		switch ( p->transform ) {
			case Standard: printf("Standard\n"); break;
			case Centered: printf("Centered\n"); break;
			case Hermitian: printf("Hermitian\n"); break;
			case CentHerm: printf("Centered hermitian\n"); break;
			default: break;
		}
	}
	printf("Color model:                    ");
	switch ( p->colormodel ) {
		case Bit: printf("Bitmap\n"); break;
		case Gray: printf("Gray scale\n"); break;
		case RGB: printf("RGB\n"); break;
		case RGBA: printf("RGBA\n"); break;
		case CMYK: printf("CMYK\n"); break;
		case Index: printf("Indexed colour map\n"); break;
		default: printf("unknown\n"); break;
	}
	printf("Voxel units/sampling:           %g %g %g A/voxel\n", p->ux, p->uy, p->uz);
	printf("Min, max, ave, std:             %g %g %g %g\n", 
			p->min, p->max, p->avg, p->std);
//	printf("Display min & max:              %g %g\n\n", p->smin, p->smax);
	
    printf("Text label:                     (length = %ld)\n%s\n", 
			strlen(p->label), p->label);
	
	UnitCell		unit_cell = {p->ua, p->ub, p->uc, p->alf, p->bet, p->gam};
	float*          frac_mat = NULL;
	if ( p->spacegroup > 0 ) {
		frac_mat = calc_skew_matrix(unit_cell);
		printf("Unit cell dimensions:           %g %g %g A\n", 
				p->ua, p->ub, p->uc);
		printf("Unit cell angles:               %g %g %g degrees\n", 
				p->alf*180/PI, p->bet*180/PI, p->gam*180/PI);
		printf("Fractionalization matrix:       %10.6f %10.6f %10.6f\n", 
				frac_mat[0], frac_mat[1], frac_mat[2]);
		printf("                                %10.6f %10.6f %10.6f\n", 
				frac_mat[3], frac_mat[4], frac_mat[5]);
		printf("                                %10.6f %10.6f %10.6f\n", 
				frac_mat[6], frac_mat[7], frac_mat[8]);
		printf("Fractional translation:         %10.6f %10.6f %10.6f\n", 
				frac_mat[9], frac_mat[10], frac_mat[11]);
	    printf("Space group:                    %d\n", p->spacegroup);
		bfree(frac_mat, 12*sizeof(float));
	}
	
	if ( p->spacelabel ) {
	    if ( strlen(p->spacelabel) )
			printf("Symmetry label:                 %s\n", p->spacelabel);
	}
	
	if ( p->fomflag && p->fom ) {
	    printf("Maximum FOM:                    %g\n", p->fommax);
	}
	
	printf("\n");
	
	return(0);
}

/************************************************************************
@Function: sub_img_info
@Description:
	Prints out header information for all sub-images.
@Algorithm:
	.
@Arguments:
	Bimage* p			the image.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int		sub_img_info(Bimage* p)
{
	if ( !p ) return(-1);
	
	unsigned int	i;
	
	printf("Image\tBackground\tOrigin\t\t\tView\n");
	for ( i=0; i<p->n; i++ )
		printf("%7d\t%12g\t%6g %6g %6g\t%7.4f %7.4f %7.4f %7.2f\n", i+1, p->image[i].background, 
				p->image[i].ox, p->image[i].oy, p->image[i].oz,
				p->image[i].vx, p->image[i].vy, p->image[i].vz, p->image[i].angle*180.0/PI);
	
	printf("\n");
	
	return(0);
}
