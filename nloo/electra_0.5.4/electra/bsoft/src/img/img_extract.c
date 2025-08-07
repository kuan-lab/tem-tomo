/*
	img_extract.c
	Library routines to extract parts of images
	Author: Bernard Heymann
	Created: 19990321 	Modified: 20050112
*/

#include "rwimg.h"
#include "img_edit.h"
#include "img_background.h"
#include "img_datatypes.h"
#include "img_util.h"
#include "spline.h"
#include "matrix.h"
//#include "colour.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated

// Internal function prototype
Bimage*		img_extract_bitmap(Bimage* p, int img_number, Vector3 origin, VectorInt3 size);

/************************************************************************
@Function: img_extract
@Description:
	Extracts part of an image into a new image.
@Algorithm:
	A single image (only one sub-image) is extracted from a given sub-image
	in an image structure, starting at a specified origin and with a 
	specified size.
	The old data is not affected.
	Statistics for the extracted image are calculated.
@Arguments:
	Bimage* p			image (not altered).
	int img_number		image number from which to extract, -1 indicates all images.
	Vector3 origin		three-value origin vector = origin of extracted image.
	VectorInt3 size		three-vector size of extracted image.
@Returns:
	Bimage* p			extracted image.
**************************************************************************/
Bimage*		img_extract(Bimage* p, int img_number, Vector3 origin, VectorInt3 size)
{
	if ( !p->dataflag ) return(NULL);
	
	if ( size.x*size.y*size.z < 1 ) {
		error_show("img_extract", __FILE__, __LINE__);
		fprintf(stderr, "Error: Size of subimage is less than one!\n");
		return(NULL);
	}
	
	if ( p->colormodel == Bit ) {
		Bimage*		pex = img_extract_bitmap(p, img_number, origin, size);
		return(pex);
	}
	
	unsigned long   i, j;
	unsigned long	x, y, z, n, new_nimg;
	long			oldn, oldx, oldy, oldz;
	
	if ( img_number < 0 ) {
		new_nimg = p->n;
		img_number = 0;
	} else {
		new_nimg = 1;
		if ( img_number >= p->n ) img_number = (int) p->n - 1;
	}
	
	Bimage*		pex = copy_img_header(p, (int) new_nimg);
	
	pex->x = pex->px = size.x;
	pex->y = pex->py = size.y;
	pex->z = pex->pz = size.z;
	pex->c = p->c;
	
	unsigned long	elementsize = pex->c*gettypesize(pex->datatype);
	unsigned long   datasize = (unsigned long)pex->n*pex->x*pex->y*pex->z;
	char*			data = (char *) p->data;
	char*			newdata = (char *) balloc(datasize*elementsize);
	pex->dataflag = 1;
	pex->data = (char *) newdata;

	if ( verbose & VERB_PROCESS ) {
		if ( new_nimg == 1 )
			printf("Extracting subimage:            %d\n", img_number);
		else
			printf("Extracting subimages:			%d - %ld\n", 0, pex->n - 1);
		printf("Size:                           %ld %ld %ld\n", pex->x, pex->y, pex->z);
		printf("Origin:                         %g %g %g\n", origin.x, origin.y, origin.z);
		printf("\n");
	} else if ( verbose & VERB_LABEL ) {
		if ( new_nimg == 1 )
			printf("Extracting subimage:            %d\n\n", img_number);
		else
			printf("Extracting subimages:			%d - %ld\n\n", 0, pex->n - 1);
	}
	
	unsigned char*  fptr = set_value_datatype(p->avg, p->datatype, (int) p->c);
	for ( i=0; i<datasize; i++ ) memcpy(newdata+i*elementsize, fptr, elementsize);
	bfree(fptr, elementsize);
	
	for ( n=0; n<new_nimg; n++ ) {
		pex->image[n].ox = origin.x;
		pex->image[n].oy = origin.y;
		pex->image[n].oz = origin.z;
		oldn = n + img_number;
		for ( z=0; z<pex->z; z++ ) {
			oldz = (long) (z + pex->image[n].oz);
			if ( oldz >= 0 && oldz < p->z ) {
				for ( y=0; y<pex->y; y++ ) {
					oldy = (long) (y + pex->image[n].oy);
					if ( oldy >= 0 && oldy < p->y ) {
						for ( x=0; x<pex->x; x++ ) {
							oldx = (long) (x + pex->image[n].ox);
							if ( oldx >= 0 && oldx < p->x ) {
								i = (((n*pex->z + z)*pex->y + y)*pex->x + x)*elementsize;
								j = (((oldn*p->z + oldz)*p->y + oldy)*p->x + oldx)*elementsize;
								memcpy(newdata+i, data+j, elementsize);
							}
						}
					}
				}
			}
		}
	}
	
	img_stats(pex);
	
	return(pex);
}

/*
@Function: img_extract_bitmap
@Description:
	Extracts part of a bitmap image into a new bitmap image.
@Algorithm:
	A single image (only one sub-image) is extracted from a given sub-image
	in an image structure, starting at a specified origin and with a 
	specified size.
	The old data is not affected.
	Statistics for the extracted image are calculated.
@Arguments:
	Bimage* p			image (not altered).
	int img_number		image number from which to extract, -1 indicates all images.
	Vector3 origin		three-value origin vector = origin of extracted image.
	VectorInt3 size		three-vector size of extracted image.
@Returns:
	Bimage* p		extracted image.
**************************************************************************/
Bimage*		img_extract_bitmap(Bimage* p, int img_number, Vector3 origin, VectorInt3 size)
{
	if ( !p->dataflag ) return(NULL);
	
	if ( p->colormodel != Bit ) {
		error_show("img_extract", __FILE__, __LINE__);
		fprintf(stderr, "Error: Cannot extract a bitmap from a non-bitmap image!\n");
		return(NULL);
	}
	
	if ( size.x*size.y*size.z < 1 ) {
		error_show("img_extract_bitmap", __FILE__, __LINE__);
		fprintf(stderr, "Error: Size of subimage is less than one!\n");
		return(NULL);
	}
	
	unsigned long	i, j, y, z, n, new_nimg;
	long			oldn, oldx, oldy, oldz;
	
	if ( img_number < 0 ) {
		new_nimg = p->n;
		img_number = 0;
	} else {
		new_nimg = 1;
		if ( img_number >= p->n ) img_number = (int) p->n - 1;
	}
	
	Bimage*		pex = copy_img_header(p, (int) new_nimg);
	
	pex->x = pex->px = size.x;
	pex->y = pex->py = size.y;
	pex->z = pex->pz = size.z;
	
	if ( verbose & VERB_PROCESS ) {
		printf("Extracting a bitmap subimage:\n");
		if ( new_nimg == 1 )
			printf("Image number:                   %d\n", img_number);
		else
			printf("Image numbers:                  %d - %ld\n", 0, pex->n - 1);
		printf("Size:                           %ld %ld %ld\n", pex->x, pex->y, pex->z);
		printf("Origin:                         %g %g %g\n", origin.x, origin.y, origin.z);
		printf("\n");
	} else if ( verbose & VERB_LABEL )
		printf("Extracting a bitmap subimage\n\n");
	
	unsigned long	elementsize = pex->c*gettypesize(pex->datatype);
	unsigned long	linesize = (pex->x*elementsize)/8;
	char*			data = (char *) p->data;
	char*			newdata = (char *) balloc(pex->n*pex->z*pex->y*linesize);
	pex->data = (char *) newdata;
	pex->dataflag = 1;

	for ( n=0; n<new_nimg; n++ ) {
		pex->image[n].ox = origin.x;
		pex->image[n].oy = origin.y;
		pex->image[n].oz = origin.z;
		oldn = n + img_number;
		for ( z=0; z<pex->z; z++ ) {
			oldz = (int) (z + pex->image[n].oz);
			for ( y=0; y<pex->y; y++ ) {
				oldy = (int) (y + pex->image[n].oy);
				oldx = (int) (pex->image[n].ox);
				i = ((n*pex->z + z)*pex->y + y)*linesize;
				j = ((((oldn*p->z + oldz)*p->y + oldy)*p->x + oldx)*elementsize)/8;
				memcpy(newdata+i, data+j, linesize);
			}
		}
	}

	img_stats(pex);
	
	return(pex);
}

/************************************************************************
@Function: img_extract_tiles
@Description:
	Extracts a set of tiles from an image into a new image.
@Algorithm:
	A set of tiles of specified size are extracted from a given sub-image 
	in an image structure, starting from the origin pixel and generating
	as many tiles as would fit into the image size given. If the image size
	given is zero, then the whole image is used. The origins of the tiles
	are inserted into the sub-image origin fields.
	The old data is not affected.
@Arguments:
	Bimage* p				image (not altered).
	int img_number			image number from which to extract, -1 indicates all images.
	VectorInt3 origin 		3-vector origin for first tile to be extracted.
	VectorInt3 size			3-vector size of part of image to be extracted (NULL = whole image).
	VectorInt3 tile_size	3-vector size of extracted image.
@Returns:
	Bimage* p				extracted set of sub-images.
**************************************************************************/
/*Bimage*		img_extract_tiles(Bimage* p, int img_number, VectorInt3 origin,
				VectorInt3 size, VectorInt3 tile_size)
{
	if ( !p->dataflag ) return(NULL);
	
	if ( vectorint3_length(tile_size) < 1 ) {
		error_show("img_extract_tiles", __FILE__, __LINE__);
		fprintf(stderr, "Error: Size of subimage is less than one!\n");
		return(NULL);
	}
	
	unsigned long   i, j;
	unsigned long	x, y, z, n, new_nimg, nx, ny, nz;
	long			oldx, oldy, oldz;
	
	if ( vectorint3_length(size) < 2 ) {
		size.x = p->x;
		size.y = p->y;
		size.z = p->z;
	}
	
	if ( origin.x > p->x ) origin.x = 0;
	if ( origin.y > p->y ) origin.y = 0;
	if ( origin.z > p->z ) origin.z = 0;
	
	if ( size.x > p->x - origin.x ) size.x = p->x - origin.x;
	if ( size.y > p->y - origin.y ) size.y = p->y - origin.y;
	if ( size.z > p->z - origin.z ) size.z = p->z - origin.z;
	
	if ( tile_size.x > size.x ) tile_size.x = size.x;
	if ( tile_size.y > size.y ) tile_size.y = size.y;
	if ( tile_size.z > size.z ) tile_size.z = size.z;
	
	if ( img_number < 0 ) img_number = 0;
	
	// Count the number of tiles to be extracted
	new_nimg = 0;
	for ( z=0; z<=size.z-tile_size.z; z+=tile_size.z )
		for ( y=0; y<=size.y-tile_size.y; y+=tile_size.y )
			for ( x=0; x<=size.x-tile_size.x; x+=tile_size.x )
				new_nimg++;

	Bimage*		pex = copy_img_header(p, new_nimg);
	
	pex->x = pex->px = tile_size.x;
	pex->y = pex->py = tile_size.y;
	pex->z = pex->pz = tile_size.z;
	
	if ( verbose & VERB_PROCESS ) {
		printf("Extracting %ld tiles:\n", pex->n);
		printf("Origin:                         %d %d %d\n", 
				origin.x, origin.y, origin.z);
		printf("Tile size:                      %ld %ld %ld\n", 
				pex->x, pex->y, pex->z);
		printf("\n");
	} else if ( verbose & VERB_LABEL )
		printf("Extracting tiles\n\n");
	
	unsigned long	elementsize = pex->c*gettypesize(pex->datatype);
	char*			data = (char *) p->data;
	char*			newdata = (char *) balloc(pex->x*pex->y*pex->z*
							pex->n*elementsize);
	pex->data = (char *) newdata;
	pex->dataflag = 1;
	
	n = 0;
	for ( z=0; z<=size.z-pex->z; z+=pex->z ) {
		for ( y=0; y<=size.y-pex->y; y+=pex->y ) {
			for ( x=0; x<=size.x-pex->x; x+=pex->x ) {
				for ( nz=0; nz<pex->z; nz++ ) {
					oldz = z + nz + origin.z;
					for ( ny=0; ny<pex->y; ny++ ) {
						oldy = y + ny + origin.y;
						for ( nx=0; nx<pex->x; nx++ ) {
							oldx = x + nx + origin.x;
							i = ((n*pex->z+nz)*pex->y + ny)*pex->x + nx;
							j = ((img_number*p->z + oldz)*p->y + oldy)*p->x + oldx;
							memcpy(newdata+i*elementsize, data+j*elementsize, elementsize);
						}
					}
				}
				pex->image[n].ox = x + origin.x;
				pex->image[n].oy = y + origin.y;
				pex->image[n].oz = z + origin.z;
				n++;
			}
		}
	}
	
	img_stats(pex);
	
	return(pex);
}
*/
/************************************************************************
@Function: img_extract_tiles_at_coords
@Description:
	Extracts a set of tiles at specified positions from an image into a new image.
@Algorithm:
	A set of tiles of specified size are extracted from a given sub-image 
	in an image structure, using given coordinates for the tiles, which have
	to fit into the image size. The origins of the tiles are inserted into 
	the sub-image origin fields.
	The old data is not affected.
@Arguments:
	Bimage* p				image (not altered).
	int img_number			image number from which to extract, -1 indicates all images.
	int ntiles 				number of tiles to extract.
	VectorInt3* coords 		coordinates for the tile origins.
	VectorInt3 tile_size	3-vector size of extracted image.
@Returns:
	Bimage* p				extracted set of sub-images.
**************************************************************************/
/*Bimage*		img_extract_tiles_at_coords(Bimage* p, int img_number,
				int ntiles, VectorInt3* coords, VectorInt3 tile_size)
{
	if ( !p->dataflag ) return(NULL);
	
	if ( vectorint3_length(tile_size) < 1 ) {
		error_show("img_extract_tiles_at_coords", __FILE__, __LINE__);
		fprintf(stderr, "Error: Size of subimage is less than one!\n");
		return(NULL);
	}
	
	unsigned long	i, j, x, y, z, n;
	long			oldx, oldy, oldz;
	
	if ( tile_size.x > p->x ) tile_size.x = p->x;
	if ( tile_size.y > p->y ) tile_size.y = p->y;
	if ( tile_size.z > p->z ) tile_size.z = p->z;
	
	if ( img_number < 0 ) img_number = 0;
	
	Bimage*		pex = copy_img_header(p, ntiles);
	
	pex->x = pex->px = tile_size.x;
	pex->y = pex->py = tile_size.y;
	pex->z = pex->pz = tile_size.z;
	
	if ( verbose & VERB_PROCESS ) {
		printf("Extracting %ld tiles using origin coordinates:\n", pex->n);
		printf("Tile size:                      %ld %ld %ld\n", 
				pex->x, pex->y, pex->z);
		printf("\n");
	} else if ( verbose & VERB_LABEL )
		printf("Extracting tiles using origin coordinates\n\n");
	
	unsigned long	elementsize = pex->c*gettypesize(pex->datatype);
	unsigned long   datasize = (unsigned long)pex->n*pex->x*pex->y*pex->z*p->c;
	char*			data = (char *) p->data;
	char*			newdata = (char *) balloc(pex->x*pex->y*pex->z*
							pex->n*elementsize);
	pex->data = (char *) newdata;
	pex->dataflag = 1;
	
	unsigned char*  fptr = set_value_datatype(p->avg, p->datatype, p->c);
//	printf("fill = %d elementsize = %ld\n", *((unsigned short *)fptr), elementsize);
	for ( i=0; i<datasize; i++ ) memcpy(newdata+i*elementsize, fptr, elementsize);
	bfree(fptr, elementsize);
	
	for ( n=0; n<pex->n; n++ ) {
		for ( z=0; z<pex->z; z++ ) {
			oldz = z + coords[n].z;
			if ( oldz >= 0 && oldz < p->z ) for ( y=0; y<pex->y; y++ ) {
				oldy = y + coords[n].y;
				if ( oldy >= 0 && oldy < p->y ) for ( x=0; x<pex->x; x++ ) {
					oldx = x + coords[n].x;
					if ( oldx >= 0 && oldx < p->x ) {
						i = ((n*pex->z+z)*pex->y + y)*pex->x + x;
						j = ((img_number*p->z + oldz)*p->y + oldy)*p->x + oldx;
						memcpy(newdata+i*elementsize, data+j*elementsize, elementsize);
					}
				}
			}
		}
		pex->image[n].ox = coords[n].x;
		pex->image[n].oy = coords[n].y;
		pex->image[n].oz = coords[n].z;
	}
	
	img_stats(pex);
		
	return(pex);
}
*/
/************************************************************************
@Function: img_extract_for_show
@Description:
	Converts a slice from an image to a 2D plane for display.
@Algorithm:
	Only the desired slice is extracted from the original image and
	resized based on the given scale argument.
	The dynamic range is rescaled using the image display minimum and 
	maximum:
		new_data = data*255/(max-min)
	with truncation of the data below 0 and above 255.
	Bit data are converted to 0 (black) and 255 (white).
	RGB data are rescaled as for gray scale images.
	Complex data types are converted to polar form and the intensities returned.
@Arguments:
	Bimage* p		image.
	int img_num 	image number to extract.
	int z			slice number to extract.
	float scale 	dimensional scaling.
@Returns:
	Bimage*			extracted slice.
**************************************************************************/
/*Bimage*		img_extract_for_show(Bimage* p, int img_num, int z, float scale)
{
	if ( p->dataflag < 1 ) return(NULL);
	
	int 				threshold = 0;
	if ( fabs(p->smax - p->smin) < 1e-37 ) threshold = 1;
	
	unsigned long		i, j, iy, jy, iz, x, y, xo, yo, c;
	float				iscalex = 1.0/scale;
	float				iscaley = 1.0/scale;
	
    Bimage*     		pshow = init_img_with_parameters(UChar, p->c, 
							(int) (scale*p->x), (int) (scale*p->y), 1, 1);
	pshow->colormodel = p->colormodel;
	if ( pshow->colormodel == Bit ) pshow->colormodel = Gray;
	if ( pshow->colormodel == CMYK ) pshow->colormodel = RGBA;
	
	unsigned char*		data = (unsigned char *) pshow->data;

    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    complex_short*  	csdata = (complex_short *) p->data;
    complex_int*  		cidata = (complex_int *) p->data;
    complex_float*  	cfdata = (complex_float *) p->data;
    
	double*				value = (double *) balloc(p->c*sizeof(double));
    
    float				dscale = 255.0/(p->smax - p->smin);
	
	// Note that the y-axis is flipped
	iz = (img_num*p->z + z)*p->y;
	for ( y=0; y<pshow->y; y++ ) {
		yo = (int) ((pshow->y - y - 1)*iscaley);
		if ( yo >= p->y ) yo = p->y - 1;
		iy = y*pshow->x;
		jy = (iz + yo)*p->x;
		for ( x=0; x<pshow->x; x++ ) {
			xo = (int) (x*iscalex);
			if ( xo >= p->x ) xo = p->x - 1;
			i = (iy + x)*pshow->c;
			j = (jy + xo)*p->c;
			if ( p->colormodel == Bit ) {
				if ( 0x80 & ( udata[j/8] << j%8 ) ) data[i] = 255;
				if ( dscale < 0 ) data[i] = 255 - data[i];
			} else {
				switch ( p->datatype ) {
    				case UChar:
						for ( c=0; c<p->c; c++ )
							value[c] = udata[j+c];
						break;
    				case SChar:
						for ( c=0; c<p->c; c++ )
							value[c] = cdata[j+c];
    	    			break;
    				case Short:
						for ( c=0; c<p->c; c++ )
							value[c] = sdata[j+c];
						break;
    				case UShort:
						for ( c=0; c<p->c; c++ )
							value[c] = usdata[j+c];
						break;
    				case Int:
						for ( c=0; c<p->c; c++ )
							value[c] = idata[j+c];
						break;
    				case Float:
						value[0] = fdata[j];
						break;
    				case ComplexShort:
						value[0] = sqrt(1.0*csdata[j].re*csdata[j].re
								+ csdata[j].im*csdata[j].im);
						break;
    				case ComplexInt:
						value[0] = sqrt(1.0*cidata[j].re*cidata[j].re
								+ cidata[j].im*cidata[j].im);
						break;
    				case ComplexFloat:
						value[0] = sqrt(cfdata[j].re*cfdata[j].re
								+cfdata[j].im*cfdata[j].im);
						break;
    				default:
						value[0] = 0;
						break;
				}
				if ( threshold ) {
					for ( c=0; c<pshow->c; c++ ) {
						if ( value[c] > p->smin ) data[i+c] = 255;
						else data[i+c] = 0;
					}
				} else {
					for ( c=0; c<pshow->c; c++ ) {
						value[c] = floor(dscale*(value[c]-p->smin)+0.5);
						if ( value[c] > 255 ) data[i+c] = 255;
						else if ( value[c] < 0 ) data[i+c] = 0;
						else data[i+c] = (unsigned char) value[c];
					}
				}
				if ( p->colormodel == CMYK ) CMYK_to_RGBA(&data[i]);
			}
		}
    }
	
	bfree(value, p->c*sizeof(double));
	
    return(pshow);
}
*/
/************************************************************************
@Function: img_extract_for_show2
@Description:
	Converts a slice from an image to a 2D plane for display.
@Algorithm:
	Only the desired slice is extracted from the original image and
	resized based on the given scale argument.
	The dynamic range is rescaled using the image display minimum and 
	maximum:
		new_data = data*255/(max-min)
	with truncation of the data below 0 and above 255.
	Bit data are converted to 0 (black) and 255 (white).
	RGB data are rescaled as for gray scale images.
	Complex data types are converted to polar form and the intensities returned.
@Arguments:
	Bimage* p		image.
	int img_num 	image number to extract.
	int z			slice number to extract.
	float scale 	dimensional scaling.
@Returns:
	Bimage*			extracted slice.
**************************************************************************/
/*Bimage*		img_extract_for_show2(Bimage* p, int img_num, int z, float scale)
{
	if ( scale >= 1 ) return(img_extract_for_show(p, img_num,z, scale));
	
	if ( p->dataflag < 1 ) return(NULL);
	
	int 				threshold = 0;
	if ( fabs(p->smax - p->smin) < 1e-37 ) threshold = 1;
	
	unsigned long		i, j, iz, x, y, c;
	long				xo, yo;
	
    Bimage*     		pshow = init_img_with_parameters(UChar, p->c, 
							(int) (scale*p->x), (int) (scale*p->y), 1, 1);
	pshow->colormodel = p->colormodel;
	if ( pshow->colormodel == Bit ) pshow->colormodel = Gray;
	if ( pshow->colormodel == CMYK ) pshow->colormodel = RGBA;
	
	unsigned char*		data = (unsigned char *) pshow->data;

    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    complex_short*  	csdata = (complex_short *) p->data;
    complex_int*  		cidata = (complex_int *) p->data;
    complex_float*  	cfdata = (complex_float *) p->data;
    
	int					showsize = pshow->x*pshow->y*pshow->c;
	int*				num = (int *) balloc(showsize*sizeof(int));
	float*				value = (float *) balloc(showsize*sizeof(float));
    
    float				dscale = 255.0/(p->smax - p->smin);

	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_extract_for_show2: pshow: x=%ld y=%ld\n", pshow->x, pshow->y);
	
	// Note that the y-axis is flipped
	iz = (img_num*p->z + z)*p->y;
	for ( y=0; y<p->y; y++ ) {
		yo = (long) ((p->y - y - 1)*scale - 0.5);
		if ( yo < 0 ) yo = 0;
		if ( yo >= pshow->y ) yo = pshow->y - 1;
		for ( x=0; x<p->x; x++ ) {
			xo = (long) (x*scale + 0.5);
			if ( xo < 0 ) xo = 0;
			if ( xo >= pshow->x ) xo = pshow->x - 1;
			i = (pshow->x*yo + xo)*pshow->c;
			j = ((iz + y)*p->x + x)*p->c;
			for ( c=0; c<p->c; c++ ) num[i+c]++;
			if ( p->colormodel == Bit ) {
				if ( 0x80 & ( udata[j/8] << j%8 ) ) value[i] += 1;
			} else {
				switch ( p->datatype ) {
    				case UChar:
						for ( c=0; c<pshow->c; c++ )
							value[i+c] += udata[j+c];
						break;
    				case SChar:
						for ( c=0; c<pshow->c; c++ )
							value[i+c] += cdata[j+c];
    	    			break;
    				case Short:
						for ( c=0; c<pshow->c; c++ )
							value[i+c] += sdata[j+c];
						break;
    				case UShort:
						for ( c=0; c<pshow->c; c++ )
							value[i+c] += usdata[j+c];
						break;
    				case Int:
						for ( c=0; c<pshow->c; c++ )
							value[i+c] += idata[j+c];
						break;
    				case Float:
						value[i] += fdata[j];
						break;
    				case ComplexShort:
						value[i] += sqrt(1.0*csdata[j].re*csdata[j].re
								+ csdata[j].im*csdata[j].im);
						break;
    				case ComplexInt:
						value[i] += sqrt(1.0*cidata[j].re*cidata[j].re
								+ cidata[j].im*cidata[j].im);
						break;
    				case ComplexFloat:
						value[i] += sqrt(cfdata[j].re*cfdata[j].re
								+cfdata[j].im*cfdata[j].im);
						break;
    				default:
						break;
				}
			}
		}
	}

	for ( y=0; y<pshow->y; y++ ) {
		for ( x=0; x<pshow->x; x++ ) {
			i = (y*pshow->x + x)*pshow->c;
			for ( c=0; c<pshow->c; c++ ) if ( num[i+c] ) value[i+c] /= num[i+c];
			if ( p->colormodel == Bit ) {
				if ( value[i] > 0.5 ) data[i] = 255;
				if ( dscale < 0 ) data[i] = 255 - data[i];
			} else {
				if ( threshold ) {
					for ( c=0; c<pshow->c; c++ ) {
						if ( value[i+c] > p->smin ) data[i+c] = 255;
						else data[i+c] = 0;
					}
				} else {
					for ( c=0; c<pshow->c; c++ ) {
						value[i+c] = floor(dscale*(value[i+c]-p->smin)+0.5);
						if ( value[i+c] > 255 ) data[i+c] = 255;
						else if ( value[i+c] < 0 ) data[i+c] = 0;
						else data[i+c] = (unsigned char) value[i+c];
					}
				}
				if ( p->colormodel == CMYK ) CMYK_to_RGBA(&data[i]);
			}
		}
    }
	
	bfree(num, showsize*sizeof(int));
	bfree(value, showsize*sizeof(float));
    
    return(pshow);
}
*/
/************************************************************************
@Function: img_extract_particles
@Description:
	Extracts particles defined by a series of coordinates.
@Algorithm:
	.
@Arguments:
	Bimage* p			image.
	float radius		radius of box to extract. 
	float bad_radius	radius of bad areas.
	int ncoord			number of coordinates.
	int nbad			number of bad area coordinates.
	float* coords		coordinates (sets of x,y,z,n) including bad areas.
@Returns:
	Bimage*			extracted particle images.
**************************************************************************/
/*Bimage*		img_extract_particles(Bimage* p, float radius, float bad_radius,
				int ncoord, int nbad, float* coords)
{
	if ( p->dataflag < 1 ) return(NULL);
	if ( ncoord < 1 ) return(NULL);
	
	unsigned long		i, j, n, x, y, z;
	long				oldx, oldy, oldz, oldn;
	
	Bimage*				pmask = NULL;
	char*				mask = NULL;
	Vector3				origin;
	if ( nbad && bad_radius > 1 ) {
		pmask = init_img_with_parameters(UChar, p->c, p->x, p->y, p->z, p->n);
		for ( i=ncoord; i<ncoord+nbad; i++ ) {
			origin.x = coords[4*i];
			origin.y = coords[4*i+1];
			origin.z = coords[4*i+2];
			img_delete_within_radius(pmask, origin, bad_radius, 0, 1);
		}
		mask = pmask->data;
	}
	
	VectorInt3          box = {(int) (2.0*radius), (int) (2.0*radius), (int) (2.0*radius)}; // Box size
	if ( box.x > p->x ) box.x = p->x;
	if ( box.y > p->y ) box.y = p->y;
	if ( box.z > p->z ) box.z = p->z;
	VectorInt3          boxrad = {box.x/2, box.y/2, box.z/2};       // Box halfsize
	
    Bimage*     		part = init_img_with_parameters(p->datatype, p->c, box.x, box.y, box.z, ncoord);
	part->ux = p->ux;
	part->uy = p->uy;
	part->uz = p->uz;
	float*				fom = img_init_fom(part);
	
	unsigned long		elementsize = part->c*gettypesize(part->datatype);
	char*               data = (char *) p->data;
	char*               newdata = (char *) part->data;

	for ( n=0; n<part->n; n++ ) {
        oldn = (int) coords[4*n+3];
		part->image[n].ox = boxrad.x;   // Set the origin to the nominal center
		part->image[n].oy = boxrad.y;
		part->image[n].oz = boxrad.z;
		for ( z=0; z<part->z; z++ ) {
			oldz = z - boxrad.z + (int) coords[4*n+2];
			for ( y=0; y<part->y; y++ ) {
				oldy = y - boxrad.y + (int) coords[4*n+1];
				for ( x=0; x<part->x; x++ ) {
					oldx = x - boxrad.x + (int) coords[4*n];
					i = ((n*part->z + z)*part->y + y)*part->x + x;
                    if ( oldx >= 0 && oldx < p->x && oldy >= 0 && oldy < p->y && oldz >= 0 && oldz < p->z) {
                        j = ((oldn*p->z + oldz)*p->y + oldy)*p->x + oldx;
						if ( !mask || (mask && mask[j] < 1) ) {
							memcpy(newdata+i*elementsize, data+j*elementsize, elementsize);
							fom[i] = 1;
						}
					}
				}
			}
		}
	}
	
	kill_img(pmask);
	
	img_correct_background(part);
		
	img_stats(part);
    
	return(part);
}
*/
/************************************************************************
@Function: img_extract_filament
@Description:
	Extracts a filament defined by a series of coordinates.
@Algorithm:
	.
@Arguments:
	Bimage* p		image.
	int img_num 	image number to extract.
	double width	width of filament image to extract. 
	int ncoord		number of coordinates.
	Vector3* coords	coordinates.
@Returns:
	Bimage*			extracted slice.
**************************************************************************/
/*Bimage*		img_extract_filament(Bimage* p, int img_num, double width,
				int ncoord, Vector3* coords)
{
	if ( p->dataflag < 1 ) return(NULL);
	if ( img_num < 0 ) return(NULL);
	if ( ncoord < 1 ) return(NULL);
	
	unsigned long   i, j, n, x, y, z, xo, yo, zo;
	int 			nspline, nx, ny, nz, xx, yy, zz, yw, zw;
	float			xf[2], yf[2], zf[2];
	float			halfwidth = width/2;
	Vector3 		vref = {1,0,0}, v, vf = {0,0,0};
	Matrix3			mat;
	
	Vector3*		spline = vector3_catmull_rom_spline(ncoord, coords, &nspline);
	
	VectorInt3		size = {nspline, (int)width, (int)width};
	if ( p->z < 2 ) size.z = 1;

	nx = ny = nz = 2;		// Number of voxels to interpolate over
	if ( nx > p->x ) nx = 1;
	if ( ny > p->y ) ny = 1;
	if ( nz > p->z ) nz = 1;
	
    Bimage*     		pfil = init_img_with_parameters(Float, 1, size.x, size.y, size.z, 1);
	float*				fildata = (float *) pfil->data;

    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
	
	for ( n=0; n<pfil->n; n++ ) {
		for ( z=0; z<pfil->z; z++ ) {
			vf.z = z - halfwidth;
			if ( pfil->z < 2 ) vf.z = 0;
			for ( y=0; y<pfil->y; y++ ) {
				vf.y = y - halfwidth;
				for ( x=0; x<pfil->x; x++ ) {
					if ( x==0 )
						v = vector3_subtract(spline[x+1], spline[x]);
					else
						v = vector3_subtract(spline[x], spline[x-1]);
					mat = matrix3_from_two_vectors3(vref, v);
					v = vector3_matrix3_multiply(mat, vf);
					v = vector3_add(spline[x], v);
					if ( v.z >= 0 && v.z < p->z ) {
						if ( v.y >= 0 && v.y < p->y ) {
							if ( v.x >= 0 && v.x < p->x ) {
								i = ((n*pfil->z + z)*pfil->y + y)*pfil->x + x;
								xo = (int) v.x;
								yo = (int) v.y;
								zo = (int) v.z;
								zf[1] = v.z - zo;
								zf[0] = 1.0 - zf[1];
								yf[1] = v.y - yo;
								yf[0] = 1.0 - yf[1];
								xf[1] = v.x - xo;
								xf[0] = 1.0 - xf[1];
								for ( zz=0; zz<nz; zz++ ) {
									zw = (img_num*p->z + zo + zz)*p->y;
									for ( yy=0; yy<ny; yy++ ) {
										yw = (zw + yo + yy)*p->x;
										for ( xx=0; xx<nx; xx++ ) {
											j = yw + xo + xx;
											switch ( p->datatype ) {
    											case UChar:
													fildata[i] += udata[j]*xf[xx]*yf[yy]*zf[zz];
													break;
    											case SChar:
													fildata[i] += cdata[j]*xf[xx]*yf[yy]*zf[zz];
    	    										break;
    											case Short:
													fildata[i] += sdata[j]*xf[xx]*yf[yy]*zf[zz];
													break;
    											case UShort:
													fildata[i] += usdata[j]*xf[xx]*yf[yy]*zf[zz];
													break;
    											case Int:
													fildata[i] += idata[j]*xf[xx]*yf[yy]*zf[zz];
													break;
    											case Float:
													fildata[i] += fdata[j]*xf[xx]*yf[yy]*zf[zz];
													break;
    											default:
													break;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	bfree(spline, nspline*sizeof(Vector3));
	
	return(pfil);
}
*/
/************************************************************************
@Function: img_extract_panels
@Description:
	Extracts panels defined by a series of coordinates.
@Algorithm:
	The image is assumed to be the result of segmentation with panels,
	each defined by a single value, typically an integer.
	The new image is a mask (0,1) of type single-channel unsigned char 
	and the same size as the input image.
@Arguments:
	Bimage* p		image.
	int ncoord		number of coordinates.
	float* coords	coordinates (sets of x,y,z,n).
@Returns:
	Bimage*			mask image with all extracted panels.
**************************************************************************/
/*Bimage*		img_extract_panels(Bimage* p, int ncoord, float* coords)
{
	if ( p->dataflag < 1 ) return(NULL);
	if ( ncoord < 1 ) return(NULL);
	
	unsigned long       i, j, k, x, y, z, n;
	double				diff, tol = 1e-30;
	
    Bimage*     		pan = init_img_with_parameters(UChar, 1, p->x, p->y, p->z, p->n);
	for ( i=0; i<p->n; i++ ) pan->image[i] = p->image[i];
	
	unsigned long		datasize = (unsigned long) p->c*p->x*p->y*p->z*p->n;
	
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
	
	char*               mask = (char *) pan->data;
	
	float*				v = (float *) balloc(ncoord*sizeof(float));

	for ( j=k=0; k<ncoord; j+=4, k++ ) {
		x = (unsigned long) (coords[j] + 0.5);
		y = (unsigned long) (coords[j+1] + 0.5);
		z = (unsigned long) (coords[j+2] + 0.5);
		n = (unsigned long) (coords[j+3] + 0.5);
		i = ((n*p->z + z)*p->y + y)*p->x + x;
		switch ( p->datatype ) {
			case UChar:
				if ( p->colormodel == Bit ) v[k] = 0x80 & ( udata[i/8] << i%8 );
				else v[k] = udata[i];
				break;
			case SChar: v[k] = cdata[i]; break;
			case UShort: v[k] = usdata[i]; break;
			case Short: v[k] = sdata[i]; break;
			case Int: v[k] = idata[i]; break;
			case Float: v[k] = fdata[i]; break;
			default: break;
		}
    }
    
	for ( i=0; i<datasize; i++ ) {
		for ( k=0; k<ncoord; k++ ) {
			switch ( p->datatype ) {
				case UChar:
					if ( p->colormodel == Bit ) {
						diff = v[k] - (0x80 & ( udata[i/8] << i%8 ));
						if ( fabs(diff) < tol ) mask[i] = 1;
					} else if ( fabs(v[k] - udata[i]) < tol ) {
						mask[i] = 1;
					}
					break;
				case SChar: if ( fabs(v[k] - cdata[i]) < tol ) mask[i] = 1; break;
				case UShort: if ( fabs(v[k] - usdata[i]) < tol ) mask[i] = 1; break;
				case Short: if ( fabs(v[k] - sdata[i]) < tol ) mask[i] = 1; break;
				case Int: if ( fabs(v[k] - idata[i]) < tol ) mask[i] = 1; break;
				case Float: if ( fabs(v[k] - fdata[i]) < tol ) mask[i] = 1; break;
				default: break;
			}
		}
	}
	
	bfree(v, ncoord*sizeof(float));
	
	img_stats(pan);
    
	return(pan);
}
*/

