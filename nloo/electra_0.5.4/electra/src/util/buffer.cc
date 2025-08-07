/*
	buffer.cc
	Functions for handling buffered images
	Author: Giovanni Cardone
	Created: 20050621 	Modified: 
*/

#include "buffer.h"

/************************************************************************
@Function: buffer_init
@Description:
	General buffer structure initialization.
@Algorithm:
	This function allocates memory and sets up a number of defaults.
@Arguments:
	long long	size		size of the buffer (bytes)
@Returns:
	Buff*				new buffer structure, NULL if initialization failed.
**************************************************************************/
Buff* 	buff_init(long long size)
{
	// Allocate memory for the image parameter structure
	Buff* 	bf = (Buff *) balloc(sizeof(Buff));
	if ( !bf ) return(bf);

	// Set parameter defaults
	bf->first = 0;
	bf->last = 0;
	bf->size = size;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG buff_init: Buffer setup done\n");
	
	return(bf);
}

/************************************************************************
@Function: buff_kill
@Description:
	General buffer structure destruction.
@Algorithm:
	This function deallocates all memory associated to the structure
@Arguments:
	Buff*				buffer structure
@Returns:
	int					error code (<0 means failure)
**************************************************************************/
int 	buff_kill(Buff* bf)
{

	bfree(bf, sizeof(Buff));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG buff_kill memory = %ld\n", memory);
	
	return(0);
}

/************************************************************************
@Function: buff_img_read_data
@Description:
	Read image data in buffered mode.
@Algorithm:
	If the buffer is not initialized, the data array of the image is
	allocated. The index and size of the buffer determine the portion of the
	whole file to load.
	Selection of a single image from a file is not yet supported.
	The data is read in the largest blocks possible for efficiency.
	Any interspersed padding and page sizes not matching the data size
	contribute to inefficiency in reading.
	Swapping:	swap = 1:	swap bytes the size of the data type
				swap > 1:	swap these number of bytes regardless of the data type
@Arguments:
	FILE* fimg		file descriptor: file opened in calling function.
	Bimage* p		image structure: header parameters filled by calling function.
	int select_img	image selection: if -1, all images, if >= 0, one image.
	int swap			flag activates byte swapping.
	int vax 			indicate vax style floating point - activates conversion.
	int pad 			any interspersed control or separation bytes.
	Buff* bf			buffer structure
@Returns:
	char*				data in one big contiguous block, NULL if reading failed.
**************************************************************************/
char*	buff_img_read_data(FILE* fimg, Bimage* p, int select_img, int swap, int vax, int pad, Buff* bf)
{
	if ( !p ) return(NULL);
	
	if ( p->dataflag < 1 ) return(NULL);

	if (select_img>=0) {
		fprintf(stderr,"buff_img_read_data: selection of images not supported!\n");
		exit(-1);
	}	
	if ( p->px < 1 ) p->px = p->x;
	if ( p->py < 1 ) p->py = p->y;
	if ( p->pz < 1 ) p->pz = p->z;
	
	// If only half of a transform is stored, it need to be handled
	long long		xstore = p->x;
	long long		xpage = p->px;
	if ( p->transform == Hermitian || p->transform == CentHerm ) {
		xstore = p->x/2 + 1;
		if ( p->px > xstore ) xpage = xstore;
	}
	
	unsigned long		i;
	unsigned long		datatypesize = gettypesize(p->datatype);
	unsigned long		valuesize = datatypesize*p->c;
	long long			datasize = (long long) p->n*xstore*p->y*p->z*valuesize;
	long long			pagesize = (long long) xpage*p->py*p->pz*valuesize;

	long long		bf_size = datasize;
	if (bf && bf->size < datasize) bf_size = bf->size;
	else if (bf && datasize < bf->size ) bf->size = datasize;

	long long		readsize = bf_size;
	if (bf && bf->first*valuesize+readsize > datasize)
		readsize = datasize - bf->first*valuesize;

	char* data = NULL;
	if (bf && bf->last!=0) data = p->data;
	else data = (char *) balloc(bf_size*sizeof(char));

	
	long long bf_offset = 0;
	if (bf) {
		bf->last = bf->first - 1 + readsize/valuesize;
		bf_offset = bf->first*valuesize;
	}
		
	char*				page = NULL;

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG buff_img_read_data: Pointer = %p  Data type size = %ld\n", p, datatypesize);
		printf("DEBUG buff_img_read_data: Data size: %ld %ld %ld %ld %ld (%lld)\n", p->x, p->y, p->z, p->c, p->n, datasize);
		printf("DEBUG buff_img_read_data: Page size: %ld %ld %ld (%lld)\n", p->px, p->py, p->pz, pagesize);
		if (bf) printf("DEBUG buff_img_read_data: Buffer size: %lld , data read: %lld (from %lld to %lld)\n", bf_size, readsize, bf->first*valuesize, (bf->last+1)*valuesize-1);
		if (bf) printf("DEBUG buff_img_read_data: Buffer size: %lld , data read: %lld (from %lld to %lld)\n", bf_size, readsize, bf->first, bf->last);
		printf("DEBUG buff_img_read_data: Swap = %d  Vax = %d  Pad = %d  Offset = %ld\n", swap, vax, pad, p->offset);
	}
	
	if ( !pad && p->x == p->px && p->y == p->py && p->z == p->pz ) { // 3D block
		if ( verbose & VERB_DEBUG )
			printf("DEBUG buff_img_read_data: Reading 3D blocks\n");
		fseeko( fimg, (off_t) p->offset + bf_offset, SEEK_SET );
		if ( pad ) page = (char *) balloc(pad*sizeof(char));
		fread( data, readsize, 1, fimg );
		if ( page ) bfree(page, pad*sizeof(char));
	} else { printf("ERROR buff_img_read_data: can not read data.\n"); exit(-1);}

	if ( verbose & VERB_DEBUG )
		printf("DEBUG buff_img_read_data: Converting image data\n");

	// Convert vax format and swap bytes if required
    if ( vax && swap < 2 && ( p->datatype == Float ) )
		for ( i=0; i<readsize; i+=4 )
    	    vax2ieee(data+i, 1-swap);
    else if ( swap == 1 )
		for ( i=0; i<readsize; i+=datatypesize ) 
			swapbytes(data+i, (unsigned int) datatypesize);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG buff_img_read_data: Finished reading and converting data\n");

	return(data);
}
