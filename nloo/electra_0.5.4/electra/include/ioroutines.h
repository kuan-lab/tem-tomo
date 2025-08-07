/*
	ioroutines.h
	Header file to read/write images
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004
*/

#ifndef _IOROUTINES_H__
#define _IOROUTINES_H__

#include "bsoft.h"
#include "proj.h"
#include "buffer.h"


// Constants
const int   IO_HEADER = -3;
const int   IO_ZEROS = -2;
const int   IO_ALL = -1;

	
// Function prototypes
Bimage* io_read_image(char* filename, int readdataflag, int select_img);
Bimage* io_read_image_buff(Bimage* p, char* filename, int select_img, Buff* bf);
int	 	io_readMRC_header(FILE *fimg, Bimage* p, int select_img, int *swap, int *vax);
int		io_readMRC(Bimage* p, int select_img);
int 		io_readMRC_buff(Bimage* p, Buff* bf);
int 		io_pwrite(char* filename, Proj* p, int writeflag, int nzp);
int 		io_update_stats(char* filename);
int 		io_update_projfile_stats(char* filename);
long long io_get_index(int x, int y, int z, Bimage* p, Buff* bf);
char* 	io_read_data(FILE* fimg, Bimage* p, int select_img, int swap, int vax, int pad);
int		io_write_image(char* filename, Bimage* p);

#endif  // #ifndef _IOROUTINES_H__
