/*
	buffer.h
	Header file for reading files in buffer mode
	Author: Giovanni Cardone
	Created: 20050621 	Modified:
*/

#ifndef _BUFFER_H__
#define _BUFFER_H__

#include "bsoft.h"

/************************************************************************
@Object: struct Buff
@Description:
	Buffer read mode structure.
@Features:
	Buffer mode structure.
*************************************************************************/
struct Buff {
	long long	first;		// Index of first element in buffer
	long long	last;		// Index of last element in buffer
							// if equal to zero, the data buffer has not been created yet
	long long	size;		// Size of buffer (given in bytes)
} ;


// Function prototypes
Buff* 	buff_init(long long size);
int 		buff_kill(Buff* bf);
char*	buff_img_read_data(FILE* fimg, Bimage* p, int select_img, int swap, int vax, int pad, Buff* bf);

#endif  // #ifndef _BUFFER_H__

