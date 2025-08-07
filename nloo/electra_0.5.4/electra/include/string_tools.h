/*
	string_tools.h
	Header file for string manipulation routines
	Author: Giovanni Cardone
	Created: 20050125	Modified:
*/

#ifndef _STRING_TOOLS_H__
#define _STRING_TOOLS_H__

#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "bsoft.h"

/*----------Constants-----------------------*/

/*----------Function Prototypes-------------*/
int str_get_list(int* il, int maxil, char* sbuf);
int str_get_minmax(int* vmin, int* vmax, char* sbuf);

#endif  // #ifndef _STRING_TOOLS_H__
