/* 
	rwmatrix.h 
	Header file for reading (and writing) matrices.
	Author: Bernard Heymann 
	Created: 20010723  	    Modified: 20030208
*/ 
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "string_util.h"

#ifndef _Bmatrix_
/************************************************************************
@Object: struct Bmatrix
@Description:
	General square matrix structure.
@Features:
	An nxn matrix with n column(row) labels.
	Intended to be used for pairwise comparisons.
*************************************************************************/
struct Bmatrix {
	int size;					// Matrix size = number of columns = number of rows
	char** label;				// Column labels
	float* data; 				// Contents
};
#define _Bmatrix_
#endif

// Function prototypes
Bmatrix*	read_matrix(char* filename);
int 		write_matrix(char* filename, Bmatrix* matrix);
int 		kill_matrix(Bmatrix* matrix);
