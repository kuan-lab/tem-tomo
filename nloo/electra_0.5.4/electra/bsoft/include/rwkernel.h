/*
	rwkernel.h
	Header file for functions to read and write kernel files
	Author: Bernard Heymann
	Created: 20031102 	Modified: 20040513
*/

#include "matrix.h"

// Function prototypes
float* 		read_kernel_file(char* filename, VectorInt3* kernelsize);
int 		write_kernel_file(char* filename, VectorInt3* kernelsize, float* kernel);
