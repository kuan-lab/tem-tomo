/*
	rwsymop.h
	Header file for reading and writing symmetry operators
	Author: Bernard Heymann
	Created: 19990509 	Modified: 20030316
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "rwimg.h"
#include "matrix.h"

// Function prototypes in symmetry.c and rwsymop.c
float* 		read_symat(char* filename, int spacegroup, int* nsym);
char* 		read_symop(char* filename, int spacegroup, int* nsym);
int 		write_symat(char* filename, int spacegroup);
int 		write_pointgroup(char* filename, char* symmetry_string, View ref_view);

