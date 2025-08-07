/*
	rwatomprop.h:
	Header file for reading atom property files
	Author: Bernard Heymann
	Created: 19980822 	Modified: 20030705
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>

#ifndef _Batomtype_
/* Atom type structure */
struct Batom_type {
	Batom_type*	next;		// Link to next atom type if not NULL
	int 		z;			// Atomic number
	char		el[4];		// Element
	char		name[8];	// Atom name
	float		mass;		// Atomic mass
	float		oxid;		// Oxidation state
	float		bond;		// Bond length
	float		vdw;		// VdW radius
	char		color[4];	// RGBA colour
	float		sfa[5];		// Cromer-Mann coefficients
	float		sfb[5];
	float		sfc;
} ;

#define _Batomtype_
#endif

/* Function prototypes */
Batom_type*		get_atom_properties(char *filename);
int 			write_atom_properties(char* filename, Batom_type* at);


