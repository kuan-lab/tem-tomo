/*
	md.h
	Header for molecular dynamics
	Author: Bernard Heymann
	Created: 20010828  	    Modified: 20050202
*/

#include "rwmd.h"
#include "rwmolecule.h"
#include "mol_util.h"
#include "utilities.h"

// Function prototypes
float		md_leapfrog(Bmolgroup* molgroup, Bangle* anglelist,
					Bmd* md, int max_iter, float velocitylimit);
double		md_bond_forces(Bmolgroup* molgroup, float Kbond, int wrap);
double		md_angular_forces(Bmolgroup* molgroup, Bangle* anglelist, float Kangle, int wrap);
double		md_nonbonded_forces(Bmolgroup* molgroup, Bmd* md);
int			atom_nonbonded_forces(Batom* atom1, Batom* atom2, Bmd* md, Vector3 box);



