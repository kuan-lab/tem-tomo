/* 
	imd_clip.cc
	Functions to launch clip program by imod 
	Author: Giovanni Cardone
	Created: 20050201	Modified:
*/ 

#include "imod.h"

/************************************************************************
@Function: imd_clip_flip
@Description:
	Flip a volume.
@Algorithm:
@Arguments:
	char*	fmode		flip mode (xy, yz, etc.)
	char* 	in_file		input file
	char*	out_file		output file
@Returns:
	int				<0 if execution failed
*************************************************************************/
int imd_clip_flip(char* fmode, char* in_file, char* out_file) {

	char* cargs = (char *) balloc((13+strlen(in_file)+strlen(out_file))*sizeof(char));
	sprintf(cargs, "clip flip%s %s %s", fmode, in_file, out_file);
	if ( system(cargs) <0 ) return(-1);

//	char* cargs = (char *) balloc((10+strlen(in_file)+strlen(out_file))*sizeof(char));
//	sprintf(cargs, "flip%s %s %s", fmode, in_file, out_file);
//	if ( ext_execute("clip",cargs) <0 ) return(-1);

	bfree_string(cargs);
	
	return (0);
}
