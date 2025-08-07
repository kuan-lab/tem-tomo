/*
	string_tools.cc
	Routines for manipulating strings
	Author: Giovanni Cardone
	Created: 20050125 	Modified:
*/

#include "string_tools.h"

/************************************************************************
@Function: str_get_list
@Description:
	Get a list of integers from a string.
@Algorithm:
	Single values as well as range values are read by scanning through the
	string. Range values are expanded into all the values in between.
@Arguments:
	int*		il		integer list array (output)
	int		maxil	size of list integer array
	char*	sbuf		input string buffer
@Returns:
	int 	 			number of values read (<0 if error)
**************************************************************************/
int str_get_list(int* il, int maxil, char* sbuf) {

	if ( sbuf==NULL) return(-1);

	int nv = 0;
	unsigned int sbufln= (unsigned int) strlen(sbuf);

	unsigned int i = 0;
	int tmin, tmax, is, iss;
	char si[3];

	while ( i< sbufln && !isdigit(sbuf[i]) ) i++;
	while ( i < sbufln ) {
		tmin = 0;
		tmax = 0;
		// read first number
		si[0]=' ';si[1]=' ';si[2]=' ';
		is = 0;
		while ( isdigit(sbuf[i]) ) si[is++] = sbuf[i++];
		// tmin = atoi(si);
		tmin = (int) strtol(si, (char **)NULL, 10);
		//scan and check if a range is defined
		while ( i< sbufln && isspace(sbuf[i]) ) i++;	
		if ( i<sbufln && sbuf[i] == '-' ) {
			si[0]=' ';si[1]=' ';si[2]=' ';
			iss = 0;
			while ( i< sbufln && !isdigit(sbuf[i]) ) i++;
			while ( isdigit(sbuf[i]) ) si[iss++] = sbuf[i++];
			tmax = (int) strtol(si, (char **)NULL, 10);
		}
		if ( tmax == 0) tmax = tmin;
		if(nv+tmax-tmin+1 >= maxil) {
			fprintf(stderr,"ERROR(str_get_list): number of input items exceeds size of output array!\n");
			return(-1);
		}
		for ( int it = tmin; it <= tmax; it++)
			il[nv++] = it;
		// skip to next value
		while ( i< sbufln && !isdigit(sbuf[i]) ) i++;
	}

	return(nv);
}

/************************************************************************
@Function: str_get_minmax
@Description:
	Get the minimum and maximum value in a list of integers from a string.
@Algorithm:
	It also returns the number of values read, which can be used for check, if only
	a continuous range is required.
@Arguments:
	int*		vmin		minimum value (output)
	int*		vmax  	maximum value (output)
	char*	sbuf		input string buffer
@Returns:
	int 	 			number of values read (<0 if error)
**************************************************************************/
int str_get_minmax(int* vmin, int* vmax, char* sbuf)
{

	const int IMAXBUFFER = 500;

	int ilist[IMAXBUFFER];
	
	int	in = str_get_list(ilist, IMAXBUFFER, sbuf);

	int imin = 0;
	int imax = 0;
	
	if ( in != 0 ) imin = imax = ilist[0];

	for (int i = 1; i < in; i++) {
		if (imin > ilist[i]) imin = ilist[i];
		if (imax < ilist[i]) imax = ilist[i];
	}

	*vmin = imin;
	*vmax = imax;
	
	return(in);
}
