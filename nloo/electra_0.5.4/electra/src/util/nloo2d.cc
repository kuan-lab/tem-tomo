/* 
	nloo2d.cc
	Functions to handle nloo2d curves
	Author: Giovanni Cardone
	Created: 20050317	Modified:
*/ 

#include "nloo2d.h"

/************************************************************************
@Function: n2d_init
@Description:
	initialize a nloo2d structure.
@Algorithm:
@Arguments:
	int    nc	number of curves
	int    ns	number of samples per curve
@Returns:
	N2d_curve*  nloo2d initialized
**************************************************************************/
N2d_curve* n2d_init( int nc, int ns)
{
	
	N2d_curve*  n2d = (N2d_curve *) balloc(sizeof(N2d_curve));

	n2d->nc = nc;
	n2d->ns = ns;
	n2d->tlt_angs = (float *) balloc(nc*sizeof(float));
	n2d->axis = (float *) balloc(ns*sizeof(float));
	n2d->res = (float *) balloc(ns*sizeof(float));
	n2d->data = (float *) balloc(nc*ns*sizeof(float));
		
	return(n2d);
}

/************************************************************************
@Function: n2d_kill
@Description:
	General nloo2d structure destruction.
@Algorithm:
	This function deallocates all memory associated to the structure
@Arguments:
	N2d_curve*  nloo2d structure
@Returns:
	int		error code (<0 means failure)
**************************************************************************/
int n2d_kill( N2d_curve* n2d)
{
	
	if ( n2d == NULL ) return(0);
	
	bfree(n2d->tlt_angs,n2d->nc*sizeof(float)); 
	bfree(n2d->axis,n2d->ns*sizeof(float)); 
	bfree(n2d->res,n2d->ns*sizeof(float)); 
	bfree(n2d->data,n2d->ns*n2d->nc*sizeof(float));

	bfree(n2d, sizeof(N2d_curve));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG n2d_kill memory = %ld\n", memory);

	return(0);
}

/************************************************************************
@Function: n2d_read_file
@Description:
	Reads from file the nloo2d curves.
@Algorithm:
	Check size and number of curves stored into the file, then init the
	related structure and load the curves into it.
@Arguments:
	char* filename		file name.
@Returns:
	N2d_curve*			structure containing the read curves.
*************************************************************************/
N2d_curve*	n2d_read_file(char* filename)
{
	const 	int NCHARS = 250;
	char  	fbuf[NCHARS];
//	char		auxstr[80];
//	char		auxch;
	FILE*	fd;

	int ns=0, nc=0, tmpns = 0;
	int i, auxi;
	float auxf, auxa;
	
	if ( verbose & VERB_PROCESS )
		printf("Reading nloo2d curves from file %s\n", filename);

	if ( ( fd = fopen(filename, "r") ) == NULL ) {
		error_show(filename, __FILE__, __LINE__);
		return(NULL);
	}

	// determine size and number of curves
	
	while(!feof(fd)) {

		fgets(fbuf, NCHARS, fd);
		if ( fbuf[0] !='P') tmpns++;
		else {
			nc++;
			if (nc>2 && tmpns!=ns) {
				fprintf(stderr,"Error(n2d_read_file): size of curves not consistent! (%d vs %d)\n",tmpns,ns);
				return(NULL);
			} else if (nc==2) ns = tmpns;
			fgets(fbuf, NCHARS, fd);  // skip comment line
			tmpns = 0;				// reset counter
		}			
	}
	
	N2d_curve* n2dc = n2d_init(nc,ns);

	rewind(fd);

	int is=0 ,ic=0;

	// read values
	while(!feof(fd)) {

		fgets(fbuf, NCHARS, fd);
		if ( fbuf[0] !='P') {
			i = (ic-1)*ns+is;
			if (ic==1) 
				sscanf(fbuf,"%d\t%f\t%f\t%f\n",&auxi,&n2dc->axis[i],&n2dc->res[i],&n2dc->data[i]);
			else
				sscanf(fbuf,"%d\t%f\t%f\t%f\n",&auxi,&auxa,&auxf,&n2dc->data[i]);
			is++;
		}
		else {
			// read tilt angle
			sscanf(fbuf,"Projection # %d - Tilt angle: %f\n",&auxi,&n2dc->tlt_angs[ic]);
			fgets(fbuf, NCHARS, fd);  // skip comment line
			ic++;
			is = 0;
		}			
	}
	
    fclose(fd);

	return(n2dc);
}
