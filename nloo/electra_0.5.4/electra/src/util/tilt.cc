/* 
	tilt.cc
	Functions to handle tilt angles 
	Author: Giovanni Cardone
	Created: 2004	Modified: 2004
*/ 

#include "tilt.h"


/************************************************************************
@Function: tlt_load_tilt_angles
@Description:
	Reads from file the value of title angles.
@Algorithm:
	A single value is read from each line in the file.
@Arguments:
	char* filename		file name.
	int* ntlt_angs		number of angles read.
@Returns:
	float*				series of tilt angles.
*************************************************************************/
float*			tlt_load_tilt_angles(char* filename, int *ntlt_angs)
{
	int			n = 0;
	FILE*		fd;
	float*		tlt_ang;
	float*		tmparray = (float *) balloc(MAXNANGLES*sizeof(float));

	int i;
	
	if ( verbose & VERB_LABEL )
		printf("Reading tilt angles from file %s\n", filename);

	if ( ( fd = fopen(filename, "r+") ) == NULL ) {
		error_show(filename, __FILE__, __LINE__);
		return(NULL);
	}
	
	while(!feof(fd) && n<=MAXNANGLES) {       /* loop through and store the numbers into the array */
		if(fscanf(fd, "%f\n", &tmparray[n]) < 1) {
			fprintf(stderr, "Error: verify %dth value in file %s!\n", n+1,filename);
			exit(-1);
		}
		n++;
    }

	if ( !feof(fd) )
		printf("WARNING: only the first %d values have been read from file %s!\n\n", MAXNANGLES,filename);

    fclose(fd);
	
	tlt_ang = (float *) balloc(n*(sizeof(float)));

	for(i=0 ; i<n ; i++)
		tlt_ang[i] = tmparray[i]*M_PI/180.0;
	
	bfree(tmparray,MAXNANGLES*sizeof(float));

	if ( verbose & VERB_FULL ) {
		printf("Number of angles read: %d\n\n", n);
		printf("The angles are:\n");

		for(i=0 ; i<n ; i++)
			printf("%f\n", tlt_ang[i]*180.0/M_PI);
		printf("\n");
	}

	*ntlt_angs = n;

	return(tlt_ang);
}

/************************************************************************
@Function: tlt_create_series_from_range
@Description:
	Create tilt series from given range and step.
@Algorithm:
@Arguments:
	float amin			minimum angle.
	float amax			maximum angle.
	float astep			incremental step.
	int* ntlt_angs		number of angles created.
@Returns:
	float*				series of tilt angles.
*************************************************************************/
float*			tlt_create_series_from_range(float amin, float amax, float astep, int* ntlt_angs)
{
	int i;
	float fsign;
	int			n = (int) ((amax-amin)/astep) + 1;
	float*		tlt_ang = (float *) balloc(n*sizeof(float));
	
	for ( i=0 ; i<n ; i++)
		tlt_ang[i] = (amin + astep*i);
	
	for ( i=0 ; i<n ; i++) {
//		printf("%d %f\n",i,tlt_ang[i]);
		fsign = tlt_ang[i]>-0.5?1.:-1.;
		tlt_ang[i] = 0.1 * (int) ( tlt_ang[i] * 10.0 + fsign*0.5);
//		printf("%d %f\n",i,tlt_ang[i]);
	}
	
	for ( i=0 ; i<n ; i++)
		tlt_ang[i] *= M_PI/180.0;
		
	if ( verbose & VERB_PROCESS ) {
		printf("Number of angles: %d\n\n", n);
		printf("The angles are:\n");

		for ( i=0 ; i<n ; i++)
			printf("%f\n", tlt_ang[i]*180.0/M_PI);
		printf("\n");
	}

	*ntlt_angs = n;

	return(tlt_ang);
}

/************************************************************************
@Function: tlt_save_angles_to_file
@Description:
	Save to file the values of title angles.
@Algorithm:
	The extension .tlt is appended to the input filename.
	Each angle is written to a separate line.
@Arguments:
	float* tlt_angs		tilt series.
	int ntlt_angs		number of angles to write.
	char* filename		file name.
@Returns:
	int					<0 if not successfull
*************************************************************************/
int			tlt_save_angles_to_file(float* tlt_angs, int ntlt_angs, char* filename)
{
	FILE*		fd;
	char *		tilt_name = NULL;
	int i;
	float fsign;
	
	tilt_name = copystring(filename);
	
	// change extension of file into tlt, if needed
	if( strcmp(extension(tilt_name),"tlt") )
		strcpy(tilt_name, filename_change_type(tilt_name,"tlt")); 

	if ( verbose & VERB_LABEL )
		printf("Writing tilt angles to file %s\n", tilt_name);

	if ( ( fd = fopen(tilt_name, "w") ) == NULL ) {
		error_show(tilt_name, __FILE__, __LINE__);
		return(-1);
	}
	
	for(i=0 ; i<ntlt_angs ; i++) {
		fsign = tlt_angs[i]*180.0/M_PI>-0.5?1.:-1.;
		fprintf(fd,"%5.1f\n", 0.1 * (int) ( tlt_angs[i] *180.0/M_PI * 10.0 + fsign*0.5 ));
	}
	
    fclose(fd);

	bfree_string(tilt_name);

	return(0);
}

/************************************************************************
@Function: tlt_views_from_angles
@Description:
	Initializes a set of views tilted around the y axis in a reference
	system determined through a given view with respect to an ideal
	reference system.
@Algorithm:
	Calculates the quaternion representing a rotation from the reference
	to the given view.
	Then, for each tilt angle its quaternion is evaluated from the
	corresponding view with respect to the local reference system.
	The requested view is finally obtained from the quaternion 
	resulting from concatenating the two quaternions/rotations.
@Arguments:
	int nviews			number of views returned.
	View tilt0_view		Reference view for tilt angles.
	float* angs			series of tilt angles (radians).
@Returns:
	View* 				a set of 4-value views.
**************************************************************************/
View* 		tlt_views_from_angles(int nviews, View tlt0_view, float* tlt_angs)
{
	int				i;
	Quaternion		q, qti;
	View			vti;
	View*			view = (View *) balloc(nviews*sizeof(View));

	if ( verbose & VERB_LABEL )
		printf("Calculating %d views for tilt angles given with respect to view (%f,%f,%f,%f)\n", 
			nviews, tlt0_view.x, tlt0_view.y, tlt0_view.z, tlt0_view.a*180.0/M_PI);

	Quaternion qt0 = quaternion_from_view(tlt0_view);
	
	for ( i=0; i<nviews; i++ ) {
//		vti.x = sin(tlt_angs[i])*cos(xyrot);
//		vti.y = sin(tlt_angs[i])*sin(xyrot);
		vti.x = sin(tlt_angs[i]);
		vti.y = 0.;
		vti.z = cos(tlt_angs[i]);
		vti.a = 0.0;
		vti.next = NULL;

		qti = quaternion_from_view(vti);
//		q = quaternion_multiply(qti,qt0);
		q = quaternion_multiply(qt0,qti);
		
		view[i] = view_from_quaternion(q);
	}

	if ( verbose & VERB_FULL ) {
		printf("The tilt angles\t/\t views are:\n");

		for(i=0 ; i<nviews ; i++)
			printf("%f\t/\t(%f,%f,%f,%f)\n",
				tlt_angs[i]*180.0/M_PI, view[i].x, view[i].y, view[i].z, view[i].a*180.0/M_PI);
		printf("\n");
	}

	return(view);
}
