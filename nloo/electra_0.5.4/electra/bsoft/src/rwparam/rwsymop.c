/*
	rwsymop.c
	Library routines to read and write symmetry operators
	Author: Bernard Heymann
	Created: 19991225 	Modified: 20030609
*/

#include "rwsymop.h"
//#include "rwstar.h"
#include "symmetry.h"
//#include "star_tags.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

// Internal function prototypes
//char* 		read_symop_star(char* filename, int spacegroup, int* nsym);
char* 		read_symop_lib(char* filename, int spacegroup, int* nsym);
//int 		write_symop_star(char* filename, int spacegroup, int nsym, char* symop, int line_len);
//int 		write_star_pointgroup(char* filename, char* symmetry_string, View ref_view);

/************************************************************************
@Function: read_symat
@Description:
	Reading crystallographic symmetry operators.
@Algorithm:
	.
@Arguments:
	char* filename		file name.
	int spacegroup		crystal space group number.
	int* nsym			number of symmetry operators.
@Returns:
	float* 				set of 12-value symmetry matrices.
**************************************************************************/
/*float* 		read_symat(char* filename, int spacegroup, int* nsym)
{
	char* 		symop = read_symop(filename, spacegroup, nsym);
	
	// Set up the symmetry operator matrix
	float*		mat = sym_matrices_from_text_list(*nsym, symop, 80);
	
	bfree(symop, (*nsym)*80*sizeof(char));
	
	return(mat);
}
*/
char* 		read_symop(char* filename, int spacegroup, int* nsym)
{
	// No space group has more than 192 operators
	if ( spacegroup < 1 || spacegroup > 100000 ) return(0);
	
	int			default_filename_used = 0;
	if ( !filename ) {
		filename = parameter_file_path("symop.star");
		default_filename_used = 1;
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG read_symop: Reading symmetry operator file %s\n", filename);

	char*		ext = extension(filename);
	char*		symop = NULL;
	
	if ( spacegroup > 1 ) {
//		if ( strstr(ext, ".star") || strstr(ext, ".cif") ) {
//			symop = read_symop_star(filename, spacegroup, nsym);
//		} else {
			symop = read_symop_lib(filename, spacegroup, nsym);
//		}
		if ( !symop )
			fprintf(stderr, "Error: No symmetry operator file read!\n");
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG read_symop: Symmetry operator file read: %s\n", filename);
	
	if ( default_filename_used ) bfree_string(filename);
	bfree_string(ext);
	
	return(symop);
}

/************************************************************************
@Function: write_symat
@Description:
	Writing crystallographic symmetry operators.
@Algorithm:
	.
@Arguments:
	char* filename		file name.
	int spacegroup		crystal space group number.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
/*int 		write_symat(char* filename, int spacegroup)
{
	int			nsym = 0, err = 0;
	char* 		symop = read_symop(NULL, spacegroup, &nsym);
	if ( !symop || !nsym ) return(-1);
	
	err = write_symop_star(filename, spacegroup, nsym, symop, 80);
	
	bfree(symop, nsym*80*sizeof(char));
	
	return(err);
}
*/
/************************************************************************
@Function: write_pointgroup
@Description:
	Writing crystallographic symmetry operators.
@Algorithm:
	.
@Arguments:
	char* filename		file name.
	int spacegroup		crystal space group number.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
/*int 		write_pointgroup(char* filename, char* symmetry_string, View ref_view)
{
	int			err = 0;
	
	if ( strstr(filename, ".star") )
		err = write_star_pointgroup(filename, symmetry_string, ref_view);
	else {
		fprintf(stderr, "Error: File type for %s not supported!\n", filename);
		err = -1;
	}

	return(err);
}
*/
/*
// Find space group label and operators in a STAR format file
char* 		read_symop_star(char* filename, int spacegroup, int* nsym)
{
	int 		j, notfound = 1, igroup = 0;
	char*		aptr;
	
 	Bstar*		star = init_star();
	
 	if ( read_star(filename, star) < 0 ) {
		kill_star(star);
		return(NULL);
	}
	
	int 		ngroups = item_get_number(star, SYMMETRY_NUMBER, -1);
	int*		groups = item_get_integer_array(star, SYMMETRY_NUMBER, -1);

	// Find the block with the spacegroup
	while ( notfound && igroup < ngroups ) {
		if ( spacegroup == groups[igroup] ) notfound = 0;
		else igroup++;
	}

	bfree(groups, ngroups*sizeof(int));
	
	if ( notfound || igroup >= ngroups ) {
		fprintf(stderr, "Error: Star symmetry operators not found!\n");
		kill_star(star);
		return(NULL);
	}
	
	// Find the number of symmetry operators
	int			n = item_get_number(star, SYMMETRY_EQUIVXYZ, igroup);
	
	// Allocate memory for the symmetry operators
	char*		symop = (char *) balloc(n*80*sizeof(char));
	
	// Find the symmetry operators
	char*		xyz = item_get_string_array(star, SYMMETRY_EQUIVXYZ, igroup);
	aptr = xyz;
	for ( j=0; j<n; j++ ) {
		strcpy(symop+80*j, aptr);
		aptr += strlen(aptr) + 1;
	}
	
	kill_star(star);
	bfree(xyz, aptr - xyz);
	
	*nsym = n;
	
	return(symop);
}
*/
// Find space group label and operators in "symop.lib" or similar file
char* 		read_symop_lib(char* filename, int spacegroup, int* nsym)
{
	FILE*		fsym;
	int 		i, j, k = 0, notfound = 1, number = 0, nlines;
	char		aline[80], symall[4000];
	memset(symall, '\0', 4000);
	
	if ( verbose & VERB_DEBUG )
		printf("Symmetry operator library: %s\n\n", filename);

	if ( strlen(filename) < 1 || ( fsym = fopen(filename, "r") ) == NULL ) {
		if ( getenv("BPARAM") ) {
			strcpy(filename, getenv("BPARAM"));
			if ( filename[strlen(filename)-1] != '/' ) strcat(filename, "/");
		} else
			strcpy(filename, "~/bsoft/parameters/");
		strcat(filename, "symop.lib");
		if ( ( fsym = fopen(filename, "r") ) == NULL ) 
			return(NULL);
	}
	
	if ( verbose & VERB_DEBUG )
		printf("Symmetry operator library: %s\n\n", filename);
	
	while ( notfound && fgets( aline, 80, fsym ) ) {
		sscanf( aline, "%d %d", &number, &nlines );
		if ( number == spacegroup ) notfound = 0;
	}
	
	if ( notfound ) {
		fclose( fsym );
		return(NULL);
	}
	
	int			n = 0;
	for ( i=0; i<nlines; i++ ) {
		fgets( aline, 80, fsym );
		for ( j=0; j<(int)strlen(aline); j++ ) {
			if ( aline[j] != ' ' && aline[j] != '\n' ) {
				symall[k] = tolower(aline[j]);
				k++;
			}
			if ( aline[j] == '*' ) n++;
		}
		symall[k] = '*';
		k++;
		n++;
	}

	if ( verbose & VERB_DEBUG )
		printf("Spacegroup = %d,  Noperators = %d\n", number, n);

	char*		symop = (char *) balloc(n*80*sizeof(char));
	
	j = 0;
	for ( i=0; i<n; i++ ) {
		k = 0;
		while( symall[j] != '*' && k<80 ) {
			symop[i*80+k] = symall[j];
			j++;
			k++;
		}
		j++;
		if ( verbose & VERB_DEBUG )
			printf("%-80s\n", &symop[80*i]);
	}
	
	fclose( fsym );
	
	*nsym = n;
	
	return(symop);
}
/*
int 		write_symop_star(char* filename, int spacegroup, int nsym, char* symop, int line_len)
{
	int			err = 0;
	
 	Bstar*		star = init_star();
	
	int 		j;
	int*		number = (int *) balloc(nsym*sizeof(int));
	char*		astring = (char *) balloc(nsym*line_len);
	char*		aptr = astring;
	
	// Convert the symmetry operators
	item_put_integer(star, SYMMETRY_NUMBER, 0, spacegroup);
	
	for ( j=0; j<nsym; j++ ) {
		number[j] = j;
		strcpy(aptr, symop+line_len*j);
		aptr += strlen(aptr) + 1;
	}
	item_put_integer_array(star, SYMMETRY_EQUIVID, 0, nsym, number);
	item_put_string_array(star, SYMMETRY_EQUIVXYZ, 0, nsym, astring);
			
	err = write_star(filename, star);
	
	if ( err < 0 ) j = err;
	
	kill_star(star);
	bfree(astring, nsym*line_len);
	bfree(number, nsym*sizeof(int));
	
	return(j);
}
*/
/*
int 		write_star_pointgroup(char* filename, char* symmetry_string, View ref_view)
{
	int				err =0;
	
	Bsymmetry*		sym = init_point_group_symmetry(symmetry_string);
	
	Bstar*			star = init_star();
	
	star->comment = copystring("# Symmetry from bsym\n\n");
	
	item_put_string(star, SYMMETRY_POINT_GROUP, 0, sym->label);
	item_put_integer(star, SYMMETRY_PG_NUMBER, 0, sym->point);
	
	int 			i;
	int*			order = (int *) balloc(sym->n*sizeof(int));
	float*			x = (float *) balloc(sym->n*sizeof(float));
	float*			y = (float *) balloc(sym->n*sizeof(float));
	float*			z = (float *) balloc(sym->n*sizeof(float));
	
	Matrix3			mat = matrix3_from_view(ref_view);
	Vector3			new_axis;
	
	for ( i=0; i<sym->n; i++ ) {
		new_axis = vector3_matrix3_multiply(mat, sym->op[i].axis);
		order[i] = sym->op[i].order;
		x[i] = new_axis.x;
		y[i] = new_axis.y;
		z[i] = new_axis.z;
	}
	
	item_put_integer_array(star, SYMMETRY_AXIS_ORDER, 0, sym->n, order);
	item_put_float_array(star, SYMMETRY_AXIS_X, 0, sym->n, x);
	item_put_float_array(star, SYMMETRY_AXIS_Y, 0, sym->n, y);
	item_put_float_array(star, SYMMETRY_AXIS_Z, 0, sym->n, z);
	
	bfree(order, sym->n*sizeof(int));
	bfree(x, sym->n*sizeof(float));
	bfree(y, sym->n*sizeof(float));
	bfree(z, sym->n*sizeof(float));
		
	err = write_star(filename, star);
		
	kill_star(star);
	
	kill_symmetry(sym);
	
	return(err);
}
*/
