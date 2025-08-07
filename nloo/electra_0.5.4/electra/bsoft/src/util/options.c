/* 
	options.c 
	Functions to handle general options 
	Author: Bernard Heymann 
	Created: 20010613	Modified: 20041031
*/ 
 
#include "options.h" 
#include "img_datatypes.h"
#include "matrix.h" 
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

/************************************************************************
@Function: get_option_list
@Description:
	Parses command line arguments based on a template.
@Algorithm:
	The usage list is parsed to find the desired option tag.
	Each option in the usage list must have the following format:
		-<tag> <value,value,...>
	The first character on the line must be '-', the tag may not
	extend beyond the 17'th character and the value must start
	before the 19'th character.
	A partial input tag is tolerated as long as it is unambiguous.
	An ambiguous or unknown tag or a tag requiring a value but
	without one causes program abortion.
	Special tags:
		-verbose 3      sets the verbosity level for all programs.
		-help           returns the usage information and quits.
@Arguments:
	char* use[]			usage list of strings, the template.
	int argc			number of command line arguments.
	char* argv[]		array of command line argument strings.
	int* optind			first argument after option list.
@Returns:
	Option* 			linked list of tag-value pairs.
**************************************************************************/
Option*		get_option_list(char* use[], int argc, char* argv[], int* optind)
{
	if ( argc < 2 ) {
		usage(use);
		exit(-1);
	}
	
	int			maxtaglength = 20;
	int			i, n, v = 1;
	char		usetag[20], tag[256] = "", value[256] = "";
	memset(usetag, 0, maxtaglength);
	Option*		option = NULL;
	Option*		curropt = NULL;
	
	for ( i=1; i<argc; i++ ) {
		if ( strncmp(argv[i], "-verbose", strlen(argv[i])) == 0 && i<argc-1 ) {
			if ( argv[i+1][0] != '-' )
				if ( sscanf(argv[i+1], "%s", value) == 1 )
					verbose = get_option_verbose(value);
		}
		if ( strncmp(argv[i], "-help", strlen(argv[i])) == 0 ) {
			usage(use);
			exit(0);
		}
	}
	
	*optind = 1;
	for ( i=1; i<argc && v; i+=v ) {
		if ( argv[i][0] == '-' ) {
			if ( !option ) option = curropt = (Option *) balloc(sizeof(Option));
			else {
				curropt->next = (Option *) balloc(sizeof(Option));
				curropt = curropt->next;
			}
			for ( n=0; use[n] != NULL; n++ ) {
				strncpy(usetag, use[n], maxtaglength-1);
				if ( usetag[0] == '-' ) {
					if ( strncmp(usetag, argv[i], strlen(argv[i])) == 0 ) {
						if ( strlen(tag) ) {
							fprintf(stderr, "Error: Ambiguous option %s!\n", argv[i]);
							exit(-1);
						} else {
							v = sscanf(usetag, "-%s %s", tag, value);
							if ( verbose & VERB_DEBUG )
								printf("DEBUG get_option_list: usetag = %s  value = %s\n", tag, value);
						}
					}
				}
			}
			if ( strlen(tag) < 1 ) {
				fprintf(stderr, "Error: Option %s is not defined!\n", argv[i]);
				exit(-1);
			}
			strncpy(curropt->tag, tag, 255);
			if ( v > 1 ) {
				if ( i+1 >= argc || strlen(argv[i+1]) < 1 ) {
					fprintf(stderr, "Error: Option -%s must have a value!\n", tag);
					exit(-1);
				}
				strncpy(curropt->value, argv[i+1], 255);
			}
			tag[0] = value[0] = 0;
			*optind = i + v;
			if ( verbose & VERB_DEBUG ) {
				if ( v == 1 ) printf("DEBUG get_option_list: tag = %s\n", curropt->tag);
				else printf("DEBUG get_option_list: tag = %s  value = %s\n", curropt->tag, curropt->value);
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG get_option_list: optind = %d  arg = %s\n", *optind, argv[*optind]);
	
	return(option);
}

/************************************************************************
@Function: option_kill
@Description:
	Deallocates a linked list of option structures.
@Algorithm:
	.
@Arguments:
	Option* option		linked list of tag-value pairs.
@Returns:
	int					0.
**************************************************************************/
int			option_kill(Option* option)
{
	Option		*curropt, *curropt2;
	
	for ( curropt = option; curropt; ) {
		curropt2 = curropt->next;
		bfree(curropt, sizeof(Option));
		curropt = curropt2;
	}
	
	return(0);
}

/************************************************************************
@Function: get_option_verbose
@Description:
	Sets the verbose option.
@Algorithm:
	The verbosity level is defined by the following constants setting 
	particular bits:
		VERB_NONE		0		No output
		VERB_RESULT		1		Program results
		VERB_LABEL		2		Function information
		VERB_PROCESS	4		Selected processing information
		VERB_STATS		8		Statistical information on objects
		VERB_FULL		16		All processing information
		VERB_TIME		32		Timing information
		VERB_MEMORY 	64		Memory allocation and freeing
		VERB_DEBUG		256 	Debugging information
		VERB_DEBUG_STAR	512 	STAR code debugging information
@Arguments:
	char* optarg		verbosity level.
@Returns:
	int 			.
**************************************************************************/
int 		get_option_verbose(char* optarg)
{
	unsigned int 		i;
	
	if ( isdigit(optarg[0]) ) {
		if ( sscanf(optarg, "%d", &verbose) < 1 ) {
			fprintf(stderr, "-verbose: A verbosity level (0-255) must be specified!");
			exit(-1);
		}
	}
	
	for ( i=0; i<strlen(optarg); i++ ) optarg[i] = tolower(optarg[i]);
	if ( strstr(optarg, "non") ) verbose |= VERB_NONE;
	if ( strstr(optarg, "res") ) verbose |= VERB_RESULT;
	if ( strstr(optarg, "lab") ) verbose |= VERB_LABEL;
	if ( strstr(optarg, "pro") ) verbose |= VERB_PROCESS;
	if ( strstr(optarg, "sta") ) verbose |= VERB_STATS;
	if ( strstr(optarg, "ful") ) verbose |= VERB_FULL;
	if ( strstr(optarg, "tim") ) verbose |= VERB_TIME;
	if ( strstr(optarg, "mem") ) verbose |= VERB_MEMORY;
	if ( strstr(optarg, "deb") ) verbose |= VERB_DEBUG;
	if ( strstr(optarg, "star") ) verbose |= VERB_DEBUG_STAR;
 
	return(verbose); 
} 
 
/************************************************************************
@Function: get_option_filename
@Description:
	Sets the desired output data type for images.
@Algorithm:
	Option character:	A capital letter
@Arguments:
	char* optarg		a file name string.
@Returns:
	char*				new file name.
**************************************************************************/
char* 		get_option_filename(char* optarg)
{
	if ( strlen(optarg) < 1 ) {
        fprintf(stderr, "Error: A valid file name must be specified!\n");
		exit(-1);
	}
	
	char*		filename = (char *) balloc((strlen(optarg) + 1)*sizeof(char));
	
	strncpy(filename, optarg, strlen(optarg));
	
	return(filename); 
} 
 
/************************************************************************
@Function: get_option_datatype
@Description:
	Sets the desired output data type for images.
@Algorithm:
	.
@Arguments:
	char* optarg		data type.
@Returns:
	DataType 			new data type.
**************************************************************************/
DataType 	get_option_datatype(char* optarg)
{ 
	DataType 	newdatatype = Unknown_Type;	// Conversion to new type
	char		newtype[20];
	
	if ( sscanf(optarg, "%s", newtype) < 1 ) {
        fprintf(stderr, "-datatype: A data type must be specified!\n");
		exit(-1);
	}
		
	newdatatype = getdatatype(newtype[0]);
	
	return(newdatatype); 
} 

/************************************************************************
@Function: get_option_size
@Description:
	Sets the 3D size for an image.
@Algorithm:
	The default dimension is voxel.
@Arguments:
	char* optarg		1-3 size values.
@Returns:
	VectorInt3 			3 value vector with the voxel dimensions in xyz.
**************************************************************************/
VectorInt3	get_option_size(char* optarg)
{
	int 		i;
	Vector3		fsize = {1,1,1};
	VectorInt3	size = {1,1,1};

	if ( ( i = sscanf(optarg, "%f,%f,%f", &fsize.x, 
			&fsize.y, &fsize.z) ) < 1 ) {
		fprintf(stderr, "-size: At least one size value must be specified!\n");
		exit(-1);
	}
	
	if ( i < 2 ) fsize.z = fsize.y = fsize.x;
	
	size.x = (int) (fsize.x + 0.5);
	size.y = (int) (fsize.y + 0.5);
	size.z = (int) (fsize.z + 0.5);
	
	return(size); 
} 
 
/************************************************************************
@Function: get_option_voxel
@Description:
	Gets a set of coordinates for a voxel in an image.
@Algorithm:
	.
@Arguments:
	char* optarg		3 values.
@Returns:
	int* 				3 value vector with the voxel coordinates in xyz.
**************************************************************************/
int*		get_option_voxel(char* optarg)
{
	int 		i;
	float		fvoxel[3] = {1,1,1};
	int*		voxel = (int *) balloc(3*sizeof(int));

	if ( ( i = sscanf(optarg, "%f,%f,%f", &fvoxel[0], 
			&fvoxel[1], &fvoxel[2]) ) < 1 ) {
		fprintf(stderr, "-voxel: Three coordinates must be specified!\n");
		exit(-1);
	}
	
	for ( i=0; i<3; i++ ) voxel[i] = (int) (fvoxel[i] + 0.5);
	
	return(voxel); 
} 
 
/************************************************************************
@Function: get_option_origin
@Description:
	Gets the 3D origin for an image.
@Algorithm:
	The default dimension is voxel.
@Arguments:
	char* optarg		1-3 size values.
@Returns:
	Vector3 			3 value vector with the origin in xyz.
**************************************************************************/
Vector3		get_option_origin(char* optarg)
{
	Vector3		origin = {0,0,0};

	if ( sscanf(optarg, "%f,%f,%f", &origin.x, 
			&origin.y, &origin.z) < 1 ) {
		fprintf(stderr, "-origin: At least one origin value must be specified!\n");
		exit(-1);
	}
	
	return(origin); 
} 
 
/************************************************************************
@Function: get_option_sampling
@Description:
	Sets the sampling or voxel size for an image.
@Algorithm:
	Sampling is assumed to be in angstrom/voxel edge.
@Arguments:
	char* optarg		1-3 sampling values.
@Returns:
	Vector3 			3 value vector with the sampling values for xyz.
**************************************************************************/
Vector3		get_option_sampling(char* optarg)
{
	int 		i;
	Vector3		sampling = {1,1,1};

	if ( ( i = sscanf(optarg, "%f,%f,%f", &sampling.x, &sampling.y, 
			&sampling.z) ) < 1 ) {
		fprintf(stderr, "-sampling: At least one sampling value must be specified!\n");
		exit(-1);
	}
	
	if ( i < 2 ) sampling.z = sampling.y = sampling.x;
	
	return(sampling); 
} 
 
/************************************************************************
@Function: get_option_vector3
@Description:
	Gets a 3-value vector.
@Algorithm:
	The vector will be normalized.
@Arguments:
	char* optarg		2-4 values.
@Returns:
	Vector3				3-value vector.
**************************************************************************/
Vector3		get_option_vector3(char* optarg)
{
	double		sum = 0;
	Vector3		vec = {0,0,0};

	if ( sscanf(optarg, "%f,%f,%f", &vec.x, &vec.y, &vec.z) < 2 ) {
		fprintf(stderr, "-vector: At least 2 values for the vector must be specified!\n");
		exit(-1);
	}
	
	sum = vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
	
	if ( sum <= 1e-30 ) {
		vec.x = vec.y = 0;
		vec.z = 1;
	} else {
		sum = sqrt(sum);
		vec.x /= sum;
		vec.y /= sum;
		vec.z /= sum;
	}
	
	return(vec); 
} 
 
/************************************************************************
@Function: get_option_view
@Description:
	Gets a view vector and the rotation angle around it.
@Algorithm:
	The vector will be normalized and the angle converted from
	degrees to radians.
@Arguments:
	char* optarg		2-4 values.
@Returns:
	View 				4-value view: 3-value vector and rotation angle.
**************************************************************************/
View		get_option_view(char* optarg)
{
	View		view = {NULL,0,0,0,0};

	if ( sscanf(optarg, "%f,%f,%f,%f", &view.x, &view.y, 
			&view.z, &view.a) < 2 ) {
		fprintf(stderr, "-View: At least 2 values for the vector must be specified!\n");
		exit(-1);
	}
	
	view.a *= PI/180;
	while ( view.a <= -PI ) view.a += 2*PI;
	while ( view.a >   PI ) view.a -= 2*PI;
	
	view_normalize(&view);
	
	return(view); 
} 
 
/************************************************************************
@Function: get_option_fill_value
@Description:
	Sets the fill value for various operations involving adding or replacing data.
@Algorithm:
	The fill value can be given as a value, or as the image average or background.
	Valid type designators are:
		back(ground) or bg or bkg
		av(erage) or avg
@Arguments:
	char* optarg		a number or fill type designator.
	int* fill_type		enumerated fill type.
@Returns:
	float 				fill value, default zero.
**************************************************************************/
float		get_option_fill_value(char* optarg, int* fill_type)
{
	unsigned int	i;
	float			fill = 0;
	
	*fill_type = FILL_USER;

	if ( sscanf(optarg, "%f", &fill) > 0 ) return(fill);

	for ( i=0; i<strlen(optarg); i++ ) optarg[i] = tolower(optarg[i]);
	if ( optarg[0] == 'a' ) *fill_type = FILL_AVERAGE;
	if ( optarg[0] == 'b' ) *fill_type = FILL_BACKGROUND;
	if ( *fill_type == FILL_USER ) {
		fprintf(stderr, "-fill: A fill value or type must be specified!\n");
		exit(-1);
	}
	
	return(fill); 
} 
 
/************************************************************************
@Function: get_option_symmetry
@Description:
	Sets the symmetry point group string.
@Algorithm:
	Valid symmetry strings are:
		C1, C2, C3, ... 	1, 2, 3, ...
		D3, D4, D5, ... 	32, 422, 52, ...
		T					23
		O					432
		I					523
@Arguments:
	char* optarg		symmetry.
@Returns:
	char* 				symmetry point group string.
**************************************************************************/
char*		get_option_symmetry(char* optarg)
{
	if ( strlen(optarg) < 1 ) {
		fprintf(stderr, "-symmetry: A symmetry must be specified!\n");
		exit(-1);
	}
	
	char		pgtemp[20] = "";
	
	strncpy(pgtemp, optarg, 10);
	
	symmetry_clean_pointgroup(pgtemp);
	
	char*		symmetry_string = (char *) balloc((strlen(pgtemp)+1)*sizeof(char));
	
	strcpy(symmetry_string, pgtemp);
	
	return(symmetry_string);
} 
 
/************************************************************************
@Function: get_option_unit_cell
@Description:
	Gets crystal unit cell parameters.
@Algorithm:
	Unit cell angles in the argument are assumed to be in degrees and
	will be converted to radians in the structure returned.
@Arguments:
	char* optarg		6 values.
@Returns:
	UnitCell 			4-value view: 3-value vector and rotation angle.
**************************************************************************/
UnitCell		get_option_unit_cell(char* optarg)
{
	UnitCell		unit_cell = {0,0,0,PI/2,PI/2,PI/2};

	if ( sscanf(optarg, "%f,%f,%f,%f,%f,%f", &unit_cell.a, &unit_cell.b, 
			&unit_cell.c, &unit_cell.alf, &unit_cell.bet, &unit_cell.gam) < 6 ) {
		fprintf(stderr, "-unitcell: All 6 values for the unit cell must be specified!\n");
		exit(-1);
	}
	
	unit_cell.alf *= PI/180;
	unit_cell.bet *= PI/180;
	unit_cell.gam *= PI/180;
	unit_cell.alf = angle_set_negPI_to_PI(unit_cell.alf);
	unit_cell.bet = angle_set_negPI_to_PI(unit_cell.bet);
	unit_cell.gam = angle_set_negPI_to_PI(unit_cell.gam);
	
	return(unit_cell); 
} 
 
