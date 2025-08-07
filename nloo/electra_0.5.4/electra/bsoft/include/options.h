/* 
	options.h 
	Header file for option handlers 
	Author: Bernard Heymann 
	Created: 20010613 	Modified: 20030520
*/

#include "rwimg.h"
#include "matrix.h"
#include "symmetry.h" 
 
#define FILL_USER		0
#define FILL_AVERAGE	1
#define FILL_BACKGROUND 2

struct Option {
	Option*	next;
	char	tag[256];
	char	value[256];
} ;
	 
// Function prototypes 
Option*		get_option_list(char* use[], int argc, char* argv[], int* optind);
int			option_kill(Option* option);
int 		get_option_verbose(char* optarg);
char* 		get_option_filename(char* optarg);
DataType 	get_option_datatype(char* optarg);
VectorInt3	get_option_size(char* optarg);
int*		get_option_voxel(char* optarg);
Vector3		get_option_origin(char* optarg);
Vector3		get_option_sampling(char* optarg);
Vector3		get_option_vector3(char* optarg);
View		get_option_view(char* optarg);
float		get_option_fill_value(char* optarg, int* fill_type);
char*		get_option_symmetry(char* optarg);
UnitCell	get_option_unit_cell(char* optarg);
