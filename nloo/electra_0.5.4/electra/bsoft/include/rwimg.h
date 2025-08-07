/*
	rwimg.h
	Header file for 2D and 3D image I/O
	Author: Bernard Heymann
	Created: 19990321 	Modified: 20050203
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>

#include "matrix.h"
#include "string_util.h"

#define SWAPTRIG	65535	// Threshold file z size above which bytes are swapped

#ifndef _transform_
/************************************************************************
@Object: enum Transform
@Description:
	Fourier transform format specifier.
@Features:
	Transforms are classified according to where the origin is and whether
	all or only the hermitian half is stored in a file:
*************************************************************************/
enum Transform {
	NoTransform = 0, 	// No transform
	Standard = 1,		// Standard transform: origin = (0,0,0)
	Centered = 2,		// Centered transform: origin = (nx/2,ny/2,nz/2)
	Hermitian = 3,		// Hermitian half: origin = (0,0,0)
	CentHerm = 4		// Centered hermitian: origin = (0,ny/2,nz/2)
} ;
#define _transform_
#endif

#ifndef _Bimage_
/************************************************************************
@Object: enum DataType
@Description:
	Data type specifier.
@Features:
	This determines what datatype is used in an image.
*************************************************************************/
enum DataType {
	Unknown_Type = 0,	// Undefined data type
	UChar = 1,			// Unsigned character or byte type
	SChar = 2,			// Signed character (for CCP4)
	UShort = 3,			// Unsigned two-byte
	Short = 4,			// Signed two-byte
	Int = 5,			// Signed integer (4-byte)
	Long = 6,			// Signed integer (4 or 8 byte, depending on system)
	Float = 7,			// Floating point (4-byte)
	Double = 8,			// Double precision floating point (8-byte)
	ComplexShort = 9,	// Complex two-byte integer (4-byte)
	ComplexInt = 10,	// Complex integer (8-byte)
	ComplexFloat = 11,	// Complex floating point (8-byte)
	ComplexDouble = 12,	// Complex floating point (8-byte)
	Polar = 13			// Polar floating point (8-byte)
} ;
	
/************************************************************************
@Object: enum ColorModel
@Description:
	Color model specifier.
@Features:
	This determines what is packed into the channels of every image.
*************************************************************************/
enum ColorModel {
	Bit = 1,			// Bitmap
    Gray = 2,			// Gray scale
    RGB = 3, 			// Red-green-blue interleaved
	RGBA = 4,			// Red-green-blue-alpha interleaved
	CMYK = 5,			// Cyan-magenta-yellow-black interleaved
	Index = 6			// Indexed color table
} ;
	
/************************************************************************
@Object: struct Bsub_image
@Description:
	General sub-image parameter structure.
@Features:
	This contains all the information pertinent to a single image in a
	multi-image file.
*************************************************************************/
struct Bsub_image {
	float		background; 	// Image background
	float		ox, oy, oz; 	// Origin (position of object origin in voxel coordinates)
	float		vx, vy, vz; 	// View orientation unit vector
	float		angle;			// Rotation around view vector in radians
} ;

/************************************************************************
@Object: struct Bimage
@Description:
	General image parameter structure.
@Features:
	All floating point coordinates are in angstroms and converted 
	using the voxel units parameters.
*************************************************************************/
struct Bimage {
	Bimage*			next;			// Pointer for linked lists
	char			filename[256];	// File name
	time_t			time;			// Time in seconds since 00:00:00 January 1, 1970, (UTC)
	unsigned int	dataflag;		// Flag to force reading of the data
	unsigned int	fomflag;		// Flag to indicate presence of FOM block
	unsigned long	x, y, z, c; 	// Dimensions, xyz and channels
	unsigned long	n, i;			// Number of images and image number (may be > n)
	unsigned long	px, py, pz; 	// Page dimensions
	unsigned long	offset; 		// Data offset
	DataType		datatype;		// Data type
	Transform		transform;  	// Transform type
	ColorModel		colormodel; 	// Gray, RGB
	unsigned int	colors; 		// Number of colours in map
	char*			colormap;		// Colour map for indexed images
	float			min, max;		// Limits
	float			avg, std;		// Average and standard deviation
	float			smin, smax; 	// Limits for display
	float			scale;			// Scale of last conversion operation
	float			resolution; 	// Resolution limit of data - used for low-pass filtering
	float			em_volt;		// Electron microscope acceleration voltage (V)
	float			em_Cs;			// Electron microscope spherical aberration coefficient (angstrom)
	float			em_amp; 		// Electron microscope amplitude contrast contribution (fraction)
	float			em_def_avg; 	// Electron microscope defocus average (angstrom)
	float			em_def_dev; 	// Electron microscope defocus deviation (angstrom)
	float			em_ast_ang; 	// Electron microscope astigmatism angle (radian)
	unsigned int	em_basetype;	// Electron microscope baseline type (1=poly, 2=double_gauss, 3=eman)
	float			em_base[5]; 	// Electron microscope baseline (up to 5 coefficients)
	float			em_env[5];		// Electron microscope envelope (up to 5 coefficients)
	float			ux, uy, uz;		// Voxel units (angstrom/pixel edge)
	float			ua, ub, uc; 	// Unit cell dimensions (angstrom)
	float			alf, bet, gam;	// Unit cell angles (radian)
	unsigned int	spacegroup;		// Space group
	char			spacelabel[80];	// Space group label
	char*			data;			// Pointer to the data
    float			fommax; 		// Maximum FOM
	float*			fom;			// Figures of merit
	Bsub_image*		image;			// Sub-images
    char*			label;			// Label block
} ;
#define _Bimage_
#endif

// Function prototypes
Bimage* 	init_img();
Bimage* 	init_img_header(DataType datatype, unsigned long c, unsigned long x, 
				unsigned long y, unsigned long z, unsigned long n);
Bimage* 	init_img_with_parameters(DataType datatype, unsigned long c, unsigned long x, 
				unsigned long y, unsigned long z, unsigned long n);
float*		img_init_fom(Bimage* p);
Bimage* 	read_img(char* filename, int readdataflag, int select);
char*		img_read_data(FILE* fimg, Bimage* p, int select, int swap, int vax, int pad);
int 		write_img(char* filename, Bimage* p);
int			write_fom_as_image(char* filename, Bimage* p);
int			write_data_block_as_image(char* filename, char* data, DataType datatype,
				unsigned long c, unsigned long x, unsigned long y, unsigned long z, unsigned long n);
int			img_write_data(FILE* fimg, char* data, unsigned long datasize, 
				unsigned long typesize, int swap);
Bimage* 	copy_img_header(Bimage* p, int new_nimg);
Bimage* 	copy_img (Bimage* p);
int 		move_img (Bimage* pfrom, Bimage* pto);
int 		img_slices_to_images(Bimage* p);
int 		img_images_to_slices(Bimage* p);
int 		kill_img(Bimage* p);
int 		kill_all_but_main_img(Bimage* p);
int 		img_clear_data(Bimage* p);
int 		img_check_param(Bimage* p);
int 		img_stats(Bimage *p);
int 		img_stats_within_radii(Bimage *p, float rad_min, float rad_max);
int  		img_set_origin(Bimage* p, Vector3 origin);
int  		img_set_sampling(Bimage* p, Vector3 sampling, int uc_flag);
int  		img_set_isotropic_sampling(Bimage* p, float sampling, int uc_flag);
float		img_max_radius(Bimage* p);
int			img_convert_fourier(Bimage* p, Transform newtransform);
int			img_info(Bimage* p);
int			sub_img_info(Bimage* p);
