/* 
	utilities.h 
	Header file for general utilities 
	Author: Bernard Heymann 
	Created: 19990722 	Modified: 20041224
*/ 

#define BVERSION "1.2.6-20050309"
 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <ctype.h> 
#include <unistd.h> 
#include <limits.h> 
#include <time.h>
#include <math.h> 
#ifdef SUNOS
#include <ieeefp.h>
#endif

#define SMALLFLOAT  1e-37
 
//#define ERROR error_show((char *)__FUNCTION__, __FILE__, __LINE__)
//#define ERRET return(ERROR)

/************************************************************************
@Constant: PI
@Description:
	Mathematical PI and 2*PI.
*************************************************************************/
#ifndef PI 
#define PI		3.14159265358979323846 
#endif 
#ifndef TWOPI 
#define TWOPI	6.283185307179586
#endif 
#ifndef MIN2PI 
#define MIN2PI	-6.283185307179586
#endif
//end

/************************************************************************
@Constant: RHO
@Description:
	Protein density.
*************************************************************************/
#ifndef RHO 
#define RHO 		0.74    	// Protein density in Da/A3
#endif 
//end

#define IMGLABELSIZE	1024	// Size of image label in the Bimage struct

/************************************************************************
@Constants: VERB_
@Description:
	Verbosity designations - controlling output to stdout.
@Features:
	The notion is that any program should not automatically generate output,
	so that it could be used within a script. Increasing levels of verbosity
	then results in increasing output.
*************************************************************************/
#define VERB_NONE		0
#define VERB_RESULT		1
#define VERB_LABEL		2
#define VERB_PROCESS	4
#define VERB_STATS		8
#define VERB_FULL		16
#define VERB_TIME		32
#define VERB_MEMORY 	64
#define VERB_DEBUG		256
#define VERB_DEBUG_STAR	512
//end

#ifndef _systype_ 
/************************************************************************
@Object: enum SysType
@Description:
	System type enumeration.
@Features:
	BigIEEE must be the lowest number type 
	LittleIEEE must be the little-endian type with the lowest number 
*************************************************************************/
enum SysType { 
	Unknown_System = 0,	// Indeterminate 
	BigIEEE = 1,		// Big-endian IEEE (unix: MIPS, Motorola PPC) 
	BigOther = 5,		// Big-endian systems other than IEEE 
	LittleIEEE = 6, 	// Little-endian IEEE (unix, Intel, VMS Alpha) 
	LittleAlpha = 7,	// Little-endian VMS non-IEEE (Alpha) 
	LittleVAX = 8,		// Little-endian VMS (VAX) 
	LittleOther = 10	// Little-endian other systems 
} ; 
#define _systype_ 
#endif 
 
#ifndef _complex_
/************************************************************************
@Object: struct complex_short
@Description:
	Complex number of two two-byte signed integers.
*************************************************************************/
struct complex_short  { short re,im; } ;

/************************************************************************
@Object: struct complex_int
@Description:
	Complex number of two four-byte signed integers.
*************************************************************************/
struct complex_int    { int   re,im; } ;

/************************************************************************
@Object: struct complex_float
@Description:
	Complex number of two four-byte floating point numbers.
*************************************************************************/
struct complex_float  { float re,im; } ;

/************************************************************************
@Object: struct complex_double
@Description:
	Complex number of two eight-byte floating point numbers.
*************************************************************************/
struct complex_double { double re,im; } ;

/************************************************************************
@Object: struct polar
@Description:
	Polar number with two four-byte floating point numbers: amplitude and phase.
*************************************************************************/
struct polar          { float amp,phi; } ;
#define _complex_
#endif

// Function prototypes 
int 		get_cmd_line(char* cmd_line, int argc, char **argv);
void		usage(char** use);
SysType 	systype(int show); 
char*		balloc(unsigned long size); 
int 		bfree(void* ptr, unsigned long size);
int			bexit(int value);
int 		bfree_string(char* string);
int 		bfree_stringlist(char* string, int n);
int 		bfree_stringarray(char** array, long n);
int 		bfreenull(void **ptr, unsigned long size);
long 		get_rand_max();
int 		findNextPowerOf(int startNumber, int powerOf);
int 		swap_integers(int* i1, int* i2);
int 		swap_unsigned_integers(unsigned int* i1, unsigned int* i2);
int 		swap_longs(long* i1, long* i2);
int 		swap_unsigned_longs(unsigned long* i1, unsigned long* i2);
int 		swap_floats(float* f1, float* f2);
int 		swap_doubles(double* f1, double* f2);
float		angle_set_negPI_to_PI(float angle) ;
float		remove_negative_zeros(float value0, float threshold);
void		swapbytes(char* v, int n);
void		vax2ieee(char* v, int swap); 
int			detect_and_fix_carriage_return(char* filename);
int			error_show(char* message, char* file, int line);

