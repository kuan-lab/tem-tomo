/* 
	utilities.c 
	Library functions useful in all the package 
	Author: Bernard Heymann 
	Created: 19990722	Modified: 20041225
*/ 

#include <errno.h>
#include "utilities.h" 
 
// Definition of the global variables 
int 	verbose = 0;		// Default: No output to the screen
long 	memory = 0;			// Total memory allocated 

struct MemBlock {
	MemBlock*		next;
	char*			ptr;
	unsigned long	size;
} ;

// Globals in the local scope to follow memory allocations and deallocations
MemBlock*	mb = NULL;
MemBlock*	mlist = NULL;

/************************************************************************
@Function: get_cmd_line
@Description:
	Returns the command line as a string.
@Algorithm:
	Concatenates the command line arguments into one string. The string
	receiving the command line must be at least 1024 bytes long.
@Arguments:
	char* cmd_line	the string receiving the command line.
	int argc		the number of command line arguments.
	char **argv 	the command line arguments.
@Returns:
	int 			the number of characters in the command line string.
**************************************************************************/
int 		get_cmd_line(char* cmd_line, int argc, char **argv)
{ 
	int			i; 
	
	memset(cmd_line, 0, strlen(cmd_line));
	for ( i=0; i<argc && strlen(cmd_line)+strlen(argv[i]) < IMGLABELSIZE; i++ ) {
		 strcat(cmd_line, argv[i]);
		 strcat(cmd_line, " ");
	}
	
	if ( verbose & VERB_DEBUG )
		printf("%s\n", cmd_line);
	 
	return((int)strlen(cmd_line)); 
} 
 
/************************************************************************
@Function: usage
@Description:
	Prints usage information.
@Algorithm:
	The usage information must be written into an array of srings, with each
	string a line and following a specific convention for the Bsoft package.
	The first line with non-space characters must start with "Usage:" 
	followed by the command-line syntax. The next lines should describe the
	program. Finally, a line starting with "Options:" indicates the start
	of a list of command-line options and their descriptions.
@Arguments:
	char** use		the string array.
@Returns:
	void			-.
**************************************************************************/
void		usage(char** use)
{
    int     i;

    for ( i = 0; use[i] != NULL; i++ )
	    printf("%s\n", use[i]);

	printf("\nBsoft version %s (%ld bit)\n\n", BVERSION, 8*sizeof(long));
	
    exit(-1);
}

/************************************************************************
@Function: systype
@Description:
	Finds the system type - mostly just for byte order.
@Algorithm:
	Test the byte order of an arbitrary byte sequence by interpreting it as
	an integer or a floating point number.
@Arguments:
	int show		a flag to indicate if the result should be shown.
@Returns:
	SysType			an enumerated type.
**************************************************************************/
SysType 	systype(int show)
{ 
	SysType 	type = Unknown_System; 
	char*		test = (char *) balloc(12); 
	int*		itest = (int *) test; 
	float*		ftest = (float *) test; 
	 
	memcpy(test, "jbh     ", 8); 
	 
	if ( itest[0] == 1784834080 && fabs(ftest[0] - 6.84272e+25) < 1e+21 ) 
		type = BigIEEE; 
	 
	if ( itest[0] == 543711850 && fabs(ftest[0] - 1.96837e-19) < 1e-23 ) 
		type = LittleIEEE; 
	 
	if ( ( ( verbose > VERB_PROCESS ) && show ) || ( verbose & VERB_DEBUG ) ) { 
		switch ( type ) { 
			case BigIEEE: 
				printf("System type:                    Big-endian IEEE (%ld)\n\n",
						8*sizeof(long)); 
				break; 
			case LittleIEEE: 
				printf("System type:                    Little-endian IEEE (%ld)\n\n", 
						8*sizeof(long)); 
				break; 
			default: 
				printf("System type:                    Unknown (%ld)\n\n", 
						8*sizeof(long)); 
				break; 
		} 
	} 
 
	if ( verbose & VERB_DEBUG )
		printf("DEBUG systype: %4s %d %g %d\n", test, itest[0], ftest[0], type); 
	 
	bfree(test, 12);
	
	return(type); 
} 
 
/************************************************************************
@Function: balloc
@Description:
	Allocates memory and clears it to zeros.
@Algorithm:
	It is called exactly like malloc, with the following enhancements:
	Successfully allocated memory is zeroed and the memory counter updated.
	If allocation of zero bytes are requested it notifies the user.
	Allocation is attempted and an error message is printed on failure.
	All failures return a NULL pointer to allow error handling from
	calling functions.
@Arguments:
	unsigned long memsize	size of memory to be allocated.
@Returns:
	char*					a pointer to the memory (NULL on failure).
**************************************************************************/
char*		balloc(unsigned long memsize) 
{ 
	char*		ptr = NULL;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG balloc: Memory requested = %ld\n", memsize);
	 
	if ( memsize == 0 ) {
		fprintf(stderr, "Error: Memory allocation size requested is zero!\n");
		return(NULL);
	}
	
//	if ( memsize < 0 ) {
//		if ( sizeof(long) < 8 ) {
//			fprintf(stderr, "Error: Memory allocation of %ld bytes not done!\n", memsize);
//			fprintf(stderr, "	The code was compiled as 32 bit and therefore cannot\n");
//			fprintf(stderr, "	allocated memory larger than 2Gb total.\n\n");
//		} else {
//			fprintf(stderr, "Error: Memory allocation of %ld bytes not done!\n", memsize);
//			fprintf(stderr, "	Reason for negative size unknown.\n\n");
//		}
//		return(NULL);
//	} 
	
	if ( memory == 0 && ( verbose & VERB_MEMORY ) )
		printf("Memory:    Allocated        Freed        Total\n");
	
	if ( ( ptr = (char *) malloc(memsize*sizeof(char)) ) == NULL ) { 
		printf("\nError: Memory allocation of %ld bytes failed!\n", memsize); 
		printf("Current total memory allocation: %12ld\n\n", memory);
		return(NULL); 
	}
	
	memory += memsize;
	 
	if ( verbose & VERB_MEMORY ) {
		printf("MEMORY: %12ld              %12ld\n", memsize, memory);
		// Add a memory block item to the memory list
		if ( mb ) {
			if ( mb->next ) {
				fprintf(stderr, "Error: Last memory block has a non-zero next pointer! (%p)\n", mb->next);
				bexit(-11);
			}
			mb->next = (MemBlock *) malloc(sizeof(MemBlock));
			mb = mb->next;
		} else {
			mb = (MemBlock *) malloc(sizeof(MemBlock));
		}
		memset(mb, 0, sizeof(MemBlock));
		if ( !mlist ) mlist = mb;
		mb->size = memsize;
		mb->ptr = ptr;
	}
	
	memset(ptr, 0, memsize); 
	 
	return(ptr); 
}

/************************************************************************
@Function: bfree
@Description:
	Frees allocated memory.
@Algorithm:
	It is called exactly like free, with the following enhancements:
	If freeing fails an error message is printed.
	The memory counter is updated on successful freeing.
@Arguments:
	void* ptr				pointer to memory to be freed.
	unsigned long memsize	size of memory block pointed to.
@Returns:
	int						0 = success, -1 = failure.
**************************************************************************/
int 		bfree(void* ptr, unsigned long memsize)
{
	if ( ptr == NULL ) return(0);
	
	if ( memsize < 1 ) {
		if ( verbose & VERB_MEMORY )
			fprintf(stderr, "Error: Freeing memory of %ld bytes failed!\n", memsize);
		return(-1);
	}
	
	free(ptr);
	
	memory -= memsize;
	
	int				error = 0;
	MemBlock		*mba, *mbb;
	if ( verbose & VERB_MEMORY ) {
		printf("MEMORY:              %12ld %12ld\n", memsize, memory);
		// Find and remove a memory block item from the memory list
		if ( !mlist ) {
			error = 1;
		} else if ( mlist->ptr == ptr ) {	// Removes the first item
			mba = mlist;
			mlist = mlist->next;
			if ( mb == mba ) mb = mlist;		// Sets the end of the list
			free(mba);
		} else {
			error = 2;
			for ( mba=mlist; mba->next; mba=mba->next ) {
				mbb = mba->next;
				if ( mbb->ptr == ptr ) {
					if ( mb == mbb ) mb = mba;	// Sets the end of the list
					mba->next = mbb->next;
					free(mbb);
					error = 0;
					break;
				}
			}
		}
		if ( error )
			fprintf(stderr, "Error: Memory block of size %ld not found for deallocation! (%d)\n", memsize, error);
	}
	
	return(0);
}

/************************************************************************
@Function: bexit
@Description:
	Exit function for cleanup and debugging.
@Algorithm:
	If the verbose memory flag is turned on, it shows all non-deallocated blocks of memory.
@Arguments:
	int value		exit value.
@Returns:
	int				given value.
**************************************************************************/
int			bexit(int value)
{
	unsigned long		n;
	
	if ( verbose & VERB_MEMORY ) {
		fprintf(stderr, "Memory blocks not deallocated:\n");
		for ( n=0, mb=mlist; mb; mb=mb->next, n++ )
//			fprintf(stderr, "%ld\n", mb->size);
			fprintf(stderr, "%p:\t%ld\n", mb, mb->size);
		fprintf(stderr, "Number:\t%ld\n", n);
	}
	
	exit(value);
	return(value);
}

/************************************************************************
@Function: bfree_string
@Description:
	Frees allocated memory for a string.
@Algorithm:
	The memory allocated for the string is assumed to be the length
	plus one.
@Arguments:
	char* string	pointer to string.
@Returns:
	int				0.
**************************************************************************/
int 		bfree_string(char* string)
{
	if ( !string ) return(0);
	
	bfree(string, (strlen(string) + 1)*sizeof(char));
	
	return(0);
}

/************************************************************************
@Function: bfree_stringlist
@Description:
	Frees allocated memory for a list of strings.
@Algorithm:
	The memory for a list of contiguous zero-terminated strings is deallocated.
	The memory allocated for each string is assumed to be the length
	plus one.
@Arguments:
	char* string	pointer to first string.
	int n			number of strings.
@Returns:
	int				0.
**************************************************************************/
int 		bfree_stringlist(char* string, int n)
{
	if ( !string ) return(0);
	
	long		i, size = 0;
	char*		aptr = string;
	
	for ( i=0; i<n; i++, aptr += strlen(aptr)+1 ) size += strlen(aptr)+1;
	
	bfree(string, size);
	
	return(0);
}

/************************************************************************
@Function: bfree_stringarray
@Description:
	Frees allocated memory for an array of strings.
@Algorithm:
	The memory for an array of strings is deallocated.
	The memory allocated for each string is assumed to be the length
	plus one.
@Arguments:
	char** array	pointer to array.
	long n			number of strings.
@Returns:
	int				0.
**************************************************************************/
int 		bfree_stringarray(char** array, long n)
{
	if ( !array ) return(0);
	
	long		i;
	
	for ( i=0; i<n; i++ )
		bfree_string(array[i]);
		
    bfree(array, n*sizeof(char*));

	return(0);
}

/************************************************************************
@Function: bfreenull
@Description:
	Frees allocated memory and set the pointer to null.
@Algorithm:
	It is called exactly like free, with the following enhancements:
	If freeing fails an error message is printed.
	The memory counter is updated on successful freeing.
	The pointer is set to NULL to avoid subsequent attempts to reference it.
@Arguments:
	void** ptr		double-referenced pointer to memory to be freed.
	unsigned long size		size of memory block pointed to.
@Returns:
	int				0 = success, -1 = failure.
**************************************************************************/
int 		bfreenull(void **ptr, unsigned long size)
{
	if ( !ptr ) return(-1);
	
	bfree(*ptr, size);
	
	*ptr = NULL;
	
	return(0);
}

/************************************************************************
@Function: get_rand_max()
@Description:
	Finds the maximum random number for a system.
@Algorithm:
	Loops through random numbers to determine if the maximum is 2 or 4 bytes.
@Returns:
	long 			the maximum random number.
**************************************************************************/
long 		get_rand_max()
{
	long 		i, rand_max = 32767;
	
	for ( i=0; i<100; i++ )
		if ( random() > rand_max ) rand_max = 2147483647;
	
	return(rand_max);
}

/************************************************************************
@Function: findNextPowerOf
Author: Dan Krainak
@Description:
	Finds the next greatest number that is a power of a given number.
@Algorithm:
	Loop through the powerOf variable, multiplying it each successive
	iteration until it is greater than the starting number.
	Eg., the next greatest power of 2 starting at 100 is 128.
@Arguments:
	int	startNumber	number to begin from.
	int	powerOf		power of this number is the number returned.
@Returns:
	int 			the next greatest power (i.e. 128)
					or 0 if error or startNumber was not positive int.
***********************************************************************/
int 		findNextPowerOf(int startNumber, int powerOf)
{
	int	rVal = powerOf;
	
	if (powerOf <= 0 ) {
		rVal = 0; 
	}
	else if (powerOf == 1) {
		rVal = startNumber;
	}
	else {
		while (rVal < startNumber) {
			rVal *= powerOf;
		}
	} 

	return rVal;
} //end findNextPowerOf

/************************************************************************
@Function: swap_integers
@Description:
	Swaps two integers.
@Algorithm:
	.
@Arguments:
	int*	i1		pointer to first integer.
	int*	i2		pointer to second integer.
@Returns:
	int 			0.
************************************************************************/
int 		swap_integers(int* i1, int* i2)
{
	int 		itemp;
	
	itemp = *i1;
	*i1 = *i2;
	*i2 = itemp;
	
	return(0);
}

/************************************************************************
@Function: swap_unsigned_integers
@Description:
	Swaps two unsigned integers.
@Algorithm:
	.
@Arguments:
	unsigned int*	i1		pointer to first integer.
	unsigned int*	i2		pointer to second integer.
@Returns:
	int 			0.
************************************************************************/
int 		swap_unsigned_integers(unsigned int* i1, unsigned int* i2)
{
	unsigned int	itemp;
	
	itemp = *i1;
	*i1 = *i2;
	*i2 = itemp;
	
	return(0);
}

/************************************************************************
@Function: swap_longs
@Description:
	Swaps two long integers.
@Algorithm:
	.
@Arguments:
	long*	i1		pointer to first integer.
	long*	i2		pointer to second integer.
@Returns:
	int 			0.
************************************************************************/
int 		swap_longs(long* i1, long* i2)
{
	long 		itemp;
	
	itemp = *i1;
	*i1 = *i2;
	*i2 = itemp;
	
	return(0);
}

/************************************************************************
@Function: swap_unsigned_longs
@Description:
	Swaps two unsigned long integers.
@Algorithm:
	.
@Arguments:
	unsigned long*	i1		pointer to first integer.
	unsigned long*	i2		pointer to second integer.
@Returns:
	int 			0.
************************************************************************/
int 		swap_unsigned_longs(unsigned long* i1, unsigned long* i2)
{
	unsigned long 		itemp;
	
	itemp = *i1;
	*i1 = *i2;
	*i2 = itemp;
	
	return(0);
}

/************************************************************************
@Function: swap_floats
@Description:
	Swaps two floating point numbers.
@Algorithm:
	.
@Arguments:
	float*	i1		pointer to first floating point.
	float*	i2		pointer to second floating point.
@Returns:
	int 			0.
************************************************************************/
int 		swap_floats(float* f1, float* f2)
{
	float 		temp;
	
	temp = *f1;
	*f1 = *f2;
	*f2 = temp;
	
	return(0);
}

/************************************************************************
@Function: swap_doubles
@Description:
	Swaps two double precision floating point numbers.
@Algorithm:
	.
@Arguments:
	double*	i1		pointer to first double precision floating point.
	double*	i2		pointer to second double precision floating point.
@Returns:
	int 			0.
************************************************************************/
int 		swap_doubles(double* f1, double* f2)
{
	double 		temp;
	
	temp = *f1;
	*f1 = *f2;
	*f2 = temp;
	
	return(0);
}

/************************************************************************
@Function: angle_set_negPI_to_PI
@Description:
	Returns an angle between -PI and PI.
@Algorithm:
	.
@Arguments:
	float angle			input angle.
@Returns:
	float				angle between -PI and PI.
**************************************************************************/
float		angle_set_negPI_to_PI(float angle) 
{
	while ( angle <= -PI ) angle += TWOPI;
	while ( angle >   PI ) angle -= TWOPI;
	
	return(angle);
}

/************************************************************************
@Function: remove_negative_zeros
@Author:   David Belnap
@Description:
	Prevent a negative sign from being placed in front of zero value in
	a text file.
@Algorithm:
	This function is intended to be used when obvious zero values are set
	to a very small negative number.  Input value is set to zero if

	              value0 > threshold  and  value0 < 0

	If so, the value is reset to zero.  Otherwise, the input value is 
	returned.
@Arguments:
	float   value0      input value to be tested
	float   threshold   a small negative number
@Returns:
	float   value       "Corrected" or input value
**************************************************************************/
float   remove_negative_zeros(float value0, float threshold)
{
	float  value;
	
	if ( (value0 > threshold) && (value0 < 0) )
		value = 0;
	else
		value = value0;

	return(value);
}

/************************************************************************
@Function: swapbytes
@Description:
	Swaps bytes.
@Algorithm:
	Byte swapping is done in place. 
@Arguments:
	char* v 			a pointer to the bytes.
	int n				the number of bytes to swap.
@Returns:
	void				-.
**************************************************************************/
void		swapbytes(char* v, int n)
{
	char		t;
	for ( int i=0; i<n/2; i++ ) {
		t = v[i];
		v[i] = v[n-1-i];
		v[n-1-i] = t;
	}
}

/************************************************************************
@Function: vax2ieee
Source: Derived from CCP4 code
@Description:
	Converts VAX floating point format to IEEE floating point format.
@Algorithm:
	Swap bytes prior to conversion if the swap flag is set.
	Handle special cases of zero, infinity, NaN or normalized values
	Otherwise assign the new byte values
@Arguments:
	char* v 		four-byte array holding the floating point value
	int	swap		flag to swap bytes before conversion
**************************************************************************/
void		vax2ieee(char* v, int swap) 
{ 
    char    	t[4]; 
    int     	shift,expnt; 
 
    if ( swap ) { 
		t[0] = v[2];  					// Swap byte order 
		t[1] = v[3]; 
		t[2] = v[0]; 
		t[3] = v[1]; 
    } else {
	    t[0] = v[1];  					// Swap for VAX format 
	    t[1] = v[0]; 
	    t[2] = v[3]; 
	    t[3] = v[2];
	}
    memcpy(v,t,4); 
    expnt = (v[0] << 1) | (v[1] >> 7);	// Get exponent 
    if ( expnt == 0 ) 
    	if ( v[0] == 0 )				// Special case zero 
    	    memset(t, '\0', 4); 
    	else {			    			// Special case infinity or NaN 
	    	t[0] = 0xff; 
	    	t[1] = v[1] | 0x80; 
    	} 
    else if ( expnt > 2 )		    	// Normalized value 
    	t[0] = v[0] - 0x01;				// Subtract 2 from the exponent 
    else {								// Denormalized value 
    	shift = 3 - expnt; 
    	t[0] = v[0] & 0x80;	    		// Keep sign, expnt = 0  
    	t[1] = ( (v[1] & 0x7f) >> shift ) | ( 0x10 << expnt ); 
		t[2] = ( v[1] << ( 8-shift ) ) | ( v[2] >> shift ); 
		t[3] = ( v[2] << ( 8-shift ) ) | ( v[3] >> shift ); 
    } 
    memcpy(v,t,4); 
}

/************************************************************************
@Function: detect_and_fix_carriage_return
@Description:
	Detects carriage returns in text files and converts them to new-lines.
@Algorithm:
	The first line is read and if any carriage returns are found, the whole
	file is scanned and carriage returns converted to new-lines.
@Arguments:
	char* filename		file name.
@Returns:
	int					0, <0 if error.
*************************************************************************/
int			detect_and_fix_carriage_return(char* filename)
{
	FILE*		fd;
	char		aline[1024], c, n='\n';
	long		loc;

	if ( ( fd = fopen(filename, "r+") ) == NULL ) {
		if ( ( fd = fopen(filename, "r") ) == NULL ) {
			error_show(filename, __FILE__, __LINE__);
		} else {
			fgets(aline, 1024, fd);
			if ( strchr(aline, '\r') ) {
				fprintf(stderr, "Warning: Carriage returns in %s! Fix with bfix.\n", filename);
				exit(-1);
			}
			fclose(fd);
		}
		return(-1);
	}
	
	fgets(aline, 1024, fd);
	
	if ( !strchr(aline, '\r') ) {
		fclose(fd);
		return(0);
	}
	
	fprintf(stderr, "Warning: Fixing carriage returns in %s!\n", filename);
	
	fseek(fd, 0L, SEEK_SET);
	for ( loc=0; ( c = fgetc(fd) ) && !EOF; loc++ ) {
		if ( c == '\r' ) {
//			fseek(fd, loc, SEEK_SET);
			fseek(fd, -1L, SEEK_CUR);
			fputc(n, fd);
//			fseek(fd, loc, SEEK_SET);
			fseek(fd, 0L, SEEK_CUR);
		}
	}
	
	fclose(fd);
	
	return(0);
}

/************************************************************************
@Function: error_show
@Description:
	Displays the error with file and line reference.
@Algorithm:
	The function uses perror() to display a message containing the source
	file and line number where it originated.
@Arguments:
	char* message		a string to be included.
	char* file			the file name (should be __FILE__).
	int line			the line number (should be __LINE__).
@Returns:
	int					error number.
*************************************************************************/
int			error_show(char* message, char* file, int line)
{
	char		string[1024];
	
	sprintf(string, "[%s:%d] %s", file, line, message);
	
	perror(string);
	
	if ( errno ) fprintf(stderr, "errno %d: ", errno);
	
	switch ( errno ) {
		case ENOENT:
			fprintf(stderr, "Make sure the file exists and the path to the file is correct\n");
			break;
		case EACCES:
			fprintf(stderr, "Make sure you have permission to write to this directory\n");
			break;
		case ENOMEM:
			fprintf(stderr, "The memory available on this computer may be insufficient\n");
			break;
		default:
			fprintf(stderr, "\n");
	}
	
	return(errno);
}

