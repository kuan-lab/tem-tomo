/* 
	string_util.c 
	Library functions dealing with string modifications 
	Author: Bernard Heymann 
	Created: 19990722	Modified: 20041130
*/ 
 
#include "string_util.h" 
 
// Definition of the global variables 
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated

/************************************************************************
@Function: extension
@Description:
	Finds the extension of a file name to identify data and parameter file types.
@Algorithm:
	The extension is in general defined as the string following the last period.
	A colon indicates that the user wants the file to be interpreted as a
	particular type - thus the string after the colon is taken as the extension.
	In a comma-delimited list only the first file name is evaluated.
	The colon or comma and string after it is stripped from the file name.
	The extension is returned as lower case.
	The string returned has a length containing the extension plus one byte.
@Arguments:
	char* filename 	file name string.
@Returns:
	char*			the lower case extension, NULL if not found.
**************************************************************************/
char*   	extension(char* filename) 
{ 
	unsigned int	i;
	int				j; 
	char			ext[256];
	
	for ( i=0, j=-1; i<1024 && filename[i] > 0 && filename[i] != ','; i++ ) {
		if ( i > 1022 ) {
			error_show("File name %s is too long!", __FILE__, __LINE__);
			exit(-1);
		}
		if ( j > -1 ) {
			ext[j] = filename[i];
			j++;
		}
		if ( filename[i] == '.' || filename[i] == ':' ) j = 0;
	}
	
	ext[j] = 0;
 
	if ( strlen(ext) < 1 ) 
		return(NULL); 		// There is no extension!
 
	char*		the_ext = (char *) balloc((strlen(ext) + 1)*sizeof(char));
	
	for ( i=0; i<strlen(ext); i++ ) 
		the_ext[i] = tolower(ext[i]); 			// All in lower case
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG extension: %s (%ld)\n", the_ext, strlen(the_ext));
		
	return(the_ext); 
} 

/************************************************************************
@Function: filename_clean_copy
@Description:
	Copies a file name to a new string, except any part after ":".
@Algorithm:
	The colon and string after it is not copied to the new file name.
@Arguments:
	char* old_filename	original file name.
	char* new_filename	string to copy file name to.
	int length			maximum length.
@Returns:
	int					<0 for an error.
**************************************************************************/
int		   	filename_clean_copy(char* old_filename, char* new_filename, int length) 
{ 
	int 		i; 
	
	for ( i=0; i<length && old_filename[i] > 0 && old_filename[i] != ':' && old_filename[i] != '#'; i++ ) {
		if ( i > length - 2 ) {
			fprintf(stderr, "Error: File name %s is too long!\n", old_filename);
			return(-1);
		}
		new_filename[i] = old_filename[i];
	}
	
	new_filename[i] = 0;
 		
	return(0); 
} 

/************************************************************************
@Function: filename_base
@Description:
	Finds the base of a file name.
@Algorithm:
	Any leading path and extension is removed and a new string returned.
@Arguments:
	char* filename	the file name string.
@Returns:
	char*			the base.
**************************************************************************/
char*   	filename_base(char* filename) 
{
	char*		bstart = NULL;
	char*		bend = NULL;
	if ( ( bstart = strrchr(filename, '/') ) ) bstart += 1;
	else bstart = filename;
	if ( ( bend = strrchr(filename, '.') ) == NULL ) bend = bstart + strlen(filename);

	int 		size = (int) (bend - bstart);
	char*		base = (char *) balloc((size+1)*sizeof(char));
	strncpy(base, bstart, size);
	base[size] = 0;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG filename_base: %s (%ld)\n", base, strlen(base));
		
	return(base); 
} 

/************************************************************************
@Function: filename_directory
@Description:
	Finds the directory of a file name.
@Algorithm:
	The string before '/' is returned if it exists, otherwise './' is returned.
@Arguments:
	char* filename	the file name and path string.
@Returns:
	char*			the directory.
**************************************************************************/
char*   	filename_directory(char* filename) 
{
	char*		slash = strrchr(filename, '/');
	char*		dir = NULL;
	long		n;
	
	if ( slash ) {
		n = slash - filename + 1;
		dir = (char *) balloc((n+1)*sizeof(char));
		strncpy(dir, filename, n);
	} else {
		dir = (char *) balloc(3*sizeof(char));
		strcpy(dir, "./");
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG filename_directory: %s (%ld)\n", dir, strlen(dir));
		
	return(dir); 
} 

/************************************************************************
@Function: filename_change_type
@Description:
	Changes the type of a file.
@Algorithm:
	Any leading path and extension is removed and a new string returned.
@Arguments:
	char* filename	the file name string.
	char* type		the new type.
@Returns:
	char*			the new file name string.
**************************************************************************/
char*   	filename_change_type(char* filename, char* type) 
{
	char*		bstart = NULL;
	char*		bend = NULL;
	if ( ( bstart = strrchr(filename, '/') ) ) bstart += 1;
	else bstart = filename;
	if ( ( bend = strrchr(filename, '.') ) == NULL ) bend = filename + strlen(filename);
	
	int 		typesize = (int) strlen(type);
	if ( typesize > 20 ) typesize = 20;

	int 		basesize = (int) (bend - bstart);
	int 		size = basesize + 1 + typesize;
	
	char*		newname = (char *) balloc((size+1)*sizeof(char));
	strncpy(newname, bstart, basesize);
	newname[basesize] = '.';
	strncpy(&newname[basesize+1], type, typesize);
	newname[size] = 0;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG filename_change_type: %s (%ld)\n", newname, strlen(newname));
		
	return(newname); 
} 

/************************************************************************
@Function: filename_change_path
@Description:
	Changes the path of a file.
@Algorithm:
	Any leading path is removed and the new path added.
@Arguments:
	char* filename	the file name string.
	char* path		the new path.
@Returns:
	char*			the new file name string.
**************************************************************************/
char*   	filename_change_path(char* filename, char* path) 
{
	char*		fstart = NULL;
	if ( ( fstart = strrchr(filename, '/') ) ) fstart += 1;
	else fstart = filename;
	
	int 		size = (int) (strlen(path) + 1 + strlen(fstart));
	
	char*		newname = (char *) balloc((size+1)*sizeof(char));
	sprintf(newname, "%s/%s", path, fstart);
	newname[size] = 0;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG filename_change_path: %s (%ld)\n", newname, strlen(newname));
		
	return(newname); 
} 

/************************************************************************
@Function: number_filename
@Description:
	Modifies a file name with a number before the last period.
@Algorithm:
	The input filename is modified and copied to a new string with
	a length of the old string plus the number of digits plus one for
	the underscore. The old string is not affected.
@Arguments:
	char* filename		the base file name.
	int n				the number to append (< 9999).
	int digits			the length of the number in digits.
@Returns:
	char* 				the numbered file name.
**************************************************************************/
char* 		number_filename(char* filename, int n, int digits)
{
	char*		extension = strrchr(filename, '.');
	char*		newname = (char *) balloc(strlen(filename)+digits+2);
	char*		number = (char *) balloc(10);
	int 		baselength = (int) (extension - filename);
	char		format[20];
	
	sprintf(format, "_%%0%dd", digits);

	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG number_filename: filename = %s   extension = %s   baselength = %d\n", 
			filename, extension, baselength);
		printf("DEBUG number_filename: digits = %d   format = %s   total length = %ld\n", 
			digits, format, strlen(filename)+digits+2);
	}
		
	strcpy(newname, filename);
	newname[baselength] = 0;
	sprintf(number, format, n);
	strcat(newname, number);
	strcat(newname, extension);
	
	bfree(number, 10);
	
	return(newname);
}

/************************************************************************
@Function: insert_in_filename
@Description:
	Modifies a file name with a given string before the last period.
@Algorithm:
	The input filename is modified and copied to a new string with
	a length of the old string plus the number of digits of the insert.
	An optional delimiting character can be added to the beginning 
	of the insert.
	The old string is not affected.
@Arguments:
	char* filename		the base file name.
	char* insert		string to be inserted.
	char delim			delimiting character before insert, 0 if none.
@Returns:
	char* 				the numbered file name.
**************************************************************************/
char* 		insert_in_filename(char* filename, char* insert, char delim)
{
	if ( verbose & VERB_DEBUG )
		printf("DEBUG insert_in_filename: filename = %s   insert = %s\n", 
			filename, insert);
	
	long		size = strlen(filename) + strlen(insert) + 1;
	if ( delim > 0 ) size++;
	
	char*		ext = strrchr(filename, '.');
	char*		newname = (char *) balloc(size);
	int 		baselength = (int) (ext - filename);
	
	strcpy(newname, filename);
	if ( delim > 0 ) newname[baselength++] = delim;
	newname[baselength] = 0;
	strcat(newname, insert);
	strcat(newname, ext);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG insert_in_filename: newname = %s\n", newname);
	
	return(newname);
}

/************************************************************************
@Function: parameter_file_path
@Description:
	Finds the parameter file path.
@Algorithm:
	The parameter file path should primarily be defined in the environmental
	variable "BPARAM".
	Otherwise a default is used which may or may not be valid.
	The string returned has a length containing the path plus one byte.
@Arguments:
	char* filename	the parameter filename.
@Returns:
	char*			the full path and file name.
**************************************************************************/
char*		parameter_file_path(char* filename)
{
	char	path[1024] = "";
	
	if ( getenv("BPARAM") ) {
		strcpy(path, getenv("BPARAM"));
 		if ( path[strlen(path)-1] != ']'  &&
 				path[strlen(path)-1] != '/')
 			strcat(path, "/");
	} else
		strcpy(path, "/usr/local/bsoft/parameters/");
	strcat(path, filename);
	
	char*	file_path = (char *) balloc((strlen(path) + 1)*sizeof(char));
	strcpy(file_path, path);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG parameter_file_path: %s\n", file_path);
	
	return(file_path);
}
 
/************************************************************************
@Function: find_char_index
@Description:
	Finds the index of the first occurrence of a character in a string.
@Algorithm:
	Searches the string until it finds the desired character.
@Arguments:
	char* string	the string to be searched.
	char c			the character to be found.
@Returns:
	int 			the index of the first occurrence of the character.
					-1 if not found.
**************************************************************************/
int 		find_char_index(char* string, char c)
{
	int			i = 0;
	
	for ( i = 0; string[i] != c && i < (int)strlen(string) ; i++ ) ;
	
	if ( string[i] != c ) i = -1;
	
	return(i);
}

/************************************************************************
@Function: char_replace
@Description:
	Replaces every occurrence of one character with another in a string.
@Algorithm:
	Searches the string until it finds the desired character and
	replaces it with the new one.
@Arguments:
	char* string	the string to be modified.
	char c1			the character to be replaced.
	char c2			the new character.
@Returns:
	int 			the number of characters replaced.
**************************************************************************/
int 		char_replace(char* string, char c1, char c2)
{
	int			i = 0;
	
	for ( ; string; string++ ) if ( *string == c1 ) {
		*string = c2;
		i++;
	}
	
	return(i);
}

/************************************************************************
@Function: copystring
@Description:
	Copy a string to a new character array with memory allocation.
@Algorithm:
	Allocates memory based on the size of the given string.
	Copies the given string into the new string.
@Arguments:
	char* oldstring	the string to be copied.
@Returns:
	char*			the copied string.
**************************************************************************/
char*		copystring(char* oldstring)
{
	if ( !oldstring ) return(NULL);
	
	if ( strlen(oldstring) < 1 ) return(NULL);
	
	char*		string = (char *) balloc((strlen(oldstring)+1)*sizeof(char));
	
	memcpy(string, oldstring, strlen(oldstring));
	
	return(string);
}

/************************************************************************
@Function: catenate2strings
@Description:
	Catenates two strings with memory allocation.
@Algorithm:
	Allocates memory based on the sizes of the two strings.
	Copies the first string into the new string.
	Catenates the second string to the new string.
	Frees the memory for the old first string.
@Arguments:
	char* string1		first string - pointer reassigned to new string.
	char* string2		second string - unchanged.
@Returns:
	char*				concatenated string.
**************************************************************************/
char* 		catenate2strings(char* string1, char* string2)
{
	char*		oldstring = string1;
	char*		string = NULL;
	
	if ( string1 && string2 ) {
		string = (char *) balloc((strlen(string1)+strlen(string2)+1)*sizeof(char));
		if ( strlen(string1) > 0 )
			memcpy(string, string1, strlen(string1));
		if ( strlen(string2) > 0 )
			memcpy(string + strlen(string1), string2, strlen(string2));
		bfree(oldstring, (strlen(oldstring)+1)*sizeof(char));
	} else if ( string1 ) {
		string = string1;
	} else if ( string2 ) {
		string = (char *) balloc((strlen(string2)+1)*sizeof(char));
		if ( strlen(string2) > 0 )
			memcpy(string, string2, strlen(string2));
	}
		
	string1 = string;
	
	return(string);
}

/************************************************************************
@Function: string_add
@Description:
	Adds a string to a linked list.
@Algorithm:
	.
@Arguments:
	Bstring** list 		the string linked list.
	char* string		string - unchanged.
@Returns:
	Bstring*				new string structure.
**************************************************************************/
Bstring*	string_add(Bstring** list, char* string)
{
	Bstring* 		this_string = *list;
	
	if ( !string ) return(this_string);
	if ( strlen(string) < 1 ) return(this_string);
	
	Bstring* 		new_string = (Bstring *) balloc(sizeof(Bstring));
	
	new_string->str = (char *) balloc(strlen(string)+1);
	
	strcpy(new_string->str, string);
	
	if ( !this_string )
		*list = new_string;
	else {
		while ( this_string->next ) this_string = this_string->next;
		this_string->next = new_string;
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG string_add: %s\n", new_string->str);
	
	return(new_string);
}

/************************************************************************
@Function: string_kill
@Description:
	Kills a string linked list.
@Algorithm:
	.
@Arguments:
	Bstring* list 		the string linked list.
@Returns:
	int					0.
**************************************************************************/
int			string_kill(Bstring* list)
{
	if ( !list ) return(0);
	
	Bstring		*s, *s2;
	
	for ( s = list; s; ) {
		s2 = s->next;
		bfree_string(s->str);
		bfree(s, sizeof(Bstring));
		s = s2;
	}
	
	return(0);
}

/************************************************************************
@Function: string_to_list
@Description:
	Converts a linked list of strings to a null-delimited list.
@Algorithm:
	.
@Arguments:
	Bstring* list 		the string linked list.
	int* size			size of the new null-delimited list.
@Returns:
	char*				new null-delimited list of strings.
**************************************************************************/
char*		string_to_list(Bstring* list, int* size)
{
	unsigned long   new_size = 0;
	Bstring*		s;
	
	for ( s = list; s; s = s->next ) new_size += strlen(s->str) + 1;
	
	char*			array = (char *) balloc(new_size);
	char*			aptr = array;
	
	for ( s = list; s; s = s->next ) {
		strcpy(aptr, s->str);
		aptr += strlen(s->str) + 1;
	}
	
	*size = (int) new_size;
	
	return(array);
}

/************************************************************************
@Function: add_string_with_comma
@Description:
	Adds a string to a comma-delimited list in a string.
@Algorithm:
	The array consists of a number of strings, delimited by commas.
	Allocates memory based on the sizes of the array and string.
	Copies the array into the new array.
	Copies the string to the end of the new array.
	Free the memory for the old array.
@Arguments:
	char* list	 		the list as a simple string - deallocated.
	char* string		string to add - unchanged.
@Returns:
	char*				concatenated string.
**************************************************************************/
char* 		add_string_with_comma(char* list, char* string)
{
	if ( !string ) return(list);
	
	int 		i, new_size = 0;
	char*		new_list = NULL;

	if ( !list ) {
		new_size = (int) (strlen(string) + 1);
		new_list = (char *) balloc(new_size);
		memcpy(new_list, string, new_size);
	} else {
		new_size = (int) (strlen(list) + strlen(string) + 2);
		new_list = (char *) balloc(new_size);
		i = (int) strlen(list);
		memcpy(new_list, list, i);
		new_list[i] = ',';
		memcpy(new_list+i+1, string, strlen(string));
	}
	
	bfree_string(list);
	
	return(new_list);
}

/******************************************************************************
@Function: split_string
@Description:	
	Splits a string into an array of strings with any given delimiter.
@Algorithm:
	The string is parsed to find all the delimiter strings to determine the
	number of string pointers to allocate. The string is parsed a second 
	time to allocate memory and copy each substring into the new string array.
	The number of substrings found is returned as well.
@Arguments:
	char* string		the line to be converted
	char* delimiter 	a delimiter or separator string
	int* number 		the number of substrings returned
@Returns:
	char** 				an array of substrings
*******************************************************************************/
char** 		split_string(char* string, char* delimiter, int* number)
{
	int 		i, n;
	int		l = (int) strlen(delimiter);
	char		*iptr, *aptr;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG split_string: splitting %s at %s\n", string, delimiter);
	
	for ( n=1, iptr = string; iptr[0]; iptr++ )
		if ( strncmp(iptr, delimiter, l) == 0 && iptr[l] > 0 ) n++;
	
    char**		stringarray = (char **) balloc(n*sizeof(char*));
	
	for ( i=0, aptr = iptr = string; i < n; iptr++ ) {
		if ( strncmp(iptr, delimiter, l) == 0 || iptr[0] == 0 ) {
			stringarray[i] = (char *) balloc((iptr-aptr+1)*sizeof(char));
			strncpy(stringarray[i], aptr, iptr-aptr);
			if ( verbose & VERB_DEBUG )
				printf("DEBUG split_string: %s\n", stringarray[i]);
			aptr = iptr + l;
			i++;
		}
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG split_string: string %s split at %s\n", string, delimiter);
	
	*number = n;
	
	return(stringarray);
}

/******************************************************************************
@Function: split_string_whitespace
@Description:	
	Splits a string into an array of strings at whitespace.
@Algorithm:
	The string is parsed to find all the whitespace to determine the
	number of string pointers to allocate. The string is parsed a second 
	time to allocate memory and copy each substring into the new string array.
	The number of substrings found is returned as well.
@Arguments:
	char* string		the line to be converted
	int* number 		the number of substrings returned
@Returns:
	char** 				an array of substrings
*******************************************************************************/
char** 		split_string_whitespace(char* string, int* number)
{
	int 		n, s;
	char		*iptr, *aptr;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG split_string_whitespace: %s\n", string);
	
	for ( n=s=0, iptr = string; iptr[0]; iptr++ ) {
		if ( s ) {
			if ( isspace(iptr[0]) ) {
				n++;
				s = 0;
			}
		} else if ( !isspace(iptr[0]) ) {
			s = 1;
		}
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG split_string_whitespace: split into %d\n", n);
	
    char**		stringarray = (char **) balloc(n*sizeof(char*));
	
	for ( n=s=0, aptr = iptr = string; iptr[0]; iptr++ ) {
		if ( s ) {
			if ( isspace(iptr[0]) ) {
				stringarray[n] = (char *) balloc((iptr-aptr+1)*sizeof(char));
				strncpy(stringarray[n], aptr, iptr-aptr);
				if ( verbose & VERB_DEBUG )
					printf("DEBUG split_string_whitespace: %s\n", stringarray[n]);
				n++;
				s = 0;
			}
		} else if ( !isspace(iptr[0]) ) {
			aptr = iptr;
			s = 1;
		}
	}
	
	*number = n;
	
	return(stringarray);
}

/******************************************************************************
@Function: split_stringlist
@Description:	
	Splits a list of null-delimited strings into an array of strings.
@Algorithm:
	.
@Arguments:
	char* string		list of null-delimited strings
	int number 			number of strings
@Returns:
	char** 				array of strings
*******************************************************************************/
char** 		split_stringlist(char* string, int number)
{
	int 		i, length;
    char**		stringarray = (char **) balloc(number*sizeof(char*));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG split_stringlist: splitting %s into %d strings\n", string, number);
	
	for ( i=0; i<number; i++ ) {
		length = (int) (strlen(string) + 1);
    	stringarray[i] = (char *) balloc(length*sizeof(char));
		strcpy(stringarray[i], string);
		string += length;
	}
	
	return(stringarray);
}

/******************************************************************************
@Function: string_array_to_list
@Description:	
	Constructs a null-delimited list of strings from an array of strings.
@Algorithm:
	.
@Arguments:
	char** array		array of strings
	int number 			number of strings
@Returns:
	char* 				list of null-delimited strings
*******************************************************************************/
char* 		string_array_to_list(char** array, int number)
{
	int 		i, length = 0;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG string_array_to_list: converting %d strings to a list\n", number);
	
	for ( i=0; i<number; i++ ) length += strlen(array[i]) + 1;
	
    char*		string = (char *) balloc(length*sizeof(char));
	
	char*		aptr = string;
	
	for ( i=0; i<number; i++ ) {
		strncpy(aptr, array[i], strlen(array[i]));
		aptr += strlen(array[i]) + 1;
	}
	
	return(string);
}

/******************************************************************************
@Function: split_string_into_floats
@Description:	
	Splits a string into an array of floating point numbers with any given delimiter.
@Algorithm:
	The string token function is used to find each substring.
	Each substring is read into a temporary character array.
	Memory is allocated for one substring at a time.
	The number of substrings found is counted and returned as well.
	The substrings are converted to floating point numbers.
@Arguments:
	char* string		the line to be converted.
	char* delimiter 	a delimiter or separator string (NULL => whitespace).
	int* number 		the number of substrings returned.
@Returns:
	float* 				an array of floating point numbers.
*******************************************************************************/
float* 		split_string_into_floats(char* string, char* delimiter, int* number)
{
    char**		stringarray = NULL;
	
	if ( !delimiter )
		stringarray = split_string_whitespace(string, number);
	else
		stringarray = split_string(string, delimiter, number);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG split_string_into_floats: splitting %s at %s\n", string, delimiter);
	
	int 		i, n = *number;
	
	float*		array = (float *) balloc(n*sizeof(float));
	
	for ( i=0; i<n; i++ ) {
		sscanf(stringarray[i], "%f", &array[i]);
		bfree(stringarray[i], (strlen(stringarray[i])+1)*sizeof(char));
	}
	
	bfree(stringarray, n*sizeof(char*));
	
	return(array);
}

/******************************************************************************
@Function: countTokens
Author: Dan Krainak
@Description:	
	Counts the number of tokens in a string line 
@Algorithm:
	The string token function is used to find each substring
	The delimiter is a space
@Arguments:
	char* aString	the line to be analyzed
@Returns:
	int 			the number of tokens in the line.
*******************************************************************************/
int			countTokens(char* aString)
{
	char *tokenPtr;
	int i = 0;

	tokenPtr = strtok(aString, " " );
	
	while (tokenPtr != NULL ) {
		i++;
		tokenPtr = strtok ( NULL, " " );
	}
	
	return(i);	/* return the number of tokens found in the string */
}

