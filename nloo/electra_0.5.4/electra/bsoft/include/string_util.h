/* 
	string_util.h 
	Header file for functions dealing with string modifications 
	Author: Bernard Heymann 
	Created: 19990722 	Modified: 20041130
*/ 

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <ctype.h> 
#include <unistd.h> 
#include <time.h>
#include <math.h> 

#include "utilities.h"

#ifndef _Bstring_
#define _Bstring_
struct Bstring {
	Bstring* next;
	char* str;
} ;
#endif

// Function prototypes 
char*   	extension(char* filename); 
int		   	filename_clean_copy(char* old_filename, char* new_filename, int length);
char*   	filename_base(char* filename);
char*   	filename_directory(char* filename);
char*   	filename_change_type(char* filename, char* type);
char*   	filename_change_path(char* filename, char* path);
char* 		number_filename(char* filename, int n, int digits);
char* 		insert_in_filename(char* filename, char* insert, char delim);
char*		parameter_file_path(char* filename);
int 		find_char_index(char* string, char c);
int 		char_replace(char* string, char c1, char c2);
char*		copystring(char* oldstring);
char* 		catenate2strings(char* string1, char* string2);
Bstring*	string_add(Bstring** list, char* string);
int			string_kill(Bstring* list);
char*		string_to_list(Bstring* list, int* size);
char* 		add_string_with_comma(char* list, char* string);
char** 		split_string(char* string, char* delimiter, int* number);
char** 		split_string_whitespace(char* string, int* number);
char** 		split_stringlist(char* string, int number);
char* 		string_array_to_list(char** array, int number);
float* 		split_string_into_floats(char* string, char* delimiter, int* number);
int			countTokens(char* aString);
