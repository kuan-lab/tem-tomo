/*
	linked_list.h
	Header file for generalized linked list functions
	Author: Bernard Heymann
	Created: 20031203 	Modified: 20041030
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//Function prototypes
char*		add_item(char** list, long size);
char*		remove_item(char** list, char* item, long size);
long 		kill_list(char* list, long size);
long 		count_list(char* list);

