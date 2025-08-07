/*
	linked_list.c
	Generalized linked list functions
	Author: Bernard Heymann
	Created: 20031203 	Modified: 20031203
*/

#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated

/************************************************************************
@Function: add_item
@Description:
	Adds an item to a linked list.
@Algorithm:
	If the list is not defined, the new item becomes the first in the list.
	Otherwise, the list is traversed to find the end and the new item appended.
	Any structure with a pointer to itself as a first element can be used in a
	linked list. However, a linked list can only consist of one type of structure.
@Arguments:
	char** list			pointer to first item in the list.
	long size			size of item.
@Returns:
	char* 				new item.
*************************************************************************/
char*		add_item(char** list, long size)
{
	char**		curr_item = (char **) *list;
	char*		new_item = (char *) balloc(size);
	char*		next_item = NULL;
	
//	printf("%p %p %p\n", list, curr_item, new_item);
	
	if ( !curr_item )
		*list = new_item;
	else {
		while ( ( next_item = *curr_item ) ) curr_item = (char **) next_item;
		*curr_item = new_item;
	}
	
	return(new_item);
}

/************************************************************************
@Function: remove_item
@Description:
	Finds the given item and deletes it from the linked list.
@Algorithm:
	If the item is the first in the list, the list pointer is set to point
	to the next item.
	Otherwise, the list is traversed to find the item, the previous item's
	pointer is set to the next item, and the current item deallocated.
@Arguments:
	char** list			pointer to first item in the linked list.
	char* item			item to be deleted.
	long size			size of item.
@Returns:
	char* 				item after the one removed.
*************************************************************************/
char*		remove_item(char** list, char* item, long size)
{
	if ( !(*list) || !item ) return(NULL);
	
	char**		curr_item = (char **) *list;
	char**		next_item;
	
//	printf("%p %p\n", list, curr_item);
	
	if ( *list == item ) {
		*list = *curr_item;
	} else {
		for ( ; curr_item && *curr_item != item; curr_item = (char **) *curr_item ) ;
		next_item = (char **) item;
		*curr_item = *next_item;
	}
	
//	*item = NULL;
	bfree(item, size);
	
	item = *curr_item;
	
	return(item);
}

/************************************************************************
@Function: kill_list
@Description:
	Frees all the items in a linked list.
@Algorithm:
	The list is traversed, setting a pointer to the next item before
	deallocating the current item.
@Arguments:
	char* list			first item in the linked list.
	long size			size of item.
@Returns:
	long 				number of items deallocated.
*************************************************************************/
long 		kill_list(char* list, long size)
{
	if ( !list ) return(0);
	
	long		n = 0;
	char**		item;
	
	for ( item = (char **) list; item; n++ ) {
		list = *item;
		bfree(item, size);
		item = (char **) list;
	}
			
	return(n);
}

/************************************************************************
@Function: count_list
@Description:
	Counts the number of items in a linked list.
@Algorithm:
	.
@Arguments:
	char* list			first item in the linked list.
@Returns:
	long 				number of items in the list.
*************************************************************************/
long 		count_list(char* list)
{
	if ( !list ) return(0);
	
	long		n = 0;
	char**		item;
	
	for ( item = (char **) list; item; item = (char **) *item ) n++;
			
	return(n);
}


