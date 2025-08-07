/*
	ext_proc.h
	Header file for execution of external programs
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004
*/

#ifndef _EXT_PROC_H__
#define _EXT_PROC_H__

#include <unistd.h>
#include <sys/wait.h>
#include "bsoft.h"

/*----------Function Prototypes-------------*/
int ext_execute(char* prog, char* parg);
int ext_rm(char* fname);
int ext_file_exist(char * fname);


#endif  // #ifndef _EXT_PROC_H__
