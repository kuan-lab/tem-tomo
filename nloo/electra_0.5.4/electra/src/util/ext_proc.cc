/*
	ext_proc.cc
	Routines for executing external programs
	Author: Giovanni Cardone
	Created: 2004 	Modified: 2004
*/

#include "ext_proc.h"

/************************************************************************
@Function: ext_execute
@Description:
	Execute an external program.
@Algorithm:
	Create a child process with fork command, then waits its completion.
@Arguments:
	char * prog		program name
	char * parg		argument to pass to program
@Returns:
	int				0 if execution successfull
**************************************************************************/
int ext_execute(char* prog, char* parg) {

	int child_status;

	fflush( (FILE *) NULL);
	
	if ( verbose & VERB_DEBUG )
		printf("original: pid %d\n", getpid());

	switch ( fork() ) {
	case -1:
		fprintf(stderr,"ERROR: fork for external program %s %s failed!",prog,parg);
		return -1;
	case 0:    /* child */
		if ( verbose & VERB_DEBUG )
			printf("child: pid %d\n", getpid());
		execlp(prog, prog, parg, NULL);
		break;
	default:	/* parent */
		sleep(1);
		wait(&child_status);	/* wait */
		if (WIFEXITED(child_status)) {
			if ( verbose & VERB_DEBUG )
				printf("parent: %d, exit %d\n", getpid(), WEXITSTATUS(child_status));
		} else {
			fprintf(stderr,"ERROR: program %s %s terminated abnormally (%d - status %d)\n",
				prog,parg,WIFEXITED(child_status),WEXITSTATUS(child_status));
			return -1;
		}
		if (WEXITSTATUS(child_status) != 0) {
			if ( verbose & (VERB_PROCESS | VERB_DEBUG) )
				printf("WARNING: program %s %s terminated with status %d\n",prog,parg,WEXITSTATUS(child_status)); 
			return -1;
		}
		break;
	}

	return 0;
}

/************************************************************************
@Function: ext_file_exist
@Description:
	Check if given file already present in the directory.
@Algorithm:
@Arguments:
	chart*			filename
@Returns:
	int				>0 if file exists
*************************************************************************/
int ext_file_exist(char * fname)
{
	FILE * fd;
	if ( ( fd = fopen(fname, "r" ) ) == NULL) {
		return(0);
	} else {
		fclose(fd);
		return(1);
	}

}
