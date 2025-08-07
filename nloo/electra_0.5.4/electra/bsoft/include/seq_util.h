/* 
	seq_util.h 
	Header file for sequence utilities 
	Author: Bernard Heymann 
	Created: 20001029 	Modified: 20040112 
*/ 
 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <ctype.h> 
#include <unistd.h> 
#include <math.h> 
#include <time.h>

#include "rwmolecule.h"
  
// Function prototypes 
int			getcode3(char c, char* cod);
char		getcode1(char* acode); 
int 		complement_all_sequences(Bmolgroup* molgroup);
int 		complement_sequence(char* nucseq);
char		get_complement(char nuc);
int 		translate_all_sequences(Bmolgroup* molgroup, int frame, char* gcname);
char* 		sequence_translate(char* nucseq, char* gencode);

