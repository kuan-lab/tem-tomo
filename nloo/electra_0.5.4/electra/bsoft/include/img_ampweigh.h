/*
	img_ampweigh.h
	Header file for image utilities
	Author: Bernard Heymann
	Created: 20000430 	Modified: 20041210
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "rwimg.h"

// Function prototypes
int 		img_weigh_amplitudes_by_fom(Bimage* p);
int 		img_weigh_amplitudes_with_C_curve(Bimage* p) ;
int 		img_weigh_amplitudes(Bimage* p, Bimage* pref, Bimage* pmask);

