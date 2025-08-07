/* 
	timer.h 
	Utilities for timing functions 
	Author: Bernard Heymann 
	Created: 20010316	Modified: 20010316 
*/ 
 
#include <stdio.h> 
#include <time.h>

// Function prototypes 
time_t		timer_start();
time_t		timer_report(time_t lasttime);

