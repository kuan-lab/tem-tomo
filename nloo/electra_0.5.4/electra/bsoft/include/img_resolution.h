/*
	img_resolution.h
	Library routines to estimate resolution
	Author: Bernard Heymann
	Created: 20000611	Modified: 20050121
*/

// Function prototypes
float* 		img_resolution(Bimage* p, Bimage* pmod, Bimage* pmask, float sampling_ratio, char* psfile);
int 		img_fourier_shell_stats(Bimage* p);
//int 		img_spectral_signal_to_power(Bimage* p);
int 		myimg_spectral_signal_to_power(Bimage* p);
int 		myimg_spectral_signal_to_power_fact(Bimage* p, float* wfrac);
