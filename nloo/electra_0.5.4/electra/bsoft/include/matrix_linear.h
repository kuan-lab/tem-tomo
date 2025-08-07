/* 
	matrix_linear.h 
	Header file for matrix functions 
	Author: Bernard Heymann 
	Created: 20000501 	Modified: 20030419
*/ 
 
// Function prototypes 
double		linear_least_squares(int n1, int n2, double *x, double *y, double *a, double *b);
double		fit_polynomial(int n, double* x, double* y, int order, double* coeff);
double 		matrix_LU_decomposition(int n, double* a, double* b);
