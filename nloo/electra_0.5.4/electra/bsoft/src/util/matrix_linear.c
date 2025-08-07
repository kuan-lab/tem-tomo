/* 
	matrix_linear.c
	Solving sets of linear equations through matrix algebra 
	Author: Bernard Heymann 
	Created: 20000501	Modified: 20040311
*/ 
 
#include "matrix.h" 
#include "matrix_linear.h" 
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

/************************************************************************
@Function: linear_least_squares
@Description:
	Does a linear least squares fit between two vectors.
@Algorithm:
	The two input vectors must have elements between indices n1 and n2.
@Arguments:
	int n1			the starting index in each vector (default 0)
	int n2			the final index in each vector
	double* x		x vector (at least n2+1 elements)
	double* y		y vector (at least n2+1 elements)
	double* a		the intercept
	double* b		the slope
@Returns:
	double			the correlation index.
**************************************************************************/
double		linear_least_squares(int n1, int n2, double *x, double *y, double *a, double *b)
{
    int     	i, n;
    double  	sx, sx2, sy, sxy, sd, dy, denom;
	
	// Initial values in case the function returns prematurely
	*a = 0;
	*b = 1;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG linear_least_squares: Starting fit\n");
    
	if ( n2 < n1 ) {
		n = n1;
		n1 = n2;
		n2 = n;
	}
	if ( n1 < 0 ) n1 = 0;
	n = n2 - n1 + 1;
	
    sx = 0;
    sx2 = 0;
    sy = 0;
    sxy = 0;
    for ( i=n1; i<=n2; i++ ) {
    	sx += x[i];
    	sx2 += x[i]*x[i];
    	sy += y[i];
    	sxy += x[i]*y[i];
    }
	
	denom = n*sx2 - sx*sx;
	if ( fabs(denom) < 1e-30 ) return(0);
	
    *a = (sx2*sy-sx*sxy)/denom;
    *b = (n*sxy-sx*sy)/denom;
	
    sy = sy/n;
    sd = 0;
    dy = 0;
    for ( i=n1; i<=n2; i++ ) {
    	sd += (y[i]-sy)*(y[i]-sy);
    	dy += (*a + (*b*x[i]) - y[i])*(*a + (*b*x[i]) - y[i]);
    } 
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG linear_least_squares: a = %g  b = %g\n", *a, *b);
    
    return(1-dy/sd);
}

/************************************************************************
@Function: fit_polynomial
@Description:
	Fits a data set to a polynomial function.
@Algorithm:
	A polynomial of any order is fitted to the data using a least squares.
	The polynomial is defined as:
		f(x) = a0 + a1*x + a2*x^2 + ...
	The number of coefficients returned is the order plus one.
	The deviation is defined as:
		R = sqrt(sum(y - f(x))^2/n)
@Arguments:
	int n			number of data points.
	double* x		x array (at least order+1 values).
	double* y		y array (at least order+1 values).
	int order		polynomial order.
	double* coeff	array in which coefficients are returned (order+1 values)
	                (if NULL, no coefficients returned).
@Returns:
	double			the deviation.
**************************************************************************/
double		fit_polynomial(int n, double* x, double* y, int order, double* coeff)
{
	int				i, j, k;
	int				nt = order + 1;
	double*			a = (double *) balloc(nt*nt*sizeof(double));
	double*			b = (double *) balloc(nt*sizeof(double));
	double*			v = (double *) balloc(nt*sizeof(double));
	double			f, df, R=0;
	
	v[0] = 1;
	for ( i=0; i<n; i++ ) {
		for ( j=1; j<nt; j++ ) v[j] = v[j-1]*x[i];
		for ( j=0; j<nt; j++ ) b[j] += v[j]*y[i];
		for ( j=0; j<nt; j++ )
			for ( k=0; k<=j; k++ ) a[nt*j+k] += v[j]*v[k];
	}
	for ( j=0; j<nt-1; j++ )
		for ( k=j+1; k<nt; k++ ) a[nt*j+k] = a[nt*k+j];

	matrix_LU_decomposition(nt, a, b);
	
	for ( i=0, R=0; i<n; i++ ) {
		for ( j=1; j<nt; j++ ) v[j] = v[j-1]*x[i];
		for ( j=0, f=0; j<nt; j++ ) f += v[j]*b[j];
		df = f - y[i];
		R += df*df;
		if ( verbose & VERB_DEBUG )
			printf("%g\t%g\t%g\t%g\n", x[i], y[i], f, df);
	}
	
	R = sqrt(R/n);
	
	if ( coeff ) for ( i=0; i<nt; i++ ) coeff[i] = b[i];
	
	if ( verbose & VERB_FULL ) {
		printf("Polynomial: f(x) = %g", b[0]);
		for ( i=1; i<nt; i++ ) printf(" + %g x^%d", b[i], i);
		printf("\nR = %g\n\n", R);
	}

	bfree(a, nt*nt*sizeof(double));
	bfree(b, nt*sizeof(double));
	bfree(v, nt*sizeof(double));
	
	return(R);
}


/************************************************************************
@Function: matrix_LU_decomposition
@Description:
	Matrix inversion and optionally solves a set of linear equations.
@Algorithm:
	This solves the equation, A*x = b, and inverts matrix A by LU decomposition.
	The matrix A must be square and is converted to and replaced by its inverse.
	If the vector b is allocated and initialize, it is replaced by the solution, x.
@Arguments:
	int n			size of square matrix.
	double* a		square matrix A of size n x n.
	double* b		vector of size n, may be the NULL pointer.
@Returns:
	double	 		determinant.
**************************************************************************/
double 		matrix_LU_decomposition(int n, double* a, double* b)
{
	int 		i, j, k, ibig;
	double		big, det = 1;
	int*		ind = (int *) balloc(n*sizeof(int));
	double*		x = (double *) balloc(n*sizeof(double));
	double*		ai = (double *) balloc(n*n*sizeof(double));
	
	// Find the largest element in each row for scaling and pivoting
	// Additionally test for singularity
	for ( i=0; i<n; i++ ) {
		ind[i] = i; 		// Initialize the indices for tracking pivoting
		big = 0;
		for ( j=0; j<n; j++ )
			if ( big < fabs(a[i*n+j]) ) big = fabs(a[i*n+j]);
		if ( big < 1e-37 )
			printf("Error: Singular matrix in matrix_LU_decomposition\n");
		else
			x[i] = 1.0/big;
	}
	
	// The LU decomposition
	for ( j=0; j<n; j++ ) {
		for ( i=0; i<j; i++ ) {			// Calculating the upper triangle
			for ( k=0; k<i; k++ )		//  not the diagonal
				a[i*n+j] -= a[i*n+k]*a[k*n+j];
		}
		big = 0;			// Largest element for pivoting
		ibig = -1;
		for ( i=j; i<n; i++ ) { 		// Calculating the lower triangle
			for ( k=0; k<j; k++ )		//  and diagonal (i==j)
				a[i*n+j] -= a[i*n+k]*a[k*n+j];
			if ( x[i]*fabs(a[i*n+j]) >= big ) { // Get the biggest element
				big = x[i]*fabs(a[i*n+j]);		//   for pivoting
				ibig = i;
			}
		}
		if ( j != ibig ) {						// Interchange rows if necessary
			for ( k=0; k<n; k++ ) swap_doubles(&a[ibig*n+k], &a[j*n+k]);
			swap_doubles(&x[ibig], &x[j]); 		// Switch scales
			swap_integers(&ind[ibig], &ind[j]);	// Switch indices
			det = -det; 						// Change the sign of the determinant
		}
		if ( a[j*n+j] == 0 ) a[j*n+j] = 1e-37;
		for ( i=j+1; i<n; i++ )	a[i*n+j] /= a[j*n+j];	// Divide by pivot element
	}
	
	// Calculate the determinant
	for ( i=0; i<n; i++ ) det *= a[i*n+i];
	
	// Invert matrix
	for ( k=0; k<n; k++ ) {
		memset(x, 0, n*sizeof(double));
		x[k] = 1;
		for ( i=0; i<n; i++ )			// Forward substitution
			for ( j=0; j<i; j++ ) x[i] -= a[i*n+j]*x[j];
		for ( i=n-1; i>=0; i-- ) {		// Backward substitution
			for ( j=i+1; j<n; j++ ) x[i] -= a[i*n+j]*x[j];
			x[i] /= a[i*n+i];
			ai[i*n+k] = x[i];
		}
	}
	
	// Reorder the matrix back to the original row order
	for ( i=0; i<n; i++ ) {
		j = ind[i]; 					// Find the old row at this position
		for ( k=0; k<n; k++ ) 			// Pack it back in order into matrix A
			a[k*n+j] = ai[k*n+i];
	}
	
	// Multiply with the input vector to get the solution
//	if ( b ) {
//		for ( i=0; i<n; i++ ) {
//			x[i] = b[i];
//			b[i] = 0;
//		}
//		for ( i=0; i<n; i++ )
//			for ( j=0; j<n; j++ ) b[i] += a[i*n+j]*x[j];
//	}
	if ( b ) matrix_vector_multiply_in_place(n, a, b);
	
	bfree(ind, n*sizeof(int));
	bfree(x, n*sizeof(double));
	bfree(ai, n*n*sizeof(double));
	
	return(det);
}

