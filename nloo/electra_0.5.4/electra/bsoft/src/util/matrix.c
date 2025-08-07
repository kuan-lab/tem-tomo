/* 
	matrix.c
	Matrix manipulation functions 
	Author: Bernard Heymann 
	Created: 20000501	Modified: 20050201
*/ 
 
#include "matrix.h" 
#include "matrix_linear.h"
#include "math_util.h"
#include "random_numbers.h"
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

/*
	General matrix functions:
	------------------------
*/

/************************************************************************
@Function: init_unit_matrix
@Description:
	Initializes a unit matrix with any number of dimensions.
@Algorithm:
	.
@Arguments:
	int dim 		number of dimensions.
	int size		size of each dimension.
@Returns:
	float* 			unit matrix.
**************************************************************************/
float*		init_unit_matrix(int dim, int size) 
{ 
	int 		i, j, ind, mat_size = 1;
	
	for ( i=0; i<dim; i++ ) mat_size *= size;
	
	float*		mat = (float *) balloc(mat_size*sizeof(float)); 
	 
	for ( i=0; i<size; i++ ) {
		ind = i;
		for ( j=1; j<dim; j++ )
			ind = ind*size + i;
		mat[ind] = 1;
	}
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG init_unit_matrix: dimension = %d,  size = %d:\n",
				dim, size);
		for ( i=0; i<size; i++ ) {
			if ( dim == 2 ) {
				for ( j=0; j<size; j++ )
					printf("\t%g", mat[i*size+j]);
				printf("\n");
			}
		}
	}
	
 
	return(mat); 
} 
 
/************************************************************************
@Function: matrix_determinant
@Description:
	Calculates the determinant of any size square matrix.
@Algorithm:
	.
@Arguments:
	int n			matrix size, n.
	float*	 		a nxn matrix.
@Returns:
	double 			determinant.
**************************************************************************/
double		matrix_determinant(int n, float* mat) 
{
	int			i;
	double*		mat_copy = (double *) balloc(n*n*sizeof(double));
	
	for ( i=0; i<n*n; i++ ) mat_copy[i] = mat[i];
	
	double		det = matrix_LU_decomposition(n, mat_copy, NULL);
	
	bfree(mat_copy, n*n*sizeof(double));
	
	return(det); 
} 
 
/************************************************************************
@Function: matrix_product
@Description:
	Multiplies two matrices.
@Algorithm:
	mat(i,j) = sum(m1(i,k)*m2(k,j)).
@Arguments:
	int m			rows in first matrix.
	int n			columns in first matrix = rows in second matrix.
	int p			columns in second matrix.
	float* m1		first matrix.
	float* m2		second matrix.
@Returns:
	float*			resultant matrix of dimension m x p.
**************************************************************************/
float*		matrix_product(int m, int n, int p, float* m1, float* m2)
{
	int 		i, j, k;
	float*		mat = (float *) balloc(m*p*sizeof(float));
	
	for ( i=0; i<m; i++ )
		for ( j=0; j<p; j++ )
			for ( k=0; k<n; k++ ) mat[i*p+j] += m1[i*n+k]*m2[k*p+j];
	
	return(mat);
}

/************************************************************************
@Function: matrix_transpose
@Description:
	Transpose a matrix in place.
@Algorithm:
	transpose(j,i) = mat(i,j).
@Arguments:
	int m			rows in matrix.
	int n			columns matrix.
	float* mat		matrix (replaced by transpose).
@Returns:
	float*			resultant matrix of dimension n x m.
**************************************************************************/
float*		matrix_transpose(int m, int n, float* mat)
{
	int 		i, j;
	
	for ( i=1; i<m; i++ )
		for ( j=0; j<i; j++ )
			swap_floats(&mat[j*m+i], &mat[i*n+j]);
	
	return(mat);
}

/************************************************************************
@Function: show_matrix
@Description:
	Print out a matrix.
@Algorithm:
	Matrix elements less than 1e-6 are set to zero for display.
@Arguments:
	int m			rows in matrix.
	int n			columns matrix.
	float* mat		matrix.
@Returns:
	int 			0.
**************************************************************************/
int 		show_matrix(int m, int n, float* mat)
{
	int 		i, j;
	float		val;
	
	for ( i=0; i<m; i++ ) {			// Rows
		for ( j=0; j<n; j++ ) {		// Columns
			val = mat[i*n+j];
			if ( fabs(val) < 1e-6 ) val = 0;
			printf("%10.6g ", val);
		}
		printf("\n");
	}
	printf("\n");
	
	return(0);
}

/************************************************************************
@Function: show_matrix
@Description:
	Print out a matrix.
@Algorithm:
	Matrix elements less than 1e-6 are set to zero for display.
@Arguments:
	int m			rows in matrix.
	int n			columns matrix.
	double* mat		matrix.
@Returns:
	int 			0.
**************************************************************************/
int 		show_matrix(int m, int n, double* mat)
{
	int 		i, j;
	double		val;
	
	for ( i=0; i<m; i++ ) {			// Rows
		for ( j=0; j<n; j++ ) {		// Columns
			val = mat[i*n+j];
			if ( fabs(val) < 1e-6 ) val = 0;
			printf("%10.6g ", val);
		}
		printf("\n");
	}
	printf("\n");
	
	return(0);
}

/************************************************************************
@Function: matrix_vector_multiply
@Description:
	Multiplies a square matrix by a vector.
@Algorithm:
	new_vector[i] = sum_j(mat[i][j]*vector[j]).
@Arguments:
	int n			size of vector and size of matrix side
	float* mat 		nxn matrix.
	float* vector 	n-value vector.
@Returns:
	float* 			resultant n-value vector.
**************************************************************************/
float*		matrix_vector_multiply(int n, float* mat, float* vector) 
{
	int			i, j;
	float*		new_vec = (float *) balloc(n*sizeof(float)); 
	 
	for ( i=0; i<n; i++ )
		for ( j=0; j<n; j++ )
			new_vec[i] += mat[n*i+j]*vector[j]; 
	 
	return(new_vec); 
} 

/************************************************************************
@Function: matrix_vector_multiply_in_place
@Description:
	Multiplies a square matrix by a vector.
@Algorithm:
	new_vector[i] = sum_j(mat[i][j]*vector[j]).
@Arguments:
	int n			size of vector and size of matrix side
	double* mat		nxn matrix.
	double* vec		n-value vector, replaced by result vector.
@Returns:
	int				0.
**************************************************************************/
int			matrix_vector_multiply_in_place(int n, double* mat, double* vec)
{
	int			i, j;
	
	double*		t = (double *) balloc(n*sizeof(double));
	
	for ( i=0; i<n; i++ )
		for ( j=0; j<n; j++ ) t[i] += mat[i*n+j]*vec[j];

	for ( i=0; i<n; i++ ) vec[i] = t[i];
	
	bfree(t, n*sizeof(double));
	
	return(0);
}

/************************************************************************
@Function: vector_normalize
@Description:
	Normalizes a vector of any size.
@Algorithm:
	The vector length is calculated and all elements of the vector divided by it.
	If the length is zero, the last element in the vector is set to 1.
@Arguments:
	int n			vector size (number of elements).
	float* vector	n-value vector.
@Returns:
	double 			vector length.
**************************************************************************/
double 		vector_normalize(int n, float* vector)
{
	int 		i;
	double		length = 0;
	
	for ( i=0; i<n; i++ )
		length += vector[i]*vector[i];
	
	if ( length > 0 ) {
		length = sqrt(length);
		for ( i=0; i<n; i++ )
			vector[i] /= length;
	} else {
		vector[n-1] = 1;
	}
	
	return(length);
}

double 		vector_normalize(int n, double* vector)
{
	int 		i;
	double		length = 0;
	
	for ( i=0; i<n; i++ )
		length += vector[i]*vector[i];
	
	if ( length > 0 ) {
		length = sqrt(length);
		for ( i=0; i<n; i++ )
			vector[i] /= length;
	} else {
		vector[n-1] = 1;
	}
	
	return(length);
}

/************************************************************************
@Function: vector_scalar_product
@Description:
	Calculates the sum of the inner product from two vectors.
@Algorithm:
	Scalar product:
		s = sum_i(v1[i]*v2[i])
@Arguments:
	int n			vector size.
	float* vector1	n-value first vector.
	float* vector2	n-value second vector.
@Returns:
	double 			scalar result.
**************************************************************************/
double 		vector_scalar_product(int n, float* vector1, float* vector2)
{
	int 		i;
	float		scalar = 0;
	
	for ( i=0; i<n; i++ )
		scalar += vector1[i]*vector2[i];
	
	return(scalar);
}

/************************************************************************
@Function: vector_vector_product
@Description:
	Calculates the vector product from two vectors.
@Algorithm:
	Vector product:
		x = y1*z2 - z1*y2
		y = z1*x2 - x1*z2
		z = x1*y2 - y1*x2
@Arguments:
	int n			vector size.
	float* vector1	n-value first vector.
	float* vector2	n-value second vector.
@Returns:
	float* 			n-value result vector.
**************************************************************************/
float* 		vector_vector_product(int n, float* vector1, float* vector2)
{
	float*		vector = (float *) balloc(n*sizeof(float));
	
	vector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
	vector[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
	vector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
	
	return(vector);
}

/************************************************************************
@Function: show_view
@Description:
	Displays a view.
@Algorithm:
	.
@Arguments:
	View v				the view.
@Returns:
	int					0.
**************************************************************************/
int			show_view(View v)
{
	printf("%7.4f\t%7.4f\t%7.4f\t%7.2f\n", v.x, v.y, v.z, v.a*180.0/PI);
	
	return(0);
}

/************************************************************************
@Function: show_view
@Description:
	Displays a linked list of views.
@Algorithm:
	.
@Arguments:
	View* v				the linked list of views.
@Returns:
	int					number of views.
**************************************************************************/
int			show_views(View* v)
{
	int			n = 0;
	
	for ( ; v; v = v->next, n++ )
		printf("%7.4f\t%7.4f\t%7.4f\t%7.2f\n", v->x, v->y, v->z, v->a*180.0/PI);
	
	return(n);
}

/************************************************************************
@Function: view_from_4_values
@Description:
	Returns a new view.
@Algorithm:
	A view is constructed from 4 floating point values.
@Arguments:
	double x			x axis value.
	double y			y axis value.
	double z			z axis value.
	double a			rotation angle (radians).
@Returns:
	View 				4-value view vector and roation angle.
**************************************************************************/
View		view_from_4_values(double x, double y, double z, double a)
{
	View		v = {NULL, x, y, z, a};
	
	view_normalize(&v);
	
	return(v);
}

/************************************************************************
@Function: view_vector
@Description:
	Returns the vector part of a view.
@Algorithm:
	The vector is not normalized.
@Arguments:
	View view			4-value view.
@Returns:
	float* 				3-value view vector.
**************************************************************************/
float*		view_vector(View view)
{
	float*		v = (float *) balloc(3*sizeof(float));
	
	v[0] = view.x;
	v[1] = view.y;
	v[2] = view.z;
	
	return(v);
}

/************************************************************************
@Function: view_vector_size
@Description:
	Returns the size of the vector part of a view.
@Algorithm:
	The vector is not normalized.
@Arguments:
	View view			4-value view.
@Returns:
	double 				vector size or length.
**************************************************************************/
double		view_vector_size(View view)
{
	double		size = sqrt(view.x*view.x + view.y*view.y + view.z*view.z);
	
	return(size);
}

/************************************************************************
@Function: view_normalize
@Description:
	Normalizes the vector part of a view.
@Algorithm:
	.
@Arguments:
	View* view			4-value view.
@Returns:
	double 				view vector size.
**************************************************************************/
double		view_normalize(View* view)
{
	double		size = view->x*view->x + view->y*view->y + view->z*view->z;
	
	if ( size < 1e-30 ) {
		view->x = view->y = 0;
		view->z = 1;
		size = 1;
	} else {
		size = sqrt(size);
		view->x /= size;
		view->y /= size;
		view->z /= size;
	}
	
	return(size);
}

/************************************************************************
@Function: view_negate
@Description:
	Inverts the view.
@Algorithm:
	.
@Arguments:
	View view			4-value view.
@Returns:
	View 				inverted/negated view.
**************************************************************************/
View		view_negate(View view)
{
	View		negview = {view.next, -view.x, -view.y, -view.z, -view.a};
	
	return(negview);
}

/************************************************************************
@Function: view_difference
@Description:
	Calculates the difference angle between two view vectors.
@Algorithm:
	Angle = arccos(x1*x2 + y1*y2 + z1*z2)
@Arguments:
	View view1		4-value view1: 3-value vector and angle.
	View view2		4-value view2: 3-value vector and angle.
@Returns:
	double 			difference angle.
**************************************************************************/
double		view_difference(View view1, View view2)
{
	double		dot = (double)view1.x*view2.x + (double)view1.y*view2.y + (double)view1.z*view2.z;
	
	if ( dot > 1 ) dot = 1;
	if ( dot < -1 ) dot = -1;
	
	double		angle = acos(dot);
	
	return(angle);
}

/*
	3-value vector functions:
	These use the Vector3 structure
*/

/************************************************************************
@Function: vector3_from_vectorint3
@Description:
	Converts a 3-value integer vector to a 3-value floating point vector.
@Algorithm:
	.
@Arguments:
	VectorInt3			3-value integer vector.
@Returns:
	Vector3 			3-value floating point vector.
**************************************************************************/
Vector3 	vector3_from_vectorint3(VectorInt3 intvec)
{
	Vector3		vector = {intvec.x, intvec.y, intvec.z};
	
	return(vector);
}

/************************************************************************
@Function: vectorint3_from_vector3
@Description:
	Converts a 3-value floating point vector to a 3-value integer vector.
@Algorithm:
	Each element is rounded to the nearest integer.
@Arguments:
	Vector3				3-value vector.
@Returns:
	VectorInt3 			3-value integer vector.
**************************************************************************/
VectorInt3 	vectorint3_from_vector3(Vector3 floatvec)
{
	VectorInt3		vector;
	
	vector.x = (int) (floatvec.x + 0.5);
	vector.y = (int) (floatvec.y + 0.5);
	vector.z = (int) (floatvec.z + 0.5);
	
	return(vector);
}

/************************************************************************
@Function: vector3_zero
@Description:
	Returns a 3-value zero vector.
@Algorithm:
	A new vector is created and all elements are set to zero.
@Arguments:
	
@Returns:
	Vector3 			{0,0,0}.
**************************************************************************/
Vector3 	vector3_zero()
{
	Vector3		vector = {0,0,0};
	
	return(vector);
}

/************************************************************************
@Function: vector3_from_3_values
@Description:
	Composes a 3-value vector from 3 floating point numbers.
@Algorithm:
	A new vector is created from the given 3 values.
@Arguments:
	double x			x-value.
	double y			y-value.
	double z			z-vlaue.
@Returns:
	Vector3 			new vector.
**************************************************************************/
Vector3 	vector3_from_3_values(double x, double y, double z)
{
	Vector3		vector = {x, y, z};
	
	return(vector);
}

/************************************************************************
@Function: vector3_random
@Description:
	Generates a 3-value vector with uniform random numbers.
@Algorithm:
	Each element is set to a random number as:
		r = random * (max - min) + min.
@Arguments:
	double min			lower bound.
	double max			upper bound.
@Returns:
	Vector3 			new vector.
**************************************************************************/
Vector3 	vector3_random(double min, double max)
{
	Vector3		vec;
	double		range = (max - min)/get_rand_max();
	
	vec.x = random()*range + min;
	vec.y = random()*range + min;
	vec.z = random()*range + min;
	
	return(vec);
}

/************************************************************************
@Function: vector3_random_gaussian
@Description:
	Generates a 3-value vector with a random gaussian amplitude.
@Algorithm:
	A random unit vector is generated using a uniform distribution from
	-1 to 1. An amplitude based on a random gaussian distribution is 
	generated and the unit vector scaled to it.
@Arguments:
	double avg 			average.
	double std 			standard deviation.
@Returns:
	Vector3 			new vector.
**************************************************************************/
Vector3 	vector3_random_gaussian(double avg, double std)
{
	Vector3		vec = vector3_random(-1, 1);
	double		amp = random_gaussian(avg, std);
	
	vec = vector3_normalize(vec);
	vec = vector3_scale(vec, amp);
	
	return(vec);
}

/************************************************************************
@Function: vector3_xy_random_gaussian
@Description:
	Generates a 3-value vector in the xy plane with a random gaussian amplitude.
@Algorithm:
	A random unit vector is generated using a uniform distribution from
	-1 to 1. An amplitude based on a random gaussian distribution is 
	generated and the unit vector scaled to it.
@Arguments:
	double avg 			average.
	double std 			standard deviation.
@Returns:
	Vector3 			new vector.
**************************************************************************/
Vector3 	vector3_xy_random_gaussian(double avg, double std)
{
	Vector3		vec = vector3_random(-1, 1);
	double		amp = random_gaussian(avg, std);
	
	vec.z = 0;
	vec = vector3_normalize(vec);
	vec = vector3_scale(vec, amp);
	
	return(vec);
}

/************************************************************************
@Function: vector3_normalize
@Description:
	Normalizes a 3-value vector.
@Algorithm:
	If the vector size is zero, the reference vector {0,0,1} is returned.
@Arguments:
	Vector3 vector		3-value vector.
@Returns:
	Vector3 			normalized 3-value vector.
**************************************************************************/
Vector3 	vector3_normalize(Vector3 vector)
{
	double		size = vector.x*vector.x + vector.y*vector.y + vector.z*vector.z;
	
	if ( size < 1e-30 ) {
		vector.x = vector.y = 0;
		vector.z = 1;
		size = 1;
	} else {
		size = sqrt(size);
		vector.x /= size;
		vector.y /= size;
		vector.z /= size;
	}
	
	return(vector);
}

/************************************************************************
@Function: vector3_floor
@Description:
	Gets a 3-value vector truncated to the given number of decimal places.
@Algorithm:
	.
@Arguments:
	Vector3 vector		3-value vector.
	int places			number of decimal places.
@Returns:
	Vector3 			truncated 3-value vector.
**************************************************************************/
Vector3 	vector3_floor(Vector3 vector, int places)
{
	vector.x = bfloor(vector.x, places);
	vector.y = bfloor(vector.y, places);
	vector.z = bfloor(vector.z, places);
	
	return(vector);
}

/************************************************************************
@Function: vector3_round
@Description:
	Rounds a 3-value vector to the given number of decimal places.
@Algorithm:
	.
@Arguments:
	Vector3 vector		3-value vector.
	int places			number of decimal places.
@Returns:
	Vector3 			rounded 3-value vector.
**************************************************************************/
Vector3 	vector3_round(Vector3 vector, int places)
{
	vector.x = bround(vector.x, places);
	vector.y = bround(vector.y, places);
	vector.z = bround(vector.z, places);
	
	return(vector);
}

/************************************************************************
@Function: vector3_remainder
@Description:
	Gets the remainder of a 3-value vector divided by a given value.
@Algorithm:
	.
@Arguments:
	Vector3 vector		3-value vector.
	int divisor			divisor value.
@Returns:
	Vector3 			remainder 3-value vector.
**************************************************************************/
Vector3 	vector3_remainder(Vector3 vector, int divisor)
{
	vector.x = fmod(vector.x, divisor);
	vector.y = fmod(vector.y, divisor);
	vector.z = fmod(vector.z, divisor);
	
	return(vector);
}

/************************************************************************
@Function: vector3_from_view
@Description:
	Gets the vector from a view.
@Algorithm:
	.
@Arguments:
	View view			4-value view.
@Returns:
	Vector3 			3-value vector.
**************************************************************************/
Vector3 	vector3_from_view(View view)
{
	Vector3 		vector;
	
	vector.x = view.x;
	vector.y = view.y;
	vector.z = view.z;
	
	return(vector);
}

/************************************************************************
@Function: vector3_length
@Description:
	Calculates the length of a 3-value vector.
@Algorithm:
	length = sqrt(x1*x2 + y1*y2 + z1*z2)
@Arguments:
	Vector3 vector		3-value vector.
@Returns:
	double	 			vector length.
**************************************************************************/
double		vector3_length(Vector3 vector)
{
  return sqrt( vector3_length2(vector) );
}

/************************************************************************
@Function: vector3_length2
@Description:
	Calculates the squared length of a 3-value vector.
@Algorithm:
	length^2 = x1*x2 + y1*y2 + z1*z2
@Arguments:
	Vector3 vector		3-value vector.
@Returns:
	double	 			vector squared length.
**************************************************************************/
double		vector3_length2(Vector3 vector)
{
  return (vector.x*vector.x + vector.y*vector.y + vector.z*vector.z);
}

/************************************************************************
@Function: vector3_scale
@Description:
	Multiplies a 3-value vector with a scalar value.
@Algorithm:
	vector = scalar*vector
@Arguments:
	Vector3 vector		3-value vector.
	double scale		scalar value.
@Returns:
	Vector3	 			scaled 3-value vector.
**************************************************************************/
Vector3		vector3_scale(Vector3 vector, double scale)
{
	Vector3 	new_vector = {scale*vector.x, scale*vector.y, scale*vector.z};
	
	return(new_vector);
}

/************************************************************************
@Function: vector3_multiply
@Description:
	Multiplies two 3-value vectors, element by element.
@Algorithm:
	vector.x = vec1.x*vec2.x
@Arguments:
	Vector3 vec1		3-value vector.
	Vector3 vec2		3-value vector.
@Returns:
	Vector3	 			new 3-value vector.
**************************************************************************/
Vector3		vector3_multiply(Vector3 vec1, Vector3 vec2)
{
	Vector3 	vector = {vec1.x*vec2.x, vec1.y*vec2.y, vec1.z*vec2.z};
	
	return(vector);
}

/************************************************************************
@Function: vector3_divide
@Description:
	Divides a 3-value vector by a another, element by element.
@Algorithm:
	vector.x = vec1.x/vec2.x
	All elements of the second vector must be non-zero.
@Arguments:
	Vector3 vec1		3-value vector.
	Vector3 vec2		3-value vector (divisor).
@Returns:
	Vector3	 			new 3-value vector.
**************************************************************************/
Vector3		vector3_divide(Vector3 vec1, Vector3 vec2)
{
	Vector3 	vector = {vec1.x/vec2.x, vec1.y/vec2.y, vec1.z/vec2.z};
	
	return(vector);
}

/************************************************************************
@Function: vector3_square
@Description:
	Squares each value of a 3-value vector.
@Algorithm:
	vector = vector^2
@Arguments:
	Vector3 vector		3-value vector.
@Returns:
	Vector3	 			squared vector.
**************************************************************************/
Vector3		vector3_square(Vector3 vector)
{
	Vector3 	new_vector = {vector.x*vector.x, vector.y*vector.y, vector.z*vector.z};
	
	return(new_vector);
}

/************************************************************************
@Function: vector3_square_root
@Description:
	Returns the square root of each value of a 3-value vector.
@Algorithm:
	vector = sqrt(vector)
@Arguments:
	Vector3 vector		3-value vector.
@Returns:
	Vector3	 			square root vector.
**************************************************************************/
Vector3		vector3_square_root(Vector3 vector)
{
	Vector3 	new_vector = {sqrt(vector.x), sqrt(vector.y), sqrt(vector.z)};
	
	return(new_vector);
}

/************************************************************************
@Function: vector3_negate
@Description:
	Negates 3-value vectors.
@Algorithm:
	vector = -vector
@Arguments:
	Vector3 vec			3-value vector.
@Returns:
	Vector3	 			the negated vector.
**************************************************************************/
Vector3		vector3_negate(Vector3 vec)
{
	Vector3 	negvec = {-vec.x, -vec.y, -vec.z};
	
	return(negvec);
}

/************************************************************************
@Function: vector3_add
@Description:
	Adds two 3-value vectors.
@Algorithm:
	vector = vector1 + vector2
@Arguments:
	Vector3 vec1		3-value first vector.
	Vector3 vec2		3-value second vector.
@Returns:
	Vector3	 			the vector sum.
**************************************************************************/
Vector3		vector3_add(Vector3 vec1, Vector3 vec2)
{
	Vector3 	vec = {vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z};
	
	return(vec);
}

/************************************************************************
@Function: vector3_subtract
@Description:
	Subtracts one 3-value vector from another.
@Algorithm:
	vector = vector1 - vector2
@Arguments:
	Vector3 vec1		3-value first vector.
	Vector3 vec2		3-value second vector.
@Returns:
	Vector3	 			the vector difference.
**************************************************************************/
Vector3		vector3_subtract(Vector3 vec1, Vector3 vec2)
{
	Vector3 	vec = {vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z};
	
	return(vec);
}

/************************************************************************
@Function: vector3_add_scalar
@Description:
	Adds a scalar value to a 3-value vector.
@Algorithm:
	vector = vector + scalar
@Arguments:
	Vector3 vec			3-value vector.
	double s				scalar.
@Returns:
	Vector3	 			the new vector.
**************************************************************************/
Vector3		vector3_add_scalar(Vector3 vec, double s)
{
	Vector3 	newvec = {vec.x + s, vec.y + s, vec.z + s};
	
	return(newvec);
}

/************************************************************************
@Function: vector3_scalar_product
@Description:
	Calculates the 3-value vector product from two 3-value vectors.
@Algorithm:
	Scalar product:
		s = x1*x2 + y1*y2 + z1*z2
@Arguments:
	Vector3 v1			3-value first vector.
	Vector3 v2			3-value second vector.
@Returns:
	double	 			scalar result.
**************************************************************************/
double	 	vector3_scalar_product(Vector3 v1, Vector3 v2)
{
	double		scalar = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	
	return(scalar);
}

/************************************************************************
@Function: vector3_vector_product
@Description:
	Calculates the 3-value vector product from two 3-value vectors.
@Algorithm:
	Vector product:
		x = y1*z2 - z1*y2
		y = z1*x2 - x1*z2
		z = x1*y2 - y1*x2
@Arguments:
	Vector3 v1			3-value first vector.
	Vector3 v2			3-value second vector.
@Returns:
	Vector3 			3-value result vector.
**************************************************************************/
Vector3 	vector3_vector_product(Vector3 v1, Vector3 v2)
{
	Vector3		v;
	
	v.x = v1.y*v2.z - v1.z*v2.y;
	v.y = v1.z*v2.x - v1.x*v2.z;
	v.z = v1.x*v2.y - v1.y*v2.x;
	
	return(v);
}

/************************************************************************
@Function: vector3_angle
@Description:
	Calculates the angle between two 3-value vectors.
@Algorithm:
	The arc cosine of the scalar product of the normalized vectors:
		a = acos(x1*x2 + y1*y2 + z1*z2)
@Arguments:
	Vector3 v1			3-value first vector.
	Vector3 v2			3-value second vector.
@Returns:
	double	 			angle between vectors.
**************************************************************************/
double	 	vector3_angle(Vector3 v1, Vector3 v2)
{
	v1 = vector3_normalize(v1);
	v2 = vector3_normalize(v2);
	
	double		prod = vector3_scalar_product(v1, v2);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG vector3_angle: prod = %g\n", prod);
	
	if ( prod > 1 ) prod = 1;
	if ( prod < -1 ) prod = -1;
	
	return(acos(prod));
}

/************************************************************************
@Function: vector3_matrix3_multiply
@Description:
	Multiplies a 3x3 matrix by a 3-value vector.
@Algorithm:
	.
@Arguments:
	Matrix3 mat 		3x3 matrix.
	Vector3 vec 		3-value vector.
@Returns:
	Vector3 			resultant 3-value vector.
**************************************************************************/
Vector3		vector3_matrix3_multiply(Matrix3 mat, Vector3 vec) 
{ 
	Vector3		new_vec;
	
	new_vec.x = mat.r00*vec.x + mat.r01*vec.y + mat.r02*vec.z;
	new_vec.y = mat.r10*vec.x + mat.r11*vec.y + mat.r12*vec.z;
	new_vec.z = mat.r20*vec.x + mat.r21*vec.y + mat.r22*vec.z;
	 
	return(new_vec); 
} 

/************************************************************************
@Function: vector3_calculate_normal
@Description:
	Calculates the normal to a plane defined by 3 vectors or points.
@Algorithm:
	The normal to a plane is calculated as:
		edge1 = vector1 - vector2
		edge2 = vector1 - vector3
		normal = edge1 X edge2.
	The resultant vector is normalized.
@Arguments:
	Vector3 v1	 		3-value vector.
	Vector3 v2	 		3-value vector.
	Vector3 v3	 		3-value vector.
@Returns:
	Vector3 			resultant 3-value normal vector.
**************************************************************************/
Vector3 	vector3_calculate_normal(Vector3 v1, Vector3 v2, Vector3 v3)
{
	Vector3 	edge1 = vector3_subtract(v1, v2);
	Vector3 	edge2 = vector3_subtract(v1, v3);

	Vector3 	normal = vector3_vector_product(edge1, edge2);
	
	normal = vector3_normalize(normal);
	
	return(normal);
}

/************************************************************************
@Function: vector3_min
@Description:
	Compares each pair of elements of two vectors and returns the minimum.
@Algorithm:
	.
@Arguments:
	Vector3 vector1		3-value vector.
	Vector3 vector2		3-value vector.
@Returns:
	Vector3	 			resultant vector.
**************************************************************************/
Vector3		vector3_min(Vector3 vector1, Vector3 vector2)
{
	Vector3		newvector = vector1;
	
	if ( vector1.x > vector2.x ) newvector.x = vector2.x;
	if ( vector1.y > vector2.y ) newvector.y = vector2.y;
	if ( vector1.z > vector2.z ) newvector.z = vector2.z;
	
	return (newvector);
}

/************************************************************************
@Function: vector3_max
@Description:
	Compares each pair of elements of two vectors and returns the maximum.
@Algorithm:
	.
@Arguments:
	Vector3 vector1		3-value vector.
	Vector3 vector2		3-value vector.
@Returns:
	Vector3	 			resultant vector.
**************************************************************************/
Vector3		vector3_max(Vector3 vector1, Vector3 vector2)
{
	Vector3		newvector = vector1;
	
	if ( vector1.x < vector2.x ) newvector.x = vector2.x;
	if ( vector1.y < vector2.y ) newvector.y = vector2.y;
	if ( vector1.z < vector2.z ) newvector.z = vector2.z;
	
	return (newvector);
}

/************************************************************************
@Function: vector3_scalar_min
@Description:
	Compares each element to the scalar and sets it to the minimum.
@Algorithm:
	.
@Arguments:
	Vector3 vector		3-value vector.
	double scalar		scalar value to compare with
@Returns:
	Vector3	 			resultant vector.
**************************************************************************/
Vector3		vector3_scalar_min(Vector3 vector, double scalar)
{
	Vector3		newvector = {scalar, scalar, scalar};
	
	if ( vector.x < scalar ) newvector.x = vector.x;
	if ( vector.y < scalar ) newvector.y = vector.y;
	if ( vector.z < scalar ) newvector.z = vector.z;
	
	return (newvector);
}

/************************************************************************
@Function: vector3_scalar_max
@Description:
	Compares each element to the scalar and sets it to the maximum.
@Algorithm:
	.
@Arguments:
	Vector3 vector		3-value vector.
	double scalar		scalar value to compare with
@Returns:
	Vector3	 			resultant vector.
**************************************************************************/
Vector3		vector3_scalar_max(Vector3 vector, double scalar)
{
	Vector3		newvector = {scalar, scalar, scalar};
	
	if ( vector.x > scalar ) newvector.x = vector.x;
	if ( vector.y > scalar ) newvector.y = vector.y;
	if ( vector.z > scalar ) newvector.z = vector.z;
	
	return (newvector);
}

/************************************************************************
@Function: vector3_scalar_range
@Description:
	Sets each element to within the scalar range.
@Algorithm:
	If an element is smaller than the minimum, it is set to the minimum.
	If an element is larger than the maximum, it is set to the maximum.
@Arguments:
	Vector3 vector		3-value vector.
	double min			scalar minimum.
	double max			scalar maximum.
@Returns:
	Vector3	 			resultant vector.
**************************************************************************/
Vector3		vector3_scalar_range(Vector3 vec, double min, double max)
{
	if ( max < min ) return(vec);
	
	if ( vec.x < min ) vec.x = min;
	if ( vec.x > max ) vec.x = max;
	if ( vec.y < min ) vec.y = min;
	if ( vec.y > max ) vec.y = max;
	if ( vec.z < min ) vec.z = min;
	if ( vec.z > max ) vec.z = max;
	
	return(vec);
}

/*
	Periodic boundary vector functions:
*/

/************************************************************************
@Function: vector3_set_PBC
@Description:
	Limits coordinates to within a given box by wrapping.
@Algorithm:
	If a coordinate lies outside the box, it is modified by adding or
	subtracting the box size until it lies within the box.
	The box origin is (0,0,0).
	Note: If the coordinates are far out (several orders of magnitude)
		this function may take a lot of time!!!
@Arguments:
	Vector3 coord		coordinate vector.
	Vector3 box			size vector.
@Returns:
	Vector3	 			resultant vector.
**************************************************************************/
Vector3		vector3_set_PBC(Vector3 coord, Vector3 box)
{
	while ( coord.x < 0 ) coord.x += box.x;
	while ( coord.y < 0 ) coord.y += box.y;
	while ( coord.z < 0 ) coord.z += box.z;
	while ( coord.x >= box.x ) coord.x -= box.x;
	while ( coord.y >= box.y ) coord.y -= box.y;
	while ( coord.z >= box.z ) coord.z -= box.z;
	
	return(coord);
}

/************************************************************************
@Function: vector3_difference_PBC
@Description:
	Finds the difference between coordinates across periodic boundaries.
@Algorithm:
	The shortest distance between two vectors is determined by comparing
	the direct difference with that considering periodic boundaries.
	The origin of the periodic box is assumed to be (0,0,0).
@Arguments:
	Vector3 v1			first set of coordinates.
	Vector3 v2			second set of coordinates.
	Vector3 box			periodic box size.
@Returns:
	Vector3	 			resultant vector.
**************************************************************************/
Vector3		vector3_difference_PBC(Vector3 v1, Vector3 v2, Vector3 box)
{
	Vector3			d = vector3_subtract(v1, v2);
	
	if ( box.x - d.x < d.x ) d.x -= box.x;
	if ( box.x + d.x < -d.x ) d.x += box.x;
	if ( box.y - d.y < d.y ) d.y -= box.y;
	if ( box.y + d.y < -d.y ) d.y += box.y;
	if ( box.z - d.z < d.z ) d.z -= box.z;
	if ( box.z + d.z < -d.z ) d.z += box.z;
	
	return(d);
}

/*
	3-value integer vector functions:
	These use the VectorInt3 structure
*/

/************************************************************************
@Function: vectorint3_from_3_values
@Description:
	Composes an integer 3-value vector from 3 integer numbers.
@Algorithm:
	A new vector is created from the given 3 values.
@Arguments:
	int x				x-value.
	int y				y-value.
	int z				z-vlaue.
@Returns:
	VectorInt3 			new integer vector.
**************************************************************************/
VectorInt3 	vectorint3_from_3_values(int x, int y, int z)
{
	VectorInt3		vector = {x, y, z};
	
	return(vector);
}

/************************************************************************
@Function: vectorint3_length
@Description:
	Calculates the length of a 3-value integer vector.
@Algorithm:
	length = sqrt(x1*x2 + y1*y2 + z1*z2)
@Arguments:
	VectorInt3 vector	3-value integer vector.
@Returns:
	double	 			vector length.
**************************************************************************/
double		vectorint3_length(VectorInt3 vector)
{
	return sqrt( vectorint3_length2(vector) );
}

/************************************************************************
@Function: vectorint3_length2
@Description:
	Calculates the squared length of a 3-value integer vector.
@Algorithm:
	length^2 = x1*x2 + y1*y2 + z1*z2
@Arguments:
	VectorInt3 vector	3-value integer vector.
@Returns:
	double	 			vector squared length.
**************************************************************************/
double		vectorint3_length2(VectorInt3 vector)
{
	return (vector.x*vector.x + vector.y*vector.y + vector.z*vector.z);
}

/************************************************************************
@Function: vectorint3_add
@Description:
	Adds two 3-value integer vectors.
@Algorithm:
	vector = vector1 + vector2
@Arguments:
	VectorInt3 vec1		3-value first integer vector.
	VectorInt3 vec2		3-value second integer vector.
@Returns:
	VectorInt3			the vector sum.
**************************************************************************/
VectorInt3  vectorint3_add(VectorInt3 vec1, VectorInt3 vec2)
{
	VectorInt3 	vec = {vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z};
	
	return(vec);
}

/************************************************************************
@Function: vectorint3_subtract
@Description:
	Subtracts one 3-value integer vector from another.
@Algorithm:
	vector = vector1 - vector2
@Arguments:
	VectorInt3 vec1		3-value first vector.
	VectorInt3 vec2		3-value second vector.
@Returns:
	VectorInt3	 		the vector difference.
**************************************************************************/
VectorInt3	vectorint3_subtract(VectorInt3 vec1, VectorInt3 vec2)
{
	VectorInt3 	vec = {vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z};
	
	return(vec);
}

/************************************************************************
@Function: vectorint3_add_scalar
@Description:
	Adds a scalar value to a 3-value integer vector.
@Algorithm:
	vector = vector + scalar
@Arguments:
	VectorInt3 vec		3-value integer vector.
	int s				scalar.
@Returns:
	VectorInt3			the new vector.
**************************************************************************/
VectorInt3  vectorint3_add_scalar(VectorInt3 vec, int s)
{
	VectorInt3 	newvec = {vec.x + s, vec.y + s, vec.z + s};
	
	return(newvec);
}

/************************************************************************
@Function: vectorint3_scale
@Description:
	Multiplies 3-value integer vector with a scalar value.
@Algorithm:
	new_vector = scale*vector
	Each element is rounded off to the nearest integer.
@Arguments:
	VectorInt3 vec		3-value first vector.
	double scale			multiplier.
@Returns:
	VectorInt3	 		the vector difference.
**************************************************************************/
VectorInt3	vectorint3_scale(VectorInt3 vec, double scale)
{
	VectorInt3 	newvec = {(int) (scale*vec.x+0.5), 
		(int) (scale*vec.y+0.5), (int) (scale*vec.z+0.5)};
	
	return(newvec);
}

/************************************************************************
@Function: vectorint3_square
@Description:
	Squares each value of a 3-value integer vector.
@Algorithm:
	vector = vector^2
@Arguments:
	VectorInt3 vector		3-value vector.
@Returns:
	VectorInt3	 			squared vector.
**************************************************************************/
VectorInt3	vectorint3_square(VectorInt3 vector)
{
	VectorInt3 	new_vector = {vector.x*vector.x, vector.y*vector.y, vector.z*vector.z};
	
	return(new_vector);
}

/************************************************************************
@Function: vectorint3_min
@Description:
	Compares each pair of elements of two integer vectors and returns the minimum.
@Algorithm:
	.
@Arguments:
	VectorInt3 vector1		3-value integer vector.
	VectorInt3 vector2		3-value integer vector.
@Returns:
	VectorInt3	 			resultant integer vector.
**************************************************************************/
VectorInt3	vectorint3_min(VectorInt3 vector1, VectorInt3 vector2)
{
	VectorInt3		newvector = vector1;
	
	if ( vector1.x > vector2.x ) newvector.x = vector2.x;
	if ( vector1.y > vector2.y ) newvector.y = vector2.y;
	if ( vector1.z > vector2.z ) newvector.z = vector2.z;
	
	return (newvector);
}

/************************************************************************
@Function: vectorint3_max
@Description:
	Compares each pair of elements of two integer vectors and returns the maximum.
@Algorithm:
	.
@Arguments:
	VectorInt3 vector1		3-value integer vector.
	VectorInt3 vector2		3-value integer vector.
@Returns:
	VectorInt3	 			resultant integer vector.
**************************************************************************/
VectorInt3	vectorint3_max(VectorInt3 vector1, VectorInt3 vector2)
{
	VectorInt3		newvector = vector1;
	
	if ( vector1.x < vector2.x ) newvector.x = vector2.x;
	if ( vector1.y < vector2.y ) newvector.y = vector2.y;
	if ( vector1.z < vector2.z ) newvector.z = vector2.z;
	
	return (newvector);
}

/************************************************************************
@Function: vectorint3_scalar_min
@Description:
	Compares each element to the scalar and sets it to the minimum.
@Algorithm:
	.
@Arguments:
	VectorInt3 vector	3-value integer vector.
	int scalar			scalar value to compare with
@Returns:
	VectorInt3	 		resultant integer vector.
**************************************************************************/
VectorInt3	vectorint3_scalar_min(VectorInt3 vector, int scalar)
{
	VectorInt3	newvector = {scalar, scalar, scalar};
	
	if ( vector.x < scalar ) newvector.x = vector.x;
	if ( vector.y < scalar ) newvector.y = vector.y;
	if ( vector.z < scalar ) newvector.z = vector.z;
	
	return (newvector);
}

/************************************************************************
@Function: vectorint3_scalar_max
@Description:
	Compares each element to the scalar and sets it to the maximum.
@Algorithm:
	.
@Arguments:
	VectorInt3 vector	3-value integer vector.
	int scalar			scalar value to compare with
@Returns:
	VectorInt3	 		resultant integer vector.
**************************************************************************/
VectorInt3	vectorint3_scalar_max(VectorInt3 vector, int scalar)
{
	VectorInt3	newvector = {scalar, scalar, scalar};
	
	if ( vector.x > scalar ) newvector.x = vector.x;
	if ( vector.y > scalar ) newvector.y = vector.y;
	if ( vector.z > scalar ) newvector.z = vector.z;
	
	return (newvector);
}

/************************************************************************
@Function: vectorlong3_from_3_values
@Description:
	Composes a long integer 3-value vector from 3 long integer numbers.
@Algorithm:
	A new vector is created from the given 3 values.
@Arguments:
	long x				x-value.
	long y				y-value.
	long z				z-vlaue.
@Returns:
	VectorLong3 			new integer vector.
**************************************************************************/
VectorLong3 vectorlong3_from_3_values(long x, long y, long z)
{
	VectorLong3		vector = {x, y, z};
	
	return(vector);
}


/*
	3x3 Matrix functions:
	These use the Matrix3 structure
*/

/************************************************************************
@Function: matrix3_show
@Description:
	Print out a matrix.
@Algorithm:
	Matrix elements less than 1e-6 are set to zero for display.
@Arguments:
	Matrix3			matrix.
@Returns:
	int 			0.
**************************************************************************/
int 		matrix3_show(Matrix3 mat)
{
	printf("%10.4f %10.4f %10.4f\n", mat.r00, mat.r01, mat.r02);
	printf("%10.4f %10.4f %10.4f\n", mat.r10, mat.r11, mat.r12);
	printf("%10.4f %10.4f %10.4f\n", mat.r20, mat.r21, mat.r22);
	
	return(0);
}

/************************************************************************
@Function: matrix3_init
@Description:
	Initializes a unit rotation matrix.
@Algorithm:
	.
@Arguments:
	.
@Returns:
	Matrix3 			3x3 unit matrix.
**************************************************************************/
Matrix3		matrix3_init() 
{
	Matrix3 	mat = {1,0,0,0,1,0,0,0,1};
 
	return(mat); 
} 
 
/************************************************************************
@Function: matrix3_determinant
@Description:
	Calculates a 3x3 matrix determinant.
@Algorithm:
	.
@Arguments:
	Matrix3 		a 3x3 matrix.
@Returns:
	double 			determinant.
**************************************************************************/
double		matrix3_determinant(Matrix3 mat) 
{ 
	double		det = 		mat.r00*(mat.r11*mat.r22-mat.r12*mat.r21) 
						- 	mat.r01*(mat.r10*mat.r22-mat.r12*mat.r20) 
						+	mat.r02*(mat.r10*mat.r21-mat.r11*mat.r20); 
 
	return(det); 
} 
 
/************************************************************************
@Function: matrix3_transpose
@Description:
	Calculates the transpose of a 3x3 rotation matrix.
@Algorithm:
	The transpose of a unit matrix is also its inverse.
@Arguments:
	Matrix3	mat			3x3 matrix.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_transpose(Matrix3 mat) 
{
	Matrix3 	tmat;
	
	tmat.r00 = mat.r00;
	tmat.r01 = mat.r10;
	tmat.r02 = mat.r20;
	tmat.r10 = mat.r01;
	tmat.r11 = mat.r11;
	tmat.r12 = mat.r21;
	tmat.r20 = mat.r02;
	tmat.r21 = mat.r12;
	tmat.r22 = mat.r22;
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_transpose:\n");
		matrix3_show(tmat);
	}
 
	return(tmat); 
}	 

/************************************************************************
@Function: matrix3_scalar_multiply
@Description:
	Multiplies a 3x3 rotation matrix by a scalar value.
@Algorithm:
	The matrix is multiplied row by row:
		new_mat(i,j) = scalar*mat(i,j)
@Arguments:
	Matrix3	mat			3x3 matrix.
	double scalar 		single scalar value.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_scalar_multiply(Matrix3 mat, double scalar) 
{
	Matrix3 	new_mat;
	
	new_mat.r00 = scalar*mat.r00;
	new_mat.r01 = scalar*mat.r01;
	new_mat.r02 = scalar*mat.r02;
	new_mat.r10 = scalar*mat.r10;
	new_mat.r11 = scalar*mat.r11;
	new_mat.r12 = scalar*mat.r12;
	new_mat.r20 = scalar*mat.r20;
	new_mat.r21 = scalar*mat.r21;
	new_mat.r22 = scalar*mat.r22;
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_scalar_multiply:\n");
		matrix3_show(new_mat);
	}
 
	return(new_mat); 
}	 
 
/************************************************************************
@Function: matrix3_vector3_multiply
@Description:
	Multiplies a 3x3 rotation matrix by a 3-value vector.
@Algorithm:
	The matrix is multiplied row by row:
		new_mat(i,j) = mat(i,j)*vec(i)
@Arguments:
	Matrix3	mat			3x3 matrix.
	Vector3 vec 		3-value vector.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_vector3_multiply(Matrix3 mat, Vector3 vec) 
{
	Matrix3 	new_mat;
	
	new_mat.r00 = mat.r00*vec.x;
	new_mat.r01 = mat.r01*vec.x;
	new_mat.r02 = mat.r02*vec.x;
	new_mat.r10 = mat.r10*vec.y;
	new_mat.r11 = mat.r11*vec.y;
	new_mat.r12 = mat.r12*vec.y;
	new_mat.r20 = mat.r20*vec.z;
	new_mat.r21 = mat.r21*vec.z;
	new_mat.r22 = mat.r22*vec.z;
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_vector3_multiply:\n");
		matrix3_show(new_mat);
	}
 
	return(new_mat); 
}	 
 
/************************************************************************
@Function: matrix3_divide
@Description:
	Divides a 3x3 rotation matrix by a 3-value vector.
@Algorithm:
	The matrix is divided row by row:
		new_mat(i,j) = mat(i,j)/vec(i)
@Arguments:
	Matrix3	mat			3x3 matrix.
	Vector3 vec 		3-value vector.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_divide(Matrix3 mat, Vector3 vec) 
{
	Matrix3 	new_mat;
	
	new_mat.r00 = mat.r00/vec.x;
	new_mat.r01 = mat.r01/vec.x;
	new_mat.r02 = mat.r02/vec.x;
	new_mat.r10 = mat.r10/vec.y;
	new_mat.r11 = mat.r11/vec.y;
	new_mat.r12 = mat.r12/vec.y;
	new_mat.r20 = mat.r20/vec.z;
	new_mat.r21 = mat.r21/vec.z;
	new_mat.r22 = mat.r22/vec.z;
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_divide:\n");
		matrix3_show(new_mat);
	}
 
	return(new_mat); 
}	 
 
/************************************************************************
@Function: matrix3_from_view
@Description:
	Calculates a 3x3 rotation matrix from a view.
@Algorithm:
	The view vector is converted to a quaternion and converted to a 3x3 matrix.
@Arguments:
	View view	 		4-value view.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_from_view(View view) 
{
	Quaternion	q = quaternion_from_view(view);
	
	// Convert to a rotation matrix 
	Matrix3 	mat = matrix3_from_quaternion(q);
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_from_vector_and_angle:\n");
		matrix3_show(mat);
	}
 
	return(mat); 
}	 
 
/************************************************************************
@Function: matrix3_from_angle_and_axis
@Description:
	Calculates a 3x3 rotation matrix from a given axis and angle.
@Algorithm:
	The axis vector is normalized first, then converted to a quaternion 
	representing the rotation around the axis by the given angle.
	The quaternion is then converted to a 3x3 matrix.
@Arguments:
	double angle 		rotation angle.
	float* axis 		three-value axis of rotation.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_from_angle_and_axis(double angle, float* axis) 
{ 
	Quaternion	q = quaternion_from_angle_and_axis(angle, axis);
	
	// Convert to a rotation matrix 
	Matrix3 	mat = matrix3_from_quaternion(q);
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_from_angle_and_axis:\n");
		matrix3_show(mat);
	}
 
	return(mat); 
}	 
 
/************************************************************************
@Function: matrix3_from_angle_and_axis3
@Description:
	Calculates a 3x3 rotation matrix from a given axis and angle.
@Algorithm:
	The axis vector is normalized first, then converted to a quaternion 
	representing the rotation around the axis by the given angle.
	The quaternion is then converted to a 3x3 matrix.
@Arguments:
	double angle 		rotation angle.
	Vector3 axis 		3-value axis of rotation.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_from_angle_and_axis3(double angle, Vector3 axis) 
{ 
	Quaternion	q = quaternion_from_angle_and_axis3(angle, axis);
	
	// Convert to a rotation matrix 
	Matrix3 	mat = matrix3_from_quaternion(q);
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_from_angle_and_axis3:\n");
		matrix3_show(mat);
	}
 
	return(mat); 
}	 
 
/************************************************************************
@Function: matrix3_from_tilt_and_axis
@Description:
	Calculates a 3x3 rotation matrix from given tilt and axis angles (2D).
@Algorithm:
	The tilt of a 2D plane around an axis defined by an angle with respect to
	the x-axis is converted to quaternions:
		q = {cos(tilt/2),cos(axis)*sin(tilt/2),sin(axis)*sin(tilt/2),0}
	The quaternions are then converted to a 3x3 matrix.
@Arguments:
	double tilt		tilt angle of 2D plane.
	double axis		tilt axis angle with respect to the x-axis.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_from_tilt_and_axis(double tilt, double axis) 
{ 
	Quaternion	q; 
	 
	// Convert to quaternions 
    q.s = cos(tilt/2.0L); 
    q.x = cos(axis)*sin(tilt/2.0L); 
    q.y = sin(axis)*sin(tilt/2.0L); 
    q.z = 0; 
	 
	// Convert to a rotation matrix 
	Matrix3 	mat = matrix3_from_quaternion(q);
	 
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_from_tilt_and_axis:\n");
		matrix3_show(mat);
	}
 
	return(mat); 
} 
 
/************************************************************************
@Function: matrix3_from_two_vectors
@Description:
	Calculates a 3x3 rotation matrix rotating from one vector to another.
@Algorithm:
	The vectors are normalized first.
	The rotation angle is derived from the scalar product:
		angle = arccos(scalar);
	The vector product gives the rotation axis.
	The axis vector and angle are converted to quaternions and then to a 
			3x3 matrix.
@Arguments:
	float* from_vec 	3-value initial vector.
	float* to_vec		3-value final vector.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_from_two_vectors(float* from_vec, float* to_vec) 
{
	// Convert to quaternion
	Quaternion	q = quaternion_from_two_vectors(from_vec, to_vec);
	
	// Convert to a rotation matrix 
	Matrix3 	mat = matrix3_from_quaternion(q);
	 
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_from_two_vectors:\n");
		matrix3_show(mat);
	}
 
	return(mat); 
}	 
 
/************************************************************************
@Function: matrix3_from_two_vectors3
@Description:
	Calculates a 3x3 rotation matrix rotating from one vector to another.
@Algorithm:
	The vectors are normalized first.
	The rotation angle is derived from the scalar product:
		angle = arccos(scalar);
	The vector product gives the rotation axis.
	The axis vector and angle are converted to quaternions and then to a 
			3x3 matrix.
@Arguments:
	Vector3 from_vec 	3-value initial vector.
	Vector3 to_vec		3-value final vector.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_from_two_vectors3(Vector3 from_vec, Vector3 to_vec) 
{
	Quaternion	q = quaternion_from_two_vectors3(from_vec, to_vec);
	
	Matrix3 	mat = matrix3_from_quaternion(q);
	 
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_from_two_vectors3:\n");
		matrix3_show(mat);
	}
 
	return(mat); 
}	 
 
/************************************************************************
@Function: matrix3_from_euler
@Description:
	Calculates a 3x3 rotation matrix from three Euler angles.
@Algorithm:
	The definitions of angles are based on SPIDER conventions: 
	1.	Rotation around z with an angle of phi 
	2.	Rotation around y with an angle of theta 
	3.	Rotation around z with an angle of psi 
	Convert each individual angle to quaternions.
	Multiply the quaternions to get the combined set of quaternions.
	The quaternions are then converted to a 3x3 matrix.
@Arguments:
	Euler euler			Euler angels: psi, theta, phi.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_from_euler(Euler euler) 
{
	return(matrix3_from_psi_theta_phi(euler.psi, euler.theta, euler.phi));
}

/************************************************************************
@Function: matrix3_from_psi_theta_phi
@Description:
	Calculates a 3x3 rotation matrix from three Euler angles.
@Algorithm:
	The definitions of angles are based on SPIDER conventions: 
	1.	Rotation around z with an angle of phi 
	2.	Rotation around y with an angle of theta 
	3.	Rotation around z with an angle of psi 
	Convert each individual angle to quaternions.
	Multiply the quaternions to get the combined set of quaternions.
	The quaternions are then converted to a 3x3 matrix.
@Arguments:
	double psi	 	rotation around z with an angle of psi.
	double theta 	rotation around y with an angle of theta.
	double phi		rotation around z with an angle of phi.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_from_psi_theta_phi(double psi, double theta, double phi) 
{ 
//	Quaternion	q2, q3;
		 
	// Convert phi to quaternions = rotation around z
	Quaternion	q1 = {cos(phi/2.0L),0,0,sin(phi/2.0L)};
	
	// Convert theta to quaternions = rotation around y 
	Quaternion	q2 = {cos(theta/2.0L),0,sin(theta/2.0L),0};
	 	 
	// Convert psi to quaternions = rotation around z 
	Quaternion	q3 = {cos(psi/2.0L),0,0,sin(psi/2.0L)};
	 
	Quaternion	q0 = quaternion_multiply(q1, q2); 
	Quaternion	q = quaternion_multiply(q0, q3); 
	 
	// Convert to a rotation matrix 
	Matrix3 	mat = matrix3_from_quaternion(q);
	 
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_from_psi_theta_phi:\n");
		matrix3_show(mat);
	}
 
	return(mat); 
} 
 
/************************************************************************
@Function: matrix3_from_theta_phi_omega
@Description:
	Calculates a 3x3 rotation matrix from three Euler angles.
@Algorithm:
	The definitions of angles are based on Tim Baker's MAP_PROJECT function: 
	1.	Rotation around z with an angle of phi 
	2.	Rotation around y with an angle of theta 
	3.	Rotation around z with an angle of -omega 
	Convert each individual angle to quaternions.
	Multiply the quaternions to get the combined set of quaternions.
	The quaternions are then converted to a 3x3 matrix.
@Arguments:
	double theta 	rotation around y with an angle of theta.
	double phi		rotation around z with an angle of phi.
	double omega 	rotation around z with an angle of -omega.
@Returns:
	Matrix3 			3x3 matrix.
**************************************************************************/
Matrix3		matrix3_from_theta_phi_omega(double theta, double phi, double omega) 
{ 
	Matrix3		mat = matrix3_from_psi_theta_phi(-omega, theta, phi);
	 
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_from_theta_phi_omega:\n");
		matrix3_show(mat);
	}
 
	return(mat); 
} 
  
/************************************************************************
@Function: matrix3_multiply
@Description:
	Multiplies two 3x3 matrices.
@Algorithm:
	.
@Arguments:
	Matrix3 mat1 	first 3x3 matrix.
	Matrix3 mat2 	second 3x3 matrix.
@Returns:
	Matrix3 			resultant 3x3 matrix.
**************************************************************************/
Matrix3		matrix3_multiply(Matrix3 mat1, Matrix3 mat2) 
{
	Matrix3 	mat = {0,0,0,0,0,0,0,0,0};
	
	double*		r1 = (double *) &mat1; 
	double*		r2 = (double *) &mat2; 
	double*		r  = (double *) &mat; 
	 
	for ( int i=0; i<3; i++ ) { 
		for ( int j=0; j<3; j++ ) { 
			for ( int k=0; k<3; k++ ) { 
				r[3*i+j] += r1[3*i+k]*r2[3*k+j]; 
			} 
		} 
	} 
	 
	return(mat); 
} 
 
/************************************************************************
@Function: matrix3_from_quaternion
@Description:
	Converts a quaternion to a 3x3 rotation matrix.
@Algorithm:
	The rotation matrix is calculated from a quaternion, (s,x,y,z), as:
    	r00 = s*s + x*x - y*y - z*z;
    	r01 = 2*x*y - 2*s*z;
    	r02 = 2*x*z + 2*s*y;
    	r10 = 2*x*y + 2*s*z;
    	r11 = s*s - x*x + y*y - z*z;
    	r12 = 2*y*z - 2*s*x;
    	r20 = 2*x*z - 2*s*y;
    	r21 = 2*y*z + 2*s*x;
    	r22 = s*s - x*x - y*y + z*z;
@Arguments:
	Quaternion q 		4-value quaternion.
@Returns:
	Matrix3 			3x3 rotation matrix.
**************************************************************************/
Matrix3 	matrix3_from_quaternion(Quaternion q)
{
	int 		i;
	Matrix3		mat; 
	double*		r = (double *) &mat; 
	 
	if ( verbose & VERB_FULL ) {
	    printf("Quaternion       :              %g %g %g %g (|q|=%g)\n",
				q.s, q.x, q.y, q.z, quaternion_size(q)); 
	} 
	 
	// Calculate the rotation matrix 
    mat.r00 = q.s*q.s + q.x*q.x - q.y*q.y - q.z*q.z; 
    mat.r01 = 2*q.x*q.y - 2*q.s*q.z; 
    mat.r02 = 2*q.x*q.z + 2*q.s*q.y; 
    mat.r10 = 2*q.x*q.y + 2*q.s*q.z; 
    mat.r11 = q.s*q.s - q.x*q.x + q.y*q.y - q.z*q.z; 
    mat.r12 = 2*q.y*q.z - 2*q.s*q.x; 
    mat.r20 = 2*q.x*q.z - 2*q.s*q.y; 
    mat.r21 = 2*q.y*q.z + 2*q.s*q.x; 
    mat.r22 = q.s*q.s - q.x*q.x - q.y*q.y + q.z*q.z; 
	 
	// Set values very close to zero to zero to avoid integer round-off errors 
	for ( i=0; i<9; i++ ) { 
		if ( fabs(r[i]) < 1e-30 ) r[i] = 0; 
		if ( r[i] > 1 ) r[i] = 1; 
		if ( r[i] < -1 ) r[i] = -1; 
	} 
	 
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG matrix3_from_quaternion:\n");
		matrix3_show(mat);
	}
	
	return(mat); 
}

/************************************************************************
@Function: matrix3_random_reslice
@Description:
	Generates a random reslicing 3x3 rotation matrix.
@Algorithm:
	The 3x3 matrix represents any one or more 90 degree rotations,
	randomly chosen.
@Arguments:
	.
@Returns:
	Matrix3 			3x3 rotation matrix.
**************************************************************************/
Matrix3		matrix3_random_reslice()
{
	double		irm = 1.0/get_rand_max();
	View		v = {NULL,0,0,0,0};
	
	int			rv = (int) (5.999*irm*random());
	int			ra = (int) (3.999*irm*random());
	switch ( rv ) {
		case 0: v = view_from_4_values(1,0,0,ra*PI/2.0); break;
		case 1: v = view_from_4_values(0,1,0,ra*PI/2.0); break;
		case 2: v = view_from_4_values(0,0,1,ra*PI/2.0); break;
		case 3: v = view_from_4_values(-1,0,0,ra*PI/2.0); break;
		case 4: v = view_from_4_values(0,-1,0,ra*PI/2.0); break;
		case 5: v = view_from_4_values(0,0,-1,ra*PI/2.0); break;
	}
	
	return(matrix3_from_view(v));
}


/*
	Quaternion section:
	------------------
*/

/************************************************************************
@Function: quaternion_from_4_values
@Description:
	Composes a quaternion from four floating point values.
@Algorithm:
	The quaternion is not normalized.
@Arguments:
	double s 			scalar value.
	double x 			first value of vector.
	double y 			second value of vector.
	double z 			third value of vector.
@Returns:
	Quaternion 			4-value quaternion.
**************************************************************************/
Quaternion	quaternion_from_4_values(double s, double x, double y, double z)
{
	Quaternion	q;
	
	q.s = s;
	q.x = x;
	q.y = y;
	q.z = z;
	
	return(q);
}

/************************************************************************
@Function: quaternion_vector
@Description:
	Returns the vector part of a quaternion.
@Algorithm:
	The vector is not normalized.
@Arguments:
	Quaternion q		4-value quaternion.
@Returns:
	double* 				quaternion vector.
**************************************************************************/
double*		quaternion_vector(Quaternion q)
{
	double*		v = (double *) balloc(3*sizeof(double));
	
	v[0] = q.x;
	v[1] = q.y;
	v[2] = q.z;
	
	return(v);
}

/************************************************************************
@Function: quaternion_show
@Description:
	Displays a quaternion.
@Algorithm:
	.
@Arguments:
	Quaternion q		4-value quaternion.
@Returns:
	int 				0.
**************************************************************************/
int 		quaternion_show(Quaternion q)
{
	printf("Quaternion = {%g,%g,%g,%g}\n",
			q.s, q.x, q.y, q.z);
	
	return(0);
}

/************************************************************************
@Function: quaternion_size
@Description:
	Quaternion size.
@Algorithm:
	The quaternion, (s, x, y, z), size is calculated as:
		|q| = sqrt(s^2 + x^2 + y^2 + z^2)
@Arguments:
	Quaternion q		4-value quaternion.
@Returns:
	double 				quaternion size.
**************************************************************************/
double 		quaternion_size(Quaternion q)
{
	int 		i;
	double		sqsum = 0, size = 0;
	double*		qv = (double *) &q;
	
	for ( i=0; i<4; i++ ) sqsum += qv[i]*qv[i];
	
	if ( sqsum ) size = sqrt(sqsum);
	
	return(size);
}

/************************************************************************
@Function: quaternion_normalize
@Description:
	Quaternion normalized.
@Algorithm:
	The quaternion size is set to 1:
		|q| = 1
@Arguments:
	Quaternion* q		4-value quaternion.
@Returns:
	double		 		quaternion size.
**************************************************************************/
double 		quaternion_normalize(Quaternion* q)
{
	double		size = quaternion_size(*q);
	
	if ( size < 1e-37 ) {
		q->s = q->x = q->y = 0;
		q->z = 1;
		size = 1;
	} else {
		q->s /= size;
		q->x /= size;
		q->y /= size;
		q->z /= size;
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_normalize: %g %g %g %g\n", 
				q->s, q->x, q->y, q->z);
 
	return(size);
}

/************************************************************************
@Function: quaternion_conjugate
@Description:
	Quaternion conjugate.
@Algorithm:
	If the input quaternion is represented as a scalar and vector:
		q = (s, v)
	The quaternion conjugate is given by:
		q* = (s, -v)
@Arguments:
	Quaternion q		4-value quaternion.
@Returns:
	Quaternion	 		4-value quaternion.
**************************************************************************/
Quaternion	quaternion_conjugate(Quaternion q)
{
	Quaternion	qnew;
	
	qnew.s =  q.s;
	qnew.x = -q.x;
	qnew.y = -q.y;
	qnew.z = -q.z;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_conjugate: %g %g %g %g\n", 
				qnew.s, qnew.x, qnew.y, qnew.z);
 
	return(qnew);
}

/************************************************************************
@Function: quaternion_inverse
@Description:
	Quaternion inverse.
@Algorithm:
	The quaternion inverse is given by:
		q-1 = q* / (|q|^2)
	where q* is the quaternion conjugate and |q| is the quaternion size.
@Arguments:
	Quaternion q		4-value quaternion.
@Returns:
	Quaternion	 		4-value quaternion.
**************************************************************************/
Quaternion	quaternion_inverse(Quaternion q)
{
	Quaternion qnew;
	
	double		size = quaternion_size(q);
	size *= size;
	
	if ( size < 1e-37 ) {
		fprintf(stderr, "Error: Size of quaternion (%g) is too small!\n", size);
		exit(-1);
	}
	
	qnew.s =  q.s/size;
	qnew.x = -q.x/size;
	qnew.y = -q.y/size;
	qnew.z = -q.z/size;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_inverse: %g %g %g %g\n", 
				qnew.s, qnew.x, qnew.y, qnew.z);
 
	return(qnew);
}

/************************************************************************
@Function: quaternion_multiply
@Description:
	Quaternion product.
@Algorithm:
	Note that the order of the input quaternions are significant.
@Arguments:
	Quaternion q1		4-value quaternion.
	Quaternion q2		4-value quaternion.
@Returns:
	Quaternion	 		4-value quaternion.
**************************************************************************/
Quaternion	quaternion_multiply(Quaternion q1, Quaternion q2)
{
	Quaternion	qnew;
	
	qnew.s = q1.s*q2.s - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;
	qnew.x = q1.s*q2.x + q1.x*q2.s + q1.y*q2.z - q1.z*q2.y;
	qnew.y = q1.s*q2.y + q1.y*q2.s + q1.z*q2.x - q1.x*q2.z;
	qnew.z = q1.s*q2.z + q1.z*q2.s + q1.x*q2.y - q1.y*q2.x;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_multiply: %g %g %g %g\n", 
				qnew.s, qnew.x, qnew.y, qnew.z);
 
	return(qnew);
}

/************************************************************************
@Function: quaternion_rotate
@Description:
	Rotates a point in 3D space using a quaternion.
@Algorithm:
	A point represented as a quaternion, p = (0,x,y,z), is rotated using 
	a quaternion:
		qnew = q * p * q-1
	where q-1 is the inverse of q.
@Arguments:
	Quaternion q		4-value quaternion representing the rotation.
	Quaternion p		4-value quaternion representing the point.
@Returns:
	Quaternion 			4-value quaternion representing the new point.
**************************************************************************/
Quaternion	quaternion_rotate(Quaternion q, Quaternion p)
{
	Quaternion	qi = quaternion_inverse(q);
	
	Quaternion	q1 = quaternion_multiply(q, p);
	
	Quaternion	qnew = quaternion_multiply(q1, qi);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_rotate: %g %g %g %g\n", 
				qnew.s, qnew.x, qnew.y, qnew.z);
 
	return(qnew);
}

/************************************************************************
@Function: quaternion_from_view_vector
@Description:
	Quaternion representing rotation from the reference to a view vector.
@Algorithm:
	Calculates a quaternion representing rotation from the reference
	vector (0,0,1) to a given view vector, v = (x,y,z).
	The vector is normalized first.
	The quaternion is given by:
		q0 = sqrt((1 + z)/2)
		if x = 0 and y = 0
			q1 = sqrt((1 - z)/2)
			q2 = 0
		else
			q1 = -y*sqrt((1 - z)/(2*(x^2+y^2)))
			q2 = x*sqrt((1 - z)/(2*(x^2+y^2)))
		q3 = 0
@Arguments:
	View view			4-value view.
@Returns:
	Quaternion			4-value quaternion.
**************************************************************************/
Quaternion	quaternion_from_view_vector(View view)
{
	Vector3 	vec = vector3_from_view(view);

	vec = vector3_normalize(vec);
	
	Quaternion	q = {1,0,0,0};
	
	double		factor, axlen;
	
	if ( vec.z < 1 ) {
		q.s =  sqrt((1.0L + vec.z)/2.0L);
		axlen = 2*(vec.x*vec.x + vec.y*vec.y);
		if ( axlen > 1e-37 ) {
			factor = sqrt((1.0L - vec.z)/axlen);
			q.x = -vec.y*factor;
			q.y =  vec.x*factor;
		} else {
			q.x =  sqrt((1.0L - vec.z)/2.0L);
		}
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_from_view_vector: %g %g %g %g\n", 
				q.s, q.x, q.y, q.z);
 
	return(q);
}

/************************************************************************
@Function: quaternion_from_two_vectors
@Description:
	Calculates a quaternion representing rotating from one vector to another.
@Algorithm:
	The vectors are normalized first.
	The rotation angle is derived from the scalar product:
		angle = arccos(scalar);
	The vector product gives the rotation axis.
	The axis vector and angle are converted to a quaternion.
@Arguments:
	float* from_vec 	3-value initial vector.
	float* to_vec		3-value final vector.
@Returns:
	Quaternion	 		4-value quaternion.
**************************************************************************/
Quaternion	quaternion_from_two_vectors(float* from_vec, float* to_vec) 
{
	float		angle = 0; 
	 
	// Convert the vectors to a unit vectors 
	vector_normalize(3, from_vec);
	vector_normalize(3, to_vec);
	
	// Calculate the vector product to get the axis of rotation
	float*		axis = vector_vector_product(3, from_vec, to_vec);
	double		scalar = vector_scalar_product(3, from_vec, to_vec);
	angle = acos(scalar);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_from_two_vectors: axis = {%g,%g,%g} angle = %g\n",
				axis[0], axis[1], axis[2], angle);
	 
	// Convert to quaternions 
	Quaternion	q = quaternion_from_angle_and_axis(angle, axis);
	 
	bfree(axis, 3*sizeof(float));
 
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_two_vectors: %g %g %g %g\n", 
				q.s, q.x, q.y, q.z);
 
	return(q); 
}	

/************************************************************************
@Function: quaternion_from_two_vectors3
@Description:
	Calculates a quaternion representing rotating from one vector to another.
@Algorithm:
	The vectors are normalized first.
	The rotation angle is derived from the scalar product:
		angle = arccos(scalar);
	The vector product gives the rotation axis.
	The axis vector and angle are converted to a quaternion.
@Arguments:
	Vector3 from_vec 	3-value initial vector.
	Vector3 to_vec		3-value final vector.
@Returns:
	Quaternion	 		4-value quaternion.
**************************************************************************/
Quaternion	quaternion_from_two_vectors3(Vector3 from_vec, Vector3 to_vec) 
{
	double		angle = 0; 
	 
	from_vec = vector3_normalize(from_vec);
	to_vec = vector3_normalize(to_vec);
	
	// Calculate the vector product to get the axis of rotation
	Vector3		axis = vector3_vector_product(from_vec, to_vec);
	double		scalar = vector3_scalar_product(from_vec, to_vec);
	angle = acos(scalar);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_from_two_vectors3: axis = {%g,%g,%g} angle = %g\n",
				axis.x, axis.y, axis.z, angle);
	 
	// Convert to quaternions 
	Quaternion	q = quaternion_from_angle_and_axis3(angle, axis);
	 
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_two_vectors3: %g %g %g %g\n", 
				q.s, q.x, q.y, q.z);
 
	return(q); 
}	

/************************************************************************
@Function: quaternion_from_angle_and_axis
@Description:
	Calculates a quaternion from a given vector and angle.
@Algorithm:
	The axis is normalized first, then converted to a quaternion:
		q = {cos(ang/2),ax_x*sin(ang/2),ax_y*sin(ang/2),ax_z*sin(ang/2)
@Arguments:
	double angle 		rotation angle.
	float* axis 		3-value axis of rotation.
@Returns:
	Quaternion	 		4-value quaternion.
**************************************************************************/
Quaternion	quaternion_from_angle_and_axis(double angle, float* axis) 
{ 
	Quaternion	q; 
	double		a = angle/2.0L;
	double		sin_a = sin(a);
	 
	vector_normalize(3, axis);
	 
	q.s = cos(a); 
	q.x = axis[0]*sin_a; 
	q.y = axis[1]*sin_a; 
	q.z = axis[2]*sin_a; 
	 
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_from_angle_and_axis: %g %g %g %g\n", 
				q.s, q.x, q.y, q.z);
 
	return(q); 
}	 
 
/************************************************************************
@Function: quaternion_from_angle_and_axis3
@Description:
	Calculates a quaternion from a given vector and angle.
@Algorithm:
	The axis is normalized first, then converted to a quaternion:
		q = {cos(ang/2),ax_x*sin(ang/2),ax_y*sin(ang/2),ax_z*sin(ang/2)
@Arguments:
	double angle 		rotation angle.
	Vector3 axis 		3-value axis of rotation.
@Returns:
	Quaternion	 		4-value quaternion.
**************************************************************************/
Quaternion	quaternion_from_angle_and_axis3(double angle, Vector3 axis) 
{ 
	Quaternion	q;
	double		a = angle/2.0L;
	double		sin_a = sin(a);
	 
	axis = vector3_normalize(axis);
	 
	q.s = cos(a); 
	q.x = axis.x*sin_a; 
	q.y = axis.y*sin_a; 
	q.z = axis.z*sin_a; 
	 
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_from_angle_and_axis3: %g %g %g %g\n", 
				q.s, q.x, q.y, q.z);
 
	return(q); 
}	 
 
/************************************************************************
@Function: quaternion_from_matrix3
@Description:
	Quaternion representing rotation based on a 3x3 matrix.
@Algorithm:
	Calculates a quaternions froma 3x3 matrix.
		s = sqrt(1+mat[1,1]+mat[2,2]+mat[3,3])/2
		x = (mat[3,2] - mat[2,3])/(4*s)
		y = (mat[1,3] - mat[3,1])/(4*s)
		z = (mat[2,1] - mat[1,2])/(4*s)
	If the trace of the rotation matrix is less or equal to zero,
	a more robust calculation is done.
@Arguments:
	Matrix3 mat 		3x3 matrix.
@Returns:
	Quaternion	 		4-value quaternion.
**************************************************************************/
Quaternion	quaternion_from_matrix3(Matrix3 mat)
{
	Quaternion	q;
	double		t = 1 + mat.r00 + mat.r11 + mat.r22;
	if ( t > 1 ) {
		t = 0.5*sqrt(t);
		q.s = t;
		t = 0.25/t;
		q.x = (mat.r21 - mat.r12)*t;
		q.y = (mat.r02 - mat.r20)*t;
		q.z = (mat.r10 - mat.r01)*t;
//		puts("** 1 **");
	} else if ( mat.r00 > mat.r11 && mat.r00 > mat.r22 ) {
		t = 0.5*sqrt(1 + mat.r00 - mat.r11 - mat.r22);
		q.x = t;
		t = 0.25/t;
		q.s = -(mat.r12 - mat.r21)*t;
		q.y = (mat.r01 + mat.r10)*t;
		q.z = (mat.r02 + mat.r20)*t;
//		puts("** 2 **");
	} else if ( mat.r11 > mat.r22 ) {
		t = 0.5*sqrt(1 - mat.r00 + mat.r11 - mat.r22);
		q.y = t;
		t = 0.25/t;
		q.s = -(mat.r02 - mat.r20)*t;
		q.x = (mat.r01 + mat.r10)*t;
		q.z = (mat.r12 + mat.r21)*t;
//		puts("** 3 **");
	} else {
		t = 0.5*sqrt(1 - mat.r00 - mat.r11 + mat.r22);
		q.z = t;
		t = 0.25/t;
		q.s = -(mat.r01 - mat.r10)*t;
		q.x = (mat.r02 + mat.r20)*t;
		q.y = (mat.r12 + mat.r21)*t;
//		puts("** 4 **");
	}
	
	quaternion_normalize(&q);
	
	if ( verbose & VERB_DEBUG ) {
		matrix3_show(mat);
		printf("DEBUG quaternion_from_matrix3: %g %g %g %g\n", 
				q.s, q.x, q.y, q.z);
  }

	return(q);
}

/************************************************************************
@Function: quaternion_from_view
@Description:
	Quaternion representing rotation from the reference to a given view.
@Algorithm:
	Calculates two quaternions, the first representing rotation from 
	the reference vector (0,0,1) to the given view vector, folowed by a
	rotation around the view vector by the view angle.
@Arguments:
	View view 			4-value view: 3-value vector and rotation angle.
@Returns:
	Quaternion	 		4-value quaternion.
**************************************************************************/
Quaternion	quaternion_from_view(View view)
{
	Quaternion	qv = quaternion_from_view_vector(view);
	
	Quaternion	qa;
	double		sin_a = sin(view.a/2.0L);
	
	qa.s = cos(view.a/2.0L);
	qa.x = view.x*sin_a;
	qa.y = view.y*sin_a;
	qa.z = view.z*sin_a;
	
	Quaternion	q = quaternion_multiply(qa, qv);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG quaternion_from_view: %g %g %g %g\n", 
				q.s, q.x, q.y, q.z);
 
	return(q);
}

/************************************************************************
@Function: view_from_quaternion
@Description:
	View vector and angle calculated from a quaternion.
@Algorithm:
	The reference vector (0,0,1) is first rotated using the given
	quaternion to obtain the view vector.
	A new quaternion is calculated from the view vector.
	The difference between the input quaternion and the new quaternion,
	is the rotation angle around the view vector.
	The rotation angle is calculated from the quaternion product between
	the inverse of the input quaternion and the new quaternion.
@Arguments:
	Quaternion q		4-value quaternion.
@Returns:
	View				4-value view: 3-value vector and rotation angle.
**************************************************************************/
View		view_from_quaternion(Quaternion q)
{
	Quaternion	pref = {0,0,0,1};
	
	Quaternion	qp = quaternion_rotate(q, pref);	// The view quaternion
	
	View		view = {NULL, qp.x, qp.y, qp.z, 0};
	view_normalize(&view);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG view_from_quaternion: vector quaternion:\n");
	Quaternion	qr = quaternion_from_view_vector(view);
	
	Quaternion	qi = quaternion_inverse(q);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG view_from_quaternion: rotation quaternion:\n");
	Quaternion	qa = quaternion_multiply(qi, qr);	// The rotation quaternion
	
	if ( qa.s < -1 ) qa.s = -1;
	if ( qa.s > 1 ) qa.s = 1;
	view.a = angle_set_negPI_to_PI(2*acos(qa.s));
	if ( qa.z > 0 ) view.a = -view.a;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG view_from_quaternion: %g %g %g %g\n", 
				view.x, view.y, view.z, view.a);
 
	return(view);
}

/************************************************************************
@Function: view_from_matrix3
@Description:
	View vector and angle calculated from a rotation matrix.
@Algorithm:
	The view vector is given by the last column in the rotation matrix.
	The view angle is given by:
		if view_x = view_y = 0 and view_z = +-1
			angle = atan(mat[[2,1]/mat[1,1])
		else
			angle = atan(mat[2,3]/mat[1,3]) + atan(-mat[3,2]/mat[3,1])
@Arguments:
	Matrix3 mat			9-value rotation matrix.
@Returns:
	View				4-value view: 3-value vector and rotation angle.
**************************************************************************/
View		view_from_matrix3(Matrix3 mat)
{
	View		view = {NULL, mat.r02, mat.r12, mat.r22, 0};
	
	view.a = view.z*atan2(mat.r10, mat.r00);
	
	if ( view.x || view.y )
		view.a = atan2(mat.r12, mat.r02) + atan2(mat.r21, -mat.r20);
	
	view.a = angle_set_negPI_to_PI(view.a);
	
	view_normalize(&view);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG view_from_matrix3: %g %g %g %g\n", 
				view.x, view.y, view.z, view.a);
 
	return(view);
}

/************************************************************************
@Function: view_from_angle_and_axis3
@Description:
	Calculates a view from a given axis and angle.
@Algorithm:
	The axis vector is normalized first, then converted to a quaternion 
	representing the rotation around the axis by the given angle.
	The quaternion is then converted to a view.
@Arguments:
	double angle 		rotation angle.
	Vector3 axis 		3-value axis of rotation.
@Returns:
	View				the view.
**************************************************************************/
View		view_from_angle_and_axis3(double angle, Vector3 axis) 
{ 
	Quaternion	q = quaternion_from_angle_and_axis3(angle, axis);
	
	// Convert to a view
	View		view = view_from_quaternion(q);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG view_from_angle_and_axis3: %g %g %g %g\n", 
				view.x, view.y, view.z, view.a);
 
	return(view); 
}	 
 
/************************************************************************
@Function: rotation_from_quaternion
@Description:
	Rotation axis and angle calculated from a quaternion.
@Algorithm:
	The rotation axis is given by the quaternion vector.
	The rotation angle is calculated from the quaternion scalar:
		angle = 2*acos(s).
@Arguments:
	Quaternion q		4-value quaternion.
@Returns:
	Rotation			7-value rotation: origin, axis and angle.
**************************************************************************/
Rotation	rotation_from_quaternion(Quaternion q)
{
	Rotation		rot = {{0,0,0},{0,0,0},0};
	
	rot.axis.x = q.x;
	rot.axis.y = q.y;
	rot.axis.z = q.z;
	rot.angle = 2*acos(q.s);
	rot.axis = vector3_normalize(rot.axis);
	
	return(rot);
}

/************************************************************************
@Function: rotation_from_matrix3
@Description:
	Rotation axis and angle calculated from a rotation matrix.
@Algorithm:
	The matrix is first converted toa quaternion.
	The rotation axis is given by the quaternion vector.
	The rotation angle is calculated from the quaternion scalar:
		angle = 2*acos(s).
@Arguments:
	Matrix3 mat			9-value rotation matrix.
@Returns:
	Rotation			7-value rotation: origin, axis and angle.
**************************************************************************/
Rotation	rotation_from_matrix3(Matrix3 mat)
{
	return(rotation_from_quaternion(quaternion_from_matrix3(mat)));
}

/*
	Euler angle conversions:
	-----------------------
*/

/************************************************************************
@Function: euler_from_three_angles
@Description:
	Converts three angles to an Euler structure.
@Algorithm:
	The angles are set to the following ranges:
		psi:	-PI to PI
		theta:  0 to PI
		phi:	-PI to PI
@Arguments:
	double psi		rotation around z ( in-plane).
	double theta		rotation around y (from z-axis).
	double phi		rotation around z (from x-axis).
@Returns:
	Euler 			3 Euler angles, psi, theta, phi.
**************************************************************************/
Euler		euler_from_three_angles(double psi, double theta, double phi)
{
	Euler		euler = {psi, theta, phi};
	
	euler.psi = angle_set_negPI_to_PI(psi); 
	euler.theta = angle_set_negPI_to_PI(theta); 
	euler.phi = angle_set_negPI_to_PI(phi); 
	
	return( euler );
}

/************************************************************************
@Function: euler_from_vector_and_angle
@Description:
	Calculates the Euler angles from a vector and a rotation around the vector.
@Algorithm:
	phi = arctan( vector_y / vector_x )
		If arccos(vector_z) is zero, the vector falls on the z-axis and phi
		does not have a meaning - it is then set to zero.
	theta = arccos( vector_z )
	psi = angle - phi
@Arguments:
	float* vector	3-value vector.
	double angle 	angle of rotation around the vector.
@Returns:
	Euler 			3 Euler angles, psi, theta, phi.
**************************************************************************/
Euler		euler_from_vector_and_angle(float* vector, double angle)
{
	Euler		euler;
	
	if ( acos(vector[2]) >  1e-30 ) euler.phi = atan2(vector[1],vector[0]);
	euler.theta = acos(vector[2]);
	euler.psi = angle - euler.phi;
	if ( euler.psi > PI )   euler.psi -= 2*PI;
	if ( euler.psi <= -PI ) euler.psi += 2*PI;
	
	return( euler );
}

/************************************************************************
@Function: euler_from_vector3_and_angle
@Description:
	Calculates the Euler angles from a vector and a rotation around the vector.
@Algorithm:
	phi = arctan( vector_y / vector_x )
		If arccos(vector_z) is zero, the vector falls on the z-axis and phi
		does not have a meaning - it is then set to zero.
	theta = arccos( vector_z )
	psi = angle - phi
@Arguments:
	Vector3 vector	3-value vector.
	double angle 	angle of rotation around the vector.
@Returns:
	Euler 			3 Euler angles, psi, theta, phi.
**************************************************************************/
Euler		euler_from_vector3_and_angle(Vector3 vector, double angle)
{
	Euler		euler;
	
	if ( acos(vector.z) >  1e-30 ) euler.phi = atan2(vector.y,vector.x);
	euler.theta = acos(vector.z);
	euler.psi = angle - euler.phi;
	if ( euler.psi > PI )   euler.psi -= 2*PI;
	if ( euler.psi <= -PI ) euler.psi += 2*PI;
	
	return( euler );
}

/************************************************************************
@Function: euler_from_view
@Description:
	Calculates the Euler angles from a view vector and a rotation angle.
@Algorithm:
	phi = arctan( vector_y / vector_x )
		If arccos(vector_z) is zero, the vector falls on the z-axis and phi
		does not have a meaning - it is then set to zero.
	theta = arccos( vector_z )
	psi = angle - phi
@Arguments:
	View view		4-value view: 3-value vector and angle.
@Returns:
	Euler 			3 Euler angles, psi, theta, phi.
**************************************************************************/
Euler		euler_from_view(View view)
{
	Euler		euler = {0,0,0};
	
	if ( acos(view.z) >  1e-30 ) euler.phi = atan2(view.y, view.x);
	euler.theta = acos(view.z);
	euler.psi = angle_set_negPI_to_PI(view.a - euler.phi);
	
	return( euler );
}

/************************************************************************
@Function: vector3_from_euler
@Description:
	Calculates the view vector from two Euler angles.
@Algorithm:
	A unit vector is calculated from the three Euler angles as:
		vector_x = cos(phi) * sin(theta)
		vector_y = sin(phi) * sin(theta)
		vector_z = cos(theta)
@Arguments:
	double theta		Euler angle theta.
	double phi		Euler angle phi.
@Returns:
	View 			3-value view vector.
**************************************************************************/
Vector3		vector3_from_euler(double theta, double phi)
{
	Vector3 	v;
	
	v.x = cos(phi)*sin(theta);
	v.y = sin(phi)*sin(theta);
	v.z = cos(theta);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG vector3_from_euler: %g %g %g\n", 
				v.x, v.y, v.z);
 
	return(v);
}

/************************************************************************
@Function: view_from_euler
@Description:
	Calculates the view from Euler angles.
@Algorithm:
	A unit vector is calculated from the three Euler angles as:
		vector_x = cos(phi) * sin(theta)
		vector_y = sin(phi) * sin(theta)
		vector_z = cos(theta)
	To complete the orientation description, the angle of rotation around
	this vector is calculated as:
		angle = psi + phi
@Arguments:
	double psi		Euler angle psi.
	double theta		Euler angle theta.
	double phi		Euler angle phi.
@Returns:
	View 			4-value view: 3-value vector and rotation angle.
**************************************************************************/
View		view_from_euler(double psi, double theta, double phi)
{
	View		view = {NULL,0,0,0,0};
	
	view.x = cos(phi)*sin(theta);
	view.y = sin(phi)*sin(theta);
	view.z = cos(theta);
	view.a = angle_set_negPI_to_PI(psi + phi);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG view_from_euler: %g %g %g %g\n", 
				view.x, view.y, view.z, view.a);
 
	return(view);
}

/************************************************************************
@Function: euler_from_rotation_matrix
@Author: David Belnap
@Description:
	Computes Euler angles (phi, theta, psi) from a rotation matrix.
@Algorithm:
	The Euler angles are calculated as:
		theta = arccos(mat[3,3])
		phi   = arctan(mat[2,3]/mat[1,3])
		psi   = arctan(mat[3,2]/mat[3,1])
	If theta = 0 or PI, then:
		phi   = 0
		psi   = arctan(mat[2,1]/mat[1,1])
	Based on Penczek et al. (1994) Ultramicroscopy 53, 251-270.
@Arguments:
	double* rot_matrix	3x3 rotation matrix
@Returns:
	Euler           	3 Euler angles as an array: phi, theta, psi.
**************************************************************************/
double*		euler_from_rotation_matrix(float* rot_matrix)
{
	int      i;
	double*  angles = (double *) balloc(3*sizeof(double));

	angles[1] = acos(rot_matrix[8]);						// theta

	if ( (angles[1] == 0) || (angles[1] == PI) )  { 		// phi is undefined when theta=0,Pi
		angles[0] = 0;										// but is set to 0
		angles[2] = atan2(rot_matrix[3],rot_matrix[0]); 	// psi
	} else {
		angles[0] = atan2(rot_matrix[5],rot_matrix[2]); 	// phi
		angles[2] = atan2(rot_matrix[7],-rot_matrix[6]);	// psi
	}

	if ( verbose & VERB_FULL )  {
		printf("Euler angles:                  ");
		for (i=0; i<3; i++)  printf(" %g",angles[i]*180/PI);
		printf("\n");
	}

	return(angles);
}

/************************************************************************
@Function: euler_from_matrix3
@Description:
	Computes Euler angles (phi, theta, psi) from a 3x3 rotation matrix.
@Algorithm:
	The Euler angles are calculated as:
		theta = arccos(mat[3,3])
		phi   = arctan(mat[2,3]/mat[1,3])
		psi   = arctan(mat[3,2]/mat[3,1])
	If theta = 0 or PI, then:
		phi   = 0
		psi   = arctan(mat[2,1]/mat[1,1])
	Based on Penczek et al. (1994) Ultramicroscopy 53, 251-270.
@Arguments:
	Matrix3 mat 		3x3 rotation matrix.
@Returns:
	Euler  				3 Euler angles: psi, theta and phi.
**************************************************************************/
Euler		euler_from_matrix3(Matrix3 mat)
{
	Euler		euler;

	euler.theta = acos(mat.r22);						// theta

	if ( ( euler.theta < 1e-37 ) || ( fabs(euler.theta-PI) < 1e-37 ) )  {
		euler.phi = 0;									// phi undefined, set to 0
		euler.psi = atan2(mat.r10, mat.r00); 			// psi
	} else {
		euler.phi = atan2(mat.r12, mat.r02); 			// phi
		euler.psi = atan2(mat.r21, -mat.r20);			// psi
	}

	if ( verbose & VERB_FULL )
		printf("Euler angles:                  %g %g %g\n",
				euler.psi, euler.theta, euler.phi);

	return(euler);
}


