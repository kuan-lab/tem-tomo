/* 
	matrix.h 
	Header file for matrix functions 
	Author: Bernard Heymann 
	Created: 20000501 	Modified: 20050201
*/ 
 
// Floating point vector structures
#ifndef _Vector_
#define _Vector_
struct Vector2 {
	float x,y;		// The vector elements
} ;

struct Vector3 {
	float x,y,z;	// The vector elements
} ;

struct Vector4 {
	float x,y,z,n;	// The vector elements
} ;
#endif

// Double precision floating point vector structures
#ifndef _VectorDouble_
#define _VectorDouble_
struct VectorDouble2 {
	double x,y;		// The vector elements
} ;

struct VectorDouble3 {
	double x,y,z;	// The vector elements
} ;

struct VectorDouble4 {
	double x,y,z,n;	// The vector elements
} ;
#endif

// Integer vector structures
#ifndef _VectorInt_
#define _VectorInt_
struct VectorInt2 {
	int x,y;		// The vector elements
} ;

struct VectorInt3 {
	int x,y,z;		// The vector elements
} ;

struct VectorInt4 {
	int x,y,z,n;	// The vector elements
} ;
#endif

// Integer long vector structures
#ifndef _VectorLong_
#define _VectorLong_
struct VectorLong2 {
	long x,y;		// The vector elements
} ;

struct VectorLong3 {
	long x,y,z;		// The vector elements
} ;

struct VectorLong4 {
	long x,y,z,n;	// The vector elements
} ;
#endif

// 9-value matrix structure
#ifndef _Matrix3_
#define _Matrix3_
struct Matrix3 {
	double r00, r01, r02;	// First row
	double r10, r11, r12;	// Second row
	double r20, r21, r22;	// Third row
} ;
#endif

// View structure
#ifndef _View_
#define _View_
struct View {
	View* next;		// Pointer for linked lists
	float x,y,z;	// The view vector - line of sight
	float a;		// The view rotation angle - rotation perpendicular to view vector
} ;
#endif

// Rotation structure
#ifndef _Rotation_
#define _Rotation_
struct Rotation {
	Vector3 origin;		// Rotation origin
	Vector3 axis;		// Rotation axis
	double angle;		// Rotation angle
} ;
#endif

// Quaternion structure
#ifndef _Quaternion_
#define _Quaternion_
struct Quaternion {
	double s;		// The scalar
	double x,y,z;	// The vector
} ;
#endif

// Euler structure
#ifndef _Euler_
#define _Euler_
struct Euler {
	double psi;		// view.angle - phi
	double theta;	// arccos(view.z)
	double phi;		// arctan(view.y/view.x)
} ;
#endif

// Function prototypes 
float*		init_unit_matrix(int dim, int size);
double		matrix_determinant(int n, float* mat) ;
float*		matrix_product(int m, int n, int p, float* m1, float* m2);
float*		matrix_transpose(int m, int n, float* mat);
int 		show_matrix(int m, int n, float* mat);
int 		show_matrix(int m, int n, double* mat);
float*		matrix_vector_multiply(int n, float* mat, float* vector);
int			matrix_vector_multiply_in_place(int n, double* mat, double* vec);
double 		vector_normalize(int n, float* vector);
double 		vector_normalize(int n, double* vector);
double 		vector_scalar_product(int n, float* vector1, float* vector2);
float* 		vector_vector_product(int n, float* vector1, float* vector2);
int			show_view(View v);
int			show_views(View* v);
View		view_from_4_values(double x, double y, double z, double a);
float*		view_vector(View view);
double		view_vector_size(View view);
double		view_normalize(View* view);
View		view_negate(View view);
double		view_difference(View view1, View view2);
Vector3 	vector3_from_vectorint3(VectorInt3 intvec);
VectorInt3 	vectorint3_from_vector3(Vector3 floatvec);
Vector3 	vector3_zero();
Vector3 	vector3_from_3_values(double x, double y, double z);
Vector3 	vector3_random(double min, double max);
Vector3 	vector3_random_gaussian(double avg, double std);
Vector3 	vector3_xy_random_gaussian(double avg, double std);
Vector3 	vector3_normalize(Vector3 vector);
Vector3 	vector3_floor(Vector3 vector, int places);
Vector3 	vector3_round(Vector3 vector, int places);
Vector3 	vector3_remainder(Vector3 vector, int divisor);
Vector3 	vector3_from_view(View view);
double		vector3_length(Vector3 vector);
double		vector3_length2(Vector3 vector);
Vector3		vector3_scale(Vector3 vector, double scale);
Vector3		vector3_multiply(Vector3 vec1, Vector3 vec2);
Vector3		vector3_divide(Vector3 vec1, Vector3 vec2);
Vector3		vector3_square(Vector3 vector);
Vector3		vector3_square_root(Vector3 vector);
Vector3		vector3_negate(Vector3 vec);
Vector3		vector3_add(Vector3 vec1, Vector3 vec2);
Vector3		vector3_subtract(Vector3 vec1, Vector3 vec2);
Vector3		vector3_add_scalar(Vector3 vec, double s);
double	 	vector3_scalar_product(Vector3 v1, Vector3 v2);
Vector3 	vector3_vector_product(Vector3 v1, Vector3 v2);
double	 	vector3_angle(Vector3 v1, Vector3 v2);
Vector3		vector3_matrix3_multiply(Matrix3 mat, Vector3 vec);
Vector3 	vector3_calculate_normal(Vector3 v1, Vector3 v2, Vector3 v3);
Vector3		vector3_min(Vector3 vector1, Vector3 vector);
Vector3		vector3_max(Vector3 vector1, Vector3 vector);
Vector3		vector3_scalar_min(Vector3 vector, double scalar);
Vector3		vector3_scalar_max(Vector3 vector, double scalar);
Vector3		vector3_scalar_range(Vector3 vec, double min, double max);
Vector3		vector3_set_PBC(Vector3 coord, Vector3 box);
Vector3		vector3_difference_PBC(Vector3 v1, Vector3 v2, Vector3 box);
VectorInt3 	vectorint3_from_3_values(int x, int y, int z);
double		vectorint3_length(VectorInt3 vector);
double		vectorint3_length2(VectorInt3 vector);
VectorInt3  vectorint3_add(VectorInt3 vec1, VectorInt3 vec2);
VectorInt3	vectorint3_subtract(VectorInt3 vec1, VectorInt3 vec2);
VectorInt3  vectorint3_add_scalar(VectorInt3 vec, int s);
VectorInt3	vectorint3_scale(VectorInt3 vec1, double scale);
VectorInt3	vectorint3_square(VectorInt3 vector);
VectorInt3	vectorint3_min(VectorInt3 vector1, VectorInt3 vector2);
VectorInt3	vectorint3_max(VectorInt3 vector1, VectorInt3 vector2);
VectorInt3	vectorint3_scalar_min(VectorInt3 vector, int scalar);
VectorInt3	vectorint3_scalar_max(VectorInt3 vector, int scalar);
VectorLong3 vectorlong3_from_3_values(long x, long y, long z);
int 		matrix3_show(Matrix3 mat);
Matrix3		matrix3_init();
double		matrix3_determinant(Matrix3 mat); 
Matrix3		matrix3_transpose(Matrix3 mat);
Matrix3		matrix3_scalar_multiply(Matrix3 mat, double scalar);
Matrix3		matrix3_vector3_multiply(Matrix3 mat, Vector3 vec);
Matrix3		matrix3_divide(Matrix3 mat, Vector3 vec);
Matrix3		matrix3_from_view(View view);
Matrix3		matrix3_from_angle_and_axis(double angle, float* axis); 
Matrix3		matrix3_from_angle_and_axis3(double angle, Vector3 axis);
Matrix3		matrix3_from_tilt_and_axis(double tilt, double axis); 
Matrix3		matrix3_from_two_vectors(float* vector1, float* vector2);
Matrix3		matrix3_from_two_vectors3(Vector3 from_vec, Vector3 to_vec);
Matrix3		matrix3_from_euler(Euler euler);
Matrix3		matrix3_from_psi_theta_phi(double psi, double theta, double phi);
Matrix3		matrix3_from_theta_phi_omega(double theta, double phi, double omega); 
Matrix3		matrix3_multiply(Matrix3 mat1, Matrix3 mat2); 
Matrix3 	matrix3_from_quaternion(Quaternion q);
Matrix3		matrix3_random_reslice();
Quaternion	quaternion_from_4_values(double s, double x, double y, double z);
double*		quaternion_vector(Quaternion q);
int 		quaternion_show(Quaternion q);
double 		quaternion_size(Quaternion q);
double 		quaternion_normalize(Quaternion* q);
Quaternion	quaternion_conjugate(Quaternion q);
Quaternion	quaternion_inverse(Quaternion q);
Quaternion	quaternion_multiply(Quaternion q1, Quaternion q2);
Quaternion	quaternion_rotate(Quaternion q, Quaternion p);
Quaternion	quaternion_from_view_vector(View view);
Quaternion	quaternion_from_two_vectors(float* from_vec, float* to_vec);
Quaternion	quaternion_from_two_vectors3(Vector3 from_vec, Vector3 to_vec);
Quaternion	quaternion_from_angle_and_axis(double angle, float* axis);
Quaternion	quaternion_from_angle_and_axis3(double angle, Vector3 axis);
Quaternion	quaternion_from_matrix3(Matrix3 mat);
Quaternion	quaternion_from_view(View view);
View		view_from_quaternion(Quaternion q);
View		view_from_matrix3(Matrix3 mat);
View		view_from_angle_and_axis3(double angle, Vector3 axis);
Rotation	rotation_from_quaternion(Quaternion q);
Rotation	rotation_from_matrix3(Matrix3 mat);
Euler		euler_from_three_angles(double psi, double theta, double phi);
Euler		euler_from_vector_and_angle(float* vector, double angle);
Euler		euler_from_vector3_and_angle(Vector3 vector, double angle);
Euler		euler_from_view(View view);
Vector3		vector3_from_euler(double theta, double phi);
View		view_from_euler(double psi, double theta, double phi);
double*		euler_from_rotation_matrix(float* rot_matrix);
Euler		euler_from_matrix3(Matrix3 mat);
