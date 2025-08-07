/*
	img_interpolate.c
	Library routines for interpolating images
	Author: Bernard Heymann
	Created: 19990904	Modified: 20041209
*/

#include "img_datatypes.h"
#include "img_interpolate.h"
#include "img_extract.h"
#include "utilities.h"
#include "matrix.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

/************************************************************************
@Function: img_transform
@Description:
	Transforms an image by translation, rotation, scaling and skewing.
@Algorithm:
	A number of transformation types are combined in this function for
	efficiency. The basic operation is the conversion of the image data 
	from the original data block to a new data block, applying a 
	transformation matrix and translation vector with a specified origin 
	followed by an optional translation:
		y = R*x + t
	where	R is the rotation/skewing/scaling matrix.
			t is the translation vector.
	Note: The image is converted to floating point to improve interpolation.
@Arguments:
	Bimage* p			image (converted to floating point).
	VectorInt3 newsize	3-value new image size.
	Vector3 scale		3-value scale factor vector to apply.
	Vector3 origin		3-value origin for rotation and skewing.
	Vector3 translate	3-value vector for translation after transformation.
	Matrix3 mat			3x3 rotation or skewing matrix.
	float fill			value to fill in empty regions.
@Returns:
	int 				0.
**************************************************************************/
int 		img_transform(Bimage* p, VectorInt3 newsize, Vector3 scale, 
				Vector3 origin, Vector3 translate, Matrix3 mat, float fill)
{
	if ( p->dataflag < 1 ) return(-1);
	
	if ( p->datatype > Float ) {
		error_show("img_transform", __FILE__, __LINE__);
		printf("Error: Interpolation not done on complex data types\n");
		return(-1);
	}
	
	
	if ( scale.x < 1e-30 || p->x < 2 ) scale.x = 1;
	if ( scale.y < 1e-30 || p->y < 2 ) scale.y = 1;
	if ( scale.z < 1e-30 || p->z < 2 ) scale.z = 1;
	
	if ( newsize.x < 1 ) newsize.x = (int) (scale.x*p->x);
	if ( newsize.y < 1 ) newsize.y = (int) (scale.y*p->y);
	if ( newsize.z < 1 ) newsize.z = (int) (scale.z*p->z);
	
	unsigned long   i, j, n, x, y, z, nx, ny, nz;
	long			xo, yo, zo, xx, yy, zz, yw, zw;
	float			oldx, oldy, oldz, newx, newy, newz, xf[2], yf[2], zf[2];
	Vector3			oldorigin, neworigin;
	
	img_to_float(p);
	
	nx = ny = nz = 2;		// Number of voxels to interpolate over
	if ( nx > p->x ) nx = 1;
	if ( ny > p->y ) ny = 1;
	if ( nz > p->z ) nz = 1;
	
	oldorigin.x = neworigin.x = origin.x;
	oldorigin.y = neworigin.y = origin.y;
	oldorigin.z = neworigin.z = origin.z;

	neworigin.x += translate.x;
	neworigin.y += translate.y;
	neworigin.z += translate.z;
	
	mat = matrix3_divide(mat, scale);
	
	if ( verbose & VERB_PROCESS ) {
		printf("\nGeometric transformation:\n");
		printf("Transformation origin:          %g %g %g\n",
				oldorigin.x, oldorigin.y, oldorigin.z);
		printf("New origin:                     %g %g %g\n",
				neworigin.x, neworigin.y, neworigin.z);
		printf("Scales:                         %g %g %g\n",
				scale.x, scale.y, scale.z);
		printf("New size:                       %d %d %d voxels\n",
				newsize.x, newsize.y, newsize.z);
		printf("                                %g %g %g A\n",
				newsize.x*p->ux/scale.x, newsize.y*p->uy/scale.y, 
				newsize.z*p->uz/scale.z);
		matrix3_show(mat);
		printf("Fill value:                     %g\n", fill);
		printf("\n");
	} else if ( verbose & VERB_LABEL )
		printf("\nGeometric transformation\n\n");

	unsigned long   olddatasize = (unsigned long) p->n*p->x*p->y*p->z*sizeof(float);
	unsigned long   newdatasize = (unsigned long) p->n*newsize.x*newsize.y*newsize.z*sizeof(float);
	float*			data = (float *) p->data;
	float*			newdata = (float *) balloc(newdatasize);
	
	// Note: the matrix is used in a back-calculation of old coordinates corresponding to new
		for ( n=0; n<p->n; n++ ) {
			for ( z=0; z<newsize.z; z++ ) {
				newz = z - neworigin.z;
				for ( y=0; y<newsize.y; y++ ) {
					newy = y - neworigin.y;
					for ( x=0; x<newsize.x; x++ ) {
						newx = x - neworigin.x;
						i = ((n*newsize.z + z)*newsize.y + y)*newsize.x + x;
//						oldz =  newx*mat.r20 + newy*mat.r21 + newz*mat.r22 + oldorigin.z;
						oldz =  newx*mat.r02 + newy*mat.r12 + newz*mat.r22 + oldorigin.z;
						zo = (long) oldz;
						if ( ( zo >= 0 ) && ( zo <= p->z - nz ) ) {
//							oldy =  newx*mat.r10 + newy*mat.r11 + newz*mat.r12 + oldorigin.y;
							oldy =  newx*mat.r01 + newy*mat.r11 + newz*mat.r21 + oldorigin.y;
							yo = (long) oldy;
							if ( ( yo >= 0 ) && ( yo <= p->y - ny ) ) {
//								oldx =  newx*mat.r00 + newy*mat.r01 + newz*mat.r02 + oldorigin.x;
								oldx =  newx*mat.r00 + newy*mat.r10 + newz*mat.r20 + oldorigin.x;
								xo = (long) oldx;
								if ( ( xo >= 0 ) && ( xo <= p->x - nx ) ) {
									zf[1] = oldz - zo;
									zf[0] = 1.0 - zf[1];
									yf[1] = oldy - yo;
									yf[0] = 1.0 - yf[1];
									xf[1] = oldx - xo;
									xf[0] = 1.0 - xf[1];
									for ( zz=0; zz<nz; zz++ ) {
										zw = (n*p->z + zo + zz)*p->y;
										for ( yy=0; yy<ny; yy++ ) {
											yw = (zw + yo + yy)*p->x;
											for ( xx=0; xx<nx; xx++ ) {
												j = yw + xo + xx;
												newdata[i] += data[j]*xf[xx]*yf[yy]*zf[zz];
											}
										}
									}
								} else newdata[i] += fill;
							} else newdata[i] += fill;
						} else newdata[i] += fill;
					}
				}
			}
			p->image[n].ox = neworigin.x;
			p->image[n].oy = neworigin.y;
			p->image[n].oz = neworigin.z;
		}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_transform: Interpolation done\n");
	
	p->data = (char *) newdata;
	bfree(data, olddatasize);
	
	p->x = p->px = newsize.x;
	p->y = p->py = newsize.y;
	p->z = p->pz = newsize.z;
	p->ux /= scale.x;
	p->uy /= scale.y;
	p->uz /= scale.z;
	
	img_stats(p);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG img_transform: all done\n");
	
	return(0);
}

/************************************************************************
@Function: img_rotate_to_view
@Description:
	Rotates an image to a specified view.
@Algorithm:
	An image is first copied and then rotated according to the input 
	view vector and angle. The input image origin is the rotation center.
	There is no translation or rescaling.  Only in-plane rotations are
	allowed for 2D images.
@Arguments:
	Bimage* p				the input image.
	View* view				View vector and angle.
@Returns:
	Bimage*					The rotated image.
**************************************************************************/
Bimage*		img_rotate_to_view(Bimage* p, View view)
{
	if ( p->n > 1 ) {
		fprintf(stderr, "Rotated views for multi-image files are not supported.\n\n");
		exit(-1);
	}
	if ( p->z == 1 && ((view.x != 0) || (view.y != 0) || (view.z != 1)) )  {
		fprintf(stderr, "View vector {x, y, z} = {%g, %g, %g}\n", view.x, view.y, view.z);
		fprintf(stderr, "Error: Only in-plane rotations (e.g. cyclic symmetry) can be applied to 2D images!\n");
		exit(-1);
	}
	
	img_to_float(p);
	
	Vector3      translate = {0,0,0};
	VectorInt3   size      = {(int) p->x, (int) p->y, (int) p->z};
	Vector3      scale     = {1,1,1};
	Vector3      origin    = {p->image[0].ox, p->image[0].oy, p->image[0].oz};
	
	if ( verbose & VERB_PROCESS )
		printf("Rotating image to view:         %g %g %g %g\n", 
			view.x, view.y, view.z, view.a*180/PI);

	Bimage*		pcopy = copy_img(p);

	Matrix3		mat = matrix3_from_view(view);
//	matrix3_transpose(mat);

	img_transform(pcopy, size, scale, origin, translate, mat, p->image[0].background);
	
	pcopy->image[0].vx = view.x;
	pcopy->image[0].vy = view.y;
	pcopy->image[0].vz = view.z;
	pcopy->image[0].angle = view.a;

	img_stats(pcopy);
	
	return(pcopy);
}

/************************************************************************
@Function: img_rotate_using_euler_angles
@Description:
	Rotates an image to a specified view.
@Algorithm:
	An image is first copied and then rotated according to the input 
	euler angles.  The input image origin is the rotation center.
	There is no translation or rescaling.  Only in-plane rotations are
	allowed for 2D images.
@Arguments:
	Bimage* p				the input image.
	Euler euler				Euler angles.
@Returns:
	int						0.
**************************************************************************/
int			img_rotate_using_euler_angles(Bimage* p, Euler euler)
{
	if ( p->n > 1 ) {
		error_show("img_rotate_using_euler_angles", __FILE__, __LINE__);
		fprintf(stderr, "Rotated views for multi-image files are not supported.\n\n");
		exit(-1);
	}
	if ( p->z == 1 && ((euler.theta != 0) || (euler.phi != 0)) )  {
		error_show("img_rotate_using_euler_angles", __FILE__, __LINE__);
		fprintf(stderr, "Euler angles {psi,theta,phi} = {%g, %g, %g}\n", 
			euler.psi*180/PI, euler.theta*180/PI, euler.phi*180/PI);
		fprintf(stderr, "Error: Only in-plane rotations (e.g. cyclic symmetry) can be applied to 2D images!\n");
		exit(-1);
	}
	
	img_to_float(p);
	
	Vector3      translate = {0,0,0};
	VectorInt3   size      = {(int) p->x, (int) p->y, (int) p->z};
	Vector3      scale     = {1,1,1};
	Vector3      origin    = {p->image[0].ox, p->image[0].oy, p->image[0].oz};
	
	if ( verbose & VERB_PROCESS )
		printf("Rotating image to Euler angles: %g %g %g\n", 
			euler.psi*180/PI, euler.theta*180/PI, euler.phi*180/PI);

	Matrix3		mat = matrix3_from_euler(euler);
//	matrix3_transpose(mat);

	img_transform(p, size, scale, origin, translate, mat, p->image[0].background);
	
	View		view = view_from_euler(euler.psi, euler.theta, euler.phi);
	p->image[0].vx = view.x;
	p->image[0].vy = view.y;
	p->image[0].vz = view.z;
	p->image[0].angle = view.a;

	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_shift
@Description:
	Shifts each image as defined in individual shift vectors.
@Algorithm:
	Each image in a Bimage structure is shifted by an amount specified
	in an array of shift vectors.
@Arguments:
	Bimage* p			image to be shifted (modified to floating point).
	Vector3* shift		array of 3-value shift vectors, one for each image.
@Returns:
	int                 number of images.
**************************************************************************/
int			img_shift(Bimage* p, Vector3* shift)
{
	Vector3			origin = {0,0,0};
	Vector3			scale = {1,1,1};
	VectorInt3		size = {(int) p->x, (int) p->y, (int) p->z};
	Matrix3			mat = {1,0,0,0,1,0,0,0,1};
	Bimage*			p1;
	
	img_to_float(p);
	unsigned long   n, datasize = (unsigned long) p->x*p->y*p->z*sizeof(float);
	
	for ( n=0; n<p->n; n++ ) {
		p1 = img_extract(p, (int) n, origin, size);
		img_transform(p1, size, scale, origin, shift[n], mat, p->image[n].background);
		memcpy(p->data + datasize*n, p1->data, datasize);
		kill_img(p1);
	}
	
	return((int) n);
}

/************************************************************************
@Function: img_unique_shift_global_rotate
@Author:  Bernard Heymann and David Belnap
@Description:
	Shifts and rotates each image as defined in individual shift vectors.
@Algorithm:
	Each image in a Bimage structure is shifted by an unique amount but 
	rotated by the same angle, no scaling or resizing is done.
	Intended for use in merging single particle images from a focal, or
	other, series of particles from micrographs.
@Arguments:
	Bimage* p			image(s) to be rotated and shifted (converted to floating point).
	Vector3* origin		array of 3-value origin vectors, one for each image.
	Vector3* shift		array of 3-value shift vectors, one for each image.
	float angle			global rotation angle to apply to all images.
@Returns:
	int                 number of images.
**************************************************************************/
int			img_unique_shift_global_rotate(Bimage* p, Vector3* origin, Vector3* shift, float angle)
{
	Vector3			xorigin = {0,0,0};
	Vector3			scale = {1,1,1}, axis = {0,0,1};
	VectorInt3		size = {(int) p->x, (int) p->y, (int) p->z};
	Matrix3			mat = matrix3_from_angle_and_axis3(angle, axis);
	Bimage*			p1;
	
	img_to_float(p);
	unsigned long   n, datasize = (unsigned long) p->x*p->y*p->z*sizeof(float);
	
	for ( n=0; n<p->n; n++ ) {
		p1 = img_extract(p, (int) n, xorigin, size);
		img_transform(p1, size, scale, origin[n], shift[n], mat, p->image[n].background);
		memcpy(p->data + datasize*n, p1->data, datasize);
		kill_img(p1);
	}
	
	return((int) n);
}


/************************************************************************
@Function: img_integer_interpolation
@Description:
	Interpolates by an integer scale with a density-preserving 
	overlapping kernel.
@Algorithm:
	An image is interpolated by integer scaling (i.e., 2, 3, 4-fold or 
	more) sometimes referred to as a form of binning.  A kernel is used
	such that it overlaps with its neighbouring positions.  Voxels where
	neighbouring kernel positions overlap contribute to 2, 4 or 8 new
	voxels based on the number of overlapping kernel positions.  Only
	the central voxel is unique to a kernel position.  The kernel is 
	calculated as:
		w(i,j,k) = (1/s^2n)*(s-|s-1-i|)*(s-|s-1-j|)*(s-|s-1-k|)
	where	s is the integer interpolation factor.
			n is the number of dimensions (1D, 2D or 3D).
@Arguments:
	Bimage* p			image.
	int integer_factor	integer interpolation factor.
@Returns:
	int 				0.
**************************************************************************/
int 		img_integer_interpolation(Bimage* p, int integer_factor)
{
	if ( integer_factor < 2 ) return(0);
	
	if ( p->colormodel > Gray ) {
		fprintf(stderr, "Error: Interpolation of colour data sets not supported!\n\n");
		return(-1);
	}
		
	// Calculate the new size and kernel weights
	unsigned long   i, j;
	long			n, x, y, z, xx, yy, zz;
	VectorInt3		newsize = {1,1,1};
	VectorInt3		oldnomori = {(int) (p->x/2),(int) (p->y/2),(int) (p->z/2)};
	VectorInt3		newnomori = {0,0,0};
	int				dimensions = 0;
	VectorInt3		kernelsize = {1,1,1};
	
	if ( p->x > 1 ) {
		kernelsize.x = 2*integer_factor - 1;
		newsize.x = (int) ((p->x + 1)/integer_factor);
		if ( 2*(newsize.x/2) == newsize.x ) newsize.x--; // Size must be odd
		newnomori.x = newsize.x/2;
		p->ux *= integer_factor;
		dimensions++;
	}
	if ( p->y > 1 ) {
		kernelsize.y = 2*integer_factor - 1;
		newsize.y = (int) ((p->y + 1)/integer_factor);
		if ( 2*(newsize.y/2) == newsize.y ) newsize.y--; // Size must be odd
		newnomori.y = newsize.y/2;
		p->uy *= integer_factor;
		dimensions++;
	}
	if ( p->z > 1 ) {
		kernelsize.z = 2*integer_factor - 1;
		newsize.z = (int) ((p->z + 1)/integer_factor);
		if ( 2*(newsize.z/2) == newsize.z ) newsize.z--; // Size must be odd
		newnomori.z = newsize.z/2;
		p->uz *= integer_factor;
		dimensions++;
	}
	double		kernel_scale = 1.0/pow(1.0*integer_factor, 2.0*dimensions);
	float*		kernel = (float *) 
			balloc(kernelsize.x*kernelsize.y*kernelsize.z*sizeof(float));
	
	for ( zz=0; zz<kernelsize.z; zz++ ) {
		z = (int) (integer_factor - abs(integer_factor - 1 - zz));
		for ( yy=0; yy<kernelsize.y; yy++ ) {
			y = (int) (integer_factor - abs(integer_factor - 1 - yy));
			for ( xx=0; xx<kernelsize.x; xx++ ) {
				x = (int) (integer_factor - abs(integer_factor - 1 - xx));
				i = (zz*kernelsize.y + yy)*kernelsize.x + xx;
				kernel[i] = x*y*z*kernel_scale;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS ) {
		printf("Special integer interpolation:\n");
		printf("Integer factor:                 %d\n", integer_factor);
		printf("Interpolation kernel size:      %d x %d x %d\n", 
				kernelsize.x, kernelsize.y, kernelsize.z);
		printf("New image size:                 %d %d %d\n\n",
				newsize.x, newsize.y, newsize.z);
	} else if ( verbose & VERB_LABEL )
		printf("Special integer interpolation\n\n");

	if ( verbose & VERB_FULL ) {
		printf("Kernel interpolation: Kernel\n");
		for ( zz=0; zz<kernelsize.z; zz++ ) {
			for ( yy=0; yy<kernelsize.y; yy++ ) {
				for ( xx=0; xx<kernelsize.x; xx++ ) {
					i = (zz*kernelsize.y + yy)*kernelsize.x + xx;
					printf(" %g", kernel[i]);
				}
				printf("\n");
			}
		}
		printf("\n");
	}
	
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    
	// Do the kernel interpolation
	int 		oldx, oldy, oldz, xkernel, ykernel, zkernel;
	int 		zlo, zhi, ylo, yhi, xlo, xhi;
	float		sum;
	for ( n=0; n<p->n; n++ ) {
		for ( z=0; z<newsize.z; z++ ) {
			oldz = (int) ((z - newnomori.z)*integer_factor + oldnomori.z - kernelsize.z/2);
			zlo = oldz;
			if ( zlo < 0 ) zlo = 0;
			zhi = oldz + kernelsize.z;
			if ( zhi > p->z ) zhi = (int) p->z;
			for ( y=0; y<newsize.y; y++ ) {
				oldy = (int) ((y - newnomori.y)*integer_factor + oldnomori.y - kernelsize.y/2);
				ylo = oldy;
				if ( ylo < 0 ) ylo = 0;
				yhi = oldy + kernelsize.y;
				if ( yhi > p->y ) yhi = (int) p->y;
				for ( x=0; x<newsize.x; x++ ) {
					oldx = (int) ((x - newnomori.x)*integer_factor + oldnomori.x - kernelsize.x/2);
					xlo = oldx;
					if ( xlo < 0 ) xlo = 0;
					xhi = oldx + kernelsize.x;
					if ( xhi > p->x ) xhi = (int) p->x;
					sum = 0;
					for ( zz=zlo; zz<zhi; zz++ ) {
						zkernel = (int) (zz - oldz);
						for ( yy=ylo; yy<yhi; yy++ ) {
							ykernel = (int) (yy - oldy);
							for ( xx=xlo; xx<xhi; xx++ ) {
								xkernel = (int) (xx - oldx);
								i = ((n*p->z + zz)*p->y + yy)*p->x + xx;
								j = (zkernel*kernelsize.y + ykernel)*kernelsize.x + xkernel;
			    				switch ( p->datatype ) {
  									case UChar:
										sum += udata[i]*kernel[j];
			    	    				break;
	    							case SChar:
										sum += cdata[i]*kernel[j];
							    	    break;
							    	case UShort:
										sum += usdata[i]*kernel[j];
						    		    break;
							    	case Short:
										sum += sdata[i]*kernel[j];
			    						break;
									case Int:
										sum += idata[i]*kernel[j];
			    						break;
									case Float:
										sum += fdata[i]*kernel[j];
			    						break;
    								default: break;
								}
							}
						}
					}
					i = ((n*newsize.z + z)*newsize.y + y)*newsize.x + x;
	    			switch ( p->datatype ) {
  						case UChar:
							udata[i] = (unsigned char) sum;
    	    				break;
    					case SChar:
							cdata[i] = (signed char) sum;
			    			break;
				    	case UShort:
							usdata[i] = (unsigned short) sum;
    						break;
				    	case Short:
							sdata[i] = (short) sum;
				    	    break;
						case Int:
							idata[i] = (int) sum;
    						break;
						case Float:
							fdata[i] = sum;
    						break;
	    				default: break;
					}
				}
			}
		}
		p->image[n].ox = (p->image[n].ox - oldnomori.x)/integer_factor + newnomori.x;
		p->image[n].oy = (p->image[n].oy - oldnomori.y)/integer_factor + newnomori.y;
		p->image[n].oz = (p->image[n].oz - oldnomori.z)/integer_factor + newnomori.z;
		if ( verbose & VERB_STATS )
			printf("Image %ld origin:                %g %g %g\n", n+1, 
				p->image[n].ox, p->image[n].oy, p->image[n].oz);
	}

	if ( verbose & VERB_DEBUG )
		printf("Kernel interpolation: Sums done\n");
	
	// Transfer the data to a new block and free the old
	unsigned long   datasize = (unsigned long) p->n*newsize.x*newsize.y*newsize.z*p->c*gettypesize(p->datatype);
	p->data = (char *) balloc(datasize);
	memcpy(p->data, udata, datasize);
	bfree(udata, p->n*p->x*p->y*p->z*p->c*gettypesize(p->datatype));
	bfree(kernel, kernelsize.x*kernelsize.y*kernelsize.z*sizeof(float));
	
	// Set new image parameters
	p->x = p->px = newsize.x;
	p->y = p->py = newsize.y;
	p->z = p->pz = newsize.z;	
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_bin
@Description:
	Bins by an integer size.
@Algorithm:
	An image is binned by an integer size, square in 2D and cubic in 3D.
@Arguments:
	Bimage* p		image.
	VectorInt3 bin	3-value vector of integer bin factors.
@Returns:
	int 			0.
**************************************************************************/
int 		img_bin(Bimage* p, VectorInt3 bin)
{
	if ( p->dataflag < 1 ) return(-1);
	
	if ( bin.x < 1 ) bin.x = 1;
	if ( bin.y < 1 ) bin.y = 1;
	if ( bin.z < 1 ) bin.z = 1;
	if ( p->x < bin.x ) bin.x = (int) p->x;
	if ( p->y < bin.y ) bin.y = (int) p->y;
	if ( p->z < bin.z ) bin.z = (int) p->z;
	
	int				binsize = bin.x*bin.y*bin.z;
	
	if ( binsize < 2 ) return(-1);
	
	if ( p->colormodel > Gray ) {
		error_show("img_bin", __FILE__, __LINE__);
		fprintf(stderr, "Error: Interpolation of colour data sets not supported!\n\n");
		return(-1);
	}
	
	unsigned long   i, n, x, y, z, xx, yy, zz;
	VectorInt3		newsize;
	float			binsum;
	
	newsize.x = (int) ((p->x + bin.x - 1)/bin.x);
	newsize.y = (int) ((p->y + bin.y - 1)/bin.y);
	newsize.z = (int) ((p->z + bin.z - 1)/bin.z);
	
	if ( verbose & VERB_PROCESS ) {
		printf("Binning:\n");
		printf("Bin size:                       %d x %d x %d = %d\n", 
				bin.x, bin.y, bin.z, binsize);
		printf("New image size:                 %d %d %d\n\n",
				newsize.x, newsize.y, newsize.z);
	} else if ( verbose & VERB_LABEL )
		printf("Binning\n\n");

    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
	
	// Add all required voxels into a floating point bin
	for ( n=0; n<p->n; n++ ) {
		for ( z=0; z<newsize.z; z++ ) {
			for ( y=0; y<newsize.y; y++ ) {
				for ( x=0; x<newsize.x; x++ ) {
					binsum = 0;
					for ( zz=bin.z*z; zz<bin.z*(z+1) && zz<p->z; zz++ ) {
						for ( yy=bin.y*y; yy<bin.y*(y+1) && yy<p->y; yy++ ) {
							for ( xx=bin.x*x; xx<bin.x*(x+1) && xx<p->x; xx++ ) {
								i = ((n*p->z + zz)*p->y + yy)*p->x + xx;
			    				switch ( p->datatype ) {
  									case UChar:
										binsum += udata[i];
			    	    				break;
	    							case SChar:
										binsum += cdata[i];
							    	    break;
							    	case UShort:
										binsum += usdata[i];
						    		    break;
							    	case Short:
										binsum += sdata[i];
			    						break;
									case Int:
										binsum += idata[i];
			    						break;
									case Float:
										binsum += fdata[i];
			    						break;
    								default: break;
								}
							}
						}
					}
					binsum /= binsize;
					i = ((n*newsize.z + z)*newsize.y + y)*newsize.x + x;
	    			switch ( p->datatype ) {
  						case UChar:
							udata[i] = (unsigned char) binsum;
    	    				break;
    					case SChar:
							cdata[i] = (signed char) binsum;
			    			break;
				    	case UShort:
							usdata[i] = (unsigned short) binsum;
    						break;
				    	case Short:
							sdata[i] = (short) binsum;
				    	    break;
						case Int:
							idata[i] = (int) binsum;
    						break;
						case Float:
							fdata[i] = binsum;
    						break;
	    				default: break;
					}
				}
			}
		}
		p->image[n].ox /= bin.x;
		p->image[n].oy /= bin.y;
		p->image[n].oz /= bin.z;
		if ( verbose & VERB_STATS )
			printf("Image %ld origin:                %g %g %g\n", n+1, 
				p->image[n].ox, p->image[n].oy, p->image[n].oz);
	}

	if ( verbose & VERB_DEBUG )
		printf("Binning: Sums done\n");
	
	// Transfer the data to a new block and free the old
	unsigned long   datasize = (unsigned long) p->n*newsize.x*newsize.y*newsize.z*p->c*gettypesize(p->datatype);
	p->data = (char *) balloc(datasize);
	memcpy(p->data, udata, datasize);
	bfree(udata, p->n*p->x*p->y*p->z*p->c*gettypesize(p->datatype));
	
	// Set new image parameters
	p->x = p->px = newsize.x;
	p->y = p->py = newsize.y;
	p->z = p->pz = newsize.z;	
	p->ux *= bin.x;
	p->uy *= bin.y;
	p->uz *= bin.z;
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_median_bin
@Description:
	Bins by an integer size.
@Algorithm:
	An image is binned by an integer size, square in 2D and cubic in 3D.
@Arguments:
	Bimage* p		image.
	int binning 	integer bin factor.
@Returns:
	int 			0.
**************************************************************************/
int 		img_median_bin(Bimage* p, int binning)
{
	if ( binning < 2 ) return(-1);
	
	if ( p->colormodel > Gray ) {
		error_show("img_median_bin", __FILE__, __LINE__);
		fprintf(stderr, "Error: Interpolation of colour data sets not supported!\n\n");
		return(-1);
	}
	
	// Calculate the new size and set up
	unsigned long   i, j;
	long 			n, x, y, z, xx, yy, zz, dimensions = 0;
	VectorInt3		newsize = {1,1,1};
	VectorInt3 		shift = {0,0,0};
	VectorInt3 		blocksize = {1,1,1};
	if ( p->x > 1 ) {
		blocksize.x = 2*binning - 1;
		newsize.x = (int) ((p->x + 1)/binning);
		if ( 2*(newsize.x/2) == newsize.x ) newsize.x--; // Size must be odd
		shift.x = (int) (newsize.x - p->x/binning);
		p->ux *= binning;
		dimensions++;
	}
	if ( p->y > 1 ) {
		blocksize.y = 2*binning - 1;
		newsize.y = (int) ((p->y + 1)/binning);
		if ( 2*(newsize.y/2) == newsize.y ) newsize.y--; // Size must be odd
		shift.y = (int) (newsize.y - p->y/binning);
		p->uy *= binning;
		dimensions++;
	}
	if ( p->z > 1 ) {
		blocksize.z = 2*binning - 1;
		newsize.z = (int) ((p->z + 1)/binning);
		if ( 2*(newsize.z/2) == newsize.z ) newsize.z--; // Size must be odd
		shift.z = (int) (newsize.z - p->z/binning);
		p->uz *= binning;
		dimensions++;
	}
	unsigned long	blocktotalsize = blocksize.x*blocksize.y*blocksize.z;
	int*			maskblock = (int *) balloc(blocktotalsize*sizeof(int));
	float*			block = (float *) balloc(blocktotalsize*sizeof(float));
	
	if ( verbose & VERB_PROCESS ) {
		printf("Binning using a median filter:\n");
		printf("Block size:                     %d x %d x %d = %ld\n", 
				blocksize.x, blocksize.y, blocksize.z, blocktotalsize);
		printf("New image size:                 %d %d %d\n\n",
				newsize.x, newsize.y, newsize.z);
	} else if ( verbose & VERB_LABEL )
		printf("Binning using a median filter\n\n");

    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    
	// Do the median binning interpolation
	long 		oldx, oldy, oldz, xblock, yblock, zblock, nblock, imedian;
	long 		zlo, zhi, ylo, yhi, xlo, xhi;
	for ( n=0; n<p->n; n++ ) {
		for ( z=0; z<newsize.z; z++ ) {
			oldz = z*binning - shift.z; 		// Origin of the block
			zlo = oldz;
			if ( zlo < 0 ) zlo = 0;
			zhi = oldz + blocksize.z;
			if ( zhi > p->z ) zhi = p->z;
			for ( y=0; y<newsize.y; y++ ) {
				oldy = y*binning - shift.y;
				ylo = oldy;
				if ( ylo < 0 ) ylo = 0;
				yhi = oldy + blocksize.y;
				if ( yhi > p->y ) yhi = p->y;
				for ( x=0; x<newsize.x; x++ ) {
					oldx = x*binning - shift.x;
					xlo = oldx;
					if ( xlo < 0 ) xlo = 0;
					xhi = oldx + blocksize.x;
					if ( xhi > p->x ) xhi = p->x;
					nblock = 0;
					memset(maskblock, 0, blocktotalsize);
					memset(block, 0, blocktotalsize);
					for ( zz=zlo; zz<zhi; zz++ ) {
						zblock = zz - oldz;
						for ( yy=ylo; yy<yhi; yy++ ) {
							yblock = yy - oldy;
							for ( xx=xlo; xx<xhi; xx++ ) {
								xblock = xx - oldx;
								i = ((n*p->z + zz)*p->y + yy)*p->x + xx;
								j = (zblock*blocksize.y + yblock)*blocksize.x + xblock;
			    				switch ( p->datatype ) {
  									case UChar:
										block[j] = udata[i];
			    	    				break;
	    							case SChar:
										block[j] = cdata[i];
							    	    break;
							    	case UShort:
										block[j] = usdata[i];
						    		    break;
							    	case Short:
										block[j] = sdata[i];
			    						break;
									case Int:
										block[j] = idata[i];
			    						break;
									case Float:
										block[j] = fdata[i];
			    						break;
    								default: break;
								}
								maskblock[j] = 1;
								nblock++;
							}
						}
					}
					for ( i=1; i<blocktotalsize; i++ ) {
						for ( j=0; j<i; j++ ) {
							if ( maskblock[j] < maskblock[i] ) {
								swap_integers(&maskblock[i], &maskblock[j]);
								swap_floats(&block[i], &block[j]);
							} else if ( maskblock[i] && block[j] < block[i] ) {
								swap_integers(&maskblock[i], &maskblock[j]);
								swap_floats(&block[i], &block[j]);
							}
						}
					}
					imedian = nblock/2;
//					for ( i=0; i<blocktotalsize; i++ ) printf("%g\t", block[i]);
//					printf("%d\n", imedian);
					i = ((n*newsize.z + z)*newsize.y + y)*newsize.x + x;
	    			switch ( p->datatype ) {
  						case UChar:
							udata[i] = (unsigned char) block[imedian];
    	    				break;
    					case SChar:
							cdata[i] = (signed char) block[imedian];
			    			break;
				    	case UShort:
							usdata[i] = (unsigned short) block[imedian];
    						break;
				    	case Short:
							sdata[i] = (short) block[imedian];
				    	    break;
						case Int:
							idata[i] = (int) block[imedian];
    						break;
						case Float:
							fdata[i] = block[imedian];
    						break;
	    				default: break;
					}
				}
			}
		}
		p->image[n].ox /= binning;
		p->image[n].oy /= binning;
		p->image[n].oz /= binning;
		if ( verbose & VERB_STATS )
			printf("Image %ld origin:                %g %g %g\n", n+1, 
				p->image[n].ox, p->image[n].oy, p->image[n].oz);
	}

	if ( verbose & VERB_DEBUG )
		printf("Median binning interpolation done\n");
	
	// Transfer the data to a new block and free the old
	unsigned long   datasize = (unsigned long) p->n*newsize.x*newsize.y*newsize.z*p->c*gettypesize(p->datatype);
	p->data = (char *) balloc(datasize);
	memcpy(p->data, udata, datasize);
	bfree(udata, p->n*p->x*p->y*p->z*p->c*gettypesize(p->datatype));
	bfree(maskblock, blocktotalsize*sizeof(int));
	bfree(block, blocktotalsize*sizeof(float));
	
	// Set new image parameters
	p->x = p->px = newsize.x;
	p->y = p->py = newsize.y;
	p->z = p->pz = newsize.z;	
	
	img_stats(p);
	
	return(0);
}

