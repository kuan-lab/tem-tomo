/*
	project.cc
	Functions for projections from volumes
	Author: Giovanni Cardone
	Created: 20050119 	Modified: 20050608
*/

#include "project.h"

/************************************************************************
@Function: prjct_generate_one
@Description:
	Calculates a single projection from a 3D density map.
@Algorithm:
	A single projection is calculated, with angle taken from
	an array of view vectors.
	The model is rotated and then projected.
(	The rotation angle around each view vector is taken as zero. )
@Arguments:
	Bimage* p			the 3D map.
	int iv				index of view to generate.
	View* view			array of 4-value views.
	int taper			tapering flag
	int bground_type		background evaluation method
	int exp_tag			exp relationship between intensity and density
	Buff* bf				buffer structure associate to the 3D map
@Returns:
	Proj* 				calculated projection.
**************************************************************************/
Proj* 	prjct_generate_one(Bimage* b, int iv, View* view, Shape3d* volut, int taper,
							int bground_type, int exp_tag, Buff* bf)
{
	
	Matrix3		rotat_mat;
	Vector3		translat_vec = {0.,0.,0.};   // need to become an input argument

	if ( verbose & VERB_PROCESS )
		printf("Calculating projection # %d:\n", iv+1);

	Proj* p = proj_init(1);
	proj_cpy_img_sett(p, b);
	Bimage* 	pb = init_img_with_parameters(Float, 1, b->x, b->y, 1, 1);
	pb->ux = b->ux;
	pb->uy = b->uy;
	pb->image[0].ox = 0.5 * ( b->x - 1.0 );
	pb->image[0].oy = 0.5 * ( b->y - 1.0 );
	pb->image[0].oz = 0.;
	proj_putimage(p, pb, 0, 0);
	proj_putview(p, view[iv], 0);
	

	rotat_mat = matrix3_from_view(view[iv]);
	prjct_rotate_and_project(p, b, 0, rotat_mat, translat_vec, volut, exp_tag, bf); 

	Vector3 origin = vector3_from_3_values(b->image[0].ox,b->image[0].oy,b->image[0].oz);
	proj_init_footprint(p, origin, volut, taper);

	proj_stats(p, bground_type);
	
	proj_mask(p);

	return(p);
}

/************************************************************************
@Function: prjct_rotate_and_project
@Description:
	Rotates a 3D map and projects it along the z-axis.
@Algorithm:
	The 3D map is rotated around its center as origin and the data 
	integrated along the z-direction after subtraction of the
	background (which must be calculated before this function).
	The resultant 2D image is translated if the translation vector is 
	non-zero.
	The rotation origin is obtained from the map origin.
	If a mask given, only the map contained in the mask is rotated and projected
@Arguments:
	Bimage*	p			3D image.
	Proj*	proj			2D output projection (must be allocated).
	int		ip			projection index
	Matrix3	mat			3x3 rotation or skewing matrix.
	Vector3	translat_vec 3-value vector for translation after transformation.
	Shape3d*	volut		volume selection for projection
	int 		exp_tag		exp relationship between intensity and density
	Buff* bf				buffer structure associate to the 3D map
@Returns:
	int					0.
**************************************************************************/
int		prjct_rotate_and_project(Proj* proj, Bimage* p, int ip, Matrix3 rotat_mat,
							 Vector3 translat_vec, Shape3d* volut, int exp_tag, Buff* bf)
{
	if ( !p->data ) return(-1);
	if ( !(proj->img[ip])->data ) return(-1);
	
	if (!bf && p->datatype!=UChar && p->datatype!=SChar && p->datatype!=Short && p->datatype!=UShort && p->datatype!=Int && p->datatype!=Float)
		img_to_float(p);

	img_to_float(proj->img[ip]);
	
	unsigned long i;
	long			j, xo, yo, xx, yy, yw, xincr, yincr;
	int          x, y, z;
	double 		oldx, oldy, oldz;
	double		newx, newy, xf, yf, w;
	Vector3		oldorigin, neworigin;
	
	
//	oldorigin.x = neworigin.x = p->image[0].ox;
//	oldorigin.y = neworigin.y = p->image[0].oy;
//	oldorigin.z = neworigin.z = p->image[0].oz;
	oldorigin.x = p->image[0].ox;
	oldorigin.y = p->image[0].oy;
	oldorigin.z = p->image[0].oz;

	neworigin.x = (proj->img[ip])->image[0].ox;
	neworigin.y = (proj->img[ip])->image[0].oy;
	neworigin.z = (proj->img[ip])->image[0].oz;

	neworigin.x += translat_vec.x;
	neworigin.y += translat_vec.y;
	neworigin.z += translat_vec.z;
	// The rotation matrix is transposed because the orientation parameters
	// give the orientation as rotated from the reference view
	rotat_mat = matrix3_transpose(rotat_mat);

	if ( verbose & VERB_FULL ) {
		printf("\nRotating around the map center and projecting:\n");
		printf("Rotation origin:                %g %g %g pixels\n",
				oldorigin.x, oldorigin.y, oldorigin.z);
		matrix3_show(rotat_mat);
		if ( translat_vec.x != 0 || translat_vec.y != 0 || translat_vec.z != 0 )
			printf("Translation:                    %g %g %g pixels\n",
					translat_vec.x, translat_vec.y, translat_vec.z);
		printf("\n");
	} else if ( verbose & ( VERB_LABEL | VERB_PROCESS) )
		printf("\nRotating around the map center and projecting\n\n");

    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    		fdata = (float *) p->data;
	float*		newdata = (float *) (proj->img[ip])->data;
	float*		num = (float *) balloc(p->x*p->y*sizeof(float));
	double		background = p->image[0].background;
	double		value = 0.;
	
	for ( z=0; z< (int) p->z; z++ ) {
		oldz = z - oldorigin.z;
		for ( y=0; y< (int) p->y; y++ ) {
			oldy = y - oldorigin.y;
			for ( x=0; x< (int) p->x; x++ ) {
				if ( !shape3d_point_inside(volut,x,y,z))  continue;
				oldx = x - oldorigin.x;
//				i = (z*p->y + y)*p->x + x;
				i = io_get_index(x,y,z,p,bf);
				newy = oldx*rotat_mat.r10 + oldy*rotat_mat.r11 + oldz*rotat_mat.r12 + neworigin.y;
				yo = (int) newy;
				if ( ( yo >= 0 ) && ( yo <  (int) p->y ) ) {
					newx = oldx*rotat_mat.r00 + oldy*rotat_mat.r01 + oldz*rotat_mat.r02 + neworigin.x;
					xo = (int) newx;
					if ( ( xo >= 0 ) && ( xo <  (int) p->x ) ) {
						yf = newy - yo;
						xf = newx - xo;
 						switch ( p->datatype ) {
    						case UChar:
							value = udata[i] - background;
						break;
    						case SChar:
							value = cdata[i] - background;
						break;
    						case UShort:
							value = usdata[i] - background;
						break;
    						case Short:
							value = sdata[i] - background;
						break;
    						case Int:
							value = idata[i] - background;
						break;
    						case Float:
							value = fdata[i] - background;
						break;
   						default: break;
					    }
						
						yincr = yo== (int) p->y-1?1:2;
						xincr = xo== (int) p->x-1?1:2;
						for ( yy=yo; yy<yo+yincr; yy++ ) {
							yw = yy*p->x;
							yf = 1 - yf;
							for ( xx=xo; xx<xo+xincr; xx++ ) {
								j = yw + xx;
								xf = 1 - xf;
								w = xf*yf;
								newdata[j] += value*w;
								num[j] += w;
							}
						}
					}
				}
			}
		}
	}
	
	// normalize with respect to interpolation coefficients
	for ( y=0; y< (int) p->y; y++ ) {
		for ( x=0; x< (int) p->x; x++ ) {
			j = y*p->x + x;
			if ( num[j] ) {
				newdata[j] = newdata[j] / num[j] + background;
				if (exp_tag) newdata[j] = exp((newdata[j]-p->min)/(p->max-p->min))*(p->max-p->min);
			}
		}
	}
	
	bfree(num, p->x*p->y*sizeof(float));
	
	(proj->img[ip])->fomflag = 0;
	(proj->img[ip])->fom = NULL;
	(proj->img[ip])->image[0].ox = neworigin.x;
	(proj->img[ip])->image[0].oy = neworigin.y;
	(proj->img[ip])->image[0].oz = 0;

	return(0);
}

/************************************************************************
@Function: prjct_locate_point
@Description:
	Calculates projection of a point in the 3D space.
@Algorithm:
@Arguments:
	Vector3 p			point in the volumetric space.
	Vector3 origin		origin of the 3D system.
	View view			reference views.
@Returns:
	Vector3 				calculated projection
**************************************************************************/
Vector3 	prjct_locate_point(Vector3 p, Vector3 origin, View view)
{

	double 		oldx, oldy, oldz;
	double		newx, newy;

	Matrix3	rotat_mat = matrix3_transpose(matrix3_from_view(view));

	if ( verbose & VERB_FULL ) {
		printf("\nProjecting 3D point:\n");
		printf("Rotation origin:                %g %g %g pixels\n",
				origin.x, origin.y, origin.z);
		matrix3_show(rotat_mat);
		printf("\n");
	}
 	
	oldz = p.z - origin.z;
	oldy = p.y - origin.y;
	oldx = p.x - origin.x;
	newy = oldx*rotat_mat.r10 + oldy*rotat_mat.r11 + oldz*rotat_mat.r12 + origin.y;
	newx = oldx*rotat_mat.r00 + oldy*rotat_mat.r01 + oldz*rotat_mat.r02 + origin.x;

	Vector3 pp = {newx,newy,0.};
	
	return(pp);
}
