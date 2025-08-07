/*
	symmetry.c
	General symmetry functions
	Author: Bernard Heymann
	Created: 20010420 	Modified: 20041229
*/

#include "symmetry.h"
#include "matrix.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated

/************************************************************************
@Function: init_point_group_symmetry
@Description:
	Sets up a symmetry structure from a string identifier.
@Algorithm:
	The point group symmetries are identified by the following strings:
		C<n>			cyclic point group of order n.
		D<n>			dihedral point group of order n.
		T				tetrahedral point group.
		O				octahedral/cubic point group.
		I				icosahedral/dodecahedral point group.
		H<r>,<a>,<d>	helical symmetry with rise r, rise angle a and dyad d (1/2).
	If the point group string is empty (NULL), the default is C1 (asymmetric).
@Arguments:
	char* symmetry_string	string containing point group identifier.
@Returns:
	Bsymmetry* sym			new symmetry structure.
**************************************************************************/
Bsymmetry* 	init_point_group_symmetry(char* symmetry_string)
{
	int 			i, theorder, asym_flag = 0, helix_dyad = 1;
	double			helix_rise = 1, helix_angle = PI;
	double			golden = (1.0 + sqrt(5.0))/2.0;
	Matrix3 		mat = {0,-1,0,1,0,0,0,0,1}; // 90 degree rotation for I90
	
	Bsymmetry*		sym = (Bsymmetry *) balloc(sizeof(Bsymmetry));
	
	if ( !symmetry_string ) asym_flag = 1;
	else if ( strlen(symmetry_string) < 1 ) asym_flag = 1;
	else {
		symmetry_clean_pointgroup(symmetry_string);
		if ( sscanf(symmetry_string, "C%d", &theorder) && theorder == 1 ) asym_flag = 1;
	}
	
	if ( asym_flag ) {
		strcpy(sym->label, "C1");
		sym->n = 1;
		sym->point = 101;
		sym->op = (Bsymop *) balloc(sizeof(Bsymop));
		sym->op[0].axis.z = 1;
		sym->op[0].order = 1;
		return(sym);
	}
	
	symmetry_clean_pointgroup(symmetry_string);
	
	if ( verbose & VERB_PROCESS )
		printf("\nGetting symmetry operators for %s\n\n", symmetry_string);
	
	strcpy(sym->label, symmetry_string);
		
	if ( sscanf(symmetry_string, "C%d", &theorder) > 0 ) { 		// Cyclic point groups
		sym->n = 1;
		sym->point = 100 + theorder;
		sym->op = (Bsymop *) balloc(sym->n*sizeof(Bsymop));
		sym->op[0].axis.z = 1;
		sym->op[0].order = theorder;
		sym->op[0].angle = PI*2.0/theorder;
	} else if ( sscanf(symmetry_string, "D%d", &theorder) > 0 ) {	// Dihedral point groups
		sym->n = 2;
		sym->point = 200 + theorder;
		sym->op = (Bsymop *) balloc(sym->n*sizeof(Bsymop));
		sym->op[0].axis.z = 1;
		sym->op[0].order = theorder;
		sym->op[0].angle = PI*2.0/theorder;
		sym->op[1].axis.x = 1;										// Rotate around 2-fold x-axis
		sym->op[1].order = 2;
		sym->op[1].angle = PI;
	} else if ( symmetry_string[0] == 'T' ) {						// Tetrahedral point group
		sym->n = 3;
		sym->point = 320;
		sym->op = (Bsymop *) balloc(sym->n*sizeof(Bsymop));
		sym->op[0].axis.x = sym->op[0].axis.y = 
				sym->op[0].axis.z = 1;								// Rotate around 3-fold (1,1,1)
		sym->op[0].order = 3;
		sym->op[0].angle = PI*2.0/3.0;
		sym->op[1].axis.z = 1;										// Rotate around 2-fold z-axis
		sym->op[1].order = 2;
		sym->op[1].angle = PI;
		sym->op[2].axis.x = 1;										// Rotate around 2-fold x-axis
		sym->op[2].order = 2;
		sym->op[2].angle = PI;
	} else if ( symmetry_string[0] == 'O' ) {						// Octahedral point group
		sym->n = 3;
		sym->point = 432;
		sym->op = (Bsymop *) balloc(sym->n*sizeof(Bsymop));
		sym->op[0].axis.x = sym->op[0].axis.y = 
				sym->op[0].axis.z = 1;								// Rotate around 3-fold (1,1,1)
		sym->op[0].order = 3;
		sym->op[0].angle = PI*2.0/3.0;
		sym->op[1].axis.z = 1;										// Rotate around 4-fold z-axis
		sym->op[1].order = 4;
		sym->op[1].angle = PI/2.0;
		sym->op[2].axis.x = 1;										// Rotate around 2-fold x-axis
		sym->op[2].order = 2;
		sym->op[2].angle = PI;
	} else if ( symmetry_string[0] == 'I' ) {						// Icosahedral point group
		sym->n = 4;
		sym->point = 532;
		sym->op = (Bsymop *) balloc(sym->n*sizeof(Bsymop));
		sym->op[0].axis.x = sym->op[0].axis.y = 
				sym->op[0].axis.z = 1;								// Rotate around 3-fold (1,1,1)
		sym->op[0].order = 3;
		sym->op[0].angle = PI*2.0/3.0;
		sym->op[1].axis.x = 1;										// Rotate around 2-fold
		sym->op[1].axis.y = 1/golden;
		sym->op[1].axis.z = golden;
		sym->op[1].order = 2;
		sym->op[1].angle = PI;
		sym->op[2].axis.x = 1/golden;								// Rotate around 5-fold
		sym->op[2].axis.y = 1;
		sym->op[2].order = 5;
		sym->op[2].angle = PI*2.0/5.0;
		sym->op[3].axis.z = 1;										// Rotate around 2-fold z-axis
		sym->op[3].order = 2;
		sym->op[3].angle = PI;
		if ( strstr(symmetry_string, "I90") )
			for ( i=0; i<sym->n; i++ )
				sym->op[i].axis = vector3_matrix3_multiply(mat, sym->op[i].axis);
	} else if ( symmetry_string[0] == 'H' ) {						// Helical symmetry
		sscanf(symmetry_string, "H%lf,%lf,%d", &helix_rise, &helix_angle, &helix_dyad);
		sym->n = 1;
		if ( helix_dyad == 2 ) sym->n = 2;
		sym->op = (Bsymop *) balloc(sym->n*sizeof(Bsymop));
		sym->point = 600;
		sym->op[0].axis.x = sym->op[0].axis.y = 0;
		sym->op[0].axis.z = 1;										// Helical axis
		sym->op[0].shift = helix_rise;								// Rise in angstrom
		sym->op[0].angle = helix_angle*PI/180.0;
		if ( helix_dyad == 2 ) {
			sym->op[1].axis.y = sym->op[0].axis.z = 0;
			sym->op[1].axis.x = 1;									// Dyad axis
			sym->op[1].order = 2;
			sym->op[1].angle = PI;
		}
	} else {
		printf("Error: This is not a valid symmetry designation: %s\n",
				symmetry_string);
		exit(-1);
	}
	
	for ( i=0; i<sym->n; i++ )
		sym->op[i].axis = vector3_normalize(sym->op[i].axis);
	
	return(sym);
}

/************************************************************************
@Function: kill_symmetry
@Description:
	Deallocates a symmetry structure.
@Algorithm:
	.
@Arguments:
	Bsymmetry* sym		symmetry structure.
@Returns:
	int					0.
**************************************************************************/
int 	 	kill_symmetry(Bsymmetry* sym)
{
	if ( !sym ) return(0);
	
	if ( sym->n > 0 && sym->op ) bfree(sym->op, sym->n*sizeof(Bsymop));
	
	bfree(sym, sizeof(Bsymmetry));
	
	return(0);
}

/************************************************************************
@Function: symmetry_clean_pointgroup
@Description:
	Corrects an existing point group string.
@Algorithm:
	A point group string is converted to .
@Arguments:
	char* symmetry_string	string containing point group identifier.
@Returns:
	int 					0.
**************************************************************************/
int			symmetry_clean_pointgroup(char* symmetry_string)
{
	unsigned int 		i, j;
	
	if ( !symmetry_string ) return(-1);
	
	// Remove leading blanks and convert to upper case
	for ( i=0, j=0; i<strlen(symmetry_string); i++ )
		if ( isalnum(symmetry_string[i]) )
			symmetry_string[j++] = toupper(symmetry_string[i]);
	symmetry_string[j] = 0;
	
	// Alternate nomenclature for point groups
	if ( isdigit(symmetry_string[0]) ) {
		sscanf(symmetry_string, "%d", &j);
		if ( strstr(symmetry_string, "532") && strstr(symmetry_string, "90") )
			strcpy(symmetry_string, "I90");
		else if ( strstr(symmetry_string, "532") )
			strcpy(symmetry_string, "I");
		else if ( strstr(symmetry_string, "432") )
			strcpy(symmetry_string, "O");
		else if ( strstr(symmetry_string, "23") )
			strcpy(symmetry_string, "T");
		else if ( j > 1 && symmetry_string[1] == 2 )
			sprintf(symmetry_string, "D%d", j);
		else if ( j > 0 )
			sprintf(symmetry_string, "C%d", j);
		else
			strcpy(symmetry_string, "C1");
	}
	
	if ( sscanf(symmetry_string, "D%d", &i) && i < 2 )
		strcpy(symmetry_string, "C1");

	return(0);
}

/************************************************************************
@Function: symmetry_get_all_views
@Description:
	Get all symmetry-related views of one given view.
@Algorithm:
	The number of views generated for a point group symmetry is
	calculated as the product of the order fields in the symmetry
	structure.
@Arguments:
	Bsymmetry* sym		symmetry structure.
	View asu_view		asymmetric unit vector and rotation angle;
@Returns:
	View*				array of n views: vectors and rotation angles.
	View*				linked list of views.
**************************************************************************/
View*		symmetry_get_all_views(Bsymmetry* sym, View asu_view)
{
	int				i, j, k;
	double			angle;
	Quaternion		q, qv; 
	
	View*			view = NULL;
	View*			v;
	View*			vn = (View *) add_item((char **) &view, sizeof(View));
	*vn = asu_view;
	vn->next = NULL;
	
	int				nview = 1;

	for ( i=0; i<sym->n; i++ ) {
		for ( j=1; j<sym->op[i].order; j++ ) {
			angle = j*TWOPI/sym->op[i].order;
			q = quaternion_from_angle_and_axis3(angle, sym->op[i].axis);
			for ( k=0, v=view; k<nview; k++, v=v->next ) {
				qv = quaternion_from_view(*v);
				qv = quaternion_multiply(q, qv);
				vn = (View *) add_item((char **) &view, sizeof(View));
				*vn = view_from_quaternion(qv);
			}
		}
		nview *= sym->op[i].order;
	}
	
	if ( verbose & VERB_FULL ) {
		for ( i=1, v=view; v; v=v->next, i++ ) {
			printf("View %5d:\t", i);
			show_view(*v);
		}
	}
	
	return(view);
}

/************************************************************************
@Function: asymmetric_unit_views
@Description:
	Initializes a well-distributed set of views in an asymmetric unit.
@Algorithm:
	A set of views is calculated with tesselation within each asymmetric
	unit such that the views are well-distributed.
	If the full flag is set, both halves of the asymmetric unit are covered.
@Arguments:
	Bsymmetry* sym		symmetry structure.
	float theta_step	angular step size from primary symmetry axis (radians).
	float phi_step		angular step size around primary symmetry axis (radians).
	int full			flag for generating a full asymmetric unit (default 0=half).
@Returns:
	View* 				a linked list of views.
**************************************************************************/
View* 		asymmetric_unit_views(Bsymmetry* sym, float theta_step, float phi_step, int full)
{
	int 			i, j, n = 0, ntheta, nphi, nrphi, istart;
	double			theta, phi, max_theta, theta_start, phi_start, phi_end;
	double			g = (1 + sqrt(5.0))/2;
	
	max_theta = PI/2.0 + 1e-6;			// Most symmetries limited to upper half
	ntheta = (int) (max_theta/theta_step + 0.5);
	nphi = (int) (PI*2.0/(sym->op[0].order*phi_step) + 0.5);
	
	if ( sym->point < 200 ) {
		if ( full ) max_theta = PI + 1e-6;
		ntheta += ntheta;
	} else if ( sym->point == 432 ) {
		max_theta = PI/4.0 + 1e-6;
		ntheta = (int) (max_theta/theta_step + 0.5);
	} else if ( sym->point == 532 ) {
		max_theta = atan(1/(g+1)) + 1e-6;
		ntheta = (int) (max_theta/theta_step + 2.5);
		nphi = (int) (atan(1/g)/phi_step + 1.5);
	} else if ( sym->point == 600 ) {
		ntheta = 1;
		nphi = (int) (PI*2.0/phi_step + 1);
	}
	
	if ( ntheta < 2 ) ntheta = 2;
	if ( nphi < 2 ) nphi = 2;
		
	if ( verbose & (VERB_PROCESS | VERB_LABEL) ) {
		if ( full )
			printf("Getting all the asymmetric unit views for symmetry %s\n", 
					sym->label);
		else
			printf("Getting half the asymmetric unit views for symmetry %s\n", 
					sym->label);
		printf("Theta and phi step sizes:       %g %g degrees\n", 
				theta_step*180/PI, phi_step*180/PI);
		printf("Theta and phi steps:            %d %d\n", 
				ntheta, nphi);
	}
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG asymmetric_unit_views: Views allocated = %d\n", 2*ntheta*nphi);
	
	View*			view = NULL;	
	View*			v = NULL;
	
	if ( sym->point < 500 ) {						// Top view for all but icosahedral and helical symmetry
		v = (View *) add_item((char **) &view, sizeof(View));
		v->z = 1;
		n = 1;
		if ( verbose & VERB_DEBUG )
			printf("DEBUG asymmetric_unit_views: Adding the top view\n");
	}
	
	if ( sym->point > 100 && sym->point < 200 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			nrphi = (int) (nphi*sin(theta)/2 + 0.5);	// Number of views at this radius and theta
			phi = 0;
			if ( nrphi ) phi = PI*1.0/(sym->op[0].order*nrphi);
			for ( j=-nrphi; j<=nrphi; j++ ) {
				v = (View *) add_item((char **) &view, sizeof(View));
				v->x = sin(theta)*cos(j*phi);
				v->y = sin(theta)*sin(j*phi);
				v->z = cos(theta);
				if ( fabs(v->z) < 1e-6 ) v->z = 0;
				v->a = j*phi;
				n++;
			}
		}
		if ( full && v->z > -0.999999 ) {
			v = (View *) add_item((char **) &view, sizeof(View));
			v->z = -1;
			n++;
		}
	} else if ( sym->point > 200 && sym->point < 300 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			nrphi = (int) (nphi*sin(theta)/2 + 0.5);	// Number of views at this radius and theta
			phi = 0;
			if ( nrphi ) phi = PI*1.0/(sym->op[0].order*nrphi);
			istart = 0;
			if ( full ) istart = -nrphi;
			for ( j=istart; j<=nrphi; j++ ) {
				v = (View *) add_item((char **) &view, sizeof(View));
				v->x = sin(theta)*cos(j*phi);
				v->y = sin(theta)*sin(j*phi);
				v->z = cos(theta);
				if ( v->z < 0 ) v->z = 0;
				v->a = j*phi;
				n++;
			}
		}
	} else if ( sym->point == 320 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			phi_start = 0;
			phi_end = theta;
			if ( theta > PI/4 ) phi_end = PI/2 - theta;
			if ( full ) phi_start = -theta;
			if ( full && theta > PI/4 ) phi_start = theta - PI/2;
			for ( phi=phi_start; phi<=phi_end; phi+=phi_step ) {
				v = (View *) add_item((char **) &view, sizeof(View));
				v->x = sin(theta);
				v->y = sin(phi);
				v->z = cos(theta);
				n++;
			}
		}
	} else if ( sym->point == 432 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			phi_start = 0;
			if ( full ) phi_start = -theta;
			for ( phi=phi_start; phi<=theta; phi+=phi_step ) {
				v = (View *) add_item((char **) &view, sizeof(View));
				v->x = sin(theta);
				v->y = sin(phi);
				v->z = cos(theta);
				n++;
			}
		}
	} else if ( sym->point == 532 ) {
		if ( strstr(sym->label, "I90") ) {
			theta_start = 0;
			if ( full ) theta_start = -ntheta*theta_step;
			for ( theta=theta_start; theta<=max_theta; theta+=theta_step ) {
				for ( phi=0; phi<=g*(max_theta - fabs(theta)); phi+=phi_step ) {
					v = (View *) add_item((char **) &view, sizeof(View));
					v->x = tan(phi);
					v->y = tan(theta);
					v->z = cos(theta);
					n++;
				}
			}
		} else {
			for ( theta=0; theta<=max_theta; theta+=theta_step ) {
				nrphi = (int) (g*(max_theta - theta)/phi_step);
				phi_start = 0;
				if ( full ) phi_start = -nrphi*phi_step;
				for ( phi=phi_start; phi<=g*(max_theta - theta); phi+=phi_step ) {
					v = (View *) add_item((char **) &view, sizeof(View));
					v->x = tan(theta);
					v->y = tan(phi);
					v->z = cos(theta);
					n++;
				}
			}
		}
	} else if ( sym->point == 600 ) {					// Helical symmetry
		for ( phi = 0; phi < 2.0*PI-0.001; phi += phi_step ) {
			v = (View *) add_item((char **) &view, sizeof(View));
			v->x = cos(phi);
			v->y = sin(phi);
			v->z = 0;
			v->a = phi;
			n++;
		}
	} else {
		fprintf(stderr, "Warning: Symmetry type %d not supported!\n", sym->point);
	}

	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) )
		printf("Views generated:                %d\n", n);

	for ( v=view; v; v = v->next ) view_normalize(v);
	
	if ( verbose & VERB_PROCESS ) {
		printf("View\tx\ty\tz\ta\n");
		for ( v=view, i=1; v; v = v->next, i++ ) {
			printf("%d\t", i);
			show_view(*v);
		}
		printf("\n");
	}
	
	return(view);
}

/************************************************************************
@Function: change_views_to_asymmetric_unit
@Description:
	Change the views to those in the asymmetric unit.
@Algorithm:
	The view is replaced with the one in the standard asymmetric unit.
@Arguments:
	Bsymmetry* sym		symmetry structure.
	View* view			list of views (replaced).
@Returns:
	int					0.
**************************************************************************/
int 		change_views_to_asymmetric_unit(Bsymmetry* sym, View* view)
{
	if ( sym->point == 101 ) return(0);
	
	View*			v;

	for ( v=view; v; v=v->next )
		find_asymmetric_unit_view(sym, v);
	
	return(0);
}

/************************************************************************
@Function: find_asymmetric_unit_view
@Description:
	Finds the corresponding view in the asymmetric unit.
@Algorithm:
	The view is replaced with the one in the standard asymmetric unit.
@Arguments:
	Bsymmetry* sym		symmetry structure.
	View* v				view (replaced).
@Returns:
	int					0.
**************************************************************************/
int 		find_asymmetric_unit_view(Bsymmetry* sym, View* theview)
{
	if ( sym->point == 101 ) return(0);
	
	View*			v, bv = {NULL,0,0,0,0};
	View*			view = symmetry_get_all_views(sym, *theview);
//	int			nview = (int) count_list((char *) view);
	
	double			tol = 1e-6;
	double			golden = (1.0 + sqrt(5.0))/2.0;
	double			lim = tan(PI*1.0/sym->op[0].order);
	
	for ( v=view; v; v=v->next ) if ( v->x + tol >= 0 ) {
		if ( sym->point < 200 ) {
			if ( sym->op[0].order == 2 ) {
				bv = *v;
			} else {
				if ( fabs(v->y) <= v->x*lim + tol ) bv = *v;
			}
		} else if ( sym->point > 200 && sym->point < 300 ) {
			if ( v->z + tol >= 0 ) {
				if ( sym->op[0].order == 2 ) {
					bv = *v;
				} else {
					if ( fabs(v->y) <= v->x*lim + tol ) bv = *v;
				}
			}
		} else if ( sym->point == 320 ) {
			if ( fabs(v->y) <= v->x + tol && fabs(v->y) <= v->z + tol ) bv = *v;
		} else if ( sym->point == 432 ) {
			if ( fabs(v->y) <= v->x + tol && v->x <= v->z + tol ) bv = *v;
		} else if ( sym->point == 532 ) {
			if ( strstr(sym->label, "I90") ) {
				if ( v->z + tol >= 0 && fabs(v->y) <= (v->z/golden - v->x)/golden + tol ) bv = *v;
			} else {
				if ( v->z + tol >= 0 && fabs(v->y) <= v->z/golden - v->x*golden + tol ) bv = *v;
			}
		} else if ( sym->point == 600 ) {	// What is it for helical symmetry?
		}
	}
	
	kill_list((char *) view, sizeof(View));
	
	if ( view_vector_size(bv) < 0.9 ) {
		fprintf(stderr, "ASU view not found: %g %g %g %g\n",
			theview->x, theview->y, theview->z, theview->a*180.0/PI);
	} else {
		if ( verbose & VERB_FULL ) {
			printf("Old view:\t");
			show_view(*theview);
			printf("New view:\t");
			show_view(bv);
		}
		theview->x = bv.x;
		theview->y = bv.y;
		theview->z = bv.z;
		theview->a = bv.a;
	}
	
	return(0);
}

/************************************************************************
@Function: calc_skew_matrix
@Description:
	Calculates the skew matrix from the unit cell parameters.
@Algorithm:
	Derived from the X-plor source rotate.s.
	New coordinates are obtained by r'(i)=sum_j matrix(i,j)*r(j)
	The convention to setup the matrices is as follows:
		a lies on the x-axis, 
		b lies in the x,y plane. 
@Arguments:
	UnitCell unit_cell	6-value unit cell parameters.
@Returns:
	float* 				the 3x3 + 3 skew matrix.
**************************************************************************/
float*  	calc_skew_matrix(UnitCell unit_cell)
{
	if ( verbose & VERB_DEBUG )
		printf("DEBUG: Calculating the fractionalization matrix\n");
	
	float*		skew = (float *) balloc(12*sizeof(float));
    
    if ( unit_cell.alf < 0 ) unit_cell.alf = -unit_cell.alf;
    if ( unit_cell.bet < 0 ) unit_cell.bet = -unit_cell.bet;
    if ( unit_cell.gam < 0 ) unit_cell.gam = -unit_cell.gam;
    if ( unit_cell.alf > PI ) unit_cell.alf *= PI/180;
    if ( unit_cell.bet > PI ) unit_cell.bet *= PI/180;
    if ( unit_cell.gam > PI ) unit_cell.gam *= PI/180;

    float   	cterm = (cos(unit_cell.bet)*cos(unit_cell.gam)-
							cos(unit_cell.alf))/
    	    	    	(sin(unit_cell.bet)*sin(unit_cell.gam));
    float   	sterm = sqrt(1.0 - cterm*cterm);
    
    skew[0] = 1.0/unit_cell.a;
    skew[1] = -cos(unit_cell.gam)/(sin(unit_cell.gam)*unit_cell.a);
    skew[2] = -(cos(unit_cell.gam)*sin(unit_cell.bet)*cterm+
				cos(unit_cell.bet)*sin(unit_cell.gam))/
     			(sin(unit_cell.bet)*sterm*sin(unit_cell.gam)*unit_cell.a);
    skew[3] = 0;
    skew[4] = 1.0/(sin(unit_cell.gam)*unit_cell.b);
    skew[5] = cterm/(sterm*sin(unit_cell.gam)*unit_cell.b);
    skew[6] = 0;
    skew[7] = 0;
    skew[8] = 1.0/(sin(unit_cell.bet)*sterm*unit_cell.c);
	
	for ( int i=0; i<9; i++ )
		if ( skew[i] < 1e-6 ) skew[i] = 0;
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG: Fractionalization matrix: %10.4f %10.4f %10.4f\n", skew[0], skew[1], skew[2]);
		printf("                                 %10.4f %10.4f %10.4f\n", skew[3], skew[4], skew[5]);
		printf("                                 %10.4f %10.4f %10.4f\n", skew[6], skew[7], skew[8]);
	}
	
    return(skew);
}

/************************************************************************
@Function: sym_matrices_from_text_list
@Description:
	Calculates symmetry matrices from a list of strings.
@Algorithm:
	The list of strings is expected to be packed into a single character
	array with a fixed length for each string. Each string encodes a
	symmetry operation in terms of x, y and z operations in reciprocal
	space. 
@Arguments:
	int nsym			number of symmetry operators.
	char* symop			array of symmetry operator lines.
	int line_len		length of text line in the array.
@Returns:
	float* 				a set of 12-value symmetry matrices.
**************************************************************************/
float*		sym_matrices_from_text_list(int nsym, char* symop, int line_len)
{
	// Set up the symmetry operator matrix
	int 		i, j, k, l;
	float*		mat = (float *) balloc(nsym*12*sizeof(float));

	char		op[200];
	for ( i=0; i<nsym; i++ ) {
		if ( verbose & VERB_DEBUG )
			printf("DEBUG sym_matrices_from_text_list: Symmetry operator %d:\n", i+1);
		k = 0;
		for ( j=0; j<3; j++ ) {
			l = 0;
			memset(op, '\0', line_len);
			while ( k<line_len && symop[i*line_len+k] != ',' ) {
				op[l] = tolower(symop[i*line_len+k]);
				k++;
				l++;
			}
			op[l] = 0;
			if ( strstr(op,"-x") ) mat[i*12+4*j] = -1;
			else if ( strstr(op,"x") ) mat[i*12+4*j] = 1;
			if ( strstr(op,"-y") ) mat[i*12+4*j+1] = -1;
			else if ( strstr(op,"y") ) mat[i*12+4*j+1] = 1;
			if ( strstr(op,"-z") ) mat[i*12+4*j+2] = -1;
			else if ( strstr(op,"z") ) mat[i*12+4*j+2] = 1;
			if ( strstr(op,"1/2") ) mat[i*12+4*j+3] = 0.5;
			if ( strstr(op,"1/4") ) mat[i*12+4*j+3] = 0.25;
			if ( strstr(op,"3/4") ) mat[i*12+4*j+3] = 0.75;
			if ( strstr(op,"1/3") ) mat[i*12+4*j+3] = 1.0/3.0;
			if ( strstr(op,"2/3") ) mat[i*12+4*j+3] = 2.0/3.0;
			if ( strstr(op,"1/6") ) mat[i*12+4*j+3] = 1.0/6.0;
			if ( strstr(op,"5/6") ) mat[i*12+4*j+3] = 5.0/6.0;
			k++;
			if ( verbose & VERB_DEBUG )
				printf("|%2.0f %2.0f %2.0f|   |%5.3f|\n", mat[i*12+4*j],
					mat[i*12+4*j+1], mat[i*12+4*j+2], mat[i*12+4*j+3]);
		}
	}
	
	return(0);
}

/************************************************************************
@Function: sym_show_matrices
@Description:
	Show symmetry matrices.
@Algorithm:
	. 
@Arguments:
	Bsymmetry* sym		symmetry structure.
@Returns:
	int 				number of symmetry matrices.
**************************************************************************/
int			sym_show_matrices(Bsymmetry* sym)
{
	int			i, nview = 1;
	for ( i=0; i<sym->n; i++ ) nview *= sym->op[i].order;

	View		ref = {NULL,0,0,1,0};
	View*		view = symmetry_get_all_views(sym, ref);
	Matrix3		mat;
	
	printf("\nSymmetry matrices:\n");
	for ( i=0; i<nview; i++ ) {
		mat = matrix3_from_view(view[i]);
		printf("Matrix %d:\n", i+1);
		matrix3_show(mat);
	}
	printf("\n");
	
	return(nview);
}
