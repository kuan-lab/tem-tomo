/*
	rwCCP4.c
	Functions for reading and writing CCP4 files
	Author: Bernard Heymann
	Created: 19990410 	Modified: 20030916
*/

#include "rwCCP4.h"
#include "rwsymop.h"
#include "symmetry.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

/************************************************************************
@Function: set_CCP4_machine_stamp
@Description:
	Setting a CCP4 style machine stamp.
@Algorithm:
	The 4-byte machine stamp:
				1	Big-endian IEEE 	17 17 00 00
                2	VAX 				34 65 00 00
				3	Cray				-
                4	Little-endian IEEE	68 65 00 00
                5	Convex				85 17 00 00
				6	Fijitsu VP			-
				(Note: not always implemented - so unreliable)
@Arguments:
	char* machine_stamp	machine stamp string.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int			set_CCP4_machine_stamp(char* machine_stamp)
{
	switch ( systype(0) ) {
		case BigIEEE:
			machine_stamp[0] = machine_stamp[1] = 17;
			break;
		case LittleIEEE:
			machine_stamp[0] = 68;
			machine_stamp[1] = 65;
			break;
		case LittleVAX:
			machine_stamp[0] = 34;
			machine_stamp[1] = 65;
			break;
		default:
			break;
	}
	
	return(0);
}

/************************************************************************
@Function: readCCP4
@Description:
	Reading a CCP4 map image file format.
@Algorithm:
	A 3D image format used in X-ray crystallography.
	Header size:				1024 bytes followed by the symmetry 
		operator table which is composed of 80 character lines, each 
		line for a symmetry operator.
	File format extensions:  	.map, .ccp, .ccp4
	The identifier is a 4-byte machine stamp:
				1	Big-endian IEEE 	17 17 00 00
                2	VAX 				34 65 00 00
				3	Cray				-
                4	Little-endian IEEE	68 65 00 00
                5	Convex				85 17 00 00
				6	Fijitsu VP			-
				(Note: not always implemented - so unreliable)
	Byte order determination:	Data type and third dimension values
								must be less than 256*256.
	Data types: 				0 = signed char, 1 = short, 2 = float,
								3 = complex short, 4 = complex float.
	Transform type: 			Centered hermitian
								The x-dimension contains the x-size
								of the full transform
@Arguments:
	Bimage*	p			the image structure.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 	readCCP4(Bimage* p)
{
    FILE        *fimg;
    if ( ( fimg = fopen(p->filename, "r") ) == NULL ) return(-1);

	CCP4head*	header = (CCP4head *) balloc(sizeof(CCP4head));
	if ( fread( header, CCP4SIZE, 1, fimg ) < 1 ) return(-2);
	
    // Determine byte order and swap bytes if from little-endian machine
    char*   	b = (char *) header;
    int     	i, swap = 0, vax = 0;
    if ( ( abs( header->mode ) > SWAPTRIG ) || ( abs(header->nz) > SWAPTRIG ) ) {
    	if ( verbose & VERB_PROCESS )
	    	fprintf(stderr, "Warning: Swapping header byte order for 4-byte types\n");
    	swap = 1;
		int 	extent = CCP4SIZE - 800; // exclude labels from swapping
    	for ( i=0; i<extent; i+=4 ) 
			if ( i < 208 || i > 212 ) swapbytes(b+i, 4);
    }
    
	if ( header->nx < 1 || header->ny < 1 || header->nz < 1 || header->mode < 0 ) {
		fprintf(stderr, "Error readCCP4: The file %s is not a valid CCP4 image!\n", p->filename);
		fclose(fimg);
		return(-3);
	}
	
    // Convert VAX floating point types if necessary
    // Machine stamp for vax is: 34 65 00 00
    if ( ( header->machst[0] == 34 ) || ( header->amin > header->amax ) ) {
    	if ( verbose & VERB_PROCESS )
	    	fprintf(stderr, "Warning: Converting VAX floating point\n");
    	vax = 1;
		for ( i=40; i<64; i+=4 )		// Unit cell parameters
    	    vax2ieee(b+i, swap);
		for ( i=76; i<88; i+=4 )		// Min, max, mean
    	    vax2ieee(b+i, swap);
		for ( i=100; i<204; i+=4 )		// Skew matrix and extra
    	    vax2ieee(b+i, swap);
    	vax2ieee(b+216, swap); 			// RMS
    }
    
	// Map the parameters
	p->x = header->nx;
	p->y = header->ny;
	p->z = header->nz;
	p->c = 1;
	p->colormodel = Gray;
	switch ( header->mode ) {
		case 0: p->datatype = SChar; break;
		case 1: p->datatype = Short; break;
		case 2: p->datatype = Float; break;
		case 3: p->datatype = ComplexShort; break;
		case 4: p->datatype = ComplexFloat; break;
		default: p->datatype = SChar; break;
	}
	p->offset = CCP4SIZE + header->nsymbt;
	if ( header->mode > 2 ) {
		p->transform = CentHerm;
		fseek(fimg, 0, SEEK_END);
		if ( ftell(fimg) > p->offset + 0.8*p->x*p->y*p->z*gettypesize(p->datatype) )
			p->x = 2*(p->x - 1);
		if ( header->mx%2 == 1 ) p->x += 1;	// Quick fix for odd x-size maps
	}
	p->min = header->amin;
	p->max = header->amax;
	p->avg = header->amean;
	p->std = header->arms;
	p->ua = header->a;
	p->ub = header->b;
	p->uc = header->c;
	p->alf = header->alpha;
	p->bet = header->beta;
	p->gam = header->gamma;
	if ( header->mx ) p->ux = header->a/header->mx;
	if ( header->my ) p->uy = header->b/header->my;
	if ( header->mz ) p->uz = header->c/header->mz;
	p->spacegroup = header->ispg;
	memcpy(p->label, header->labels, 800);
	i = 800;
	while ( --i > 0 && ( p->label[i] == ' ' || p->label[i] == 0 ) )
		p->label[i] = 0;
	
	// Allocating the single sub-image and setting its origin
	p->image = (Bsub_image *) balloc(sizeof(Bsub_image));
	p->image[0].ox = -header->nxStart;
	p->image[0].oy = -header->nyStart;
	p->image[0].oz = -header->nzStart;

	bfree(header, sizeof(CCP4head));
		
	if ( verbose & VERB_DEBUG )
		printf("DEBUG rwCCP4: Memory freed\n");

	p->data = img_read_data( fimg, p, -1, swap, vax, 0 );
	
	fclose(fimg);
	    
	if ( verbose & VERB_DEBUG )
		printf("DEBUG rwCCP4: Data read\n");

	return(0);
}

/************************************************************************
@Function: writeCCP4
@Description:
	Writing a CCP4 map image file format.
@Algorithm:
	A 3D image format used in X-ray crystallography.
@Arguments:
	Bimage*				the image structure.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int 	writeCCP4(Bimage* p)
{
	img_RGB2gray(p);
	
	long 			i, j;
	
    switch ( p->datatype ) {
    	case UChar: img_to_signed_byte(p); break;
    	case UShort: img_to_short(p); break;
    	case Int: img_to_float(p); break;
    	case ComplexInt: img_to_complex_float(p); break;
    	case Polar: img_polar2complex(p); break;
    	default: break;
    }
	
	if ( p->transform != NoTransform )
		img_convert_fourier(p, CentHerm);
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG rwCCP4: Writing CCP4 file\n");

	CCP4head*	header = (CCP4head *) balloc(sizeof(CCP4head));
	
	// Map the parameters
	strncpy(header->map, "MAP ", 4);
	set_CCP4_machine_stamp(header->machst);
	header->nx = (unsigned int) p->x;
	header->ny = (unsigned int) p->y;
	header->nz = (unsigned int) p->z;
	if ( p->transform == CentHerm ) header->nx = (int) (p->x/2 + 1);	// If a transform, physical storage is nx/2 + 1
	switch ( p->datatype ) {
		case SChar: header->mode = 0; break;
		case Short: header->mode = 1; break;
		case Float: header->mode = 2; break;
		case ComplexShort: header->mode = 3; break;
		case ComplexFloat: header->mode = 4; break;
		default: header->mode = 0; break;
	}
	header->nxStart = (int) (-p->image[0].ox - 0.5);
	header->nyStart = (int) (-p->image[0].oy - 0.5);
	header->nzStart = (int) (-p->image[0].oz - 0.5);
	header->mx = (int) (p->ua/p->ux + 0.5);
	header->my = (int) (p->ub/p->uy + 0.5);
	header->mz = (int) (p->uc/p->uz + 0.5);
	header->mapc = 1;
	header->mapr = 2;
	header->maps = 3;
	header->amin = p->min;
	header->amax = p->max;
	header->amean = p->avg;
	header->arms = p->std;
	header->a = p->ua;
	header->b = p->ub;
	header->c = p->uc;
	
	// This is a band-aid to overcome the limitations of the image format
	if ( fabs(p->ua - p->ux*header->mx) > 0.001 || fabs(p->ub - p->uy*header->my) > 0.001 ||
			fabs(p->uc - p->uz*header->mz) > 0.001 ) {
		header->a = p->ux*header->mx;
		header->b = p->uy*header->my;
		header->c = p->uz*header->mz;
		if ( verbose )
			fprintf(stderr, "Warning: Resetting the unit cell to: %g %g %g A\n", header->a, header->b, header->c );
	}
	
	UnitCell		unit_cell = {p->ua, p->ub, p->uc, p->alf, p->bet, p->gam};
	float*          frac_mat = calc_skew_matrix(unit_cell);
	header->alpha = p->alf*180/PI;
	header->beta = p->bet*180/PI;
	header->gamma = p->gam*180/PI;
	for ( i=0; i<9; i++ ) header->skwmat[i] = frac_mat[i];
	for ( i=0; i<3; i++ ) header->skwtrn[i] = frac_mat[i+9];
	header->ispg = p->spacegroup;
	bfree(frac_mat, 12*sizeof(float));
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG rwCCP4: Writing labels\n");

	i = 0;
	while ( ( j = strlen(p->label + i) ) ) {
		
		if ( verbose & VERB_DEBUG )
			printf("DEBUG rwCCP4: i=%ld  j=%ld  %s\n", i, j, p->label + i);

		if ( j > 80 ) j = 80;
		memcpy(header->labels + i, p->label + i, j);
		i += 80;
	}
	header->nlabl = (unsigned int) i/80;
	
	if ( verbose & VERB_DEBUG ) {
		printf("DEBUG rwCCP4: Nlabels = %d\n", header->nlabl);
		for ( i=0; i<header->nlabl; i++ )
			printf("%-80s\n", &header->labels[i*80]);
	}
	
	int			nsym = 0;
	char* 		symop = read_symop(NULL, p->spacegroup, &nsym);

	header->nsymbt = nsym*80;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG rwCCP4: Nsymbt = %d\n", header->nsymbt);
	
	p->offset = CCP4SIZE + header->nsymbt;
	
	unsigned long	datatypesize = gettypesize(p->datatype);
	unsigned long   datasize = (unsigned long) header->nx*header->ny*header->nz*datatypesize;
	
	if ( verbose & VERB_DEBUG )
		printf("DEBUG rwCCP4: Offset = %ld,  Typesize = %ld,  Datasize = %ld\n", 
				p->offset, datatypesize, datasize);
	
    FILE        *fimg;
    if ( ( fimg = fopen(p->filename, "w") ) == NULL ) return(-1);
	
	fwrite( header, CCP4SIZE, 1, fimg );
	fwrite( symop, header->nsymbt, 1, fimg );
	fwrite( p->data, datasize, 1, fimg );
	
	fclose(fimg);
	
	if ( symop ) bfree(symop, header->nsymbt*sizeof(char));
	bfree(header, sizeof(CCP4head));
		
	if ( verbose & VERB_DEBUG )
		printf("DEBUG rwCCP4: Done writing\n");
	
	return(0);
}

