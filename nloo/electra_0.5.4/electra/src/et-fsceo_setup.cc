/*
	et-fsceo_setup.cc
	Prepare data for even/odd FSC analysis
	Author: Giovanni Cardone
	Created: 20041215 	Modified: 20050329
*/

#include "bsoft.h"
#include "proj.h"
#include "tilt.h"
#include "ioroutines.h"
#include "imod.h"

// Usage assistance
char* use[] = {
" ",
"Usage: et-fsceo_setup -imod par.file [options] [output.vol] ",
"------------------------------------------------------------------",
"Generate the two reconstructions, each from half the projections,",
"for the even/odd FSC analysis.",
"If a name for the output volumes is given, this is used instead of that in the imod",
"parameter file. The name of the reconstructions from the odd and the even",
"projections are determined by appending the strings '_0' and '_1', respectively,",
"to the given name of the output volume.",
" ",
"Note: at the moment only imod is supported.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output",
"-imod tilt.com           imod reconstruction code selected, with related input file",
"-task 1                  selection of subvolume to generate",
"                         (default: both",
"                           1=only volume from odd projections",
"                           2=only volume from even projections )",
"-datatype u              Force writing of a new data type (byte, short and float supported).",
" ",
"Output:",
"-intermediate            keep intermediate files",
" ",
NULL
};

int 		main(int argc, char **argv)
{

	// Initialize variables
	
	DataType 	datatype = Unknown_Type;		// Conversion to new type
	char*		in_pfile = NULL;				// input projections filename
	char*		out_pfile = NULL;			// output projections filename
	char*		out_pfx_vfile = NULL;			// output prefix volume filename
	char*		out_vfile = NULL;			// output volume filename
	char*		tlt_file = NULL;				// tilt angles file
	char*		out_tlt_file = NULL;			// output prefix tilt angles file
	char*		imod_parfile = NULL;			// imod parameters file
	
	Bimage*		bimg = NULL;					// individual image
	Proj*        p = NULL;					// projections
	
	char*        fname_imdtlt = NULL;			// imod parameters filename
	Imod_tilt*   imd_tpars = NULL;			// imod tilt parameters
	Imod_tilt*   imd_tpars_out = NULL;		// imod tilt parameters

	int			ntlt_angs = 0;				// number of tilt angles
	float*		tlt_angs = NULL;				// tilt series
	
	int			intermediate = 0;				// intermediate files flag
	int			itask = 0;						// volumes to generate:
												// 0=from both odd and even projections;
												// 1=from odd projections only;
												// 2=from even projections only;
	
	int			np, np_half, mnp, pnp;
	int			optind;

	Option*		option = get_option_list(use, argc, argv, &optind);
	Option*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
 		if ( strcmp(curropt->tag, "intermediate") == 0 )
			intermediate = 1;
 		if ( strcmp(curropt->tag, "imod") == 0 )
			imod_parfile = get_option_filename(curropt->value);
 		if ( strcmp(curropt->tag, "datatype") == 0 ) {
			datatype = get_option_datatype(curropt->value);
 			if ( datatype != UChar && datatype != Short && datatype != Float ) {
 				fprintf(stderr, "-datatype: only byte, short and float are supported\n");
 				exit(-1);
 			}				
 		}
 		if ( strcmp(curropt->tag, "task") == 0 )
			if ( sscanf(curropt->value, "%d", &itask) < 1 || (itask<1 || itask>2))
				fprintf(stderr, "-task: a value between 1 and 2 must be specified");
    }
	option_kill(option);

	if ( verbose ) {
		printf("\n-----------------------------------------------\n");
		printf("                   et-fsceo_setup\n");
		printf("    Prepare data for even/odd FSC analysis\n");
		printf("-----------------------------------------------\n\n");
	}

	if ( optind < argc ) out_pfx_vfile = copystring(argv[optind]);

	// check if selected imod
	if ( imod_parfile == NULL ) {
		fprintf(stderr, "Error: imod selection needed!\n");
		exit(-1);
	}

	if ( verbose & VERB_PROCESS )
		printf("Reading imod paramter file %s\n",imod_parfile);

	imd_tpars = imd_tilt_read(imod_parfile);
	if ( imd_tpars == NULL ) exit(-1);
	imd_tilt_checkvol(imd_tpars);
	imd_tilt_setparallel(imd_tpars);
	if ( datatype != Unknown_Type ) imd_tilt_putdatatype(imd_tpars, datatype);

	// get file names and initialize
	tlt_file = imd_tilt_gettiltfile(imd_tpars);
	if ( verbose & VERB_PROCESS )
		printf("Reading tilt angles file %s\n",tlt_file);
	tlt_angs = tlt_load_tilt_angles(tlt_file,&ntlt_angs);

	if ( out_pfx_vfile == NULL) out_pfx_vfile = imd_tilt_getvolfile(imd_tpars);

	in_pfile = imd_tilt_getprojfile(imd_tpars);
	bimg = io_read_image(in_pfile, 0, -1);
	
	if ( datatype == Unknown_Type ) datatype = bimg->datatype;
	
	// check consistency
	if( ntlt_angs != (int) ((bimg->n==1&&bimg->z>1)?bimg->z:bimg->n)) {
		fprintf(stderr,"Error: number of tilt angles in %s and number of projections in %s are not consistant!\n",
			tlt_file, in_pfile);
		exit(-1);
	}
	
	np = ntlt_angs;
	
	kill_img(bimg);
	
	
	for ( int j = 0; j < 2; j++) {

		// skip unrequested reconstructions
		if (itask!=0 && itask!=(j+1) ) continue;

		np_half = ( np + 1 - j) / 2;
		p = proj_init(np_half);
 
		pnp = 0;
		mnp = 0;
		for ( int l=j; l<np; l+=2) {
		
			if (imd_tilt_angle_missing(imd_tpars,l+1)) {
				mnp++;
				continue;
			}

			bimg = io_read_image(in_pfile,1,l);
			
			if (pnp==0)
				proj_cpy_img_sett(p,bimg);

			proj_puttilt(p,tlt_angs[l],pnp);
			proj_putimage(p,bimg,pnp,pnp);
		
			pnp++;
		}

		if ((pnp+mnp)!=np_half) {
			fprintf(stderr,"Error while generating half projections (%d of 2)\n",j);
			exit(-1);
		}

		if ( verbose & VERB_FULL ) {
			printf("Collected %d images for set of projections # %d of 2\n",np_half,j+1);
			proj_display(p);
		}

		// initialize output file names
		fname_imdtlt = insert_in_filename("tilt.com",filename_base(number_filename(in_pfile,j,1)),'_');
		out_pfile = number_filename(in_pfile,j,1);
		out_vfile = number_filename(out_pfx_vfile,j,1);
		out_tlt_file = number_filename(tlt_file,j,1);

		// write tilt angles and projections to file
		proj_write_angles(out_tlt_file, p, 0);
		io_pwrite(out_pfile,p,IO_ALL,0);
		
		// initialize reconstruction and run
		imd_tpars_out = (Imod_tilt *) balloc(sizeof(Imod_tilt));
		imd_tilt_cpy_and_modify(imd_tpars_out, imd_tpars, fname_imdtlt,
			out_pfile, out_vfile, -1, -1, -1, 100., out_tlt_file) ;
		imd_tilt_rmxcldlist(imd_tpars_out);
		if ( imd_tilt_writepars(imd_tpars_out) < 0 ) exit(-1);
		if ( imd_tilt_execute(imd_tpars_out) < 0 ) exit(-1);

		if(!intermediate) {
			// delete intermediate files
			if (imd_tilt_rmallfiles_butvol(imd_tpars_out)<0) exit(-1);
		}

		bfree_string(out_tlt_file);
		bfree_string(out_vfile);
		bfree_string(out_pfile);
		bfree_string(fname_imdtlt);
		bfree(imd_tpars_out,sizeof(Imod_tilt));
		
		proj_kill(p);
	}

	bfree(imd_tpars,sizeof(Imod_tilt));
	bfree_string(in_pfile);
	bfree_string(out_pfx_vfile);
	bfree(tlt_angs,ntlt_angs*sizeof(float));
	bfree_string(tlt_file);

	if ( verbose ) {
		printf("-----------------------------------------------\n");
		printf("       et-fsceo_setup concluded successfully!\n");
		printf("-----------------------------------------------\n");
	}
	return(0);
}

