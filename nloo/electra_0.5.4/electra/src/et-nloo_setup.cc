/*
	et-nloo_setup.cc
	Prepare data for Noise-compensated Leave-One-Out analysis
	Author: Giovanni Cardone
	Created: 20050118	Modified: 20070718
*/

#include "bsoft.h"
#include "proj.h"
#include "shape.h"
#include "project.h"
#include "tilt.h"
#include "ioroutines.h"
#include "imod.h"
#include "img_op.h"
#include "buffer.h"

// Usage assistance
char* use[] = {
" ",
"Usage: et-nloo_setup -imod par.file [options] [output_prefix] ",
"------------------------------------------------------------------",
"Generate the set of projections/reprojections needed",
"for the Noise-compensated Leave-One-Out resolution analysis.",
"Note: at the moment only imod is supported.",
"Output notation: if a prefix for the output files is not specified, then the",
"                 prefix is derived from the input projection file name.",
"                 The output file names are so obtained:",
"                  {output_prefix}_input.mrc   masked input projections;",
"                  {output_prefix}_full.mrc    masked reprojections from full reconstruction;",
"                  {output_prefix}_miss.mrc    masked reprojections from partial reconstructions;",
" ",
"Actions:",
"-box 1,1,1,50,50,30      Lower corner coords (first three values) and size of the resolution analysis box (pixels)",
"-sphere 100,100,100,50   Coords of center (first three values) and radius of the resolution analysis sphere (pixels)",
"-maxsphere               Perform analysis on maximum centered sphere included in the volume",
"-ibox 10,10,10,25,25,25  Lower corner coords (first three values) and size of the resolution analysis inner box (pixels)",
"-isphere 50,50,50,20     Coords of center (first three values) and radius of the resolution analysis inner sphere (pixels)",
" ",
"Parameters:",
"-verbose 7               Verbosity of output",
"-imod tilt.com           Imod reconstruction code selected, with related input file",
"-background  frame       Region, with respect to the reprojection from the selected volume,",
"                         where calculating the background to apply:",
"                           1=frame (default, frame along the edge of the projection)",
"                           2=inside (all the projected area)",
"-proc 3,25               Distribution of projection generation among multiple processors",
"                         First vaue is the processor number, second value is the total number of processors to use",
"                           +processor #1 needs to be the first to run in order to initialize ",
"						     all the output files",
"                           +input/output files need to be in a directory shared by all the computers",
"                           +total number of processor can not exceed the value set in -fast (if the option is used)",
"-fast 25                 Approximate mode: total number of partial reconstructions from which to",
"                         generate all the needed reprojections (default: a partial reconstruction",
"                         for each input projection, i.e. no approximation)",
"-undersamp 5             Approximate mode: undersampling factor reducing the number",
"                         of projections/reprojections generated (default: 1, i.e. no approximation)",
"-buffer 512              Size of buffer created for loading the tomogram (MB, default 256).",
"                         If the buffer size is greater than the size of the tomogram, then",
"                         all the data are loaded into memory.",
"-no_overwrite            Do not overwrite output projection files if already existent (useful for distributed processing)",
" ",
NULL
};

int 		main(int argc, char **argv)
{

	// Initialize variables

	long long 	buff_size = 256 * 1024 * 1024; // 256 MB
	Buff* 		bfm = NULL, *bf=NULL;
	DataType 	datatype = Float;      		// Conversion to new type
	Vector3 		vol_origin = {0.,0.,0.}; 		// Origin coords of volume
	Vector3 		proj_origin = {0.,0.,0.}; 	// Origin coords of projections
	int			log_tag = 0;					// exponential relationship between intensity and projected density
	float		log_base = 1.;				// base value to add before calculating the logarithm
	char*		out_pfx_pfile = NULL;			// output prefix projection filename
	char*		fullproj_file = NULL;			// projections from full reconstruction filename
	char*		refproj_file = NULL;			// masked input projections filename
	char*		missproj_file = NULL;			// projections from missing reconstruction filename
	char*		flip_vol_file = NULL;			// flipped volume file
	char*		in_vol_file = NULL;			// input reconstruction filename
	char*		in_proj_file = NULL;			// input projection filename
	char*		bin_proj_file = NULL;			// input projection filename (for bsoft read routines)
	char*		imod_pfx = NULL;				// output prefix for imod-generated files
	char*		missvol_pfx_file = NULL;		// partial reconstruction prefix name
	char*		missvol_file = NULL;			// partial reconstruction filename
	char*		tlt_file = NULL;				// tilt angles file
	char*		imod_parfile = NULL;			// imod parameters file

	Bimage*		b = NULL;					// generic image
	Bimage*		bb = NULL;					// test image
	Bimage*		miss_vol = NULL;				// partial reconstruction volume
	Bimage*		full_vol = NULL;				// original reconstruction volume
	Proj*		p = NULL;					// generic projection
	Proj*		miss_proj = NULL;			// reprojection from partial reconstruction
	Proj*		full_proj = NULL;			// reprojection from original reconstruction
	Proj*		ref_proj = NULL;				// input projection


	Imod_tilt*   imd_tpars = NULL;			// imod tilt parameters
	Imod_tilt*   imd_miss = NULL;  			// partial recpnstruction imod tilt parameters

	int			nangs = 0;					// number of tilt angles
	float*		tlt_angs = NULL;				// tilt series
	int			nviews = 0;					// number of tilt views
	View*		tlt_views = NULL;			// tilt series expressed as views
	float		psi0 = 0.0;					// Euler angles
	float		theta0 = 0.0;				//   corresponding to the reference view
	float		phi0 = 0.0;					//   of the ideal planar acquisition
	View			tlt0_view = {0,0,1,0};		// Reference view of ideal planar acquisition
	float		xtilt = 0.0;					// default x tilt angle

	int			np;						// number of projections
	int			taper = 0;				// projection tapering flag

	Shape3d*	volut = NULL;     		// volume under test
	int			vol_type = VALL;		// resolution space constraint 0=all; 1=box
										// 2=sphere ; 3=maxsphere
	int			ivol_type = VNONE;		// resolution inner space constraint 0=none;
										// 1=box ; 2=sphere
	float		ox=-1., oy=-1., oz=-1.;
	float		sx=-1., sy=-1., sz=-1.;
	float		iox=-1., ioy=-1., ioz=-1.;
	float		isx=-1., isy=-1., isz=-1.;

	int			bground_type = SFRAME;	// how to evaluate the background:
										// 1 = frame = along a frame internal to the selected projection area
										// 2 = inside = internal to the selected projection area

	int			tusamp = 1;				// tilt series undersampling factor

	int			proc_id = 1, Nproc = 1;
	int			set_distribute = 0;		// task distribution flag
	int			tfirst, tlast;			// first and last task to locally perform
	int			fast_step = -1;			// fast mode, i.e. number of partial reconstructions
	                                    // to calculate: it also corresponds to the
	                                    // step between reprojections
										// evaluated from the same partial volume
	int			ifirst, ilast, nrun;
	int			tt, it, itt, nl, nl_max;
	int			overwrite_flag = 1;		// ovewrite already existent output projection files

	int			optind;

	Option*		option = get_option_list(use, argc, argv, &optind);
	Option*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
 		if ( strcmp(curropt->tag, "imod") == 0 )
			imod_parfile = get_option_filename(curropt->value);
 		if ( strcmp(curropt->tag, "box") == 0 ) {
			if ( sscanf(curropt->value, "%f,%f,%f,%f,%f,%f", &ox,&oy,&oz,&sx,&sy,&sz) < 6 )
				fprintf(stderr, "-box: lower corner and size (6 values) must be specified");
		    else if (vol_type != VALL)
				fprintf(stderr, "Error: selection of only one type of volume allowed!");
		    else vol_type = VBOX;
 		}
 		if ( strcmp(curropt->tag, "sphere") == 0 ) {
			if ( sscanf(curropt->value, "%f,%f,%f,%f", &ox,&oy,&oz,&sx) < 4 )
				fprintf(stderr, "-sphere: center and radius (4 values) must be specified");
		    else if (vol_type != VALL)
				fprintf(stderr, "Error: selection of only one type of volume allowed!");
		    else vol_type = VSPHERE;
 		}
 		if ( strcmp(curropt->tag, "maxsphere") == 0 ) {
		    if (vol_type != VALL)
				fprintf(stderr, "Error: selection of only one type of volume allowed!");
		    else vol_type = VMAXSPHERE;
 		}
 		if ( strcmp(curropt->tag, "ibox") == 0 ) {
			if ( sscanf(curropt->value, "%f,%f,%f,%f,%f,%f", &iox,&ioy,&ioz,&isx,&isy,&isz) < 6 )
				fprintf(stderr, "-ibox: lower corner and size (6 values) must be specified");
		    else if (ivol_type != VNONE)
				fprintf(stderr, "Error: selection of only one type of inner volume allowed!");
		    else ivol_type = VBOX;
 		}
 		if ( strcmp(curropt->tag, "isphere") == 0 ) {
			if ( sscanf(curropt->value, "%f,%f,%f,%f", &iox,&ioy,&ioz,&isx) < 1 )
				fprintf(stderr, "-isphere: center and radius (4 values) must be specified");
		    else if (ivol_type != VNONE)
				fprintf(stderr, "Error: selection of only one type of inner volume allowed!");
		    else ivol_type = VSPHERE;
 		}
		if ( strcmp(curropt->tag, "background") == 0 ) {
			if ( sscanf(curropt->value, "%d", &bground_type) < 1 ) {
				for ( int i=0; i<(int)strlen(curropt->value); i++ )
					curropt->value[i] = tolower(curropt->value[i]);
				if ( curropt->value[0] == 'f' ) {
					bground_type = SFRAME;
				} else if ( curropt->value[0] == 'i' ) {
					bground_type =SINSIDE;
				}
			}
			if ( bground_type < 1 || bground_type > 2 ) bground_type = SFRAME;
		}
 		if ( strcmp(curropt->tag, "proc") == 0 )
			if ( sscanf(curropt->value, "%d,%d", &proc_id, &Nproc) < 2  )
				fprintf(stderr, "-proc: a pair of values must be specified");
 			else
 				set_distribute = 1;
 				if ( proc_id < 1 || Nproc < 1) {
					fprintf(stderr, "-proc: values less than 1 are not allowed\n");
					exit(-1);
 				}
  		if ( strcmp(curropt->tag, "fast") == 0 )
			if ( sscanf(curropt->value, "%d", &fast_step) < 1 )
				fprintf(stderr, "-fast: a number of reconstructions must be specified");
			else if (fast_step<3) {
				fprintf(stderr, "-fast: values less than 3 are not allowed\n");
				exit(-1);
			}
 		if ( strcmp(curropt->tag, "undersamp") == 0 )
			if ( sscanf(curropt->value, "%d", &tusamp) < 1 )
				fprintf(stderr, "-undersamp: an undersampling factor must be specified");
  		if ( strcmp(curropt->tag, "buffer") == 0 )
			if (sscanf(curropt->value, "%lld", &buff_size) < 1 )
				fprintf(stderr, "-buffer: a valid size of the buffer must be specified");
 			else
 				buff_size *= 1024*1024;
  		if ( strcmp(curropt->tag, "no_overwrite") == 0 )
			overwrite_flag = 0;
    }
	option_kill(option);

	if ( verbose ) {
		printf("-----------------------------------------------\n");
		printf("                   et-nloo_setup\n");
		printf("    Prepare data for NLOO resolution analysis\n");
		printf("-----------------------------------------------\n");
	}

	if ( optind < argc ) out_pfx_pfile = copystring(argv[optind]);

	// check if selected imod
	if ( imod_parfile == NULL ) {
		fprintf(stderr, "Error: imod selection needed!\n");
		exit(-1);
	}

	// only one approximation mode is allowed
	if ( fast_step>0 && tusamp!=1) {
		fprintf(stderr, "Error: can not select both -fast and -undersamp options!\n");
		exit(-1);
	}

/////////////////////////////////////////////////
//		initialization
/////////////////////////////////////////////////

	if ( verbose & VERB_PROCESS )
		printf("Reading imod parameter file %s\n",imod_parfile);

	imd_tpars = imd_tilt_read(imod_parfile);
	if ( imd_tpars == NULL ) exit(-1);
	if ( verbose & VERB_PROCESS ) printf("Parameters read. Checking orientation of tomogram %s\n",imd_tilt_getvolfile(imd_tpars));
	imd_tilt_checkvol(imd_tpars);
	if (!imd_tilt_isparallel(imd_tpars)) {
		if (set_distribute > 0) {
			// distribution of computation requires flipping from outside the program
			fprintf(stderr,"\nERROR: the tomogram has been generated by tilt from IMOD in the PERPENDICULAR mode!\n");
			fprintf(stderr,"       Please flip the tomogram with the following command from IMOD:\n");
			fprintf(stderr,"           clip flipyz [old_file] [new_file]\n");
			fprintf(stderr,"       Then substitute the name of [old_file] with [new_file] in %s\n", imod_parfile);
			exit(-1);
		} else {
			if (verbose && VERB_LABEL) {
				printf("Attention: because of internal requirements,");
				printf(" a copy of %s, with y and z axes flipped, is being created!\n",imd_tilt_getvolfile(imd_tpars));
				printf("This action could take a while, depending on the size of the volume\n");
			}
			flip_vol_file = insert_in_filename(imd_tilt_getvolfile(imd_tpars),"tmp",'_');
			if (imd_clip_flip("yz", imd_tilt_getvolfile(imd_tpars), flip_vol_file) <0) exit(-1);;
			imd_tilt_putvolfile(imd_tpars, flip_vol_file);
			imd_tilt_setparallel(imd_tpars);
			imd_tilt_checkvol(imd_tpars);
		}
	}

	// initialize input and output file names
	in_vol_file = catenate2strings(imd_tilt_getvolfile(imd_tpars),":mrc");
	bin_proj_file = catenate2strings(imd_tilt_getprojfile(imd_tpars),":mrc");
	in_proj_file = imd_tilt_getprojfile(imd_tpars);
	if (out_pfx_pfile == NULL)
		out_pfx_pfile = filename_change_type(in_proj_file,"mrc");
	else if (extension(out_pfx_pfile))
		out_pfx_pfile = copystring(filename_change_type(out_pfx_pfile,"mrc"));
	else
		out_pfx_pfile = catenate2strings(out_pfx_pfile, ".mrc");

	imod_pfx = filename_base(out_pfx_pfile);

	missproj_file = insert_in_filename(out_pfx_pfile,"miss",'_');
	fullproj_file = insert_in_filename(out_pfx_pfile,"full",'_');
	refproj_file = insert_in_filename(out_pfx_pfile,"input",'_');

	missvol_pfx_file = insert_in_filename(filename_change_type(imod_pfx, "mrc"),"miss",'_');

	log_tag = imd_tilt_getlog(imd_tpars);
	if (log_tag) {
		log_base = imd_tilt_getlogbase(imd_tpars);
		b = io_read_image(bin_proj_file, 0, -1);
		if (b->min + log_base < 1. ) {
				if (verbose) {
					printf("Warning: input projections have negative or not correct values\n");
					printf(" changing the base value to add before calculating the logarithm from %f to %f\n",log_base, 1.-b->min);
				}
				log_base = 1. - b->min;
		}
		kill_img(b);
	}

	// projection tilt series
	b = io_read_image(bin_proj_file, 0, -1);

	proj_origin.x = 0.5 * ( b->x - 1.0 );
	proj_origin.y = 0.5 * ( b->y - 1.0 );
	proj_origin.z = 0.;

	kill_img(b);

	// volume under resolution analysis
	b = io_read_image(in_vol_file, 0, -1);

	vol_origin.x = 0.5 * ( b->x - 1.0 );
	vol_origin.y = 0.5 * ( b->y - 1.0 );
	vol_origin.z = 0.5 * ( b->z - 1.0 );
	if (imd_tpars->setshift) {
		vol_origin.x += imd_tpars->shift_x;
		vol_origin.z += imd_tpars->shift_z;
	}

	if ( vol_type == VALL || vol_type == VMAXSPHERE ) {
		volut = shape3d_maxinit( volut, b, vol_type );
	} else
		volut = shape3d_init( volut, vol_type, ox-1, oy-1, oz-1, sx, sy, sz, 0);

	if ( ivol_type != VNONE )
		volut = shape3d_init( volut, ivol_type, iox-1, ioy-1, ioz-1, isx, isy, isz, 1);

	if(shape3d_check(b,volut)<0) exit(-1);

	kill_img(b);

	shape3d_taper_init(volut, 0);

	tlt_file = imd_tilt_gettiltfile(imd_tpars);
	if ( verbose & VERB_PROCESS )
		printf("Reading tilt angles file %s\n",tlt_file);
	tlt_angs = tlt_load_tilt_angles(tlt_file,&nangs);

	if ( verbose & VERB_PROCESS )
		printf("Loaded %d tilt angles\n",nangs);

	// check consistency between angles and projections
	b = io_read_image(bin_proj_file, 0, -1);

	if( nangs != (int) ((b->n==1&&b->z>1)?b->z:b->n)) {
		fprintf(stderr,"Error: number of tilt angles in %s and number of projections in %s are not consistant!\n",
			tlt_file, in_proj_file);
		exit(-1);
	}

	np = nangs;

	// initialize speed-up mechanism
	if (fast_step < 2 || fast_step > np)
		fast_step = np;

	if ( fast_step < np && (verbose & VERB_PROCESS) )
			printf("Distance between projections evaluated from the same partial reconstruction: %d\n",fast_step);


	// initialize distribution of processes
	if ( Nproc > fast_step ) {
		Nproc = fast_step;
		if (verbose)
			if ( fast_step < np )
				printf("Warning: number of total processors reduced to %d (can not exceed the fast factor)\n",Nproc);
			else
				printf("Warning: number of total processors reduced to %d (can not exceed the number of projections)\n",Nproc);
	}

	if ( proc_id > Nproc ) {
		fprintf(stderr,"Warning: run of processor %d aborted because its number is greater than the number of runs (%d) required!\n",proc_id,Nproc);
		exit(0);
	}

	if (set_distribute) {
		div_t ntsk = div(fast_step,Nproc);
		if ( proc_id > (Nproc - ntsk.rem) ) {
			tfirst = (Nproc - ntsk.rem)*ntsk.quot + ( proc_id - Nproc + ntsk.rem -1)*(ntsk.quot+1) + 1;
			tlast = tfirst + ntsk.quot;
		} else {
			tfirst = (proc_id-1)*ntsk.quot + 1;
			tlast = tfirst + ntsk.quot - 1;
		}

	} else {
		tfirst = 1;
		tlast = fast_step;
	}

	ifirst = tfirst-1;
	ilast = tlast-1;


    // max number of projections to evaluate from each partial reconstruction
	nl_max = np / fast_step;
	while ( (np-nl_max*fast_step)>0 ) nl_max++;
	int* lprojs = (int *) balloc(nl_max*sizeof(int));
	nrun = fast_step;


	// generates projections filled with zeros
	// (if task mode is chosen, only the first process is in charge)
	if (tfirst<2) {
		p = proj_init(1);
		proj_cpy_img_sett(p, b);
		p->img[0] = init_img();

		p->datatype = datatype;
		get_cmd_line(p->img[0]->label,argc,argv);
		if ( !overwrite_flag ) bb = io_read_image(missproj_file, 0, -1);
		if ( bb==NULL )
			io_pwrite(missproj_file, p, IO_ZEROS, np);
		kill_img(bb);
		if ( !overwrite_flag ) bb = io_read_image(fullproj_file, 0, -1);
		if ( bb==NULL )
			io_pwrite(fullproj_file, p, IO_ZEROS, np);
		kill_img(bb);
		if ( !overwrite_flag ) bb = io_read_image(refproj_file, 0, -1);
		if ( bb==NULL )
			io_pwrite(refproj_file, p, IO_ZEROS, np);

		proj_kill(p);
	}

	kill_img(b);

	// initialize tilt views
	xtilt = imd_tilt_getxtilt(imd_tpars);
	theta0 = xtilt;
	phi0 = 0.5*M_PI;
	psi0 = -0.5*M_PI;
	tlt0_view = view_from_euler(psi0, theta0, phi0);

	nviews = nangs;

	tlt_views = tlt_views_from_angles(nviews,tlt0_view,tlt_angs);


/////////////////////////////////////////////////
// reprojection from volumes with missing tilt(s)
/////////////////////////////////////////////////
	if ( verbose & VERB_LABEL )
		printf("\n\nStarting generation of reprojections from partial reconstructions\n\n");

	for ( it = ifirst; it<=ilast; it++) {

		if ( it%tusamp ) continue;

		nl = imd_tilt_makelist(it, np, fast_step, imd_tpars, lprojs);

		// skip cycle if all projections not utilized for reconstruction
		if ( nl == 0 ) continue;

		if ( verbose & VERB_FULL )
			printf("Iteration # %d of %d\n",it-ifirst+1, ilast-ifirst+1);

		// generate map volume after excluding one or more tilt angles
		if ( verbose & VERB_PROCESS )
			printf("Launching partial reconstruction\n");

		missvol_file = number_filename(number_filename(missvol_pfx_file,fast_step,3),it+1,3);
		imd_miss = (Imod_tilt *) balloc(sizeof(Imod_tilt));
		imd_tilt_cpy_and_modify(imd_miss, imd_tpars,
				number_filename(number_filename(filename_change_type(imod_pfx, "com"),fast_step,3),it+1,3),
				in_proj_file, missvol_file, -1, -1, -1, 100., (char *) NULL) ;

		for (tt=0; tt<nl; tt++)
			imd_tilt_add2xcldlist(imd_miss,lprojs[tt]+1);

		if ( imd_tilt_writepars(imd_miss) < 0 ) exit(-1);
		if ( imd_tilt_execute(imd_miss) < 0 ) exit(-1);

		bfm = buff_init(buff_size);
		miss_vol = io_read_image_buff((Bimage *) NULL, missvol_file, -1, bfm);
		miss_vol->image[0].background = miss_vol->min;

		if ( verbose & VERB_PROCESS )
			printf("Reprojecting from volume generated without corresponding input projection\n");

		// force origin to the center
		img_set_origin(miss_vol, vol_origin);

		// generate and save projection(s) at the given angle
		for ( tt = 0; tt<nl; tt++) {
			itt = lprojs[tt];
			miss_proj = prjct_generate_one(miss_vol, itt, tlt_views, volut, taper, bground_type, 0, bfm);

			proj_to_datatype(miss_proj, datatype);
			io_pwrite(missproj_file, miss_proj, itt, 0);

			proj_kill(miss_proj);
		}
		kill_img(miss_vol);
		buff_kill(bfm);

		imd_tilt_rmalloutfiles(imd_miss);
		imd_tilt_rmparfile(imd_miss);
		bfree(imd_miss,sizeof(Imod_tilt));
		bfree_string(missvol_file);

	}

/////////////////////////////////////////////////
//	reprojection from full volume
/////////////////////////////////////////////////
	if ( verbose & VERB_LABEL )
		printf("Reading original reconstruction\n");

	bf = buff_init(buff_size);
	full_vol = io_read_image_buff((Bimage *) NULL, in_vol_file, -1, bf);
	full_vol->image[0].background = full_vol->min;

	// force origin to the center
	img_set_origin(full_vol, vol_origin);

	if ( verbose & VERB_LABEL )
		printf("\n\nStarting generation of reprojections from original reconstruction\n\n");

	for ( it = ifirst; it<=ilast; it++) {

		if ( it%tusamp ) continue;

		nl = imd_tilt_makelist(it, np, fast_step, imd_tpars, lprojs);

		// skip cycle if all projections not utilized for reconstruction
		if ( nl == 0 ) continue;

		if ( verbose & VERB_FULL )
			printf("Iteration # %d of %d\n",it-ifirst+1, ilast-ifirst+1);

		// generate and save projection(s) at the given angle
		for ( tt = 0; tt<nl; tt++) {
			itt = lprojs[tt];
			full_proj = prjct_generate_one(full_vol, itt, tlt_views, volut, taper, bground_type, 0, bf);

			proj_to_datatype(full_proj, datatype);
			io_pwrite(fullproj_file, full_proj, itt, 0);

			proj_kill(full_proj);
		}
	}
	kill_img(full_vol);
	buff_kill(bf);

////////////////////////////////////////////////////////
//	Generation of a masked version of input projections
////////////////////////////////////////////////////////

	if ( verbose & VERB_LABEL )
		printf("\n\nStarting generation of masked versions of the input projections\n\n");

	for ( it = ifirst; it<=ilast; it++) {

		if ( it%tusamp ) continue;

		nl = imd_tilt_makelist(it, np, fast_step, imd_tpars, lprojs);

		// skip cycle if all projections not utilized for reconstruction
		if ( nl == 0 ) continue;

		if ( verbose & VERB_FULL )
			printf("Iteration # %d of %d\n",it-ifirst+1, ilast-ifirst+1);

		// read, mask and save projection(s) at the given angle
		for ( tt = 0; tt<nl; tt++) {
			itt = lprojs[tt];
			b = io_read_image(bin_proj_file, 1, itt);
			b->image[0].ox = proj_origin.x;
			b->image[0].oy = proj_origin.y;
			b->image[0].oz = 0.;
			if (log_tag) if (img_op_ln(b, log_base) < 0) exit(-1);
			ref_proj = proj_init(1);
			proj_cpy_img_sett(ref_proj, b);
			proj_putimage(ref_proj, b, 0, 0);
			proj_putview(ref_proj, tlt_views[itt], 0);

			proj_init_footprint(ref_proj, vol_origin, volut, taper);
			proj_stats(ref_proj, bground_type);
			proj_mask(ref_proj);

			proj_to_datatype(ref_proj, datatype);
			io_pwrite(refproj_file, ref_proj, itt, 0);

			proj_kill(ref_proj);
		}
	}


	// the task generating the last useful projection is in charge of
	// updating the statistics in the header of the files
	int ilast_true = nrun-1;
	while(imd_tilt_makelist(ilast_true, np, fast_step, imd_tpars, lprojs)==0)
		ilast_true--;
	if (ilast_true>=ifirst && ilast_true<=ilast) {
		io_update_stats(missproj_file);
		io_update_stats(fullproj_file);
		io_update_stats(refproj_file);
	}

	if (flip_vol_file) {
		remove(flip_vol_file);
		bfree_string(flip_vol_file);
	}

	// deallocation
	bfree(lprojs,nl_max*sizeof(int));
	bfree(tlt_angs,nangs*sizeof(float));
	bfree(tlt_views,nviews*sizeof(View));
	bfree_string(tlt_file);
	shape3d_kill(volut);
	bfree_string(missvol_pfx_file);
	bfree_string(imod_pfx);
	bfree_string(refproj_file);
	bfree_string(fullproj_file);
	bfree_string(missproj_file);
	bfree_string(out_pfx_pfile);
	bfree_string(in_proj_file);
	bfree_string(bin_proj_file);
	bfree_string(in_vol_file);
	bfree(imd_tpars,sizeof(Imod_tilt));
	bfree_string(imod_parfile);


	if ( verbose ) {
		printf("-----------------------------------------------\n");
		printf("       et-nloo_setup concluded successfully!\n");
		printf("-----------------------------------------------\n");
	}
	return(0);
}
