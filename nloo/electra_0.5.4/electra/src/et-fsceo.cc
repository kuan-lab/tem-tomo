/*
	et-fsceo.cc
	Evaluate even/odd FSC
	Author: Giovanni Cardone
	Created: 20050111 	Modified: 20050329
*/

#include "bsoft.h"
#include "ioroutines.h"
#include "shape.h"
#include "wedge.h"
#include "norm.h"
#include "resol.h"

// Usage assistance
char* use[] = {
" ",
"Usage: et-fsceo [options] input1.vol input2.vol ",
"------------------------------------------------------------------",
"Evaluate even/odd Fourier Shell Correlation between two input volumes.",
"Alternatively, FSC is evaluated by setting the -fsc option.",
" ",
"Actions:",
"-taper                   Apply tapering function to borders of volume of interest.",
"-box 1,1,1,50,50,30      Lower corner coords (first three values) and size of the resolution analysis box (pixels)",
"-sphere 100,100,100,50   Coords of center (first three values) and radius of the resolution analysis sphere (pixels)",
"-maxsphere               Perform analysis on maximum centered sphere included in the volume",
"-ibox 10,10,10,25,25,25  Lower corner coords (first three values) and size of the resolution analysis inner box (pixels)",
"-isphere 50,50,50,20     Coords of center (first three values) and radius of the resolution analysis inner sphere (pixels)",
"-rescale  ls             Method for rescaling pairs of volumes:",
"                           0=avgstd (default, same average and standard deviation)",
"                           1=ls (least-squares linear regression)",
" ",
"Parameters:",
"-verbose 7               Verbosity of output",
"-sampling 1.5,1.5,1.5    Sampling (angstrom/pixel, a single value sets all three)",
"-resolution 15           Resolution limit for output (angstrom)",
"-threshold 0.3,0.5,0.7   Thresholds at which resolution is evaluated (up to three values, default 0.3,0.5)",
"-bin 4                   Width of annulus region in fourier space (pixels - default 2)",
"-tiltrange -80,80        Minimum and maximum tilt angle to consider (degrees - default none)", 
"-xtilt 2.4               Rotation around x axis by which reconstruction was rotated (degrees - default 0)", 
"-background  frame       Region, with respect to the selected volume, where calculating the background to apply:",
"                           0=outside (external to the volume)",
"                           1=frame (frame along the edge of the selected volume)",
"                           2=inside (default, all the selected volume)",
"-fsc                     Evaluate FSC instead of FSCeo",
" ",
"Input:",
"-mask mask.pif           Reciprocal space mask (same size as input volumes).",
" ",
"Output:",
"-output  res_plot.txt    resolution curve (text file)",
"-intermediate            save to file the masked volumes ( _mask suffix appended to original names)",
" ",
NULL
};

int		main(int argc, char* argv[])
{
	// Initialize all settings
	char* 		map_file = NULL;			// Compared map
	char*		ref_file = NULL;			// Reference map
	char*		curve_file = NULL;		// Output resolution curve file
	char*		mask_file = NULL;		// Input mask file
	Vector3		sampling = {0,0,0};		// Pixel size
	float		res = 0;
	int			bin_size = 2;			// annulus width for resolution analysis
	float		res_thr[RES_NTHRESHOLDS] = {0.3, 0.5, 0.}; // resolution criteria
	int			nrt = 2;						// number of resolution thresholds
	float*		res_est; 
	
	float		tilt_min = 90.;
	float		tilt_max = 90.;
	float		xtilt = 0.;
	float		xyrot = 0.;
		 
	int			fsc_flag = 0;			// fsc evaluation flag
	int			taper = 0;				// volumes tapering flag
	Shape3d*		volut = NULL;			// volume under test
	int			vol_type = VALL;		    // resolution space constraint 0=all; 1=box
										// 2=sphere ; 3=maxsphere 
	int			ivol_type = VNONE;		// resolution inner space constraint 0=none;
										// 1=box ; 2=sphere 
	float		ox=-1., oy=-1., oz=-1.;
	float		sx=-1., sy=-1., sz=-1.;
	float		iox=-1., ioy=-1., ioz=-1.;
	float		isx=-1., isy=-1., isz=-1.;

	int			intermediate = 0;		// write processed volumes before analysis
	int			rescale_type = 0;		// how to normalize pairs of images:
										// 0 = avgstd = same average and standard deviation
										// 1 = ls = least-squares linear regression
	int			bground_type = 2;		// how to evaluate the background:
										// 0 = outside = external to the selected volume
										// 1 = frame = along a frame internal to the volume
										// 2 = inside = internal to the selected volume
											
	Bimage *     b = NULL;
	Bimage *		pmask = NULL;
	    
	int			optind;
	Option*		option = get_option_list(use, argc, argv, &optind);
	Option*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
 		if ( strcmp(curropt->tag, "intermediate") == 0 )
 			intermediate = 1;
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
 		if ( strcmp(curropt->tag, "fsc") == 0 ) {
       	    fsc_flag = 1;
		}
 		if ( strcmp(curropt->tag, "taper") == 0 ) {
       	    taper = 1;
		}
		if ( strcmp(curropt->tag, "sampling") == 0 )
			sampling = get_option_sampling(curropt->value);
		if ( strcmp(curropt->tag, "mask") == 0 )
			mask_file = get_option_filename(curropt->value);
		if ( strcmp(curropt->tag, "output") == 0 )
			curve_file = get_option_filename(curropt->value);
		if ( strcmp(curropt->tag, "resolution") == 0 )
			if ( sscanf(curropt->value, "%f", &res) < 1 )
				fprintf(stderr, "-resolution: A resolution must be specified!\n");
 		if ( strcmp(curropt->tag, "threshold") == 0 )
			if ( (nrt = sscanf(curropt->value, "%f,%f,%f",&res_thr[0],&res_thr[1],&res_thr[2])) < 1 )
				fprintf(stderr, "-threshold: a list of maximum three thresholds must be specified");
		if ( strcmp(curropt->tag, "bin") == 0 )
			if ( sscanf(curropt->value, "%d", &bin_size) < 1 )
				fprintf(stderr, "-bin: a bin size must be specified");
 		if ( strcmp(curropt->tag, "tiltrange") == 0 )
			if ( sscanf(curropt->value, "%f,%f", &tilt_min,&tilt_max) < 2 )
				fprintf(stderr, "-tiltrange: two tilt angles needed");
 		if ( strcmp(curropt->tag, "xtilt") == 0 )
			if ( sscanf(curropt->value, "%f", &xtilt) < 1 )
				fprintf(stderr, "-xtilt: a rotation around x must be specified");
		if ( strcmp(curropt->tag, "background") == 0 ) {
			if ( sscanf(curropt->value, "%d", &bground_type) < 1 ) {
				for ( int i=0; i<(int)strlen(curropt->value); i++ )
					curropt->value[i] = tolower(curropt->value[i]);
				if ( curropt->value[0] == 'e' ) {
					bground_type = VOUTSIDE;
				} else if ( curropt->value[0] == 'f' ) {
					bground_type = VFRAME;
				} else if ( curropt->value[0] == 'i' ) {
					bground_type =VINSIDE;
				}
			}
			if ( bground_type < 0 || bground_type > 2 ) bground_type = VINSIDE;
		}
		if ( strcmp(curropt->tag, "rescale") == 0 ) {
			if ( sscanf(curropt->value, "%d", &rescale_type) < 1 ) {
				for ( int i=0; i<(int)strlen(curropt->value); i++ )
					curropt->value[i] = tolower(curropt->value[i]);
				if ( curropt->value[0] == 'a' ) {
					rescale_type = 0;
				} else if ( curropt->value[0] == 'l' ) {
					rescale_type = 1;
				} else {
					if ( rescale_type < 0 || rescale_type > 1 ) rescale_type = 0;
				}
			}
		}
    }
	option_kill(option);
    
	if ( verbose ) {
		printf("\n---------------------------------------------------------\n");
		printf("                      et-fsceo\n");
		printf("Resolution analysis by even/odd Fourier Shell Correlation\n");
		printf("---------------------------------------------------------\n\n");
	}

	// check if given input files
	if ( optind > argc - 2 ) {
		fprintf(stderr, "Error: two input files needed!\n");
		exit(-1);
	}
	else {
		ref_file = copystring(argv[optind++]);
		map_file = copystring(argv[optind]);
	}

	// initializations
	
	for (int i = nrt; i<RES_NTHRESHOLDS; i++) res_thr[i]=0.;

	if (fsc_flag && verbose && VERB_LABEL)
		printf("Evaluating Fourier Shell Correlation instead of even/odd Fourier Shell Correlation.\n");

	// missing wedge
	Wedge* mss_wdg = wedge_init(tilt_min, tilt_max, xtilt, xyrot);
	
	// volume under resolution analysis	
	b = io_read_image(ref_file, 0, -1);
	if ( vol_type == VALL || vol_type == VMAXSPHERE ) {
		volut = shape3d_maxinit( volut, b, vol_type );
	} else
		volut = shape3d_init( volut, vol_type, ox-1, oy-1, oz-1, sx, sy, sz, 0);
		
	if ( ivol_type != VNONE )
		volut = shape3d_init( volut, ivol_type, iox-1, ioy-1, ioz-1, isx, isy, isz, 1);  

	if(shape3d_check(b,volut)<0) exit(-1);
	
	kill_img(b);

	if (taper == 1 || bground_type == VFRAME)
		shape3d_taper_init( volut, taper); 

	if (bin_size <1 ) bin_size = 1;

	Bimage* 		pref = io_read_image(ref_file, 1, -1);
	img_to_float(pref);
	Bimage* 		pmap = io_read_image(map_file, 1, -1);
	img_to_float(pmap);

	if ( pmap->z == pref->z && pmap->n == pref->n ) {

		if ( mask_file ) pmask = io_read_image(mask_file, 1, -1);

		if ( res > 0 )
			pmap->resolution = pref->resolution = res;
		if ( sampling.x*sampling.y*sampling.z > 0 ) {
			pmap->ux = pref->ux = sampling.x;
			pmap->uy = pref->uy = sampling.y;
			pmap->uz = pref->uz = sampling.z;
		}


		shape3d_image_stats(pmap, volut, bground_type);
		shape3d_image_stats(pref, volut, bground_type);

		shape3d_mask_image(pref, volut, 0);
		
		if(rescale_type==0)
			norm3d_rescale_to_avg_std(pmap, pref->avg, pref->std, &(pref->image[0].background), volut);
		else if(rescale_type==1)
			norm3d_rescale_linregression(pmap, pref, volut);

		if (taper==1) {
			shape3d_taper_image(pmap, volut);
			shape3d_taper_image(pref, volut);
		}
		
		img_stats(pref);
		img_stats(pmap);
		
		if (intermediate) {
			img_check_param(pref);
			get_cmd_line(pref->label,argc,argv);
			io_write_image(insert_in_filename(ref_file,"masked",'_'),pref);
			img_check_param(pmap);
			get_cmd_line(pmap->label,argc,argv);
			io_write_image(insert_in_filename(map_file,"masked",'_'),pmap);
		}
		
		res_est = res_fsc_estimate(pmap, pref, bin_size, pmask, curve_file, res_thr, mss_wdg, !fsc_flag, (Bimage *) NULL);

		if ( verbose & VERB_PROCESS ) {
			printf("Resolution estimate (A):\nn");
			for (int i = 0; i<nrt; i++) {
				if (fsc_flag)
					printf("\tFSC");
				else
					printf("\tFSCeo");
				printf("(%3.1f)",res_thr[i]);
			}
			printf("\n");
			for (int j = 0; j< (int) pmap->n; j++) {
				printf("%d",j+1);
				int nn = RES_NTHRESHOLDS*j;
				for (int ii = 0; ii<nrt; ii++)
					printf("\t%8.3f",res_est[nn+ii]);
				printf("\n");
			}
		}
		
		if (curve_file) {
			FILE* fpt = fopen(insert_in_filename(curve_file,"estimate",'_'), "w");
			fprintf(fpt,"n");
			for (int i = 0; i<nrt; i++) {
				if (fsc_flag)
					fprintf(fpt,"\tFSC");
				else
					fprintf(fpt,"\tFSCeo");
				fprintf(fpt,"(%3.1f)",res_thr[i]);
			}
			fprintf(fpt,"\n");
			for (int j = 0; j< (int) pmap->n; j++) {
				fprintf(fpt,"%d",j+1);
				int nn = RES_NTHRESHOLDS*j;
				for (int ii = 0; ii<nrt; ii++)
					fprintf(fpt,"\t%8.3f",res_est[nn+ii]);
				fprintf(fpt,"\n");
			}
			fclose(fpt);
		}

		bfree(res_est,pmap->n*RES_NTHRESHOLDS*sizeof(float));

	} else
		fprintf(stderr,"Resolution analysis not performed: data sizes not consistant\n");
	
	if (mask_file) kill_img(pmask);
	kill_img(pref);
	kill_img(pmap);
	
	shape3d_kill(volut);
	wedge_kill(mss_wdg);
	
	bfree_string(map_file);
	bfree_string(ref_file);
	bfree_string(mask_file);

	if ( verbose ) {
		printf("-----------------------------------------------\n");
		printf("        et-fsceo concluded successfully!\n");
		printf("-----------------------------------------------\n");
	}
	return(0);
}

