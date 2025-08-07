/*
	et-findbeads.cc
	Determine positions of beads in a tomogram
	Author: Giovanni Cardone
	Created: 20050123	Modified: 20070716
*/

#include "bsoft.h"
#include "ioroutines.h"
#include "imod.h"
#include "resol.h"
#include "project.h"

// Usage assistance
char* use[] = {
" ",
"Usage: et-findbeads [options] [output_prefix] ",
"-------------------------------------------------------------------------------------",
"Locate gold beads in the tomogram and remove them, if requested.",
"If removal of beads is selected, a prefix for the output file can be specified.",
"The coordinates corresponding to the beads identified can be printed out by setting '-verbose result'.",
" ",
"Actions:",
"-removefrom projs        remove located beads",
"                           0=none  (do not remove beads)",
"                           1=projs (default, remove from projections)",
"-paint random            paint gold particles",
"                           1=random (assign random value, default)",
"                           2=constant  (default, assign constant value obtained from the adjacent pixels)",
" ",
"Parameters:",
"-verbose 7               Verbosity of output",
"-map vol.map             reconstructed tomogram",
"-projections proj.mrc    input projections",
"-tilt     series.tlt     tilt angles",
"-xtilt 2.4               Rotation around x axis by which reconstruction was rotated (degrees - default 0)",
"-imod tilt.com           Imod reconstruction parameter file (alternative to the previous four parameters)",
"-diameter 20             gold particles diameter (nm, default 10)",
"-ccfactor 3.5            coefficient for determining the cross-correlation peaks",
"                           (threshold calculated as cc_average+ccfactor*cc_standard_deviation, default 10)",
"-ccmap cc.mrc            load a cross-correlation map previously obtained (skip the estimation process)",
" ",
"Output:",
"-save_ccmap cc.mrc       Save cross-correlation map",
" ",
NULL
};

int	ring_image_stats(Bimage* p, VectorInt3 pc, float iradius, float oradius, float* avg, float* std);
int	ring_image_randomize(Bimage* p, VectorInt3 pc, float iradius, float oradius, float avg, float std);

int 		main(int argc, char **argv)
{

	// Initialize variables
	
	char*		in_vol_file = NULL;			// tomogram filename
	char*		in_proj_file = NULL;		// projections filename
	char*		out_pfx_file = NULL;		// output prefix filename
	char*		tlt_file = NULL;			// tilt angles file
	char*		ccmap_file = NULL;			// cross-correlation map file

	Bimage*		b = NULL;				    // input volume
	Bimage*		beadvol = NULL;			    // volume containing a bead template
	Bimage*		p = NULL;				    // single image
	Proj*		pp = NULL;				    // single projection

	char*		imod_parfile = NULL;		// imod parameters file
	Imod_tilt*  imd_tpars = NULL;			// imod tilt parameters

	int 		paint_flag = 2;             // paint gold beads
											// 1 = with random value
											// 2 = with constant value
	int 		load_ccmap = 0;				// load cc map previously saved
	int 		remove_flag = 1;
	float		gdiameter = 10.;			// gold beads diameter (nm)

	int			icheck = 0;				// input consistency check
	int			set_xtilt = 0;

	int			bin_factor = 2;
	int			bin_set = 0;

	Vector3 	vol_origin = {0.,0.,0.}; 	// Origin coords of volume
	Vector3		pcvol, pcproj;
	VectorInt3	ipcproj;

	Wedge*		mss_wdg = NULL;
	int			nangs = 0;					// number of tilt angles
	float*		tlt_angs = NULL;			// tilt series
	float		tmin, tmax;
	int			nviews = 0;					// number of tilt views
	View*		tlt_views = NULL;			// tilt series expressed as views
	float		psi0 = 0.0;					// Euler angles
	float		theta0 = 0.0;				//   corresponding to the reference view
	float		phi0 = 0.0;					//   of the ideal planar acquisition
	View		tlt0_view = {0,0,1,0};		// Reference view of ideal planar acquisition
	float		xtilt = 0.0;				// default x tilt angle
	float		cc_factor = 10.0;			// multiplication factor for determining cc threshold
	
	int			optind;

	Option*		option = get_option_list(use, argc, argv, &optind);
	Option*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
 		if ( strcmp(curropt->tag, "map") == 0 ) {
			in_vol_file = get_option_filename(curropt->value);
 			icheck++;
 		}	
 		if ( strcmp(curropt->tag, "projections") == 0 ) {
			in_proj_file = get_option_filename(curropt->value);
 			icheck++;
 		}	
 		if ( strcmp(curropt->tag, "tilt") == 0 ) {
			tlt_file = get_option_filename(curropt->value);
 			icheck++;
 		}
 		if ( strcmp(curropt->tag, "xtilt") == 0 )
			if ( sscanf(curropt->value, "%f", &xtilt) < 1 )
				fprintf(stderr, "-xtilt: a rotation around x must be specified");
 			else {
 				xtilt = xtilt*M_PI/180.0;
 				set_xtilt = 1;
 			}
 		if ( strcmp(curropt->tag, "imod") == 0 )
			imod_parfile = get_option_filename(curropt->value);
 		if ( strcmp(curropt->tag, "removefrom") == 0 ) {
			if ( sscanf(curropt->value, "%d", &remove_flag) < 1 ) {
				for ( int i=0; i<(int)strlen(curropt->value); i++ )
					curropt->value[i] = tolower(curropt->value[i]);
				if ( curropt->value[0] == 'n' ) {
					remove_flag = 0;
				} else if ( curropt->value[0] == 'p' ) {
					remove_flag = 1;
				}
			}
			if ( remove_flag < 0 || remove_flag > 1 ) remove_flag = 1;
		}
 		if ( strcmp(curropt->tag, "paint") == 0 ) {
			if ( sscanf(curropt->value, "%d", &paint_flag) < 1 ) {
				for ( int i=0; i<(int)strlen(curropt->value); i++ )
					curropt->value[i] = tolower(curropt->value[i]);
				if ( curropt->value[0] == 'r' ) {
					paint_flag = 1;
				} else if ( curropt->value[0] == 'c' ) {
					paint_flag = 2;
				}
			}
			if ( paint_flag < 0 || paint_flag > 2 ) paint_flag = 2;
		}
		if ( strcmp(curropt->tag, "diameter") == 0 )
			if ( sscanf(curropt->value, "%f", &gdiameter) < 1 )
				fprintf(stderr, "-diameter: diameter of gold beads must be specified");
		if ( strcmp(curropt->tag, "ccfactor") == 0 )
			if ( sscanf(curropt->value, "%f", &cc_factor) < 1 )
				fprintf(stderr, "-ccfactor: a crosscorrelation factor must be specified");
 		if ( strcmp(curropt->tag, "save_ccmap") == 0 )
			ccmap_file = get_option_filename(curropt->value);
 		if ( strcmp(curropt->tag, "ccmap") == 0 ) {
			ccmap_file = get_option_filename(curropt->value);
			load_ccmap = 1;
 		}
      }
     
	option_kill(option);

	// check if given input files
	if ( optind == argc - 1 ) {
		out_pfx_file = copystring(argv[optind]);
	}

	if ( verbose & VERB_LABEL) {
		printf("-----------------------------------------------------\n");
		printf("                      et-findbeads\n");
		printf("Determination of position of gold beads in a tomogram\n");
		printf("-----------------------------------------------------\n");
	}

	// check if selected imod or all parameters given
	if ( imod_parfile == NULL && icheck!=3) {
		fprintf(stderr, "ERROR: some of the following mandatory input parameters are missing:\n");
		fprintf(stderr, "   map, projections, tilt.\n");
		exit(-1);
	} else if (imod_parfile != NULL && ( icheck>0 || set_xtilt)) {
		fprintf(stderr, "ERROR: some input parameters are overlapping with those read from the imod file!\n");
		fprintf(stderr, "Parameters read from imod file: map, projections, tilt, xtilt\n");
		exit(-1);
	}		

/////////////////////////////////////////////////
//		initialization
/////////////////////////////////////////////////

	gdiameter *= 10.;
	
	if ( imod_parfile ) {
		if ( verbose & VERB_PROCESS )
			printf("Reading imod parameter file %s\n",imod_parfile);

		imd_tpars = imd_tilt_read(imod_parfile);
		if ( imd_tpars == NULL ) exit(-1);
		imd_tilt_checkvol(imd_tpars);
		if (!imd_tilt_isparallel(imd_tpars)) {
		// distribution of computation requires flipping from outside the program
			fprintf(stderr,"\nERROR: the tomogram has been generated by tilt from IMOD in the PERPENDICULAR mode!\n");
			fprintf(stderr,"       Please flip the tomogram with the following command from IMOD:\n");
			fprintf(stderr,"           clip flipyz [old_file] [new_file]\n");
			fprintf(stderr,"       Then substitute the name of [old_file] with [new_file] in %s\n", imod_parfile);
			exit(-1);
		}

	// initialize input and output file names, and other parameters
		in_vol_file = imd_tilt_getvolfile(imd_tpars);
		in_proj_file = imd_tilt_getprojfile(imd_tpars);

		tlt_file = imd_tilt_gettiltfile(imd_tpars);
		xtilt = imd_tilt_getxtilt(imd_tpars);
	}

    // load tilt angles
	tlt_angs = tlt_load_tilt_angles(tlt_file,&nangs);
	
	// generate wedge mask object
	if ( imod_parfile ) {
		mss_wdg = wedge_init_from_imod(imd_tpars);
	} else {
		tmin = 100.;
		tmax = -100.;
		for (int j=0; j<nangs; j++) {
			if (tlt_angs[j]<tmin) tmin = tlt_angs[j];
			if (tlt_angs[j]>tmax) tmax = tlt_angs[j];
		}
		mss_wdg = wedge_init(tmin*180./M_PI, tmax*180./M_PI, xtilt, 0.0);
	}	
	
	if (remove_flag) {
		if (out_pfx_file == NULL) {
			if (remove_flag == 1)			
				out_pfx_file = copystring(in_proj_file);
			out_pfx_file = insert_in_filename(out_pfx_file,"nobeads",'_');
		}
		if (extension(out_pfx_file))
			out_pfx_file = copystring(filename_change_type(out_pfx_file,"mrc")); 
		else 
			out_pfx_file = catenate2strings(out_pfx_file,".mrc"); 
	}

	in_vol_file = catenate2strings(in_vol_file,":mrc");
	in_proj_file = catenate2strings(in_proj_file,":mrc");


	if (!load_ccmap) {
	// volume to be filtered
		b = io_read_image(in_vol_file, 1, -1);
		img_to_float(b);
	} else  // only header needed
		b = read_img(in_vol_file, 0, -1);

	if ( b->transform != NoTransform ) {
		fprintf(stderr, "Error: File %s must be a real space map!\n", b->filename); 
		exit(-1); 
	}
	
	vol_origin.x = 0.5 * ( b->x - 1.0 );
	vol_origin.y = 0.5 * ( b->y - 1.0 );
	vol_origin.z = 0.5 * ( b->z - 1.0 );

	if (gdiameter/b->ux > 5.) {
		bin_set = 1;
		bin_factor = ((1 + (int) (gdiameter/(5.*b->ux)) ) / 2) * 2;
		VectorInt3	bin = {bin_factor,bin_factor,bin_factor};
		img_bin(b,bin);
	}

//define the gold bead as a sphere
//	int gx = 2*(int)(gdiameter/b->ux);
//	int gy = 2*(int)(gdiameter/b->uy);
//	int gz = 2*(int)(gdiameter/b->uz);
//	float fgx = gx*b->ux;
//	float fgy = gy*b->uy;
//	float fgz = gz*b->uz;
	int gdx = (int)(gdiameter/b->ux);
	int gdy = (int)(gdiameter/b->uy);
	int gdz = (int)(gdiameter/b->uz);
	
	unsigned long	i, n;
	int    x, y, z, xx, yy, zz;
	float  fx, fy, fz;

	if (!load_ccmap) {
		beadvol = copy_img(b);
		img_clear_data(beadvol);
		img_to_float(beadvol);
	
		img_fft(FFTW_FORWARD, b);
		img_stats(b);
	
		if ( verbose & VERB_DEBUG ) {
			printf("DEBUG: FFT done\n");
			printf("DEBUG: Average = %g, StDev = %g\n\n", b->avg, b->std);
		}


		float* data = (float *) beadvol->data;

		for ( n=0; n<b->n; n++) {
			for ( z=-gdz; z<gdz; z++ ) {
				fz = (float) z * b->uz;
				zz = z;
				if (zz<0) zz+=b->z; 
				for ( y=-gdy; y< gdy; y++ ) {
					fy = (float) y * b->uy;
					yy = y;
					if (yy<0) yy+=b->y; 
					for ( x=-gdx; x< gdx; x++ ) {
						fx = (float) x * b->ux;
						xx = x;
						if (xx<0) xx+=b->x; 
					
						if ( fx*fx+fy*fy+fz*fz <= 0.25*gdiameter*gdiameter ) {
							i = ((n*b->z + zz)*b->y + yy)*b->x + xx;
							data[i] = -1.;
						}
					} 
				}
			}
		}	

	
		img_fft(FFTW_FORWARD, beadvol);
		img_stats(beadvol);

		wedge_mask_fourier_vol(beadvol,mss_wdg);

		img_complex_product(beadvol, b);

		wedge_kill(mss_wdg);
		kill_img(b);

		img_fft(FFTW_BACKWARD, beadvol);
		img_stats(beadvol);

		if (ccmap_file != NULL) {
			ccmap_file = copystring(filename_change_type(ccmap_file, "mrc"));
			get_cmd_line(beadvol->label, argc, argv);
			write_img(ccmap_file, beadvol);		
		}

	} else {
		beadvol = io_read_image(ccmap_file, 1, -1);
	}

	unsigned long ncoord = 0;
	double threshold = beadvol->avg+cc_factor*beadvol->std;
	VectorInt3* pcoords = img_find_particles_in_ccmap(beadvol, gdiameter/beadvol->ux, 
					&threshold, &ncoord);

	if (verbose & VERB_FULL)
		printf("Number of gold beads found: %ld (threshold=%lf)\n",ncoord,threshold);
			
	if (bin_set==1) {
		for (n=0; n<ncoord; n++) {		
			pcoords[n].x = pcoords[n].x*bin_factor + (bin_factor-1)/2;
			pcoords[n].y = pcoords[n].y*bin_factor + (bin_factor-1)/2;
			pcoords[n].z = pcoords[n].z*bin_factor + (bin_factor-1)/2;
		}
	}

	if (verbose & VERB_RESULT) {
		printf("Beads coordinates:\n");
		for (n=0; n<ncoord; n++) {		
			printf("%ld %d %d %d\n",n+1,pcoords[n].x,pcoords[n].y,pcoords[n].z);
		}
	}
	
	kill_img(beadvol);

	if (remove_flag==1 && ncoord>0 ) {


		// initialize tilt views
		theta0 = xtilt;//*M_PI/180.0;
		phi0 = 0.5*M_PI;
		psi0 = -0.5*M_PI;
		tlt0_view = view_from_euler(psi0, theta0, phi0);

		nviews = nangs;

		tlt_views = tlt_views_from_angles(nviews,tlt0_view,tlt_angs);

		p = io_read_image(in_proj_file, 0, -1);
		get_cmd_line(p->label, argc, argv);
		pp = proj_init(1);
		proj_cpy_img_sett(pp, p);
	
		io_pwrite(out_pfx_file, pp, IO_ZEROS, nviews);

		proj_kill(pp);
		kill_img(p);

		float ravg, rstd;
		for ( int it=0; it<nviews; it++) {

			if ( verbose & VERB_FULL )
				printf("Adjust projection # %d of %d\n",it+1, nviews);

			p = io_read_image(in_proj_file, 1, it);
			pp = proj_init(1);
			proj_cpy_img_sett(pp, p);
			proj_putimage(pp, p, 0, 0);

			for ( int ig=0; ig<(int)ncoord; ig++) {
				pcvol = vector3_from_vectorint3(pcoords[ig]);
				pcproj = prjct_locate_point(pcvol, vol_origin, tlt_views[it]);
				ipcproj = vectorint3_from_vector3(pcproj);
				ring_image_stats(p, ipcproj, 0.5*gdiameter/p->ux, gdiameter/p->ux, &ravg, &rstd);
				if (paint_flag == 2) rstd = 0.;
				ring_image_randomize(p, ipcproj, 0., 1.2*(0.5*gdiameter/p->ux), ravg, rstd);
			}
			io_pwrite(out_pfx_file, pp, it, 0);

			proj_kill(pp);
		}

		io_update_stats(out_pfx_file);
	}
	

	bfree(tlt_angs,nangs*sizeof(float));	
	bfree(tlt_views,nviews*sizeof(View));	
	bfree_string(tlt_file);
	bfree(pcoords, ncoord*sizeof(VectorInt3));
	bfree_string(imod_parfile);
	bfree(imd_tpars,sizeof(Imod_tilt));	
	if(out_pfx_file)	bfree_string(out_pfx_file);
	bfree_string(in_vol_file);
	bfree_string(in_proj_file);

	if ( verbose  & VERB_LABEL) {
		printf("-----------------------------------------------\n");
		printf("       et-findbeads concluded successfully!\n");
		printf("-----------------------------------------------\n");
	}
	return(0);
}

int	ring_image_stats(Bimage* p, VectorInt3 pc, float iradius, float oradius, float* avg, float* std)
{
	if ( !p->data ) return(-1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short*	usdata = (unsigned short *) p->data;
    short* 		    	sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    		fdata = (float *) p->data;
    complex_short* 	csdata = (complex_short *) p->data;
    complex_int* 	cidata = (complex_int *) p->data;
    complex_float*	cfdata = (complex_float *) p->data;
    polar*  			pdata = (polar *) p->data;
    
	double		sum, sum2;
	float		amp, dist2;
    long     	i, number;
	unsigned int x, y, z;
	float oradius2 = oradius*oradius;
	float iradius2 = iradius*iradius;

	sum = sum2 = 0.;
	number = 0;
	
	for ( unsigned int n = 0; n < p->n; n++) {
	
		for ( z=0; z<p->z; z++ ) {
			if (fabs((float) pc.z-z)>oradius) continue;
			for ( y=0; y<p->y; y++ ) {
				if (fabs((float) pc.y-y)>oradius) continue;
				for ( x=0; x<p->x; x++ ) {
					if (fabs((float) pc.x-x)>oradius) continue;
					i = ((n*p->z+z)*p->y+y)*p->x+x;
					dist2 = (pc.x-x)*(pc.x-x)+(pc.y-y)*(pc.y-y)+(pc.z-z)*(pc.z-z);
					// average and standard deviation
					if ( dist2<oradius2 && dist2>=iradius2 ) {
						switch ( p->datatype ) {
							case UChar:
								sum += udata[i];
								sum2 += udata[i]*1.0*udata[i];
								break;
							case SChar:
								sum += cdata[i];
								sum2 += cdata[i]*1.0*cdata[i];
								break;
							case UShort:
								sum += usdata[i];
								sum2 += usdata[i]*1.0*usdata[i];
								break;
							case Short:
								sum += sdata[i];
								sum2 += sdata[i]*1.0*sdata[i];
								break;
							case Int:
								sum += idata[i];
								sum2 += idata[i]*1.0*idata[i];
								break;
							case Float:
								sum += fdata[i];
								sum2 += fdata[i]*fdata[i];
								break;
							case ComplexShort:
								amp = sqrt(1.0*csdata[i].re*csdata[i].re +
    	    								csdata[i].im*csdata[i].im);
								sum += amp;
								sum2 += amp*amp;
								break;
							case ComplexInt:
								amp = sqrt(1.0*cidata[i].re*cidata[i].re +
    	    								cidata[i].im*cidata[i].im);
								sum += amp;
								sum2 += amp*amp;
								break;
							case ComplexFloat:
								amp = sqrt(cfdata[i].re*cfdata[i].re +
    	    								cfdata[i].im*cfdata[i].im);
								sum += amp;
								sum2 += amp*amp;
								break;
							case Polar:
								sum += pdata[i].amp;
								sum2 += pdata[i].amp*pdata[i].amp;
								break;
							default: break;
						}
						number++;
					}
				}
			}
		}
	}
	if (number) {
		*avg = sum/number;
		*std = sum2/number - *avg * *avg;
		if ( *std > 0 ) *std = sqrt(*std);
		else *std = 0.;
	}

	return(0);
}


int	ring_image_randomize(Bimage* p, VectorInt3 pc, float iradius, float oradius, float avg, float std)
{
	if ( !p->data ) return(-1);
	
    unsigned char* 	udata = (unsigned char *) p->data;
    signed char* 	cdata = (signed char *) p->data;
    unsigned short*	usdata = (unsigned short *) p->data;
    short* 		    	sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    		fdata = (float *) p->data;
    complex_short* 	csdata = (complex_short *) p->data;
    complex_int* 	cidata = (complex_int *) p->data;
    complex_float*	cfdata = (complex_float *) p->data;
    polar*  			pdata = (polar *) p->data;
    
	float		rval, dist2;
    long     	i, number;
	long         x, y, z;
	float oradius2 = oradius*oradius;
	float iradius2 = iradius*iradius;

	
	for ( unsigned int n = 0; n < p->n; n++) {
		for ( z=0; z<p->z; z++ ) {
			if (fabs((float) pc.z-z)>oradius) continue;
			for ( y=0; y<p->y; y++ ) {
//printf("* %d - %d -> %f vs %f*",pc.y,y,fabs(pc.y-y), oradius);
				if (fabs((float) pc.y-y)>oradius) continue;
				for ( x=0; x<p->x; x++ ) {
					if (fabs((float) pc.x-x)>oradius) continue;
					dist2 = (pc.x-x)*(pc.x-x)+(pc.y-y)*(pc.y-y)+(pc.z-z)*(pc.z-z);
					// average and standard deviation
					if ( dist2<oradius2 && dist2>=iradius2 ) {
//					if ( z==pc.z && y==pc.y && x==pc.x) {
						i = ((n*p->z+z)*p->y+y)*p->x+x;
						rval = random_gaussian(avg, std);
//						rval = 500.;
						switch ( p->datatype ) {
							case UChar:
								udata[i] = (unsigned char) rval;
								break;
							case SChar:
								cdata[i] = (signed char) rval;
								break;
							case UShort:
								usdata[i] = (unsigned short) rval;;
								break;
							case Short:
								sdata[i] = (short) rval;
								break;
							case Int:
								idata[i] = (int) rval;
								break;
							case Float:
								fdata[i] = rval;
								break;
							case ComplexShort:
								csdata[i].re = (short) rval;
								csdata[i].im = (short) random_gaussian(avg, std);
								break;
							case ComplexInt:
								cidata[i].re = (int) rval;
								cidata[i].im = (int) random_gaussian(avg, std);
								break;
							case ComplexFloat:
								cfdata[i].re = rval;
								cfdata[i].im = random_gaussian(avg, std);
								break;
							case Polar:
								pdata[i].amp = rval;
								break;
							default: break;
						}
						number++;
					}
				}
			}
		}
	}
	return(0);
}
