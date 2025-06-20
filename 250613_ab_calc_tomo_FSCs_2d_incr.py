from FSC import *
import mrcfile
import time
import pandas as pd
import os
from resolution_measure_2D_mrc import *
from glob import glob

####### Edit these params
num_cores = 16
sub_region = [-1, 1890, 1890]
sub_sampling_zxy = [32, 8, 8] # z is slice step
cutout_size = 45 # used by resolution measure 2d 
cube_size = 45 # only used by z_clip (to match legacy 3d fsc FOVs)

# For const increment, e.g. rev1 Fig. 1g
#const_incr = True
const_incr = False
num_angs = [41,35,31,25,21,15,11]
max_angs = [60,51,45,36,30,21,15]

# For 2d array, e.g. rev1 Fig 1e,f
num_angs = [121, 33, 21, 17, 11, 5]
max_angs = [60,50,40,30,20,10]
 
plane = 'beam'
output_dir = 'results/250617_11k_2D_incr_abHalfBit315_45pix_4x'
#slice_step = 16
#fake = True
fake = False
overwrite = False
#overwrite = True
###########

# Working with file structure to analyze multiple datasets
home_dir = '/home/atk13/repos/tem-tomo'
#tomo_lst = 'tomo_lists/tomograms_lst - double_tilt_tomos_3.3k.csv'
tomo_lst = 'tomo_lists/tomograms_lst - double_tilt_tomos_11k.csv'
df = pd.read_csv(tomo_lst)

for index,row in df.iterrows():
    proj = 'microscopy_%i' % int(row['MPID'])
    tomo = row['Tomogram']
    thickness = row['Thickness']
    pixel_size = row['Pixel Size bin 4 (nm)']
    ccdbprod  = row['ccdbprod']
    pid = row['PID']
    mag = row['Magnification']
    z_min = int(row['z_min'])
    z_max = int(row['z_max'])

    # calculate z clips to avoid top and bottom artifacts
    z_center = math.floor((z_max + z_min) / 2)
    if mag == '3.3k':
        if thickness == "250nm":
            pixel_z_clip = cube_size*1 
        elif thickness == '500nm':
            pixel_z_clip = cube_size*2	
        elif thickness == '750nm':
            pixel_z_clip = cube_size*3
        elif thickness == '1000nm':
            pixel_z_clip = cube_size*4
    elif mag == '11k':
        if thickness == "250nm":
                pixel_z_clip = cube_size*3
        elif thickness == '500nm':
                pixel_z_clip = cube_size*6
        elif thickness == '750nm':
                pixel_z_clip = cube_size*9
        elif thickness == '1000nm':
                pixel_z_clip = cube_size*12
    z_clip = (math.floor(z_center - pixel_z_clip/2)+1, math.floor(z_center + pixel_z_clip/2))
    print('z_clip : %i, %i' % z_clip)
    data_path = '/ccdbprod/%s/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_%s/' % (ccdbprod, pid)
    tomo_path = os.sep.join([data_path, proj, 'processed_data',tomo,'txbr-backprojection','limited-bin4'])
		
	# Match equal incrmements
    for idx,num_ang in enumerate(num_angs):
        if const_incr:
            loop_max_angs = [max_angs[idx]]
        else: 
            loop_max_angs = max_angs
        for max_ang in loop_max_angs:
            if num_ang == 121 and max_ang == 10: # we don't have these data
                continue

            recon_dir = os.sep.join([tomo_path,'%i-limited[%.1f_-%.1f]' % (num_ang,max_ang,max_ang)])
            os.chdir(recon_dir)
            
            a_paths = glob(tomo+'a_z_-*0.out')
            if len(a_paths) != 1:
                print('Problem dir for a recon: %s' % recon_dir)
                print(a_path)
                fake = True
            else:
                a_path = os.sep.join([recon_dir, a_paths[0]])
				
            b_paths = glob(tomo+'b_z_-*0.out')
            if len(b_paths) != 1:
                print('Problem dir for b recon : %s' % recon_dir)
                print(b_path)
                fake = True
            else:
                b_path = os.sep.join([recon_dir, b_paths[0]])
				
            ofn = os.sep.join([output_dir, 'FSC2D_%s_%s_%i-limited[%.1f_-%.1f].csv' % (thickness, tomo,num_ang,max_ang,max_ang)])
            os.chdir(home_dir)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            print('Calculating FSC for %s' % ofn)
            if not fake:
                if overwrite or not os.path.isfile(ofn):		
                    resolution_measure_2D(a_path, b_path, num_cores, cutout_size = cutout_size, snrt = 0.2071, \
                    pixel_size = pixel_size, z_clip = z_clip, sub_region = sub_region, \
                    sub_sampling_zxy = sub_sampling_zxy, plane = plane, ofn=ofn)
                    # half biit threshold for ab
