from FSC import *
import mrcfile
import time
import pandas as pd
import os
from resolution_measure_2D_mrc import *
from glob import glob

####### Edit these params
num_cores = 24
cutout_size = -1
sub_region = -1
num_angs = [121, 33, 21, 5]
max_angs = [60, 40, 20]
plane = 'beam'
output_dir = '240207_dtFSC2D_' + plane
slice_step = 2
#fake = True
fake = False
#overwrite = False
overwrite = True
###########

# Working with file structure to analyze multiple datasets
home_dir = '/home/atk13/repos/tem-tomo'
#data_path = '/Users/atk42/OneDrive - Yale University/Lab/Projects/TEM_tomo/tomo_data'
#data_path = '/home/atk13/new_project_20471'
data_path = '/ccdbprod/ccdbprod29/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_20471/'
#tomo_lst = 'tomograms_lst - Local Tomograms for FSC.csv'
tomo_lst = 'tomograms_lst - double_tilt_tomos.csv'

df = pd.read_csv(tomo_lst)

for index,row in df.iterrows():
	proj = 'microscopy_%i' % int(row['MPID'])
	tomo = row['Tomogram']
	thickness = row['Thickness']
	pixel_size = row['Pixel Size bin 4 (nm)']

	tomo_path = os.sep.join([data_path, proj, 'processed_data',tomo,'txbr-backprojection','limited-bin4'])
	for num_ang in num_angs:
		for max_ang in max_angs:
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
			ofn = os.sep.join([output_dir, 'dtFSC2D_%s_%s_%i-limited[%.1f_-%.1f].csv' % (thickness, tomo,num_ang,max_ang,max_ang)])
			os.chdir(home_dir)

			if not os.path.exists(output_dir):
				os.makedirs(output_dir)
			print('Calculating FSC for %s' % ofn)
			if not fake:
				if overwrite or not os.path.isfile(ofn):		
					resolution_measure_2D(a_path, b_path, num_cores, cutout_size = cutout_size, \
					pixel_size = pixel_size, sub_region = sub_region, plane = plane, \
					slice_step = slice_step, ofn=ofn)
