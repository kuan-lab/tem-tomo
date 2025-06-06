from FSC import *
import mrcfile
import time
import pandas as pd
import os
from resolution_measure_2D_mrc import *

####### Edit these params
num_cores = 24
cutout_size = -1
sub_region = -1
num_angs = [33, 17, 5]
max_angs = [60, 40, 20, 10]
plane = 'beam'
output_dir = '240126_FSC2D_' + plane
slice_step = 2
fake = True
#fake = False
#overwrite = False
overwrite = True
###########

# Working with file structure to analyze multiple datasets

#data_path = '/Users/atk42/OneDrive - Yale University/Lab/Projects/TEM_tomo/tomo_data'
data_path = '/home/atk13/new_project_20471'
tomo_lst = 'tomograms_lst - Local Tomograms for FSC.csv'

df = pd.read_csv(tomo_lst)

for index,row in df.iterrows():
	proj = 'microscopy_%i' % int(row['MPID'])
	tomo = row['Tomogram']
	thickness = row['Thickness']
	pixel_size = row['Pixel Size bin 4 (nm)']

	tomo_path = os.sep.join([data_path, proj, 'processed_data',tomo,'txbr-backprojection','limited-bin4'])
	for num_ang in num_angs:
		for max_ang in max_angs:
			a_dir = os.sep.join([tomo_path,'%i-limited[%.1f_-%.1f]_fsc-a' % (num_ang,max_ang,max_ang)])
			if len(os.listdir(a_dir)) == 1:
    				a_path = os.sep.join([a_dir,os.listdir(a_dir)[0]])
			else:  				
				print('Problem dir: %s' % os.listdir(a_dir))
				fake = True
				# sys.exit('tomo does not have exactly 1 output file')	
			b_dir = os.sep.join([tomo_path,'%i-limited[%.1f_-%.1f]_fsc-b' % (num_ang,max_ang,max_ang)])
			if len(os.listdir(b_dir)) == 1:
    				b_path = os.sep.join([b_dir,os.listdir(a_dir)[0]])
			else:
				print('Problem dir: %s' % os.listdir(b_dir))
				#sys.exit('tomo dir does not have exactly 1 output file')
				fake = True
			ofn = os.sep.join([output_dir, 'FSC3D_%s_%s_%i-limited[%.1f_-%.1f].csv' % (thickness, tomo,num_ang,max_ang,max_ang)])
			if not os.path.exists(output_dir):
				os.makedirs(output_dir)
			print('Calculating FSC for %s' % ofn)
			if not fake:
				if overwrite or not os.path.isfile(ofn):		
					resolution_measure_2D(a_path, b_path, num_cores, cutout_size = cutout_size, \
					pixel_size = pixel_size, sub_region = sub_region, plane = plane, \
					slice_step = slice_step, ofn=ofn)
