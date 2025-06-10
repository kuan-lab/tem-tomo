#!/usr/bin/env python
# 240125 modified ATK to calculate 2D FSC

from __future__ import division, print_function
#import zarr
from scipy import *
import numpy as np
import os, sys
import multiprocessing
#import h5py
import tqdm
import json
from scipy.ndimage import fourier_shift
from scipy.special import erf
from time import time
from FSC import *
import mrcfile
import argparse

# this will be the function the worker threads run to analyze one block
def parallel_FSC_worker(a):
    with mrcfile.open(a['fn1']) as mrc:
        v1 = mrc.data
    with mrcfile.open(a['fn2']) as mrc:
        v2 = mrc.data
    plane = a['plane']
    s = a['cutout_size']
    z,x,y = a['top_left']

    if s > 0:
        if plane == 'beam': # This is plane with the full size images
            a['center'] = (z,x+s//2,y+s//2)
            v1c=v1[z,x:x+s,y:y+s]
            v2c=v2[z,x:x+s,y:y+s]
        elif plane == 'tilt': # This is the plane that will have missing wedge artifacts
            a['center'] = (z+s//2,x,y+s//2)
            v1c=v1[z:z+s,x,y:y+s]
            v2c=v2[z:z+s,x,y:y+s]
        elif plane == 'edge': 
            a['center'] = (z+s//2,x+s//2,y)
            v1c=v1[z:z+s,x:x+s,y]
            v2c=v2[z:z+s,x:x+s,y]
    else: # use whole image
        if plane == 'beam': # This is plane with the full size images
            a['center'] = (z,x+s//2,y+s//2)
            v1c=v1[z,:,:]
            v2c=v2[z,:,:]
        elif plane == 'tilt': # This is the plane that will have missing wedge artifacts
            a['center'] = (z+s//2,x,y+s//2)
            v1c=v1[:,x,:]
            v2c=v2[:,x,:]
        elif plane == 'edge': 
            a['center'] = (z+s//2,x+s//2,y)
            v1c=v1[:,:,y]
            v2c=v2[:,:,y]

    a['mean_pix'] = (v1c.mean() + v2c.mean())/2.
    if( v1c.max() > 0 and v2c.max() > 0):
        FSCvol = FSCPlot(v1c,v2c,a['snrt'],a['rt'],a['rad_apod'],a['ax_apod'])
        a['resolution'] = a['pixel_size']/FSCvol.get_intersect()
        if a['savefig']:
            FSCvol.plot()
            FSCvol.save_fig(a['prefix'])
    else:
        a['resolution'] = -1
    my_corr_coef = np.corrcoef(np.ndarray.flatten(v1c),np.ndarray.flatten(v2c))   
    a['corr'] = my_corr_coef[0,1]
    mrc.close()
    return a.copy()
    

def makedir(dn):
    if os.path.exists(dn):
        if os.path.isdir(dn):
            return
        else:
            raise Exception("Dir %s exists as a file"%dn)
    else:
        os.makedirs(dn)

def resolution_measure_2D(vol1, vol2, num_cores, cutout_size=-1, \
    project_name='FSC', sub_region=-1, sub_sampling_zxy = [1,1,1], \
    use_json=False, slice_step = 1, \
    snrt = 0.2071, pixel_size = 1, z_clip = None,\
    ofn = None, plane = 'beam'):
    makedir(project_name)
    
    z_st,x_st,y_st = (0,0,0)
    with mrcfile.open(vol1) as mrc:
        z_size,x_size,y_size = mrc.data.shape
    if z_clip is not None:
        z_st = z_clip[0]
        z_size = z_clip[1] - z_clip[0] + 1 

    # Note: sub_region z overrides z_clip
    if len(sub_region) == 3:
        if sub_region[0] > 0:
            z_st = (z_size - sub_region[0])//2
            z_size = sub_region[0]
        if sub_region[1] > 0:        
            x_st = (x_size - sub_region[1])//2
            x_size = sub_region[1] 
        if sub_region[2] > 0:       
            y_st = (y_size - sub_region[2])//2
            y_size = sub_region[2]
	
    elif len(sub_region == 1):
        if sub_region > 0:
             z_st = (z_size - sub_region)//2
             x_st = (x_size - sub_region)//2
             y_st = (y_size - sub_region)//2
             z_size = sub_region
             x_size = sub_region
             y_size = sub_region        
    mrc.close()
    tmp = dict()
    
    if os.path.exists("%s/default.json"%project_name) and use_json:
        tmp = json.load(open("%s/default.json"%project_name))
    else:
        tmp['fn1']=vol1
        tmp['fn2']=vol2
        tmp['plane'] = plane
        tmp['cutout_size'] = cutout_size
        tmp['snrt'] = snrt
        tmp['rt'] = 6
        tmp['rad_apod'] = 60   
        tmp['ax_apod'] = 60
        tmp['pixel_size'] = pixel_size # nm
        tmp['savefig']=False
        tmp['prefix'] =""
        json.dump(tmp,open("%s/default.json"%project_name,'w'))
    
    print("Estimating the resolution by FSC...")
    startfsc = time()
    
    # prepare the pool
    pool = multiprocessing.Pool(num_cores)
    par_args = []
    
    print("Base arguments: %s"%tmp) 
    if cutout_size > 0:
        if plane == 'beam':
            for k in range(z_st,z_st + z_size, sub_sampling_zxy[0]):
                for i in range( x_st, x_st + x_size - cutout_size, cutout_size*sub_sampling_zxy[1]):
                    for j in range( y_st, y_st + y_size - cutout_size, cutout_size*sub_sampling_zxy[2]):
                        tmp['top_left'] = (k,i,j) 
                        par_args.append(tmp.copy())
        elif plane == 'tilt': # Warning - tilt and edge have no been updated to work with subsampling
            for i in range( x_st, x_st + x_size, slice_step):
                for j in range( y_st, y_st + y_size - cutout_size, cutout_size):
                    for k in range(z_st,z_st + z_size - cutout_size, cutout_size):
                        tmp['top_left'] = (k,i,j) 
                        par_args.append(tmp.copy())
        elif plane == 'edge':
            for j in range( y_st, y_st + y_size, slice_step):
                for i in range( x_st, x_st + x_size - cutout_size, cutout_size):
                    for k in range(z_st,z_st + z_size-cutout_size,cutout_size):
                        tmp['top_left'] = (k,i,j) 
                        par_args.append(tmp.copy())

    else: # Use whole image
        if plane == 'beam':
            for k in range(z_st,z_st + z_size,slice_step):
                i = x_st
                j = y_st
                tmp['top_left'] = (k,i,j) 
                par_args.append(tmp.copy())
        elif plane == 'tilt': 
            for i in range( x_st, x_st + x_size, slice_step):
                j = y_st
                k = z_st 
                tmp['top_left'] = (k,i,j) 
                par_args.append(tmp.copy())
        elif plane == 'edge':
            for j in range( y_st, y_st + y_size, slice_step):
                i = x_st
                k = z_st
                tmp['top_left'] = (k,i,j) 
                par_args.append(tmp.copy())
    #run 
    print("Running across %s cores"%num_cores)
    print(len(par_args))
    ret = list(tqdm.tqdm(pool.imap(parallel_FSC_worker, par_args), total=len(par_args)))

    if ofn is None:
        ofn = "%s/FSC_%s.csv"%(project_name,cutout_size)
    print("Outputting to %s"%ofn)
    of = open(ofn,'w')
    for r in ret:
        tl = "%s, %s, %s"%r['center']
        of.write("%s, %s, %s\n"%(tl,r['resolution'],r['mean_pix']))
    of.close()
    pool.close()
