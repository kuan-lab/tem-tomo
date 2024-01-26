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
    project_name='FSC', sub_region=-1, use_json=False, slice_step = 1, \
    snrt = 0.2071, pixel_size = 1, \
    ofn = None, plane = 'beam'):
    makedir(project_name)
    
    z_st,x_st,y_st = (0,0,0)
    with mrcfile.open(vol1) as mrc:
        z_size,x_size,y_size = mrc.data.shape

    if sub_region > 0:
        z_st = (z_size - sub_region)//2
        x_st = (x_size - sub_region)//2
        y_st = (y_size - sub_region)//2
        z_size = sub_region
        x_size = sub_region
        y_size = sub_region
    
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
        if plane == 'tilt': # This is plane with the full size images
            for i in range( x_st, x_st + x_size + 1, slice_step):
                for j in range( y_st, y_st + y_size - cutout_size+1, cutout_size):
                    for k in range(z_st,z_st + z_size - cutout_size+1, cutout_size):
                        tmp['top_left'] = (k,i,j) 
                        par_args.append(tmp.copy())
        elif plane == 'edge':
            for j in range( y_st, y_st + y_size +1, slice_step):
                for i in range( x_st, x_st + x_size - cutout_size+1, cutout_size):
                    for k in range(z_st,z_st + z_size-cutout_size+1,cutout_size):
                        tmp['top_left'] = (k,i,j) 
                        par_args.append(tmp.copy())
        elif plane == 'beam':
            for k in range(z_st,z_st + z_size+1,slice_step):
                for i in range( x_st, x_st + x_size - cutout_size+1, cutout_size):
                    for j in range( y_st, y_st + y_size - cutout_size+1, cutout_size):
                        tmp['top_left'] = (k,i,j) 
                        par_args.append(tmp.copy())
    else: # Use whole image
        if plane == 'tilt': # This is plane with the full size images
            for i in range( x_st, x_st + x_size + 1, slice_step):
                j = y_st
                k = z_st 
                tmp['top_left'] = (k,i,j) 
                par_args.append(tmp.copy())
        elif plane == 'edge':
            for j in range( y_st, y_st + y_size +1, slice_step):
                i = x_st
                k = z_st
                tmp['top_left'] = (k,i,j) 
                par_args.append(tmp.copy())
        elif plane == 'beam':
            for k in range(z_st,z_st + z_size+1,slice_step):
                i = x_st
                j = y_st
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

######### Haven't edited command line version yet
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("vol1",help="fn of first mrc volume")
    parser.add_argument("vol2",help="fn of second mrc volume")
    parser.add_argument("num_cores",default=1,type=int,help="number of cores to par over")
    parser.add_argument("cutout_size",default=200,type=int,help="len of a side of the cube to partiton dataset")
    parser.add_argument("-ps","--pixel_size",default=50,type=int,help="resolution in nm of dataset")
    parser.add_argument("-sn","--snrt",default=0.2071,type=float,help="snrt value to use, .2071 is default however, 1/7 is 0.143") 
    parser.add_argument("--param_sweep",action="store_true",help="flag to do an initial parameter sweep")
    parser.add_argument("--use_json",action="store_true",help="flag to set using the default.json")
    parser.add_argument("-sub","--sub_region",default=-1, type=int, help="size of the subvolume to run on, will take a centered cube of side length specified")
    args = parser.parse_args()

    project_name = ""
    for ind in range(len(args.vol1)):
        if args.vol1[ind] == args.vol2[ind]:
            project_name += args.vol1[ind]
        else:
            if project_name == "":
                project_name = "resolution_test"
            break
    
    project_name = "%s_resolution"%project_name
    print("Project name: %s"%project_name)
    makedir(project_name)
    
    
    
    #confirm these are h5py files before launching par jobs
    #f1 = zarr.open(args.vol1,'r')
    #f2 = zarr.open(args.vol2,'r')
    
    z_st,x_st,y_st = (0,0,0)
    with mrcfile.open(args.vol1) as mrc:
        z_size,x_size,y_size = mrc.data.shape
    #z_size,x_size,y_size = f1['data'].shape
    if args.sub_region > 0:
        z_st = (z_size - args.sub_region)//2
        x_st = (x_size - args.sub_region)//2
        y_st = (y_size - args.sub_region)//2
        z_size = args.sub_region
        x_size = args.sub_region
        y_size = args.sub_region
    
    
    # TODO assert here if f2 dosnt have same shape?
    # prepare the dictionary to be passed to worker threads
    
    
    tmp = dict()
    
    if os.path.exists("%s/default.json"%project_name) and args.use_json:
        tmp = json.load(open("%s/default.json"%project_name))
    else:
        tmp['fn1']=args.vol1
        tmp['fn2']=args.vol2
        tmp['cutout_size'] = args.cutout_size
        tmp['snrt'] = args.snrt
        tmp['rt'] = 6
        tmp['rad_apod'] = 60   
        tmp['ax_apod'] = 60
        tmp['pixel_size'] = args.pixel_size # nm
        tmp['savefig']=False
        tmp['prefix'] =""
        json.dump(tmp,open("%s/default.json"%project_name,'w'))
    
    print("Estimating the resolution by FSC...")
    startfsc = time()
    
    
    # prepare the pool
    pool = multiprocessing.Pool(args.num_cores)
    par_args = []
    
    if not args.param_sweep:
        print("Base arguments: %s"%tmp)
        for i in range( x_st, x_st + x_size - args.cutout_size+1, args.cutout_size):
            for j in range( y_st, y_st + y_size - args.cutout_size+1, args.cutout_size):
                for k in range(z_st,z_st + z_size-args.cutout_size+1,args.cutout_size):
                    tmp['top_left'] = (k,i,j) 
                    par_args.append(tmp.copy())
    
        #run 
        print("Running across %s cores"%args.num_cores)
    
        ret = list(tqdm.tqdm(pool.imap(parallel_FSC_worker, par_args), total=len(par_args)))
    
        ofn = "%s/FSC_%s.csv"%(project_name,args.cutout_size)
        print("Outputting to %s"%ofn)
        of = open(ofn,'w')
        for r in ret:
            tl = "%s, %s, %s"%r['center']
            of.write("%s, %s %s\n"%(tl,r['resolution'],r['mean_pix']))
    else:
        print("Running parameter sweep")
        par_args = []
        par_args.extend(sweep_param('cutout_size',range(50,800,100),tmp,project_name))
        par_args.extend(sweep_param(['ax_apod','rad_apod'],range(10,1000,100),tmp,project_name))
        par_args.extend(sweep_param('rt',range(2,14,2),tmp,project_name))
        print("Running across %s cores"%args.num_cores)
        pool.map(parallel_FSC_worker, par_args)


