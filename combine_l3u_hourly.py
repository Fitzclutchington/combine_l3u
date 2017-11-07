import glob
import sys
import os
import datetime

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.io as sio
from scipy import interpolate
import netCDF4

import utils
import window_functions

def interpolate_quantized_matrix(mat, interp_mat):
    diff = np.diff(mat)
    end = mat.shape[1]
    num_rows = mat.shape[0]
    x_inds = np.arange(end)
    for i in range(num_rows):
        edges = np.where(diff[i,:] != 0)[0]
        mid_ind = np.zeros(edges.shape[0] + 1,dtype=np.uint16)
        mid_ind[-1] = end-1
        mid_ind[1:-1] = np.round((edges[0:-1]+edges[1:])/2)
        f = interpolate.interp1d(mid_ind, mat[i,mid_ind])
        interp_mat[i,:] = f(x_inds)

def matlab_style_gauss2D(w, sigma):
    """
    2D gaussian mask - should give the same result as MATLAB's
    fspecial('gaussian',[shape],[sigma])
    """
    shape = (2*w+1,2*w+1)
    m,n = [(ss-1.)/2. for ss in shape]
    y,x = np.ogrid[-m:m+1,-n:n+1]
    h = np.exp( -(x*x + y*y) / (2.*sigma*sigma) )
    h[ h < np.finfo(h.dtype).eps*h.max() ] = 0
    sumh = h.sum()
    if sumh != 0:
        h /= sumh
    return h

if len(sys.argv) != 5:
    print "usage:python combine_l3u_hourly.py input_l3u_folder input_l2p_folder output_concatenated_l3u_folder modis_flag"
    sys.exit()
    
save_folder = sys.argv[3]
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

modis_flag = bool(int(sys.argv[4]))

folders = { 
            'l3u' : sys.argv[1],
            'l2p' : sys.argv[2]
        }

files = {
         'l3u' : sorted(glob.glob(folders['l3u']+'*.nc')),
         'l2p' : sorted(glob.glob(folders['l2p']+"*.nc")) 
        }

sets = { 
         'l3u' : set(files['l3u']),
         'l2p' : set(files['l2p'])
        }

endings = { 
            'l3u': utils.compute_ending(files['l3u'][0]),
            'l2p': utils.compute_ending(files['l2p'][0]) 
        }

dimensions = { 
               'l3u' : (utils.get_dimensions(files['l3u'][0],'lat','lon')),
            }


dL = 0.02
p = [-dL,   90+dL/2] 
q = [ dL,  -180-dL/2]
N = 180/dL
M = 360/dL

w=5
f = matlab_style_gauss2D(w, w/2)
f_sum = f.sum()

current_day_string = files['l3u'][0].split('/')[-1].split('-')[0][0:8]

current_time = datetime.datetime.strptime(current_day_string+'0000', '%Y%m%d%H%M%S')
if modis_flag:
    time_delta = datetime.timedelta(minutes=5)
else:
    time_delta = datetime.timedelta(minutes=10)

next_day = current_time + datetime.timedelta(days=1)
next_hour = current_time + datetime.timedelta(hours=1)

l3u_sst = np.full(dimensions['l3u'],np.nan)
l3u_sza = np.full(dimensions['l3u'],np.nan)
l3u_time = np.full(dimensions['l3u'],np.nan).astype(np.float32)
l3u_day = np.full(dimensions['l3u'],-1).astype(np.int8)

H = np.zeros(l3u_sst.shape).astype(np.float32)
S = np.zeros(l3u_sst.shape).astype(np.float32)
T = np.zeros(l3u_sst.shape).astype(np.float32)

while current_time != next_day:
#while current_time != next_hour:
    date_time_string = current_time.strftime('%Y%m%d%H%M%S')
    l3u_filename = folders['l3u'] + date_time_string + endings['l3u']
    l2p_filename = folders['l2p'] + date_time_string + endings['l2p']

    # append sst
    if l3u_filename in sets['l3u']:
        current_sst = utils.read_var(l3u_filename,'sea_surface_temperature')
        QL_flags = utils.read_var(l3u_filename, 'quality_level')
        l2p_flags = utils.read_var(l3u_filename, 'l2p_flags')
        day_mask = np.bitwise_and(l2p_flags,512).astype(bool)
        QL_mask = QL_flags >= 0
        

        l3u_sst[QL_mask] = current_sst[QL_mask]
        l3u_day[day_mask] = 1
        l3u_day[~day_mask & QL_mask] = 0
        # compute time and sza
        if l2p_filename in sets['l2p']:

            sza = utils.read_var(l2p_filename,'satellite_zenith_angle')
            sza_interp = np.zeros(sza.shape).astype(np.float32)
            t_utc = np.zeros(sza.shape).astype(np.float32)
            interpolate_quantized_matrix(sza, sza_interp)
            
            d_time = utils.read_var(l2p_filename,'sst_dtime')
            current_time_val = current_time.hour + current_time.minute/60.0
            t_utc = current_time_val + d_time/(60*60)

            lat_L2 = utils.read_var(l2p_filename,'lat')
            lon_L2 = utils.read_var(l2p_filename,'lon')
            ii = np.round((lat_L2 - p[1])/p[0]).astype(np.int32) - 1
            jj = np.round((lon_L2 - q[1])/q[0]).astype(np.int32) - 1
            

            window_functions.remap(ii, jj, H,  T,  S, t_utc, sza_interp )

            N = l3u_sza.shape[0]
            M = l3u_sza.shape[1]          

            
            window_functions.smoothing(QL_flags, S, T, H, l3u_sza, l3u_time, f, f_sum)

            print "concatenated", current_time

        if (current_time + time_delta).minute == 0:
            if not modis_flag:
                save_file = (current_time - datetime.timedelta(minutes=50)).strftime('%Y%m%d%H%M%S') + ".nc"
            else:
                save_file = (current_time - datetime.timedelta(minutes=55)).strftime('%Y%m%d%H%M%S') + ".nc"
            utils.save_to_netCDF(l3u_sst, l3u_sza, l3u_time, l3u_day, save_folder+save_file)
            print "saved",(current_time - datetime.timedelta(minutes=50)).strftime('%Y%m%d%H%M%S')
            
            l3u_sst.fill(np.nan)
            l3u_time.fill(np.nan)
            l3u_sza.fill(np.nan)
            l3u_day.fill(-1)
            H.fill(0)
            T.fill(0)
            S.fill(0)
    
        current_time = current_time + time_delta

    else:

        if (current_time + time_delta).minute == 0:
            if not modis_flag:
                save_file = (current_time - datetime.timedelta(minutes=50)).strftime('%Y%m%d%H%M%S') + ".nc"
            else:
                save_file = (current_time - datetime.timedelta(minutes=55)).strftime('%Y%m%d%H%M%S') + ".nc"
            utils.save_to_netCDF(l3u_sst, l3u_sza, l3u_time, l3u_day, save_folder+save_file)
            print "saved",(current_time - datetime.timedelta(minutes=50)).strftime('%Y%m%d%H%M%S')
            
            l3u_sst.fill(np.nan)
            l3u_time.fill(np.nan)
            l3u_sza.fill(np.nan)
            l3u_day.fill(-1)
            H.fill(0)
            T.fill(0)
            S.fill(0)
            
        current_time = current_time + time_delta
        continue

