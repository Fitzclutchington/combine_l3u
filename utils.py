import netCDF4
import numpy as np

def read_var(filepath, variable):
    cdf = netCDF4.Dataset(filepath)
    data = np.squeeze(cdf[variable][:])
    cdf.close()
    return data

def read_pixels(point_file):
    points = []
    with open(point_file, 'r') as f:
        for line in f:
            point = map(int,line.strip().strip('()').split(','))
            points.append((point[0],point[1]))
    return points

def get_dimensions(granule_file, height_id, width_id):
    cdf = netCDF4.Dataset(granule_file)
    height = cdf.dimensions[height_id].size
    width = cdf.dimensions[width_id].size
    return height, width

def get_crops(crop_file):
    crops = []
    with open(crop_file,'r') as f:
        for line in f:
            coords = line.split(',')
            ys = map(int,coords[0].split(':'))
            xs = map(int,coords[1].split(':'))
            crops.append(np.s_[ys[0]:ys[1], xs[0]:xs[1]])
    return crops

def get_files(files):
    filenames = []
    with open(files,'r') as f:
        for line in f:
            filenames.append(line.strip())
    return filenames

def get_hour(filepath):
    filename = filepath.split('/')[-1]
    hour = filename[8:10]
    return hour

def compute_ending(filepath):
    filename = filepath.split('/')[-1]
    ending = filename.split('-')
    return '-'+'-'.join(ending[1:])

def save_to_netCDF(sst, sza, time, day, save_loc):
    variable_names = ["sea_surface_temperature", "solar_zenith_angle",
                      "retrieval_time", "day_mask"]
    variable_data = [sst, sza, time, day]
    l2p_name = "l2p_flags"

    h, w = sst.shape
    
    rootgrp = netCDF4.Dataset(save_loc, "w", format="NETCDF4")
    height = rootgrp.createDimension("height", h)
    width = rootgrp.createDimension("width", w)

    for (variable_name, data) in zip(variable_names,variable_data):
        var = rootgrp.createVariable(variable_name,"f4",("height","width"),zlib=True)
        var[:] = data

    rootgrp.close()
    
