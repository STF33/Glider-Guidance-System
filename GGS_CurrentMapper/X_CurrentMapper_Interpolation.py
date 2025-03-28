"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import numpy as np
import os
from scipy.interpolate import interp1d
import xarray as xr

from X_CurrentMapper_Functions import format_save_datetime, print_starttime, print_endtime, print_runtime

# =========================

### FUNCTION:
def interpolation_model(u, v, depths, max_bins, config_bins):

    '''
    Interpolate model data to depth-averaged and bin-averaged values.

    Args:
    - u (np.array): Zonal velocity.
    - v (np.array): Meridional velocity.
    - depths (np.array): Depth values.
    - max_bins (int): Maximum number of bins.
    - config_bins (int): Number of bins based on the configuration.

    Returns:
    - u_bin_avg (np.array): Zonal velocity bin-averaged.
    - v_bin_avg (np.array): Meridional velocity bin-averaged.
    - mag_bin_avg (np.array): Magnitude bin-averaged.
    - dir_bin_avg (np.array): Direction bin-averaged.
    - u_depth_avg (float): Zonal velocity depth-averaged.
    - v_depth_avg (float): Meridional velocity depth-averaged.
    - mag_depth_avg (float): Magnitude depth-averaged.
    - dir_depth_avg (float): Direction depth-averaged.
    '''

    bins = int(np.ceil(max_bins))
    u_bin_avg = np.full(bins, np.nan)
    v_bin_avg = np.full(bins, np.nan)
    mag_bin_avg = np.full(bins, np.nan)
    dir_bin_avg = np.full(bins, np.nan)
    
    valid_mask = ~np.isnan(u) & ~np.isnan(v)
    if np.any(valid_mask):
        valid_depths = depths[valid_mask]
        u_valid, v_valid = u[valid_mask], v[valid_mask]

        valid_bins = int(np.ceil(valid_depths.max())) + 1
        target_bins = np.arange(min(valid_bins, config_bins))

        u_interp = interp1d(valid_depths, u_valid, bounds_error=False, fill_value='extrapolate')
        v_interp = interp1d(valid_depths, v_valid, bounds_error=False, fill_value='extrapolate')

        u_bin_avg[target_bins] = u_interp(target_bins)
        v_bin_avg[target_bins] = v_interp(target_bins)
        mag_bin_avg[target_bins] = np.sqrt(u_bin_avg[target_bins]**2 + v_bin_avg[target_bins]**2)
        dir_bin_avg[target_bins] = (np.degrees(np.arctan2(v_bin_avg[target_bins], u_bin_avg[target_bins])) + 360) % 360

        valid = ~np.isnan(u_bin_avg[target_bins])
        counts_valid = valid.sum()
        u_depth_avg = np.nansum(u_bin_avg[target_bins][valid]) / counts_valid
        v_depth_avg = np.nansum(v_bin_avg[target_bins][valid]) / counts_valid
        mag_depth_avg = np.nansum(mag_bin_avg[target_bins][valid]) / counts_valid
        dir_depth_avg = np.nansum(dir_bin_avg[target_bins][valid]) / counts_valid
    else:
        u_depth_avg = np.nan
        v_depth_avg = np.nan
        mag_depth_avg = np.nan
        dir_depth_avg = np.nan
    return (u_bin_avg, v_bin_avg, mag_bin_avg, dir_bin_avg, u_depth_avg, v_depth_avg, mag_depth_avg, dir_depth_avg)

### FUNCTION:
def interpolate_rtofs(config, directory, model_data, chunk=False, save_depth_average=True, save_bin_average=False):
    
    '''
    Compute depth-averaged values using xarray's apply_ufunc for the entire dataset.
    Optimized for RTOFS model datasets.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - model_data (xarray.Dataset): Ocean model dataset.
    - chunk (bool): Chunk the data.
        - default: 'False'
    - save_depth_average (bool): Save the depth average data.
        - default: 'True'
    - save_bin_average (bool): Save the bin average data.
        - default: 'False'

    Returns:
    - model_depth_average (xarray.Dataset): Depth average data.
    - model_bin_average (xarray.Dataset): Bin average data.
    '''

    print("\n### INTERPOLATING RTOFS MODEL DATA ###\n")
    start_time = print_starttime()

    if chunk:
        model_data = model_data.chunk({'Y': 'auto', 'X': 'auto'})

    model_data = model_data.load()

    config_depth = config['MISSION']['max_depth']
    config_bins = config_depth + 1
    max_depth = model_data.Depth.max().item()
    max_bins = max_depth + 1

    results = xr.apply_ufunc(
        interpolation_model,
        model_data['u'],
        model_data['v'],
        model_data['Depth'],
        kwargs={'max_bins': max_bins, 'config_bins': config_bins},
        input_core_dims=[['Depth'], ['Depth'], ['Depth']],
        output_core_dims=[['bin'], ['bin'], ['bin'], ['bin'], [], [], [], []],
        output_dtypes=[float, float, float, float, float, float, float, float],
        dask_gufunc_kwargs={'output_sizes': {'bin': max_depth}},
        vectorize=True
    )

    u_bin_avg, v_bin_avg, mag_bin_avg, dir_bin_avg, u_depth_avg, v_depth_avg, mag_depth_avg, dir_depth_avg = results

    model_depth_average = xr.Dataset({
        'u_depth_avg': (('Y', 'X'), u_depth_avg.data),
        'v_depth_avg': (('Y', 'X'), v_depth_avg.data),
        'mag_depth_avg': (('Y', 'X'), mag_depth_avg.data),
        'dir_depth_avg': (('Y', 'X'), dir_depth_avg.data)
    }, coords={'lat': model_data.Latitude, 'lon': model_data.Longitude})
    model_depth_average = model_depth_average.expand_dims('time')
    model_depth_average.attrs['model_datetime'] = model_data.attrs['model_datetime']
    model_depth_average.attrs['model_name'] = model_data.attrs['model_name']

    model_bin_average = xr.Dataset({
        'u_bin_avg': (('Y', 'X', 'bin'), u_bin_avg.data),
        'v_bin_avg': (('Y', 'X', 'bin'), v_bin_avg.data),
        'mag_bin_avg': (('Y', 'X', 'bin'), mag_bin_avg.data),
        'dir_bin_avg': (('Y', 'X', 'bin'), dir_bin_avg.data)
    }, coords={'lat': model_data.Latitude, 'lon': model_data.Longitude, 'bin': np.arange(max_bins)})
    model_bin_average = model_bin_average.expand_dims('time')
    model_bin_average.attrs['model_datetime'] = model_data.attrs['model_datetime']
    model_bin_average.attrs['model_name'] = model_data.attrs['model_name']

    if save_depth_average:
        model_datetime = model_data.attrs['model_datetime']
        file_datetime = format_save_datetime(model_datetime)
        model_depth_average.to_netcdf(os.path.join(directory, f"{config['MISSION'].get('mission_name', 'UnknownMission')}_RTOFS_DepthAverage_{file_datetime}.nc"), unlimited_dims=['time'])
    if save_bin_average:
        model_datetime = model_data.attrs['model_datetime']
        file_datetime = format_save_datetime(model_datetime)
        model_bin_average.to_netcdf(os.path.join(directory, f"{config['MISSION'].get('mission_name', 'UnknownMission')}_RTOFS_BinAverage_{file_datetime}.nc"), unlimited_dims=['time'])
    
    end_time = print_endtime()
    print_runtime(start_time, end_time)
    
    return model_depth_average, model_bin_average

### FUNCTION:
def interpolate_cmems(config, directory, model_data, chunk=False, save_depth_average=True, save_bin_average=False):

    '''
    Compute depth-averaged values using xarray's apply_ufunc for the entire dataset.
    Optimized for CMEMS model datasets.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - model_data (xarray.Dataset): Ocean model dataset.
    - chunk (bool): Chunk the data.
        - default: 'False'
    - save_depth_average (bool): Save the depth average data.
        - default: 'True'
    - save_bin_average (bool): Save the bin average data.
        - default: 'False'

    Returns:
    - model_depth_average (xarray.Dataset): Depth average data.
    - model_bin_average (xarray.Dataset): Bin average data.
    '''

    print("\n### INTERPOLATING CMEMS MODEL DATA ###\n")
    start_time = print_starttime()

    if chunk:
        model_data = model_data.chunk({'lat': 'auto', 'lon': 'auto'})

    model_data = model_data.load()

    config_depth = config['MISSION']['max_depth']
    config_bins = config_depth + 1
    max_depth = model_data.depth.max().item()
    max_bins = max_depth + 1

    results = xr.apply_ufunc(
        interpolation_model,
        model_data['u'],
        model_data['v'],
        model_data['depth'],
        kwargs={'max_bins': max_bins, 'config_bins': config_bins},
        input_core_dims=[['depth'], ['depth'], ['depth']],
        output_core_dims=[['bin'], ['bin'], ['bin'], ['bin'], [], [], [], []],
        output_dtypes=[float, float, float, float, float, float, float, float],
        dask_gufunc_kwargs={'output_sizes': {'bin': max_depth}},
        vectorize=True
    )

    u_bin_avg, v_bin_avg, mag_bin_avg, dir_bin_avg, u_depth_avg, v_depth_avg, mag_depth_avg, dir_depth_avg = results

    u_depth_avg = u_depth_avg.data.squeeze()
    v_depth_avg = v_depth_avg.data.squeeze()
    mag_depth_avg = mag_depth_avg.data.squeeze()
    dir_depth_avg = dir_depth_avg.data.squeeze()
    u_bin_avg = u_bin_avg.data.squeeze()
    v_bin_avg = v_bin_avg.data.squeeze()
    mag_bin_avg = mag_bin_avg.data.squeeze()
    dir_bin_avg = dir_bin_avg.data.squeeze()

    model_depth_average = xr.Dataset({
        'u_depth_avg': (('lat', 'lon'), u_depth_avg.data),
        'v_depth_avg': (('lat', 'lon'), v_depth_avg.data),
        'mag_depth_avg': (('lat', 'lon'), mag_depth_avg.data),
        'dir_depth_avg': (('lat', 'lon'), dir_depth_avg.data)
    }, coords={'lat': model_data.lat, 'lon': model_data.lon})
    model_depth_average = model_depth_average.expand_dims('time')
    model_depth_average.attrs['model_datetime'] = model_data.attrs['model_datetime']
    model_depth_average.attrs['model_name'] = model_data.attrs['model_name']

    model_bin_average = xr.Dataset({
        'u_bin_avg': (('lat', 'lon', 'bin'), u_bin_avg.data),
        'v_bin_avg': (('lat', 'lon', 'bin'), v_bin_avg.data),
        'mag_bin_avg': (('lat', 'lon', 'bin'), mag_bin_avg.data),
        'dir_bin_avg': (('lat', 'lon', 'bin'), dir_bin_avg.data)
    }, coords={'lat': model_data.lat, 'lon': model_data.lon, 'bin': np.arange(max_bins)})
    model_bin_average = model_bin_average.expand_dims('time')
    model_bin_average.attrs['model_datetime'] = model_data.attrs['model_datetime']
    model_bin_average.attrs['model_name'] = model_data.attrs['model_name']
    
    if save_depth_average:
        model_datetime = model_data.attrs['model_datetime']
        file_datetime = format_save_datetime(model_datetime)
        model_depth_average.to_netcdf(os.path.join(directory, f"{config['MISSION'].get('mission_name', 'UnknownMission')}_CMEMS_DepthAverage_{file_datetime}.nc"), unlimited_dims=['time'])
    if save_bin_average:
        model_datetime = model_data.attrs['model_datetime']
        file_datetime = format_save_datetime(model_datetime)
        model_bin_average.to_netcdf(os.path.join(directory, f"{config['MISSION'].get('mission_name', 'UnknownMission')}_CMEMS_BinAverage_{file_datetime}.nc"), unlimited_dims=['time'])
    
    end_time = print_endtime()
    print_runtime(start_time, end_time)
    
    return model_depth_average, model_bin_average

### FUNCTION:
def interpolate_gofs(config, directory, model_data, chunk=False, save_depth_average=True, save_bin_average=False):
    
    '''
    Compute depth-averaged values using xarray's apply_ufunc for the entire dataset.
    Optimized for GOFS model datasets.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - model_data (xarray.Dataset): Ocean model dataset.
    - chunk (bool): Whether to chunk the data for dask processing.
        - default: 'False'
    - save_depth_average (bool): Whether to save the depth-averaged data.
        - default: 'True'
    - save_bin_average (bool): Whether to save the bin-averaged data.
        - default: 'False'

    Returns:
    - model_depth_average (xarray.Dataset): The dataset containing depth-averaged data.
    - model_bin_average (xarray.Dataset): The dataset containing bin-averaged data.
    '''

    print("\n### INTERPOLATING GOFS MODEL DATA ###\n")
    start_time = print_starttime()

    if chunk:
        model_data = model_data.chunk({'lat': 'auto', 'lon': 'auto'})

    model_data.load()

    config_depth = config['MISSION']['max_depth']
    config_bins = config_depth + 1
    max_depth = model_data.depth.max().item()
    max_bins = max_depth + 1

    results = xr.apply_ufunc(
        interpolation_model,
        model_data['u'],
        model_data['v'],
        model_data['depth'],
        kwargs={'max_bins': max_bins, 'config_bins': config_bins},
        input_core_dims=[['depth'], ['depth'], ['depth']],
        output_core_dims=[['bin'], ['bin'], ['bin'], ['bin'], [], [], [], []],
        output_dtypes=[float, float, float, float, float, float, float, float],
        dask_gufunc_kwargs={'output_sizes': {'bin': max_bins}},
        vectorize=True
    )

    u_bin_avg, v_bin_avg, mag_bin_avg, dir_bin_avg, u_depth_avg, v_depth_avg, mag_depth_avg, dir_depth_avg = results

    model_depth_average = xr.Dataset({
        'u_depth_avg': (('lat', 'lon'), u_depth_avg.data),
        'v_depth_avg': (('lat', 'lon'), v_depth_avg.data),
        'mag_depth_avg': (('lat', 'lon'), mag_depth_avg.data),
        'dir_depth_avg': (('lat', 'lon'), dir_depth_avg.data)
    }, coords={'lat': model_data.lat, 'lon': model_data.lon})
    model_depth_average = model_depth_average.expand_dims('time')
    model_depth_average.attrs['model_datetime'] = model_data.attrs['model_datetime']
    model_depth_average.attrs['model_name'] = model_data.attrs['model_name']
    
    model_bin_average = xr.Dataset({
        'u_bin_avg': (('lat', 'lon', 'bin'), u_bin_avg.data),
        'v_bin_avg': (('lat', 'lon', 'bin'), v_bin_avg.data),
        'mag_bin_avg': (('lat', 'lon', 'bin'), mag_bin_avg.data),
        'dir_bin_avg': (('lat', 'lon', 'bin'), dir_bin_avg.data)
    }, coords={'lat': model_data.lat, 'lon': model_data.lon, 'bin': np.arange(max_bins)})
    model_bin_average = model_bin_average.expand_dims('time')
    model_bin_average.attrs['model_datetime'] = model_data.attrs['model_datetime']
    model_bin_average.attrs['model_name'] = model_data.attrs['model_name']
    
    if save_depth_average:
        model_datetime = model_data.attrs['model_datetime']
        file_datetime = format_save_datetime(model_datetime)
        model_depth_average.to_netcdf(os.path.join(directory, f"{config['MISSION'].get('mission_name', 'UnknownMission')}_GOFS_DepthAverage_{file_datetime}.nc"), unlimited_dims=['time'])
    if save_bin_average:
        model_datetime = model_data.attrs['model_datetime']
        file_datetime = format_save_datetime(model_datetime)
        model_bin_average.to_netcdf(os.path.join(directory, f"{config['MISSION'].get('mission_name', 'UnknownMission')}_GOFS_BinAverage_{file_datetime}.nc"), unlimited_dims=['time'])
    
    end_time = print_endtime()
    print_runtime(start_time, end_time)

    return model_depth_average, model_bin_average
