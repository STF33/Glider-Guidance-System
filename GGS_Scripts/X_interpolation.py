# =========================
# X - IMPORTS
# =========================

import numpy as np
import os
from scipy.interpolate import interp1d
import xarray as xr

# =========================
# INTERPOLATION CALCULATION
# =========================

### FUNCTION:
def interpolation_function(u, v, temp, sal, depths, max_depth):
    
    '''
    Interpolate and process calculations for bin averages and depth average data for a single XY grid point.

    Args:
    - u, v, temp, sal (xr.DataArray): DataArrays of ocean model variables over depth.
    - depths (xr.DataArray): DataArray of depth values.
    - max_depth (int): Maximum depth for bin calculation.

    Returns:
    - Tuple of arrays containing bin averages and depth averages for u, v, temperature, and salinity.
    '''

    # Initialize bin output arrays
    num_bins = max_depth
    u_bin_avg = np.full(num_bins, np.nan)
    v_bin_avg = np.full(num_bins, np.nan)
    temp_bin_avg = np.full(num_bins, np.nan)
    salt_bin_avg = np.full(num_bins, np.nan)
    mag_bin_avg = np.full(num_bins, np.nan)
    dir_bin_avg = np.full(num_bins, np.nan)

    # Check for valid data
    valid_mask = ~np.isnan(u) & ~np.isnan(v) & ~np.isnan(temp) & ~np.isnan(sal)
    if np.any(valid_mask):
        valid_depths = depths[valid_mask]

        # Create interpolation functions using valid data
        u_interp = interp1d(valid_depths, u[valid_mask], bounds_error=False, fill_value="extrapolate")
        v_interp = interp1d(valid_depths, v[valid_mask], bounds_error=False, fill_value="extrapolate")
        temp_interp = interp1d(valid_depths, temp[valid_mask], bounds_error=False, fill_value="extrapolate")
        salt_interp = interp1d(valid_depths, sal[valid_mask], bounds_error=False, fill_value="extrapolate")

        # Interpolating and averaging over each bin
        for bin_idx in range(num_bins):
            if bin_idx < valid_depths.max():
                u_bin_avg[bin_idx] = u_interp(bin_idx)
                v_bin_avg[bin_idx] = v_interp(bin_idx)
                temp_bin_avg[bin_idx] = temp_interp(bin_idx)
                salt_bin_avg[bin_idx] = salt_interp(bin_idx)

                # Calculate magnitude and direction
                magnitude = np.sqrt(u_bin_avg[bin_idx]**2 + v_bin_avg[bin_idx]**2)
                direction = (270 - np.rad2deg(np.arctan2(v_bin_avg[bin_idx], u_bin_avg[bin_idx])) + 180) % 360
                mag_bin_avg[bin_idx] = magnitude
                dir_bin_avg[bin_idx] = direction

        # Compute depth-averaged values
        u_depth_avg = np.nansum(u_bin_avg[:num_bins]) / valid_depths.max()
        v_depth_avg = np.nansum(v_bin_avg[:num_bins]) / valid_depths.max()
        temp_depth_avg = np.nansum(temp_bin_avg[:num_bins]) / valid_depths.max()
        salt_depth_avg = np.nansum(salt_bin_avg[:num_bins]) / valid_depths.max()
        mag_depth_avg = np.nansum(mag_bin_avg[:num_bins]) / valid_depths.max()
        dir_depth_avg = np.nansum(dir_bin_avg[:num_bins]) / valid_depths.max()
    else:
        # If no valid data, use the full range of bins and set depth averages to NaN
        u_depth_avg = np.nan
        v_depth_avg = np.nan
        temp_depth_avg = np.nan
        salt_depth_avg = np.nan
        mag_depth_avg = np.nan
        dir_depth_avg = np.nan

    return (u_bin_avg, v_bin_avg, temp_bin_avg, salt_bin_avg, mag_bin_avg, dir_bin_avg, 
            u_depth_avg, v_depth_avg, temp_depth_avg, salt_depth_avg, mag_depth_avg, dir_depth_avg)

### FUNCTION:
def interpolation_model(config, directory, model_data):
    
    '''
    Compute depth-averaged values using xarray's apply_ufunc for the entire dataset.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - model_data (xarray.Dataset): Ocean model dataset.

    Returns:
    - Depth average data and bin average data as xarray.Dataset.
    '''
    print("\n")
    print("### INTERPOLATING MODEL DATA ###")
    print("\n")

    # Chunk the data for optimal parallel processing
    chunked_data = model_data.chunk({'y': 50, 'x': 50})
    max_depth = config['max_depth']
    depth_data = chunked_data['depth']

    # Apply interpolation_function to each grid point
    results = xr.apply_ufunc(
        interpolation_function,
        chunked_data['u'],
        chunked_data['v'],
        chunked_data['temperature'],
        chunked_data['salinity'],
        depth_data,
        kwargs={'max_depth': max_depth},
        input_core_dims=[['depth'], ['depth'], ['depth'], ['depth'], ['depth']],
        output_core_dims=[['bin'], ['bin'], ['bin'], ['bin'], ['bin'], ['bin'], [], [], [], [], [], []],
        dask_gufunc_kwargs={'output_sizes': {'bin': max_depth}},
        vectorize=True,
        dask='parallelized')
    
    # Unpack results
    u_bin_avg, v_bin_avg, temp_bin_avg, salt_bin_avg, mag_bin_avg, dir_bin_avg, u_depth_avg, v_depth_avg, temp_depth_avg, salt_depth_avg, mag_depth_avg, dir_depth_avg = [da.data for da in results]

    # Create datasets
    depth_average_data = xr.Dataset({
        'u_depth_avg': (('y', 'x'), u_depth_avg),
        'v_depth_avg': (('y', 'x'), v_depth_avg),
        'temp_depth_avg': (('y', 'x'), temp_depth_avg),
        'salt_depth_avg': (('y', 'x'), salt_depth_avg),
        'mag_depth_avg': (('y', 'x'), mag_depth_avg),
        'dir_depth_avg': (('y', 'x'), dir_depth_avg)
    }, coords={'lat': model_data.lat, 'lon': model_data.lon})
    bin_average_data = xr.Dataset({
        'u_bin_avg': (('y', 'x', 'bin'), u_bin_avg),
        'v_bin_avg': (('y', 'x', 'bin'), v_bin_avg),
        'temp_bin_avg': (('y', 'x', 'bin'), temp_bin_avg),
        'salt_bin_avg': (('y', 'x', 'bin'), salt_bin_avg),
        'mag_bin_avg': (('y', 'x', 'bin'), mag_bin_avg),
        'dir_bin_avg': (('y', 'x', 'bin'), dir_bin_avg)
    }, coords={'lat': model_data.lat, 'lon': model_data.lon, 'bin': np.arange(max_depth)})

    # Save datasets as NetCDF files
    depth_average_data.to_netcdf(os.path.join(directory, f"{config['glider_name']}_DepthAverageData.nc"))
    bin_average_data.to_netcdf(os.path.join(directory, f"{config['glider_name']}_BinAverageData.nc"))

    print("INTERPOLATED DATASETS CREATED")

    return depth_average_data, bin_average_data
