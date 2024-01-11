# =========================
# X - IMPORTS
# =========================

import numpy as np
import os
import pandas as pd
from scipy.interpolate import interp1d
import xarray as xr
from X_functions import format_datetime, get_filename_datetime

# =========================
# [RTOFS] DATA PROCESSING
# =========================

### CLASS:
class RTOFS():
    
    '''
    Class for handling RTOFS data.

    Attributes:
    - data_origin (xarray.Dataset): Original RTOFS data.
    - data (xarray.Dataset): Main RTOFS data.
    - x (np.array): X-coordinates of the RTOFS grid.
    - y (np.array): Y-coordinates of the RTOFS grid.
    - grid_lons (np.array): Longitudes of the RTOFS grid.
    - grid_lats (np.array): Latitudes of the RTOFS grid.
    - rtofs_qc (xarray.Dataset): RTOFS data for quality control.

    Methods:
    - __init__: Initialize the RTOFS instance.
    - rtofs_load: Fetch the RTOFS data from the given URL source and set its coordinates.
    - rtofs_subset: Subset the RTOFS data based on the GGS mission extent.
    - rtofs_save: Save the subset RTOFS data as a NetCDF file.
    '''

    ### FUNCTION:
    def __init__(self, datetime_index=None) -> None:
        
        '''
        Initialize the RTOFS instance.

        Args:
        - datetime_index (datetime.datetime): The datetime index for data loading.
        '''

        if datetime_index is not None:
            self.data_origin = self.rtofs_load(datetime_index)
        else:
            self.data_origin = xr.Dataset()

        self.data_origin = self.data_origin.set_coords(['lat', 'lon'])
        self.data = self.data_origin.copy()

        self.x = self.data.x.values
        self.y = self.data.y.values
        self.grid_lons = self.data.lon.values[0,:]
        self.grid_lats = self.data.lat.values[:,0]

        self.rtofs_qc = self.data_origin.copy()

    ### FUNCTION:
    def rtofs_load(self, datetime_index):
        
        '''
        Fetch the RTOFS data from the given URL source and set its coordinates.

        Args:
        - datetime_index (str): Index of the datetime to fetch.

        Returns:
        - rtofs_raw (xarray.Dataset): RTOFS data
        '''

        rtofs_access = "https://tds.marine.rutgers.edu/thredds/dodsC/cool/rtofs/rtofs_us_east_scraped"

        try:
            rtofs_raw = xr.open_dataset(rtofs_access)
            
            datetime = pd.Timestamp(datetime_index).tz_localize(None)
            time_values = rtofs_raw.time.values
            time_index = np.argmin(np.abs(time_values - np.datetime64(datetime)))
            rtofs_raw = rtofs_raw.isel(time=time_index)
            rtofs_raw.attrs['requested_datetime'] = datetime_index
            rtofs_raw.attrs['acquired_datetime'] = str(rtofs_raw.time.values)

            print("\n### RTOFS data acquired ###\n")
            print("Requested datetime:", datetime_index)
            print("Nearest datetime index in the dataset:", time_index)
            print("Acquired datetime:", rtofs_raw.attrs['acquired_datetime'])

            return rtofs_raw

        except Exception as e:
            print(f"Error fetching RTOFS data: {e}")
            return None
    
    ### FUNCTION:
    def rtofs_subset(self, config, buffer=0, subset=True):

        '''
        Subset the RTOFS data based on the GGS mission extent.

        Args:
        - config (dict): Glider Guidance System mission configuration.
        - buffer (float): Buffer (in degrees) to add to the bounding box.
            - default: 0
        - subset (bool): Subset the data.
            - default: 'True'
            - if True, the data is subset based on GGS mission extent.
            - if False, the entire RTOFS grid is used (no subsetting).

        Returns:
        - None
        '''

        if subset:
            lats, lons = zip(*config['extent'])

            min_lon, max_lon = min(lons) - buffer, max(lons) + buffer
            min_lat, max_lat = min(lats) - buffer, max(lats) + buffer

            lons_ind = np.interp([min_lon, max_lon], self.grid_lons, self.x)
            lats_ind = np.interp([min_lat, max_lat], self.grid_lats, self.y)

            extent = [
                np.floor(lons_ind[0]).astype(int),
                np.ceil(lons_ind[1]).astype(int),
                np.floor(lats_ind[0]).astype(int),
                np.ceil(lats_ind[1]).astype(int)
            ]

            self.data = self.data_origin.isel(
                x=slice(extent[0], extent[1]),
                y=slice(extent[2], extent[3])
            )

            self.data_lons = self.data.lon.values[0,:]
            self.data_lats = self.data.lat.values[:,0]
            
            self.data = self.data.where(self.data['depth'] <= config["max_depth"], drop=True)
            self.rtofs_qc = self.rtofs_qc.where(self.rtofs_qc['depth'] <= config["max_depth"], drop=True)

        else:
            self.data = self.data_origin

            self.data_lons = self.data.lon.values[0, :]
            self.data_lats = self.data.lat.values[:, 0]

            self.data = self.data.where(self.data['depth'] <= config["max_depth"], drop=True)
            self.rtofs_qc = self.rtofs_qc.where(self.rtofs_qc['depth'] <= config["max_depth"], drop=True)

    ### FUNCTION:
    def rtofs_save(self, config, directory):
        
        '''
        Save the subset RTOFS data as a NetCDF file.

        Args:
        - config (dict): Glider Guidance System mission configuration.
        - directory (str): Glider Guidance System mission directory.

        Returns:
        - None
        '''

        requested_datetime = self.data.attrs.get('requested_datetime', 'unknown_datetime')
        formatted_datetime = format_datetime(requested_datetime)

        rtofs_data_file = f"{config['glider_name']}_RTOFS_{config['max_depth']}m.nc"
        rtofs_data_path = os.path.join(directory, rtofs_data_file)
        self.data.to_netcdf(rtofs_data_path)

### FUNCTION:
def interp_depth_average(config, directory, model_data):
    
    '''
    Compute depth-averaged values, direction, and magnitude for the model data.
    Save the computational output and 1-meter bin averages.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - model_data (xarray.Dataset): Ocean model dataset.

    Returns:
    - calculated_data (xarray.Dataset): Dataset with the computed variables and layer information.
    - bin_data (xarray.Dataset): Dataset with the bin averages.
    '''

    # Extract the RTOFS dataset date
    model_datetime = model_data.attrs.get('requested_datetime', 'unknown_datetime')

    # Extracting data from the model dataset
    u_currents = model_data['u']
    v_currents = model_data['v']
    temperature = model_data['temperature']
    salinity = model_data['salinity']
    depths = model_data['depth'].values

    # Initialize arrays for depth-averaged and bin-averaged variables
    averaged_u = np.full((len(model_data.y), len(model_data.x)), np.nan)
    averaged_v = np.full_like(averaged_u, np.nan)
    averaged_temp = np.full_like(averaged_u, np.nan)
    averaged_sal = np.full_like(averaged_u, np.nan)
    averaged_direction = np.full_like(averaged_u, np.nan)
    averaged_magnitude = np.full_like(averaged_u, np.nan)

    # Determine the maximum number of bins across all grid points
    max_num_bins = int(np.nanmax([valid_depths.max() for valid_depths in depths if valid_depths.size > 0])) + 1

    # Initialize arrays to store bin averages for each variable at each depth and position
    bin_avg_u = np.full((len(model_data.y), len(model_data.x), max_num_bins), np.nan)
    bin_avg_v = np.full_like(bin_avg_u, np.nan)
    bin_avg_temp = np.full_like(bin_avg_u, np.nan)
    bin_avg_sal = np.full_like(bin_avg_u, np.nan)
    bin_avg_direction = np.full_like(bin_avg_u, np.nan)
    bin_avg_magnitude = np.full_like(bin_avg_u, np.nan)

    # Looping over each spatial point to interpolate and average the current data
    for y in range(len(model_data.y)):
        for x in range(len(model_data.x)):
            
            # Extract data values at the current point
            u_vals = u_currents.isel(y=y, x=x).values
            v_vals = v_currents.isel(y=y, x=x).values
            temp_vals = temperature.isel(y=y, x=x).values
            sal_vals = salinity.isel(y=y, x=x).values
            
            # Identifying valid depth indices
            valid_indices = ~np.isnan(u_vals) & ~np.isnan(v_vals)
            valid_depths = depths[valid_indices]

            # If there are valid depths, interpolate and average the data
            if valid_depths.size > 0:

                # Creating bins for each depth
                max_valid_depth = valid_depths.max()
                bins = np.arange(0, max_valid_depth + 1, 1)
                
                # Interpolation and averaging for each variable
                u_interp = interp1d(valid_depths, u_vals[valid_indices], bounds_error=False, fill_value="extrapolate")
                v_interp = interp1d(valid_depths, v_vals[valid_indices], bounds_error=False, fill_value="extrapolate")
                temp_interp = interp1d(valid_depths, temp_vals[valid_indices], bounds_error=False, fill_value="extrapolate")
                sal_interp = interp1d(valid_depths, sal_vals[valid_indices], bounds_error=False, fill_value="extrapolate")

                # Compute direction and magnitude at each depth bin
                for bin_idx, depth in enumerate(bins[:-1]):
                    u_at_depth = u_interp(depth)
                    v_at_depth = v_interp(depth)
                    direction_at_depth = (270 - np.rad2deg(np.arctan2(v_at_depth, u_at_depth)) + 180) % 360
                    bin_avg_direction[y, x, bin_idx] = direction_at_depth
                    magnitude_at_depth = np.sqrt(u_at_depth**2 + v_at_depth**2)
                    bin_avg_magnitude[y, x, bin_idx] = magnitude_at_depth

                # Storing bin averages for each variable
                for bin_idx, depth in enumerate(bins[:-1]):
                    bin_avg_u[y, x, bin_idx] = u_interp(depth).mean()
                    bin_avg_v[y, x, bin_idx] = v_interp(depth).mean()
                    bin_avg_temp[y, x, bin_idx] = temp_interp(depth).mean()
                    bin_avg_sal[y, x, bin_idx] = sal_interp(depth).mean()

                # Computing depth-averaged values for each variable
                total_depth = valid_depths.max()
                averaged_u[y, x] = np.nansum(bin_avg_u[y, x, :len(bins)-1]) / total_depth
                averaged_v[y, x] = np.nansum(bin_avg_v[y, x, :len(bins)-1]) / total_depth
                averaged_temp[y, x] = np.nansum(bin_avg_temp[y, x, :len(bins)-1]) / total_depth
                averaged_sal[y, x] = np.nansum(bin_avg_sal[y, x, :len(bins)-1]) / total_depth
                averaged_direction[y, x] = np.nansum(bin_avg_direction[y, x, :len(bins)-1]) / total_depth
                averaged_magnitude[y, x] = np.nansum(bin_avg_magnitude[y, x, :len(bins)-1]) / total_depth

    # Creating 'calculated_data' dataset
    depth_average_data = xr.Dataset({
        'u_avg': (('y', 'x'), averaged_u),
        'v_avg': (('y', 'x'), averaged_v),
        'temperature_avg': (('y', 'x'), averaged_temp),
        'salinity_avg': (('y', 'x'), averaged_sal),
        'magnitude_avg': (('y', 'x'), averaged_magnitude),
        'direction_avg': (('y', 'x'), averaged_direction)
    }, coords={'lat': model_data.lat, 'lon': model_data.lon})

    # Creating 'bin_data' dataset
    bin_average_data = xr.Dataset({
        'bin_avg_u': (('y', 'x', 'bin'), bin_avg_u),
        'bin_avg_v': (('y', 'x', 'bin'), bin_avg_v),
        'bin_avg_temp': (('y', 'x', 'bin'), bin_avg_temp),
        'bin_avg_sal': (('y', 'x', 'bin'), bin_avg_sal),
        'bin_avg_magnitude': (('y', 'x', 'bin'), bin_avg_magnitude),
        'bin_avg_direction': (('y', 'x', 'bin'), bin_avg_direction)
    }, coords={'lat': model_data.lat, 'lon': model_data.lon, 'bin': np.arange(max_num_bins)})

    # Get the filename datetime
    filename_datetime = get_filename_datetime(model_data)

    # Save 'depth_average_data' dataset as a NetCDF file
    depth_average_data_file = f"{config['glider_name']}_DepthAverageData.nc"
    depth_average_data_path = os.path.join(directory, depth_average_data_file)
    depth_average_data.to_netcdf(depth_average_data_path)

    # Save 'bin_average_data' dataset as a NetCDF file
    bin_average_data_file = f"{config['glider_name']}_BinAverageData.nc"
    bin_average_data_path = os.path.join(directory, bin_average_data_file)
    bin_average_data.to_netcdf(bin_average_data_path)

    return depth_average_data, bin_average_data
