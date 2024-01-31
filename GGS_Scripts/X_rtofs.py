# =========================
# X - IMPORTS
# =========================

import numpy as np
import os
import pandas as pd
import xarray as xr
from X_functions import datetime_format

# =========================
# [RTOFS] DATA PROCESSING
# =========================

### CLASS:
# @profile
class RTOFS():
    
    '''
    Class for handling RTOFS data.

    Attributes:
    - data_origin (xarray.Dataset): Original RTOFS data.
    - data (xarray.Dataset): Subset RTOFS data.
    - x (np.ndarray): RTOFS x-axis grid.
    - y (np.ndarray): RTOFS y-axis grid.
    - grid_lons (np.ndarray): RTOFS longitude grid.
    - grid_lats (np.ndarray): RTOFS latitude grid.
    - rtofs_qc (xarray.Dataset): RTOFS data for quality control.

    Methods:
    - __init__ (datetime_index=None): Initialize the RTOFS instance.
    - rtofs_load (datetime_index): Fetch the RTOFS data from the given URL source and set its coordinates.
    - rtofs_subset (config, buffer=0, subset=True): Subset the RTOFS data based on the GGS mission extent.
    - rtofs_save (config, directory): Save the subset RTOFS data as a NetCDF file.
    '''

    ### FUNCTION:
    def __init__(self, datetime_index=None) -> None:
        
        '''
        Initialize the RTOFS instance.

        Args:
        - datetime_index (datetime.datetime): The datetime index for data loading.

        Returns:
        - None
        '''

        print("\n")
        print("### MODEL DATA: [RTOFS] ###")
        print("\n")

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

            print("MODEL DATA ACQUIRED")
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
            
            self.data = self.data.sel(depth=slice(None, config["max_depth"]))
            self.rtofs_qc = self.rtofs_qc.sel(depth=slice(None, config["max_depth"]))

        else:

            self.data = self.data_origin

            self.data_lons = self.data.lon.values[0, :]
            self.data_lats = self.data.lat.values[:, 0]

            self.data = self.data.sel(depth=slice(None, config["max_depth"]))
            self.rtofs_qc = self.rtofs_qc.sel(depth=slice(None, config["max_depth"]))

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
        formatted_datetime = datetime_format(requested_datetime)

        rtofs_data_file = f"{config['glider_name']}_RTOFS_{config['max_depth']}m.nc"
        rtofs_data_path = os.path.join(directory, rtofs_data_file)
        self.data.to_netcdf(rtofs_data_path)
