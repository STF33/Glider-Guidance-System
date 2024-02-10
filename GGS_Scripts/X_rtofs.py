# =========================
# X - IMPORTS
# =========================

import numpy as np
import os
import pandas as pd
import xarray as xr

# =========================

### CLASS:
class RTOFS():
    
    '''
    Class for handling RTOFS data.

    Attributes:
    - data_origin (xarray.Dataset): Original RTOFS data.
    - data (xarray.Dataset): Subset or full RTOFS data.
    - x (np.ndarray): RTOFS x-axis grid.
    - y (np.ndarray): RTOFS y-axis grid.
    - grid_lons (np.ndarray): RTOFS longitude grid.
    - grid_lats (np.ndarray): RTOFS latitude grid.
    - rtofs_qc (xarray.Dataset): RTOFS data for quality control.
    '''

    ### FUNCTION:
    def __init__(self) -> None:
        
        '''
        Initialize the RTOFS instance.

        Args:
        - datetime_index (datetime.datetime): The datetime index for data loading.

        Returns:
        - None
        '''

        print("\n### MODEL DATA: [RTOFS] ###\n")

        self.data_origin = None
        self.x = None
        self.y = None
        self.grid_lons = None
        self.grid_lats = None

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
            if rtofs_raw:
                datetime = pd.Timestamp(datetime_index).tz_localize(None)
                time_values = rtofs_raw.time.values
                time_index = np.argmin(np.abs(time_values - np.datetime64(datetime)))
                rtofs_raw = rtofs_raw.isel(time=time_index)

                self.data_origin = rtofs_raw
                self.x = self.data_origin.x.values
                self.y = self.data_origin.y.values
                self.grid_lons = self.data_origin.lon.values[0,:]
                self.grid_lats = self.data_origin.lat.values[:,0]
                self.data_origin.attrs['model_datetime'] = str(rtofs_raw.time.values)

                print("MODEL DATA ACQUIRED")
                print(f"Requested datetime: {datetime_index}")
                print(f"Nearest datetime index in the dataset: {time_index}")
                print(f"Acquired datetime: {rtofs_raw.attrs['model_datetime']}")
            else:
                print("Failed to load RTOFS data. Dataset is None.")
        except Exception as e:
            print(f"Error fetching RTOFS data: {e}")
    
    ### FUNCTION:
    def rtofs_subset(self, config, subset=True):

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
        
        max_depth = config["max_depth"]
        depth_indices = self.data_origin.depth.values
        target_depth_index = depth_indices[depth_indices >= max_depth][0]

        if subset:
            lats, lons = zip(*config['extent'])
            min_lon, max_lon = min(lons), max(lons)
            min_lat, max_lat = min(lats), max(lats)

            lons_idx = np.interp([min_lon, max_lon], self.grid_lons, self.x)
            lats_idx = np.interp([min_lat, max_lat], self.grid_lats, self.y)

            extent = [
                np.floor(lons_idx[0]).astype(int),
                np.ceil(lons_idx[1]).astype(int),
                np.floor(lats_idx[0]).astype(int),
                np.ceil(lats_idx[1]).astype(int)
            ]

            self.data_origin = self.data_origin.isel(
                x=slice(extent[0], extent[1]),
                y=slice(extent[2], extent[3])
            ).sel(depth=slice(0, target_depth_index))

        else:
            self.data_origin = self.data_origin.sel(depth=slice(0, target_depth_index))
        
    ### FUNCTION:
    def rtofs_save(self, config, directory, save_data=True, save_qc=True):
        
        '''
        Save the subset RTOFS data as a NetCDF file.

        Args:
        - config (dict): Glider Guidance System mission configuration.
        - directory (str): Glider Guidance System mission directory.

        Returns:
        - None
        '''

        self.data = self.data_origin.copy()
        self.qc = self.data_origin.copy()

        if save_data:
            rtofs_data_file = f"{config['glider_name']}_RTOFS_Data_{config['max_depth']}m.nc"
            rtofs_data_path = os.path.join(directory, rtofs_data_file)
            self.data.to_netcdf(rtofs_data_path)
            print(f"RTOFS Data saved to: {rtofs_data_path}")
        
        if save_qc:
            rtofs_qc_file = f"{config['glider_name']}_RTOFS_QC_{config['max_depth']}m.nc"
            rtofs_qc_path = os.path.join(directory, rtofs_qc_file)
            self.qc.to_netcdf(rtofs_qc_path)
            print(f"RTOFS qc saved to: {rtofs_qc_path}")
