# =========================
# IMPORTS
# =========================

import copernicusmarine as cm
from dateutil import parser
import numpy as np
import os
import pandas as pd
import xarray as xr

# =========================

### CLASS:
class RTOFS():
    
    '''
    Class for handling RTOFS data.
    '''

    ### FUNCTION:
    def __init__(self) -> None:
        
        '''
        Initialize the RTOFS instance.

        Args:
        - None

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
    def rtofs_load(self, config, datetime_index, subset_extent=True, subset_depth=True):
        
        '''
        Fetch the RTOFS ocean model data and standardize it.

        Args:
        - config (dict): Glider Guidance System mission configuration.
        - datetime_index (str): Index of the datetime to fetch.
        - subset_extent (bool): Subset the data based on the mission extent.
            - default: 'True'
        - subset_depth (bool): Subset the data based on the mission depth.
            - default: 'True'

        Returns:
        - None
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
                self.data_origin.attrs['model_name'] = 'RTOFS'

                if subset_extent and config:
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
                    )

                if subset_depth and config:
                    max_depth = config["max_depth"]
                    depth_indices = self.data_origin.depth.values
                    target_depth_index = depth_indices[depth_indices >= max_depth][0]
                    self.data_origin = self.data_origin.sel(depth=slice(0, target_depth_index))
            else:
                print("Failed to load RTOFS data. Dataset is None.")
        except Exception as e:
            print(f"Error fetching RTOFS data: {e}")
    
    ### FUNCTION:
    def rtofs_save(self, config, directory, save_data=True, save_qc=True):
        
        '''
        Save the subset RTOFS data as a NetCDF file.

        Args:
        - config (dict): Glider Guidance System mission configuration.
        - directory (str): Glider Guidance System mission directory.
        - save_data (bool): Save the data file.
            - default: 'True'
        - save_qc (bool): Save the quality control file.
            - default: 'True'

        Returns:
        - None
        '''

        self.data = self.data_origin.copy()
        self.qc = self.data_origin.copy()

        if save_data:
            rtofs_data_file = f"RTOFS_Data_{config['max_depth']}m.nc"
            rtofs_data_path = os.path.join(directory, rtofs_data_file)
            self.data.to_netcdf(rtofs_data_path)
            print(f"RTOFS Data saved to: {rtofs_data_path}")
        
        if save_qc:
            rtofs_qc_file = f"RTOFS_QC_{config['max_depth']}m.nc"
            rtofs_qc_path = os.path.join(directory, rtofs_qc_file)
            self.qc.to_netcdf(rtofs_qc_path)
            print(f"RTOFS QC saved to: {rtofs_qc_path}")

### CLASS:
class CMEMS:
    
    '''
    Class for handling CMEMS data.
    '''
    
    ### FUNCTION:
    def __init__(self, username, password) -> None:
        
        '''
        Initialize the CMEMS instance.

        Args:
        - username (str): CMEMS username.
        - password (str): CMEMS password.

        Returns:
        - None
        '''

        print("\n### MODEL DATA: [CMEMS] ###\n")

        self.username = username
        self.password = password
        self.data_origin = None

    ### FUNCTION:
    def cmems_load(self, config, datetime_index, subset_extent=True, subset_depth=True):
        
        '''
        Fetch the CMEMS ocean model data and standardize it.

        Args:
        - config (dict): Glider Guidance System mission configuration.
        - datetime_index (str): Index of the datetime to fetch.
        - subset_extent (bool): Subset the data based on the mission extent.
            - default: 'True'
        - subset_depth (bool): Subset the data based on the mission depth.
            - default: 'True'

        Returns:
        - None
        '''
        
        datetime_index = parser.parse(datetime_index)
        formatted_datetime_index = datetime_index.strftime('%Y-%m-%dT%H:%M:%S')
        start_datetime = end_datetime = formatted_datetime_index

        lats, lons = zip(*config['extent'])
        min_lon, max_lon = min(lons), max(lons)
        min_lat, max_lat = min(lats), max(lats)

        self.data_origin = cm.open_dataset(
            dataset_id="cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i",
            minimum_longitude=min_lon if subset_extent else None,
            maximum_longitude=max_lon if subset_extent else None,
            minimum_latitude=min_lat if subset_extent else None,
            maximum_latitude=max_lat if subset_extent else None,
            start_datetime=start_datetime,
            end_datetime=end_datetime,
            variables=["uo", "vo"],
            username=self.username,
            password=self.password
        )

        self.data_origin.attrs['model_datetime'] = datetime_index.strftime('%Y-%m-%dT%H:%M:%S')
        self.data_origin.attrs['model_name'] = 'CMEMS'

        if subset_depth and 'depth' in self.data_origin.coords:
            max_depth = config["max_depth"]
            depth_indices = self.data_origin.depth.values
            target_depth_index = np.searchsorted(depth_indices, max_depth, side='right') - 1
            self.data_origin = self.data_origin.isel(depth=slice(None, target_depth_index + 1))
        
        rename_dict = {'uo': 'u', 'vo': 'v', 'latitude': 'lat', 'longitude': 'lon'}
        existing_vars = set(self.data_origin.variables.keys()) & set(rename_dict.keys())
        final_rename_dict = {k: rename_dict[k] for k in existing_vars}
        self.data_origin = self.data_origin.rename(final_rename_dict)

    ### FUNCTION:
    def cmems_save(self, config, directory, save_data=True, save_qc=True):
        
        '''
        Save the subset CMEMS data as a NetCDF file.

        Args:
        - config (dict): Glider Guidance System mission configuration.
        - directory (str): Glider Guidance System mission directory.
        - save_data (bool): Save the data file.
            - default: 'True'
        - save_qc (bool): Save the quality control file.
            - default: 'True'

        Returns:
        - None
        '''
        
        self.data = self.data_origin.copy()
        self.qc = self.data_origin.copy()

        if save_data:
            cmems_data_file = f"CMEMS_Data_{config['max_depth']}m.nc"
            cmems_data_path = os.path.join(directory, cmems_data_file)
            self.data.to_netcdf(cmems_data_path)
            print(f"CMEMS Data saved to: {cmems_data_path}")
        
        if save_qc:
            cmems_qc_file = f"CMEMS_QC_{config['max_depth']}m.nc"
            cmems_qc_path = os.path.join(directory, cmems_qc_file)
            self.qc.to_netcdf(cmems_qc_path)
            print(f"CMEMS QC saved to: {cmems_qc_path}")

### CLASS:
class GOFS:

    '''
    Class for handling GOFS data.
    '''
    
    ### FUNCTION:
    def __init__(self) -> None:
        
        '''
        Initialize the GOFS instance.

        Args:
        - None

        Returns:
        - None
        '''

        print("\n### MODEL DATA: [GOFS] ###\n")
        
        self.data_origin = None

    ### FUNCTION:
    def gofs_load(self, config, datetime_index, subset_extent=True, subset_depth=True):
        
        '''
        Fetch the GOFS ocean model data and standardize it.

        Args:
        - config (dict): Glider Guidance System mission configuration.
        - datetime_index (str): Index of the datetime to fetch.
        - subset_extent (bool): Subset the data based on the mission extent.
            - default: 'True'
        - subset_depth (bool): Subset the data based on the mission depth.
            - default: 'True'

        Returns:
        - None
        '''

        datetime_index = parser.parse(datetime_index)
        datetime_index = datetime_index.replace(tzinfo=None)

        gofs_access = "https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0"

        try:
            gofs_raw = xr.open_dataset(gofs_access, drop_variables="tau")
            if gofs_raw:
                gofs_raw = gofs_raw.sel(time=datetime_index, method='nearest')

                gofs_raw['lon'] = ((gofs_raw['lon'] + 180) % 360) - 180
                gofs_raw = gofs_raw.sortby(gofs_raw['lon'])

                self.data_origin = gofs_raw

                self.data_origin.attrs['model_datetime'] = datetime_index.strftime('%Y-%m-%dT%H:%M:%S')
                self.data_origin.attrs['model_name'] = 'GOFS'

                if subset_extent:
                    lats, lons = zip(*config['extent'])
                    min_lon, max_lon = min(lons), max(lons)
                    min_lat, max_lat = min(lats), max(lats)

                    self.data_origin = self.data_origin.sel(lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))

                if subset_depth and 'depth' in self.data_origin.coords:
                    max_depth = config["max_depth"]
                    depth_indices = self.data_origin.depth.values
                    target_depth_index = np.searchsorted(depth_indices, max_depth, side='right') - 1
                    
                    self.data_origin = self.data_origin.isel(depth=slice(None, target_depth_index + 1))

                rename_dict = {
                    "surf_el": "sea_surface_height",
                    "water_temp": "temperature",
                    "water_u": "u",
                    "water_v": "v"
                }
                existing_vars = set(self.data_origin.variables.keys()) & set(rename_dict.keys())
                final_rename_dict = {k: rename_dict[k] for k in existing_vars}
                self.data_origin = self.data_origin.rename(final_rename_dict)

            else:
                print("Failed to load GOFS data. Dataset is None.")
        except Exception as e:
            print(f"Error fetching GOFS data: {e}")

    ### FUNCTION:
    def gofs_save(self, config, directory, save_data=True, save_qc=True):

        '''
        Save the subset GOFS data as a NetCDF file.

        Args:
        - config (dict): Glider Guidance System mission configuration.
        - directory (str): Glider Guidance System mission directory.
        - save_data (bool): Save the data file.
            - default: 'True'
        - save_qc (bool): Save the quality control file.
            - default: 'True'
        
        Returns:
        - None
        '''   
        
        self.data = self.data_origin.copy()
        self.qc = self.data_origin.copy()
        
        if save_data:
            gofs_data_file = f"GOFS_Data_{config['max_depth']}m.nc"
            gofs_data_path = os.path.join(directory, gofs_data_file)
            self.data.to_netcdf(gofs_data_path)
            print(f"GOFS Data saved to: {gofs_data_path}")
        
        if save_qc:
            gofs_qc_file = f"GOFS_QC_{config['max_depth']}m.nc"
            gofs_qc_path = os.path.join(directory, gofs_qc_file)
            self.qc.to_netcdf(gofs_qc_path)
            print(f"GOFS QC saved to: {gofs_qc_path}")
