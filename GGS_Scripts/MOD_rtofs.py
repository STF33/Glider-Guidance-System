# =========================
# X - IMPORTS
# =========================

### /// RTOFS ///
import numpy as np
import os
from scipy.interpolate import interp1d
import xarray as xr

# =========================
# [RTOFS] DATA PROCESSING
# =========================

### CLASS:
class RTOFS():
    
    '''
    Class for handling RTOFS data.

    Attributes:
    - data_orig (xarray.Dataset): original RTOFS data
    - data (xarray.Dataset): RTOFS data
    - x (np.array): x-coordinates of the RTOFS grid
    - y (np.array): y-coordinates of the RTOFS grid
    - grid_lons (np.array): longitudes of the RTOFS grid
    - grid_lats (np.array): latitudes of the RTOFS grid
    - rtofs_qc (xarray.Dataset): RTOFS data for quality control

    Methods:
    - __init__ (initialize the RTOFS instance)
    - rtofs_load (fetch the RTOFS data from the given URL and set its coordinates)
    - rtofs_subset (subset the RTOFS data based on the bounding box created by the given points)
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

        self._data_orig = self.rtofs_load()
        self._data_orig = self._data_orig.set_coords(['lat', 'lon'])
        self.data = self._data_orig.copy()

        self.x = self.data.x.values
        self.y = self.data.y.values
        self.grid_lons = self.data.lon.values[0,:]
        self.grid_lats = self.data.lat.values[:,0]

        self.rtofs_qc = self._data_orig.copy()

    ### FUNCTION:
    def rtofs_load(self):
        
        '''
        Fetch the RTOFS data from the given URL and set its coordinates.

        Args:
        - None

        Returns:
        - rtofs_raw (xarray.Dataset): RTOFS data
        '''

        rtofs_access = "https://tds.marine.rutgers.edu/thredds/dodsC/cool/rtofs/rtofs_us_east_scraped"

        try:
            rtofs_raw = xr.open_dataset(rtofs_access).set_coords(['lon', 'lat'])
            rtofs_raw.attrs['model'] = 'RTOFS'
            rtofs_raw = rtofs_raw.isel(time=-1)
            return rtofs_raw
        except Exception as e:
            print(f"Error fetching RTOFS data: {e}")
            return None
    
    ### FUNCTION:
    def rtofs_subset(self, config, waypoints, buffer=0.5, subset=True):

        '''
        Subset the RTOFS data based on the bounding box created by the given points.

        Args:
        - config (dict): Glider Guidance System configuration
        - waypoints (list): list of waypoints
        - buffer (float): buffer (in degrees) to add to the bounding box
            - default: 0.5
        - subset (bool): whether or not to subset the data
            - default: True
            - if True, the data is subset based on the bounding box of the waypoints
            - if False, the entire RTOFS grid is used (no subsetting)

        Returns:
        - None
        '''

        if subset:
            lats, lons = zip(*waypoints)

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

            self.data = self._data_orig.isel(
                x=slice(extent[0], extent[1]),
                y=slice(extent[2], extent[3])
            )

            self.data_lons = self.data.lon.values[0,:]
            self.data_lats = self.data.lat.values[:,0]
            
            self.data = self.data.where(self.data['depth'] <= config["max_depth"], drop=True)
            self.rtofs_qc = self.rtofs_qc.where(self.rtofs_qc['depth'] <= config["max_depth"], drop=True)

        else:
            self.data = self._data_orig

            self.data_lons = self.data.lon.values[0, :]
            self.data_lats = self.data.lat.values[:, 0]

            self.data = self.data.where(self.data['depth'] <= config["max_depth"], drop=True)
            self.rtofs_qc = self.rtofs_qc.where(self.rtofs_qc['depth'] <= config["max_depth"], drop=True)

    ### FUNCTION:
    def rtofs_save(self, config, directory):
        
        '''
        Save the subset RTOFS data as a NetCDF file.

        Args:
        - config (dict): Glider Guidance System configuration
        - directory (str): directory to save the file

        Returns:
        - None
        '''

        rtofs_data_file = f"{config['glider_name']}_rtofs_{config['max_depth']}m_.nc"
        rtofs_data_path = os.path.join(directory, rtofs_data_file)
        self.data.to_netcdf(rtofs_data_path)

### FUNCTION:
def interp_average(config, directory, model_data):
    
    '''
    Compute depth-averaged ocean currents, temperature, and salinity, 
    and store 1-meter bin averages for the input model data.

    Args:
    - config (dict): Glider Guidance System configuration
    - directory (str): directory to save the file
    - model_data (xarray.Dataset): Ocean model data

    Returns:
    - calculated_data (xarray.Dataset): dataset with the computed variables and layer information
    - bin_data (xarray.Dataset): dataset with the bin averages
    '''

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

    # Determine the maximum number of bins across all grid points
    max_num_bins = int(np.nanmax([valid_depths.max() for valid_depths in depths if valid_depths.size > 0])) + 1

    # Initialize arrays to store bin averages for each variable at each depth and position
    bin_avg_u = np.full((len(model_data.y), len(model_data.x), max_num_bins), np.nan)
    bin_avg_v = np.full((len(model_data.y), len(model_data.x), max_num_bins), np.nan)
    bin_avg_temp = np.full((len(model_data.y), len(model_data.x), max_num_bins), np.nan)
    bin_avg_sal = np.full((len(model_data.y), len(model_data.x), max_num_bins), np.nan)

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
                
                # Computing bin averages for each variable
                bin_averages_u = [u_interp(depth).mean() for depth in bins[:-1]]
                bin_averages_v = [v_interp(depth).mean() for depth in bins[:-1]]
                bin_averages_temp = [temp_interp(depth).mean() for depth in bins[:-1]]
                bin_averages_sal = [sal_interp(depth).mean() for depth in bins[:-1]]
                
                # Computing depth-averaged values at the XY gridpoint
                averaged_u[y, x] = np.mean(bin_averages_u)
                averaged_v[y, x] = np.mean(bin_averages_v)
                averaged_temp[y, x] = np.mean(bin_averages_temp)
                averaged_sal[y, x] = np.mean(bin_averages_sal)

                # Storing bin averages for each variable
                for i in range(len(bins) - 1):
                    bin_avg_u[y, x, i] = bin_averages_u[i]
                    bin_avg_v[y, x, i] = bin_averages_v[i]
                    bin_avg_temp[y, x, i] = bin_averages_temp[i]
                    bin_avg_sal[y, x, i] = bin_averages_sal[i]

    # Compute the depth-averaged current magnitude
    magnitude = np.sqrt(averaged_u ** 2 + averaged_v ** 2)

    # Creating 'calculated_data' dataset
    calculated_data = xr.Dataset({
        'u_avg': (('y', 'x'), averaged_u),
        'v_avg': (('y', 'x'), averaged_v),
        'temperature_avg': (('y', 'x'), averaged_temp),
        'salinity_avg': (('y', 'x'), averaged_sal),
        'magnitude': (('y', 'x'), magnitude)
    }, coords={'lat': model_data.lat, 'lon': model_data.lon})

    # Creating 'bin_data' dataset
    bin_data = xr.Dataset({
        'bin_avg_u': (('y', 'x', 'bin'), bin_avg_u),
        'bin_avg_v': (('y', 'x', 'bin'), bin_avg_v),
        'bin_avg_temp': (('y', 'x', 'bin'), bin_avg_temp),
        'bin_avg_sal': (('y', 'x', 'bin'), bin_avg_sal),
    }, coords={'lat': model_data.lat, 'lon': model_data.lon, 'bin': np.arange(max_num_bins)})

    # Save 'calculated_data' dataset as a NetCDF file
    calculated_data_file = f"{config['glider_name']}_calculated_data.nc"
    calculated_data_path = os.path.join(directory, calculated_data_file)
    calculated_data.to_netcdf(calculated_data_path)

    # Save 'bin_data' dataset as a NetCDF file
    bin_data_file = f"{config['glider_name']}_bin_data.nc"
    bin_data_path = os.path.join(directory, bin_data_file)
    bin_data.to_netcdf(bin_data_path)

    return calculated_data, bin_data

# =========================
#
# rtofs = RTOFS()
# rtofs.rtofs_subset(config, waypoints, subset=True)
# rtofs_data = rtofs.data
# rtofs_qc = rtofs.rtofs_qc
# rtofs.rtofs_save(config, directory)
#
# calculated_data, bin_data = interp_average(config, directory, rtofs_data)
#
# =========================
