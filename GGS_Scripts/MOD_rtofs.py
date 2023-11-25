# =========================
# X - IMPORTS
# =========================

### /// RTOFS ///
import numpy as np
import xarray as xr

# =========================
# [RTOFS] DATA PROCESSING
# =========================

### CLASS:
class RTOFS():

    ### FUNCTION:
    def __init__(self) -> None:
        
        '''
        Initialize the RTOFS instance and fetch initial RTOFS data.

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
        # rtofs_access = "rtofs_us_east_scraped.nc"

        try:
            rtofs_raw = xr.open_dataset(rtofs_access).set_coords(['lon', 'lat'])
            rtofs_raw.attrs['model'] = 'RTOFS'
            rtofs_raw= rtofs_raw.isel(time=-1)
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
            self.rtofs_qc = self.data.copy()

        else:
            self.data = self._data_orig

            self.data_lons = self.data.lon.values[0, :]
            self.data_lats = self.data.lat.values[:, 0]

            self.data = self.data.where(self.data['depth'] <= config["max_depth"], drop=True)
            self.rtofs_qc = self.data.copy()

### FUNCTION:
def model_currents(model_data, mask=False, mask_value=0.5):
    
    '''
    Calculate the depth-averaged currents for the given depth configuration.

    Args:
    - model_data (xarray.Dataset): ocean model data
    - mask (bool): whether or not to mask the data
        - default: False
        - if True, the data is masked based on the given mask_value
        - if False, the data is not masked
    - mask_value (float): the value to mask the data with
        - default: 0.5 (m/s)

    Returns:
    - u_avg (xarray.DataArray): depth-averaged u currents
    - v_avg (xarray.DataArray): depth-averaged v currents
    '''
    
    u_currents = model_data.data['u']
    v_currents = model_data.data['v']
    depths = model_data.data['depth']
    
    layer_thick = depths.diff(dim='depth', label='upper')

    n1_layer_thick = layer_thick.isel(depth=0)
    layer_thick = xr.concat([n1_layer_thick, layer_thick], dim='depth')
    layer_thick['depth'].values[0] = 0

    depth_weighted_u = u_currents * layer_thick
    depth_weighted_v = v_currents * layer_thick

    total_weighted_u = depth_weighted_u.sum(dim='depth')
    total_weighted_v = depth_weighted_v.sum(dim='depth')

    u_avg = total_weighted_u / layer_thick.sum(dim='depth')
    v_avg = total_weighted_v / layer_thick.sum(dim='depth')

    magnitude = np.sqrt(u_avg ** 2 + v_avg ** 2)

    if mask:
        u_avg = np.where(magnitude > mask_value, u_avg, np.nan)
        v_avg = np.where(magnitude > mask_value, v_avg, np.nan)
    
    return u_avg, v_avg, magnitude

# =========================
# X - MAIN
# =========================
# rtofs_data = RTOFS()
# rtofs_data.rtofs_subset(config, waypoints, subset=True)
# u_avg, v_avg, magnitude = model_currents(rtofs_data, mask=False)
# =========================

