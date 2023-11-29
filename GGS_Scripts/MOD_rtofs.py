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
def compute_currents(model_data):
    
    '''
    Calculate the thickness-weighted depth-averaged currents for the model data.
    Create an xarray dataset with the computed variables and layer information.

    Args:
    - model_data (xarray.Dataset): Ocean model data

    Returns:
    - u_avg (xarray.DataArray): depth-averaged zonal currents
    - v_avg (xarray.DataArray): depth-averaged meridional currents
    - magnitude (xarray.DataArray): depth-averaged current magnitude
    - currents_data (xarray.Dataset): dataset with the computed variables and layer information
    '''

    depths = model_data['depth']
    u_currents = model_data['u']
    v_currents = model_data['v']
    
    layer_bins = depths.diff(dim='depth', label='upper')
    layer_n1 = layer_bins.isel(depth=0)
    layer_bins = xr.concat([layer_n1, layer_bins], dim='depth')
    layer_bins['depth'].values[0] = 0

    z_weighted_u = u_currents * layer_bins
    z_weighted_v = v_currents * layer_bins

    sum_weighted_u = z_weighted_u.sum(dim='depth')
    sum_weighted_v = z_weighted_v.sum(dim='depth')
    sum_layer_bins = layer_bins.sum(dim='depth')

    u_avg = sum_weighted_u / sum_layer_bins
    v_avg = sum_weighted_v / sum_layer_bins
    magnitude = np.sqrt(u_avg ** 2 + v_avg ** 2)

    currents_data = xr.Dataset({
        'u_avg': u_avg,
        'v_avg': v_avg,
        'magnitude': magnitude,
        'layer_thickness': layer_bins,
        'z_weighted_u': z_weighted_u.sum(dim='depth'),
        'z_weighted_v': z_weighted_v.sum(dim='depth')
        }, coords=model_data.coords)

    return u_avg, v_avg, magnitude, currents_data

# =========================
# X - MAIN
# =========================
# rtofs = RTOFS()
# rtofs.rtofs_subset(config, waypoints, subset=True)
#
# rtofs_data = rtofs.data
# rtofs_qc = rtofs.rtofs_qc
#
# u_avg, v_avg, magnitude, currents_data = compute_currents(rtofs_data)
# =========================
