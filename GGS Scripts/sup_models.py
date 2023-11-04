# =========================
# X - Imports
# =========================

import numpy as np
import xarray as xr

# =========================
# RTOFS
# =========================

def rtofs():
    url = "https://tds.marine.rutgers.edu/thredds/dodsC/cool/rtofs/rtofs_us_east_scraped"
    ds = xr.open_dataset(url).set_coords(['lon', 'lat'])
    ds.attrs['model'] = 'RTOFS'
    return ds

class RTOFS():
    def __init__(self) -> None:
        self._data_orig = self.load()
        self._data_orig = self._data_orig.set_coords(['u', 'v'])
        self.data = self._data_orig.copy()
        self.x = self.data.x.values
        self.y = self.data.y.values
        self.grid_lons = self.data.lon.values[0,:]
        self.grid_lats = self.data.lat.values[:,0]

    def load(self):
        url = "https://tds.marine.rutgers.edu/thredds/dodsC/cool/rtofs/rtofs_us_east_scraped"
        ds = xr.open_dataset(url).set_coords(['lon', 'lat'])
        ds.attrs['model'] = 'RTOFS'
        return ds

    def subset(self, extent):
         # Find x, y indexes of the area we want to subset
        lons_ind = np.interp(extent[:2], self.grid_lons, self.x)
        lats_ind = np.interp(extent[2:], self.grid_lats, self.y)

        # Use np.floor on the first index and np.ceiling on  the second index
        # of each slice in order to widen the area of the extent slightly.
        extent = [
            np.floor(lons_ind[0]).astype(int),
            np.ceil(lons_ind[1]).astype(int),
            np.floor(lats_ind[0]).astype(int),
            np.ceil(lats_ind[1]).astype(int)
        ]

        # Use the xarray .isel selector on x/y 
        # since we know the exact indexes we want to slice
        self.data = self._data_orig.isel(
            x=slice(extent[0], extent[1]),
            y=slice(extent[2], extent[3])
            )

    def transect(self, start, end, grid_spacing=5000):
        """
        # Return a transect 

        Args:
            start (array): Start of transect (lon, lat)
            end (array): End of transect (lon, lat)
            grid_spacing (int, optional): Distance (meters) between each point along transect. Defaults to 5000.

        Returns:
            _type_: _description_
        """
        return calculate_transect(start, end, grid_spacing)

    def profile(self, lon, lat, method='nearest'):
        # Find x, y indexes of the area we want to subset
        lons_ind = np.interp(lon, self.grid_lons, self.x)
        lats_ind = np.interp(lat, self.grid_lats, self.y)

        if method == 'nearest':
            pass
        elif method == 'interp':
            pass

# =========================
# CMEMS
# =========================

def cmems(rename=False):
    username = 'maristizabalvar'
    password = 'MariaCMEMS2018'
    from pydap.client import open_url
    from pydap.cas.get_cookies import setup_session
    cas_url = 'https://cmems-cas.cls.fr/cas/login'
    session = setup_session(cas_url, username, password)
    session.cookies.set("CASTGC", session.cookies.get_dict()['CASTGC'])

    try:
        url = f'https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i'
        data_store1 = xr.backends.PydapDataStore(open_url(url, session=session, user_charset='utf-8')) # needs PyDAP >= v3.3.0 see https://github.com/pydap/pydap/pull/223/commits 
    except:
        url = f'https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i'
        data_store1 = xr.backends.PydapDataStore(open_url(url, session=session, user_charset='utf-8')) # needs PyDAP >= v3.3.0 see https://github.com/pydap/pydap/pull/223/commits

    try:
        url = f'https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i'
        data_store2 = xr.backends.PydapDataStore(open_url(url, session=session, user_charset='utf-8')) # needs PyDAP >= v3.3.0 see https://github.com/pydap/pydap/pull/223/commits 
    except:
        url = f'https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i'
        data_store2 = xr.backends.PydapDataStore(open_url(url, session=session, user_charset='utf-8')) # needs PyDAP >= v3.3.0 see https://github.com/pydap/pydap/pull/223/commits

    try:
        url = f'https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy-so_anfc_0.083deg_PT6H-i'
        data_store3 = xr.backends.PydapDataStore(open_url(url, session=session, user_charset='utf-8')) # needs PyDAP >= v3.3.0 see https://github.com/pydap/pydap/pull/223/commits 
    except:
        url = f'https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy-so_anfc_0.083deg_PT6H-i'
        data_store3 = xr.backends.PydapDataStore(open_url(url, session=session, user_charset='utf-8')) # needs PyDAP >= v3.3.0 see https://github.com/pydap/pydap/pull/223/commits
        

    # # Downloading and reading Copernicus grid
    ds1 = xr.open_dataset(data_store1, drop_variables='tau')
    ds2 = xr.open_dataset(data_store2, drop_variables='tau')
    ds3 = xr.open_dataset(data_store3, drop_variables='tau')
    ds = xr.merge([ds1, ds2, ds3])

    ds.attrs['model'] = 'CMEMS'

    if rename:
        ds = ds.rename(
            {
                'thetao': 'temperature', 
                'so': 'salinity',
                'latitude': 'lat',
                'longitude': 'lon',
                'uo': 'u',
                'vo': 'v'
                }
            )
    return ds