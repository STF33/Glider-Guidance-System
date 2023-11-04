import datetime as dt
import inspect
import multiprocessing
from collections import namedtuple
from pprint import pprint
# urllib.error.HTTPError
from urllib.error import HTTPError as uHTTPError
from urllib.error import URLError

import pandas as pd
from erddapy import ERDDAP
from joblib import Parallel, delayed
from numpy import isin
from requests.exceptions import HTTPError as rHTTPError
import requests 

Argo = namedtuple('Argo', ['name', 'lon', 'lat'])
Glider = namedtuple('Glider', ['name', 'lon', 'lat'])
time_formatter = '%Y-%m-%dT%H:%M:%SZ'

rename_gliders = {}
# rename_gliders["time (UTC)"] = "time"
rename_gliders["longitude"] = "lon"
rename_gliders["latitude"] = "lat"

rename_argo = {}
rename_argo["platform_number"] = "argo"
rename_argo["time (UTC)"] = "time"
rename_argo["longitude (degrees_east)"] = "lon"
rename_argo["latitude (degrees_north)"] = "lat"


def get_argo_floats_by_time(bbox=(-100, -45, 5, 46),
                            time_start=None, time_end=dt.date.today(),
                            wmo_id=None, variables=None):
    """_summary_

    Args:
        bbox (_type_, optional): _description_. Defaults to None.
        time_start (_type_, optional): Start time. Defaults to None.
        time_end (_type_, optional): End time. Defaults to dt.date.today().
        floats (_type_, optional): _description_. Defaults to None.
        add_vars (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    # Accept both tuples and lists, but turn tuple into a list if passed a tuple
    if isinstance(variables, tuple):
        variables = list(variables)
        
    time_start = time_start or (time_end - dt.timedelta(days=1))
    
    default_variables = ['platform_number', 'time', 'longitude', 'latitude']
    
    constraints = {
        'time>=': str(time_start),
        'time<=': str(time_end),
    }

    if bbox:
        constraints['longitude>='] = bbox[0]
        constraints['longitude<='] = bbox[1]
        constraints['latitude>='] = bbox[2]
        constraints['latitude<='] = bbox[3]

    if wmo_id:
        if isinstance(wmo_id, int) or isinstance(wmo_id, float):
            wmo_id = str(wmo_id)
            
        constraints['platform_number='] = wmo_id

    if variables:
        default_variables = default_variables + variables
        default_variables = list(set(default_variables)) # remove duplicates
        
    e = ERDDAP(
        server='IFREMER',
        protocol='tabledap',
        response='csv'
    )

    e.dataset_id = 'ArgoFloats'
    e.constraints = constraints
    e.variables = default_variables

    try:
        df = e.to_pandas(
            index_col="time (UTC)",
            parse_dates=True,
        ).dropna().tz_localize(None)
        df = df.reset_index().rename(rename_argo, axis=1)
        df = df.set_index(["argo", "time"]).sort_index()
    except rHTTPError:
        df = pd.DataFrame()
    return df


def get_active_gliders(bbox=None, t0=None, t1=dt.date.today(), variables=None, 
                       timeout=5, parallel=False):
    variables = variables or ['time', 'latitude', 'longitude']
    bbox = bbox or [-100, -40, 18, 60]
    t0 = t0 or (t1 - dt.timedelta(days=1))

    # Convert dates to strings
    t0 = t0.strftime('%Y-%m-%dT%H:%M:%SZ')
    t1 = t1.strftime('%Y-%m-%dT%H:%M:%SZ')

    # Initialize GliderDAC Object
    e = ERDDAP(server='NGDAC')

    # Set timeout (seconds)
    e.requests_kwargs['timeout'] = timeout

    # Grab every dataset available
    # Search constraints
    kw = dict()
    kw['min_time'] = t0
    kw['max_time'] = t1

    if bbox:
        kw['min_lon'] = bbox[0]
        kw['max_lon'] = bbox[1]
        kw['min_lat'] = bbox[2]
        kw['max_lat'] = bbox[3]

    search_url = e.get_search_url(search_for=None, response='csv', **kw)

    try:
        # Grab the results
        search = pd.read_csv(search_url)
    except uHTTPError as error:
        print(f"{inspect.currentframe().f_code.co_name} - Error: {error}")
        # return empty dataframe if there are no results
        return pd.DataFrame()
    except URLError as e:
        print(f"{inspect.currentframe().f_code.co_name} - Error: {error}")
        # return empty dataframe if there are no results
        return pd.DataFrame()

    # Extract the IDs
    gliders = search['Dataset ID'].values

    msg = f"Found {len(gliders)} Glider Datasets: "
    pprint(msg + ', '.join(gliders.tolist()))

    # Setting constraints
    constraints = {
            'time>=': t0,
            'time<=': t1,
            # 'longitude>=': bbox[0],
            # 'longitude<=': bbox[1],
            # 'latitude>=': bbox[2],
            # 'latitude<=': bbox[3],
            }
    

    def request_multi(dataset_id, protocol="tabledap"): # variables=None):
        # variables = variables or ['depth', 'latitude', 'longitude']
        e.constraints = constraints
        e.protocol = protocol
        # e.variables = ['time', 'latitude', 'longitude']
        e.variables = variables
        e.dataset_id = dataset_id
        
        # Drop units in the first line and Nans
        try:
            df = e.to_pandas(
                response="csv", 
                index_col="time",
                parse_dates=True,
                skiprows=(1,)
                ).dropna().tz_localize(None)
        except rHTTPError:
            df = pd.DataFrame()
        return (dataset_id, df)

    # If we want to take advantage of parallel processing,
    # parellel = True as an optional argument input.
    if parallel:
        if isinstance(parallel, bool):
            num_cores = multiprocessing.cpu_count()
        else:
            num_cores = parallel
        
        downloads = Parallel(n_jobs=num_cores)(
            delayed(request_multi)(dataset_id) for dataset_id in gliders
        )
        dfs = {glider: df for (glider, df) in downloads}
    else:
        dfs = {glider: df for (glider, df) in [request_multi(id) for id in gliders]}

    try:
        df = pd.concat(dfs)
        df.index.names = ["glider", "time"]
        df = df.sort_index().rename(rename_gliders, axis=1)
    except ValueError:
        df = pd.DataFrame()
    return df


def get_glider_by_id(dataset_id=None, bbox=None, start=None, end=None, vars=None):
    """_summary_

    Args:
        dataset_id (_type_, optional): _description_. Defaults to None.
        bbox (_type_, optional): _description_. Defaults to None.
        start (_type_, optional): _description_. Defaults to None.
        end (_type_, optional): _description_. Defaults to None.
        vars (_type_, optional): _description_. Defaults to None.

    Raises:
        TypeError: _description_

    Returns:
        _type_: _description_
    """
    if dataset_id is None:
        raise TypeError('glider_id must be a string containing a dataset id from https://gliders.ioos.us/')
    

    variables = [
        "time",
        "longitude",
        "latitude",
        "pressure",
        "depth",
        "temperature",
        "salinity",
        "conductivity",
        "density"
        ]
      
    e = ERDDAP(
        server='NGDAC',
        protocol="tabledap")

    constraints = {}

    if bbox:
        constraints['longitude>='] =  bbox[0]
        constraints['longitude<='] = bbox[1]
        constraints['latitude>='] = bbox[2]
        constraints['latitude<='] = bbox[3]

    if start:
        constraints["time>="] =  start
        constraints["time<="] =  end
        
    
    if constraints: 
        e.constraints = constraints
    
    if vars:
        variables = variables + vars
        
    # print('Reading ' + id)
    e.dataset_id = dataset_id
    e.variables = variables
    
    # checking data frame is not empty
    try:
        df = e.to_pandas(
            index_col='time (UTC)',
            parse_dates=True,
            skiprows=(1,)  # units information can be dropped.
        ).dropna().tz_localize(None)
    except rHTTPError:
        print("Please enter a valid dataset id")
    return df


def active_drifters(bbox=None, time_start=None, time_end=None):
    bbox = bbox or [-100, -40, 18, 60]
    time_end = time_end or dt.date.today()
    time_start = time_start or (time_end - dt.timedelta(days=1))
    t0 = time_start.strftime('%Y-%m-%dT%H:%M:%SZ')
    t1 = time_end.strftime('%Y-%m-%dT%H:%M:%SZ')

    e = ERDDAP(server='OSMC', protocol="tabledap")
    e.dataset_id = "gdp_interpolated_drifter"

    # Setting constraints
    e.constraints = {
        "time>=": t0,
        "time<=": t1,
        'longitude>=': bbox[0],
        'longitude<=': bbox[1],
        'latitude>=': bbox[2],
        'latitude<=': bbox[3],
    }

    # e.variables = [
    #     "WMO",
    #     "latitude",
    #     "longitude",
    #     "time",
    # ]

    try:
        df = e.to_pandas()
    except ValueError:
        return pd.DataFrame()

    return df


def get_ndbc(bbox=None, time_start=None, time_end=None, buoy=None):
    bbox = bbox or [-100, -45, 5, 46]
    time_end = time_end or dt.date.today()
    time_start = time_start or (time_end - dt.timedelta(days=1))
    buoy = buoy or False
    time_formatter = '%Y-%m-%dT%H:%M:%SZ'

    e = ERDDAP(
        server='CSWC',
        protocol='tabledap',
        response='csv'
    )

    e.dataset_id = 'cwwcNDBCMet'
    e.constraints = {
        'time>=': time_start.strftime(time_formatter),
        'time<=': time_end.strftime(time_formatter),
    }

    if bbox:
        e.constraints['longitude>='] = bbox[0]
        e.constraints['longitude<='] = bbox[1]
        e.constraints['latitude>='] = bbox[2]
        e.constraints['latitude<='] = bbox[3]

    e.variables = [
        "station",
        "latitude",
        "longitude",
        "time"
    ]

    if buoy:
        e.constraints['station='] = buoy

    df = e.to_pandas(
        parse_dates=['time (UTC)'],
        skiprows=(1,)  # units information can be dropped.
    ).dropna()

    stations = df.station.unique()

    # e.variables = [
    #     "station",
    #     "latitude",
    #     "longitude",
    #     "wd",
    #     "wspd",
    #     "gst",
    #     "wvht",
    #     "dpd",
    #     "apd",
    #     "mwd",
    #     "bar",
    #     "atmp",
    #     "wtmp",
    #     "dewp",
    #     # "vis",
    #     # "ptdy",
    #     # "tide",
    #     "wspu",
    #     "wspv",
    #     "time",
    # ]

    try:
        df = e.to_pandas(
            parse_dates=['time (UTC)'],
            skiprows=(1,)  # units information can be dropped.
        ).dropna()
    except rHTTPError:
        df = pd.DataFrame()

    return df


def get_bathymetry(bbox=None):
    """
    Function to select bathymetry within a bounding box.
    This function pulls GEBCO 2014 bathy data from hfr.marine.rutgers.edu 

    Args:
        bbox (list, optional): Cartopy bounding box. Defaults to None.

    Returns:
        xarray.Dataset: xarray Dataset containing bathymetry data
    """
    bbox = bbox or [-100, -45, 5, 46]

    lons = bbox[:2]
    lats = bbox[2:]

    e = ERDDAP(
        server="https://hfr.marine.rutgers.edu/erddap/",
        protocol="griddap"
    )

    e.dataset_id = "bathymetry_gebco_2014_grid"

    e.griddap_initialize()

    # Modify constraints
    e.constraints["latitude<="] = max(lats)
    e.constraints["latitude>="] = min(lats)
    e.constraints["longitude>="] = max(lons)
    e.constraints["longitude<="] = min(lons)

    # return xarray dataset
    return e.to_xarray()


def get_ohc(bbox=None, time=None):
    bbox = bbox or [-100, -45, 5, 46]

    lons = bbox[:2]
    lats = bbox[2:]

    time = time or dt.date.today()

    e = ERDDAP(
        server="https://coastwatch.noaa.gov/erddap/",
        protocol="griddap"
    )

    e.dataset_id = "noaacwOHCna"

    e.griddap_initialize()

    # Modify constraints
    e.constraints["latitude<="] = max(lats)
    e.constraints["latitude>="] = min(lats)
    e.constraints["longitude>="] = max(lons)
    e.constraints["longitude<="] = min(lons)
    e.constraints['time>='] = time.strftime(time_formatter)
    e.constraints['time<='] = time.strftime(time_formatter)
    
    # e.griddap_initialize()

    # return xarray dataset
    try:
        return e.to_xarray()
    except requests.exceptions.HTTPError:
        print("No data available for this time period.")
        return