# =========================
# X - IMPORTS
# =========================

import cartopy.crs as ccrs
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import numpy as np
import xarray as xr

# =========================
# CHECK INPUTS
# =========================

### FUNCTION:
EXIT_KEYWORD = "EXIT"
def check_abort(input_string):
    
    '''
    Check if the input string matches the exit keyword.
    If the input is the exit keyword, abort the program.
    
    Args:
    - input_string (str): The input string to be checked.

    Returns:
    - None
    '''
    
    if input_string.upper() == EXIT_KEYWORD:
        print("!!! [CONFIGURATION ABORTED] !!!")
        exit()

### FUNCTION:
def check_float(prompt_msg):
    
    '''
    Check for a valid float input.
    If invalid, prompt again with an error message.
    
    Args:
    - prompt_msg (str): The message to prompt the user with.

    Returns:
    - float: The parsed float value.
    '''
    
    while True:
        user_input = input(prompt_msg)
        check_abort(user_input)
        try:
            return float(user_input)
        except ValueError:
            prompt_msg = "[INPUT ERROR] " + prompt_msg

### FUNCTION:
def check_coordinate(coord_str):
    
    '''
    Check if the input decimal degree coordinates can be parsed into a valid latitude and longitude coordinate.
    
    Args:
    - coord_str (str): The input string containing coordinate values.

    Returns:
    - tuple: A tuple containing the latitude and longitude if valid. None otherwise.
    '''
    
    try:
        lat, lon = map(float, coord_str.split(","))
        if -90 <= lat <= 90 and -180 <= lon <= 180:
            return (lat, lon)
    except ValueError:
        return None

# =========================
# CALCULATE VARIABLES
# =========================

### FUNCTION:
def calculate_distance(coordinate_1, coordinate_2):
    
    '''
    Calculate the Haversine distance between two sets of decimal degree coordinates.
    
    Args:
    - coordinate_1 (tuple): Starting coordinate as (latitude, longitude).
    - coordinate_2 (tuple): Ending coordinate as (latitude, longitude).

    Returns:
    - float: The distance between the two coordinates in meters.
    '''
    
    R = 6371000

    lat1, lon1 = np.radians(coordinate_1)
    lat2, lon2 = np.radians(coordinate_2)
    
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    a = (np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    return R * c

### FUNCTION:
def calculate_heading(coordinate_1, coordinate_2):
    
    '''
    Calculate the compass bearing between two sets of decimal degree coordinates.
    
    Args:
    - coordinate_1 (tuple): Starting coordinate as (latitude, longitude).
    - coordinate_2 (tuple): Ending coordinate as (latitude, longitude).

    Returns:
    - float: The bearing from coordinate_1 to coordinate_2 in degrees.
    '''

    lat1, lon1 = np.radians(coordinate_1)
    lat2, lon2 = np.radians(coordinate_2)
    
    dlon = lon2 - lon1
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - (np.sin(lat1) * np.cos(lat2) * np.cos(dlon))
    
    initial_bearing = np.degrees(np.arctan2(x, y))

    return (initial_bearing + 360) % 360

### FUNCTION:
def calculate_gridpoint(model_data, target_lat, target_lon):
    
    '''
    Calculate the nearest XY gridpoint in a model dataset to the input latitude and longitude.

    Args:
    - model_data (xarray.Dataset): The model dataset.
    - target_lat (float): The target latitude.
    - target_lon (float): The target longitude.

    Returns:
    - (y_index, x_index) (tuple): The indices of the nearest point in the dataset.
    - (lat_index, lon_index) (tuple): The coordinates of the nearest point in the dataset.
    '''

    lat_diff = model_data['lat'] - float(target_lat)
    lon_diff = model_data['lon'] - float(target_lon)
    distance_square = lat_diff**2 + lon_diff**2

    y_index, x_index = np.unravel_index(distance_square.argmin(), distance_square.shape)

    lat_index = model_data['lat'].isel(y=y_index, x=x_index).values
    lon_index = model_data['lon'].isel(y=y_index, x=x_index).values

    print(f"Input Coordinates: ({target_lat}, {target_lon})")
    print(f"Dataset Indices: ({y_index}, {x_index})")
    print(f"Dataset Coordinates: ({lat_index:.3f}, {lon_index:.3f})")

    return (y_index, x_index), (lat_index, lon_index)

# =========================
# CONVERT VARIABLES
# =========================

### FUNCTION:
def DD_to_DM(DD, coord_type='longitude'):
    
    '''
    Convert a decimal degree (DD) coordinate to a degree minute (DM) coordinate.

    Args:
    - DD (float): A decimal degree coordinate.
    - coord_type (str): A string indicating whether the position is a longitude or latitude.

    Returns:
    - str: A degree minute coordinate.
    '''

    degrees = int(DD)
    minutes = abs(int((DD - degrees) * 60))

    direction = ''
    if degrees > 0:
        if coord_type == 'longitude':
            direction = 'E'
        elif coord_type == 'latitude':
            direction = 'N'
    elif degrees < 0:
        if coord_type == 'longitude':
            direction = 'W'
        elif coord_type == 'latitude':
            direction = 'S'

    return f"{abs(degrees)}°{minutes}'{direction}"

### FUNCTION:
def DD_to_DMS(DD, coord_type='longitude'):

    '''
    Convert a decimal degree (DD) coordinate to a degree minute second (DMS) coordinate.

    Args:
    - DD (float): A decimal degree coordinate.
    - coord_type (str): A string indicating whether the position is a longitude or latitude.

    Returns:
    - str: A degree minute second coordinate.
    '''

    degrees = int(DD)
    minutes = int((DD - degrees) * 60)
    seconds = (DD - degrees - minutes/60) * 3600.00
    
    direction = ''
    if degrees > 0:
        if coord_type == 'longitude':
            direction = 'E'
        elif coord_type == 'latitude':
            direction = 'N'
    elif degrees < 0:
        if coord_type == 'longitude':
            direction = 'W'
        elif coord_type == 'latitude':
            direction = 'S'
    
    return f"{abs(degrees)}°{abs(minutes)}'{abs(seconds):.2f}\"{direction}"

### FUNCTION:
def DD_to_DDM(DD):
    
    '''
    Convert a decimal degree (DD) coordinate to a degree decimal minute (DDM) coordinate.

    Args:
    - DD (float): A decimal degree coordinate.

    Returns:
    - str: A degree decimal minute coordinate.    
    '''
    
    degrees = int(DD)
    minutes = abs(DD - degrees) * 60

    return f"{degrees:02d}{minutes:05.2f}"

def get_6h_interval():
    
    '''
    Get the previous 6 hour interval in UTC time.

    Args:
    - None
    
    Returns:
    - str: The start time of the previous 6 hour interval in UTC time.
    '''

    current_time = datetime.utcnow()
    rounded_hour = current_time.hour - (current_time.hour % 6)
    previous_interval_time = current_time.replace(hour=rounded_hour, minute=0, second=0, microsecond=0)
    
    return previous_interval_time.strftime("%H:%M")

# =========================
# PLOT HELPERS
# =========================

### FUNCTION:

def set_map_ticks(ax, extent_lon, extent_lat):
    
    '''
    Set simple tick marks for a matplotlib map based on the extent.

    Args:
    - ax (matplotlib.axes): The axes to set the ticks for.
    - extent_lon (array-like): The longitude extent of the map.
    - extent_lat (array-like): The latitude extent of the map.

    Returns:
    - None
    '''

    def get_tick_interval(extent):
        
        '''
        Get the tick interval for a given extent.

        Args:
        - extent (array-like): The extent of the map.

        Returns:
        - float: The tick interval.
        '''
        
        min_extent, max_extent = np.min(extent), np.max(extent)
        range_extent = max_extent - min_extent

        intervals = [1, 2, 5, 10, 15, 30]
        for interval in intervals:
            num_ticks = int(range_extent / interval) + 1
            if 3 <= num_ticks <= 6:
                return interval
        return intervals[-1]

    lon_interval = get_tick_interval(extent_lon)
    lat_interval = get_tick_interval(extent_lat)

    lon_start = np.ceil(np.min(extent_lon) / lon_interval) * lon_interval
    lon_end = np.floor(np.max(extent_lon) / lon_interval) * lon_interval
    lat_start = np.ceil(np.min(extent_lat) / lat_interval) * lat_interval
    lat_end = np.floor(np.max(extent_lat) / lat_interval) * lat_interval

    lon_ticks = np.arange(lon_start, lon_end + lon_interval/2, lon_interval)
    lat_ticks = np.arange(lat_start, lat_end + lat_interval/2, lat_interval)

    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())

    def degree_formatter(val, pos):
        
        '''
        Format the tick labels to be in degree minute format.

        Args:
        - val (float): The tick value.
        - pos (int): The tick position.

        Returns:
        - str: The formatted tick label.
        '''
        
        return f"{int(val)}°"

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(degree_formatter))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(degree_formatter))

    # Minor ticks within the extent
    minor_lon_ticks = np.arange(np.floor(np.min(extent_lon)), np.ceil(np.max(extent_lon)) + 1, 1)
    minor_lat_ticks = np.arange(np.floor(np.min(extent_lat)), np.ceil(np.max(extent_lat)) + 1, 1)

    # Filter minor ticks to keep only those that are within the plot extent
    minor_lon_ticks = minor_lon_ticks[(minor_lon_ticks >= np.min(extent_lon)) & (minor_lon_ticks <= np.max(extent_lon))]
    minor_lat_ticks = minor_lat_ticks[(minor_lat_ticks >= np.min(extent_lat)) & (minor_lat_ticks <= np.max(extent_lat))]

    ax.set_xticks(minor_lon_ticks, minor=True, crs=ccrs.PlateCarree())
    ax.set_yticks(minor_lat_ticks, minor=True, crs=ccrs.PlateCarree())

    ax.tick_params(axis='both', which='both', direction='out', length=6, labelbottom=True, labelleft=True, labelright=False, labeltop=False, bottom=True, top=True, left=True, right=True)
    ax.tick_params(axis='both', which='minor', length=3)

### FUNCTION:
def calculate_cbar_ticks(magnitude):
    
    '''
    Calculate simple tick marks for a matplotlib colorbar based on the magnitude.

    Args:
    - magnitude (array-like): The magnitude data.

    Returns:
    - ticks (array-like): The tick values.
    '''
    
    valid_magnitude = magnitude[~np.isnan(magnitude)]
    if valid_magnitude.size == 0:
        return []

    min_val = np.min(valid_magnitude)
    max_val = np.max(valid_magnitude)

    interval_options = [0.1, 0.2]
    optimal_interval = min(interval_options, key=lambda x: len(np.arange(round(min_val / x) * x, round(max_val / x) * x + x, x)))

    ticks = np.arange(round(min_val / optimal_interval) * optimal_interval, round(max_val / optimal_interval) * optimal_interval + optimal_interval, optimal_interval)
    
    return ticks

### FUNCTION
def add_bathymetry(ax, model_data, isobath_levels=[-100, -1000], show_legend=True):
    
    '''
    Add isobath lines to a plot and fill ocean color between isobaths with varying shades of blue.

    Args:
    - ax: The axis object of the matplotlib plot.
    - model_data (xarray.Dataset): Dataset used to define the extent of the bathymetry data.
    - isobath_levels (list): Depths at which to plot the isobaths. Must be a negative value.
    - show_legend (bool): Whether to show the isobath legend.
        - default: True

    Returns:
    - None
    '''

    bathymetry_file = 'C:/Users/salfr/OneDrive/Desktop/STF-0/!-GGS/GGS_Files/GEBCO_2023_sub_ice_topo.nc'
    bathymetry = xr.open_dataset(bathymetry_file)

    data_lons = model_data.lon.values
    data_lats = model_data.lat.values
    extent = [np.min(data_lons), np.max(data_lons), np.min(data_lats), np.max(data_lats)]

    bathymetry_subset = bathymetry.sel(lon=slice(extent[0], extent[1]), lat=slice(extent[2], extent[3]))

    isobath_levels = sorted([-level for level in isobath_levels])
    min_elevation = -bathymetry_subset.elevation.max()
    max_elevation = -bathymetry_subset.elevation.min()

    isobath_levels = [lvl for lvl in isobath_levels if min_elevation <= lvl <= max_elevation]
    if not isobath_levels or isobath_levels[0] > min_elevation:
        isobath_levels.insert(0, min_elevation)
    if isobath_levels[-1] < max_elevation:
        isobath_levels.append(max_elevation)

    cmap = plt.get_cmap('Blues')
    cs = ax.contourf(bathymetry_subset.lon, bathymetry_subset.lat, -bathymetry_subset.elevation, levels=isobath_levels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())

    for level in isobath_levels:
        ax.contour(bathymetry_subset.lon, bathymetry_subset.lat, -bathymetry_subset.elevation, levels=[level], colors='dimgrey', linestyles='dashed', linewidths=0.25, zorder=50, transform=ccrs.PlateCarree())

    if show_legend:
        isobath_patches = []
        for i, level in enumerate(isobath_levels):
            color = cmap(float(i) / (len(isobath_levels) - 1))
            patch = mpatches.Patch(color=color, label=f'{-level} m')
            isobath_patches.append(patch)

        isobath_legend = ax.legend(handles=isobath_patches, loc='lower left', title='Isobath Levels', facecolor='grey', edgecolor='black', framealpha=0.75, fontsize='small', title_fontsize='medium')
        isobath_legend.set_zorder(100)
        isobath_legend.get_title().set_color('white')
        for text in isobath_legend.get_texts():
            text.set_color('white')

    bathymetry.close()
