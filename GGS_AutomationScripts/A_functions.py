# =========================
# X - IMPORTS
# =========================

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime as dt
from datetime import datetime as datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import xarray as xr

# =========================
# CONFIG FUNCTIONS
# =========================

### FUNCTION:
def get_date_list(target_datetime):
    
    '''
    Get the list of dates for the next 24 hours in 6 hour intervals, formatted as 'YYYY-MM-DDT00:00:00Z'.

    Args:
    - target_datetime (dt.datetime): The target datetime.

    Returns:
    - date_list (list of str): The list of formatted dates for the next 24 hours in 6 hour intervals.
    '''

    date_start = target_datetime.replace(hour=0, minute=0, second=0, microsecond=0)
    date_end = target_datetime.replace(hour=0, minute=0, second=0, microsecond=0) + dt.timedelta(days=1)
    freq = '6H'

    date_range = pd.date_range(date_start, date_end, freq=freq)
    date_list = [date.strftime('%Y-%m-%dT%H:%M:%SZ') for date in date_range]

    return date_list

# =========================
# CHECK FUNCTIONS
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
# CALCULATE FUNCTIONS
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

    print("\n")
    print(f"Input Coordinates: ({target_lat}, {target_lon})")
    print(f"Dataset Indices: ({y_index}, {x_index})")
    print(f"Dataset Coordinates: ({lat_index:.3f}, {lon_index:.3f})")

    return (y_index, x_index), (lat_index, lon_index)

### FUNCTION:
def calculate_ticks(extent, direction):
    
    '''
    Define major and minor tick locations and major tick labels

    Args:
    - extent (tuple or list): extent (x0, x1, y0, y1) of the map in the given coordinate system.
    - direction (str): 'longitude' or 'latitude'

    Returns:
    - minor_ticks (numpy.ndarray): minor tick locations
    - major_ticks (numpy.ndarray): major tick locations
    - major_tick_labels (list): major tick labels
    '''

    direction = direction.lower()

    extent = [float(x) for x in extent]

    if direction == 'longitude':
        l0 = extent[0]
        l1 = extent[1]
        o0 = extent[2]
        o1 = extent[3]
    elif direction == 'latitude':
        l0 = extent[2]
        l1 = extent[3]
        o0 = extent[0]
        o1 = extent[1]

    r = np.max([l1 - l0, o1 - o0])

    if r <= 1.5:
        minor_int = 1.0 / 12.0
        major_int = 1.0 / 4.0
    elif r <= 3.0:
        minor_int = 1.0 / 6.0
        major_int = 0.5
    elif r <= 7.0:
        minor_int = 0.25
        major_int = float(1)
    elif r <= 15:
        minor_int = 0.5
        major_int = float(2)
    elif r <= 30:
        minor_int = float(1)
        major_int = float(3)
    elif r <=50:
        minor_int = float(1)
        major_int = float(5)
    elif r <=80:
        minor_int = float(5)
        major_int = float(10)
    elif r <=120:
        minor_int = float(5)
        major_int = float(15)
    elif r <=160:
        minor_int = float(5)
        major_int = float(20)
    elif r <=250:
        minor_int = float(10)
        major_int = float(30)
    else:
        minor_int = float(15)
        major_int = float(45)

    minor_ticks = np.arange(
        np.ceil(l0 / minor_int) * minor_int, 
        np.ceil(l1 / minor_int) * minor_int + minor_int,
        minor_int)
    minor_ticks = minor_ticks[minor_ticks <= l1]
    
    major_ticks = np.arange(
        np.ceil(l0 / major_int) * major_int, 
        np.ceil(l1 / major_int) * major_int + major_int,
        major_int)
    major_ticks = major_ticks[major_ticks <= l1]

    if major_int < 1:
        degree, minute, second = DD_to_DMS(np.array(major_ticks))
        if direction == 'longitude':
            n = 'W' * sum(degree < 0)
            p = 'E' * sum(degree >= 0)
            dir = n + p
            major_tick_labels = [str(np.abs(int(degree[i]))) + u"\N{DEGREE SIGN}" + str(int(minute[i])) + "'" + dir[i] for i in range(len(degree))]
        elif direction == 'latitude':
            n = 'S' * sum(degree < 0)
            p = 'N' * sum(degree >= 0)
            dir = n + p
            major_tick_labels = [str(np.abs(int(degree[i]))) + u"\N{DEGREE SIGN}" + str(int(minute[i])) + "'" + dir[i] for i in range(len(degree))]
        else:
            major_tick_labels = [str(int(degree[i])) + u"\N{DEGREE SIGN}" + str(int(minute[i])) + "'" for i in range(len(degree))]
    else:
        degree = major_ticks
        if direction == 'longitude':
            n = 'W' * sum(degree < 0)
            p = 'E' * sum(degree >= 0)
            dir = n + p
            major_tick_labels = [str(np.abs(int(degree[i]))) + u"\N{DEGREE SIGN}" + dir[i] for i in range(len(degree))]
        elif direction == 'latitude':
            n = 'S' * sum(degree < 0)
            p = 'N' * sum(degree >= 0)
            dir = n + p
            major_tick_labels = [str(np.abs(int(degree[i]))) + u"\N{DEGREE SIGN}" + dir[i] for i in range(len(degree))]
        else:
            major_tick_labels = [str(int(degree[i])) + u"\N{DEGREE SIGN}" for i in range(len(degree))]

    return minor_ticks, major_ticks, major_tick_labels

# =========================
# CONVERSION FUNCTIONS
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
def DD_to_DMS(DD):
    
    '''
    Convert decimal degrees to degrees, minutes, and seconds.

    Args:
    - decimal_degrees (np.ndarray): Numpy array of decimal degrees.

    Returns:
    - degrees (np.ndarray): Degrees part of the DMS.
    - minutes (np.ndarray): Minutes part of the DMS.
    - seconds (np.ndarray): Seconds part of the DMS.
    '''

    negative_DD = DD < 0

    absolute_decimal_degrees = np.abs(DD)

    degree = np.floor(absolute_decimal_degrees)

    remainder = (absolute_decimal_degrees - degree) * 60
    minute = np.floor(remainder)
    second = np.round((remainder - minute) * 60)

    degree[negative_DD] *= -1

    return degree, minute, second

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

# =========================
# PLOT FUNCTIONS
# =========================

### FUNCTION:
def plot_formatted_ticks(ax, extent_lon, extent_lat, proj=ccrs.PlateCarree(), fontsize=13, label_left=True, label_right=False, label_bottom=True, label_top=False, gridlines=True):
    
    '''
    Calculate and add formatted tick marks to the map based on longitude and latitude extents.

    Args:
    - ax (matplotlib.axes._subplots.AxesSubplot): The axes to set the ticks for.
    - extent_lon (list): Longitude bounds of the map, [min_longitude, max_longitude].
    - extent_lat (list): Latitude bounds of the map, [min_latitude, max_latitude].
    - proj (cartopy.crs class, optional): Define a projected coordinate system for ticks. Defaults to ccrs.PlateCarree().
    - fontsize (int, optional): Font size of tick labels. Defaults to 13.
    - gridlines (bool, optional): Add gridlines to map. Defaults to False.
    
    Returns:
    - None
    '''
    
    if not (len(extent_lon) == 2 and len(extent_lat) == 2):
        raise ValueError("extent_lon and extent_lat must each contain exactly two elements: [min_val, max_val].")
    overall_extent = [extent_lon[0], extent_lon[1], extent_lat[0], extent_lat[1]]
    
    minor_lon_ticks, major_lon_ticks, major_lon_labels = calculate_ticks(overall_extent, 'longitude')
    ax.set_xticks(minor_lon_ticks, minor=True, crs=proj)
    ax.set_xticks(major_lon_ticks, crs=proj)
    ax.set_xticklabels(major_lon_labels, fontsize=fontsize)

    minor_lat_ticks, major_lat_ticks, major_lat_labels = calculate_ticks(overall_extent, 'latitude')
    ax.set_yticks(minor_lat_ticks, minor=True, crs=proj)
    ax.set_yticks(major_lat_ticks, crs=proj)
    ax.set_yticklabels(major_lat_labels, fontsize=fontsize)

    ax.tick_params(which='major', direction='out', bottom=True, top=True, labelbottom=label_bottom, labeltop=label_top, left=True, right=True, labelleft=label_left, labelright=label_right, length=5, width=2)

    ax.tick_params(which='minor', direction='out', bottom=True, top=True, left=True, right=True, width=1)

    if gridlines:
        gl = ax.gridlines(draw_labels=False, linewidth=.25, color='black', alpha=0.1, linestyle='--', crs=proj, zorder=1000) # zorder = [1000]
        gl.xlocator = mticker.FixedLocator(minor_lon_ticks)
        gl.ylocator = mticker.FixedLocator(minor_lat_ticks)

### FUNCTION:
def plot_contour_cbar(magnitude, max_levels=10):
    
    '''
    Dynamically calculate levels and ticks for a matplotlib colorbar based on the magnitude.
    Merge the uppermost level with the second-highest level if it contains a very small percentage of data.

    Args:
    - magnitude (array-like): The magnitude data.
    - max_levels (int): Maximum number of levels. Default is 10.

    Returns:
    - levels (array-like): The level values for contour plot.
    - ticks (array-like): The tick values for the colorbar, aligned with levels.
    '''

    valid_magnitude = magnitude[~np.isnan(magnitude)]
    if valid_magnitude.size == 0:
        return [], []

    min_val = np.min(valid_magnitude)
    max_val = np.max(valid_magnitude)
    magnitude_range = max_val - min_val

    interval_options = [0.1, 0.2, 0.5, 1.0]
    best_interval = min(interval_options, key=lambda x: abs(max_levels - np.ceil(magnitude_range / x)))

    upper_threshold = np.ceil(max_val / best_interval) * best_interval
    while upper_threshold - min_val > magnitude_range:
        upper_threshold -= best_interval

    levels = np.arange(min_val, upper_threshold + best_interval, best_interval)
    levels = levels[levels <= max_val + best_interval * 0.1]

    if len(levels) > 2:
        upper_bound = levels[-1]
        lower_bound = levels[-2]
        upper_interval_count = np.sum((valid_magnitude > lower_bound) & (valid_magnitude <= upper_bound))
        total_count = len(valid_magnitude)
        upper_interval_percent = upper_interval_count / total_count

        if upper_interval_percent <= 0.01:
            magnitude[magnitude > lower_bound] = lower_bound
            levels = levels[:-1]

    ticks = levels

    return levels, ticks

### FUNCTION
def plot_bathymetry(ax, model_data, isobath1=-100, isobath2=-1000, show_legend=False):
    
    '''
    Add bathymetry to a plot.

    Args:
    - ax (matplotlib.axes._subplots.AxesSubplot): Matplotlib subplot.
    - model_data (xarray.core.dataset.Dataset): Model data.
    - isobath1 (int): First isobath level.
        - default: -100
    - isobath2 (int): Second isobath level.
        - default: -1000
    - show_legend (bool): Show legend.
        - default: False

    Returns:
    - none
    '''

    bathymetry_file = 'C:/Users/sal_f/OneDrive/Desktop/STF-0/!-GGS/GGS_Files/GEBCO_2023_sub_ice_topo.nc'
    bathy_data = xr.open_dataset(bathymetry_file)

    subset_bathy = bathy_data.sel(lat=slice(model_data.lat.min(), model_data.lat.max()), lon=slice(model_data.lon.min(), model_data.lon.max()))

    isobath_levels = sorted([isobath1, isobath2])
    depth_intervals = [-np.inf] + isobath_levels + [0]

    subset_bathy['elevation'] = subset_bathy['elevation'].fillna(0).where(~np.isinf(subset_bathy['elevation']), 0)

    cornflowerblue = mcolors.to_rgba('cornflowerblue')
    water = cfeature.COLORS['water']
    lightsteelblue = mcolors.to_rgba('lightsteelblue')

    colors = [cornflowerblue, water, lightsteelblue]
    
    ax.contour(subset_bathy.lon, subset_bathy.lat, subset_bathy.elevation, levels=isobath_levels, colors='dimgrey', linestyles='dashed', linewidths=0.25, zorder=50)

    for i in range(len(depth_intervals) - 1):
        ax.contourf(subset_bathy.lon, subset_bathy.lat, subset_bathy.elevation,
                    levels=[depth_intervals[i], depth_intervals[i + 1]], colors=[colors[i]])

    if show_legend:
        isobath1 = -isobath1
        isobath2 = -isobath2
        legend_colors = [lightsteelblue, water, cornflowerblue]
        legend_labels = [f'0m - {isobath1}m', f'{isobath1}m - {isobath2}m', f'> {isobath2}m']
        patches = [plt.plot([], [], marker="o", ms=10, ls="", mec=None, color=color, label=label)[0] for color, label in zip(legend_colors, legend_labels)]
        legend = ax.legend(handles=patches, loc='upper left', facecolor='lightgrey', edgecolor='black', framealpha=0.75, fontsize='x-small', markerscale=0.75)
        legend.set_zorder(1000)
        for text in legend.get_texts():
            text.set_color('black')

# =========================
# DATETIME FUNCTIONS
# =========================

### FUNCTION:
def datetime_format(input_datetime):
    
    '''
    Format a datetime string from 'YYYY-MM-DDTHH:MM:SSZ' to 'YYYYMMDDTHH'.

    Args:
    - input_datetime (str): The datetime string to format.

    Returns:
    - str: The formatted datetime string.
    '''
    
    if input_datetime != 'unknown_datetime':
        datetime_obj = datetime.strptime(input_datetime, "%Y-%m-%dT%H:%M:%SZ")
        formatted_datetime = datetime_obj.strftime("%Y%m%dT%H")
    else:
        formatted_datetime = 'unknown_datetime'
    
    return formatted_datetime

### FUNCTION:
def datetime_filename(model_data):

    '''
    Format a datetime string from 'YYYY-MM-DDTHH:MM:SSZ' to 'YYYYMMDDTHH'.

    Args:
    - model_data (xarray.core.dataset.Dataset): Model data.

    Returns:
    - str: The datetime formatted as 'YYYYMMDDTHH'.
    '''

    model_datetime = model_data.attrs.get('requested_datetime', 'unknown_datetime')

    if model_datetime != 'unknown_datetime':
        try:
            datetime_obj = datetime.strptime(model_datetime, "%Y-%m-%dT%H:%M:%SZ")
            filename_datetime = datetime_obj.strftime("%Y%m%dT%H")
        except ValueError:
            filename_datetime = 'invalid_datetime'
    else:
        filename_datetime = 'unknown_datetime'
    
    return filename_datetime

### FUNCTION:
def datetime_title(model_data):
    
    '''
    Format a datetime string from 'YYYY-MM-DDTHH:MM:SSZ' to 'YYYY-MM-DD HH:MM'.

    Args:
    - model_data (xarray.core.dataset.Dataset): Model data.

    Returns:
    - str: The datetime formatted as 'YYYY-MM-DD HH:MM'.
    '''
    
    model_datetime = model_data.attrs.get('requested_datetime', 'unknown_datetime')

    if model_datetime != 'unknown_datetime':
        try:
            datetime_obj = datetime.strptime(model_datetime, "%Y-%m-%dT%H:%M:%SZ")
            formatted_datetime = datetime_obj.strftime("%Y-%m-%d %H:%M")
        except ValueError:
            formatted_datetime = 'invalid_datetime'
    else:
        formatted_datetime = 'unknown_datetime'

    return formatted_datetime