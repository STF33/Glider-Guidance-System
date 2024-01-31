# =========================
# X - IMPORTS
# =========================

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import dask.array
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
# CALCULATE FUNCTIONS
# =========================

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
def plot_contour_cbar(magnitude, max_levels=10, extend_max=True):
    
    '''
    Calculate levels for a matplotlib colorbar based on the magnitude. Optionally extend the maximum color level.

    Args:
    - magnitude (array-like or xarray.DataArray): The magnitude data.
    - max_levels (int): Maximum number of levels. Default is 10.
    - extend_max (bool): Extend the maximum color level to indicate values exceeding the set levels.
        - default: False

    Returns:
    - levels (array-like): The level values for contour plot.
    - ticks (array-like): The tick values for the colorbar, aligned with levels.
    - extend (str): String indicating if the colorbar should be extended. 'max', 'neither', or 'both'.
    '''

    if isinstance(magnitude, xr.DataArray):
        if isinstance(magnitude.data, dask.array.Array):
            magnitude = magnitude.compute()
        magnitude = magnitude.values

    valid_magnitude = magnitude[~np.isnan(magnitude)]
    if valid_magnitude.size == 0:
        return [], []

    min_val = np.min(valid_magnitude)
    max_val = np.max(valid_magnitude)
    magnitude_range = max_val - min_val

    interval_options = [0.1, 0.2, 0.5, 1.0]
    best_interval = min(interval_options, key=lambda x: abs(max_levels - np.ceil(magnitude_range / x)))

    levels = np.arange(min_val, max_val, best_interval)

    if extend_max and max_val > levels[-1]:
        extend = 'max'
    else:
        extend = 'neither'

    ticks = levels

    return levels, ticks, extend

### FUNCTION:
def plot_bathymetry(ax, config, model_data, isobath1=-100, isobath2=-1000, show_legend=False):
    
    '''
    Add bathymetry to a plot.

    Args:
    - ax (matplotlib.axes._subplots.AxesSubplot): Matplotlib subplot.
    - config (dict): Glider Guidance System mission configuration.
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

    bathymetry_path = config["bathymetry_path"]
    bathy_data = xr.open_dataset(bathymetry_path)

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
        legend = ax.legend(handles=patches, loc='upper left', facecolor='white', edgecolor='black', framealpha=0.75, fontsize='x-small', markerscale=0.75)
        legend.set_zorder(1000)
        for text in legend.get_texts():
            text.set_color('black')

### FUNCTION:
def plot_profile_thresholds(ax, data, threshold, color):
    
    '''
    Apply shading to regions where data exceeds the threshold and update the legend.
    
    Args:
    - ax (matplotlib.axes.Axes): The axes object to apply shading to.
    - data (array): Data array for the plot.
    - threshold (float): Threshold value for shading.
    - color (str): Color for the shaded region.
    '''

    depth_values = np.arange(len(data))
    regions = []
    start = None

    for i, value in enumerate(data):
        if value > threshold:
            if start is None:
                start = i
        else:
            if start is not None:
                regions.append((start, i))
                start = None

    if start is not None:
        regions.append((start, len(data)))

    for start, end in regions:
        ax.fill_betweenx(depth_values[start:end], ax.get_xlim()[0], ax.get_xlim()[1], color=color, alpha=0.25)

    ax.plot([], [], color=color, alpha=0.5, linewidth=10, label=f'Above Threshold = [{threshold}]')

# =========================
# CONVERSION FUNCTIONS
# =========================

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

# =========================
# PRINT FUNCTIONS
# =========================

### FUNCTION:
def print_starttime():

    '''
    Print the start time of the program.

    Args:
    - None

    Returns:
    - start_time (datetime.datetime): The start time of the program.
    '''

    start_time = dt.datetime.utcnow()
    print(f"Start time (UTC): {start_time}")

    return start_time

### FUNCTION:
def print_endtime():

    '''
    Print the end time of the program.

    Args:
    - None

    Returns:
    - end_time (datetime.datetime): The end time of the program.
    '''

    end_time = dt.datetime.utcnow()
    print(f"End time (UTC): {end_time}")

    return end_time

### FUNCTION:
def print_runtime(start_time, end_time):

    '''
    Print the runtime of the program.

    Args:
    - None

    Returns:
    - None
    '''

    print(f"Run time: {end_time - start_time}")
    