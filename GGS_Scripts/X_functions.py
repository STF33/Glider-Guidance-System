# =========================
# IMPORTS
# =========================

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import dask.array
import datetime as dt
from datetime import datetime as datetime
from dateutil import parser
from erddapy import ERDDAP
from joblib import Parallel, delayed
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import multiprocessing
import numpy as np
import os
import pandas as pd
import xarray as xr

# =========================

# OPTIMIAL LOGIC FUNCTIONS

### FUNCTION:
def optimal_workers(power=1.0):
    
    '''
    Calculate the optimal number of workers for parallel processing based on the available CPU cores and a power factor.

    Args:
    - power (float): The percentage of available resources to use in processing.
        - default: 1.0
    
    Returns:
    - num_workers (int): The optimal number of workers for parallel processing.
    '''

    print(f"\n### ALLOCATING RESOURCES ###\n")

    if not 0 <= power <= 1:
        raise ValueError("Power must be between 0 and 1.")
    
    total_cores = os.cpu_count()
    
    if total_cores is None:
        total_cores = 4
    
    num_workers = max(1, math.floor(total_cores * power))
    print(f"Number of workers: {num_workers}")

    return num_workers

# DATA ACQUISITION FUNCTIONS

### FUNCTION:
def acquire_gliders(extent=None, target_date=dt.datetime.now(), date_delta=dt.timedelta(days=1), requested_variables=["time", "longitude", "latitude", "profile_id", "depth"], print_vars=False, target="all", request_timeout=5, enable_parallel=False):
    
    '''
    Fetches active glider datasets from the IOOS Glider DAC ERDDAP server, focusing on specified targets within a given spatial extent and time frame.
    
    Args:
    - extent (list): The spatial extent to search for gliders. Format: [min_lon, max_lon, min_lat, max_lat].
    - target_date (datetime.datetime): The target date for the search window.
        - default: dt.datetime.now()
    - date_delta (datetime.timedelta): The duration of the search window.
        - default: dt.timedelta(days=1)
    - requested_variables (list): The variables to request from the glider datasets.
        - default: ["time", "longitude", "latitude", "profile_id", "depth"]
    - print_vars (bool): Print available variables for each glider dataset.
        - default: False
    - target (str or list): The target glider(s) to fetch datasets for. Accepts "all" or specific glider IDs.
        - default: "all"
    - request_timeout (int): The timeout for the ERDDAP requests.
        - default: 5
    - enable_parallel (bool or int): Enables parallel downloads if True or an integer specifying the number of cores.
        - default: False
    
    Returns:
    - glider_dataframes (pandas.DataFrame): The concatenated glider datasets.
    '''

    print(f"\n### ACQUIRING GLIDER DATASETS ###\n")
    start_time = print_starttime()

    if extent is None:
        extent = [-180, 180, -90, 90]

    start_date = target_date - date_delta
    end_date = target_date

    formatted_start_date = start_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    formatted_end_date = end_date.strftime('%Y-%m-%dT%H:%M:%SZ')

    erddap_server = ERDDAP(server='https://data.ioos.us/gliders/erddap')
    erddap_server.requests_kwargs = {'timeout': request_timeout}

    search_params = {
        'min_time': formatted_start_date,
        'max_time': formatted_end_date,
        'min_lon': extent[0],
        'max_lon': extent[1],
        'min_lat': extent[2],
        'max_lat': extent[3],
    }

    search_url = erddap_server.get_search_url(search_for="gliders", response='csv', **search_params)

    try:
        search_results = pd.read_csv(search_url)
    except Exception as error:
        print(f"Error during initial glider search: {error}")
        return pd.DataFrame()

    glider_ids = search_results['Dataset ID'].values
    print(f"Found {len(glider_ids)} Glider Datasets within search window.")
    print("Glider Indexes found:", ", ".join(glider_ids))

    target_glider_ids = []
    if target != "all":
        target_glider_ids = [target] if isinstance(target, str) else target

    def available_variables(glider_id, print_vars):
        info_url = erddap_server.get_info_url(dataset_id=glider_id, response="csv")
        try:
            info = pd.read_csv(info_url)
            available_vars = info[info['Variable Name'].notnull()]['Variable Name'].tolist()
            if print_vars:
                print(f"Available variables for {glider_id}: {available_vars}")
            return available_vars
        except Exception as error:
            print(f"Error fetching available variables for glider {glider_id}: {error}")
            return []

    def acquire_glider_dataset(glider_id):
        if target_glider_ids and glider_id not in target_glider_ids:
            return glider_id, pd.DataFrame()

        available_vars = available_variables(glider_id, print_vars)
        vars_to_request = [var for var in requested_variables if var in available_vars]

        missing_vars = set(requested_variables) - set(vars_to_request)
        if missing_vars:
            print(f"Variables {missing_vars} not available for glider {glider_id} and will be skipped.")

        if not vars_to_request:
            print(f"No requested variables are available for glider {glider_id}.")
            return glider_id, pd.DataFrame()

        erddap_server.constraints = {}
        erddap_server.protocol = 'tabledap'
        erddap_server.variables = vars_to_request
        erddap_server.dataset_id = glider_id

        try:
            dataset_df = erddap_server.to_pandas(
                response="csv",
                index_col="time",
                parse_dates=True,
                skiprows=(1,)
            ).tz_localize(None)
            
            if 'lon' in dataset_df.columns and 'lat' in dataset_df.columns:
                dataset_df.rename(columns={'lon': 'longitude', 'lat': 'latitude'}, inplace=True)
            elif 'longitude' not in dataset_df.columns or 'latitude' not in dataset_df.columns:
                print(f"Missing 'longitude' or 'latitude' columns in dataset for glider {glider_id}. Requested variables might not be available.")
                return glider_id, pd.DataFrame()
            
            return glider_id, dataset_df
        except Exception as error:
            print(f"Error fetching dataset for glider {glider_id} with requested variables: {error}")
            return glider_id, pd.DataFrame()

    if enable_parallel:
        num_cores = multiprocessing.cpu_count() if isinstance(enable_parallel, bool) else enable_parallel
        parallel_downloads = Parallel(n_jobs=num_cores)(delayed(acquire_glider_dataset)(glider_id) for glider_id in glider_ids)
        glider_datasets = {glider: df for glider, df in parallel_downloads}
    else:
        glider_datasets = {glider: df for glider, df in (acquire_glider_dataset(glider_id) for glider_id in glider_ids)}
    
    non_empty_glider_datasets = {glider: df for glider, df in glider_datasets.items() if not df.empty}

    try:
        if non_empty_glider_datasets:
            glider_dataframes = pd.concat(non_empty_glider_datasets, names=["glider", "time"]).sort_index()
        else:
            print("No non-empty datasets found for concatenation.")
            glider_dataframes = pd.DataFrame()
    except ValueError as e:
        print(f"Error during DataFrame concatenation: {e}")
        glider_dataframes = pd.DataFrame()
        
    end_time = print_endtime()
    print_runtime(start_time, end_time)

    return glider_dataframes

# CALCULATE FUNCTIONS

### FUNCTION:
def calculate_gridpoint(model_data, target_lat, target_lon):
    
    '''
    Calculate the nearest XY gridpoint in a model dataset to the input latitude and longitude, accommodating both 1D and 2D lat/lon arrays.

    Args:
    - model_data (xarray.Dataset): The model dataset.
    - target_lat (float): The target latitude.
    - target_lon (float): The target longitude.

    Returns:
    - (y_index, x_index) (tuple): The indices of the nearest point in the dataset.
    - (lat_index, lon_index) (tuple): The coordinates of the nearest point in the dataset.
    '''
    
    if 'lat' in model_data.dims and 'lon' in model_data.dims:
        lat = model_data['lat'].values
        lon = model_data['lon'].values
        lon_grid, lat_grid = np.meshgrid(lon, lat)
        squared_dist = (lat_grid - target_lat)**2 + (lon_grid - target_lon)**2
    else:
        lat_diff = model_data['lat'] - target_lat
        lon_diff = model_data['lon'] - target_lon
        squared_dist = lat_diff**2 + lon_diff**2
    
    y_index, x_index = np.unravel_index(squared_dist.argmin(), squared_dist.shape)
    
    if 'lat' in model_data.dims and 'lon' in model_data.dims:
        lat_index = lat[y_index]
        lon_index = lon[x_index]
    else:
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

# PLOT FUNCTIONS

### FUNCTION:
def plot_formatted_ticks(ax, extent_lon, extent_lat, proj=ccrs.PlateCarree(), fontsize=10, label_left=True, label_right=False, label_bottom=True, label_top=False, gridlines=True):
    
    '''
    Calculate and add formatted tick marks to the map based on longitude and latitude extents.

    Args:
    - ax (matplotlib.axes._subplots.AxesSubplot): The axes to set the ticks for.
    - extent_lon (list): Longitude bounds of the map, [min_longitude, max_longitude].
    - extent_lat (list): Latitude bounds of the map, [min_latitude, max_latitude].
    - proj (cartopy.crs class, optional): Define a projected coordinate system for ticks.
        - default: ccrs.PlateCarree()
    - fontsize (int, optional): Font size of tick labels.
        - default: 10
    - label_left (bool, optional): Label the left side of the map.
        - default: True
    - label_right (bool, optional): Label the right side of the map.
        - default: False
    - label_bottom (bool, optional): Label the bottom side of the map.  
        - default: True
    - label_top (bool, optional): Label the top side of the map.
        - default: False
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
    - max_levels (int): Maximum number of levels.
        - default: 10
    - extend_max (bool): Extend the maximum color level to indicate values exceeding the set levels.
        - default: True

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
def plot_threshold_legend(ax, mag2, mag3, mag4, mag5):
    
    '''
    Creates and adds a custom legend to the plot indicating threshold zones for current magnitudes.

    Args:
    - ax (matplotlib.axes.Axes): The axes object to add the legend to.
    - mag2 (float): Second threshold magnitude.
    - mag3 (float): Third threshold magnitude.
    - mag4 (float): Fourth threshold magnitude.
    - mag5 (float): Fifth threshold magnitude.

    Returns:
    - None
    '''

    patches = [
        mpatches.Patch(facecolor='yellow', label=f'{mag2} - {mag3} m/s'),
        mpatches.Patch(facecolor='orange', label=f'{mag3} - {mag4} m/s'),
        mpatches.Patch(facecolor='orangered', label=f'{mag4} - {mag5} m/s'),
        mpatches.Patch(facecolor='maroon', label=f'> {mag5} m/s')
    ]

    threshold_legend = ax.legend(handles=patches, loc='upper right', facecolor='white', edgecolor='black', fontsize='x-small')
    threshold_legend.set_zorder(10000)
    for text in threshold_legend.get_texts():
        text.set_color('black')

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
        bathymetry_legend = ax.legend(handles=patches, loc='upper left', facecolor='white', edgecolor='black', fontsize='x-small', markerscale=0.75)
        bathymetry_legend.set_zorder(10000)
        for text in bathymetry_legend.get_texts():
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

    Returns:
    - None
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

### FUNCTION:
def plot_add_gliders(ax, glider_data_frame, legend=True):
    
    '''
    Adds glider tracks to a cartopy map.

    Args:
    - ax (cartopy.mpl.geoaxes.GeoAxesSubplot): The cartopy map to add the glider tracks to.
    - glider_data_frame (pd.DataFrame): The glider data as a pandas DataFrame.
    - legend (bool): Show legend.
        - default: True
        
    Returns:
    - None
    '''

    num_gliders = len(glider_data_frame.groupby(level=0))
    colors = plt.cm.jet(np.linspace(0, 1, num_gliders))
    
    legend_handles = []

    for (glider_id, glider_data), color in zip(glider_data_frame.groupby(level=0), colors):
        glider_name = glider_id.split('-2')[0]
        latest_observation = glider_data.iloc[-1]
        ax.plot(glider_data['longitude'], glider_data['latitude'], '-', color=color, linewidth=1.5, transform=ccrs.PlateCarree(), zorder=1000)
        ax.plot(latest_observation['longitude'], latest_observation['latitude'], '^', color=color, markeredgecolor='black', markersize=8.5, transform=ccrs.PlateCarree(), zorder=1000)

        custom_handle = mlines.Line2D([], [], color=color, markeredgecolor='black', markersize=8.5, label=glider_name, linestyle='-', marker='^')
        legend_handles.append(custom_handle)

    if legend:
        glider_legend = ax.legend(handles=legend_handles, title="Gliders", loc='lower right', facecolor='white', edgecolor='black', fontsize='x-small', markerscale=0.75)
        glider_legend.set_zorder(10000)
        for text in glider_legend.get_texts():
            text.set_color('black')

# FORMATTING FUNCTIONS

### FUNCTION:
def format_colorbar(ax, cbar):
    
    '''
    Adjusts the colorbar position to match the height of the map object (ax) it is plotted next to.

    Args:
    - ax (matplotlib.axes.Axes): The axes object containing the map.
    - cbar (matplotlib.colorbar.Colorbar): The colorbar object to adjust.

    Returns:
    - None
    '''
    
    ax_pos = ax.get_position()
    cbar_ax = cbar.ax
    cbar_pos = cbar_ax.get_position()
    
    cbar_height = ax_pos.height
    cbar_bottom = ax_pos.y0
    cbar_left = cbar_pos.x0
    cbar_width = cbar_pos.width

    cbar_ax.set_position([cbar_left, cbar_bottom, cbar_width, cbar_height])

### FUNCTION:
def format_figure_titles(ax, fig, config, datetime_index, model_name, title):
    
    '''
    Sets the main title, subtitle, and suptitle for a plot with correct positioning.

    Args:
    - ax (matplotlib.axes.Axes): The axes object to set the main title on.
    - fig (matplotlib.figure.Figure): The figure object to set subtitles and suptitles on.
    - config (dict): Configuration dictionary containing plot settings and metadata.
    - datetime_index (str): Datetime index used to derive the subtitle.
    - model_name (str): Model name used to derive the subtitle.
    - title (str): Custom title text for the main title of the plot.

    Returns:
    - None
    '''

    ax_pos = ax.get_position()
    top = ax_pos.y1
    bottom = ax_pos.y0
    
    title_distance_top = 0.1
    title_position = top + title_distance_top
    
    subtitle_distance_title = 0.05
    subtitle_position = title_position - subtitle_distance_title
    
    suptitle_distance_bottom = 0.1
    suptitle_position = bottom - suptitle_distance_bottom

    title_datetime = format_title_datetime(datetime_index)
    
    fig.text(0.5, title_position, title, fontsize=14, fontweight='bold', ha='center', va='bottom')

    subtitle_text = f"{model_name} {title_datetime} UTC"
    fig.text(0.5, subtitle_position, subtitle_text, fontsize=10, ha='center', va='bottom')

    suptitle_text = f"Generated by the Glider Guidance System (GGS) - {config['glider_name']}"
    fig.text(0.5, suptitle_position, suptitle_text, fontsize='smaller', fontweight='bold', ha='center', va='top', color='gray')

### FUNCTION:
def format_subplot_titles(fig, config, datetime, title):
    
    '''
    Sets the main title, subtitle, and suptitle for a figure.

    Args:
    - fig (matplotlib.figure.Figure): The figure object to set subtitles and suptitles on.
    - config (dict): Configuration dictionary containing plot settings and metadata.
    - datetime (str): Datetime index used to derive the subtitle.
    - title (str): Custom title text for the main title of the plot.

    Returns:
    - None
    '''

    fig.suptitle(title, fontsize=26, fontweight='bold', ha='center', y=0.95)
    
    title_datetime = format_title_datetime(datetime)
    subtitle_text = f"{title_datetime} UTC"
    fig.text(0.5, 0.92, subtitle_text, fontsize=22, ha='center', va='top')
    
    suptitle_text = f"Generated by the Glider Guidance System (GGS) - {config['glider_name']}"
    fig.text(0.5, 0.05, suptitle_text, fontsize=20, ha='center', va='top', color='gray')

### FUNCTION:
def format_subplot_headers(axes, fig, subplot_titles):

    '''
    Adds subplot headers as bolded titles above each subplot row figure.

    Args:
    - axes (numpy.ndarray): Array of axes objects to add headers to.
    - fig (matplotlib.figure.Figure): Figure object to add text to.
    - subplot_titles (list): List of strings to use as header titles.

    Returns:
    - None
    '''

    if axes.ndim == 1:
        axes = np.array([axes])

    for ax, model_name in zip(axes[:, 0], subplot_titles):
        bbox = ax.get_position()
        text_x = bbox.x0
        text_y = bbox.y1 + 0.02
        
        fig.text(text_x, text_y, model_name, fontsize=14, fontweight='bold', ha='left')

### FUNCTION:
def format_title_datetime(datetime):
    
    '''
    Format a datetime string from any recognized format to '%Y-%m-%d %H:%M'.

    Args:
    - datetime (str): The datetime string to format.
    
    Returns:
    - str: The datetime formatted as '%Y-%m-%d %H:%M'.
    '''
    

    if datetime:
        try:
            datetime_obj = parser.parse(datetime)
            formatted_datetime = datetime_obj.strftime('%Y-%m-%d %H:%M')
        except ValueError:
            formatted_datetime = 'unknown_datetime'
    
    return formatted_datetime

### FUNCTION:
def format_save_datetime(datetime):
    
    '''
    Format a datetime string from any recognized format to '%Y%m%dT%HZ'.

    Args:
    - datetime (str): The datetime string to format.

    Returns:
    - str: The formatted datetime string in the '%Y%m%dT%HZ' format.
    '''

    if datetime:
        try:
            datetime_obj = parser.parse(datetime)
            formatted_datetime = datetime_obj.strftime('%Y%m%dT%HZ')
        except ValueError:
            formatted_datetime = 'unknown_datetime'
        
    return formatted_datetime

# CONVERSION FUNCTIONS

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

# PRINT FUNCTIONS

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
    - start_time (datetime.datetime): The start time of the program.
    - end_time (datetime.datetime): The end time of the program.

    Returns:
    - None
    '''

    print(f"Run time: {end_time - start_time}")
