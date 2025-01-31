"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature
import cmocean.cm as cmo
import csv
import dask.array
import datetime as dt
from datetime import datetime as datetime
from dateutil import parser
from erddapy import ERDDAP
import glob
import heapq
from joblib import Parallel, delayed
import math
from math import radians, cos, sin, asin, sqrt
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

# OPERATIONAL FUNCTIONS

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

# ALGORITHM FUNCTIONS

### FUNCTION:
def compute_optimal_path(config, directory, model_dataset, glider_raw_speed=0.5):
    
    '''
    Calculates the optimal path between waypoints for a mission, considering the impact of ocean currents and distance.
    
    Args:
    - config (dict): A dictionary containing mission config details including waypoints.
    - directory (str): The directory path to save the output statistics file.
    - model_dataset (xarray.Dataset): An xarray dataset containing depth-averaged ocean current data.
    - glider_raw_speed (float, optional): The glider's base speed in meters per second.
        - default: 0.5

    Returns:
    - optimal_mission_path (list of tuples): A list of latitude and longitude tuples representing the optimal route.
    '''
    
    model_name = model_dataset.attrs['model_name']
    print(f"\n### COMPUTING OPTIMAL PATH [{model_name}] ###\n")
    start_time = print_starttime()

    model_name = model_dataset.attrs['model_name']
    csv_data = [("Segment Start", "Segment End", "Segment Time (s)", "Segment Distance (m)")]

    def calculate_haversine_distance(longitude1, latitude1, longitude2, latitude2):
        '''Calculates the great circle distance between two points on the earth using the Haversine formula.'''
        longitude1, latitude1, longitude2, latitude2 = map(radians, [longitude1, latitude1, longitude2, latitude2])
        delta_longitude = longitude2 - longitude1
        delta_latitude = latitude2 - latitude1
        a = sin(delta_latitude / 2)**2 + cos(latitude1) * cos(latitude2) * sin(delta_longitude / 2)**2
        distance = 2 * asin(sqrt(a)) * 6371000
        return distance

    def calculate_direct_path(start_index, end_index, glider_raw_speed):
        '''Fallback to the direct great circle path if no optimal path is found.'''
        start_lat, start_lon = convert_grid2coord(*start_index)
        end_lat, end_lon = convert_grid2coord(*end_index)
        distance = calculate_haversine_distance(start_lon, start_lat, end_lon, end_lat)
        time = distance / glider_raw_speed
        return [(start_lat, start_lon), (end_lat, end_lon)], time, distance

    def convert_coord2grid(latitude, longitude):
        '''Converts geographical latitude and longitude to the nearest index on the dataset grid.'''
        latitude_index = np.argmin(np.abs(latitude_array - latitude))
        longitude_index = np.argmin(np.abs(longitude_array - longitude))
        return latitude_index, longitude_index
    
    def convert_grid2coord(latitude_index, longitude_index):
        '''Converts dataset grid indices back to geographical latitude and longitude coordinates.'''
        latitude = latitude_array[latitude_index]
        longitude = longitude_array[longitude_index]
        return latitude, longitude
    
    def calculate_remaining_distance(current_index, goal_index):
        '''Estimates the distance from the current index to the goal using the Haversine formula as a heuristic.'''
        current_latitude, current_longitude = convert_grid2coord(*current_index)
        goal_latitude, goal_longitude = convert_grid2coord(*goal_index)
        estimated_distance = calculate_haversine_distance(current_longitude, current_latitude, goal_longitude, goal_latitude)
        return estimated_distance
    
    def calculate_movement_variables(model_dataset, start_index, end_index, glider_raw_speed):
        '''Calculates the time and distance cost of moving from one grid point to the next, considering ocean currents.'''
        start_lat, start_lon = convert_grid2coord(*start_index)
        end_lat, end_lon = convert_grid2coord(*end_index)
        if start_lat == end_lat and start_lon == end_lon:
            return 0, 0, 0
        heading_vector = np.array([end_lon - start_lon, end_lat - start_lat])
        norm = np.linalg.norm(heading_vector)
        if norm == 0:
            return 0, 0, 0
        heading_vector = heading_vector / norm
        u_current = model_dataset['u_depth_avg'].isel(lat=start_index[0], lon=start_index[1]).values.item()
        v_current = model_dataset['v_depth_avg'].isel(lat=start_index[0], lon=start_index[1]).values.item()
        current_vector = np.array([u_current, v_current])
        current_along_heading = np.dot(current_vector, heading_vector)
        net_speed = glider_raw_speed + current_along_heading
        net_speed = max(net_speed, 0.1)
        distance = calculate_haversine_distance(start_lon, start_lat, end_lon, end_lat)
        time = distance / net_speed
        return distance, time, current_along_heading
    
    def generate_adjacent_nodes(index):
        '''Generates neighboring index nodes for exploration based on the current index's position.'''
        latitude_index, longitude_index = index
        for delta_latitude in [-1, 0, 1]:
            for delta_longitude in [-1, 0, 1]:
                if delta_latitude == 0 and delta_longitude == 0:
                    continue
                new_latitude_index, new_longitude_index = latitude_index + delta_latitude, longitude_index + delta_longitude
                if 0 <= new_latitude_index < len(latitude_array) and 0 <= new_longitude_index < len(longitude_array):
                    yield (new_latitude_index, new_longitude_index)
    
    def reconstruct_path(came_from_dictionary, start_index, goal_index):
        '''Reconstructs the path from the start index to the goal index using the came_from dictionary populated by the A* algorithm.'''
        current_index = goal_index
        optimal_path = [current_index]
        while current_index != start_index:
            current_index = came_from_dictionary[current_index]
            optimal_path.append(current_index)
        optimal_path.reverse()
        optimal_path_coords = [convert_grid2coord(*index) for index in optimal_path]
        return optimal_path_coords
    
    def GGS_algorithm(model_dataset, start_index, end_index, glider_raw_speed):
        '''Executes the GGS: Current Mapper A* algorithm to find the most efficient path from the start index to the goal index.'''
        open_set = [(calculate_remaining_distance(start_index, end_index), start_index)]
        came_from = {start_index: None}
        distance_traveled = {start_index: 0}
        path_found = False
        while open_set:
            _, current_node = heapq.heappop(open_set)
            if current_node == end_index:
                path_found = True
                break
            current_distance_to_target = calculate_remaining_distance(current_node, end_index)
            for neighbor in generate_adjacent_nodes(current_node):
                neighbor_distance, neighbor_time, _ = calculate_movement_variables(model_dataset, current_node, neighbor, glider_raw_speed)
                neighbor_distance_to_target = calculate_remaining_distance(neighbor, end_index)
                distance_gain = current_distance_to_target - neighbor_distance_to_target
                distance_gain_per_time = distance_gain / neighbor_time
                tentative_distance_traveled = distance_traveled[current_node] + neighbor_distance / (glider_raw_speed + max(0.1, distance_gain_per_time))
                if tentative_distance_traveled < distance_traveled.get(neighbor, float('inf')):
                    came_from[neighbor] = current_node
                    distance_traveled[neighbor] = tentative_distance_traveled
                    heapq.heappush(open_set, (tentative_distance_traveled, neighbor))
        if path_found:
            path = reconstruct_path(came_from, start_index, end_index)
            time, distance = calculate_movement_variables(model_dataset, start_index, end_index, glider_raw_speed)[:2]
        else:
            print(f"Direct path used from {convert_grid2coord(*start_index)} to {convert_grid2coord(*end_index)}.")
            path, time, distance = calculate_direct_path(start_index, end_index, glider_raw_speed)
        return path, time, distance
    
    mission_waypoints = config['MISSION']['GPS_coords']
    mission_waypoints = [(float(lat), float(lon)) for lat, lon in mission_waypoints]
    latitude_array = model_dataset['lat'].values
    longitude_array = model_dataset['lon'].values
    
    optimal_mission_path = []
    total_time = 0
    total_distance = 0
    
    for i in range(len(mission_waypoints) - 1):
        start_index = convert_coord2grid(*mission_waypoints[i])
        end_index = convert_coord2grid(*mission_waypoints[i + 1])
        segment_path, segment_time, segment_distance = GGS_algorithm(model_dataset, start_index, end_index, glider_raw_speed)
        optimal_mission_path.extend(segment_path[:-1])
        total_time += segment_time
        total_distance += segment_distance
        
        csv_data.append((mission_waypoints[i], mission_waypoints[i+1], segment_time, segment_distance))
        print(f"Segment {i+1}: Start {mission_waypoints[i]} End {mission_waypoints[i+1]} Time {segment_time} seconds Distance {segment_distance} meters")
    
    optimal_mission_path.append(mission_waypoints[-1])
    
    print(f"Total mission time (adjusted): {total_time} seconds")
    print(f"Total mission distance: {total_distance} meters")

    csv_file_path = os.path.join(directory, f"{model_name}_mission_statistics.csv")
    with open(csv_file_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(csv_data)

    end_time = print_endtime()
    print_runtime(start_time, end_time)

    return optimal_mission_path

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
    print(f"Found {len(glider_ids)} Glider Datasets within search window: {formatted_start_date} to {formatted_end_date}")
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
def calculate_bearing(lat1, lon1, lat2, lon2):
    
    '''
    Calculate the compass bearing between two latitude and longitude coordinates.

    Args:
    - lat1 (float): Latitude of the first point.
    - lon1 (float): Longitude of the first point.
    - lat2 (float): Latitude of the second point.
    - lon2 (float): Longitude of the second point.

    Returns:
    - compass_bearing (float): The compass bearing between the two points.
    '''

    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])

    delta_lon = lon2 - lon1
    x = np.sin(delta_lon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - (np.sin(lat1) * np.cos(lat2) * np.cos(delta_lon))
    initial_bearing = np.arctan2(x, y)

    initial_bearing = np.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing

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
        degree, minute, second = convert_DD_to_DMS(np.array(major_ticks))
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
def plot_formatted_ticks(ax, extent_lon, extent_lat, proj=ccrs.Mercator(), fontsize=16, label_left=True, label_right=False, label_bottom=True, label_top=False, gridlines=True):
    
    '''
    Calculate and add formatted tick marks to the map based on longitude and latitude extents.

    Args:
    - ax (matplotlib.axes._subplots.AxesSubplot): The axes to set the ticks for.
    - extent_lon (list): Longitude bounds of the map, [min_longitude, max_longitude].
    - extent_lat (list): Latitude bounds of the map, [min_latitude, max_latitude].
    - proj (cartopy.crs class, optional): Define a projected coordinate system for ticks.
        - default: ccrs.Mercator()
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
    ax.set_xticklabels(major_lon_labels, fontsize=fontsize, rotation=45, ha='right')

    minor_lat_ticks, major_lat_ticks, major_lat_labels = calculate_ticks(overall_extent, 'latitude')
    ax.set_yticks(minor_lat_ticks, minor=True, crs=proj)
    ax.set_yticks(major_lat_ticks, crs=proj)
    ax.set_yticklabels(major_lat_labels, fontsize=fontsize)

    ax.tick_params(which='major', direction='out', bottom=True, top=True, labelbottom=label_bottom, labeltop=label_top, left=True, right=True, labelleft=label_left, labelright=label_right, length=5, width=2)

    ax.tick_params(which='minor', direction='out', bottom=True, top=True, left=True, right=True, width=1)

    if gridlines:
        gl = ax.gridlines(draw_labels=False, linewidth=0.25, color='black', alpha=0.25, linestyle='--', crs=proj, zorder=1000)
        gl.xlocator = mticker.FixedLocator(minor_lon_ticks)
        gl.ylocator = mticker.FixedLocator(minor_lat_ticks)

### FUNCTION:
def plot_bathymetry(ax, config, model_data, isobath1=-100, isobath2=-1000, downsample="auto", show_legend=False):
    
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
    - downsample (bool or str): Downsample bathymetry data.
        - default: "auto"
    - show_legend (bool): Show legend.
        - default: False

    Returns:
    - none
    '''

    bathymetry_path = config['DATA']['bathymetry_path']
    bathy_data = xr.open_dataset(bathymetry_path, chunks={'lat': 1000, 'lon': 1000})
    bathy_data = bathy_data.sel(lat=slice(model_data.lat.min(), model_data.lat.max()), lon=slice(model_data.lon.min(), model_data.lon.max()))
    
    if downsample == "auto":
        lat_min, lat_max = model_data.lat.min().item(), model_data.lat.max().item()
        lon_min, lon_max = model_data.lon.min().item(), model_data.lon.max().item()
        bathy_data = bathy_data.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))

        lat_range = lat_max - lat_min
        lon_range = lon_max - lon_min
        extent_area = lat_range * lon_range

        if extent_area > 250:
            downsample = True
        else:
            downsample = False

    if downsample:
        bathy_data = bathy_data.coarsen(lat=25, lon=25, boundary='trim').mean()
        ax.add_feature(cfeature.OCEAN, zorder=1)

    bathy_data = bathy_data.compute()

    isobath_levels = sorted([isobath1, isobath2])
    depth_intervals = [-np.inf] + isobath_levels + [0]

    bathy_data['elevation'] = bathy_data['elevation'].fillna(0).where(~np.isinf(bathy_data['elevation']), 0)

    cornflowerblue = mcolors.to_rgba('cornflowerblue')
    water = cfeature.COLORS['water']
    lightsteelblue = mcolors.to_rgba('lightsteelblue')
    colors = [cornflowerblue, water, lightsteelblue]

    ax.contour(bathy_data.lon, bathy_data.lat, bathy_data.elevation, levels=isobath_levels, colors='dimgrey', linestyles='dashed', linewidths=0.25, zorder=50, transform=ccrs.PlateCarree())

    for i in range(len(depth_intervals) - 1):
        ax.contourf(bathy_data.lon, bathy_data.lat, bathy_data.elevation, levels=[depth_intervals[i], depth_intervals[i + 1]], colors=[colors[i]], transform=ccrs.PlateCarree())
    
    if show_legend:
        isobath1 = -isobath1
        isobath2 = -isobath2
        legend_colors = [lightsteelblue, water, cornflowerblue]
        legend_labels = [f'0m - {isobath1}m', f'{isobath1}m - {isobath2}m', f'> {isobath2}m']
        patches = [plt.plot([], [], marker="o", ms=10, ls="", mec=None, color=color, label=label)[0] for color, label in zip(legend_colors, legend_labels)]
        bathymetry_legend = ax.legend(handles=patches, loc='upper left', facecolor='white', edgecolor='black', fontsize='small', markerscale=0.75)
        bathymetry_legend.set_zorder(10000)
        for text in bathymetry_legend.get_texts():
            text.set_color('black')

### FUNCTION:
def plot_streamlines(ax, longitude, latitude, u_depth_avg, v_depth_avg, density=2):
    
    '''
    Adds streamlines to the plot.
    
    Args:
    - ax (cartopy.mpl.geoaxes.GeoAxesSubplot): The cartopy map to add streamlines to.
    - longitude (array-like): Longitude values.
    - latitude (array-like): Latitude values.
    - u_depth_avg (array-like): U-component of depth-averaged currents.
    - v_depth_avg (array-like): V-component of depth-averaged currents.
    - density (int): Density of the streamlines.
        - default: 2

    Returns:
    - None
    '''

    streamplot = ax.streamplot(longitude, latitude, u_depth_avg, v_depth_avg, transform=ccrs.PlateCarree(), density=density, linewidth=0.5, color='black', zorder=10)
    streamplot.lines.set_alpha(1.0)

### FUNCTION:
def plot_magnitude_contour(ax, fig, longitude, latitude, mag_depth_avg, max_levels=10, extend_max=True):
    
    '''
    Plots a magnitude contour and adds a formatted color bar to the plot.

    Args:
    - ax (matplotlib.axes.Axes): The axes object to add the contour to.
    - fig (matplotlib.figure.Figure): The figure object for the plot.
    - longitude (array-like): Longitude values.
    - latitude (array-like): Latitude values.
    - mag_depth_avg (array-like or xarray.DataArray): Magnitude of depth-averaged current values.
    - max_levels (int): Maximum number of levels for contour.
        - default: 10
    - extend_max (bool): Whether to extend the maximum color level.
        - default: True

    Returns:
    - None
    '''

    levels, ticks, extend = format_contour_cbar(mag_depth_avg, max_levels=max_levels, extend_max=extend_max)
    
    contourf = ax.contourf(longitude, latitude, mag_depth_avg, levels=levels, cmap=cmo.speed, transform=ccrs.PlateCarree(), zorder=10, extend=extend)
    
    cbar = fig.colorbar(contourf, orientation='vertical', extend=extend, ax=ax)
    cbar.set_label('Depth Averaged Current Magnitude (m/s)', labelpad=10)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{tick:.1f}" for tick in ticks])
    format_cbar_position(ax, cbar)

### FUNCTION:
def plot_threshold_zones(ax, longitude, latitude, mag_depth_avg, mag1, mag2, mag3, mag4, mag5, threshold_legend=True):
    
    '''
    Adds threshold zones to the map.

    Args:
    - ax (cartopy.mpl.geoaxes.GeoAxesSubplot): The cartopy map to add the threshold zones to.
    - longitude (array-like): Longitude values.
    - latitude (array-like): Latitude values.
    - mag_depth_avg (array-like): Depth-averaged magnitude data.
    - mag1 (float): First threshold magnitude.
    - mag2 (float): Second threshold magnitude.
    - mag3 (float): Third threshold magnitude.
    - mag4 (float): Fourth threshold magnitude.
    - mag5 (float): Fifth threshold magnitude.
    - threshold_legend (bool): Show legend.
        - default: True
    
    Returns:
    - None
    '''

    max_mag = np.nanmax(mag_depth_avg)
    max_label = f'{max_mag:.2f}'

    levels = []
    colors = []
    labels = []

    if max_mag <= mag1:
        print("Warning: Threshold zone levels are too small to parse, no zones will be displayed")
        pass
    elif mag1 < max_mag <= mag2:
        print("Warning: Threshold zone levels are too small to parse, no zones will be displayed")
        pass
    elif mag2 < max_mag <= mag3:
        levels = [mag1, mag2, mag3]
        colors = ['none', 'yellow']
        labels = [None, f'{mag2} - {max_label} m/s']
    elif mag3 < max_mag <= mag4:
        levels = [mag1, mag2, mag3, mag4]
        colors = ['none', 'yellow', 'orange']
        labels = [None, f'{mag2} - {mag3} m/s', f'{mag3} - {max_label} m/s']
    elif mag4 < max_mag <= mag5:
        levels = [mag1, mag2, mag3, mag4, mag5]
        colors = ['none', 'yellow', 'orange', 'orangered']
        labels = [None, f'{mag2} - {mag3} m/s', f'{mag3} - {mag4} m/s', f'{mag4} - {max_label} m/s']
    else:
        levels = [mag1, mag2, mag3, mag4, mag5, max_mag]
        colors = ['none', 'yellow', 'orange', 'orangered', 'maroon', 'maroon']
        labels = [None, f'{mag2} - {mag3} m/s', f'{mag3} - {mag4} m/s', f'{mag4} - {mag5} m/s', f'{mag5} - {max_label} m/s']
    
    if levels:
        threshold_contourf = ax.contourf(longitude, latitude, mag_depth_avg, levels=levels, colors=colors[:len(levels)-1], extend='both', transform=ccrs.PlateCarree(), zorder=10)
    
    patches = []
    for color, label in zip(colors, labels):
        if label:
            patches.append(mpatches.Patch(color=color, label=label))

    if threshold_legend and labels:
        threshold_legend = ax.legend(handles=patches, loc='upper right', facecolor='white', edgecolor='black', fontsize='medium')
        threshold_legend.set_zorder(10000)
        for text in threshold_legend.get_texts():
            text.set_color('black')
        ax.add_artist(threshold_legend)

### FUNCTION:
def plot_advantage_zones(ax, config, longitude, latitude, dir_depth_avg, tolerance, advantage_legend=True):
    
    '''
    Adds advantage zones to the map.

    Args:
    - ax (cartopy.mpl.geoaxes.GeoAxesSubplot): The cartopy map to add the advantage zones to.
    - config (dict): Glider Guidance System mission configuration.
    - depth_average_data (xarray.core.dataset.Dataset): Depth-averaged model data.
    - tolerance (float): Tolerance for the advantage zones.
    - advantage_legend (bool): Show legend.
        - default: True

    Returns:
    - None
    '''

    GPS_coords = config['MISSION']['GPS_coords']
    if not GPS_coords or len(GPS_coords) < 2:
        print("Insufficient GPS route coordinates provided. Skipping advantage zone plotting.")
        return

    start_lat, start_lon = GPS_coords[0]
    end_lat, end_lon = GPS_coords[-1]

    direct_bearing = calculate_bearing(start_lat, start_lon, end_lat, end_lon)
    
    bearing_lower = (direct_bearing - tolerance) % 360
    bearing_upper = (direct_bearing + tolerance) % 360

    if bearing_lower < bearing_upper:
        mask = (dir_depth_avg >= bearing_lower) & (dir_depth_avg <= bearing_upper)
    else:
        mask = (dir_depth_avg >= bearing_lower) | (dir_depth_avg <= bearing_upper)
    
    acceptable_bearing = np.full_like(dir_depth_avg, np.nan)
    acceptable_bearing[mask] = dir_depth_avg[mask]

    advantage_contourf = ax.contourf(longitude, latitude, acceptable_bearing, levels=[bearing_lower, bearing_upper, 360], colors=['purple'], alpha=0.5, transform=ccrs.PlateCarree(), zorder=10)

    if advantage_legend:
        bearing_label_lower = round((direct_bearing - tolerance) % 360)
        bearing_label_upper = round((direct_bearing + tolerance) % 360)
        bearing_label = f"Advantage Zone: {bearing_label_lower}° to {bearing_label_upper}°"
        patches = [mpatches.Patch(color='purple', label=bearing_label)]
        advantage_legend = ax.legend(handles=patches, loc='upper left', facecolor='white', edgecolor='black', fontsize='medium')
        advantage_legend.set_zorder(10000)
        for text in advantage_legend.get_texts():
            text.set_color('black')
        ax.add_artist(advantage_legend)

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
def plot_optimal_path(ax, config, model_depth_average, optimal_path):

    '''
    Plots the optimal path for the mission between all successive waypoints given in the GPS_coords list.

    Args:
    - ax (matplotlib.axes.Axes): Matplotlib axes object.
    - config (dict): Configuration dictionary.
    - model_depth_average (xarray.Dataset): Depth-averaged ocean current dataset.
    - optimal_path (list): List of optimal path coordinates.

    Returns:
    - None
    '''

    if optimal_path:
        route_lats, route_lons = zip(*optimal_path)
        ax.plot(route_lons, route_lats, 'o-', transform=ccrs.PlateCarree(), markersize=5, linewidth=3, color='black', zorder=94)
    else:
        print("Invalid GPS waypoint list provided. Skipping optimal path plotting.")
        return

### FUNCTION:
def profile_rtofs(axs, dataset, latitude_qc, longitude_qc, threshold):
    
    '''
    Plot the RTOFS model data profiles for a given latitude and longitude.

    Args:
    - axs (list): List of matplotlib axes objects.
    - dataset (tuple): Tuple containing the RTOFS model data, depth-averaged data, and bin-averaged data.
    - latitude_qc (float): Target latitude value.
    - longitude_qc (float): Target longitude value.
    - threshold (float): Threshold value for shading.

    Returns:
    - None
    '''

    rtofs_model_data, rtofs_depth_average, rtofs_bin_average = dataset

    (y_rtofs_model, x_rtofs_model), _ = calculate_gridpoint(rtofs_model_data, latitude_qc, longitude_qc)
    (y_rtofs_avg, x_rtofs_avg), _ = calculate_gridpoint(rtofs_depth_average, latitude_qc, longitude_qc)
    (y_rtofs_bin, x_rtofs_bin), _ = calculate_gridpoint(rtofs_bin_average, latitude_qc, longitude_qc)
    
    rtofs_model_u = rtofs_model_data['u'].isel(y=y_rtofs_model, x=x_rtofs_model).values
    rtofs_avg_u = rtofs_depth_average['u_depth_avg'].isel(y=y_rtofs_avg, x=x_rtofs_avg).values
    rtofs_avg_u = rtofs_avg_u.item()
    rtofs_bin_u = rtofs_bin_average['u_bin_avg'].isel(y=y_rtofs_bin, x=x_rtofs_bin).values
    rtofs_bin_u1d = rtofs_bin_u[0]

    rtofs_model_v = rtofs_model_data['v'].isel(y=y_rtofs_model, x=x_rtofs_model).values
    rtofs_avg_v = rtofs_depth_average['v_depth_avg'].isel(y=y_rtofs_avg, x=x_rtofs_avg).values
    rtofs_avg_v = rtofs_avg_v.item()
    rtofs_bin_v = rtofs_bin_average['v_bin_avg'].isel(y=y_rtofs_bin, x=x_rtofs_bin).values
    rtofs_bin_v1d = rtofs_bin_v[0]

    rtofs_avg_mag = rtofs_depth_average['mag_depth_avg'].isel(y=y_rtofs_avg, x=x_rtofs_avg).values
    rtofs_avg_mag = rtofs_avg_mag.item()
    rtofs_bin_mag = rtofs_bin_average['mag_bin_avg'].isel(y=y_rtofs_bin, x=x_rtofs_bin).values
    rtofs_bin_mag1d = rtofs_bin_mag[0]

    rtofs_bin_dir = rtofs_bin_average['dir_bin_avg'].isel(y=y_rtofs_bin, x=x_rtofs_bin).values
    rtofs_bin_dir = np.mod(rtofs_bin_dir, 360)

    rtofs_avg_dir = rtofs_depth_average['dir_depth_avg'].isel(y=y_rtofs_avg, x=x_rtofs_avg).values
    rtofs_avg_dir = rtofs_avg_dir.item()
    rtofs_avg_dir = np.mod(rtofs_avg_dir, 360)

    rtofs_max_depth = rtofs_model_data.depth.max().item()
    rtofs_max_bins = rtofs_max_depth + 1
    rtofs_bin_depths = np.arange(rtofs_max_bins)

    ax = axs[0]
    ax.scatter(rtofs_model_u, rtofs_model_data.depth, marker='x', color='black', s=100, label='Model Datapoint', alpha=1.0, zorder=3)
    ax.scatter(rtofs_bin_u, rtofs_bin_depths, label='1m Interpolation', color='cyan', alpha=1.0, zorder=2)
    ax.axvline(x=rtofs_avg_u, label=f'Depth Average = [{rtofs_avg_u:.2f}]', color='darkcyan', linestyle='--', linewidth=2, zorder=1)
    ax.set_xlabel('u Velocity (m/s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)

    ax = axs[1]
    ax.scatter(rtofs_model_v, rtofs_model_data.depth, label='Model Datapoint', marker='x', color='black', s=100, alpha=1.0, zorder=3)
    ax.scatter(rtofs_bin_v, rtofs_bin_depths, label='1m Interpolation', color='orange', alpha=1.0, zorder=2)
    ax.axvline(x=rtofs_avg_v, label=f'Depth Avgerage = [{rtofs_avg_v:.2f}]', color='darkorange', linestyle='--', linewidth=2, zorder=1)
    ax.set_xlabel('v Velocity (m/s)', fontsize=12, fontweight='bold')
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)

    ax = axs[2]
    ax.scatter(rtofs_bin_mag, rtofs_bin_depths, label='1m Interpolation', color='green', alpha=1.0, zorder=2)
    ax.set_xlabel('Current Magnitude (m/s)', fontsize=12, fontweight='bold')
    ax.axvline(x=rtofs_avg_mag, label=f'Depth Avgerage = [{rtofs_avg_mag:.2f}]', color='darkgreen', linestyle='--', linewidth=2, zorder=1)
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)

    ax = axs[3]
    ax.scatter(rtofs_bin_dir, rtofs_bin_depths, label='1m Interpolation', color='purple', alpha=1.0, zorder=2)
    ax.axvline(x=rtofs_avg_dir, label=f'Depth Average = [{rtofs_avg_dir:.2f}]', color='darkviolet', linestyle='--', linewidth=2, zorder=1)
    ax.set_xlabel('Current Direction (degrees)', fontsize=12, fontweight='bold')
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)

    shading_colors = ['cyan', 'orange', 'green']
    for i, (data_1d, color) in enumerate(zip([rtofs_bin_u1d, rtofs_bin_v1d, rtofs_bin_mag1d], shading_colors)):
        plot_profile_thresholds(axs[i], data_1d, threshold, color)

    axs[0].set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
    for ax in axs:
        ax.invert_yaxis()
        ax.tick_params(axis='x', labelrotation=45)
        ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(handles, labels, loc='lower center', facecolor='lightgrey', edgecolor='black', framealpha=1.0, fontsize=14)

### FUNCTION:
def profile_cmems(axs, dataset, latitude_qc, longitude_qc, threshold):
    
    '''
    Plot the CMEMS model data profiles for a given latitude and longitude.

    Args:
    - axs (list): List of matplotlib axes objects.
    - dataset (tuple): Tuple containing the CMEMS model data, depth-averaged data, and bin-averaged data.
    - latitude_qc (float): Target latitude value.
    - longitude_qc (float): Target longitude value.
    - threshold (float): Threshold value for shading.

    Returns:
    - None
    '''

    cmems_model_data, cmems_depth_average, cmems_bin_average = dataset
        
    (y_cmems_model, x_cmems_model), _ = calculate_gridpoint(cmems_model_data, latitude_qc, longitude_qc)
    (y_cmems_avg, x_cmems_avg), _ = calculate_gridpoint(cmems_depth_average, latitude_qc, longitude_qc)
    (y_cmems_bin, x_cmems_bin), _ = calculate_gridpoint(cmems_bin_average, latitude_qc, longitude_qc)

    cmems_model_u = cmems_model_data['u'].isel(lat=y_cmems_model, lon=x_cmems_model).values
    cmems_avg_u = cmems_depth_average['u_depth_avg'].isel(lat=y_cmems_avg, lon=x_cmems_avg).values
    cmems_avg_u = cmems_avg_u.item()
    cmems_bin_u = cmems_bin_average['u_bin_avg'].isel(lat=y_cmems_bin, lon=x_cmems_bin).values
    cmems_bin_u1d = cmems_bin_u[0]

    cmems_model_v = cmems_model_data['v'].isel(lat=y_cmems_model, lon=x_cmems_model).values
    cmems_avg_v = cmems_depth_average['v_depth_avg'].isel(lat=y_cmems_avg, lon=x_cmems_avg).values
    cmems_avg_v = cmems_avg_v.item()
    cmems_bin_v = cmems_bin_average['v_bin_avg'].isel(lat=y_cmems_bin, lon=x_cmems_bin).values
    cmems_bin_v1d = cmems_bin_v[0]

    cmems_avg_mag = cmems_depth_average['mag_depth_avg'].isel(lat=y_cmems_avg, lon=x_cmems_avg).values
    cmems_avg_mag = cmems_avg_mag.item()
    cmems_bin_mag = cmems_bin_average['mag_bin_avg'].isel(lat=y_cmems_bin, lon=x_cmems_bin).values
    cmems_bin_mag1d = cmems_bin_mag[0]

    cmems_bin_dir = cmems_bin_average['dir_bin_avg'].isel(lat=y_cmems_bin, lon=x_cmems_bin).values
    cmems_bin_dir = np.mod(cmems_bin_dir, 360)

    cmems_avg_dir = cmems_depth_average['dir_depth_avg'].isel(lat=y_cmems_bin, lon=x_cmems_bin).values
    cmems_avg_dir = cmems_avg_dir.item()
    cmems_avg_dir = np.mod(cmems_avg_dir, 360)

    cmems_max_depth = cmems_model_data.depth.max().item()
    cmems_max_bins = cmems_max_depth + 1
    cmems_bin_depths = np.arange(cmems_max_bins)

    ax = axs[0]
    ax.scatter(cmems_model_u, cmems_model_data.depth, marker='x', color='black', s=100, label='Model Datapoint', alpha=1.0, zorder=3)
    ax.scatter(cmems_bin_u, cmems_bin_depths, label='1m Interpolation', color='cyan', alpha=1.0, zorder=2)
    ax.axvline(x=cmems_avg_u, label=f'Depth Average = [{cmems_avg_u:.2f}]', color='darkcyan', linestyle='--', linewidth=2, zorder=1)
    ax.set_xlabel('u Velocity (m/s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)

    ax = axs[1]
    ax.scatter(cmems_model_v, cmems_model_data.depth, label='Model Datapoint', marker='x', color='black', s=100, alpha=1.0, zorder=3)
    ax.scatter(cmems_bin_v, cmems_bin_depths, label='1m Interpolation', color='orange', alpha=1.0, zorder=2)
    ax.axvline(x=cmems_avg_v, label=f'Depth Avgerage = [{cmems_avg_v:.2f}]', color='darkorange', linestyle='--', linewidth=2, zorder=1)
    ax.set_xlabel('v Velocity (m/s)', fontsize=12, fontweight='bold')
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    
    ax = axs[2]
    ax.scatter(cmems_bin_mag, cmems_bin_depths, label='1m Interpolation', color='green', alpha=1.0, zorder=2)
    ax.set_xlabel('Current Magnitude (m/s)', fontsize=12, fontweight='bold')
    ax.axvline(x=cmems_avg_mag, label=f'Depth Avgerage = [{cmems_avg_mag:.2f}]', color='darkgreen', linestyle='--', linewidth=2, zorder=1)
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)

    ax = axs[3]
    ax.scatter(cmems_bin_dir, cmems_bin_depths, label='1m Interpolation', color='purple', alpha=1.0, zorder=2)
    ax.axvline(x=cmems_avg_dir, label=f'Depth Average = [{cmems_avg_dir:.2f}]', color='darkviolet', linestyle='--', linewidth=2, zorder=1)
    ax.set_xlabel('Current Direction (degrees)', fontsize=12, fontweight='bold')
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)

    shading_colors = ['cyan', 'orange', 'lawngreen']
    for i, (data_1d, color) in enumerate(zip([cmems_bin_u1d, cmems_bin_v1d, cmems_bin_mag1d], shading_colors)):
        plot_profile_thresholds(axs[i], data_1d, threshold, color)

    axs[0].set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
    for ax in axs:
        ax.invert_yaxis()
        ax.tick_params(axis='x', labelrotation=45)
        ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(handles, labels, loc='lower center', facecolor='lightgrey', edgecolor='black', framealpha=1.0, fontsize=14)

### FUNCTION:
def profile_gofs(axs, dataset, latitude_qc, longitude_qc, threshold):
    
    '''
    Plot the GOFS model data profiles for a given latitude and longitude.

    Args:
    - axs (list): List of matplotlib axes objects.
    - dataset (tuple): Tuple containing the GOFS model data, depth-averaged data, and bin-averaged data.
    - latitude_qc (float): Target latitude value.
    - longitude_qc (float): Target longitude value.
    - threshold (float): Threshold value for shading.

    Returns:
    - None
    '''

    gofs_model_data, gofs_depth_average, gofs_bin_average = dataset
        
    (y_gofs_model, x_gofs_model), _ = calculate_gridpoint(gofs_model_data, latitude_qc, longitude_qc)
    (y_gofs_avg, x_gofs_avg), _ = calculate_gridpoint(gofs_depth_average, latitude_qc, longitude_qc)
    (y_gofs_bin, x_gofs_bin), _ = calculate_gridpoint(gofs_bin_average, latitude_qc, longitude_qc)

    gofs_model_u = gofs_model_data['u'].isel(lat=y_gofs_model, lon=x_gofs_model).values
    gofs_avg_u = gofs_depth_average['u_depth_avg'].isel(lat=y_gofs_avg, lon=x_gofs_avg).values
    gofs_avg_u = gofs_avg_u.item()
    gofs_bin_u = gofs_bin_average['u_bin_avg'].isel(lat=y_gofs_bin, lon=x_gofs_bin).values
    gofs_bin_u1d = gofs_bin_u[0]

    gofs_model_v = gofs_model_data['v'].isel(lat=y_gofs_model, lon=x_gofs_model).values
    gofs_avg_v = gofs_depth_average['v_depth_avg'].isel(lat=y_gofs_avg, lon=x_gofs_avg).values
    gofs_avg_v = gofs_avg_v.item()
    gofs_bin_v = gofs_bin_average['v_bin_avg'].isel(lat=y_gofs_bin, lon=x_gofs_bin).values
    gofs_bin_v1d = gofs_bin_v[0]

    gofs_avg_mag = gofs_depth_average['mag_depth_avg'].isel(lat=y_gofs_avg, lon=x_gofs_avg).values
    gofs_avg_mag = gofs_avg_mag.item()
    gofs_bin_mag = gofs_bin_average['mag_bin_avg'].isel(lat=y_gofs_bin, lon=x_gofs_bin).values
    gofs_bin_mag1d = gofs_bin_mag[0]

    gofs_bin_dir = gofs_bin_average['dir_bin_avg'].isel(lat=y_gofs_bin, lon=x_gofs_bin).values
    gofs_bin_dir = np.mod(gofs_bin_dir, 360)

    gofs_avg_dir = gofs_depth_average['dir_depth_avg'].isel(lat=y_gofs_avg, lon=x_gofs_avg).values
    gofs_avg_dir = gofs_avg_dir.item()
    gofs_avg_dir = np.mod(gofs_avg_dir, 360)

    gofs_max_depth = gofs_model_data.depth.max().item()
    gofs_max_bins = gofs_max_depth + 1
    gofs_bin_depths = np.arange(gofs_max_bins)

    ax = axs[0]
    ax.scatter(gofs_model_u, gofs_model_data.depth, marker='x', color='black', s=100, label='Model Datapoint', alpha=1.0, zorder=3)
    ax.scatter(gofs_bin_u, gofs_bin_depths, label='1m Interpolation', color='cyan', alpha=1.0, zorder=2)
    ax.axvline(x=gofs_avg_u, label=f'Depth Average = [{gofs_avg_u:.2f}]', color='darkcyan', linestyle='--', linewidth=2, zorder=1)
    ax.set_xlabel('u Velocity (m/s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)

    ax = axs[1]
    ax.scatter(gofs_model_v, gofs_model_data.depth, label='Model Datapoint', marker='x', color='black', s=100, alpha=1.0, zorder=3)
    ax.scatter(gofs_bin_v, gofs_bin_depths, label='1m Interpolation', color='orange', alpha=1.0, zorder=2)
    ax.axvline(x=gofs_avg_v, label=f'Depth Avgerage = [{gofs_avg_v:.2f}]', color='darkorange', linestyle='--', linewidth=2, zorder=1)
    ax.set_xlabel('v Velocity (m/s)', fontsize=12, fontweight='bold')
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)

    ax = axs[2]
    ax.scatter(gofs_bin_mag, gofs_bin_depths, label='1m Interpolation', color='green', alpha=1.0, zorder=2)
    ax.set_xlabel('Current Magnitude (m/s)', fontsize=12, fontweight='bold')
    ax.axvline(x=gofs_avg_mag, label=f'Depth Avgerage = [{gofs_avg_mag:.2f}]', color='darkgreen', linestyle='--', linewidth=2, zorder=1)
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    
    ax = axs[3]
    ax.scatter(gofs_bin_dir, gofs_bin_depths, label='1m Interpolation', color='purple', alpha=1.0, zorder=2)
    ax.axvline(x=gofs_avg_dir, label=f'Depth Average = [{gofs_avg_dir:.2f}]', color='darkviolet', linestyle='--', linewidth=2, zorder=1)
    ax.set_xlabel('Current Direction (degrees)', fontsize=12, fontweight='bold')
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)

    shading_colors = ['cyan', 'orange', 'lawngreen']
    for i, (data_1d, color) in enumerate(zip([gofs_bin_u1d, gofs_bin_v1d, gofs_bin_mag1d], shading_colors)):
        plot_profile_thresholds(axs[i], data_1d, threshold, color)

    axs[0].set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
    for ax in axs:
        ax.invert_yaxis()
        ax.tick_params(axis='x', labelrotation=45)
        ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(handles, labels, loc='lower center', facecolor='lightgrey', edgecolor='black', framealpha=1.0, fontsize=14)

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
        glider_legend = ax.legend(handles=legend_handles, title="Gliders", loc='lower right', facecolor='white', edgecolor='black', fontsize='small', markerscale=0.75)
        glider_legend.set_zorder(10000)
        for text in glider_legend.get_texts():
            text.set_color('black')

### FUNCTION:
def plot_add_eez(ax, config, color='dimgrey', linewidth=3, zorder=90):
    
    '''
    Adds Exclusive Economic Zones (EEZ) to a cartopy map.

    Args:
    - ax (cartopy.mpl.geoaxes.GeoAxesSubplot): The cartopy map to add the EEZ to.
    - config (dict): Glider Guidance System mission configuration.
    - color (str): Color of the EEZ border.
        - default: 'dimgrey'
    - linewidth (int): Width of the EEZ border.
        - default: 3
    - zorder (int): Z-order of the EEZ border.
        - default: 90
    
    Returns:
    - None
    '''

    eez_path = config['DATA']['eez_path']

    eez_feature = cfeature.ShapelyFeature(
        Reader(eez_path).geometries(),
        ccrs.PlateCarree(),
        edgecolor=color,
        facecolor='none',
        linewidth=linewidth,
        linestyle='-',
        alpha=0.75
    )
    
    ax.add_feature(eez_feature, zorder=zorder)

### FUNCTION:
def plot_glider_route(ax, config):
    
    '''
    Adds the glider route to a cartopy map.

    Args:
    - ax (cartopy.mpl.geoaxes.GeoAxesSubplot): The cartopy map to add the glider route to.
    - config (dict): Glider Guidance System mission configuration.

    Returns:
    - None
    '''

    if config['MISSION']['GPS_coords'] is not None:
        GPS_coords = config['MISSION']['GPS_coords']
        num_waypoints = len(GPS_coords)

        for idx, (lat, lon) in enumerate(GPS_coords):
            if idx == 0:
                color = 'green'
                label = 'Start'
            elif idx == num_waypoints - 1:
                color = 'red'
                label = 'End'
            else:
                color = 'purple'
                label = None

            ax.scatter(lon, lat, color=color, s=100, transform=ccrs.PlateCarree(), zorder=95, label=label)

        if num_waypoints > 2:
            ax.scatter([], [], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=95, label='Intermediate')
    else:
        print("Invalid GPS waypoint list provided. Skipping glider route plotting.")
        return

# FORMATTING FUNCTIONS

### FUNCTION:
def format_contour_cbar(magnitude, max_levels=10, extend_max=True):
    
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
def format_cbar_position(ax, cbar):
    
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
    
    title_distance_top = 0.15
    title_position = top + title_distance_top
    
    subtitle_distance_title = 0.05
    subtitle_position = title_position - subtitle_distance_title
    
    suptitle_distance_bottom = 0.1
    suptitle_position = bottom - suptitle_distance_bottom

    title_datetime = format_title_datetime(datetime_index)
    
    fig.text(0.5, title_position, title, fontsize=20, fontweight='bold', ha='center', va='bottom')

    subtitle_text = f"{model_name} {title_datetime} UTC"
    fig.text(0.5, subtitle_position, subtitle_text, fontsize=18, ha='center', va='bottom')

    suptitle_text = f"Generated by the Glider Guidance System (GGS) - {config['MISSION']['mission_name']}"
    fig.text(0.5, suptitle_position, suptitle_text, fontsize=16, fontweight='bold', ha='center', va='top', color='gray')

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

    fig.suptitle(title, fontsize=26, fontweight='bold', ha='center', y=1.030)
    
    title_datetime = format_title_datetime(datetime)
    subtitle_text = f"{title_datetime} UTC"
    fig.text(0.5, 0.950, subtitle_text, fontsize=22, ha='center', va='top')
    
    suptitle_text = f"Generated by the Glider Guidance System (GGS) - {config['MISSION']['mission_name']}"
    fig.text(0.5, -0.015, suptitle_text, fontsize=20, ha='center', va='top', color='gray')

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
        text_x = (bbox.x0 + bbox.x1) / 2
        text_y = bbox.y1 + bbox.height * 0.025
        
        fig.text(text_x, text_y, model_name, fontsize=16, fontweight='bold', ha='center', va='bottom')

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
def convert_DD_to_DMS(DD):
    
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
