# =========================
# X - IMPORTS
# =========================

### /// FUNCTIONS ///
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import numpy as np

# =========================
# CHECK ABORT
# =========================

### FUNCTION:
EXIT_KEYWORD = "EXIT"
def check_abort(input_string):
    
    '''
    Check if the input string matches the exit keyword. If it does, abort the program.
    
    Args:
    - input_string (str): The input string to be checked.

    Returns:
    - None
    '''
    
    if input_string.upper() == EXIT_KEYWORD:
        print("Configuration aborted by user.")
        exit()

# =========================
# CHECK FLOAT
# =========================

### FUNCTION:
def check_float(prompt_msg):
    
    '''
    Repeatedly prompt the user until they provide a valid float input.
    
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

# =========================
# CHECK COORDINATE
# =========================

### FUNCTION:
def check_coordinate(coord_str):
    
    '''
    Validate if the provided string can be parsed into valid latitude and longitude coordinates.
    
    Args:
    - coord_str (str): The input string containing coordinates.

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
# CALCULATE DISTANCE
# =========================

### FUNCTION:
def calculate_distance(coord1, coord2):
    
    '''
    Calculate the Haversine distance between two sets of GPS coordinates.
    
    Args:
    - coord1 (tuple): Starting coordinate as (latitude, longitude).
    - coord2 (tuple): Ending coordinate as (latitude, longitude).

    Returns:
    - float: The distance between the two coordinates in meters.
    '''
    
    R = 6371000
    lat1, lon1 = np.radians(coord1)
    lat2, lon2 = np.radians(coord2)
    
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    a = (np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    return R * c

# =========================
# CALCULATE NEARPOINT
# =========================

### FUNCTION:
def calculate_nearpoint(dataset, target_lat, target_lon):
    
    '''
    Find the nearest point in the dataset to the input latitude and longitude.

    Args:
    - dataset (xarray.Dataset): The dataset to search in.
    - target_lat (float): The target latitude.
    - target_lon (float): The target longitude.

    Returns:
    - (y_index, x_index) (tuple): The indices of the nearest point in the dataset.
    - (lat_index, lon_index) (tuple): The coordinates of the nearest point in the dataset.
    '''

    lat_diff = dataset['lat'] - float(target_lat)
    lon_diff = dataset['lon'] - float(target_lon)
    distance_square = lat_diff**2 + lon_diff**2

    y_index, x_index = np.unravel_index(distance_square.argmin(), distance_square.shape)

    lat_index = dataset['lat'].isel(y=y_index, x=x_index).values
    lon_index = dataset['lon'].isel(y=y_index, x=x_index).values

    print(f"Input Coordinates: ({target_lat}, {target_lon})")
    print(f"Dataset Indices: ({y_index}, {x_index})")
    print(f"Dataset Coordinates: ({lat_index:.3f}, {lon_index:.3f})")

    return (y_index, x_index), (lat_index, lon_index)

# =========================
# CALCULATE HEADING
# =========================

### FUNCTION:
def calculate_heading(coord1, coord2):
    
    '''
    Calculate the initial compass bearing/heading from one point to another on Earth.
    
    Args:
    - coord1 (tuple): Starting coordinate as (latitude, longitude).
    - coord2 (tuple): Ending coordinate as (latitude, longitude).

    Returns:
    - float: The initial heading/bearing in degrees.
    '''

    lat1, lon1 = np.radians(coord1)
    lat2, lon2 = np.radians(coord2)
    
    dlon = lon2 - lon1
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - (np.sin(lat1) * np.cos(lat2) * np.cos(dlon))
    
    initial_bearing = np.degrees(np.arctan2(x, y))
    return (initial_bearing + 360) % 360

# =========================
# DD TO DMS
# =========================

### FUNCTION:
def DD_to_DMS(DD, coord_type='longitude'):

    '''
    Convert a decimal degree position to a degree-minute-second string.

    Args:
    - DD (float): A decimal degree position.
    - coord_type (str): A string indicating whether the position is a longitude or latitude.

    Returns:
    - str: A degree-minute-second string.
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

# =========================
# DD TO DM
# =========================

### FUNCTION:
def DD_to_DM(DD, coord_type='longitude'):
    
    '''
    Convert a decimal degree position to a degree-minute string.

    Args:
    - DD (float): A decimal degree position.
    - coord_type (str): A string indicating whether the position is a longitude or latitude.

    Returns:
    - str: A degree-minute string.
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

# =========================
# DD TO DDM
# =========================

### FUNCTION:
def DD_to_DDM(DD):
    
    '''
    Convert a decimal degree (DD) coordinate to a degree-minute (DDMM) coordinate.

    Args:
    - DD (float): A decimal degree coordinate.

    Returns:
    - str: A degree-minute coordinate.    
    '''
    
    degrees = int(DD)
    minutes = abs(DD - degrees) * 60

    return f"{degrees:02d}{minutes:05.2f}"

# =========================
# SET TICKMARKS
# =========================

### FUNCTION:
def set_ticks(ax, extent_lon, extent_lat):
    
    '''
    Set the ticks for a map based on the extent without changing the map's bounds, ensuring at least three ticks.

    Args:
    - ax (matplotlib.axes): The axes to set the ticks for.
    - extent_lon (array-like): The longitude extent of the map.
    - extent_lat (array-like): The latitude extent of the map.

    Returns:
    - None
    '''

    ### FUNCTION:
    def get_tick_interval(extent):

        '''
        Calculate the tick interval for a given extent.

        Args:
        - extent (array-like): The extent of the map.

        Returns:
        - float: The tick interval.
        '''

        min_extent, max_extent = np.min(extent), np.max(extent)
        range_extent = max_extent - min_extent

        intervals = [0.25, 0.5, 1, 2, 5, 10, 15, 30]
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

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda val, pos: DD_to_DM(val, 'longitude')))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda val, pos: DD_to_DM(val, 'latitude')))
