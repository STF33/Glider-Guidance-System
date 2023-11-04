# =========================
# 0 - Imports
# =========================

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader
from cartopy.io.img_tiles import GoogleTiles
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import cmocean
import cool_maps.plot as cplt
import dask
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
# import plotly.graph_objects as go
# import plotly.io as pio
import xarray as xr

from sup_functions import check_abort, check_float, check_coordinate
from sup_functions import calculate_distance, calculate_heading, DD_to_DDMM

from sup_slocum import goto_l10

# =========================
# 1 - GGS Config
# =========================

### FUNCTION:           DEBUG           DEBUG           DEBUG           DEBUG           DEBUG
def GGS_config_static():
    
    """
    Return a hardcoded Glider Guidance System configuration for debugging.
    
    Returns:
    - tuple: A tuple containing the hardcoded configuration dictionary and a list of all waypoints.
    """
    
    config = {
        "glider_name": "RU00", #
        "max_depth": 100.0, #
        "avg_velocity": 0.225, #
        "battery_capacity": 150, #
        "battery_drain": 5, #
        "satisfying_radius": 1000, #
        "waypoints": [(20.058417927338272, -85.59492605234463),
                      (21.96308733263867, -85.96115990509725),
                      (23.315103962610003, -87.9492865343257),
                      (22.9462365115331, -90.8791573563466),
                      (24.08161230606445, -95.23908417483008),
                      (29.231215967007177, -94.70901167378757)
                      ] # Yucatan
        # "waypoints": [(39.457283, -74.193117),
        #               (39.29305, -73.94345),
        #               (39.129067, -73.665433),
        #               (38.93165, -73.343783),
        #               (38.767433, -73.078417),
        #               (38.93165, -73.343783),
        #               (39.129067, -73.665433),
        #               (39.29305, -73.94345)
        #               ] # E-Line
    }

    waypoints = config["waypoints"]

    GGS_config_output(config)

    return config, waypoints

### FUNCTION:
def GGS_config_import():
    
    """
    Import an existing Glider Guidance System configuration.
    
    Returns:
    - tuple: A tuple containing the configuration dictionary and a list of all waypoints.
    """
    
    while True:
        file_name = input("Please enter the path to your GGS config file: ")
        check_abort(file_name)
        try:
            config = pd.read_pickle(file_name)
            waypoints = config["waypoints"]
            GGS_config_output(config)
            return config, waypoints
        except Exception as e:
            print(f"Error reading the file: {e}")
            retry = input("Invalid file. Would you like to try another file? (yes or no): ").lower()
            check_abort(retry)
            if retry != "yes":
                print("Import aborted. Please restart the configuration process.")
                return None

### FUNCTION:
def GGS_config_new():
    
    """
    Create a new Glider Guidance System configuration.

    Returns:
    - tuple: A tuple containing the configuration dictionary and a list of all waypoints.
    """

    config = {}

    prompt_glider_name = "Enter the glider's name: "
    while True:
        glider_name = input(prompt_glider_name)
        check_abort(glider_name)
        if glider_name:
            config["glider_name"] = glider_name
            break
        prompt_glider_name = "[INPUT ERROR] " + prompt_glider_name

    config["max_depth"] = check_float("Enter the maximum mission depth (meters): ")
    config["avg_velocity"] = check_float("Enter the average velocity (meters per second): ")
    config["battery_capacity"] = check_float("Enter the maximum battery capacity (amp-hours): ")
    config["battery_drain"] = check_float("Enter the average battery drain rate (amp-hours per day): ")
    config["satisfying_radius"] = check_float("Enter the satisfying radius (meters): ")

    prompt_num_waypoints = "Enter the total number of waypoints (1-100): "
    waypoints = []
    while True:
        try:
            num_waypoints = int(input(prompt_num_waypoints))
            check_abort(str(num_waypoints))
            if 1 <= num_waypoints <= 100:
                break
        except ValueError:
            pass
        prompt_num_waypoints = "[INPUT ERROR] " + prompt_num_waypoints

    for i in range(num_waypoints):
        waypoint_prompt = f"Enter waypoint {i} coordinates as 'latitude, longitude': "
        while True:
            waypoint = check_coordinate(input(waypoint_prompt))
            if waypoint:
                waypoints.append(waypoint)
                break
            waypoint_prompt = "[INPUT ERROR] " + waypoint_prompt

    config["waypoints"] = waypoints

    GGS_config_output(config)
    return config, waypoints

### FUNCTION:
def GGS_config():
    
    """
    Import an existing Glider Guidance System configuration or guide the user to create a new one.
    
    Returns:
    - tuple: A tuple containing the configuration dictionary and a list of all waypoints.
    """
    
    print(f"Glider Guidance System (GGS) Configuration Setup\n(Note: Type '{EXIT_KEYWORD}' at any time to abort the configuration.)\n")
    choice = input("Do you want to import an existing GGS config or create a new one? (Enter 'import' or 'new'): ").lower()
    check_abort(choice)
    
    if choice == "import":
        return GGS_config_import()
    elif choice == "new":
        return GGS_config_new()
    else:
        print("Invalid choice. Please enter 'import' or 'new'.")
        return GGS_config()

### FUNCTION:
def GGS_config_output(config):
    
    """
    Display the provided Glider Guidance System configuration in a formatted manner and save to a file.
    
    Args:
    - config (dict): The configuration dictionary.

    Returns:
    - directory (str): The directory path where the configuration is saved.
    """

    output_str = "\nGlider Guidance System (GGS) Configuration:\n\n"

    for key, value in config.items():
        if key == "waypoints":
            for i, coord in enumerate(value, 1):
                output_str += "Waypoint {}: {}\n".format(i, coord)
        elif key == "glider_name":
            output_str += "{}: {}\n\n".format(key.capitalize(), value)
        elif key in ["max_depth", "avg_velocity", "battery_capacity", "battery_drain", "satisfying_radius"]:
            formatted_key = key.capitalize().replace('_', ' ')
            output_str += "{}: {}\n".format(formatted_key, value)

    directory = os.path.join(os.path.expanduser("~"), "Downloads", f"GGS_{config['glider_name']}")
    os.makedirs(directory, exist_ok=True)

    config_pickle = os.path.join(directory, f"GGS_{config['glider_name']}_config.pkl")
    pd.to_pickle(config, config_pickle)

    config_text = os.path.join(directory, f"GGS_{config['glider_name']}_config.txt")
    with open(config_text, 'w') as file:
        file.write(output_str)

    print(output_str)

    return directory

# =========================
# 2 - Route Calculations
# =========================

### FUNCTION:
def route_analysis(config, waypoints):
    
    """
    Analyze each leg of the route to compute metrics such as distance, time, and battery drain.
    
    Args:
    - config (dict): The GGS configuration dictionary.
    - waypoints (list): A list of coordinates representing the mission route.

    Returns:
    - list: A list of dictionaries, each containing data for a specific leg of the route.
    """
    
    results = []
    
    for i in range(len(waypoints) - 1):
        leg_distance = calculate_distance(waypoints[i], waypoints[i + 1])
        leg_time_seconds = leg_distance / config['avg_velocity']
        leg_time_hours = leg_time_seconds / 3600
        battery_drain = config['battery_drain'] * (leg_time_seconds / 86400)
        
        leg_data = {
            "Leg Start": waypoints[i],
            "Leg End": waypoints[i + 1],
            "Leg Distance (meters)": leg_distance,
            "Time Taken (seconds)": leg_time_seconds,
            "Time Taken (hours)": leg_time_hours,
            "Expected Battery Drain (amp-hours)": battery_drain
        }
        
        results.append(leg_data)
    
    return results

### FUNCTION:
def route_analysis_output(config, directory, analysis_results):
    
    """
    Print a summary of the route analysis and save the details in a text file.
    
    Args:
    - config (dict): The GGS configuration dictionary.
    - directory (str): The directory path where the configuration is saved.
    - analysis_results (list): List of dictionaries containing data for each leg of the route.

    Returns:
    - total_distance (float): Total distance traveled in the route.
    - total_time_seconds (float): Total time taken for the route in seconds.
    - total_time_hours (float): Total time taken for the route in hours.
    - total_battery_drain (float): Total expected battery drain for the route.
    """
    
    output_str = ""

    point_names = [f"Waypoint {i+1}" for i in range(len(config['waypoints']))]
    
    total_distance = 0
    total_time_seconds = 0
    total_time_hours = 0
    total_battery_drain = 0

    for idx, result in enumerate(analysis_results):
        leg_description = f"Leg {idx+1}: {point_names[idx]} {result['Leg Start']} to {point_names[idx+1]} {result['Leg End']}"
        leg_heading = calculate_heading(result['Leg Start'], result['Leg End'])
        
        output_str += (
            f"{leg_description}\n"
            f"Heading (degrees): {leg_heading}\n"
            f"Leg Distance (meters): {result['Leg Distance (meters)']}\n"
            f"Time Taken (seconds): {result['Time Taken (seconds)']}\n"
            f"Time Taken (hours): {result['Time Taken (hours)']}\n"
            f"Expected Battery Drain (amp-hours): {result['Expected Battery Drain (amp-hours)']}\n\n"
        )
        
        total_distance += result['Leg Distance (meters)']
        total_time_seconds += result['Time Taken (seconds)']
        total_time_hours += result['Time Taken (hours)']
        total_battery_drain += result['Expected Battery Drain (amp-hours)']

    output_str += (
        "Round Trip Details:\n"
        f"Total Distance (meters): {total_distance}\n"
        f"Total Time (seconds): {total_time_seconds}\n"
        f"Total Time (hours): {total_time_hours}\n"
        f"Total Expected Battery Drain (amp-hours): {total_battery_drain}\n"
    )
    
    print(output_str)

    with open(os.path.join(directory, 'leg_analysis.txt'), 'w') as file:
        file.write(output_str)

    return total_distance, total_time_seconds, total_time_hours, total_battery_drain

# =========================
# 3 - RTOFS Data Processing
# =========================

### CLASS:
class RTOFS():

    ### FUNCTION:
    def __init__(self) -> None:
        
        """
        Initialize the RTOFS instance and fetch initial RTOFS data.
        """
        
        self._data_orig = self.rtofs_load()
        self._data_orig = self._data_orig.set_coords(['u', 'v'])
        self.data = self._data_orig.copy()
        self.x = self.data.x.values
        self.y = self.data.y.values
        self.grid_lons = self.data.lon.values[0,:]
        self.grid_lats = self.data.lat.values[:,0]
    
    ### FUNCTION:
    def rtofs_load(self):
        
        """
        Fetch the RTOFS data from the given URL and set its coordinates.
        """
        
        rtofs_access = "https://tds.marine.rutgers.edu/thredds/dodsC/cool/rtofs/rtofs_us_east_scraped"
        # rtofs_access = "rtofs_us_east_scraped.nc"
        
        try:
            rtofs_dataset = xr.open_dataset(rtofs_access).set_coords(['lon', 'lat'])
            rtofs_dataset.attrs['model'] = 'RTOFS'
            rtofs_dataset= rtofs_dataset.isel(time=-1)
            return rtofs_dataset
        except Exception as e:
            print(f"Error fetching RTOFS data: {e}")
            return None
    
    ### FUNCTION:
    def rtofs_subset(self, config, waypoints, buffer=0.5):
        
        """
        Subset the RTOFS data based on the bounding box created by the given points.
        """
        
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
    
        self.longitudes = self.data.lon.values[0,:]
        self.latitudes = self.data.lat.values[:,0]

        self.data = self.data.where(self.data['depth'] <= config["max_depth"], drop=True)


### FUNCTION:
def rtofs_depth_avg(RTOFS):
    
    """
    Calculate the depth-averaged currents for the given depth configurations.
    """
    
    u_currents = RTOFS.data['u']
    v_currents = RTOFS.data['v']
    depths = RTOFS.data['depth']
    
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

    return u_avg, v_avg

### FUNCTION:
def rtofs_current_magnitude(u_avg, v_avg):
        
    """
    Calculate the current magnitude for the given RTOFS data.
    """
    
    magnitude = np.sqrt(u_avg ** 2 + v_avg ** 2)
    
    return magnitude

# =========================
# 4 - Mission Plotting
# =========================

### FUNCTION:
def GGS_plot_route(config, waypoints, RTOFS_class, u_avg, v_avg, magnitude):

    """
    Plot the glider's mission route along with the depth-averaged currents and the advected path.

    Args:
    - rtofs (RTOFS): The RTOFS instance.
    - waypoints (list): A list of coordinates representing the mission route.
    - config (dict): The GGS configuration dictionary.
    - advected_route (list): A list of coordinates representing the advected route.
    """

    longitudes = RTOFS_class.data.lon.values
    latitudes = RTOFS_class.data.lat.values

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    # Add shaded relief
    terrain_tiles = GoogleTiles(style='terrain')
    ax.add_image(terrain_tiles, 8)

    contour = ax.contourf(longitudes, latitudes, magnitude, cmap=cmocean.cm.speed, transform=ccrs.PlateCarree())
    ax.streamplot(longitudes, latitudes, u_avg, v_avg, color='black', transform=ccrs.PlateCarree())
    # #
    # lats, lons = zip(*waypoints)
    # ax.plot(lons, lats, 'w-', transform=ccrs.PlateCarree(), linewidth=2.5, zorder=1)
    # ax.plot(lons, lats, 'k', transform=ccrs.PlateCarree(), linewidth=1.0, linestyle='--', alpha=0.6, zorder=1)
    
    # start_coords = config["waypoints"][0]
    # end_coords = config["waypoints"][-1]
    # ax.scatter(*start_coords[::-1], color='green', s=100, transform=ccrs.PlateCarree(), zorder=2)
    # for waypoint in config["waypoints"][1:-1]:
    #     ax.scatter(*waypoint[::-1], color='blue', s=100, transform=ccrs.PlateCarree(), zorder=2)
    # ax.scatter(*end_coords[::-1], color='red', s=100, transform=ccrs.PlateCarree(), zorder=2)
    # #        
    ax.set_extent([longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, edgecolor='black', zorder=0)
    ax.add_feature(cfeature.OCEAN, zorder=0, edgecolor='k', facecolor='lightblue')
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor='k',
                                                facecolor=cfeature.COLORS['land']))
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'lakes', '10m',
                                                edgecolor='none',
                                                facecolor=cfeature.COLORS['water']))
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '10m',
                                                edgecolor=cfeature.COLORS['water'],
                                                facecolor='none'))
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    cbar_ax = fig.add_axes([box.x0 + box.width + 0.01, box.y0, 0.03, box.height])
    fig.colorbar(contour, cax=cbar_ax, orientation='vertical', label='Depth Averaged Current Magnitude (m/s)')

    ax.set_xticks(np.linspace(longitudes.min(), longitudes.max(), 5), crs=ccrs.PlateCarree())
    ax.set_yticks(np.linspace(latitudes.min(), latitudes.max(), 5), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.set_title(f'{config["glider_name"]} Mission Route Along Depth-Averaged Currents', pad=20)
    
    mission_directory = os.path.join(os.path.expanduser("~"), "Downloads", f"GGS_{config['glider_name']}")
    fig_filename = f"{config['glider_name']}_depth_avg_currents.png"
    fig_path = os.path.join(mission_directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')

# =========================
# 5 - Slocum Files
# =========================

### FUNCTION:
def slocum_file_list(config, waypoints, directory, calculate_heading, DD_to_DDMM):
    
    """
    Offer a list of potential Slocum files and asks the user if they'd like to generate each.

    Args:
    - config (dict): The GGS configuration dictionary.
    - waypoints (list): List of waypoints.
    - directory (str): The path to the directory where the file should be saved.
    - calculate_heading (func): A function to calculate the heading.
    - DD_to_DDMM (func): A function to convert from decimal degrees to DDMM format.
    """

    slocum_functions = [
        {
            "name": "goto_l10",
            "function": goto_l10,
            "description": "generate a goto_l10.ma file",
            "args": (config, waypoints, directory, calculate_heading, DD_to_DDMM)
        },
        # {
        #     "name": "name",
        #     "function": function,
        #     "description": "description",
        #     "args": (arg1, arg2, arg3, ...)
        # }
    ]

    for function in slocum_functions:
        while True:
            user_input = input(f"Would you like to {function['description']}? (yes/no): ").lower()
            check_abort(user_input)
            if user_input == 'yes':
                function['function'](*function['args'])
                break
            elif user_input == 'no':
                break
            else:
                print("Invalid input. Please enter 'yes', 'no', or 'EXIT'.")

# =========================
# X - Execution
# =========================

### RUN:
EXIT_KEYWORD = "EXIT"
def main():
    
    """
    GGS main function.
    """

    buffer = 0

    config, waypoints = GGS_config_static() # DEBUG
    # config, waypoints = GGS_config() # MAIN
    directory = GGS_config_output(config)

    analysis_results = route_analysis(config, waypoints)
    route_analysis_output(config, directory, analysis_results)
    
    RTOFS_class = RTOFS()
    RTOFS_class.rtofs_subset(config, waypoints, buffer)

    u_avg, v_avg = rtofs_depth_avg(RTOFS_class)
    magnitude = rtofs_current_magnitude(u_avg, v_avg)

    GGS_plot_route(config, waypoints, RTOFS_class, u_avg, v_avg, magnitude)

    # slocum_file_list(config, waypoints, directory, calculate_heading, DD_to_DDMM)

if __name__ == "__main__":
    main()

# =========================
# ///// END OF SCRIPT \\\\\
# =========================