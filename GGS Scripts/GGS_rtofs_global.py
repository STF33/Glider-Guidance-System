# =========================
# 0 - Imports
# =========================

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import cmocean
import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import requests
import struct
from struct import unpack
import tarfile
from tqdm import tqdm
import xarray as xr

from sup_functions import check_abort, check_float, check_coordinate
from sup_functions import calculate_distance, calculate_heading, DD_to_DDMM
from sup_functions import read_ab, find_indx_lonlat, dx_dy, get_zz_zm
from sup_functions import dist_sphcrd

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
        "waypoints": [(21.11312876467859, -86.52437115832592),
                      (21.786507032637843, -86.61634071404619),
                      (22.031112029601204, -87.59121800468087)
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
BASE_URL = "https://noaa-nws-rtofs-pds.s3.amazonaws.com/rtofs.parallel.v2.3/"
class RTOFS_Global:

    """
    Represents the RTOFS (Real Time Ocean Forecast System) data and operations.
    """

    def __init__(self, directory):

        '''
        Initialize the RTOFS instance.

        Args:
        - directory (str): The directory where the RTOFS files should be downloaded.

        Returns:
        - RTOFS: The RTOFS instance.
        '''
        
        self.target_directory = directory

    @staticmethod
    def get_RTOFS_url():

        '''
        Get the URLs of the most recent RTOFS files.

        Returns:
        - list: A list of URLs for the most recent RTOFS files.
        '''

        desired_files = ['rtofs_glo.t00z.f12.archv.a.tgz', 'rtofs_glo.t00z.f12.archv.b']
        today = datetime.date.today()
        
        urls_today = [BASE_URL + "rtofs." + today.strftime('%Y%m%d') + "/" + file for file in desired_files]
        if all(requests.head(url).status_code == 200 for url in urls_today):
            return urls_today

        yesterday = today - datetime.timedelta(days=1)
        urls_yesterday = [BASE_URL + "rtofs." + yesterday.strftime('%Y%m%d') + "/" + file for file in desired_files]
        if all(requests.head(url).status_code == 200 for url in urls_yesterday):
            return urls_yesterday

        return None

    def download_RTOFS_files(self, file_urls):
        
        '''
        Download the specified files to the target directory.

        Args:
        - file_urls (list): A list of URLs for the files to be downloaded.

        Returns:
        - list: A list of the downloaded (and extracted) file paths.
        '''

        for file_url in file_urls:
            file_name = file_url.split('/')[-1]
            destination_path = os.path.join(self.target_directory, file_name)

            if self.check_RTOFS_files(file_name):
                print(f"File {file_name} or its extracted content already exists. Skipping download.")
                continue

            response = requests.get(file_url, stream=True)
            total_size = int(response.headers.get('content-length', 0))
            block_size = 1024
            t = tqdm(total=total_size, unit='iB', unit_scale=True)
            with open(destination_path, 'wb') as fd:
                for chunk in response.iter_content(block_size):
                    t.update(len(chunk))
                    fd.write(chunk)
            t.close()

            if destination_path.endswith('.tgz'):
                self.extract_RTOFS_files(destination_path)
                destination_path = destination_path.replace('.tgz', '')

    def check_RTOFS_files(self, file_name):
        
        """
        Check if the specified file or its extracted content already exists.
        
        Args:
        - file_name (str): Name of the file to check.

        Returns:
        - bool: True if file or its extracted content exists, False otherwise.
        """

        if os.path.exists(os.path.join(self.target_directory, file_name)):
            return True
        
        if file_name.endswith('.tgz'):
            extracted_file_name = file_name.replace('.tgz', '')
            if os.path.exists(os.path.join(self.target_directory, extracted_file_name)):
                return True
        
        return False

    def extract_RTOFS_files(self, file_path):
        
        """
        Extract the contents of a .tgz file to its current directory.
        
        Args:
        - file_path (str): Path to the .tgz file to be extracted.
        """
        
        try:
            with tarfile.open(file_path, 'r:gz') as tar_ref:
                tar_ref.extractall(self.target_directory)
            os.remove(file_path)
            print(f"Extracted and removed {file_path}")
        except Exception as e:
            print(f"Failed to extract {file_path}. Error: {e}")

    def fetch_RTOFS_files(self):
        recent_files = self.get_RTOFS_url()
        if recent_files:
            return self.download_RTOFS_files(recent_files)
        else:
            print("Could not find the recent RTOFS files.")
            return []

### READ AND PROCESS AB FILES ###





### READ AND PROCESS AB FILES ###

# =========================
# X - Execution
# =========================

### RUN:
EXIT_KEYWORD = "EXIT"
def main():
    
    """
    GGS main function.
    """

    config, waypoints = GGS_config_static() # DEBUG
    # config, waypoints = GGS_config() # MAIN
    directory = GGS_config_output(config)

    analysis_results = route_analysis(config, waypoints)
    route_analysis_output(config, directory, analysis_results)
    
    rtofs = RTOFS()
    rtofs.rtofs_data_subset(waypoints)

    GGS_plot_route(config, waypoints, rtofs)

    # slocum_file_list(config, waypoints, directory, calculate_heading, DD_to_DDMM)

if __name__ == "__main__":
    main()

# =========================
# ///// END OF SCRIPT \\\\\
# =========================