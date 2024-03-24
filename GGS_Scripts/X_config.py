# =========================
# IMPORTS
# =========================

import datetime as dt
from datetime import datetime, timezone
import ast
import os
import pandas as pd
from X_functions import acquire_gliders

# =========================

### FUNCTION:
def GGS_config_import(config_name):

    '''
    Import a Glider Guidance System mission configuration from a text file.

    Args:
    - config_name (str): Name of the configuration file to import.

    Returns:
    - config (dict): Glider Guidance System mission configuration.
    '''

    current_directory = os.path.dirname(__file__)
    config_path = os.path.join(current_directory, "config", f"{config_name}.txt")

    with open(config_path, 'r') as file:
        lines = file.readlines()

    config = {}
    for line in lines:
        key, value = line.strip().split(": ", 1)
        if key == 'extent' or key == 'GPS_coords':
            config[key] = ast.literal_eval(value)
        elif key == 'max_depth':
            config[key] = int(value)
        elif key == 'target_date':
            if value.strip('"') == "None":
                config[key] = dt.datetime.now(timezone.utc)
            else:
                config[key] = dt.datetime.fromisoformat(value.strip('"'))
        else:
            config[key] = value.strip('"')

    if config.get('glider_id', "None") != "None":
        glider_df = acquire_gliders(extent=None, target_date=config['target_date'], date_delta=dt.timedelta(days=1), requested_variables=["time", "longitude", "latitude", "profile_id", "depth"], print_vars=False, target=config['glider_id'], request_timeout=5, enable_parallel=False)
        if not glider_df.empty:
            last_lon = glider_df['longitude'].iloc[-1]
            last_lat = glider_df['latitude'].iloc[-1]
            buffer_size = 5
            config['extent'] = [(last_lat - buffer_size, last_lon - buffer_size), (last_lat + buffer_size, last_lon + buffer_size)]

    config["bathymetry_path"] = os.path.join(current_directory, "data", "bathymetry", "GEBCO_2023_sub_ice_topo.nc")
    config["eez_path"] = os.path.join(current_directory, "data", "eez", "eez_boundaries_v12.shp")

    return config

### FUNCTION:
def GGS_config_output(config, path="default"):
    
    '''
    Output and save the configured Glider Guidance System mission.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - path (str): Path for saving the mission configuration.
        - default: "default"

    Returns:
    - root_directory (str): Root directory for saving the mission configuration.
    '''

    output_str = "\n\n### Glider Guidance System (GGS) Configuration ###\n\n"

    for key, value in config.items():
        if key == "GPS_coords":
            output_str += "GPS_coords:\n"
            for coord in value:
                output_str += f"  - {coord}\n"
            output_str += "\n"
        elif key == "extent":
            output_str += "Extent:\n  Min Latitude, Min Longitude: {}\n  Max Latitude, Max Longitude: {}\n\n".format(value[0], value[1])
        else:
            formatted_value = value if not isinstance(value, datetime) else value.strftime('%Y-%m-%dT%H:%M:%S%z')
            output_str += f"{key.replace('_', ' ').capitalize()}: {formatted_value}\n"
            if key == "glider_id":
                output_str += "\n"

    mission_name = config.get('mission_name', 'UnknownMission')
    if path == "default":
        root_directory = os.path.join(os.path.expanduser("~"), "Downloads", f"GGS_{mission_name}")
    else:
        root_directory = os.path.join(path, f"GGS_{mission_name}")

    os.makedirs(root_directory, exist_ok=True)

    print(output_str)

    return root_directory
