# =========================
# IMPORTS
# =========================

import datetime as dt
from dateutil import parser
import json
import os
from X_functions import acquire_gliders

# =========================

### FUNCTION:
def GGS_config_import(config_name):
    
    '''
    Import a Glider Guidance System mission configuration from a JSON file.
    
    Args:
    - config_name (str): Name of the configuration file to import.
    
    Returns:
    - config (dict): Glider Guidance System mission configuration.
    '''

    print(f"\n### IMPORTING GGS CONFIGURATION: {config_name} ###\n")

    current_directory = os.path.dirname(__file__)
    config_path = os.path.join(current_directory, "config", f"{config_name}.json")

    try:
        with open(config_path, 'r') as file:
            config = json.load(file)
            
            mission_config = config['MISSION']
            mission_config['target_date'] = dt.datetime.now(dt.timezone.utc) if mission_config['target_date'] is None else dt.datetime.strptime(mission_config['target_date'], '%Y-%m-%d %H:%M:%S').replace(tzinfo=dt.timezone.utc)
            mission_config['GPS_coords'] = None if mission_config['GPS_coords'] is None else mission_config['GPS_coords']
            mission_config['extent'] = tuple(map(tuple, mission_config['extent']))
            
            plot_config = config['PLOT']
            plot_config['manual_extent'] = None if plot_config['manual_extent'] is None else tuple(map(tuple, plot_config['manual_extent']))

            data_config = config['DATA']
            data_config['bathymetry_path'] = os.path.join(current_directory, data_config['bathymetry_path'])
            data_config['eez_path'] = os.path.join(current_directory, data_config['eez_path'])
    
    except Exception as e:
        print(f"Error during config import: {e}")
        return None

    glider_id = config['MISSION'].get('glider_id')
    if glider_id and glider_id != "0":
        print(f"Locating glider: {glider_id}")
        try:
            glider_df = acquire_gliders(
                extent=None,
                target_date=dt.datetime.now(dt.timezone.utc),
                date_delta=dt.timedelta(days=1),
                requested_variables=["time", "longitude", "latitude", "profile_id", "depth"],
                print_vars=False,
                target=glider_id,
                request_timeout=5,
                enable_parallel=False
            )
            if not glider_df.empty:
                last_lon = glider_df['longitude'].iloc[-1]
                last_lat = glider_df['latitude'].iloc[-1]
                buffer = 5
                config['MISSION']['extent'] = [[last_lat - buffer, last_lon - buffer], [last_lat + buffer, last_lon + buffer]]
        except Exception as e:
            print(f"Error updating extent based on glider ID {glider_id}: {e}")

    print("Configuration import success!")

    return config

### FUNCTION:
def GGS_config_process(config, path="default"):
    
    '''
    Output the configured Glider Guidance System configuration.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - path (str): Path for saving the mission configuration. 'default' saves to the user's Downloads directory.

    Returns:
    - root_directory (str): The directory where the configuration was saved.
    '''

    print("\n### PROCESSING GGS CONFIGURATION ###\n")

    mission_name = config['MISSION'].get('mission_name', 'UnknownMission')
    
    if path == "default":
        root_directory = os.path.join(os.path.expanduser("~"), "Downloads", f"GGS_{mission_name}")
    else:
        root_directory = os.path.join(path, f"GGS_{mission_name}")
    
    os.makedirs(root_directory, exist_ok=True)

    output_str = "Configuration:\n"

    for section, settings in config.items():
        output_str += f"\n--> {section}\n"
        for key, value in settings.items():
            formatted_value = str(value) if not isinstance(value, dt.datetime) else value.strftime('%Y-%m-%dT%H:%M:%S%z')
            output_str += f"{key}: {formatted_value}\n"

    print(output_str)
    
    return root_directory
