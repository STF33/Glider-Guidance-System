# =========================
# IMPORTS
# =========================

import datetime as dt
from datetime import timezone
import os
import pandas as pd

# =========================

### FUNCTION:
def GGS_config_static(date=dt.datetime.now(timezone.utc)):
    
    '''
    Configure mission from a hardcoded configuration.
    
    Args:
    - date (datetime): Date of mission start.
        - default: dt.datetime.now(timezone.utc)

    Returns:
    - config (dict): Glider Guidance System mission configuration.
    '''

    execution_date = date

    current_directory = os.path.dirname(__file__)
    bathymetry_path = os.path.join(current_directory, "data", "bathymetry", "GEBCO_2023_sub_ice_topo.nc")
    eez_path = os.path.join(current_directory, "data", "eez", "eez_boundaries_v12.shp")
    
    config = {
        "glider_name": "SENTINEL2",
        "execution_date": execution_date,
        "max_depth": 1000,
        
        # UGOS (East)
        # "extent": [(15, -90), (30, -78)],
        # "GPS_coords": [(0, 0), (0, 0)],

        # Sentinel 1
        # "extent": [(10, -78), (44, -10)],
        # "GPS_coords": [(41.675, -70.522), (15.067, -23.650)],

        # Sentinel 2
        "extent": [(20, -30), (-40, 55)],
        "GPS_coords": [(15.067, -23.650), (-33.907, 18.564)],

        # GLOBAL
        # "extent": [(-80, -180), (90, 180)],
        # "GPS_coords": [(0, 0), (0, 0)],

        "bathymetry_path": bathymetry_path,
        "eez_path": eez_path
        }

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
            output_str += "\nGPS_coords:\n"
            for i, coord in enumerate(value, 1):
                output_str += "  GPS_coord {}: {}\n".format(i, coord)
        elif key == "extent":
            output_str += "\nExtent:\n"
            for i, coord in enumerate(value, 1):
                output_str += "  Boundary {}: {}\n".format(i, coord)
        elif key in ["glider_name", "execution_date"]:
            output_str += "{}: {}\n\n".format(key.capitalize().replace('_', ' '), value)
        elif key == "max_depth":
            output_str += "Max Depth: {}\n".format(value)

    formatted_date = config['execution_date'].strftime("%Y%m%d")

    if path == "default":
        root_directory = os.path.join(os.path.expanduser("~"), "Downloads", f"GGS_{config['glider_name']}")
    else:
        root_directory = os.path.join(path, f"GGS_{config['glider_name']}")

    os.makedirs(root_directory, exist_ok=True)

    config_pickle = os.path.join(root_directory, f"GGS_{config['glider_name']}_config_{formatted_date}.pkl")
    pd.to_pickle(config, config_pickle)

    config_text = os.path.join(root_directory, f"GGS_{config['glider_name']}_config_{formatted_date}.txt")
    with open(config_text, 'w') as file:
        file.write(output_str)

    print(output_str)

    return root_directory
