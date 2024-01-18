# =========================
# X - IMPORTS
# =========================

import datetime as dt
from datetime import timezone
import os
import pandas as pd
from X_functions import get_date_list

# =========================
# GGS CONFIGURATION
# =========================

### FUNCTION:
def GGS_config_static(date=dt.datetime.now(timezone.utc)):
    
    '''
    Configure mission from a hardcoded configuration.
    Intended for testing only.
    
    Args:
    - date (datetime): Date of mission start.
        - default: dt.datetime.now(timezone.utc)

    Returns:
    - config (dict): Glider Guidance System mission configuration.
    '''
    
    target_date = date
    date_list = get_date_list(date)

    current_directory = os.path.dirname(__file__)
    bathymetry_path = os.path.join(current_directory, "data", "bathymetry", "GEBCO_2023_sub_ice_topo.nc")

    config = {
        "glider_name": "Yucatan",
        "target_date": target_date,
        "date_list": date_list,
        "max_depth": 1000,
        "extent": [(15.75, -89.25), (27.00, -80.00)],  # Yucatan
        # "extent": [(17.0, -98.0), (30.5, -80.0)],  # GoM
        "GPS_coords": [(16.645, -87.880), (16.952, -87.163), (18.196, -86.727), (19.562, -86.281), (21.621, -86.099)],
        "bathymetry_path": bathymetry_path
        }

    return config

### FUNCTION:
def GGS_config_output(config, path="default"):
    
    '''
    Output and save the configured Glider Guidance System mission.
    
    Args:
    - config (dict): Glider Guidance System mission configuration.
    - path (str): Path for saving the mission configuration. 
        - options: "default", "(custom path)"

    Returns:
    - directory (str): Glider Guidance System mission directory.
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
        elif key in ["glider_name", "target_date", "date_list"]:
            output_str += "{}: {}\n\n".format(key.capitalize().replace('_', ' '), value)
        elif key == "max_depth":
            output_str += "Max Depth: {}\n".format(value)

    formatted_date = config['target_date'].strftime("%Y%m%d")

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
