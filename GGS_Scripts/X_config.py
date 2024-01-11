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

    config = {
        "glider_name": "GoM",
        "target_date": target_date,
        "date_list": date_list,
        "max_depth": 10,
        # "extent": [(20.604, -87.022), (21.321, -86.140)],  # TEST
        # "extent": [(15.25, -90.5), (25.25, -80.0)],  # Yucatan
        "extent": [(17.0, -98.0), (30.5, -80.0)],  # GoM
        "GPS_coords": [(20.375, -86.541), (21.025, -86.349), (21.506, -86.528)]
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
