# =========================
# X - IMPORTS
# =========================

import datetime as dt
from datetime import timezone
import os
import pandas as pd
from X_functions import check_abort, check_float, check_coordinate, get_date_list

# =========================
# STATIC CONFIGURATION
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
        "glider_name": "Test",
        "target_date": target_date,
        "date_list": date_list,
        "max_depth": 1000,
        "extent": [(20.604, -87.022), (21.321, -86.140)],  # Test
        # "extent": [(17.0, -98.0), (30.5, -80.0)],  # Yucatan
        # "extent": [(15.25, -90.5), (25.25, -80.0)],  # GoM
        "GPS_coords": [(20.375, -86.541), (21.025, -86.349), (21.506, -86.528)]
        }

    return config

# =========================
# CONFIGURATION PATHS
# =========================

### FUNCTION:
# def GGS_config_new():
    
#     '''
#     Configure a new mission.

#     Args:
#     - None

#     Returns:
#     - config (dict): Glider Guidance System mission configuration.
#     - GPS_coords (list): Glider Guidance System mission GPS_coords.
#     '''

#     config = {}

#     prompt_glider_name = "Enter the glider's name: "
#     while True:
#         glider_name = input(prompt_glider_name)
#         check_abort(glider_name)
#         if glider_name:
#             config["glider_name"] = glider_name
#             break
#         prompt_glider_name = "[INPUT ERROR] " + prompt_glider_name

#     config["max_depth"] = check_float("Enter the maximum mission depth (meters): ")

#     prompt_num_GPS_coords = "Enter the total number of GPS_coords (1-100): "
#     GPS_coords = []
#     while True:
#         try:
#             num_GPS_coords = int(input(prompt_num_GPS_coords))
#             check_abort(str(num_GPS_coords))
#             if 1 <= num_GPS_coords <= 100:
#                 break
#         except ValueError:
#             pass
#         prompt_num_GPS_coords = "[INPUT ERROR] " + prompt_num_GPS_coords

#     for i in range(num_GPS_coords):
#         GPS_coord_prompt = f"Enter GPS_coord {i} coordinates as 'latitude, longitude': "
#         while True:
#             GPS_coord = check_coordinate(input(GPS_coord_prompt))
#             if GPS_coord:
#                 GPS_coords.append(GPS_coord)
#                 break
#             GPS_coord_prompt = "[INPUT ERROR] " + GPS_coord_prompt

#     config["GPS_coords"] = GPS_coords

#     GGS_config_output(config)
#     return config, GPS_coords

### FUNCTION:
# def GGS_config_import():
    
#     '''
#     Configure mission from an imported file.
    
#     Args:
#     - None

#     Returns:
#     - config (dict): Glider Guidance System mission configuration.
#     - GPS_coords (list): Glider Guidance System mission GPS_coords.
#     '''
    
#     while True:
#         file_name = input("Please enter the path to your GGS config file: ")
#         check_abort(file_name)
#         try:
#             config = pd.read_pickle(file_name)
#             GPS_coords = config["GPS_coords"]
#             GGS_config_output(config)
#             return config, GPS_coords
#         except Exception as e:
#             print(f"Error reading the file: {e}")
#             retry = input("Invalid file. Would you like to try another file? (yes or no): ").lower()
#             check_abort(retry)
#             if retry != "yes":
#                 print("Import aborted. Please restart the configuration process.")
#                 return None

# =========================
# CONFIGURATION PROCESSING
# =========================

### FUNCTION:
# EXIT_KEYWORD = "EXIT"
# def GGS_config():
    
#     '''
#     Initialize the Glider Guidance System configuration protocol.
    
#     Args:
#     - None

#     Returns:
#     - config (dict): Glider Guidance System mission configuration.
#     - GPS_coords (list): Glider Guidance System mission GPS_coords.
#     '''
    
#     print(f"Glider Guidance System (GGS) Configuration Setup\n(Note: Type '{EXIT_KEYWORD}' at any time to abort the configuration.)\n")
#     choice = check_abort(input("Do you want to import an existing GGS config or create a new one? (Enter 'import' or 'new'): ").lower())
    
#     if choice == "import":
#         return GGS_config_import()
#     elif choice == "new":
#         return GGS_config_new()
#     else:
#         print("Invalid choice. Please enter 'import' or 'new'.")
#         return GGS_config()

### FUNCTION:
def GGS_config_output(config):
    
    '''
    Output and save the configured Glider Guidance System mission.
    
    Args:
    - config (dict): Glider Guidance System mission configuration.

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
    directory = os.path.join(os.path.expanduser("~"), "Downloads", f"GGS_{formatted_date}_{config['glider_name']}")
    os.makedirs(directory, exist_ok=True)

    config_pickle = os.path.join(directory, f"GGS_{formatted_date}_{config['glider_name']}_config.pkl")
    pd.to_pickle(config, config_pickle)

    config_text = os.path.join(directory, f"GGS_{formatted_date}_{config['glider_name']}_config.txt")
    with open(config_text, 'w') as file:
        file.write(output_str)

    print(output_str)

    return directory
