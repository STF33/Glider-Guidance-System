# =========================
# X - IMPORTS
# =========================

import os
import pandas as pd
from X_functions import check_abort, check_float, check_coordinate

# =========================
# STATIC CONFIGURATION
# =========================

### FUNCTION:
def GGS_config_static():
    
    '''
    Configure mission from a hardcoded configuration.
    Intended for testing only.
    
    Args:
    - None

    Returns:
    - config (dict): Glider Guidance System mission configuration.
    - waypoints (list): Glider Guidance System mission waypoints.
    '''
    
    config = {
        "glider_name": "Gulf of Mexico",
        "max_depth": 100.0,
        "avg_velocity": 0.5,
        "battery_capacity": 1000,
        "battery_drain": 5,
        "satisfying_radius": 1000,
        "waypoints": [(18.00510558149081, -97.75507520595109),
                      (30.68658969837126, -80.17002034924428)
                    ] # CURRENT: Gulf of Mexico
    }

    waypoints = config["waypoints"]

    return config, waypoints

# =========================
# CONFIGURATION PATHS
# =========================

### FUNCTION:
def GGS_config_new():
    
    '''
    Configure a new mission.

    Args:
    - None

    Returns:
    - config (dict): Glider Guidance System mission configuration.
    - waypoints (list): Glider Guidance System mission waypoints.
    '''

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
def GGS_config_import():
    
    '''
    Configure mission from an imported file.
    
    Args:
    - None

    Returns:
    - config (dict): Glider Guidance System mission configuration.
    - waypoints (list): Glider Guidance System mission waypoints.
    '''
    
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

# =========================
# CONFIGURATION PROCESSING
# =========================

### FUNCTION:
EXIT_KEYWORD = "EXIT"
def GGS_config():
    
    '''
    Initialize the Glider Guidance System configuration protocol.
    
    Args:
    - None

    Returns:
    - config (dict): Glider Guidance System mission configuration.
    - waypoints (list): Glider Guidance System mission waypoints.
    '''
    
    print(f"Glider Guidance System (GGS) Configuration Setup\n(Note: Type '{EXIT_KEYWORD}' at any time to abort the configuration.)\n")
    choice = check_abort(input("Do you want to import an existing GGS config or create a new one? (Enter 'import' or 'new'): ").lower())
    
    if choice == "import":
        return GGS_config_import()
    elif choice == "new":
        return GGS_config_new()
    else:
        print("Invalid choice. Please enter 'import' or 'new'.")
        return GGS_config()

### FUNCTION:
def GGS_config_output(config):
    
    '''
    Output and save the configured Glider Guidance System mission.
    
    Args:
    - config (dict): Glider Guidance System mission configuration.

    Returns:
    - directory (str): Glider Guidance System mission directory.
    '''

    output_str = "\nGlider Guidance System (GGS) Configuration:\n"

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
