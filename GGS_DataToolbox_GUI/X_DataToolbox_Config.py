"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import datetime as dt
from dateutil import parser
import json
import os

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
            
            glider_config = config['GLIDER']
            glider_config['glider_unit'] = glider_config.get('glider_unit', 'unknown')
            glider_config['glider_version'] = glider_config.get('glider_version', 'unknown')
            glider_config['glider_type'] = glider_config.get('glider_type', 'unknown')
           
            sensors_config = config['SENSORS']
            sensors_config['sensor_list'] = sensors_config.get('sensor_list', [])

            decompression_config = config['DECOMPRESSION']
            decompression_config['run_decompression'] = decompression_config.get('run_decompression', False)
            
            conversion_config = config['CONVERSION']
            conversion_config['run_conversion'] = conversion_config.get('run_conversion', False)
            
            data_config = config['DATA']
            data_config['run_dataframe'] = data_config.get('run_dataframe', False)

            product_config = config['PRODUCTS']
            product_config['run_plot'] = product_config.get('run_plot', False)
            product_config['run_excel'] = product_config.get('run_excel', False)
            product_config['run_data_sorter'] = product_config.get('run_data_sorter', False)
            
            advanced_config = config['ADVANCED']
            advanced_config['run_data_cleanup'] = advanced_config.get('run_data_cleanup', False)
    
    except Exception as e:
        print(f"Error during config import: {e}")
        return None

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

    try:
        glider_unit = config['GLIDER']['glider_unit']
    except:
        glider_unit = "UnknownGlider"
        print(f"Error using provided 'glider_unit'. Defaulting glider unit to 'UnknownGlider'.")
    try:
        glider_version = config['GLIDER']['glider_version']
    except:
        glider_version = "UnknownVersion"
        print(f"Error using provided 'glider_version'. Defaulting glider version to 'UnknownVersion'.")
    try:
        glider_type = config['GLIDER']['glider_type']
    except:
        glider_type = "UnknownType"
        print(f"Error using provided 'glider_type'. Defaulting glider type to 'UnknownType'.")

    if path == "default":
        root_directory = os.path.join(os.path.expanduser("~"), "Downloads", f"GGS_{glider_unit}-{glider_version}-{glider_type}")
    else:
        root_directory = os.path.join(path, f"GGS_{glider_unit}-{glider_version}-{glider_type}")
    
    os.makedirs(root_directory, exist_ok=True)

    output_str = "Configuration:\n"

    for section, settings in config.items():
        output_str += f"\n--> {section}\n"
        for key, value in settings.items():
            formatted_value = str(value) if not isinstance(value, dt.datetime) else value.strftime('%Y-%m-%dT%H:%M:%S%z')
            output_str += f"{key}: {formatted_value}\n"

    print(output_str)
    print("### END OF CONFIGURATION ###")
    print("\n")
    
    return root_directory