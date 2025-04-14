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
# CONFIG HANDLING
# =========================

def config_import(config_name):

    ''' 
    Import a Glider Guidance System mission configuration from a JSON file.
    
    Arguments:
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
            data_config['run_data_filter'] = data_config.get('run_data_filter', False)
            data_config['run_data_sorter'] = data_config.get('run_data_sorter', False)
            data_config['run_logfile_search'] = data_config.get('run_logfile_search', False)
            
            product_config = config['PRODUCTS']
            product_config['run_plot'] = product_config.get('run_plot', False)
            product_config['run_excel'] = product_config.get('run_excel', False)
            product_config['run_energy_evaluation'] = product_config.get('run_energy_evaluation', False)
            
            advanced_config = config['ADVANCED']
            advanced_config['cleanup_run'] = advanced_config.get('cleanup_run', False)
    
    except Exception as e:
        print(f"Error during config import: {e}")
        return None
    
    print("Configuration import success!")
    
    return config

def config_process(config, path="default"):

    ''' 
    Process the imported configuration and create an output directory based on glider parameters.
    
    Arguments:
    - config (dict): Glider Guidance System mission configuration.
    - path (str): Directory path to save output to. Use "default" to save to the user's Downloads directory.
      
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

# =========================
# PROGRAM CLEANUP
# =========================

def cleanup_define_directories():

    ''' 
    Define directories that contain run-specific data files.
    
    Arguments:
    - None
      
    Returns:
    - directories_to_clean (list): A list of directory paths.
    '''
    
    current_directory = os.path.dirname(__file__)
    directories_to_clean = [
        os.path.join(current_directory, 'DBD_Files'),
        os.path.join(current_directory, 'DBD_Files', 'ProcessedAscii'),
        os.path.join(current_directory, 'DBD_Files', 'Decompressed'),
        os.path.join(current_directory, 'DBD_Files', 'Logfiles')
    ]
    
    return directories_to_clean

def cleanup_run(directories_to_clean):

    ''' 
    Delete data files in the specified directories.
    
    Arguments:
    - directories_to_clean (list): A list of directory paths.
      
    Returns:
    - None
    '''

    print(f"\n### RUNNING: DATA LIBRARY CLEANUP ###\n")
    
    for directory in directories_to_clean:
        if os.path.exists(directory):
            for item in os.listdir(directory):
                item_path = os.path.join(directory, item)
                if os.path.isfile(item_path) and item != '.gitignore':
                    os.remove(item_path)
                    print(f"Deleted file: {item_path}")
                else:
                    print(f"Skipped item: {item_path}")
        else:
            print(f"Directory does not exist: {directory}")
