"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

# =========================

### FUNCTION:
def pull_sensor_rate(sensor, master_dict):

    '''
    Retrieve the data rate for a specified sensor from the master dictionary.
    
    Args:
    - sensor (str): The name of the sensor.
    - master_dict (dict): The master dictionary containing sensor data.
    
    Returns:
    - rate (float): The data rate of the specified sensor.
    '''

    key = list(master_dict.keys())[0]
    rate = master_dict[key]['data'][sensor]['rate']
    return rate

### FUNCTION:
def pull_sensor_units(sensor, master_dict):

    '''
    Retrieve the units for a specified sensor from the master dictionary.
    
    Args:
    - sensor (str): The name of the sensor.
    - master_dict (dict): The master dictionary containing sensor data.
    
    Returns:
    - units (str): The units of the specified sensor.
    '''

    key = list(master_dict.keys())[0]
    units = master_dict[key]['data'][sensor]['units']
    return units

### FUNCTION:
def pull_sensor_data(sensor, master_dict):

    '''
    Retrieve the data values for a specified sensor from the master dictionary.
    
    Args:
    - sensor (str): The name of the sensor.
    - master_dict (dict): The master dictionary containing sensor data.
    
    Returns:
    - data (list): A list of data values for the specified sensor.
    '''
    
    data = []
    for key in master_dict:
        if sensor in master_dict[key]['data']:
            active_data = master_dict[key]['data'][sensor]
            for value in active_data['values']:
                data.append(float(active_data['values'][value]))
    return data

### FUNCTION:
def pull_sensor_list_data(sensor_list, chosen_dataframe, master_dict):

    '''
    Retrieve data for a list of sensors and append it to the chosen DataFrame.
    
    Args:
    - sensor_list (list): A list of sensor names.
    - chosen_dataframe (pd.DataFrame): The DataFrame to append the sensor data to.
    - master_dict (dict): The master dictionary containing sensor data.
    
    Returns:
    - chosen_dataframe (pd.DataFrame): The updated DataFrame with appended sensor data.
    '''
    
    subframe = {}
    max_length = 0
    for sensor in sensor_list:
        try:
            sensor_data = pull_sensor_data(sensor, master_dict)
            subframe[sensor] = sensor_data
            if len(sensor_data) > max_length:
                max_length = len(sensor_data)
        except KeyError:
            print(f"Warning: Sensor '{sensor}' not found in the data. Skipping...")
            continue
    
    # Pad shorter lists with NaN values
    for sensor in subframe:
        if len(subframe[sensor]) < max_length:
            subframe[sensor] += [np.nan] * (max_length - len(subframe[sensor]))
    
    subframe = pd.DataFrame(subframe)
    chosen_dataframe = pd.concat([chosen_dataframe, subframe], ignore_index=True)
    return chosen_dataframe
