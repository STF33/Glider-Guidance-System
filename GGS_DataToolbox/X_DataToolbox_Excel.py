"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import os
import pandas as pd
from openpyxl import load_workbook

from X_DataToolbox_SensorTools import *

# =========================

### FUNCTION:
def read_ascii(directory, file_name):

    '''
    Read and parse an ASCII file into a dictionary format.
    
    Args:
    - directory (str): The directory containing the ASCII file.
    - file_name (str): The name of the ASCII file to read.
    
    Returns:
    - active_dictionary (dict): The parsed data from the ASCII file.
    - file_identifier (str): The identifier of the file (either 'filename' or 'full_filename').
    '''

    entry = 0
    active_dictionary = {}
    active_file = open(directory + '/' + file_name, 'r')
    file_content = active_file.readlines()
    line_count = 0

    for line in file_content:
        line = line.replace('\n', '')
        if line_count < 13:
            line = line.replace(' ', '')
            line = line.split(':')
            active_dictionary[str(line[0])] = line[1]
        if line_count == 14:
            active_dictionary['data'] = {}
            all_sensors = []
            line = line.split(' ')
            for sensor in line:
                active_dictionary['data'][str(sensor)] = {'units': '', 'rate': '', 'values': {}}
                all_sensors.append(sensor)
        if line_count == 15:
            line = line.split(' ')
            sens_index = 0
            for unit in line:
                sensor = all_sensors[sens_index]
                active_dictionary['data'][str(sensor)]['units'] = unit
                sens_index += 1
        if line_count == 16:
            line = line.split(' ')
            entry = 0
            sub_entry = 0
            for val in line:
                sensor = all_sensors[sub_entry]
                active_dictionary['data'][str(sensor)]['rate'] = val
                sub_entry += 1
        if line_count > 16:
            line = line.split(' ')
            index = 0
            for val in line:
                sensor = all_sensors[index]
                active_dictionary['data'][str(sensor)]['values'][str(entry)] = val
                index += 1
            entry += 1
        line_count += 1
    active_dictionary['data'].pop('', None)

    if 'filename' in active_dictionary.keys():
        return active_dictionary, active_dictionary['filename']
    else:
        return active_dictionary, active_dictionary['full_filename']

### FUNCTION:
def dbd_directory_to_dict(directory, files_read, num_read=0):

    '''
    Convert a directory of DBD files to a dictionary format.
    
    Args:
    - directory (str): The directory containing the DBD files.
    - files_read (list): A list of files that have already been read.
    - num_read (int): The number of files read so far.
    
    Returns:
    - master_dict (dict): The master dictionary containing data from the DBD files.
    - files_read (list): The updated list of files that have been read.
    - num_read (int): The updated number of files read.
    - first_file (str): The name of the first file processed in this batch.
    - last_file (str): The name of the last file processed in this batch.
    '''

    master_dict = {}
    chunk_count = 0
    for file in os.listdir(directory):
        if chunk_count >= 30:
            break
        if file not in files_read:
            if chunk_count == 0:
                first_file = file
            if file.endswith("dbd.asc") or file.endswith("sbd.asc") or file.endswith("DBD.asc") or\
                    file.endswith("SBD.asc") or file.endswith("ebd.asc"):
                sub_dict, sub_name = read_ascii(directory, file)
                master_dict[sub_name] = sub_dict
                files_read.append(file)
                last_file = file
                num_read += 1
                chunk_count += 1
        else:
            continue
    print("Files {} to {} Processed to Dict, Following Values pertain to that Chunck". format(first_file, last_file))
    return master_dict, files_read, num_read, first_file, last_file

### FUNCTION:
def initialize_data_frame(sensor_list):

    '''
    Initialize a DataFrame with specified sensor columns.
    
    Args:
    - sensor_list (list): A list of sensor names to include as columns in the DataFrame.
    
    Returns:
    - frame (pd.DataFrame): The initialized DataFrame with sensor columns.
    '''

    frame = {}

    for sensor in sensor_list:
        frame[sensor] = []
    frame = pd.DataFrame(frame)

    return frame

### FUNCTION:
def run_dataframe(input_directory, sensor_list):
    
    '''
    Process sensor data from a directory and tabulate it into a DataFrame.
    
    Args:
    - input_directory (str): The directory containing the input files.
    - sensor_list (list): A list of sensor names to include in the DataFrame.

    Returns:
    - dataframe (pd.DataFrame): The tabulated DataFrame containing sensor data.
    '''
    
    print(f"\n### RUNNING: DATAFRAME ###\n")

    files_processed = []
    tot_num_files = len(os.listdir(input_directory))
    num_read = 0
    dataframe = initialize_data_frame(sensor_list)

    while num_read < tot_num_files:
        master_dict, files_processed, num_read, first_file, last_file = dbd_directory_to_dict(input_directory, files_processed, num_read)
        dataframe = pull_sensor_list_data(sensor_list, dataframe, master_dict)
        print('{} Files Processed'.format(num_read))
        master_dict.clear()

    return dataframe

### FUNCTION:
def run_excel(root_directory, glider_info, data_frame):
    
    '''
    Save the DataFrame to an Excel file.
    
    Args:
    - data_frame (pd.DataFrame): The DataFrame containing sensor data.
    - root_directory (str): The root directory from the configuration file.
    - glider_info (tuple): A tuple containing glider unit, version, and type information.

    Returns:
    - None
    '''
    
    print(f"\n### RUNNING: EXCEL ###\n")

    glider_unit, glider_version, glider_type = glider_info
    output_file = os.path.join(root_directory, f"{glider_unit}-{glider_version}-{glider_type}_DataOutput.xlsx")
    writer = pd.ExcelWriter(output_file, engine='openpyxl')
    data_frame.to_excel(writer, index=False)
    writer.close()

    workbook = load_workbook(output_file)
    worksheet = workbook.active
    for col in worksheet.columns:
        worksheet.column_dimensions[col[0].column_letter].width = 20
    workbook.save(output_file)

    print(f"All data saved to {output_file}")