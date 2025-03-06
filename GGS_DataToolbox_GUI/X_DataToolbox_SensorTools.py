"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import os
import pandas as pd
import numpy as np

# =========================

### FUNCTION:
def read_ascii(directory, file_name):

    '''
    Read and parse an ASCII file into a dictionary format.
    
    Expected format:
      - The first N lines are metadata (where N is given by the metadata key "num_ascii_tags"; default is 13 if missing).
      - The next M lines are header lines (where M is given by the metadata key "num_label_lines"; default is 3).
          * The first header line should contain sensor names.
          * The second header line should contain sensor units.
          * The third header line should contain sensor rates.
      - The remaining lines are data rows with whitespace‐separated numeric values.
    
    Args:
      directory (str): The directory containing the ASCII file.
      file_name (str): The name of the ASCII file.
      
    Returns:
      result (dict): A dictionary with keys:
         "metadata": dict of metadata,
         "data": dict mapping each sensor name to a dict with keys "units", "rate", and "values" (a list of floats).
      file_identifier (str): The identifier for the file (from metadata if available, else the file name).
    '''

    filepath = os.path.join(directory, file_name)
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    lines = [line.strip() for line in lines if line.strip()]
    
    if len(lines) < 10:
        raise ValueError(f"File {file_name} does not contain enough lines for proper parsing.")
    
    default_metadata_lines = 13
    default_header_lines = 3
    try:
        num_ascii_tags = default_metadata_lines
        for line in lines[:default_metadata_lines+2]:
            if "num_ascii_tags" in line:
                parts = line.split(":", 1)
                if len(parts) == 2:
                    num_ascii_tags = int(parts[1].strip())
                    break
    except Exception as e:
        print(f"Warning: Unable to parse num_ascii_tags in file {file_name}. Using default {default_metadata_lines}.")
        num_ascii_tags = default_metadata_lines

    metadata = {}
    for i in range(num_ascii_tags):
        parts = lines[i].split(":", 1)
        if len(parts) == 2:
            key = parts[0].strip()
            value = parts[1].strip()
            metadata[key] = value
        else:
            metadata[f"line_{i}"] = lines[i]
    
    try:
        num_label_lines = default_header_lines
        if "num_label_lines" in metadata:
            num_label_lines = int(metadata["num_label_lines"])
    except Exception as e:
        print(f"Warning: Unable to parse num_label_lines in file {file_name}. Using default {default_header_lines}.")
        num_label_lines = default_header_lines

    header_start = num_ascii_tags
    header_end = num_ascii_tags + num_label_lines
    if len(lines) < header_end:
        raise ValueError(f"File {file_name} does not contain enough header lines. Expected {header_end} lines, found {len(lines)}.")

    if num_label_lines < 3:
        raise ValueError(f"File {file_name} does not provide enough header lines (found {num_label_lines}, expected at least 3).")
    
    sensor_names_line = lines[header_start]
    if ":" in sensor_names_line:
        sensor_names = sensor_names_line.split(":", 1)[1].split()
    else:
        sensor_names = sensor_names_line.split()
    
    sensor_units_line = lines[header_start+1]
    if ":" in sensor_units_line:
        sensor_units = sensor_units_line.split(":", 1)[1].split()
    else:
        sensor_units = sensor_units_line.split()
    
    sensor_rates_line = lines[header_start+2]
    if ":" in sensor_rates_line:
        sensor_rates = sensor_rates_line.split(":", 1)[1].split()
    else:
        sensor_rates = sensor_rates_line.split()
    
    min_length = min(len(sensor_names), len(sensor_units), len(sensor_rates))
    if not (len(sensor_names) == len(sensor_units) == len(sensor_rates)):
        print(f"Warning: Sensor header lengths do not match in file {file_name}.")
        print(f"sensor_names: {len(sensor_names)}, sensor_units: {len(sensor_units)}, sensor_rates: {len(sensor_rates)}")
        print(f"Using only the first {min_length} tokens for each header.")
        sensor_names = sensor_names[:min_length]
        sensor_units = sensor_units[:min_length]
        sensor_rates = sensor_rates[:min_length]
    
    sensor_data = {sensor: [] for sensor in sensor_names}
    
    data_start = header_end
    for row in lines[data_start:]:
        row_values = row.split()
        if len(row_values) < len(sensor_names):
            row_values += [None] * (len(sensor_names) - len(row_values))
        for sensor, value in zip(sensor_names, row_values):
            try:
                sensor_data[sensor].append(float(value))
            except:
                sensor_data[sensor].append(np.nan)
    
    result = {
        "metadata": metadata,
        "data": {}
    }
    for sensor, unit, rate in zip(sensor_names, sensor_units, sensor_rates):
        result["data"][sensor] = {
            "units": unit,
            "rate": rate,
            "values": sensor_data[sensor]
        }
    
    file_identifier = metadata.get("filename", file_name)
    return result, file_identifier

def dbd_directory_to_dict(directory, files_read=None, max_files=30):
    
    '''
    Convert a directory of DBD ASCII files into a master dictionary.
    
    For each file in the directory with an extension in the valid set,
    parse it and add the result to the master dictionary.
    
    Args:
      directory (str): The directory containing the ASCII files.
      files_read (list): Optional list of files already processed.
      max_files (int): Maximum number of files to process.
      
    Returns:
      master_dict (dict): Dictionary with keys as the base file names and values as parsed results.
      files_read (list): Updated list of processed file names.
      processed (int): Number of files processed.
      first_file (str): Name of the first processed file.
      last_file (str): Name of the last processed file.
    '''

    if files_read is None:
        files_read = []
    master_dict = {}
    processed = 0
    first_file = None
    last_file = None
    valid_extensions = ("dbd.asc", "sbd.asc", "DBD.asc", "SBD.asc", "ebd.asc")
    
    for file in os.listdir(directory):
        if processed >= max_files:
            break
        if file in files_read:
            continue
        if any(file.endswith(ext) for ext in valid_extensions):
            try:
                result, file_id = read_ascii(directory, file)
                file_key = os.path.splitext(file)[0]
                master_dict[file_key] = result
                files_read.append(file)
                if processed == 0:
                    first_file = file
                last_file = file
                processed += 1
            except Exception as e:
                print(f"Error processing {file}: {e}")
                continue
    if first_file:
        print(f"Files {first_file} to {last_file} processed.")
    else:
        print("No files processed.")
    return master_dict, files_read, processed, first_file, last_file

def pull_sensor_rate(sensor, master_dict):
    
    '''
    Retrieve the data rate for a specified sensor from the master dictionary.
    
    Args:
      sensor (str): Sensor name.
      master_dict (dict): Master dictionary of parsed files.
      
    Returns:
      rate (str): Sensor rate from the first file that contains it, or None.
    '''

    for content in master_dict.values():
        if "data" in content and sensor in content["data"]:
            return content["data"][sensor]["rate"]
    return None

def pull_sensor_units(sensor, master_dict):
    
    '''
    Retrieve the units for a specified sensor from the master dictionary.
    
    Args:
      sensor (str): Sensor name.
      master_dict (dict): Master dictionary of parsed files.
      
    Returns:
      units (str): Sensor units from the first file that contains it, or None.
    '''

    for content in master_dict.values():
        if "data" in content and sensor in content["data"]:
            return content["data"][sensor]["units"]
    return None

def pull_sensor_data(sensor, master_dict):
    
    '''
    Retrieve the concatenated data values for a specified sensor across all files.
    
    Args:
      sensor (str): Sensor name.
      master_dict (dict): Master dictionary of parsed files.
      
    Returns:
      data (list): List of all sensor values (floats) across files.
    '''

    data = []
    for content in master_dict.values():
        if "data" in content and sensor in content["data"]:
            data.extend(content["data"][sensor]["values"])
    return data

def pull_sensor_list_data(sensor_list, chosen_dataframe, master_dict):
    
    '''
    Retrieve data for a list of sensors from the master dictionary and append it to a DataFrame.
    
    This function:
      - Builds a dictionary of sensor data lists.
      - Determines the maximum length among these lists.
      - Pads any shorter lists with NaN so that all lists are equal in length.
      - Creates a DataFrame from the dictionary and concatenates it with the chosen_dataframe.
    
    Args:
      sensor_list (list): List of sensor names.
      chosen_dataframe (pd.DataFrame): DataFrame to which sensor data will be appended.
      master_dict (dict): Master dictionary of parsed files.
      
    Returns:
      chosen_dataframe (pd.DataFrame): Updated DataFrame with sensor data.
    '''

    data_dict = {}
    for sensor in sensor_list:
        try:
            sensor_values = pull_sensor_data(sensor, master_dict)
            data_dict[sensor] = sensor_values
        except Exception as e:
            print(f"Warning: Could not retrieve data for sensor '{sensor}': {e}")
            continue

    if not data_dict:
        return chosen_dataframe

    max_length = max(len(vals) for vals in data_dict.values())
    for sensor, vals in data_dict.items():
        if len(vals) < max_length:
            data_dict[sensor] = vals + [np.nan] * (max_length - len(vals))
    
    df_new = pd.DataFrame(data_dict)
    chosen_dataframe = pd.concat([chosen_dataframe, df_new], ignore_index=True)
    return chosen_dataframe

def run_dataframe(input_directory, sensor_list):
    
    '''
    Process sensor data from ASCII files in a directory and compile them into a single DataFrame.
    
    This function:
      - Scans the given directory for valid ASCII files.
      - Parses each file using read_ascii.
      - If sensor_list is empty, it builds the union of all sensor names from every file.
      - Compiles the data from all files into a single DataFrame.
    
    Args:
      input_directory (str): Directory containing ASCII files.
      sensor_list (list): List of sensors to include. If empty, all unique sensors are used.
      
    Returns:
      dataframe (pd.DataFrame): The compiled sensor data.
    '''

    files_processed = []
    master_dict, files_processed, processed, first_file, last_file = dbd_directory_to_dict(input_directory, files_processed)
    
    if processed == 0:
        print("No valid ASCII files were processed.")
        return pd.DataFrame()
    
    if not sensor_list:
        all_sensors = set()
        for content in master_dict.values():
            all_sensors.update(content["data"].keys())
        sensor_list = list(all_sensors)
        print(f"Sensor list was empty. Using union of all sensors: {sensor_list}")
    
    df = pd.DataFrame()
    df = pull_sensor_list_data(sensor_list, df, master_dict)
    return df
