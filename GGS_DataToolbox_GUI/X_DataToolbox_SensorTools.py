"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System

This module provides functions for parsing ASCII data files and compiling sensor data into a DataFrame.
"""

# =========================
# IMPORTS
# =========================

import os
import pandas as pd
import numpy as np
from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLineEdit, QComboBox, QPushButton, QLabel, QScrollArea, QWidget
from PyQt6.QtCore import Qt

# =========================

OPERATORS = {
    '<': lambda series, val: series < val,
    '<=': lambda series, val: series <= val,
    '=': lambda series, val: series == val,
    '==': lambda series, val: series == val,
    '>=': lambda series, val: series >= val,
    '>': lambda series, val: series > val
}

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
    - directory (str): The directory containing the ASCII file.
    - file_name (str): The name of the ASCII file.
      
    Returns:
    - result (dict): A dictionary with keys:
        - "metadata": dict of metadata,
        - "data": dict mapping each sensor name to a dict with keys "units", "rate", and "values" (a list of floats).
    - file_identifier (str): The identifier for the file (from metadata if available, else the file name).
    '''

    filepath = os.path.join(directory, file_name)
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    lines = [line.strip() for line in lines if line.strip()]
    
    if len(lines) < 16:
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
    sensor_names = sensor_names_line.split(":", 1)[-1].split() if ":" in sensor_names_line else sensor_names_line.split()
    
    sensor_units_line = lines[header_start+1]
    sensor_units = sensor_units_line.split(":", 1)[-1].split() if ":" in sensor_units_line else sensor_units_line.split()
    
    sensor_rates_line = lines[header_start+2]
    sensor_rates = sensor_rates_line.split(":", 1)[-1].split() if ":" in sensor_rates_line else sensor_rates_line.split()
    
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

### FUNCTION:
def dbd_directory_to_dict(directory, files_read=None, max_files=30, selected_files=None):
    
    '''
    Convert a list of ASCII files into a master dictionary.
    
    If selected_files is provided, only those files (full file names) are processed.
    
    Args:
    - directory (str): The directory containing the ASCII files.
    - files_read (list): Optional list of files already processed.
    - max_files (int): Maximum number of files to process.
    - selected_files (list): Optional list of file names to process.
      
    Returns:
    - master_dict (dict): Dictionary with keys as the base file names and values as parsed results.
    - files_read (list): Updated list of processed file names.
    - processed (int): Number of files processed.
    - first_file (str): Name of the first processed file.
    - last_file (str): Name of the last processed file.
    '''

    if files_read is None:
        files_read = []
    master_dict = {}
    processed = 0
    first_file = None
    last_file = None
    
    if selected_files is None:
        file_list = os.listdir(directory)
    else:
        file_list = selected_files

    valid_extensions = ("dbd.asc", "sbd.asc", "DBD.asc", "SBD.asc", "ebd.asc")
    
    for file in file_list:
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

### FUNCTION:
def pull_sensor_rate(sensor, master_dict):
    
    '''
    Retrieve the data rate for a specified sensor from the master dictionary.
    
    Args:
    - sensor (str): Sensor name.
    - master_dict (dict): Master dictionary of parsed files.
      
    Returns:
    - rate (str): Sensor rate from the first file that contains it, or None.
    '''

    for content in master_dict.values():
        if "data" in content and sensor in content["data"]:
            return content["data"][sensor]["rate"]
    return None

### FUNCTION:
def pull_sensor_units(sensor, master_dict):
    
    '''
    Retrieve the units for a specified sensor from the master dictionary.
    
    Args:
    - sensor (str): Sensor name.
    - master_dict (dict): Master dictionary of parsed files.
      
    Returns:
    - units (str): Sensor units from the first file that contains it, or None.
    '''

    for content in master_dict.values():
        if "data" in content and sensor in content["data"]:
            return content["data"][sensor]["units"]
    return None

### FUNCTION:
def pull_sensor_data(sensor, master_dict):
    
    '''
    Retrieve the concatenated data values for a specified sensor across all files.
    
    Args:
    - sensor (str): Sensor name.
    - master_dict (dict): Master dictionary of parsed files.
      
    Returns:
    - data (list): List of all sensor values (floats) across files.
    '''

    data = []
    for content in master_dict.values():
        if "data" in content and sensor in content["data"]:
            data.extend(content["data"][sensor]["values"])
    return data

### FUNCTION:
def pull_sensor_list_data(sensor_list, chosen_dataframe, master_dict):
    
    '''
    Retrieve data for a list of sensors from the master dictionary and append it to a DataFrame.
    
    Args:
    - sensor_list (list): List of sensor names.
    - chosen_dataframe (pd.DataFrame): DataFrame to which sensor data will be appended.
    - master_dict (dict): Master dictionary of parsed files.
      
    Returns:
    - chosen_dataframe (pd.DataFrame): Updated DataFrame with sensor data.
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

FLIGHT_EXT_PRIORITY = {'.dbd': 1, '.mbd': 2, '.sbd': 3}
SCIENCE_EXT_PRIORITY = {'.ebd': 1, '.nbd': 2, '.tbd': 3}

### FUNCTION:
def group_files_by_category(directory):
    
    '''
    Group files (ending with ".asc") into flight and science groups based on their extension.

    Args:
    - directory (str): The directory containing the files to be grouped.

    Returns:
    - flight_files (list): List of files belonging to the flight group.
    - science_files (list): List of files belonging to the science group.
    '''

    flight_files = []
    science_files = []
    for file in os.listdir(directory):
        if not file.endswith(".asc"):
            continue
        base_file = file[:-4]
        base, ext = os.path.splitext(base_file)
        ext = ext.lower()
        if ext in FLIGHT_EXT_PRIORITY:
            flight_files.append(file)
        elif ext in SCIENCE_EXT_PRIORITY:
            science_files.append(file)
    return flight_files, science_files

### FUNCTION:
def select_priority_files(file_list, priority_map):
    
    '''
    From a list of files (all ending with ".asc"), group by their base name (excluding the extension part before ".asc") and select one file per group based on the given priority mapping.
    
    Args:
    - file_list (list): List of file names.
    - priority_map (dict): Mapping from extension (e.g. ".dbd") to a numeric priority (lower value = higher priority).
    
    Returns:
    - selected_files (list): List of file names chosen based on the priority.
    '''

    groups = {}
    for file in file_list:
        base_file = file[:-4]
        base, ext = os.path.splitext(base_file)
        ext = ext.lower()
        groups.setdefault(base, []).append((file, priority_map.get(ext, float('inf'))))
    
    selected_files = []
    for base, files in groups.items():
        best_file = min(files, key=lambda x: x[1])[0]
        selected_files.append(best_file)
    return selected_files

### FUNCTION:
def process_files(file_list, directory):
    
    '''
    Process a given list of files from the specified directory using read_ascii.

    Args:
    - file_list (list): List of files to be processed.
    - directory (str): The directory containing the files.

    Returns:
    - master_dict (dict): Dictionary containing the processed file data.
    - files_read (list): List of files that were successfully read.
    - processed (int): Number of files that were processed.
    - first_file (str): The first file that was processed.
    - last_file (str): The last file that was processed.
    '''

    files_read = []
    master_dict = {}
    processed = 0
    first_file = None
    last_file = None
    for file in file_list:
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

### FUNCTION:
def get_units_dict(sensor_list, master_dict):
    
    '''
    Build a dictionary mapping each sensor in sensor_list to its unit based on the first file in master_dict that provides a value.

    Args:
    - sensor_list (list): List of sensor names.
    - master_dict (dict): Master dictionary of parsed files.

    Returns:
    - units_dict (dict): Mapping from sensor name to its unit (or an empty string if not found).
    '''

    units_dict = {}
    for sensor in sensor_list:
        unit = ""
        for content in master_dict.values():
            if "data" in content and sensor in content["data"]:
                unit = content["data"][sensor].get("units", "")
                break
        units_dict[sensor] = unit
    return units_dict

### FUNCTION:
def merge_flight_science(flight_df, science_df, flight_master, science_master, sensor_list):
    
    '''
    Merge the flight and science DataFrames by creating a new column "time" in each, then concatenates and sorts by "time". Also, build a dictionary of units for each sensor.

    For flight_df, "time" is computed from the column 'm_present_time'.
    For science_df, "time" is computed from 'sci_m_present_time'.

    Args:
    - flight_df (pd.DataFrame): DataFrame from flight files.
    - science_df (pd.DataFrame): DataFrame from science files.
    - flight_master (dict): Master dictionary from flight files.
    - science_master (dict): Master dictionary from science files.
    - sensor_list (list): List of sensor names used.

    Returns:
    - glider_df (pd.DataFrame): Merged DataFrame sorted by the new "time" column with units stored in glider_df.attrs["units"].
    '''

    if not flight_df.empty and "m_present_time" in flight_df.columns:
        flight_df = flight_df.copy()
        flight_df["time"] = pd.to_datetime(
            pd.to_numeric(flight_df["m_present_time"], errors='coerce'),
            unit='s', errors='coerce'
        )

    if not science_df.empty and "sci_m_present_time" in science_df.columns:
        science_df = science_df.copy()
        science_df["time"] = pd.to_datetime(
            pd.to_numeric(science_df["sci_m_present_time"], errors='coerce'),
            unit='s', errors='coerce'
        )
    
    glider_df = pd.concat([flight_df, science_df], ignore_index=True)
    if "time" in glider_df.columns:
        glider_df = glider_df.sort_values(by="time").reset_index(drop=True)
    
    combined_master = {}
    combined_master.update(flight_master)
    combined_master.update(science_master)
    units_dict = get_units_dict(sensor_list, combined_master)
    glider_df.attrs["units"] = units_dict

    return glider_df

### FUNCTION:
def run_dataframe(input_directory, sensor_list):
    
    '''
    Process sensor data from ASCII files in a directory and compile them into a single merged DataFrame.
    
    Args:
    - input_directory (str): Directory containing ASCII files.
    - sensor_list (list): List of sensors to include. If empty, all unique sensors are used.
      
    Returns:
    - glider_df (pd.DataFrame): The merged sensor data with a "time" column, and unit metadata in attrs.
    '''

    flight_files, science_files = group_files_by_category(input_directory)
    
    selected_flight = select_priority_files(flight_files, FLIGHT_EXT_PRIORITY)
    selected_science = select_priority_files(science_files, SCIENCE_EXT_PRIORITY)
    
    print("Processing flight files:")
    flight_master, _, flight_processed, f_first, f_last = process_files(selected_flight, input_directory)
    print("Processing science files:")
    science_master, _, science_processed, s_first, s_last = process_files(selected_science, input_directory)
    
    if not sensor_list:
        flight_sensors = set()
        for content in flight_master.values():
            flight_sensors.update(content["data"].keys())
        science_sensors = set()
        for content in science_master.values():
            science_sensors.update(content["data"].keys())
        sensor_list = list(flight_sensors.union(science_sensors))
        print(f"Sensor list was empty. Using union of all sensors: {sensor_list}")
    
    flight_df = pd.DataFrame()
    flight_df = pull_sensor_list_data(sensor_list, flight_df, flight_master)
    
    science_df = pd.DataFrame()
    science_df = pull_sensor_list_data(sensor_list, science_df, science_master)
    
    glider_df = merge_flight_science(flight_df, science_df, flight_master, science_master, sensor_list)

    return glider_df

### CLASS:
class DataFilterDialog(QDialog):

    '''
    Create a dialog for adding multiple data filter criteria.
    
    Args:
    - dataframe (pd.DataFrame): The input dataframe to be filtered.
      
    Returns:
    - None
    '''

    def __init__(self, dataframe, parent=None):
        super().__init__(parent)
        self.dataframe = dataframe.copy()
        self.filter_rows = []
        self.filtered_df = None
        self.initUI()
    
    def initUI(self):

        '''
        Initialize the filter dialog layout with a scroll area and buttons.
        
        Args:
        - None
          
        Returns:
        - None
        '''

        self.setWindowTitle("Data Filter")
        self.setMinimumSize(500, 300)
        self.layout = QVBoxLayout(self)
        
        self.filterWidget = QWidget()
        self.filterLayout = QVBoxLayout(self.filterWidget)
        self.filterWidget.setLayout(self.filterLayout)
        
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.filterWidget)
        self.layout.addWidget(self.scroll)
        
        self.addFilterButton = QPushButton("+ Add Filter")
        self.addFilterButton.clicked.connect(self.add_filter_row)
        self.layout.addWidget(self.addFilterButton)
        
        self.runFilterButton = QPushButton("Run Data Filter")
        self.runFilterButton.clicked.connect(self.apply_filters)
        self.layout.addWidget(self.runFilterButton)
    
    def add_filter_row(self):

        '''
        Add a new filter row to the dialog. Each row contains a QLineEdit for the variable, a QComboBox for the operator, and a QLineEdit for the value.
        
        Args:
        - None
          
        Returns:
        - None
        '''

        row_widget = QWidget()
        row_layout = QHBoxLayout(row_widget)
        
        var_edit = QLineEdit()
        var_edit.setPlaceholderText("Variable")
        row_layout.addWidget(var_edit)
        
        op_combo = QComboBox()
        op_combo.addItems(list(OPERATORS.keys()))
        row_layout.addWidget(op_combo)
        
        val_edit = QLineEdit()
        val_edit.setPlaceholderText("Value")
        row_layout.addWidget(val_edit)
        
        self.filterLayout.addWidget(row_widget)
        self.filter_rows.append((var_edit, op_combo, val_edit))
    
    def apply_filters(self):

        '''
        Apply all filter rows to the dataframe and close the dialog.
        
        Args:
        - None
          
        Returns:
        - None
        '''

        df_filtered = self.dataframe.copy()
        for var_edit, op_combo, val_edit in self.filter_rows:
            var = var_edit.text().strip()
            op = op_combo.currentText().strip()
            val_str = val_edit.text().strip()
            if var == "" or val_str == "":
                continue
            try:
                try:
                    val = float(val_str)
                except ValueError:
                    val = val_str
                if var not in df_filtered.columns:
                    print(f"Warning: Variable '{var}' not found in dataframe. Skipping filter.")
                    continue
                condition = OPERATORS[op](pd.to_numeric(df_filtered[var], errors='coerce'), val)
                df_filtered = df_filtered[condition]
            except Exception as e:
                print(f"Error applying filter on {var} {op} {val_str}: {e}")
        self.filtered_df = df_filtered
        self.accept()

### FUNCTION:
def run_data_filter(dataframe):

    '''
    Run the data filter dialog on the input dataframe.
    
    Args:
    - dataframe (pd.DataFrame): The input dataframe to filter.
      
    Returns:
    - filtered_df (pd.DataFrame): The filtered dataframe after applying user-selected criteria.
    '''

    dialog = DataFilterDialog(dataframe)
    if dialog.exec() == QDialog.DialogCode.Accepted:
        return dialog.filtered_df
    else:
        return dataframe
