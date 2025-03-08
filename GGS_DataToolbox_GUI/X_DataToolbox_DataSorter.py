"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import os
import shutil
import re
from datetime import datetime

# =========================

### FUNCTION:
def organize_data_files(input_data_folder, output_folder):

    '''
    Organize data files into mission-specific folders.

    Args:
      input_data_folder (str): The directory containing the input data files.
      output_folder (str): The directory where the organized data files will be saved.
      
    Returns:
      None
    '''

    file_groups = {}

    sys_log_path = os.path.join(input_data_folder, 'sys.log')
    sys_log_exists = os.path.exists(sys_log_path)

    for filename in os.listdir(input_data_folder):
        if filename.endswith('.dbd'):
            with open(os.path.join(input_data_folder, filename), 'r', encoding='utf-8', errors='ignore') as file:
                lines = file.readlines()
                mission_name_line = lines[7]
                mission_name = mission_name_line.split('mission_name:')[1].strip()
            base_name = os.path.splitext(filename)[0]
            if mission_name not in file_groups:
                file_groups[mission_name] = []
            file_groups[mission_name].append(base_name)

    for mission_name, base_names in file_groups.items():
        mission_folder = os.path.join(output_folder, mission_name)
        os.makedirs(mission_folder, exist_ok=True)
        for base_name in base_names:
            for filename in os.listdir(input_data_folder):
                if filename.startswith(base_name):
                    src_file = os.path.join(input_data_folder, filename)
                    shutil.move(src_file, mission_folder)
        
        if sys_log_exists:
            shutil.copy(sys_log_path, mission_folder)

### FUNCTION:
def organize_log_files(input_log_folder, output_folder):

    '''
    Organize log files by matching them with corresponding data files.

    Args:
      input_log_folder (str): The directory containing the input log files.
      output_folder (str): The directory where the organized log files will be saved.
      
    Returns:
      None
    '''

    for filename in os.listdir(input_log_folder):
        file_path = os.path.join(input_log_folder, filename)
        if os.path.getsize(file_path) == 0 or os.path.getsize(file_path) < 3072:
            os.remove(file_path)

    data_file_pattern = re.compile(r'\b\w+-\d{4}-\d{3}-\d+-\d+\b')

    for filename in os.listdir(input_log_folder):
        if filename.endswith('.log'):
            with open(os.path.join(input_log_folder, filename), 'r', encoding='utf-8', errors='ignore') as file:
                log_content = file.read()
                match = data_file_pattern.search(log_content)
                if match:
                    data_filename = match.group()
                    for mission_name in os.listdir(output_folder):
                        mission_folder = os.path.join(output_folder, mission_name)
                        if os.path.isdir(mission_folder):
                            for data_file in os.listdir(mission_folder):
                                if data_file.startswith(data_filename):
                                    src_file = os.path.join(input_log_folder, filename)
                                    dest_file = os.path.join(mission_folder, filename)
                                    file.close()
                                    shutil.move(src_file, dest_file)
                                    break

### FUNCTION:
def run_data_sorter(root_directory, glider_info):

    '''
    Run the data sorting process to organize data and log files.

    Args:
      root_directory (str): The root directory from the configuration file.
      glider_info (tuple): A tuple containing glider unit, version, and type information.
      
    Returns:
      None
    '''
    
    print(f"\n### RUNNING: DATA SORTER ###\n")
    datetime_string = datetime.now().strftime("%Y%m%dT%H%M%S")
    glider_unit, glider_version, glider_type = glider_info
    data_sort_directory = os.path.join(root_directory, f"GliderDataSorter")
    runtime_directory = os.path.join(data_sort_directory, f"{glider_unit}-{glider_type}-{glider_version}-{datetime_string}")
    os.makedirs(runtime_directory, exist_ok=True)
    current_directory = os.path.abspath(os.path.dirname(__file__))
    logfile_directory = os.path.join(current_directory, 'DBD_Files/Logfiles')
    data_directory = os.path.join(current_directory, 'DBD_Files')
    logfile_directory = os.path.join(current_directory, 'DBD_Files', 'Logfiles')
    organize_data_files(data_directory, runtime_directory)
    organize_log_files(logfile_directory, runtime_directory)
    print(f"\nGlider Data/Log Files saved to: {runtime_directory}\n")
