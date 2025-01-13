import os
import shutil
import re

### FUNCTION:
def organize_data_files(input_data_folder, output_folder):

    '''
    Organize data files into mission-specific folders.
    
    Args:
    - input_data_folder (str): The directory containing the input data files.
    - output_folder (str): The directory to save the organized data files.
    
    Returns:
    - None
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
    - input_log_folder (str): The directory containing the input log files.
    - output_folder (str): The directory to save the organized log files.
    
    Returns:
    - None
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
def run_data_sorter(output_directory):

    '''
    Run the data sorting process to organize data and log files.
    
    Args:
    - output_directory (str): The directory to save the organized files.
    
    Returns:
    - None
    '''

    print(f"\n### RUNNING: DATA SORTER ###\n")

    os.makedirs('DBD_Files/Logfiles', exist_ok=True)

    current_directory = os.path.dirname(__file__)
    data_directory = os.path.join(current_directory, 'DBD_Files')
    logfile_directory = os.path.join(current_directory, 'DBD_Files', 'Logfiles')

    organize_data_files(data_directory, output_directory)
    organize_log_files(logfile_directory, output_directory)
