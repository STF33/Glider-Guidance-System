# =========================
# IMPORTS
# =========================

import os

# =========================

### FUNCTION:
def define_directories_to_clean():
    
    '''
    Define directories which contain run-specific data files.
    
    Args:
    - None
    
    Returns:
    - directories_to_clean (list): A list of directory paths.
    '''
    
    current_directory = os.path.dirname(__file__)
    directories_to_clean = [
        os.path.join(current_directory, 'DBD_Files'),
        os.path.join(current_directory, 'DBD_Files', 'ProcessedAscii'),
        os.path.join(current_directory, 'DBD_Files', 'RenamedBinary'),
        os.path.join(current_directory, 'DBD_Files', 'Decompressed'),
        os.path.join(current_directory, 'DBD_Files', 'Logfiles')
    ]
    
    return directories_to_clean

### FUNCTION:
def run_data_cleanup(directories_to_clean):

    '''
    Delete data files in the data directories.
    
    Args:
    - directories_to_clean (list): A list of directory paths.
    
    Returns:
    - None
    '''

    print(f"\n### RUNNING: DATA LIBRARY CLEANUP ###\n")

    for directory in directories_to_clean:
        if os.path.exists(directory):
            for item in os.listdir(directory):
                item_path = os.path.join(directory, item)
                if os.path.isfile(item_path):
                    os.remove(item_path)
                    print(f"Deleted file: {item_path}")
                else:
                    print(f"Skipped folder: {item_path}")
        else:
            print(f"Directory does not exist: {directory}")

