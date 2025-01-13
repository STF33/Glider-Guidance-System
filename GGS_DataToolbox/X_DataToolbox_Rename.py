# =========================
# IMPORTS
# =========================

import os

# =========================

### FUNCTION:
def file_rename(file_type):

    '''
    Rename binary files of a specified type.
    
    Args:
    - file_type (str): The file extension of the binary files to rename.
    
    Returns:
    - None
    '''

    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    for file in os.listdir('DBD_Files'):
        if file.endswith(file_type):
            source = os.path.join('DBD_Files', file)
            destination = os.path.join('DBD_Files/RenamedBinary', file)
            if not os.path.exists(destination):
                os.system(f'EXE_rename_dbd_files {source} > {destination}')
                print(f'{file} converted')
            else:
                print(f'{file} already converted')
    
    print(f'Renaming {file_type} files complete')

### FUNCTION:
def run_file_rename():

    '''
    Run the renaming process for multiple binary file types.
    
    Args:
    - None
    
    Returns:
    - None
    '''
    
    print(f"\n### RUNNING: RENAME ###\n")

    os.makedirs('DBD_Files/RenamedBinary', exist_ok=True)
    file_types = ['.sbd', '.scd', '.dbd', '.dcd', '.tbd', '.tcd', '.ebd', '.ecd']
    
    for file_type in file_types:
        file_rename(file_type)
