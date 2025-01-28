"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import os


# =========================

### FUNCTION:
def ascii_converter(file_type):

    '''
    Convert binary files of a specified type to ASCII format.
    
    Args:
    - file_type (str): The file extension of the binary files to convert.
    
    Returns:
    - None
    '''

    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    for binaryfile in os.listdir('DBD_Files'):
        if binaryfile.endswith(file_type):
            source = os.path.join('DBD_Files', binaryfile)
            destination = os.path.join('DBD_Files/ProcessedAscii', binaryfile + '.asc')
            if not os.path.exists(destination):
                os.system(f'EXE_dbd2asc {source} > {destination}')
                print(f'{binaryfile} converted')
            else:
                print(f'{binaryfile} already converted')
    
    print(f'Converting {file_type} files to ASCII complete')

### FUNCTION:
def run_ascii_converter():

    '''
    Run the ASCII conversion process for multiple binary file types.
    
    Args:
    - None
    
    Returns:
    - None
    '''
    
    print(f"\n### RUNNING: CONVERTER ###\n")

    current_directory = os.path.abspath(os.path.dirname(__file__))
    ascii_directory = os.path.join(current_directory, 'DBD_Files/ProcessedAscii')
    os.makedirs(ascii_directory, exist_ok=True)

    file_types = ['.sbd', '.dbd', '.tbd', '.ebd']
    
    for file_type in file_types:
        ascii_converter(file_type)