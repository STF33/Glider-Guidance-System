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
def file_decompression(file_type, new_extension):

    '''
    Decompress files of a specified type to a new extension.
    
    Args:
      file_type (str): The file extension of the compressed files to decompress.
      new_extension (str): The new file extension after decompression.
      
    Returns:
      None
    '''

    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    for compressed_file in os.listdir('DBD_Files'):
        if compressed_file.endswith(file_type):
            source = os.path.join('DBD_Files', compressed_file)
            destination = os.path.join('DBD_Files/Decompressed', compressed_file.replace(file_type, new_extension))
            if not os.path.exists(destination):
                os.system(f'EXE_compexp.exe x {source} {destination}')
                print(f'{compressed_file} decompressed to {destination}')
            else:
                print(f'{compressed_file} already decompressed')
    
    print(f'Decompressing {file_type} files complete')

### FUNCTION:
def run_file_decompression():

    '''
    Run the decompression process for multiple file type mappings.
    
    Args:
      None
      
    Returns:
      None
    '''
    
    print(f"\n### RUNNING: DECOMPRESSION ###\n")
    current_directory = os.path.abspath(os.path.dirname(__file__))
    decompressed_directory = os.path.join(current_directory, 'DBD_Files/Decompressed')
    os.makedirs(decompressed_directory, exist_ok=True)
    
    file_mappings = {
        '.dcd': '.dbd',
        '.ecd': '.ebd',
        '.scd': '.sbd',
        '.tcd': '.tbd',
        '.mcd': '.mbd',
        '.ncd': '.nbd',
        '.mcg': '.mlg',
        '.ncg': '.nlg'
    }
    
    for file_type, new_extension in file_mappings.items():
        file_decompression(file_type, new_extension)
