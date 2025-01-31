"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

from X_GGS_DataSorter import *
from X_GGS_FlightDataEval import *
from X_GGS_LowPowerEval import *
from X_GGS_YoEval import *

from X_DataToolbox_ClearDataLibrary import *
from X_DataToolbox_Config import *
from X_DataToolbox_Converter import *
from X_DataToolbox_Decompression import *
from X_DataToolbox_Filter import *
from X_DataToolbox_Rename import *
from X_DataToolbox_SensorTools import *
from X_DataToolbox_Excel import *

# =========================

### FUNCTION:
def GGS_DataToolbox_Main(config_name=None):
    
    '''
    GGS: Data Toolbox main function.

    Args:
    - config_name (str): The name of the config file without the extension.
    - path (str): The path directory to save output to. Options: 'local' or a specific directory path.
    
    Returns:
    - None
    '''

    current_directory = os.path.abspath(os.path.dirname(__file__))

    cache_directory = os.path.join(current_directory, 'cache')
    os.makedirs(cache_directory, exist_ok=True)
    decompressed_directory = os.path.join(current_directory, 'DBD_Files/Decompressed')
    os.makedirs(decompressed_directory, exist_ok=True)
    logfile_directory = os.path.join(current_directory, 'DBD_Files/Logfiles')
    os.makedirs(logfile_directory, exist_ok=True)
    ascii_directory = os.path.join(current_directory, 'DBD_Files/ProcessedAscii')
    os.makedirs(ascii_directory, exist_ok=True)
    renamed_directory = os.path.join(current_directory, 'DBD_Files/RenamedBinary')
    os.makedirs(renamed_directory, exist_ok=True)

    config = GGS_config_import(config_name)
    root_directory = GGS_config_process(config, path="default")

    glider_config = config['GLIDER']
    glider_info = (glider_config.get('glider_unit'), glider_config.get('glider_version'), glider_config.get('glider_type'))

    sensor_config = config['SENSORS']
    sensor_list = sensor_config.get('sensor_list')

    if config is None:
        print("Failed to import configuration.")
        return
    
    if 'DECOMPRESSION' in config:
        decompression_config = config['DECOMPRESSION']
        if decompression_config.get('run_decompression'):
            run_file_decompression()
        else:
            pass
    
    if 'CONVERSION' in config:
        conversion_config = config['CONVERSION']
        if conversion_config.get('run_conversion'):
            run_ascii_converter()
        else:
            pass
    
    if 'DATA' in config:
        data_config = config['DATA']
        if data_config.get('run_dataframe'):
            current_directory = os.path.dirname(__file__)
            input_directory = os.path.join(current_directory, 'DBD_Files', 'ProcessedAscii')
            dataframe = run_dataframe(input_directory, sensor_list)
        else:
            pass

    if 'PRODUCTS' in config:
        product_config = config['PRODUCTS']
        # Generic Products
        if product_config.get('run_excel'):
            run_excel(root_directory, glider_info, dataframe)
        else:
            pass
        if product_config.get('run_data_sorter'):
            run_data_sorter(root_directory, glider_info)
        else:
            pass
        if product_config.get('run_fde'):
            run_fde(dataframe)
        else:
            pass
        # Specific Products
        if product_config.get('run_low_power_eval'):
            run_low_power_eval(root_directory, glider_info)
        elif product_config.get('run_yo_eval'):
            run_yo_eval(root_directory, glider_info)
        else:
            pass
    
    if 'ADVANCED' in config:
        advanced_config = config['ADVANCED']
        if advanced_config.get('run_data_cleanup'):
            directories_to_clean = define_directories_to_clean()
            run_data_cleanup(directories_to_clean)
        else:
            pass

if __name__ == "__main__":
    GGS_DataToolbox_Main(config_name="yo_eval")