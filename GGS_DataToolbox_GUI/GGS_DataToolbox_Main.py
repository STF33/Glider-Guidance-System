"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

from X_DataToolbox_Data import *
from X_DataToolbox_Products import *
from X_DataToolbox_Advanced import *

# =========================

def GGS_DataToolbox_Main(config_name=None):

    ''' 
    Process the glider guidance system data toolbox by running decompression, conversion,
    data compilation, and product generation based on the provided configuration.
    
    Arguments:
    - config_name (str): The name of the config file without the extension.
      
    Returns:
    - dataframe (pd.DataFrame): The merged sensor data DataFrame (if created), otherwise None.
    - root_directory (str): The output directory as defined by the configuration.
    '''

    import os
    current_directory = os.path.abspath(os.path.dirname(__file__))

    cache_directory = os.path.join(current_directory, 'cache')
    os.makedirs(cache_directory, exist_ok=True)
    decompressed_directory = os.path.join(current_directory, 'DBD_Files', 'Decompressed')
    os.makedirs(decompressed_directory, exist_ok=True)
    logfile_directory = os.path.join(current_directory, 'DBD_Files', 'Logfiles')
    os.makedirs(logfile_directory, exist_ok=True)
    ascii_directory = os.path.join(current_directory, 'DBD_Files', 'ProcessedAscii')
    os.makedirs(ascii_directory, exist_ok=True)

    config = GGS_config_import(config_name)
    root_directory = GGS_config_process(config, path="default")
    glider_config = config['GLIDER']
    glider_info = (glider_config.get('glider_unit'), glider_config.get('glider_version'), glider_config.get('glider_type'))
    sensor_config = config['SENSORS']
    sensor_list = sensor_config.get('sensor_list')

    if config is None:
        print("Failed to import configuration.")
        return None, root_directory

    if 'DECOMPRESSION' in config:
        decompression_config = config['DECOMPRESSION']
        if decompression_config.get('run_decompression'):
            run_file_decompression()
    if 'CONVERSION' in config:
        conversion_config = config['CONVERSION']
        if conversion_config.get('run_conversion'):
            run_ascii_converter()

    dataframe = None
    if 'DATA' in config:
        data_config = config['DATA']
        if data_config.get('run_dataframe'):
            input_directory = os.path.join(current_directory, 'DBD_Files', 'ProcessedAscii')
            dataframe = run_dataframe(input_directory, sensor_list)
        if data_config.get('run_data_filter') and dataframe is not None:
            dataframe = run_data_filter(dataframe)
        if data_config.get('run_data_sorter'):
            run_data_sorter(root_directory, glider_info)
        if data_config.get('run_logfile_search'):
            run_logfile_search(root_directory)

    if 'PRODUCTS' in config:
        product_config = config['PRODUCTS']
        if product_config.get('run_plot') and dataframe is not None:
            run_plot(dataframe)
        if product_config.get('run_excel') and dataframe is not None:
            run_excel(root_directory, glider_info, dataframe)
        if product_config.get('run_energy_evaluation') and dataframe is not None:
            run_energy_evaluation(root_directory, glider_info, dataframe)

    if 'ADVANCED' in config:
        advanced_config = config['ADVANCED']
        if advanced_config.get('run_data_cleanup'):
            directories_to_clean = define_directories_to_clean()
            run_data_cleanup(directories_to_clean)

    return dataframe, root_directory

if __name__ == "__main__":
    GGS_DataToolbox_Main(config_name="config")
