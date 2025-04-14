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

    config = config_import(config_name)
    root_directory = config_process(config, path="default")
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
            decompressor_run()
    if 'CONVERSION' in config:
        conversion_config = config['CONVERSION']
        if conversion_config.get('run_conversion'):
            converter_run()

    dataframe = None
    if 'DATA' in config:
        data_config = config['DATA']
        if data_config.get('run_dataframe'):
            input_directory = os.path.join(current_directory, 'DBD_Files', 'ProcessedAscii')
            dataframe = glider_dataframe_run(input_directory, sensor_list)
        if data_config.get('run_data_filter') and dataframe is not None:
            dataframe = datafilter_run(dataframe)
        if data_config.get('run_data_sorter'):
            datasorter_run(root_directory, glider_info)
        if data_config.get('run_logfile_search'):
            logsearch_run(root_directory)

    if 'PRODUCTS' in config:
        product_config = config['PRODUCTS']
        if product_config.get('run_plot') and dataframe is not None:
            plot_run(dataframe)
        if product_config.get('run_excel') and dataframe is not None:
            excel_run(root_directory, glider_info, dataframe)
        if product_config.get('run_energy_evaluation') and dataframe is not None:
            energy_run(root_directory, glider_info, dataframe)

    if 'ADVANCED' in config:
        advanced_config = config['ADVANCED']
        if advanced_config.get('run_data_cleanup'):
            directories_to_clean = cleanup_define_directories()
            cleanup_run(directories_to_clean)

    return dataframe, root_directory

