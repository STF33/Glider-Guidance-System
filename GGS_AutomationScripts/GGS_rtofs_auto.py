# =========================
# IMPORTS
# =========================

from A_functions import *
from A_config import *
from A_models import *
from A_interpolation import *
from A_plots import *

from concurrent.futures import ProcessPoolExecutor

# =========================
# MAIN
# =========================

### PARALELLIZE DATETIME PROCESSING:
def GGS_process_date(datetime_index, config, root_directory, enable_rtofs, enable_cmems, enable_gofs, latitude_qc, longitude_qc, glider_data):
    
    '''
    Process a single datetime index.

    Args:
    - datetime_index (int): Index of the date to process in date_list.
    - config (dict): Configuration dictionary.
    - root_directory (str): Root directory to save output to.
    - enable_rtofs (bool): Flag to process RTOFS data.
    - enable_cmems (bool): Flag to process CMEMS data.
    - enable_gofs (bool): Flag to process GOFS data.
    - latitude_qc (str): Latitude for QC plotting.
    - longitude_qc (str): Longitude for QC plotting.
    - glider_data (pd.DataFrame): Glider data.

    Returns:
    - None
    '''
    
    # INITIALIZATION
    sub_directory_plots = os.path.join(root_directory, "plots", ''.join(datetime_index[:10].split('-')))
    os.makedirs(sub_directory_plots, exist_ok=True)
    sub_directory_data = os.path.join(root_directory, "data", ''.join(datetime_index[:10].split('-')))
    os.makedirs(sub_directory_data, exist_ok=True)
    
    check_datetime = pd.to_datetime(datetime_index).strftime('%Y%m%dT%HZ')
    check_datafile = os.path.join(sub_directory_data, f"{config['glider_name']}_DepthAverageData_{check_datetime}.nc")

    if os.path.exists(check_datafile):
        print(f'{check_datafile} already exists. Skipping processing for datetime index: {datetime_index}')
        return
    
    # MODEL DATA PROCESSING
    model_datasets = []
    if enable_rtofs:
        try:
            rtofs = RTOFS()
            rtofs.rtofs_load(config, datetime_index, subset_extent=True, subset_depth=True)
            rtofs.rtofs_save(config, sub_directory_data, save_data=False, save_qc=False)
            rtofs_model_data = rtofs.data
            
            rtofs_depth_average, rtofs_bin_average = interpolate_rtofs(config, sub_directory_data, rtofs_model_data, chunk=False, save_depth_average=True, save_bin_average=False)

            rtofs_datasets = (rtofs_model_data, rtofs_depth_average, rtofs_bin_average)
            model_datasets.append(rtofs_datasets)
        except Exception as e:
            print(f"Error during RTOFS processing: {e}")
    if enable_cmems:
        try:
            cmems = CMEMS(username='sfricano1', password='GlobalGliders1')
            cmems.cmems_load(config, datetime_index, subset_extent=True, subset_depth=True)
            cmems.cmems_save(config, sub_directory_data, save_data=False, save_qc=False)
            cmems_model_data = cmems.data
            
            cmems_depth_average, cmems_bin_average = interpolate_cmems(config, sub_directory_data, cmems_model_data, chunk=False, save_depth_average=True, save_bin_average=False)
            
            cmems_datasets = (cmems_model_data, cmems_depth_average, cmems_bin_average)
            model_datasets.append(cmems_datasets)
        except Exception as e:
            print(f"Error during CMEMS processing: {e}")
    if enable_gofs:
        try:
            gofs = GOFS()
            gofs.gofs_load(config, datetime_index, subset_extent=True, subset_depth=True)
            gofs.gofs_save(config, sub_directory_data, save_data=False, save_qc=False)
            gofs_model_data = gofs.data

            gofs_depth_average, gofs_bin_average = interpolate_gofs(config, sub_directory_data, gofs_model_data, chunk=False, save_depth_average=True, save_bin_average=False)

            gofs_datasets = (gofs_model_data, gofs_depth_average, gofs_bin_average)
            model_datasets.append(gofs_datasets)
        except Exception as e:
            print(f"Error during GOFS processing: {e}")
    
    # PLOTTING
    GGS_plot_magnitudes(config, sub_directory_plots, datetime_index, model_datasets, latitude_qc=latitude_qc, longitude_qc=longitude_qc, density=2, gliders=glider_data, show_route=False, show_qc=False, manual_extent=None)
    GGS_plot_threshold(config, sub_directory_plots, datetime_index, model_datasets, latitude_qc=latitude_qc, longitude_qc=longitude_qc, density=2, mag1=0.0, mag2=0.2, mag3=0.3, mag4=0.4, mag5=0.5, gliders=glider_data, show_route=False, show_qc=False, manual_extent=None)

### MAIN:
def main(power=1, path="local", target_datetime=dt.datetime.now(dt.timezone.utc), single_datetime=False, enable_rtofs=True, enable_cmems=True, enable_gofs=True, latitude_qc='0', longitude_qc='0'):
    
    '''
    GGS main function.

    Args:
    - power (int): Percentage of available resources to use in processing.
        - Default: '1'
    - path (str): Path directory to save output to.
        - Options: 'local' or 'rucool'.
    - target_datetime (datetime): Target date to process.
        - Default: 'dt.datetime.now(dt.timezone.utc)'
    - single (bool): Process a single date.
        - Default: 'False'
    - enable_rtofs (bool): Flag to process RTOFS data.
        - Default: 'True'
    - enable_cmems (bool): Flag to process CMEMS data.
        - Default: 'True'
    - enable_gofs (bool): Flag to process GOFS data.
        - Default: 'True'
    - latitude_qc (str): Latitude for QC plotting.
        - Default: '0'
    - longitude_qc (str): Longitude for QC plotting.
        - Default: '0'

    Returns:
    - None
    '''

    # CONFIGURATION
    config = GGS_config_static(date=target_datetime)
    
    # TARGET DIRECTORY
    if path == "local":
        root_directory = GGS_config_output(config, path="default")
    elif path == "rucool":
        root_directory = GGS_config_output(config, path="/www/web/rucool/hurricane/model_comparisons/maps/yucatan")
    else:
        raise ValueError("Invalid root directory. Please set root to either 'local' or 'rucool'.")
    
    # TARGET DATE LIST
    if single_datetime:
        datetime_list = [target_datetime.replace(hour=0, minute=0, second=0, microsecond=0).strftime('%Y-%m-%dT%H:%M:%SZ')]
    else:
        datetime_start = target_datetime.replace(hour=0, minute=0, second=0, microsecond=0)
        datetime_end = datetime_start + dt.timedelta(days=1)
        datetime_range = pd.date_range(datetime_start, datetime_end, freq='6H').tz_localize(None)
        datetime_list = [datetime.strftime('%Y-%m-%dT%H:%M:%SZ') for datetime in datetime_range]
    
    # GLIDER DATA PROCESSING
    glider_dataframes = None
    if any([enable_rtofs, enable_cmems, enable_gofs]):
        search_extent = [config['extent'][0][1], config['extent'][1][1], config['extent'][0][0], config['extent'][1][0]]
        glider_dataframes = acquire_gliders(extent=search_extent, target_date=dt.datetime.now(), date_delta=dt.timedelta(days=1), requested_variables=["time", "longitude", "latitude", "profile_id", "depth"], print_vars=False, target="all", request_timeout=5, enable_parallel=False)

    # PARALELL PROCESSING
    config_flag = [config] * len(datetime_list)
    root_directory_flag = [root_directory] * len(datetime_list)
    enable_rtofs_flag = [enable_rtofs] * len(datetime_list)
    enable_cmems_flag = [enable_cmems] * len(datetime_list)
    enable_gofs_flag = [enable_gofs] * len(datetime_list)
    latitude_qc_flag = [latitude_qc] * len(datetime_list)
    longitude_qc_flag = [longitude_qc] * len(datetime_list)
    glider_data_flag = [glider_dataframes] * len(datetime_list)

    num_workers = optimal_workers(power=power)
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        executor.map(GGS_process_date, datetime_list, config_flag, root_directory_flag, enable_rtofs_flag, enable_cmems_flag, enable_gofs_flag, latitude_qc_flag, longitude_qc_flag, glider_data_flag)

if __name__ == "__main__":
    main(power=0.9, path="local", target_datetime=dt.datetime.now(dt.timezone.utc), single_datetime=True, enable_rtofs=True, enable_cmems=True, enable_gofs=True, latitude_qc='21', longitude_qc='-86')

# =========================
# ///// END OF SCRIPT \\\\\
# =========================