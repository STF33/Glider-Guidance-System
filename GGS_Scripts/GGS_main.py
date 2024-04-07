# =========================
# IMPORTS
# =========================

from X_functions import *
from X_config import *
from X_models import *
from X_interpolation import *
from X_products import *

from concurrent.futures import ProcessPoolExecutor

# =========================
# MAIN
# =========================

### REPROCESSOR:
def GGS_reprocessor(task):

    '''
    Reprocess all datafiles in the 'reprocess' directory.

    Args:
    - task (dict): A dictionary containing all necessary parameters for processing.

    Returns:
    - None
    '''

    current_directory = os.path.dirname(__file__)
    reprocess_path = os.path.join(current_directory, "data/reprocess")

    model_datasets = []
    model_files = glob.glob(os.path.join(reprocess_path, '*.nc'))
    for model_file in model_files:
        depth_average_dataset = xr.open_dataset(model_file)
        model_name = depth_average_dataset.attrs['model_name']
        if model_name == 'RTOFS':
            rtofs_datasets = (None, depth_average_dataset, None)
            model_datasets.append(rtofs_datasets)
        elif model_name == 'CMEMS':
            cmems_datasets = (None, depth_average_dataset, None)
            model_datasets.append(cmems_datasets)
        elif model_name == 'GOFS':
            gofs_datasets = (None, depth_average_dataset, None)
            model_datasets.append(gofs_datasets)

    datetime_index = model_datasets[0][1].attrs['model_datetime']
    config_flag = task['config_flag']
    root_directory_flag = task['root_directory_flag']
    glider_data_flag = task['glider_data_flag']

    create_magnitude_plot_flag = config_flag['PRODUCT']['create_magnitude_plot']
    create_threshold_plot_flag = config_flag['PRODUCT']['create_threshold_plot']
    create_advantage_plot_flag = config_flag['PRODUCT']['create_advantage_plot']
    create_profiles_plot_flag = config_flag['PRODUCT']['create_profiles_plot']
    create_gpkg_file_flag = config_flag['PRODUCT']['create_gpkg_file']
    latitude_qc_flag = config_flag['PRODUCT']['latitude_qc']
    longitude_qc_flag = config_flag['PRODUCT']['longitude_qc']
    density_flag = config_flag['PRODUCT']['density']
    mag1_flag = config_flag['PRODUCT']['mag1']
    mag2_flag = config_flag['PRODUCT']['mag2']
    mag3_flag = config_flag['PRODUCT']['mag3']
    mag4_flag = config_flag['PRODUCT']['mag4']
    mag5_flag = config_flag['PRODUCT']['mag5']
    tolerance_flag = config_flag['PRODUCT']['tolerance']
    show_route_flag = config_flag['PRODUCT']['show_route']
    show_eez_flag = config_flag['PRODUCT']['show_eez']
    show_qc_flag = config_flag['PRODUCT']['show_qc']
    manual_extent_flag = config_flag['PRODUCT']['manual_extent']
    
    sub_directory_plots = os.path.join(root_directory_flag, "REPROCESSED", "plots", ''.join(datetime_index[:10].split('-')))
    os.makedirs(sub_directory_plots, exist_ok=True)
    sub_directory_data = os.path.join(root_directory_flag, "REPROCESSED", "data", ''.join(datetime_index[:10].split('-')))
    os.makedirs(sub_directory_data, exist_ok=True)
    
    check_datetime = pd.to_datetime(datetime_index).strftime('%Y%m%dT%HZ')
    check_pattern = os.path.join(sub_directory_data, f"*_DepthAverageData_{check_datetime}.nc")
    check_files = glob.glob(check_pattern)
    if check_files:
        print(f"Datetime {datetime_index} already processed: {check_files[0]}, skipping task.")
    else:
        print(f"Datetime {datetime_index} unprocessed, proceeding with task.")

    if create_magnitude_plot_flag:
        GGS_plot_magnitude(
            config_flag,
            sub_directory_plots,
            datetime_index,
            model_datasets,
            latitude_qc=latitude_qc_flag, longitude_qc=longitude_qc_flag,
            density=density_flag,
            gliders=glider_data_flag,
            show_route=show_route_flag, show_eez=show_eez_flag, show_qc=show_qc_flag,
            manual_extent=manual_extent_flag
        )
    if create_threshold_plot_flag:
        GGS_plot_threshold(
            config_flag,
            sub_directory_plots,
            datetime_index,
            model_datasets,
            latitude_qc=latitude_qc_flag, longitude_qc=longitude_qc_flag,
            density=density_flag,
            mag1=mag1_flag, mag2=mag2_flag, mag3=mag3_flag, mag4=mag4_flag, mag5=mag5_flag,
            gliders=glider_data_flag,
            show_route=show_route_flag, show_eez=show_eez_flag, show_qc=show_qc_flag,
            manual_extent=manual_extent_flag
        )
    if create_advantage_plot_flag:
        GGS_plot_advantage(
            config_flag,
            sub_directory_plots,
            datetime_index,
            model_datasets,
            latitude_qc=latitude_qc_flag, longitude_qc=longitude_qc_flag,
            density=density_flag,
            tolerance=tolerance_flag,
            mag1=mag1_flag, mag2=mag2_flag, mag3=mag3_flag, mag4=mag4_flag, mag5=mag5_flag,
            gliders=glider_data_flag,
            show_route=show_route_flag, show_eez=show_eez_flag, show_qc=show_qc_flag,
            manual_extent=manual_extent_flag
        )
    if create_profiles_plot_flag:
        GGS_plot_profiles(
            config_flag,
            sub_directory_plots,
            datetime_index,
            model_datasets,
            latitude_qc=latitude_qc_flag, longitude_qc=longitude_qc_flag,
            threshold=0.5
        )
    if create_gpkg_file_flag:
        GGS_export_gpkg(
            sub_directory_data,
            datetime_index,
            model_datasets
        )

### EXECUTIONER:
def GGS_executioner(task):
    
    '''
    Process a single datetime index.

    Args:
    - task (dict): A dictionary containing all necessary parameters for processing.

    Returns:
    - None
    '''
    
    datetime_index = task['datetime_index']
    config_flag = task['config_flag']
    root_directory_flag = task['root_directory_flag']
    glider_data_flag = task['glider_data_flag']
    
    enable_rtofs_flag = config_flag['MODEL']['enable_rtofs']
    enable_cmems_flag = config_flag['MODEL']['enable_cmems']
    enable_gofs_flag = config_flag['MODEL']['enable_gofs']
    save_model_data_flag = config_flag['MODEL']['save_model_data']
    save_depth_average_flag = config_flag['MODEL']['save_depth_average']
    save_bin_average_flag = config_flag['MODEL']['save_bin_average']
    chunk_flag = config_flag['MODEL']['chunk']

    create_magnitude_plot_flag = config_flag['PRODUCT']['create_magnitude_plot']
    create_threshold_plot_flag = config_flag['PRODUCT']['create_threshold_plot']
    create_advantage_plot_flag = config_flag['PRODUCT']['create_advantage_plot']
    create_profiles_plot_flag = config_flag['PRODUCT']['create_profile_plot']
    create_gpkg_file_flag = config_flag['PRODUCT']['create_gpkg_file']
    latitude_qc_flag = config_flag['PRODUCT']['latitude_qc']
    longitude_qc_flag = config_flag['PRODUCT']['longitude_qc']
    density_flag = config_flag['PRODUCT']['density']
    mag1_flag = config_flag['PRODUCT']['mag1']
    mag2_flag = config_flag['PRODUCT']['mag2']
    mag3_flag = config_flag['PRODUCT']['mag3']
    mag4_flag = config_flag['PRODUCT']['mag4']
    mag5_flag = config_flag['PRODUCT']['mag5']
    tolerance_flag = config_flag['PRODUCT']['tolerance']
    show_route_flag = config_flag['PRODUCT']['show_route']
    show_eez_flag = config_flag['PRODUCT']['show_eez']
    show_qc_flag = config_flag['PRODUCT']['show_qc']
    manual_extent_flag = config_flag['PRODUCT']['manual_extent']
    
    sub_directory_plots = os.path.join(root_directory_flag, "plots", ''.join(datetime_index[:10].split('-')))
    os.makedirs(sub_directory_plots, exist_ok=True)
    sub_directory_data = os.path.join(root_directory_flag, "data", ''.join(datetime_index[:10].split('-')))
    os.makedirs(sub_directory_data, exist_ok=True)
    
    check_datetime = pd.to_datetime(datetime_index).strftime('%Y%m%dT%HZ')
    check_pattern = os.path.join(sub_directory_data, f"*_DepthAverageData_{check_datetime}.nc")
    check_files = glob.glob(check_pattern)
    if check_files:
        print(f"Datetime {datetime_index} already processed: {check_files[0]}, skipping task.")
    else:
        print(f"Datetime {datetime_index} unprocessed, proceeding with task.")
    
    model_datasets = []
    if enable_rtofs_flag:
        try:
            rtofs = RTOFS()
            rtofs.rtofs_load(config_flag, datetime_index)
            rtofs.rtofs_save(config_flag, sub_directory_data, save_data=save_model_data_flag)
            rtofs_model_data = rtofs.data
            
            rtofs_depth_average, rtofs_bin_average = interpolate_rtofs(config_flag, sub_directory_data, rtofs_model_data, chunk=chunk_flag, save_depth_average=save_depth_average_flag, save_bin_average=save_bin_average_flag)

            rtofs_datasets = (rtofs_model_data, rtofs_depth_average, rtofs_bin_average)
            model_datasets.append(rtofs_datasets)
        except Exception as e:
            print(f"Error during RTOFS processing: {e}")
    if enable_cmems_flag:
        try:
            cmems = CMEMS(username='sfricano1', password='GlobalGliders1')
            cmems.cmems_load(config_flag, datetime_index)
            cmems.cmems_save(config_flag, sub_directory_data, save_data=save_model_data_flag)
            cmems_model_data = cmems.data
            
            cmems_depth_average, cmems_bin_average = interpolate_cmems(config_flag, sub_directory_data, cmems_model_data, chunk=chunk_flag, save_depth_average=save_depth_average_flag, save_bin_average=save_bin_average_flag)
            
            cmems_datasets = (cmems_model_data, cmems_depth_average, cmems_bin_average)
            model_datasets.append(cmems_datasets)
        except Exception as e:
            print(f"Error during CMEMS processing: {e}")
    if enable_gofs_flag:
        try:
            gofs = GOFS()
            gofs.gofs_load(config_flag, datetime_index)
            gofs.gofs_save(config_flag, sub_directory_data, save_data=save_model_data_flag)
            gofs_model_data = gofs.data

            gofs_depth_average, gofs_bin_average = interpolate_gofs(config_flag, sub_directory_data, gofs_model_data, chunk=chunk_flag, save_depth_average=save_depth_average_flag, save_bin_average=save_bin_average_flag)

            gofs_datasets = (gofs_model_data, gofs_depth_average, gofs_bin_average)
            model_datasets.append(gofs_datasets)
        except Exception as e:
            print(f"Error during GOFS processing: {e}")
    
    if create_magnitude_plot_flag:
        GGS_plot_magnitude(
            config_flag,
            sub_directory_plots,
            datetime_index,
            model_datasets,
            latitude_qc=latitude_qc_flag, longitude_qc=longitude_qc_flag,
            density=density_flag,
            gliders=glider_data_flag,
            show_route=show_route_flag, show_eez=show_eez_flag, show_qc=show_qc_flag,
            manual_extent=manual_extent_flag
        )
    if create_threshold_plot_flag:
        GGS_plot_threshold(
            config_flag,
            sub_directory_plots,
            datetime_index,
            model_datasets,
            latitude_qc=latitude_qc_flag, longitude_qc=longitude_qc_flag,
            density=density_flag,
            mag1=mag1_flag, mag2=mag2_flag, mag3=mag3_flag, mag4=mag4_flag, mag5=mag5_flag,
            gliders=glider_data_flag,
            show_route=show_route_flag, show_eez=show_eez_flag, show_qc=show_qc_flag,
            manual_extent=manual_extent_flag
        )
    if create_advantage_plot_flag:
        GGS_plot_advantage(
            config_flag,
            sub_directory_plots,
            datetime_index,
            model_datasets,
            latitude_qc=latitude_qc_flag, longitude_qc=longitude_qc_flag,
            density=density_flag,
            tolerance=tolerance_flag,
            mag1=mag1_flag, mag2=mag2_flag, mag3=mag3_flag, mag4=mag4_flag, mag5=mag5_flag,
            gliders=glider_data_flag,
            show_route=show_route_flag, show_eez=show_eez_flag, show_qc=show_qc_flag,
            manual_extent=manual_extent_flag
        )
    if create_profiles_plot_flag:
        GGS_plot_profiles(
            config_flag,
            sub_directory_plots,
            datetime_index,
            model_datasets,
            latitude_qc=latitude_qc_flag, longitude_qc=longitude_qc_flag,
            threshold=0.5
        )
    if create_gpkg_file_flag:
        GGS_export_gpkg(
            sub_directory_data,
            datetime_index,
            model_datasets
        )

### MAIN:
def GGS_main(power=1, path="local", config_name=None):
    
    '''
    GGS main function.

    Args:
    - config_name (str): The name of the config file without the extension.
    - path (str): The path directory to save output to. Options: 'local' or a specific directory path.
    
    Returns:
    - None
    '''

    if config_name is None:
        print("No config file specified. Exiting.")
        return

    config = GGS_config_import(config_name)
    
    target_datetime = config['MISSION'].get('target_date')
    if not target_datetime:
        print("Issue with target datetime. Using current datetime.")
        target_datetime = dt.datetime.now(dt.timezone.utc)

    if path == "local":
        root_directory = GGS_config_process(config, path="default")
    elif path == "rucool":
        root_directory = GGS_config_process(config, path="/www/web/rucool/hurricane/model_comparisons/maps/yucatan")
    else:
        raise ValueError("Invalid root directory.")
    
    if config['MODEL'].get('single_datetime'):
        datetime_list = [target_datetime.replace(hour=0, minute=0, second=0, microsecond=0).strftime('%Y-%m-%dT%H:%M:%SZ')]
    else:
        datetime_start = target_datetime.replace(hour=0, minute=0, second=0, microsecond=0)
        datetime_end = datetime_start + dt.timedelta(days=1)
        datetime_range = pd.date_range(datetime_start, datetime_end, freq='6H').tz_localize(None)
        datetime_list = [datetime.strftime('%Y-%m-%dT%H:%M:%SZ') for datetime in datetime_range]

    glider_dataframes = None
    if config['PRODUCT'].get('show_gliders'):
        min_lat, min_lon = config['MISSION']['extent'][0]
        max_lat, max_lon = config['MISSION']['extent'][1]
        search_extent = [min_lon, max_lon, min_lat, max_lat]
        glider_dataframes = acquire_gliders(
            extent=search_extent,
            target_date=target_datetime,
            date_delta=dt.timedelta(days=1),
            requested_variables=["time", "longitude", "latitude", "profile_id", "depth"],
            print_vars=False,
            target="all",
            request_timeout=5,
            enable_parallel=False
        )
    
    if config['ADVANCED']['reprocess']:
        print(f"\n### !!! ALERT: REPROCESSING MODE ENABLED !!! ###\n")
        task = {
            'config_flag': config,
            'root_directory_flag': root_directory,
            'glider_data_flag': glider_dataframes
        }
        GGS_reprocessor(task)
    else:
        tasks = [{
            'datetime_index': datetime_index,
            'config_flag': config,
            'root_directory_flag': root_directory,
            'glider_data_flag': glider_dataframes
        } for datetime_index in datetime_list]

        num_workers = optimal_workers(power=power)
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            print("Starting parallel processing with the following tasks:")
            for i, task in enumerate(tasks, start=1):
                print(f"Task {i}: {task['datetime_index']}")
            executor.map(GGS_executioner, tasks)

if __name__ == "__main__":
    GGS_main(power=1, path="local", config_name="sentinel1")
