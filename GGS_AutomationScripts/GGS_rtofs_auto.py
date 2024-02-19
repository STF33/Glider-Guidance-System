# =========================
# IMPORTS
# =========================

from A_functions import *

from A_config import *

from A_models import *

from A_interpolation import *

from A_plots import *

# =========================
# MAIN
# =========================

from concurrent.futures import ProcessPoolExecutor

### PARALELLIZE:
def GGS_process_date(datetime_index, config, root_directory):
    
    '''
    Process a single datetime index.

    Args:
    - datetime_index (int): Index of the date to process in date_list.
    - config (dict): Configuration dictionary.
    - root_directory (str): Root directory to save output to.

    Returns:
    - None
    '''
    
    # INITIALIZATION
    sub_directory_plots = os.path.join(root_directory, "plots", ''.join(datetime_index[:10].split('-')))
    sub_directory_data = os.path.join(root_directory, "data", ''.join(datetime_index[:10].split('-')))
    
    os.makedirs(sub_directory_plots, exist_ok=True)
    os.makedirs(sub_directory_data, exist_ok=True)
    
    tstr = pd.to_datetime(datetime_index).strftime('%Y%m%dT%HZ')
    sname = os.path.join(sub_directory_data, f"{config['glider_name']}_DepthAverageData_{tstr}.nc")

    if os.path.exists(sname):
        print(f'{sname} already exists. Skipping processing for datetime index: {datetime_index}')
        return
    
    # [RTOFS] MODEL DATA PROCESSING
    rtofs = RTOFS()
    rtofs.rtofs_load(datetime_index)
    rtofs.rtofs_subset(config, subset=True) # [subset] OPTIONS: True = Config Extent, False = Full Extent
    rtofs.rtofs_save(config, sub_directory_data, save_data=False, save_qc=False) # [save] OPTIONS: True = Save Data File, False = Do Not Save Data File
    rtofs_data = rtofs.data
    rtofs_qc = rtofs.qc
    
    depth_average_rtofs, bin_average_rtofs = interpolate_rtofs(config, sub_directory_data, rtofs_data, chunk=False, save_depth_average=True, save_bin_average=False) # [save] OPTIONS: True = Save Data File, False = Do Not Save Data File

    # [CMEMS] MODEL DATA PROCESSING
    cmems = CMEMS(username='sfricano1', password='GlobalGliders1')
    cmems.cmems_load(config, datetime_index)
    cmems.cmems_subset(config, subset=True) # [subset] OPTIONS: True = Depth Subset, False = Full Depth Range
    cmems.cmems_save(config, sub_directory_data, save_data=False, save_qc=False) # [save] OPTIONS: True = Save Data File, False = Do Not Save Data File
    cmems_data = cmems.data
    cmems_qc = cmems.qc
    
    depth_average_cmems, bin_average_cmems = interpolate_cmems(config, sub_directory_data, cmems_data, chunk=False, save_depth_average=True, save_bin_average=False) # [save] OPTIONS: True = Save Data File, False = Do Not Save Data File

    # PLOTTING
    latitude_qc = '20.30'
    longitude_qc = '-86.50'
    GGS_plot_magnitude(config, sub_directory_plots, depth_average_rtofs, latitude_qc, longitude_qc, model_name="rtofs", density=2, show_gliders=True, show_route=False, show_qc=False, manual_extent=None)
    GGS_plot_threshold(config, sub_directory_plots, depth_average_rtofs, latitude_qc, longitude_qc, model_name="rtofs", density=2, mag1=0.0, mag2=0.2, mag3=0.3, mag4=0.4, mag5=0.5, show_gliders=True, show_route=False, show_qc=False, manual_extent=None)
    GGS_plot_magnitude(config, sub_directory_plots, depth_average_cmems, latitude_qc, longitude_qc, model_name="cmems", density=2, show_gliders=True, show_route=False, show_qc=False, manual_extent=None)
    GGS_plot_threshold(config, sub_directory_plots, depth_average_cmems, latitude_qc, longitude_qc, model_name="cmems", density=2, mag1=0.0, mag2=0.2, mag3=0.3, mag4=0.4, mag5=0.5, show_gliders=True, show_route=False, show_qc=False, manual_extent=None)
    GGS_plot_magnitude_dual(config, sub_directory_plots, depth_average_rtofs, depth_average_cmems, latitude_qc, longitude_qc, model_name_1="rtofs", model_name_2="cmems", density=2, show_gliders=True, show_route=False, show_qc=False, manual_extent=[-88, -84, 18, 22])

### MAIN:
def main(root="local", single=False):

    '''
    GGS main function.

    Args:
    - root (str): Root directory to save output to. Options: 'local' or 'rucool'.
    - single (bool): Process a single date (True or False).

    Returns:
    - None
    '''

    config = GGS_config_static(date=dt.datetime.now(timezone.utc))
    
    if root == "local":
        root_directory = GGS_config_output(config, path="default")
    elif root == "rucool":
        root_directory = GGS_config_output(config, path="/www/web/rucool/hurricane/model_comparisons/maps/yucatan/GGS_Yucatan")
    else:
        raise ValueError("Invalid root directory. Please set root to either 'local' or 'rucool'.")

    if single:
        date_list = [dt.datetime.now(timezone.utc).replace(hour=0, minute=0, second=0, microsecond=0).strftime('%Y-%m-%dT%H:%M:%SZ')]
    else:
        date_start = dt.datetime.now(timezone.utc).replace(hour=0, minute=0, second=0, microsecond=0)
        date_end = date_start + dt.timedelta(days=1)
        freq = '6H'
        date_range = pd.date_range(date_start, date_end, freq=freq).tz_localize(None)
        date_list = [date.strftime('%Y-%m-%dT%H:%M:%SZ') for date in date_range]

    with ProcessPoolExecutor(max_workers=6) as executor:
        executor.map(GGS_process_date, date_list, [config] * len(date_list), [root_directory] * len(date_list))

if __name__ == "__main__":
    main(root="local", single=False)

# =========================
# ///// END OF SCRIPT \\\\\
# =========================