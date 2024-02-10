# =========================
# X - IMPORTS
# =========================

# Functions
from A_functions import *

# Configuration
from A_config import *

# [RTOFS] Model Data Processing
from A_rtofs import *

# Interpolation
from A_interpolation import *

# Plotting
from A_plots import *

# =========================
# X - MAIN
# =========================

from concurrent.futures import ProcessPoolExecutor, as_completed

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
    
    sub_directory_plots = os.path.join(root_directory, "plots", ''.join(datetime_index[:10].split('-')))
    sub_directory_data = os.path.join(root_directory, "data", ''.join(datetime_index[:10].split('-')))
    
    os.makedirs(sub_directory_plots, exist_ok=True)
    os.makedirs(sub_directory_data, exist_ok=True)
    
    tstr = pd.to_datetime(datetime_index).strftime('%Y%m%dT%HZ')
    sname = os.path.join(sub_directory_data, f"{config['glider_name']}_DepthAverageData_{tstr}.nc")

    if os.path.exists(sname):
        print(f'{sname} already exists. Skipping processing for datetime index: {datetime_index}')
        return
    
    rtofs = RTOFS()
    rtofs.rtofs_load(datetime_index)
    rtofs.rtofs_subset(config, subset=True)
    rtofs.rtofs_save(config, sub_directory_data, save_data=True, save_qc=True)
    rtofs_data = rtofs.data
    rtofs_qc = rtofs.qc

    depth_average_data, bin_average_data = interpolation_model(config, sub_directory_data, rtofs_data, save_depth_average=True, save_bin_average=True)

    latitude_qc = '21.100'
    longitude_qc = '-86.25'
    GGS_plot_profiles(config, sub_directory_plots, rtofs_qc, depth_average_data, bin_average_data, latitude_qc, longitude_qc, threshold=0.5)
    GGS_plot_currents(config, sub_directory_plots, rtofs_data, depth_average_data, latitude_qc, longitude_qc, density=2, show_gliders=True, show_route=False, show_qc=False)
    GGS_plot_threshold(config, sub_directory_plots, rtofs_data, depth_average_data, latitude_qc, longitude_qc, density=2, mag1=0.0, mag2=0.2, mag3=0.3, mag4=0.4, mag5=0.5, show_gliders=True, show_route=False, show_qc=False)
    
### MAIN:
def main(root="local"):

    '''
    GGS main function.

    Args:
    - date_indices (list): List of datetime indices to process.

    Returns:
    - None
    '''

    config = GGS_config_static(date=dt.datetime.now(timezone.utc))
    
    if root == "local":
        root_directory = GGS_config_output(config, path="default")
    elif root == "rucool":
        root_directory = GGS_config_output(config, path="/www/web/rucool/hurricane/model_comparisons/maps/yucatan/GGS_Yucatan")
    
    date_start = dt.datetime.now(timezone.utc).replace(hour=0, minute=0, second=0, microsecond=0)
    date_end = date_start + dt.timedelta(days=1)
    freq = '6H'
    date_range = pd.date_range(date_start, date_end, freq=freq).tz_localize(None)
    date_list = [date.strftime('%Y-%m-%dT%H:%M:%SZ') for date in date_range]

    with ProcessPoolExecutor(max_workers=6) as executor:
        executor.map(GGS_process_date, date_list, [config] * len(date_list), [root_directory] * len(date_list))

if __name__ == "__main__":
    main()

# =========================
# ///// END OF SCRIPT \\\\\
# =========================