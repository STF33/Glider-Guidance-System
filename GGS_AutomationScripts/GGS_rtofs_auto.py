# =========================
# FUNCTIONS
# =========================

from A_functions import *

# =========================
# GGS CONFIGURATION
# =========================

from A_config import *

# =========================
# [RTOFS] MODEL DATA PROCESSING
# =========================

from A_rtofs import *

# =========================
# INTERPOLATION CALCULATION
# =========================

from A_interpolation import *

# =========================
# QUALITY CONTROL CHECKS
# =========================

from A_quality_control import *

# =========================
# PLOTS
# =========================

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
    
    formatted_datetime = datetime_format(datetime_index)
    sub_directory = os.path.join(root_directory, formatted_datetime)

    if not os.path.exists(sub_directory):
        os.makedirs(sub_directory, exist_ok=True)
        
        rtofs = RTOFS(datetime_index)
        rtofs.rtofs_subset(config, buffer=0, subset=True)
        rtofs_data = rtofs.data
        rtofs_qc = rtofs.rtofs_qc
        # rtofs.rtofs_save(config, sub_directory)

        depth_average_data, bin_average_data = interpolation_model(config, sub_directory, rtofs_data, save_depth_average=True, save_bin_average=True)

        qc_latitude = '21.100'
        qc_longitude = '-86.25'
        GGSS_plot_qc(config, sub_directory, rtofs_qc, depth_average_data, bin_average_data, qc_latitude, qc_longitude, threshold=0.5)
        GGS_plot_currents(config, sub_directory, rtofs_data, depth_average_data, qc_latitude, qc_longitude, density=2, show_route=False, show_qc=False)
        GGS_plot_threshold(config, sub_directory, rtofs_data, depth_average_data, qc_latitude, qc_longitude, density=2, mag1=0.0, mag2=0.2, mag3=0.3, mag4=0.4, mag5=0.5, show_route=False, show_qc=False)
    else:
        print(f'Skipping processing for datetime index: {datetime_index}')

### MAIN:
def main(date_indices=None):

    '''
    GGS main function.

    Args:
    - date_indices (list): List of datetime indices to process.

    Returns:
    - None
    '''

    config = GGS_config_static(date=dt.datetime.now(timezone.utc))
    root_directory = GGS_config_output(config, path="default")

    if date_indices is not None:
        date_list = [config['date_list'][i] for i in date_indices]
    else:
        date_list = config['date_list']

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(GGS_process_date, datetime_index, config, root_directory) for datetime_index in date_list]

        for future in as_completed(futures):
            print(future.result())

if __name__ == "__main__":
    main(date_indices=None)

# =========================
# ///// END OF SCRIPT \\\\\
# =========================