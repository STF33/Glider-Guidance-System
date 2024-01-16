# =========================
# FUNCTIONS
# =========================

from X_functions import *

# =========================
# GGS CONFIGURATION
# =========================

from X_config import *

# =========================
# [RTOFS] MODEL DATA PROCESSING
# =========================

from X_rtofs import *

# =========================
# QUALITY CONTROL CHECKS
# =========================

from X_quality_control import *

# =========================
# PLOTS
# =========================

from X_plots import *

# =========================
# X - MAIN
# =========================

### RUN:
# EXIT_KEYWORD = "EXIT"
def main(date_indices=None):
    
    '''
    GGS main function.

    Args:
    - date_indices (list of int): Indices of the dates to process in date_list.
    '''

    config = GGS_config_static(date=dt.datetime.now(timezone.utc))
    root_directory = GGS_config_output(config)

    if date_indices is not None:
        date_list = [config['date_list'][i] for i in date_indices]
    else:
        date_list = config['date_list']

    for datetime_index in date_list:

        formatted_datetime = datetime_format(datetime_index)
        sub_directory = os.path.join(root_directory, formatted_datetime)
        os.makedirs(sub_directory, exist_ok=True)

        rtofs = RTOFS(datetime_index)
        rtofs.rtofs_subset(config, subset=True)
        rtofs_data = rtofs.data
        rtofs_qc = rtofs.rtofs_qc
        rtofs.rtofs_save(config, sub_directory)

        depth_average_data, bin_average_data = interp_depth_average(config, sub_directory, rtofs_data)
        
        # depth_average_data = xr.open_dataset('C:/Users/sal_f/OneDrive/Desktop/STF-0/!-GGS/0-Demo/Yucatan_DepthAverageData.nc')
        # bin_average_data = xr.open_dataset('C:/Users/sal_f/OneDrive/Desktop/STF-0/!-GGS/0-Demo/Yucatan_BinAverageData.nc')

        qc_latitude = '21.100'
        qc_longitude = '-86.25'
        qc_uv_profile(config, sub_directory, rtofs_qc, depth_average_data, bin_average_data, qc_latitude, qc_longitude)

        GGS_plot_currents(config, sub_directory, rtofs_data, depth_average_data, qc_latitude, qc_longitude, show_route=False, show_qc=False)
        GGS_plot_threshold(config, sub_directory, rtofs_data, depth_average_data, qc_latitude, qc_longitude, mag1=0.0, mag2=0.2, mag3=0.3, mag4=0.4, mag5=0.5, show_route=False, show_qc=False)

if __name__ == "__main__":
    main(date_indices=[1])

# =========================
# ///// END OF SCRIPT \\\\\
# =========================