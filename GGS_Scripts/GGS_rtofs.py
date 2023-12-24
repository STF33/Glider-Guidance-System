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
EXIT_KEYWORD = "EXIT"
def main():
    
    '''
    GGS main function.
    '''

    config, GPS_coords = GGS_config_static() # Manual
    # config, waypoints = GGS_config() # Automatic
    directory = GGS_config_output(config)

    rtofs = RTOFS()
    rtofs.rtofs_subset(config, GPS_coords, subset=True)
    rtofs_data = rtofs.data
    rtofs_qc = rtofs.rtofs_qc
    rtofs.rtofs_save(config, directory)
    calculated_data, bin_data = interp_depth_average(config, directory, rtofs_data)

    qc_latitude = '21.100'
    qc_longitude = '-86.25'
    qc_uv_profile(config, directory, rtofs_qc, calculated_data, bin_data, qc_latitude, qc_longitude)

    GGS_plot_currents(config, directory, GPS_coords, rtofs_data, calculated_data, qc_latitude, qc_longitude, extent='data', map_lons=[0, 0], map_lats=[0, 0], show_route=False, show_qc=False)
    GGS_plot_threshold(config, directory, GPS_coords, rtofs_data, calculated_data, qc_latitude, qc_longitude, mag1=0.0, mag2=0.2, mag3=0.3, mag4=0.4, mag5=0.5, extent='data', map_lons=[0, 0], map_lats=[0, 0], show_route=False, show_qc=False)

if __name__ == "__main__":
    main()

# =========================
# ///// END OF SCRIPT \\\\\
# =========================