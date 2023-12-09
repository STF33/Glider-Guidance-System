# =========================
# FUNCTIONS
# =========================

from SUB_functions import *

# =========================
# GGS CONFIGURATION
# =========================

from SUP_config import *

# =========================
# ROUTE ANALYSIS
# =========================

from SUP_route_analysis import *

# =========================
# [RTOFS] MODEL DATA PROCESSING
# =========================

from MOD_rtofs import *

# =========================
# QUALITY CONTROL CHECKS
# =========================

from SUP_qualitycontrol import *

# =========================
# PLOTS
# =========================

from SUP_plots import *

# =========================
# X - MAIN
# =========================

### RUN:
EXIT_KEYWORD = "EXIT"
def main():
    
    '''
    GGS main function.
    '''

    config, waypoints = GGS_config_static() # Manual
    # config, waypoints = GGS_config() # Automatic
    directory = GGS_config_output(config)

    analysis_results = route_analysis(config, waypoints)
    route_analysis_output(config, directory, analysis_results)
    
    rtofs_data = xr.open_dataset('C:/Users/salfr/OneDrive/Desktop/STF-0/!-GGS/0-Demo/GGS_GoM/GoM_rtofs_1000.0m.nc')
    rtofs_qc = xr.open_dataset('C:/Users/salfr/OneDrive/Desktop/STF-0/!-GGS/0-Demo/GGS_GoM/GoM_rtofs_1000.0m.nc')
    
    calculated_data = xr.open_dataset('C:/Users/salfr/OneDrive/Desktop/STF-0/!-GGS/0-Demo/GGS_GoM/GoM_calculated_data.nc')
    bin_data = xr.open_dataset('C:/Users/salfr/OneDrive/Desktop/STF-0/!-GGS/0-Demo/GGS_GoM/GoM_bin_data.nc')

    qc_latitude = '21.5'
    qc_longitude = '-85.5'
    qc_uv_profile(config, directory, rtofs_qc, calculated_data, bin_data, qc_latitude, qc_longitude)

    GGS_plot_currents(config, directory, waypoints, rtofs_data, calculated_data, qc_latitude, qc_longitude, extent='data', map_lons=[0, 0], map_lats=[0, 0], show_route=False, show_qc=False)

if __name__ == "__main__":
    main()

# =========================
# ///// END OF SCRIPT \\\\\
# =========================