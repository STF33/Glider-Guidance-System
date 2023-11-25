# =========================
# X - Imports
# =========================

### /// CONFIGURATION ///
import os
import pandas as pd
from SUB_functions import check_abort, check_float, check_coordinate

### /// ROUTE ANALYSIS ///
import os
from SUB_functions import calculate_distance, calculate_heading

### /// RTOFS ///
import numpy as np
import xarray as xr

### /// QUALITY CONTROL CHECKS ///
import cmocean.cm as cmo
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr

### /// PLOTS ///
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import plotly.graph_objects as go
import plotly.io as pio
from SUB_functions import set_ticks
from SUP_route_analysis import route_analysis_output

# =========================
# GGS CONFIGURATION
# =========================

### IMPORT:
from SUP_config import GGS_config_static, GGS_config_import, GGS_config_new, GGS_config, GGS_config_output

### NOTES:
# N/A

# =========================
# ROUTE ANALYSIS
# =========================

### IMPORT:
from SUP_route_analysis import route_analysis, route_analysis_output

### NOTES:
# N/A

# =========================
# [RTOFS] DATA PROCESSING
# =========================

### IMPORT:
from MOD_rtofs import RTOFS, model_currents

### NOTES:
# N/A

# =========================
# QUALITY CONTROL CHECKS
# =========================

### IMPORT:
from SUP_qc_checks import qc_selection, qc_plots

### NOTES:
# N/A

# =========================
# PLOTS
# =========================

### IMPORT:
from SUP_plots import GGS_plot_currents

### NOTES:
# Gulf of Mexico: map_lons=[-95,0], map_lats=[0,50]

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
    
    rtofs_data = RTOFS()
    rtofs_data.rtofs_subset(config, waypoints, subset=True)
    rtofs_qc = rtofs_data.rtofs_qc
    u_avg, v_avg, magnitude = model_currents(rtofs_data, mask=False)

    qc_coordinates = qc_selection(rtofs_qc)
    qc_plots(directory, rtofs_qc, qc_coordinates, u_avg, v_avg)

    GGS_plot_currents(config, waypoints, directory, rtofs_data, u_avg, v_avg, magnitude, extent='data', map_lons=[-95, 0], map_lats=[0, 50], show_route=False)

if __name__ == "__main__":
    main()

# =========================
# ///// END OF SCRIPT \\\\\
# =========================