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

# =========================
# ROUTE ANALYSIS
# =========================

### IMPORT:
from SUP_route_analysis import route_analysis, route_analysis_output

# =========================
# [RTOFS] DATA PROCESSING
# =========================

### IMPORT:
from MOD_rtofs import RTOFS, model_currents

# =========================
# PLOTS
# =========================

### IMPORT:
from SUP_plots import GGS_plot_currents

# Gulf of Mexico: map_lons=[-95,0], map_lats=[0,50]

# =========================
# X - MAIN
# =========================

### RUN:
EXIT_KEYWORD = "EXIT"
def main():
    
    """
    GGS main function.
    """

    config, waypoints = GGS_config_static() # Manual
    # config, waypoints = GGS_config() # Automatic

    directory = GGS_config_output(config)

    analysis_results = route_analysis(config, waypoints)
    route_analysis_output(config, directory, analysis_results)
    
    RTOFS_class = RTOFS()
    RTOFS_class.rtofs_subset(config, waypoints, subset=True)
    u_avg, v_avg, magnitude = model_currents(RTOFS_class, mask=False)

    GGS_plot_currents(config, waypoints, directory, RTOFS_class, u_avg, v_avg, magnitude, extent='map', map_lons=[-95, 0], map_lats=[0, 50], show_route=False)

if __name__ == "__main__":
    main()

# =========================
# ///// END OF SCRIPT \\\\\
# =========================