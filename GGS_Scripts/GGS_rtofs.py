# =========================
# X - Imports
# =========================

### /// FUNCTIONS ///
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import numpy as np

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
from MOD_rtofs import compute_currents
from SUB_functions import calculate_nearpoint

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
from SUP_config import *

### NOTES:
# N/A

# =========================
# ROUTE ANALYSIS
# =========================

### IMPORT:
from SUP_route_analysis import *

### NOTES:
# N/A

# =========================
# [RTOFS] DATA PROCESSING
# =========================

### IMPORT:
from MOD_rtofs import *

### NOTES:
# Find layer depth derivation code from Dimitry hycom python files
# Function for calculating the center thickness u and v value by taking u1 and u2 at z1 and z2 and calculating the average between them

# =========================
# QUALITY CONTROL CHECKS
# =========================

### IMPORT:
from SUP_qualitycontrol import *

### NOTES:
# N/A

# =========================
# PLOTS
# =========================

### IMPORT:
from SUP_plots import *

### NOTES:
# N/A

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

    rtofs = RTOFS()
    rtofs.rtofs_subset(config, waypoints, subset=True)
    
    rtofs_data = rtofs.data
    rtofs_qc = rtofs.rtofs_qc
    
    u_avg, v_avg, magnitude, currents_data = compute_currents(rtofs_data)
    
    qc_currents_comparison(directory, rtofs_data, rtofs_qc, latitude=20.5, longitude=-86.0)
    qc_currents_profile(directory, rtofs_data, latitude=20.5, longitude=-86.0)

    GGS_plot_currents(config, waypoints, directory, rtofs_data, u_avg, v_avg, magnitude, extent='data', map_lons=[0, 0], map_lats=[0, 0], show_route=False)

if __name__ == "__main__":
    main()

# =========================
# ///// END OF SCRIPT \\\\\
# =========================