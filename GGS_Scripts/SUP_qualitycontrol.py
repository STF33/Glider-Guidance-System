# =========================
# X - IMPORTS
# =========================

### /// QUALITY CONTROL CHECKS ///
import cmocean.cm as cmo
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr
from MOD_rtofs import compute_currents
from SUB_functions import calculate_nearpoint

# =========================
# QUALITY CONTROL PLOTS
# =========================

### FUNCTION:
def qc_currents_comparison(directory, rtofs_data, rtofs_qc, latitude, longitude):
    
    '''
    Compare 'u' and 'v' profiles over depth from processed and original datasets.

    Args:
    - rtofs_data (xarray.Dataset): Processed dataset.
    - rtofs_qc (xarray.Dataset): Original dataset.
    - latitude (float): Latitude of the point of interest.
    - longitude (float): Longitude of the point of interest.

    Returns:
    - None
    '''

    (y_index, x_index), (lat_index, lon_index) = calculate_nearpoint(rtofs_data, latitude, longitude)
    (y_index_qc, x_index_qc), (lat_index_qc, lon_index_qc) = calculate_nearpoint(rtofs_qc, latitude, longitude)

    u_data = rtofs_data['u'][:, y_index, x_index]
    v_data = rtofs_data['v'][:, y_index, x_index]
    u_data_qc = rtofs_qc['u'][:, y_index_qc, x_index_qc]
    v_data_qc = rtofs_qc['v'][:, y_index_qc, x_index_qc]

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    axes[0].plot(u_data, rtofs_data['depth'], label='u')
    axes[0].plot(v_data, rtofs_data['depth'], label='v')
    axes[0].set_title('Subset Data')
    axes[0].set_xlabel('Component Velocity Magnitude(m/s)')
    axes[0].set_ylabel('Depth (m)')
    axes[0].invert_yaxis()
    axes[0].legend()

    axes[1].plot(u_data_qc, rtofs_qc['depth'], label='u')
    axes[1].plot(v_data_qc, rtofs_qc['depth'], label='v')
    axes[1].set_title('Original Data')
    axes[1].set_xlabel('Component Velocity Magnitude(m/s)')
    axes[1].invert_yaxis()
    axes[1].legend()

    plt.suptitle(f'Quality Control Comparison - Lat: {lat_index:.3f}, Lon: {lon_index:.3f}', fontsize=14)
    plt.tight_layout()

    fig_filename = "GGS_QualityControl_Comparison.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')

### FUNCTION:
def qc_currents_profile(directory, rtofs_data, latitude, longitude):
    
    '''
    Plot 'u' and 'v' magnitudes over depth and thickness-weighted depth-averaged 'u' and 'v' as vertical lines.

    Args:
    - rtofs_data (xarray.Dataset): Dataset containing 'u' and 'v' data.
    - latitude (float): Latitude of the point of interest.
    - longitude (float): Longitude of the point of interest.

    Returns:
    - None
    '''

    (y_index, x_index), (lat_index, lon_index) = calculate_nearpoint(rtofs_data, latitude, longitude)

    u_data = rtofs_data['u'][:, y_index, x_index]
    v_data = rtofs_data['v'][:, y_index, x_index]
    
    u_avg, v_avg, magnitude, currents_data = compute_currents(rtofs_data)

    avg_u = u_avg.isel(y=y_index, x=x_index).values
    avg_v = v_avg.isel(y=y_index, x=x_index).values

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    axes[0].plot(u_data, rtofs_data['depth'], label='u')
    axes[0].axvline(x=avg_u, color='r', linestyle='--', label='Weighted Depth Average u')
    axes[0].set_title('U Magnitude Over Depth')
    axes[0].set_xlabel('U Velocity (m/s)')
    axes[0].set_ylabel('Depth (m)')
    axes[0].invert_yaxis()
    axes[0].legend()

    axes[1].plot(v_data, rtofs_data['depth'], label='V')
    axes[1].axvline(x=avg_v, color='r', linestyle='--', label='Weighted Depth Average u')
    axes[1].set_title('V Magnitude Over Depth')
    axes[1].set_xlabel('V Velocity (m/s)')
    axes[1].invert_yaxis()
    axes[1].legend()

    plt.suptitle(f'Quality Control Profiles - Lat: {lat_index:.3f}, Lon: {lon_index:.3f}', fontsize=14)
    plt.tight_layout()

    fig_filename = "GGS_QualityControl_Profile.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')

# =========================
# X - MAIN
# =========================
# latitude = 20.5 # Input latitude
# longitude = -86.0 # Input longitude
#
# qc_currents_comparison(directory, rtofs_data, rtofs_qc, latitude, longitude)
#
# qc_currents_profile(directory, rtofs_data, latitude, longitude)
# =========================