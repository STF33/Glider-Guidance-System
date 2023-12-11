# =========================
# X - IMPORTS
# =========================

import matplotlib.pyplot as plt
import numpy as np
import os
from SUB_functions import calculate_gridpoint

# =========================
# QUALITY CONTROL PLOTS
# =========================

### FUNCTION:
def qc_uv_profile(config, directory, model_data, calculated_data, bin_data, qc_latitude, qc_longitude):
    
    '''
    Intake a latitude and longitude to provide Quality Control data profiles at.
    Scatter the ocean model 'u' and 'v' magnitudes over depth, Scatter the1-meter bin averages, and display the depth-averaged 'u' and 'v' component velocities as vertical lines.
    Also, plot the bin averages and original data points.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - model_data (xarray.Dataset): Model dataset containing 'u' and 'v' data.
    - currents_data (xarray.Dataset): Dataset containing depth-averaged 'u' and 'v' data.
    - bin_data (xarray.Dataset): Dataset containing bin-averaged 'u' and 'v' data.
    - latitude (float): Latitude of the point of interest.
    - longitude (float): Longitude of the point of interest.

    Returns:
    - None
    '''
    
    (y_model_index, x_model_index), (lat_model_index, lon_model_index) = calculate_gridpoint(model_data, qc_latitude, qc_longitude)
    (y_calc_index, x_calc_index), (lat_calc_index, lon_calc_index) = calculate_gridpoint(calculated_data, qc_latitude, qc_longitude)
    (y_bin_index, x_bin_index), (lat_bin_index, lon_bin_index) = calculate_gridpoint(bin_data, qc_latitude, qc_longitude)
    
    u_data = model_data['u'].isel(y=y_model_index, x=x_model_index).values
    v_data = model_data['v'].isel(y=y_model_index, x=x_model_index).values
    
    avg_u = calculated_data['u_avg'].isel(y=y_calc_index, x=x_calc_index).values
    avg_v = calculated_data['v_avg'].isel(y=y_calc_index, x=x_calc_index).values

    bin_u_data = bin_data['bin_avg_u'].isel(y=y_bin_index, x=x_bin_index).values
    bin_v_data = bin_data['bin_avg_v'].isel(y=y_bin_index, x=x_bin_index).values
    
    fig, axes = plt.subplots(1, 2, figsize=(6, 10))

    axes[0].scatter(u_data, model_data['depth'], marker='x', color='black', s=100, label='Model Datapoint', alpha=1.0, zorder=3)
    axes[0].scatter(bin_u_data, np.arange(len(bin_u_data)), color='cyan', label='1m Interpolation', alpha=1.0, zorder=2)
    axes[0].axvline(x=avg_u, color='darkcyan', linestyle='--', linewidth=2, label='Depth Average u', zorder=1)
    axes[0].set_xlabel('u Velocity (m/s)', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
    axes[0].invert_yaxis()
    axes[0].legend(loc='lower center', facecolor='lightgrey')
    
    axes[1].scatter(v_data, model_data['depth'], marker='x', color='black', s=100, label='Model Datapoint', alpha=1.0, zorder=3)
    axes[1].scatter(bin_v_data, np.arange(len(bin_v_data)), color='orange', label='1m Interpolation', alpha=1.0, zorder=2)
    axes[1].axvline(x=avg_v, color='darkorange', linestyle='--', linewidth=2, label='Depth Average v', zorder=1)
    axes[1].set_xlabel('v Velocity (m/s)', fontsize=12, fontweight='bold')
    axes[1].invert_yaxis()
    axes[1].legend(loc='lower center', facecolor='lightgrey')
    
    plt.suptitle(f'Quality Control Profiles - Lat: {lat_model_index:.3f}, Lon: {lon_model_index:.3f}', fontsize=14)
    plt.subplots_adjust(top=0.90)
    plt.tight_layout()

    fig_filename = f"{config['glider_name']}_QualityControl_uvProfile.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    
    plt.close(fig)
