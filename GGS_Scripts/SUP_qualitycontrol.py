# =========================
# X - IMPORTS
# =========================

### /// QUALITY CONTROL CHECKS ///
import matplotlib.pyplot as plt
import numpy as np
import os
from SUB_functions import calculate_nearpoint

# =========================
# QUALITY CONTROL PLOTS
# =========================

### FUNCTION:
def qc_currents_comparison(config, directory, subset_data, original_data, latitude, longitude):
    
    '''
    Compare 'u' and 'v' profiles over depth from subset and original datasets.

    Args:
    - subset_data (xarray.Dataset): Processed dataset.
    - original_data (xarray.Dataset): Original dataset.
    - latitude (float): Latitude of the point of interest.
    - longitude (float): Longitude of the point of interest.

    Returns:
    - None
    '''

    (y_index, x_index), (lat_index, lon_index) = calculate_nearpoint(subset_data, latitude, longitude)
    (y_index_qc, x_index_qc), (lat_index_qc, lon_index_qc) = calculate_nearpoint(original_data, latitude, longitude)

    u_data = subset_data['u'][:, y_index, x_index]
    v_data = subset_data['v'][:, y_index, x_index]
    u_data_qc = original_data['u'][:, y_index_qc, x_index_qc]
    v_data_qc = original_data['v'][:, y_index_qc, x_index_qc]

    magnitude_subset = np.sqrt(u_data**2 + v_data**2)
    magnitude_original = np.sqrt(u_data_qc**2 + v_data_qc**2)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    axes[0].plot(u_data, subset_data['depth'], label='u', color='blue')
    axes[0].plot(v_data, subset_data['depth'], label='v', color='orange')
    axes[0].plot(magnitude_subset, subset_data['depth'], label='Magnitude', color='black', linestyle='--')
    axes[0].set_title('Subset Data')
    axes[0].set_xlabel('Velocity (m/s)')
    axes[0].set_ylabel('Depth (m)')
    axes[0].invert_yaxis()
    axes[0].legend()

    axes[1].plot(u_data_qc, original_data['depth'], label='u', color='blue')
    axes[1].plot(v_data_qc, original_data['depth'], label='v', color='orange')
    axes[1].plot(magnitude_original, original_data['depth'], label='Magnitude', color='black', linestyle='--')
    axes[1].set_title('Original Data')
    axes[1].set_xlabel('Velocity (m/s)')
    axes[1].invert_yaxis()
    axes[1].legend()

    plt.suptitle(f'Quality Control Comparison - Lat: {lat_index:.3f}, Lon: {lon_index:.3f}', fontsize=14)
    plt.tight_layout()

    fig_filename = f"{config['glider_name']}_QualityControl_Comparison.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')

    plt.close(fig)

### FUNCTION:
def qc_currents_profile(config, directory, model_data, currents_data, bin_data, latitude, longitude):
    
    '''
    Plot the model output 'u' and 'v' magnitudes over depth and display the depth-averaged 'u' and 'v' as vertical lines.
    Also, plot the bin averages and original data points.

    Args:
    - config (dict): Configuration dictionary
    - directory (str): Directory to save the quality control plot.
    - model_data (xarray.Dataset): Dataset containing 'u' and 'v' data.
    - currents_data (xarray.Dataset): Dataset containing depth-averaged 'u' and 'v' data.
    - bin_data (xarray.Dataset): Dataset containing bin-averaged 'u' and 'v' data.
    - latitude (float): Latitude of the point of interest.
    - longitude (float): Longitude of the point of interest.

    Returns:
    - None
    '''

    (y_index, x_index), (lat_index, lon_index) = calculate_nearpoint(model_data, latitude, longitude)

    u_data = model_data['u'][:, y_index, x_index]
    v_data = model_data['v'][:, y_index, x_index]
    
    avg_u = currents_data['u_avg'].isel(y=y_index, x=x_index).values
    avg_v = currents_data['v_avg'].isel(y=y_index, x=x_index).values

    bin_u_data = bin_data['bin_avg_u'].isel(y=y_index, x=x_index).values
    bin_v_data = bin_data['bin_avg_v'].isel(y=y_index, x=x_index).values

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    axes[0].plot(u_data, model_data['depth'], label='u', color='black')
    axes[0].scatter(bin_u_data, np.arange(len(bin_u_data)), color='blue', label='Bin Average u', alpha=0.2)
    axes[0].scatter(u_data, model_data['depth'], color='black', label='u Datapoint', alpha=1.0)
    axes[0].axvline(x=avg_u, color='green', linestyle='--', label='Depth Average u')
    axes[0].set_title('U Magnitude Over Depth')
    axes[0].set_xlabel('U Velocity (m/s)')
    axes[0].set_ylabel('Depth (m)')
    axes[0].invert_yaxis()
    axes[0].legend()

    axes[1].plot(v_data, model_data['depth'], label='v', color='black')
    axes[1].scatter(bin_v_data, np.arange(len(bin_v_data)), color='orange', label='Bin Average v', alpha=0.2)
    axes[1].scatter(v_data, model_data['depth'], color='black', label='v Datapoint', alpha=1.0)
    axes[1].axvline(x=avg_v, color='green', linestyle='--', label='Depth Average v')
    axes[1].set_title('V Magnitude Over Depth')
    axes[1].set_xlabel('V Velocity (m/s)')
    axes[1].invert_yaxis()
    axes[1].legend()

    plt.suptitle(f'Quality Control Profiles - Lat: {lat_index:.3f}, Lon: {lon_index:.3f}', fontsize=14)
    plt.tight_layout()

    fig_filename = f"{config['glider_name']}_QualityControl_Profile.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    
    plt.close(fig)

# =========================
#
# latitude = '21.5'
# longitude = '-85.5'
#
# qc_currents_comparison(config, directory, rtofs_data, rtofs_qc, latitude, longitude)
#
# qc_currents_profile(config, directory, rtofs_data, currents_data, bin_data, latitude, longitude)
#
# =========================