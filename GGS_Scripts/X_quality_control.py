# =========================
# X - IMPORTS
# =========================

import matplotlib.pyplot as plt
import numpy as np
import os
from X_functions import calculate_gridpoint, plot_profile_thresholds, datetime_filename

# =========================
# QUALITY CONTROL PLOTS
# =========================

### FUNCTION:
def GGSS_plot_qc(config, directory, model_data, depth_average_data, bin_average_data, qc_latitude, qc_longitude, threshold=0.5):
    
    '''
    Produce quality control profiles for 'u', 'v', 'magnitude', and 'direction' data at the specified point of interest.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - model_data (xarray.Dataset): Model dataset containing 'u' and 'v' data.
    - depth_average_data (xarray.Dataset): Dataset containing depth-averaged 'u' and 'v' data.
    - bin_average_data (xarray.Dataset): Dataset containing bin-averaged 'u', 'v', 'magnitude', and 'direction' data.
    - qc_latitude (float): Latitude of the point of interest.
    - qc_longitude (float): Longitude of the point of interest.

    Returns:
    - None
    '''

    print("\n")
    print("### CREATING QC PLOT ###")
    print("\n")


    model_datetime = model_data.attrs.get('requested_datetime', 'unknown_datetime')

    qc_latitude = float(qc_latitude)
    qc_longitude = float(qc_longitude)

    (y_model_index, x_model_index), _ = calculate_gridpoint(model_data, qc_latitude, qc_longitude)
    (y_avg_index, x_avg_index), _ = calculate_gridpoint(depth_average_data, qc_latitude, qc_longitude)
    (y_bin_index, x_bin_index), _ = calculate_gridpoint(bin_average_data, qc_latitude, qc_longitude)

    u_model = model_data['u'].isel(y=y_model_index, x=x_model_index).values
    v_model = model_data['v'].isel(y=y_model_index, x=x_model_index).values
    u_depth_avg = depth_average_data['u_depth_avg'].isel(y=y_avg_index, x=x_avg_index).values
    v_depth_avg = depth_average_data['v_depth_avg'].isel(y=y_avg_index, x=x_avg_index).values
    mag_depth_avg = depth_average_data['mag_depth_avg'].isel(y=y_avg_index, x=x_avg_index).values
    u_bin_avg = bin_average_data['u_bin_avg'].isel(y=y_bin_index, x=x_bin_index).values
    v_bin_avg = bin_average_data['v_bin_avg'].isel(y=y_bin_index, x=x_bin_index).values
    mag_bin_avg = bin_average_data['mag_bin_avg'].isel(y=y_bin_index, x=x_bin_index).values
    dir_bin_avg = bin_average_data['dir_bin_avg'].isel(y=y_bin_index, x=x_bin_index).values

    adjusted_direction_data = np.mod(dir_bin_avg + 180, 360) - 180

    fig, axes = plt.subplots(1, 4, figsize=(15, 10))

    axes[0].scatter(u_model, model_data['depth'], marker='x', color='black', s=100, label='Model Datapoint', alpha=1.0, zorder=3)
    axes[0].scatter(u_bin_avg, np.arange(len(u_bin_avg)), color='cyan', label='1m Interpolation', alpha=1.0, zorder=2)
    axes[0].axvline(x=u_depth_avg, color='darkcyan', linestyle='--', linewidth=2, label=f'Depth Avgerage = [{u_depth_avg:.2f}]', zorder=1)
    axes[0].set_xlabel('u Velocity (m/s)', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
    axes[0].invert_yaxis()
    axes[0].grid(color='lightgrey', linestyle='-', linewidth=0.5)

    axes[1].scatter(v_model, model_data['depth'], marker='x', color='black', s=100, label='Model Datapoint', alpha=1.0, zorder=3)
    axes[1].scatter(v_bin_avg, np.arange(len(v_bin_avg)), color='orange', label='1m Interpolation', alpha=1.0, zorder=2)
    axes[1].axvline(x=v_depth_avg, color='darkorange', linestyle='--', linewidth=2, label=f'Depth Avgerage = [{v_depth_avg:.2f}]', zorder=1)
    axes[1].set_xlabel('v Velocity (m/s)', fontsize=12, fontweight='bold')
    axes[1].invert_yaxis()
    axes[1].grid(color='lightgrey', linestyle='-', linewidth=0.5)

    axes[2].scatter(mag_bin_avg, np.arange(len(mag_bin_avg)), color='green', label='1m Interpolation', alpha=1.0, zorder=2)
    axes[2].set_xlabel('Current Magnitude (m/s)', fontsize=12, fontweight='bold')
    axes[2].axvline(x=mag_depth_avg, color='darkgreen', linestyle='--', linewidth=2, label=f'Depth Avgerage = [{mag_depth_avg:.2f}]', zorder=1)
    axes[2].invert_yaxis()
    axes[2].grid(color='lightgrey', linestyle='-', linewidth=0.5)

    axes[3].scatter(adjusted_direction_data, np.arange(len(dir_bin_avg)), color='purple', label='1m Interpolation', alpha=1.0, zorder=2)
    axes[3].set_xlabel('Current Direction (degrees)', fontsize=12, fontweight='bold')
    axes[3].invert_yaxis()
    axes[3].grid(color='lightgrey', linestyle='-', linewidth=0.5)

    shading_colors = ['cyan', 'orange', 'green']
    for ax, data, color in zip(axes[:3], [u_bin_avg, v_bin_avg, mag_bin_avg], shading_colors):
        plot_profile_thresholds(ax, data, threshold, color)
    
    for ax in axes:
        ax.legend(loc='lower center', facecolor='lightgrey', edgecolor='black', framealpha=1.0)

    title_text = f"{config['glider_name']} Mission - Quality Control Profiles - Lat: {qc_latitude:.3f}, Lon: {qc_longitude:.3f}"
    fig.suptitle(title_text, fontsize=14, fontweight='bold')
        
    plt.tight_layout()

    filename_datetime = datetime_filename(model_data)
    fig_filename = f"GGS_{config['glider_name']}_QualityControl.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')

    plt.close(fig)

    print(f"### QC PLOT SAVED ###")
