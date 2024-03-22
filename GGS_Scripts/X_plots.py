# =========================
# IMPORTS
# =========================

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import os

from X_functions import calculate_gridpoint, plot_formatted_ticks, plot_contour_cbar, plot_threshold_legend, plot_bathymetry, plot_profile_thresholds, plot_add_gliders, plot_add_eez, plot_advantage_zones, plot_glider_route, format_colorbar, format_figure_titles, format_subplot_titles, format_subplot_headers, format_save_datetime, print_starttime, print_endtime, print_runtime

# =========================

### FUNCTION:
def GGS_plot_profiles(config, directory, datetime_index, model_datasets, latitude_qc=None, longitude_qc=None, threshold=0.5):

    '''
    Produce quality control profiles for 'u', 'v', 'magnitude', and 'direction' data at the specified point of interest.

    Args:
    - config (dict): The configuration dictionary.
    - directory (str): The directory to save the plot.
    - latitude_qc (float): The latitude of the point of interest.
    - longitude_qc (float): The longitude of the point of interest.
    - threshold (float): The threshold value for the shading.
    - model_datasets (tuple): Tuple containing the three model datasets.

    Returns:
    - None
    '''

    # INITIALIZATION
    print(f"\n### CREATING PROFILE PLOT ###\n")
    start_time = print_starttime()

    latitude_qc = float(latitude_qc)
    longitude_qc = float(longitude_qc)

    # DATA UNPACKING
    rtofs_datasets = None
    cmems_datasets = None
    gofs_datasets = None

    for dataset in model_datasets:
        if dataset:
            model_data, depth_average_data, bin_average_data = dataset
            
            model_name = model_data.attrs.get('model_name')
            
            if model_name == 'RTOFS':
                rtofs_datasets = (model_data, depth_average_data, bin_average_data)
            elif model_name == 'CMEMS':
                cmems_datasets = (model_data, depth_average_data, bin_average_data)
            elif model_name == 'GOFS':
                gofs_datasets = (model_data, depth_average_data, bin_average_data)
            else:
                print(f"Unknown model name: {model_name}")

    # DATASET EXTRACTION
    if rtofs_datasets is not None:
        rtofs_model_data, rtofs_depth_average, rtofs_bin_average = rtofs_datasets

        (y_rtofs_model, x_rtofs_model), _ = calculate_gridpoint(rtofs_model_data, latitude_qc, longitude_qc)
        (y_rtofs_avg, x_rtofs_avg), _ = calculate_gridpoint(rtofs_depth_average, latitude_qc, longitude_qc)
        (y_rtofs_bin, x_rtofs_bin), _ = calculate_gridpoint(rtofs_bin_average, latitude_qc, longitude_qc)
        
        rtofs_model_u = rtofs_model_data['u'].isel(y=y_rtofs_model, x=x_rtofs_model).values
        rtofs_avg_u = rtofs_depth_average['u_depth_avg'].isel(y=y_rtofs_avg, x=x_rtofs_avg).values
        rtofs_avg_u = rtofs_avg_u.item()
        rtofs_bin_u = rtofs_bin_average['u_bin_avg'].isel(y=y_rtofs_bin, x=x_rtofs_bin).values
        rtofs_bin_u1d = rtofs_bin_u[0]

        rtofs_model_v = rtofs_model_data['v'].isel(y=y_rtofs_model, x=x_rtofs_model).values
        rtofs_avg_v = rtofs_depth_average['v_depth_avg'].isel(y=y_rtofs_avg, x=x_rtofs_avg).values
        rtofs_avg_v = rtofs_avg_v.item()
        rtofs_bin_v = rtofs_bin_average['v_bin_avg'].isel(y=y_rtofs_bin, x=x_rtofs_bin).values
        rtofs_bin_v1d = rtofs_bin_v[0]

        rtofs_avg_mag = rtofs_depth_average['mag_depth_avg'].isel(y=y_rtofs_avg, x=x_rtofs_avg).values
        rtofs_avg_mag = rtofs_avg_mag.item()
        rtofs_bin_mag = rtofs_bin_average['mag_bin_avg'].isel(y=y_rtofs_bin, x=x_rtofs_bin).values
        rtofs_bin_mag1d = rtofs_bin_mag[0]

        rtofs_bin_dir = rtofs_bin_average['dir_bin_avg'].isel(y=y_rtofs_bin, x=x_rtofs_bin).values
        rtofs_bin_dir = np.mod(rtofs_bin_dir + 180, 360) - 180

        rtofs_max_depth = rtofs_model_data.depth.max().item()
        rtofs_max_bins = rtofs_max_depth + 1
        rtofs_bin_depths = np.arange(rtofs_max_bins)

    if cmems_datasets is not None:
        cmems_model_data, cmems_depth_average, cmems_bin_average = cmems_datasets
        
        (y_cmems_model, x_cmems_model), _ = calculate_gridpoint(cmems_model_data, latitude_qc, longitude_qc)
        (y_cmems_avg, x_cmems_avg), _ = calculate_gridpoint(cmems_depth_average, latitude_qc, longitude_qc)
        (y_cmems_bin, x_cmems_bin), _ = calculate_gridpoint(cmems_bin_average, latitude_qc, longitude_qc)

        cmems_model_u = cmems_model_data['u'].isel(lat=y_cmems_model, lon=x_cmems_model).values
        cmems_avg_u = cmems_depth_average['u_depth_avg'].isel(lat=y_cmems_avg, lon=x_cmems_avg).values
        cmems_avg_u = cmems_avg_u.item()
        cmems_bin_u = cmems_bin_average['u_bin_avg'].isel(lat=y_cmems_bin, lon=x_cmems_bin).values
        cmems_bin_u1d = cmems_bin_u[0]

        cmems_model_v = cmems_model_data['v'].isel(lat=y_cmems_model, lon=x_cmems_model).values
        cmems_avg_v = cmems_depth_average['v_depth_avg'].isel(lat=y_cmems_avg, lon=x_cmems_avg).values
        cmems_avg_v = cmems_avg_v.item()
        cmems_bin_v = cmems_bin_average['v_bin_avg'].isel(lat=y_cmems_bin, lon=x_cmems_bin).values
        cmems_bin_v1d = cmems_bin_v[0]

        cmems_avg_mag = cmems_depth_average['mag_depth_avg'].isel(lat=y_cmems_avg, lon=x_cmems_avg).values
        cmems_avg_mag = cmems_avg_mag.item()
        cmems_bin_mag = cmems_bin_average['mag_bin_avg'].isel(lat=y_cmems_bin, lon=x_cmems_bin).values
        cmems_bin_mag1d = cmems_bin_mag[0]

        cmems_bin_dir = cmems_bin_average['dir_bin_avg'].isel(lat=y_cmems_bin, lon=x_cmems_bin).values
        cmems_bin_dir = np.mod(cmems_bin_dir + 180, 360) - 180

        cmems_max_depth = cmems_model_data.depth.max().item()
        cmems_max_bins = cmems_max_depth + 1
        cmems_bin_depths = np.arange(cmems_max_bins)

    if gofs_datasets is not None:
        gofs_model_data, gofs_depth_average, gofs_bin_average = gofs_datasets
        
        (y_gofs_model, x_gofs_model), _ = calculate_gridpoint(gofs_model_data, latitude_qc, longitude_qc)
        (y_gofs_avg, x_gofs_avg), _ = calculate_gridpoint(gofs_depth_average, latitude_qc, longitude_qc)
        (y_gofs_bin, x_gofs_bin), _ = calculate_gridpoint(gofs_bin_average, latitude_qc, longitude_qc)

        gofs_model_u = gofs_model_data['u'].isel(lat=y_gofs_model, lon=x_gofs_model).values
        gofs_avg_u = gofs_depth_average['u_depth_avg'].isel(lat=y_gofs_avg, lon=x_gofs_avg).values
        gofs_avg_u = gofs_avg_u.item()
        gofs_bin_u = gofs_bin_average['u_bin_avg'].isel(lat=y_gofs_bin, lon=x_gofs_bin).values
        gofs_bin_u1d = gofs_bin_u[0]

        gofs_model_v = gofs_model_data['v'].isel(lat=y_gofs_model, lon=x_gofs_model).values
        gofs_avg_v = gofs_depth_average['v_depth_avg'].isel(lat=y_gofs_avg, lon=x_gofs_avg).values
        gofs_avg_v = gofs_avg_v.item()
        gofs_bin_v = gofs_bin_average['v_bin_avg'].isel(lat=y_gofs_bin, lon=x_gofs_bin).values
        gofs_bin_v1d = gofs_bin_v[0]

        gofs_avg_mag = gofs_depth_average['mag_depth_avg'].isel(lat=y_gofs_avg, lon=x_gofs_avg).values
        gofs_avg_mag = gofs_avg_mag.item()
        gofs_bin_mag = gofs_bin_average['mag_bin_avg'].isel(lat=y_gofs_bin, lon=x_gofs_bin).values
        gofs_bin_mag1d = gofs_bin_mag[0]

        gofs_bin_dir = gofs_bin_average['dir_bin_avg'].isel(lat=y_gofs_bin, lon=x_gofs_bin).values
        gofs_bin_dir = np.mod(gofs_bin_dir + 180, 360) - 180

        gofs_max_depth = gofs_model_data.depth.max().item()
        gofs_max_bins = gofs_max_depth + 1
        gofs_bin_depths = np.arange(gofs_max_bins)

    # PLOTTING SETUP
    fig, axes = plt.subplots(3, 4, figsize=(20, 25), gridspec_kw={'height_ratios': [1, 1, 1], 'hspace': 0.3})

    if rtofs_datasets is not None:
        # U-VELOCITY PROFILE
        axes[0,0].scatter(rtofs_model_u, rtofs_model_data.depth, marker='x', color='black', s=100, label='Model Datapoint', alpha=1.0, zorder=3)
        axes[0,0].scatter(rtofs_bin_u, rtofs_bin_depths, label='1m Interpolation', color='cyan', alpha=1.0, zorder=2)
        axes[0,0].axvline(x=rtofs_avg_u, label=f'Depth Average = [{rtofs_avg_u:.2f}]', color='darkcyan', linestyle='--', linewidth=2, zorder=1)
        axes[0,0].set_xlabel('u Velocity (m/s)', fontsize=12, fontweight='bold')
        axes[0,0].set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
        axes[0,0].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[0,0].invert_yaxis()
    
        # V-VELOCITY PROFILE
        axes[0,1].scatter(rtofs_model_v, rtofs_model_data.depth, label='Model Datapoint', marker='x', color='black', s=100, alpha=1.0, zorder=3)
        axes[0,1].scatter(rtofs_bin_v, rtofs_bin_depths, label='1m Interpolation', color='orange', alpha=1.0, zorder=2)
        axes[0,1].axvline(x=rtofs_avg_v, label=f'Depth Avgerage = [{rtofs_avg_v:.2f}]', color='darkorange', linestyle='--', linewidth=2, zorder=1)
        axes[0,1].set_xlabel('v Velocity (m/s)', fontsize=12, fontweight='bold')
        axes[0,1].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[0,1].invert_yaxis()

        # MAGNITUDE PROFILE
        axes[0,2].scatter(rtofs_bin_mag, rtofs_bin_depths, label='1m Interpolation', color='green', alpha=1.0, zorder=2)
        axes[0,2].set_xlabel('Current Magnitude (m/s)', fontsize=12, fontweight='bold')
        axes[0,2].axvline(x=rtofs_avg_mag, label=f'Depth Avgerage = [{rtofs_avg_mag:.2f}]', color='darkgreen', linestyle='--', linewidth=2, zorder=1)
        axes[0,2].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[0,2].invert_yaxis()

        # DIRECTION PROFILE
        axes[0,3].scatter(rtofs_bin_dir, rtofs_bin_depths, label='1m Interpolation', color='purple', alpha=1.0, zorder=2)
        axes[0,3].set_xlabel('Current Direction (degrees)', fontsize=12, fontweight='bold')
        axes[0,3].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[0,3].invert_yaxis()
    
        # THRESHOLDS
        shading_colors = ['cyan', 'orange', 'lawngreen']
        for i, (data_1d, color) in enumerate(zip([rtofs_bin_u1d, rtofs_bin_v1d, rtofs_bin_mag1d], shading_colors)):
            plot_profile_thresholds(axes[0, i], data_1d, threshold, color)

        # LEGEND
        for row in axes:
            for ax in row:
                handles, labels = ax.get_legend_handles_labels()
                if labels:
                    ax.legend(handles, labels, loc='lower center', facecolor='lightgrey', edgecolor='black', framealpha=1.0)

    if cmems_datasets is not None:
        # U-VELOCITY PROFILE
        axes[1,0].scatter(cmems_model_u, cmems_model_data.depth, marker='x', color='black', s=100, label='Model Datapoint', alpha=1.0, zorder=3)
        axes[1,0].scatter(cmems_bin_u, cmems_bin_depths, label='1m Interpolation', color='cyan', alpha=1.0, zorder=2)
        axes[1,0].axvline(x=cmems_avg_u, label=f'Depth Average = [{cmems_avg_u:.2f}]', color='darkcyan', linestyle='--', linewidth=2, zorder=1)
        axes[1,0].set_xlabel('u Velocity (m/s)', fontsize=12, fontweight='bold')
        axes[1,0].set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
        axes[1,0].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[1,0].invert_yaxis()

        # V-VELOCITY PROFILE
        axes[1,1].scatter(cmems_model_v, cmems_model_data.depth, label='Model Datapoint', marker='x', color='black', s=100, alpha=1.0, zorder=3)
        axes[1,1].scatter(cmems_bin_v, cmems_bin_depths, label='1m Interpolation', color='orange', alpha=1.0, zorder=2)
        axes[1,1].axvline(x=cmems_avg_v, label=f'Depth Avgerage = [{cmems_avg_v:.2f}]', color='darkorange', linestyle='--', linewidth=2, zorder=1)
        axes[1,1].set_xlabel('v Velocity (m/s)', fontsize=12, fontweight='bold')
        axes[1,1].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[1,1].invert_yaxis()
        
        # MAGNITUDE PROFILE
        axes[1,2].scatter(cmems_bin_mag, cmems_bin_depths, label='1m Interpolation', color='green', alpha=1.0, zorder=2)
        axes[1,2].set_xlabel('Current Magnitude (m/s)', fontsize=12, fontweight='bold')
        axes[1,2].axvline(x=cmems_avg_mag, label=f'Depth Avgerage = [{cmems_avg_mag:.2f}]', color='darkgreen', linestyle='--', linewidth=2, zorder=1)
        axes[1,2].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[1,2].invert_yaxis()

        # DIRECTION PROFILE
        axes[1,3].scatter(cmems_bin_dir, cmems_bin_depths, label='1m Interpolation', color='purple', alpha=1.0, zorder=2)
        axes[1,3].set_xlabel('Current Direction (degrees)', fontsize=12, fontweight='bold')
        axes[1,3].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[1,3].invert_yaxis()

        # THRESHOLDS
        shading_colors = ['cyan', 'orange', 'lawngreen']
        for i, (data_1d, color) in enumerate(zip([cmems_bin_u1d, cmems_bin_v1d, cmems_bin_mag1d], shading_colors)):
            plot_profile_thresholds(axes[1, i], data_1d, threshold, color)

        # LEGEND
        for row in axes:
            for ax in row:
                handles, labels = ax.get_legend_handles_labels()
                if labels:
                    ax.legend(handles, labels, loc='lower center', facecolor='lightgrey', edgecolor='black', framealpha=1.0)

    if gofs_datasets is not None:
        # U-VELOCITY PROFILE
        axes[2,0].scatter(gofs_model_u, gofs_model_data.depth, marker='x', color='black', s=100, label='Model Datapoint', alpha=1.0, zorder=3)
        axes[2,0].scatter(gofs_bin_u, gofs_bin_depths, label='1m Interpolation', color='cyan', alpha=1.0, zorder=2)
        axes[2,0].axvline(x=gofs_avg_u, label=f'Depth Average = [{gofs_avg_u:.2f}]', color='darkcyan', linestyle='--', linewidth=2, zorder=1)
        axes[2,0].set_xlabel('u Velocity (m/s)', fontsize=12, fontweight='bold')
        axes[2,0].set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
        axes[2,0].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[2,0].invert_yaxis()

        # V-VELOCITY PROFILE
        axes[2,1].scatter(gofs_model_v, gofs_model_data.depth, label='Model Datapoint', marker='x', color='black', s=100, alpha=1.0, zorder=3)
        axes[2,1].scatter(gofs_bin_v, gofs_bin_depths, label='1m Interpolation', color='orange', alpha=1.0, zorder=2)
        axes[2,1].axvline(x=gofs_avg_v, label=f'Depth Avgerage = [{gofs_avg_v:.2f}]', color='darkorange', linestyle='--', linewidth=2, zorder=1)
        axes[2,1].set_xlabel('v Velocity (m/s)', fontsize=12, fontweight='bold')
        axes[2,1].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[2,1].invert_yaxis()

        # MAGNITUDE PROFILE
        axes[2,2].scatter(gofs_bin_mag, gofs_bin_depths, label='1m Interpolation', color='green', alpha=1.0, zorder=2)
        axes[2,2].set_xlabel('Current Magnitude (m/s)', fontsize=12, fontweight='bold')
        axes[2,2].axvline(x=gofs_avg_mag, label=f'Depth Avgerage = [{gofs_avg_mag:.2f}]', color='darkgreen', linestyle='--', linewidth=2, zorder=1)
        axes[2,2].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[2,2].invert_yaxis()
        
        # DIRECTION PROFILE
        axes[2,3].scatter(gofs_bin_dir, gofs_bin_depths, label='1m Interpolation', color='purple', alpha=1.0, zorder=2)
        axes[2,3].set_xlabel('Current Direction (degrees)', fontsize=12, fontweight='bold')
        axes[2,3].grid(color='lightgrey', linestyle='-', linewidth=0.5)
        axes[2,3].invert_yaxis()

        # THRESHOLDS
        shading_colors = ['cyan', 'orange', 'lawngreen']
        for i, (data_1d, color) in enumerate(zip([gofs_bin_u1d, gofs_bin_v1d, gofs_bin_mag1d], shading_colors)):
            plot_profile_thresholds(axes[2, i], data_1d, threshold, color)

        # LEGEND
        for row in axes:
            for ax in row:
                handles, labels = ax.get_legend_handles_labels()
                if labels:
                    ax.legend(handles, labels, loc='lower center', facecolor='lightgrey', edgecolor='black', framealpha=1.0)
    
    # TITLES
    title_text = f"Vertical Profile Subplots - (Lat: {latitude_qc:.3f}, Lon: {longitude_qc:.3f})"
    format_subplot_titles(fig, config, datetime_index, title=title_text)

    subplot_titles = ['RTOFS', 'CMEMS', 'GOFS']
    format_subplot_headers(axes, fig, subplot_titles)

    # SAVE & CLOSE
    file_datetime = format_save_datetime(datetime_index)
    fig_filename = f"DepthAverageProfiles_{config['max_depth']}m_{file_datetime}.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # LOGGING
    end_time = print_endtime()
    print_runtime(start_time, end_time)

### FUNCTION:
def GGS_plot_magnitudes(config, directory, datetime_index, model_datasets, latitude_qc=None, longitude_qc=None, density=2, gliders=None, show_route=False, show_qc=False, show_eez=False, manual_extent=None):
    
    '''
    Plot the depth-averaged current fields from three datasets side by side.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Directory to save the plot.
    - datetime_index (int): Index of the datetime for the plot title.
    - model_datasets (tuple): Tuple containing the model datasets.
    - latitude_qc (float): Latitude for QC plotting.
    - longitude_qc (float): Longitude for QC plotting.
    - density (int): Density of the streamplot.
    - gliders (optional): DataFrame containing glider data for plotting.
    - show_route (bool): Flag to show the glider route.
    - show_qc (bool): Flag to show the QC sample point.
    - show_eez (bool): Flag to show the Exclusive Economic Zone (EEZ).
    - manual_extent (list or None): Manual specification of plot extent.

    Returns:
    - None
    '''

    # INITIALIZATION
    print(f"\n### CREATING MAGNITUDE PLOT ###\n")
    start_time = print_starttime()

    # DATASET EXTRACTION
    valid_datasets = [datasets for datasets in model_datasets if datasets is not None]
    num_datasets = len(valid_datasets)
    if num_datasets == 0:
        print("No datasets provided for plotting.")
        return

    # PLOTTING FUNCTION
    def plot_magnitude(ax, config, model_depth_average, latitude_qc, longitude_qc, density, gliders, show_route, show_qc, show_eez, manual_extent):
        
        # DATA EXTRACTION
        longitude = model_depth_average.lon.values.squeeze()
        latitude = model_depth_average.lat.values.squeeze()
        u_depth_avg = model_depth_average['u_depth_avg'].values.squeeze()
        v_depth_avg = model_depth_average['v_depth_avg'].values.squeeze()
        mag_depth_avg = model_depth_average['mag_depth_avg'].values.squeeze()
        
        # EXTENT SETUP
        if manual_extent == "config":
            map_extent = [config['extent'][0][1], config['extent'][1][1], config['extent'][0][0], config['extent'][1][0]]
        elif isinstance(manual_extent, list) and len(manual_extent) == 4:
            map_extent = manual_extent
        else:
            data_extent_lon = [float(longitude.min()), float(longitude.max())]
            data_extent_lat = [float(latitude.min()), float(latitude.max())]
            map_extent = data_extent_lon + data_extent_lat
        ax.set_extent(map_extent, crs=ccrs.PlateCarree())
        plot_formatted_ticks(ax, map_extent[:2], map_extent[2:], proj=ccrs.PlateCarree(), fontsize=16, label_left=True, label_right=False, label_bottom=True, label_top=False, gridlines=True)

        # PLOT ELEMENTS
        levels, ticks, extend = plot_contour_cbar(mag_depth_avg, max_levels=10, extend_max=True)
        contourf = ax.contourf(longitude, latitude, mag_depth_avg, levels=levels, cmap=cmo.speed, transform=ccrs.PlateCarree(), zorder=10, extend=extend)
        streamplot = ax.streamplot(longitude, latitude, u_depth_avg, v_depth_avg, transform=ccrs.PlateCarree(), density=density, linewidth=0.5, color='black', zorder=10)

        # GLIDERS
        if gliders is not None:
            plot_add_gliders(ax, gliders, legend=True)
            glider_legend = ax.get_legend()
            if glider_legend:
                glider_legend.get_frame().set_alpha(0.5)
                glider_legend.get_frame().set_facecolor('white')
                ax.add_artist(glider_legend)
        
        # ROUTE
        if show_route:
            plot_glider_route(ax, config)
        
        # QUALITY CONTROL
        if show_qc:
            (y_index, x_index), (lat_index, lon_index) = calculate_gridpoint(model_depth_average, latitude_qc, longitude_qc)
            qc_lon = model_depth_average['lon'].isel(x=x_index, y=y_index).values
            qc_lat = model_depth_average['lat'].isel(x=x_index, y=y_index).values
            circle = Circle((qc_lon, qc_lat), radius=0.25, edgecolor='purple', facecolor='none', linewidth=2, transform=ccrs.PlateCarree(), zorder=95)
            ax.add_patch(circle)
        
        # EEZ
        if show_eez:
            plot_add_eez(ax, config, color='dimgrey', linewidth=3, zorder=90)
        
        # FEATURES
        ax.add_feature(cfeature.GSHHSFeature(scale='full'), edgecolor="black", facecolor="tan", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.RIVERS, edgecolor="steelblue", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.LAKES, edgecolor="black", facecolor="lightsteelblue", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=0.25, zorder=90)

        plot_bathymetry(ax, config, model_depth_average, isobath1=-100, isobath2=-1000, downsample=True, show_legend=False)
        bathymetry_legend = ax.get_legend()
        if bathymetry_legend:
            bathymetry_legend.get_frame().set_alpha(0.5)
            ax.add_artist(bathymetry_legend)

        # COLORBAR
        cbar = fig.colorbar(contourf, orientation='vertical', extend=extend)
        format_colorbar(ax, cbar)
        cbar.set_label('Depth Averaged Current Magnitude (m/s)', labelpad=10)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f"{tick:.1f}" for tick in ticks])

    # PLOTTING
    fig, axs = plt.subplots(1, num_datasets, subplot_kw={'projection': ccrs.Mercator()}, figsize=(10*num_datasets, 10))
    if num_datasets == 1:
        axs = [axs]
    
    model_names = []
    for ax, (model_data, depth_average_data, bin_average_data) in zip(axs, valid_datasets):
        model_name = depth_average_data.attrs['model_name']
        model_names.append(model_name)
        plot_magnitude(ax, config, depth_average_data, latitude_qc, longitude_qc, density, gliders, show_route, show_qc, show_eez, manual_extent)
        ax.set_title(f"{model_name}", fontsize=14, fontweight='bold', pad=20)

    title_text = f"Depth Averaged Currents - Depth Range: {config['max_depth']}m"
    model_names_combined = " vs. ".join(model_names)
    format_figure_titles(axs[0], fig, config, datetime_index, model_name=model_names_combined, title=title_text)

    # SAVE & CLOSE
    file_datetime = format_save_datetime(datetime_index)
    fig_filename = f"DepthAverageMagnitude_{config['max_depth']}m_{file_datetime}.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # LOGGING
    end_time = print_endtime()
    print_runtime(start_time, end_time)

### FUNCTION:
def GGS_plot_threshold(config, directory, datetime_index, model_datasets, latitude_qc=None, longitude_qc=None, density=2, mag1=0.0, mag2=0.2, mag3=0.3, mag4=0.4, mag5=0.5, gliders=None, show_route=False, show_qc=False, show_eez=False, manual_extent=None):
    
    '''
    Plot the depth-averaged current fields from three datasets side by side.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Directory to save the plot.
    - datetime_index (int): Index of the datetime for the plot title.
    - model_datasets (tuple): Tuple containing the model datasets.
    - latitude_qc (float): Latitude for QC plotting.
    - longitude_qc (float): Longitude for QC plotting.
    - density (int): Density of the streamplot.
    - mag1 (float): Threshold for the first magnitude level.
    - mag2 (float): Threshold for the second magnitude level.
    - mag3 (float): Threshold for the third magnitude level.
    - mag4 (float): Threshold for the fourth magnitude level.
    - mag5 (float): Threshold for the fifth magnitude level.
    - gliders (optional): DataFrame containing glider data for plotting.
    - show_route (bool): Flag to show the glider route.
    - show_qc (bool): Flag to show the QC sample point.
    - show_eez (bool): Flag to show the Exclusive Economic Zone (EEZ).
    - manual_extent (list or None): Manual specification of plot extent.

    Returns:
    - None
    '''

    # INITIALIZATION
    print(f"\n### CREATING THRESHOLD PLOT ###\n")
    start_time = print_starttime()

    # DATASET EXTRACTION
    valid_datasets = [datasets for datasets in model_datasets if datasets is not None]
    num_datasets = len(valid_datasets)
    if num_datasets == 0:
        print("No datasets provided for plotting.")
        return

    # PLOTTING FUNCTION
    def plot_threshold(ax, config, model_depth_average, latitude_qc, longitude_qc, density, mag1, mag2, mag3, mag4, mag5, gliders, show_route, show_qc, show_eez, manual_extent):
        
        # DATA EXTRACTION
        longitude = model_depth_average.lon.values.squeeze()
        latitude = model_depth_average.lat.values.squeeze()
        u_depth_avg = model_depth_average['u_depth_avg'].values.squeeze()
        v_depth_avg = model_depth_average['v_depth_avg'].values.squeeze()
        mag_depth_avg = model_depth_average['mag_depth_avg'].values.squeeze()
        
        # EXTENT SETUP
        if manual_extent == "config":
            map_extent = [config['extent'][0][1], config['extent'][1][1], config['extent'][0][0], config['extent'][1][0]]
        elif isinstance(manual_extent, list) and len(manual_extent) == 4:
            map_extent = manual_extent
        else:
            data_extent_lon = [float(longitude.min()), float(longitude.max())]
            data_extent_lat = [float(latitude.min()), float(latitude.max())]
            map_extent = data_extent_lon + data_extent_lat
        ax.set_extent(map_extent, crs=ccrs.PlateCarree())
        plot_formatted_ticks(ax, map_extent[:2], map_extent[2:], proj=ccrs.PlateCarree(), fontsize=16, label_left=True, label_right=True, label_bottom=True, label_top=False, gridlines=True)

        # PLOT ELEMENTS
        levels = [mag1, mag2, mag3, mag4, mag5, np.nanmax(mag_depth_avg)]
        colors = ['none', 'yellow', 'orange', 'orangered', 'maroon']
        contourf = ax.contourf(longitude, latitude, mag_depth_avg, levels=levels, colors=colors, extend='both', transform=ccrs.PlateCarree(), zorder=10)
        streamplot = ax.streamplot(longitude, latitude, u_depth_avg, v_depth_avg, transform=ccrs.PlateCarree(), density=density, linewidth=0.5, color='black', zorder=10)
        streamplot.lines.set_alpha(1.0)
        plot_threshold_legend(ax, mag2, mag3, mag4, mag5)
        threshold_legend = ax.get_legend()
        if threshold_legend:
            threshold_legend.get_frame().set_alpha(0.75)
            threshold_legend.get_frame().set_facecolor('white')
            ax.add_artist(threshold_legend)

        # GLIDERS
        if gliders is not None:
            plot_add_gliders(ax, gliders, legend=True)
            glider_legend = ax.get_legend()
            if glider_legend:
                glider_legend.get_frame().set_alpha(0.5)
                glider_legend.get_frame().set_facecolor('white')
                ax.add_artist(glider_legend)

        # ROUTE
        if show_route:
            plot_glider_route(ax, config)
        
        # QUALITY CONTROL
        if show_qc:
            (y_index, x_index), (lat_index, lon_index) = calculate_gridpoint(model_depth_average, latitude_qc, longitude_qc)
            qc_lon = model_depth_average['lon'].isel(x=x_index, y=y_index).values
            qc_lat = model_depth_average['lat'].isel(x=x_index, y=y_index).values
            circle = Circle((qc_lon, qc_lat), radius=0.25, edgecolor='purple', facecolor='none', linewidth=2, transform=ccrs.PlateCarree(), zorder=95)
            ax.add_patch(circle)

        # EEZ
        if show_eez:
            plot_add_eez(ax, config, color='dimgrey', linewidth=3, zorder=90)
        
        # FEATURES
        ax.add_feature(cfeature.GSHHSFeature(scale='full'), edgecolor="black", facecolor="tan", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.RIVERS, edgecolor="steelblue", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.LAKES, edgecolor="black", facecolor="lightsteelblue", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=0.25, zorder=90)

        plot_bathymetry(ax, config, model_depth_average, isobath1=-100, isobath2=-1000, downsample=True, show_legend=False)
        bathymetry_legend = ax.get_legend()
        if bathymetry_legend:
            bathymetry_legend.get_frame().set_alpha(0.5)
            ax.add_artist(bathymetry_legend)
        
    # PLOT SETUP
    fig, axs = plt.subplots(1, num_datasets, subplot_kw={'projection': ccrs.Mercator()}, figsize=(10*num_datasets, 10))
    if num_datasets == 1:
        axs = [axs]

    model_names = []
    for ax, (model_data, depth_average_data, bin_average_data) in zip(axs, valid_datasets):
        model_name = depth_average_data.attrs['model_name']
        model_names.append(model_name)
        plot_threshold(ax, config, depth_average_data, latitude_qc, longitude_qc, density, mag1, mag2, mag3, mag4, mag5, gliders, show_route, show_qc, show_eez, manual_extent)
        ax.set_title(f"{model_name}", fontsize=14, fontweight='bold', pad=20)
    
    title_text = f"Depth Averaged Current Threshold Zones - Depth Range: {config['max_depth']}m"
    model_names_combined = " vs. ".join(model_names)
    format_figure_titles(axs[0], fig, config, datetime_index, model_name=model_names_combined, title=title_text)

    # SAVE & CLOSE
    file_datetime = format_save_datetime(datetime_index)
    fig_filename = f"DepthAverageThreshold_{config['max_depth']}m_{file_datetime}.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    # LOGGING
    end_time = print_endtime()
    print_runtime(start_time, end_time)

### FUNCTION:
def GGS_plot_advantage(config, directory, datetime_index, model_datasets, latitude_qc=None, longitude_qc=None, density=2, tolerance=15, mag1=0.0, mag2=0.2, mag3=0.3, mag4=0.4, mag5=0.5, gliders=None, show_route=False, show_qc=False, show_eez=False, manual_extent=None):
    
    '''
    Plot the depth-averaged current fields from three datasets side by side.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Directory to save the plot.
    - datetime_index (int): Index of the datetime for the plot title.
    - model_datasets (tuple): Tuple containing the model datasets.
    - latitude_qc (float): Latitude for QC plotting.
    - longitude_qc (float): Longitude for QC plotting.
    - density (int): Density of the streamplot.
    - tolerance (float): Tolerance for the bearing.
    - mag1 (float): Threshold for the first magnitude level.
    - mag2 (float): Threshold for the second magnitude level.
    - mag3 (float): Threshold for the third magnitude level.
    - mag4 (float): Threshold for the fourth magnitude level.
    - mag5 (float): Threshold for the fifth magnitude level.
    - gliders (optional): DataFrame containing glider data for plotting.
    - show_route (bool): Flag to show the glider route.
    - show_qc (bool): Flag to show the QC sample point.
    - show_eez (bool): Flag to show the Exclusive Economic Zone (EEZ).
    - manual_extent (list or None): Manual specification of plot extent.

    Returns:
    - None
    '''

    # INITIALIZATION
    print(f"\n### CREATING ADVANTAGE PLOT ###\n")
    start_time = print_starttime()

    # DATASET EXTRACTION
    valid_datasets = [datasets for datasets in model_datasets if datasets is not None]
    num_datasets = len(valid_datasets)
    if num_datasets == 0:
        print("No datasets provided for plotting.")
        return

    # PLOTTING FUNCTION
    def plot_advantage(ax, config, model_depth_average, latitude_qc, longitude_qc, density, tolerance, mag1, mag2, mag3, mag4, mag5, gliders, show_route, show_qc, show_eez, manual_extent):
        
        # DATA EXTRACTION
        longitude = model_depth_average.lon.values.squeeze()
        latitude = model_depth_average.lat.values.squeeze()
        u_depth_avg = model_depth_average['u_depth_avg'].values.squeeze()
        v_depth_avg = model_depth_average['v_depth_avg'].values.squeeze()
        mag_depth_avg = model_depth_average['mag_depth_avg'].values.squeeze()
        
        # EXTENT SETUP
        if manual_extent == "config":
            map_extent = [config['extent'][0][1], config['extent'][1][1], config['extent'][0][0], config['extent'][1][0]]
        elif isinstance(manual_extent, list) and len(manual_extent) == 4:
            map_extent = manual_extent
        else:
            data_extent_lon = [float(longitude.min()), float(longitude.max())]
            data_extent_lat = [float(latitude.min()), float(latitude.max())]
            map_extent = data_extent_lon + data_extent_lat
        ax.set_extent(map_extent, crs=ccrs.PlateCarree())
        plot_formatted_ticks(ax, map_extent[:2], map_extent[2:], proj=ccrs.PlateCarree(), fontsize=16, label_left=True, label_right=True, label_bottom=True, label_top=False, gridlines=True)

        # PLOT ELEMENTS
        levels = [mag1, mag2, mag3, mag4, mag5, np.nanmax(mag_depth_avg)]
        colors = ['none', 'yellow', 'orange', 'orangered', 'maroon']
        contourf = ax.contourf(longitude, latitude, mag_depth_avg, levels=levels, colors=colors, extend='both', transform=ccrs.PlateCarree(), zorder=10)
        plot_advantage_zones(ax, config, model_depth_average, tolerance)
        streamplot = ax.streamplot(longitude, latitude, u_depth_avg, v_depth_avg, transform=ccrs.PlateCarree(), density=density, linewidth=0.5, color='black', zorder=10)
        streamplot.lines.set_alpha(1.0)
        plot_threshold_legend(ax, mag2, mag3, mag4, mag5)
        threshold_legend = ax.get_legend()
        if threshold_legend:
            threshold_legend.get_frame().set_alpha(0.75)
            threshold_legend.get_frame().set_facecolor('white')
            ax.add_artist(threshold_legend)

        # GLIDERS
        if gliders is not None:
            plot_add_gliders(ax, gliders, legend=True)
            glider_legend = ax.get_legend()
            if glider_legend:
                glider_legend.get_frame().set_alpha(0.5)
                glider_legend.get_frame().set_facecolor('white')
                ax.add_artist(glider_legend)

        # ROUTE
        if show_route:
            plot_glider_route(ax, config)
        
        # QUALITY CONTROL
        if show_qc:
            (y_index, x_index), (lat_index, lon_index) = calculate_gridpoint(model_depth_average, latitude_qc, longitude_qc)
            qc_lon = model_depth_average['lon'].isel(x=x_index, y=y_index).values
            qc_lat = model_depth_average['lat'].isel(x=x_index, y=y_index).values
            circle = Circle((qc_lon, qc_lat), radius=0.25, edgecolor='purple', facecolor='none', linewidth=2, transform=ccrs.PlateCarree(), zorder=95)
            ax.add_patch(circle)

        # EEZ
        if show_eez:
            plot_add_eez(ax, config, color='dimgrey', linewidth=3, zorder=90)
        
        # FEATURES
        ax.add_feature(cfeature.GSHHSFeature(scale='full'), edgecolor="black", facecolor="tan", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.RIVERS, edgecolor="steelblue", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.LAKES, edgecolor="black", facecolor="lightsteelblue", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=0.25, zorder=90)

        plot_bathymetry(ax, config, model_depth_average, isobath1=-100, isobath2=-1000, downsample=True, show_legend=False)
        bathymetry_legend = ax.get_legend()
        if bathymetry_legend:
            bathymetry_legend.get_frame().set_alpha(0.5)
            ax.add_artist(bathymetry_legend)
        
    # PLOT SETUP
    fig, axs = plt.subplots(1, num_datasets, subplot_kw={'projection': ccrs.Mercator()}, figsize=(10*num_datasets, 10))
    if num_datasets == 1:
        axs = [axs]

    model_names = []
    for ax, (model_data, depth_average_data, bin_average_data) in zip(axs, valid_datasets):
        model_name = depth_average_data.attrs['model_name']
        model_names.append(model_name)
        plot_advantage(ax, config, depth_average_data, latitude_qc, longitude_qc, density, tolerance, mag1, mag2, mag3, mag4, mag5, gliders, show_route, show_qc, show_eez, manual_extent)
        ax.set_title(f"{model_name}", fontsize=14, fontweight='bold', pad=20)
    
    title_text = f"Depth Averaged Current Advantage Zones - Depth Range: {config['max_depth']}m"
    model_names_combined = " vs. ".join(model_names)
    format_figure_titles(axs[0], fig, config, datetime_index, model_name=model_names_combined, title=title_text)

    # SAVE & CLOSE
    file_datetime = format_save_datetime(datetime_index)
    fig_filename = f"DepthAverageAdvantage_{config['max_depth']}m_{file_datetime}.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    # LOGGING
    end_time = print_endtime()
    print_runtime(start_time, end_time)
    