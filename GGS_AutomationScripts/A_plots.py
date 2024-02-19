# =========================
# IMPORTS
# =========================

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import os
from A_functions import calculate_gridpoint, plot_profile_thresholds, plot_formatted_ticks, plot_contour_cbar, plot_bathymetry, acquire_gliders, plot_add_gliders, format_colorbar, plot_threshold_legend, format_titles, format_save_datetime, print_starttime, print_endtime, print_runtime

# =========================

### FUNCTION:
def GGS_plot_profiles(config, directory, model_data, depth_average_data, bin_average_data, latitude_qc, longitude_qc, threshold=0.5):
    
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
    - threshold (float): Threshold value for the profile plots.
        - default: 0.5

    Returns:
    - None
    '''

    # INITIALIZATION
    print("\n### CREATING PROFILE PLOT ###\n")
    start_time = print_starttime()

    latitude_qc = float(latitude_qc)
    longitude_qc = float(longitude_qc)

    # DATA EXTRACTION
    (y_model_index, x_model_index), _ = calculate_gridpoint(model_data, latitude_qc, longitude_qc)
    (y_avg_index, x_avg_index), _ = calculate_gridpoint(depth_average_data, latitude_qc, longitude_qc)
    (y_bin_index, x_bin_index), _ = calculate_gridpoint(bin_average_data, latitude_qc, longitude_qc)

    u_model = model_data['u'].isel(y=y_model_index, x=x_model_index).values
    u_depth_avg = depth_average_data['u_depth_avg'].isel(y=y_avg_index, x=x_avg_index).values
    u_depth_avg = u_depth_avg.item()
    u_bin_avg = bin_average_data['u_bin_avg'].isel(y=y_bin_index, x=x_bin_index).values
    u_bin_avg_1d = u_bin_avg[0]

    v_model = model_data['v'].isel(y=y_model_index, x=x_model_index).values
    v_depth_avg = depth_average_data['v_depth_avg'].isel(y=y_avg_index, x=x_avg_index).values
    v_depth_avg = v_depth_avg.item()
    v_bin_avg = bin_average_data['v_bin_avg'].isel(y=y_bin_index, x=x_bin_index).values
    v_bin_avg_1d = v_bin_avg[0]

    mag_depth_avg = depth_average_data['mag_depth_avg'].isel(y=y_avg_index, x=x_avg_index).values
    mag_depth_avg = mag_depth_avg.item()
    mag_bin_avg = bin_average_data['mag_bin_avg'].isel(y=y_bin_index, x=x_bin_index).values
    mag_bin_avg_1d = mag_bin_avg[0]

    dir_bin_avg = bin_average_data['dir_bin_avg'].isel(y=y_bin_index, x=x_bin_index).values
    dir_bin_avg = np.mod(dir_bin_avg + 180, 360) - 180

    max_depth = model_data.depth.max().item()
    max_bins = max_depth + 1
    bin_depths = np.arange(max_bins)

    # PLOTTING SETUP
    fig, axes = plt.subplots(1, 4, figsize=(15, 9))

    # U-VELOCITY PROFILE
    axes[0].scatter(u_model, model_data.depth, marker='x', color='black', s=100, label='Model Datapoint', alpha=1.0, zorder=3)
    axes[0].scatter(u_bin_avg, bin_depths, label='1m Interpolation', color='cyan', alpha=1.0, zorder=2)
    axes[0].axvline(x=u_depth_avg, label=f'Depth Average = [{u_depth_avg:.2f}]', color='darkcyan', linestyle='--', linewidth=2, zorder=1)
    axes[0].set_xlabel('u Velocity (m/s)', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Depth (m)', fontsize=12, fontweight='bold')
    axes[0].grid(color='lightgrey', linestyle='-', linewidth=0.5)
    axes[0].invert_yaxis()
    
    # V-VELOCITY PROFILE
    axes[1].scatter(v_model, model_data.depth, label='Model Datapoint', marker='x', color='black', s=100, alpha=1.0, zorder=3)
    axes[1].scatter(v_bin_avg, bin_depths, label='1m Interpolation', color='orange', alpha=1.0, zorder=2)
    axes[1].axvline(x=v_depth_avg, label=f'Depth Avgerage = [{v_depth_avg:.2f}]', color='darkorange', linestyle='--', linewidth=2, zorder=1)
    axes[1].set_xlabel('v Velocity (m/s)', fontsize=12, fontweight='bold')
    axes[1].grid(color='lightgrey', linestyle='-', linewidth=0.5)
    axes[1].invert_yaxis()

    # MAGNITUDE PROFILE
    axes[2].scatter(mag_bin_avg, bin_depths, label='1m Interpolation', color='green', alpha=1.0, zorder=2)
    axes[2].set_xlabel('Current Magnitude (m/s)', fontsize=12, fontweight='bold')
    axes[2].axvline(x=mag_depth_avg, label=f'Depth Avgerage = [{mag_depth_avg:.2f}]', color='darkgreen', linestyle='--', linewidth=2, zorder=1)
    axes[2].grid(color='lightgrey', linestyle='-', linewidth=0.5)
    axes[2].invert_yaxis()

    # DIRECTION PROFILE
    axes[3].scatter(dir_bin_avg, bin_depths, label='1m Interpolation', color='purple', alpha=1.0, zorder=2)
    axes[3].set_xlabel('Current Direction (degrees)', fontsize=12, fontweight='bold')
    # axes[3].axvline(x=dir_depth_avg, label=f'Depth Avgerage = [{dir_depth_avg:.2f}]', color='darkpurple', linestyle='--', linewidth=2, zorder=1)
    axes[3].grid(color='lightgrey', linestyle='-', linewidth=0.5)
    axes[3].invert_yaxis()
    
    # THRESHOLDS
    shading_colors = ['cyan', 'orange', 'green']
    for ax, data_1d, color in zip(axes[:3], [u_bin_avg_1d, v_bin_avg_1d, mag_bin_avg_1d], shading_colors):
        plot_profile_thresholds(ax, data_1d, threshold, color)

    # LEGEND
    for ax in axes:
        ax.legend(loc='lower center', facecolor='lightgrey', edgecolor='black', framealpha=1.0)

    # TITLE
    title_text = f"{config['glider_name']} Mission - Vertical Profiles - Lat: {latitude_qc:.3f}, Lon: {longitude_qc:.3f}"
    fig.suptitle(title_text, fontsize=14, fontweight='bold')
    
    # SAVE & CLOSE
    file_datetime = format_save_datetime(model_data)
    fig_filename = f"GGS_{config['glider_name']}_Profiles_{file_datetime}.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    # LOGGING
    end_time = print_endtime()
    print_runtime(start_time, end_time)

### FUNCTION:
def GGS_plot_magnitude(config, directory, depth_average_data, latitude_qc, longitude_qc, model_name=None, density=2, show_gliders=False, show_route=False, show_qc=False, manual_extent=None):
    
    '''
    Plot the depth-averaged current fields.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - depth_average_data (xarray.Dataset): Dataset with the computed variables and layer information.
    - qc_latitude (float): Latitude of the QC sample point.
    - qc_longitude (float): Longitude of the QC sample point.
    - model_name (str): Name of the model.
        - default: 'None'
    - density (int): Density of the streamplot.
        - default: 2
    - show_route (bool): Show the route on the plot.
        - default: 'False'
    - show_qc (bool): Show the QC sample point.
        - default: 'False'
    - manual_extent (list or str): Optional manual specification of plot extent as [west, east, south, north], or "config" to use extent from config.
        - default: 'None'

    Returns:
    - None
    '''

    # INITIALIZATION
    if model_name is None:
        model_name = "Unknown"
    else:
        model_name = model_name.upper()

    print(f"\n### CREATING {model_name} MAGNITUDE PLOT ###\n")
    start_time = print_starttime()

    # DATA EXTRACTION
    longitude = depth_average_data.lon.values.squeeze()
    latitude = depth_average_data.lat.values.squeeze()
    u_depth_avg = depth_average_data['u_depth_avg'].values.squeeze()
    v_depth_avg = depth_average_data['v_depth_avg'].values.squeeze()
    mag_depth_avg = depth_average_data['mag_depth_avg'].values.squeeze()
    
    # PLOT SETUP
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    
    # EXTENT SETUP
    if manual_extent == "config":
        map_extent = [config['extent'][0][1], config['extent'][1][1], config['extent'][0][0], config['extent'][1][0]]
    elif isinstance(manual_extent, list) and len(manual_extent) == 4:
        map_extent = manual_extent
    else:
        data_extent_lon = [float(depth_average_data.lon.min()), float(depth_average_data.lon.max())]
        data_extent_lat = [float(depth_average_data.lat.min()), float(depth_average_data.lat.max())]
        map_extent = data_extent_lon + data_extent_lat
    ax.set_extent(map_extent, crs=ccrs.PlateCarree())
    plot_formatted_ticks(ax, map_extent[:2], map_extent[2:], proj=ccrs.PlateCarree(), fontsize=10, label_left=True, label_right=False, label_bottom=True, label_top=False, gridlines=True)
    
    # PLOT ELEMENTS
    levels, ticks, extend = plot_contour_cbar(mag_depth_avg, max_levels=10, extend_max=True)
    contourf = ax.contourf(longitude, latitude, mag_depth_avg, levels=levels, cmap=cmo.speed, transform=ccrs.PlateCarree(), zorder=10, extend=extend)
    streamplot = ax.streamplot(longitude, latitude, u_depth_avg, v_depth_avg, density=density, linewidth=0.5, color='black', transform=ccrs.PlateCarree(), zorder=10)

    # GLIDERS
    if show_gliders:
        glider_dataframes = acquire_gliders(extent=map_extent, target_date=dt.datetime.now(), date_delta=dt.timedelta(days=1), requested_variables=["time", "longitude", "latitude", "profile_id", "depth"], print_vars=False, target="all", request_timeout=5, enable_parallel=False)
        plot_add_gliders(ax, glider_dataframes, legend=True)
        glider_legend = ax.get_legend()
        if glider_legend:
            glider_legend.get_frame().set_alpha(0.5)
            ax.add_artist(glider_legend)
    
    # ROUTE
    if show_route:
        lats, lons = zip(*config["GPS_coords"])
        ax.plot(lons, lats, 'w-', transform=ccrs.PlateCarree(), linewidth=2.5, zorder=91)
        ax.plot(lons, lats, 'k', transform=ccrs.PlateCarree(), linewidth=1.0, linestyle='--', alpha=0.6, zorder=92)
        
        start_coords = config["GPS_coords"][0]
        end_coords = config["GPS_coords"][-1]
        ax.scatter(*start_coords[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
        for GPS_coord in config["GPS_coords"][1:-1]:
            ax.scatter(*GPS_coord[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
        ax.scatter(*end_coords[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
    
    # QUALITY CONTROL
    if show_qc:
        (y_index, x_index), (lat_index, lon_index) = calculate_gridpoint(depth_average_data, latitude_qc, longitude_qc)
        qc_lon = depth_average_data['lon'].isel(x=x_index, y=y_index).values
        qc_lat = depth_average_data['lat'].isel(x=x_index, y=y_index).values
        circle = Circle((qc_lon, qc_lat), radius=0.25, edgecolor='purple', facecolor='none', linewidth=2, transform=ccrs.PlateCarree(), zorder=95)
        ax.add_patch(circle)
    
    # FEATURES
    ax.add_feature(cfeature.GSHHSFeature(scale='full'), edgecolor="black", facecolor="tan", linewidth=0.25, zorder=90)
    ax.add_feature(cfeature.RIVERS, edgecolor="steelblue", linewidth=0.25, zorder=90)
    ax.add_feature(cfeature.LAKES, edgecolor="black", facecolor="lightsteelblue", linewidth=0.25, zorder=90)
    ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=0.25, zorder=90)
    
    plot_bathymetry(ax, config, depth_average_data, isobath1=-100, isobath2=-1000, show_legend=False)
    bathymetry_legend = ax.get_legend()
    if bathymetry_legend:
        bathymetry_legend.get_frame().set_alpha(0.5)
        ax.add_artist(bathymetry_legend)

    # COLORBAR
    cbar = fig.colorbar(contourf, orientation='vertical', extend=extend)
    format_colorbar(fig, ax, cbar)
    cbar.set_label('Depth Averaged Current Magnitude (m/s)', labelpad=10)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{tick:.1f}" for tick in ticks])

    # TITLES
    title_text = f"Depth Averaged Currents - Depth Range: {config['max_depth']}m"
    format_titles(ax, fig, config, model_data=depth_average_data, model_name=model_name, title=title_text)

    # SAVE & CLOSE
    file_datetime = format_save_datetime(depth_average_data)
    fig_filename = f"GGS_{config['glider_name']}_{model_name}_Magnitude_{config['max_depth']}m_{file_datetime}.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # LOGGING
    end_time = print_endtime()
    print_runtime(start_time, end_time)

### FUNCTION:
def GGS_plot_magnitude_dual(config, directory, depth_average_data_1, depth_average_data_2, latitude_qc, longitude_qc, model_name_1=None, model_name_2=None, density=2, show_gliders=False, show_route=False, show_qc=False, manual_extent=None):
    
    '''
    Plot the depth-averaged current fields from two datasets side by side.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - depth_average_data1, depth_average_data2 (xarray.Dataset): Datasets with the computed variables and layer information for side-by-side comparison.
    - qc_latitude (float): Latitude of the QC sample point.
    - qc_longitude (float): Longitude of the QC sample point.
    - model_name (str, optional): Name of the model.
    - density (int, optional): Density of the streamplot.
    - show_route (bool, optional): Show the route on the plot.
    - show_qc (bool, optional): Show the QC sample point.
    - manual_extent (list or str, optional): Manual specification of plot extent as [west, east, south, north], or "config" to use extent from config.

    Returns:
    - None
    '''

    # INITIALIZATION
    if model_name_1 is None:
        model_name_1 = "Unknown"
    else:
        model_name_1 = model_name_1.upper()
    if model_name_2 is None:
        model_name_2 = "Unknown"
    else:
        model_name_2 = model_name_2.upper()
    
    print(f"\n### CREATING DUAL ({model_name_1}, {model_name_2}) MAGNITUDE PLOT ###\n")
    start_time = dt.datetime.now()

    # PLOTTING FUNCTION
    def plot_magnitude(ax, config, depth_average_data, latitude_qc, longitude_qc, model_name, density, show_gliders, show_route, show_qc, manual_extent):
        
        # DATA EXTRACTION
        longitude = depth_average_data.lon.values.squeeze()
        latitude = depth_average_data.lat.values.squeeze()
        u_depth_avg = depth_average_data['u_depth_avg'].values.squeeze()
        v_depth_avg = depth_average_data['v_depth_avg'].values.squeeze()
        mag_depth_avg = depth_average_data['mag_depth_avg'].values.squeeze()
        
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
        plot_formatted_ticks(ax, map_extent[:2], map_extent[2:], proj=ccrs.PlateCarree(), fontsize=10, label_left=True, label_right=False, label_bottom=True, label_top=False, gridlines=True)

        # PLOT ELEMENTS
        levels, ticks, extend = plot_contour_cbar(mag_depth_avg, max_levels=10, extend_max=True)
        contourf = ax.contourf(longitude, latitude, mag_depth_avg, levels=levels, cmap=cmo.speed, transform=ccrs.PlateCarree(), zorder=10, extend=extend)
        streamplot = ax.streamplot(longitude, latitude, u_depth_avg, v_depth_avg, density=density, linewidth=0.5, color='black', transform=ccrs.PlateCarree(), zorder=10)

        # GLIDERS
        if show_gliders:
            glider_dataframes = acquire_gliders(extent=map_extent, target_date=dt.datetime.now(), date_delta=dt.timedelta(days=1), requested_variables=["time", "longitude", "latitude", "profile_id", "depth"], print_vars=False, target="all", request_timeout=5, enable_parallel=False)
            plot_add_gliders(ax, glider_dataframes, legend=True)
            glider_legend = ax.get_legend()
            if glider_legend:
                glider_legend.get_frame().set_alpha(0.5)
                ax.add_artist(glider_legend)
        
        # ROUTE
        if show_route:
            lats, lons = zip(*config["GPS_coords"])
            ax.plot(lons, lats, 'w-', transform=ccrs.PlateCarree(), linewidth=2.5, zorder=91)
            ax.plot(lons, lats, 'k', transform=ccrs.PlateCarree(), linewidth=1.0, linestyle='--', alpha=0.6, zorder=92)
            
            start_coords = config["GPS_coords"][0]
            end_coords = config["GPS_coords"][-1]
            ax.scatter(*start_coords[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
            for GPS_coord in config["GPS_coords"][1:-1]:
                ax.scatter(*GPS_coord[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
            ax.scatter(*end_coords[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
        
        # QUALITY CONTROL
        if show_qc:
            (y_index, x_index), (lat_index, lon_index) = calculate_gridpoint(depth_average_data, latitude_qc, longitude_qc)
            qc_lon = depth_average_data['lon'].isel(x=x_index, y=y_index).values
            qc_lat = depth_average_data['lat'].isel(x=x_index, y=y_index).values
            circle = Circle((qc_lon, qc_lat), radius=0.25, edgecolor='purple', facecolor='none', linewidth=2, transform=ccrs.PlateCarree(), zorder=95)
            ax.add_patch(circle)
        
        # FEATURES
        ax.add_feature(cfeature.GSHHSFeature(scale='full'), edgecolor="black", facecolor="tan", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.RIVERS, edgecolor="steelblue", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.LAKES, edgecolor="black", facecolor="lightsteelblue", linewidth=0.25, zorder=90)
        ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=0.25, zorder=90)
        
        plot_bathymetry(ax, config, depth_average_data, isobath1=-100, isobath2=-1000, show_legend=False)
        bathymetry_legend = ax.get_legend()
        if bathymetry_legend:
            bathymetry_legend.get_frame().set_alpha(0.5)
            ax.add_artist(bathymetry_legend)

        # COLORBAR
        cbar = fig.colorbar(contourf, orientation='vertical', extend=extend)
        format_colorbar(fig, ax, cbar)
        cbar.set_label('Depth Averaged Current Magnitude (m/s)', labelpad=10)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f"{tick:.1f}" for tick in ticks])

    # PLOT SETUP
    fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(20, 10))

    # PLOT EXECUTION
    plot_magnitude(ax1, config, depth_average_data_1, latitude_qc, longitude_qc, model_name_1, density, show_gliders, show_route, show_qc, manual_extent)
    plot_magnitude(ax2, config, depth_average_data_2, latitude_qc, longitude_qc, model_name_2, density, show_gliders, show_route, show_qc, manual_extent)

    # TITLES
    ax1.set_title(f"{model_name_1}", fontsize=14, fontweight='bold', pad=20)
    ax2.set_title(f"{model_name_2}", fontsize=14, fontweight='bold', pad=20)
    title_text = f"Depth Averaged Currents - Depth Range: {config['max_depth']}m"
    model_name_dual = f"{model_name_1} vs. {model_name_2}"
    format_titles(ax1, fig, config, model_data=depth_average_data_1, model_name=model_name_dual, title=title_text)
    
    # SAVE & CLOSE
    file_datetime = format_save_datetime(depth_average_data_1)
    fig_filename = f"GGS_{config['glider_name']}_Magnitude_Dual_{config['max_depth']}m_{file_datetime}.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # LOGGING
    end_time = dt.datetime.now()
    print(f"Runtime: {end_time - start_time}")

### FUNCTION:
def GGS_plot_threshold(config, directory, depth_average_data, latitude_qc, longitude_qc, model_name=None, density=2, mag1=0.0, mag2=0.2, mag3=0.3, mag4=0.4, mag5=0.5, show_gliders=False, show_route=False, show_qc=False, manual_extent=None):
    
    '''
    Plot the depth-averaged current threshold zones.
    
    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - model_data (xarray.Dataset): Ocean model dataset.
    - depth_average_data (xarray.Dataset): Dataset with the computed variables and layer information.
    - qc_latitude (float): Latitude of the QC sample point.
    - qc_longitude (float): Longitude of the QC sample point.
    - model_name (str): Name of the model.
        - default: 'None'
    - mag1 (float): First threshold mag_depth_avg.
        - default: 0.0
    - mag2 (float): Second threshold mag_depth_avg.
        - default: 0.2
    - mag3 (float): Third threshold mag_depth_avg.
        - default: 0.3
    - mag4 (float): Fourth threshold mag_depth_avg.
        - default: 0.4
    - mag5 (float): Fifth threshold mag_depth_avg.
        - default: 0.5
    - show_route (bool): Show the route on the plot.
        - default: 'False'
    - show_qc (bool): Whow the QC sample point.
        - default: 'False'
    - manual_extent (list or str): Optional manual specification of plot extent as [west, east, south, north], or "config" to use extent from config.
        - default: 'None'

    Returns:
    - None
    '''

    # INITIALIZATION
    if model_name is None:
        model_name = "Unknown"
    else:
        model_name = model_name.upper()
    
    print(f"\n### CREATING {model_name} THRESHOLD PLOT ###\n")
    start_time = print_starttime()

    # DATA EXTRACTION
    longitude = depth_average_data.lon.values.squeeze()
    latitude = depth_average_data.lat.values.squeeze()
    u_depth_avg = depth_average_data['u_depth_avg'].values.squeeze()
    v_depth_avg = depth_average_data['v_depth_avg'].values.squeeze()
    mag_depth_avg = depth_average_data['mag_depth_avg'].values.squeeze()

    # PLOT SETUP
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

    # EXTENT SETUP
    if manual_extent == "config":
        map_extent = [config['extent'][0][1], config['extent'][1][1], config['extent'][0][0], config['extent'][1][0]]
    elif isinstance(manual_extent, list) and len(manual_extent) == 4:
        map_extent = manual_extent
    else:
        data_extent_lon = [float(depth_average_data.lon.min()), float(depth_average_data.lon.max())]
        data_extent_lat = [float(depth_average_data.lat.min()), float(depth_average_data.lat.max())]
        map_extent = data_extent_lon + data_extent_lat
    ax.set_extent(map_extent, crs=ccrs.PlateCarree())
    plot_formatted_ticks(ax, map_extent[:2], map_extent[2:], proj=ccrs.PlateCarree(), fontsize=10, label_left=True, label_right=False, label_bottom=True, label_top=False, gridlines=True)
    
    # PLOT ELEMENTS
    levels = [mag1, mag2, mag3, mag4, mag5, np.nanmax(mag_depth_avg)]
    colors = ['none', 'yellow', 'orange', 'orangered', 'maroon']
    contourf = ax.contourf(longitude, latitude, mag_depth_avg, levels=levels, colors=colors, extend='both', transform=ccrs.PlateCarree(), zorder=10)
    streamplot = ax.streamplot(longitude, latitude, u_depth_avg, v_depth_avg, transform=ccrs.PlateCarree(), density=density, linewidth=0.5, color='dimgrey', zorder=10)
    streamplot.lines.set_alpha(1.0)
    plot_threshold_legend(ax, mag2, mag3, mag4, mag5)
    threshold_legend = ax.get_legend()
    if threshold_legend:
        threshold_legend.get_frame().set_alpha(0.5)
        ax.add_artist(threshold_legend)

    # GLIDERS
    if show_gliders:
        glider_dataframes = acquire_gliders(extent=map_extent, target_date=dt.datetime.now(), date_delta=dt.timedelta(days=1), requested_variables=["time", "longitude", "latitude", "profile_id", "depth"], print_vars=False, target="all", request_timeout=5, enable_parallel=False)
        plot_add_gliders(ax, glider_dataframes, legend=True)        
        glider_legend = ax.get_legend()
        if glider_legend:
            glider_legend.get_frame().set_alpha(0.5)
            ax.add_artist(glider_legend)

    # ROUTE
    if show_route:
        lats, lons = zip(*config["GPS_coords"])
        ax.plot(lons, lats, 'w-', transform=ccrs.PlateCarree(), linewidth=2.5, zorder=91)
        ax.plot(lons, lats, 'k', transform=ccrs.PlateCarree(), linewidth=1.0, linestyle='--', alpha=0.6, zorder=92)
        
        start_coords = config["GPS_coords"][0]
        end_coords = config["GPS_coords"][-1]
        ax.scatter(*start_coords[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
        for GPS_coord in config["GPS_coords"][1:-1]:
            ax.scatter(*GPS_coord[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
        ax.scatter(*end_coords[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)

    # QUALITY CONTROL
    if show_qc:
        (y_index, x_index), (lat_index, lon_index) = calculate_gridpoint(depth_average_data, latitude_qc, longitude_qc)
        qc_lon = depth_average_data['lon'].isel(x=x_index, y=y_index).values
        qc_lat = depth_average_data['lat'].isel(x=x_index, y=y_index).values
        circle = Circle((qc_lon, qc_lat), radius=0.25, edgecolor='white', facecolor='none', linewidth=2, transform=ccrs.PlateCarree(), zorder=95)
        ax.add_patch(circle)
    
    # FEATURES
    ax.add_feature(cfeature.GSHHSFeature(scale='full'), edgecolor="black", facecolor="tan", linewidth=0.25, zorder=90)
    ax.add_feature(cfeature.RIVERS, edgecolor="steelblue", linewidth=0.25, zorder=90)
    ax.add_feature(cfeature.LAKES, edgecolor="black", facecolor="lightsteelblue", linewidth=0.25, zorder=90)
    ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=0.25, zorder=90)
    
    plot_bathymetry(ax, config, depth_average_data, isobath1=-100, isobath2=-1000, show_legend=True)
    bathymetry_legend = ax.get_legend()
    if bathymetry_legend:
        bathymetry_legend.get_frame().set_alpha(0.5)
        ax.add_artist(bathymetry_legend)

    # TITLES
    title_text = f"Depth Averaged Threshold Zones - Depth Range: {config['max_depth']}m"
    format_titles(ax, fig, config, model_data=depth_average_data, model_name=model_name, title=title_text)

    # SAVE & CLOSE
    file_datetime = format_save_datetime(depth_average_data)
    fig_filename = f"GGS_{config['glider_name']}_{model_name}_Threshold_{config['max_depth']}m_{file_datetime}.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    # LOGGING
    end_time = print_endtime()
    print_runtime(start_time, end_time)
