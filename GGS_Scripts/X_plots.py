# =========================
# X - Imports
# =========================

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
import numpy as np
import os
from X_functions import calculate_gridpoint, plot_formatted_ticks, plot_contour_cbar, plot_bathymetry, datetime_filename, datetime_title

# =========================
# PLOTTING
# =========================

### FUNCTION:
def GGS_plot_currents(config, directory, model_data, depth_average_data, qc_latitude, qc_longitude, density=2, show_route=False, show_qc=False):
    
    '''
    Plot the depth-averaged current fields.
    Optionally display the mission route and/or quality control sample location.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - model_data (xarray.Dataset): Ocean model dataset.
    - depth_average_data (xarray.Dataset): Dataset with the computed variables and layer information.
    - qc_latitude (float): Latitude of the QC sample point.
    - qc_longitude (float): Longitude of the QC sample point.
    - density (int): Density of the streamplot.
        - default: '2'
    - show_route (bool): Show the route on the plot.
        - default: 'False'
        - if True, show the route.
        - if False, do not show the route.
    - show_qc (bool): Whow the QC sample point.
        - default: 'False'
        - if True, show the QC sample point.
        - if False, do not show the QC sample point.

    Returns:
    - None
    '''

    print("\n")
    print("### CREATING MAGNITUDE PLOT ###")
    print("\n")

    u_depth_avg = depth_average_data['u_depth_avg'].values
    v_depth_avg = depth_average_data['v_depth_avg'].values
    mag_depth_avg = depth_average_data['mag_depth_avg'].values
    
    data_lons = model_data.lon.values
    data_lats = model_data.lat.values

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

    data_extent_lon = [np.min(data_lons), np.max(data_lons)]
    data_extent_lat = [np.min(data_lats), np.max(data_lats)]
    ax.set_extent(data_extent_lon + data_extent_lat, crs=ccrs.PlateCarree())
    plot_formatted_ticks(ax, data_extent_lon, data_extent_lat, proj=ccrs.PlateCarree(), fontsize=10, label_left=True, label_right=False, label_bottom=True, label_top=False, gridlines=True)
    
    levels, ticks = plot_contour_cbar(mag_depth_avg, max_levels=10)
    contour = ax.contourf(data_lons, data_lats, mag_depth_avg, levels=levels, cmap=cmo.speed, transform=ccrs.PlateCarree(), zorder=10) # zorder = [1]
    ax.streamplot(data_lons, data_lats, u_depth_avg, v_depth_avg, color='black', transform=ccrs.PlateCarree(), density=density, linewidth=0.5, zorder=10) # zorder = [1]

    if show_route: # zorder = [9]
        lats, lons = zip(*config["GPS_coords"])
        ax.plot(lons, lats, 'w-', transform=ccrs.PlateCarree(), linewidth=2.5, zorder=91)
        ax.plot(lons, lats, 'k', transform=ccrs.PlateCarree(), linewidth=1.0, linestyle='--', alpha=0.6, zorder=92)
        
        start_coords = config["GPS_coords"][0]
        end_coords = config["GPS_coords"][-1]
        ax.scatter(*start_coords[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
        for GPS_coord in config["GPS_coords"][1:-1]:
            ax.scatter(*GPS_coord[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
        ax.scatter(*end_coords[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)

    if show_qc: # zorder = [9]
        (y_index, x_index), (lat_index, lon_index) = calculate_gridpoint(model_data, qc_latitude, qc_longitude)
        qc_lon = model_data['lon'].isel(x=x_index, y=y_index).values
        qc_lat = model_data['lat'].isel(x=x_index, y=y_index).values
        circle = Circle((qc_lon, qc_lat), radius=0.25, edgecolor='purple', facecolor='none', linewidth=2, transform=ccrs.PlateCarree(), zorder=95)
        ax.add_patch(circle)
    
    ax.add_feature(cfeature.GSHHSFeature(scale='full'), edgecolor="black", facecolor="tan", linewidth=0.25, zorder=90) # zorder = [9]
    ax.add_feature(cfeature.RIVERS, edgecolor="steelblue", linewidth=0.25, zorder=90) # zorder = [9]
    ax.add_feature(cfeature.LAKES, edgecolor="black", facecolor="lightsteelblue", linewidth=0.25, zorder=90) # zorder = [9]
    ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=0.25, zorder=90) # zorder = [9]
    
    plot_bathymetry(ax, config, model_data, isobath1=-100, isobath2=-1000, show_legend=False)

    box = ax.get_position()
    left_margin = box.x0
    bottom_margin = box.y0
    plot_width = box.width
    plot_height = box.height
    ax.set_position([left_margin, bottom_margin, plot_width, plot_height])

    cbar_ax = fig.add_axes([left_margin + plot_width + 0.02, bottom_margin, 0.03, plot_height])
    cbar = fig.colorbar(contour, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Depth Averaged Current mag_depth_avg (m/s)', labelpad=10)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{tick:.1f}" for tick in ticks])

    title_datetime = datetime_title(model_data)
    title_text = f"Depth Averaged Currents - Depth Range: {config['max_depth']}m"
    ax.set_title(title_text, fontsize=14, fontweight='bold', pad=25)
    subtitle_text = f"[RTOFS] {title_datetime} UTC"
    fig.text(0.5, 0.915, subtitle_text, fontsize=10, fontweight='bold', ha='center', va='center')
    suptitle_text = f"Generated by the Glider Guidance System (GGS) - {config['glider_name']}"
    fig.suptitle(suptitle_text, fontsize='smaller', fontweight='bold', x=0.5, y=0.01, ha='center', va='bottom', color='gray')
    
    filename_datetime = datetime_filename(model_data)
    fig_filename = f"GGS_{config['glider_name']}_Currents_{config['max_depth']}m.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')

    plt.close(fig)

    print("### MAGNITUDE PLOT SAVED ###")

### FUNCTION:
def GGS_plot_threshold(config, directory, model_data, depth_average_data, qc_latitude, qc_longitude, mag1=0.0, mag2=0.2, mag3=0.3, mag4=0.4, mag5=0.5, show_route=False, show_qc=False):
    
    '''
    Plots the depth-averaged current mag_depth_avg threshold zones for the currents data.
    Optionally display the mission route and/or quality control sample location.
    
    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - model_data (xarray.Dataset): Ocean model dataset.
    - depth_average_data (xarray.Dataset): Dataset with the computed variables and layer information.
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
    - qc_latitude (float): Latitude of the QC sample point.
    - qc_longitude (float): Longitude of the QC sample point.
    - show_route (bool): Show the route on the plot.
        - default: 'False'
        - if True, show the route.
        - if False, do not show the route.
    - show_qc (bool): Whow the QC sample point.
        - default: 'False'
        - if True, show the QC sample point.
        - if False, do not show the QC sample point.

    Returns:
    - None
    '''

    print("\n")
    print("### CREATING THRESHOLD PLOT ###")
    print("\n")

    u_depth_avg = depth_average_data['u_depth_avg'].values
    v_depth_avg = depth_average_data['v_depth_avg'].values
    mag_depth_avg = depth_average_data['mag_depth_avg'].values
    mag_depth_avg = np.nan_to_num(mag_depth_avg, nan=0.0, posinf=0.0, neginf=0.0)

    data_lons = model_data.lon.values
    data_lats = model_data.lat.values

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

    data_extent_lon = [np.min(data_lons), np.max(data_lons)]
    data_extent_lat = [np.min(data_lats), np.max(data_lats)]
    ax.set_extent(data_extent_lon + data_extent_lat, crs=ccrs.PlateCarree())
    plot_formatted_ticks(ax, data_extent_lon, data_extent_lat, proj=ccrs.PlateCarree(), fontsize=10, label_left=True, label_right=True, label_bottom=True, label_top=False, gridlines=True)
    
    levels = [mag1, mag2, mag3, mag4, mag5, np.max(mag_depth_avg)]
    colors = ['none', 'yellow', 'orange', 'orangered', 'maroon']
    
    contourf = ax.contourf(data_lons, data_lats, mag_depth_avg, levels=levels, colors=colors, extend='both', transform=ccrs.PlateCarree(), zorder=10) # zorder = [1]
    streamplot = ax.streamplot(data_lons, data_lats, u_depth_avg, v_depth_avg, color='dimgrey', transform=ccrs.PlateCarree(), density=2, linewidth=0.5, zorder=11) # zorder = [1]
    streamplot.lines.set_alpha(1.0)

    if show_route: # zorder = [9]
        lats, lons = zip(*config["GPS_coords"])
        ax.plot(lons, lats, 'w-', transform=ccrs.PlateCarree(), linewidth=2.5, zorder=91)
        ax.plot(lons, lats, 'k', transform=ccrs.PlateCarree(), linewidth=1.0, linestyle='--', alpha=0.6, zorder=92)
        
        start_coords = config["GPS_coords"][0]
        end_coords = config["GPS_coords"][-1]
        ax.scatter(*start_coords[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
        for GPS_coord in config["GPS_coords"][1:-1]:
            ax.scatter(*GPS_coord[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)
        ax.scatter(*end_coords[::-1], color='purple', s=100, transform=ccrs.PlateCarree(), zorder=93)

    if show_qc: # zorder = [9]
        (y_index, x_index), (lat_index, lon_index) = calculate_gridpoint(model_data, qc_latitude, qc_longitude)
        qc_lon = model_data['lon'].isel(x=x_index, y=y_index).values
        qc_lat = model_data['lat'].isel(x=x_index, y=y_index).values
        circle = Circle((qc_lon, qc_lat), radius=0.25, edgecolor='white', facecolor='none', linewidth=2, transform=ccrs.PlateCarree(), zorder=95)
        ax.add_patch(circle)
    
    ax.add_feature(cfeature.GSHHSFeature(scale='full'), edgecolor="black", facecolor="tan", linewidth=0.25, zorder=90) # zorder = [9]
    ax.add_feature(cfeature.RIVERS, edgecolor="steelblue", linewidth=0.25, zorder=90) # zorder = [9]
    ax.add_feature(cfeature.LAKES, edgecolor="black", facecolor="lightsteelblue", linewidth=0.25, zorder=90) # zorder = [9]
    ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=0.25, zorder=90) # zorder = [9]
    
    plot_bathymetry(ax, config, model_data, isobath1=-100, isobath2=-1000, show_legend=True)
    bathymetry_legend = ax.get_legend()

    patches = [
        mpatches.Patch(facecolor='yellow', label=f'{mag2} - {mag3} m/s'),
        mpatches.Patch(facecolor='orange', label=f'{mag3} - {mag4} m/s'),
        mpatches.Patch(facecolor='orangered', label=f'{mag4} - {mag5} m/s'),
        mpatches.Patch(facecolor='maroon', label=f'> {mag5} m/s')]
    legend = ax.legend(handles=patches, loc='upper right', facecolor='white', edgecolor='black', framealpha=0.75, fontsize='x-small')
    legend.set_zorder(1000) # zorder = [100]
    legend.get_title().set_color('black')
    for text in legend.get_texts():
        text.set_color('black')
    if bathymetry_legend:
        ax.add_artist(bathymetry_legend) # zorder = [100]

    title_datetime = datetime_title(model_data)
    title_text = f"Depth Averaged Threshold Zones - Depth Range: {config['max_depth']}m"
    ax.set_title(title_text, fontsize=14, fontweight='bold', pad=25)
    subtitle_text = f"[RTOFS] {title_datetime} UTC"
    fig.text(0.5, 0.915, subtitle_text, fontsize=10, fontweight='bold', ha='center', va='center')
    suptitle_text = f"Generated by the Glider Guidance System (GGS) - {config['glider_name']}"
    fig.suptitle(suptitle_text, fontsize='smaller', fontweight='bold', x=0.5, y=0.01, ha='center', va='bottom', color='gray')
    
    filename_datetime = datetime_filename(model_data)
    fig_filename = f"GGS_{config['glider_name']}_ThresholdZones_{config['max_depth']}m.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    
    plt.close(fig)

    print("### THRESHOLD PLOT SAVED ###")
