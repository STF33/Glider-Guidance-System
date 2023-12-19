# =========================
# X - Imports
# =========================

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from X_functions import calculate_gridpoint, set_map_ticks, calculate_cbar_ticks, add_bathymetry

# =========================
# PLOTTING
# =========================

### FUNCTION:
def GGS_plot_currents(config, directory, waypoints, model_data, currents_data, qc_latitude, qc_longitude, density=2, extent='data', map_lons=[0, 0], map_lats=[0, 0], show_route=False, show_qc=False):
    
    '''
    Plot the depth-averaged current fields.
    Optionally display the mission route and/or quality control sample location.

    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - waypoints (list): Glider Guidance System mission waypoints.
    - model_data (xarray.Dataset): Ocean model dataset.
    - currents_data (xarray.Dataset): Dataset with the computed variables and layer information.
    - qc_latitude (float): Latitude of the QC sample point.
    - qc_longitude (float): Longitude of the QC sample point.
    - density (int): Density of the streamplot.
        - default: '2'
    - extent (str): Extent of the plot.
        - if 'map', use the map_lons and map_lats.
        - if 'data', use the model_data extent.
    - map_lons (list): Longitude bounds of the map, if extent='map'.
    - map_lats (list): Latitude bounds of the map, if extent='map'.
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

    u_avg = currents_data['u_avg'].values
    v_avg = currents_data['v_avg'].values
    magnitude = currents_data['magnitude_avg'].values
    
    data_lons = model_data.lon.values
    data_lats = model_data.lat.values

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

    if extent == 'map':
        ax.set_extent([min(map_lons), max(map_lons), min(map_lats), max(map_lats)], crs=ccrs.PlateCarree())
        set_map_ticks(ax, map_lons, map_lats)
    elif extent == 'data':
        ax.set_extent([np.min(data_lons), np.max(data_lons), np.min(data_lats), np.max(data_lats)], crs=ccrs.PlateCarree())
        set_map_ticks(ax, data_lons, data_lats)
    else:
        raise ValueError("Invalid extent option. Use 'map' or 'data'.")

    ticks = calculate_cbar_ticks(magnitude)
    contour = ax.contourf(data_lons, data_lats, magnitude, levels=ticks, cmap=cmo.speed, transform=ccrs.PlateCarree(), zorder=10) # zorder = [1]
    ax.streamplot(data_lons, data_lats, u_avg, v_avg, color='black', transform=ccrs.PlateCarree(), density=density, zorder=10) # zorder = [1]

    if show_route: # zorder = [2]
        lats, lons = zip(*waypoints)
        ax.plot(lons, lats, 'w-', transform=ccrs.PlateCarree(), linewidth=2.5, zorder=21)
        ax.plot(lons, lats, 'k', transform=ccrs.PlateCarree(), linewidth=1.0, linestyle='--', alpha=0.6, zorder=22)
        
        start_coords = config["waypoints"][0]
        end_coords = config["waypoints"][-1]
        ax.scatter(*start_coords[::-1], color='green', s=100, transform=ccrs.PlateCarree(), zorder=23)
        for waypoint in config["waypoints"][1:-1]:
            ax.scatter(*waypoint[::-1], color='blue', s=100, transform=ccrs.PlateCarree(), zorder=23)
        ax.scatter(*end_coords[::-1], color='red', s=100, transform=ccrs.PlateCarree(), zorder=23)

    if show_qc: # zorder = [3]
        (y_index, x_index), (lat_index, lon_index) = calculate_gridpoint(model_data, qc_latitude, qc_longitude)
        qc_lon = model_data['lon'].isel(x=x_index, y=y_index).values
        qc_lat = model_data['lat'].isel(x=x_index, y=y_index).values
        ax.scatter(qc_lon, qc_lat, color='red', s=100, transform=ccrs.PlateCarree(), zorder=35)

    add_bathymetry(ax, model_data, isobath_levels=[-100, -1000])

    ax.add_feature(cfeature.GSHHSFeature(scale='full'), edgecolor="black", facecolor="tan", zorder=90) # zorder = [9]
    ax.add_feature(cfeature.LAND, edgecolor="black", facecolor="tan", zorder=90) # zorder = [9]
    ax.add_feature(cfeature.BORDERS, edgecolor="black", zorder=90) # zorder = [9]
    ax.add_feature(cfeature.RIVERS, edgecolor="steelblue", zorder=90) # zorder = [9]
    ax.add_feature(cfeature.LAKES, edgecolor="black", facecolor="steelblue", zorder=90) # zorder = [9]

    box = ax.get_position()

    left_margin = box.x0
    bottom_margin = box.y0
    plot_width = box.width
    plot_height = box.height

    ax.set_position([left_margin, bottom_margin, plot_width, plot_height])

    cbar_ax = fig.add_axes([left_margin + plot_width + 0.02, bottom_margin, 0.03, plot_height])
    cbar = fig.colorbar(contour, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Depth-Averaged Current Magnitude (m/s)', labelpad=10)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{tick:.1f}" for tick in ticks])

    title_text = f"{config['glider_name']} Mission - Depth-Averaged Currents"
    ax.set_title(title_text, pad=25)
    subtitle_text = f"Depth Range: {config['max_depth']}m"
    fig.text(0.5, 0.915, subtitle_text, ha='center', va='center', fontsize=10)
    suptitle_text = f"Generated by the Glider Guidance System (GGS) on {datetime.utcnow().strftime('%m-%d-%Y')} at {datetime.utcnow().strftime('%H:%M')} UTC"
    fig.suptitle(suptitle_text, fontsize='smaller', x=0.5, y=0.01, ha='center', va='bottom', color='gray')
    
    fig_filename = f"{config['glider_name']}_{config['max_depth']}m_currents.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')

    plt.close(fig)

### FUNCTION:
def GGS_plot_threshold(config, directory, waypoints, model_data, currents_data, qc_latitude, qc_longitude, mag1=0.25, mag2=0.50, mag3=0.75, extent='data', map_lons=[0, 0], map_lats=[0, 0], show_route=False, show_qc=False):
    
    '''
    Plots the depth-averaged current magnitude threshold zones for the currents data.
    Optionally display the mission route and/or quality control sample location.
    
    Args:
    - config (dict): Glider Guidance System mission configuration.
    - directory (str): Glider Guidance System mission directory.
    - waypoints (list): Glider Guidance System mission waypoints.
    - model_data (xarray.Dataset): Ocean model dataset.
    - currents_data (xarray.Dataset): Dataset with the computed variables and layer information.
    - mag1 (float): First threshold magnitude.
        - default: 0.25
    - mag2 (float): Second threshold magnitude.
        - default: 0.50
    - mag3 (float): Third threshold magnitude.
        - default: 0.75
    - qc_latitude (float): Latitude of the QC sample point.
    - qc_longitude (float): Longitude of the QC sample point.
    - extent (str): Extent of the plot.
        - if 'map', use the map_lons and map_lats.
        - if 'data', use the model_data extent.
    - map_lons (list): Longitude bounds of the map, if extent='map'.
    - map_lats (list): Latitude bounds of the map, if extent='map'.
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

    u_avg = currents_data['u_avg'].values
    v_avg = currents_data['v_avg'].values
    magnitude = currents_data['magnitude_avg'].values
    magnitude = np.nan_to_num(magnitude, nan=0.0, posinf=0.0, neginf=0.0)

    data_lons = model_data.lon.values
    data_lats = model_data.lat.values

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

    if extent == 'map':
        ax.set_extent([min(map_lons), max(map_lons), min(map_lats), max(map_lats)], crs=ccrs.PlateCarree())
        set_map_ticks(ax, map_lons, map_lats)
    elif extent == 'data':
        ax.set_extent([np.min(data_lons), np.max(data_lons), np.min(data_lats), np.max(data_lats)], crs=ccrs.PlateCarree())
        set_map_ticks(ax, data_lons, data_lats)
    else:
        raise ValueError("Invalid extent option. Use 'map' or 'data'.")

    levels = [0, mag1, mag2, mag3, np.max(magnitude)]
    colors = ['none', 'yellow', 'orange', 'firebrick']
    
    contourf = ax.contourf(data_lons, data_lats, magnitude, levels=levels, colors=colors, extend='both', transform=ccrs.PlateCarree(), zorder=10) # zorder = [1]
    ax.streamplot(data_lons, data_lats, u_avg, v_avg, color='black', transform=ccrs.PlateCarree(), density=2, zorder=11) # zorder = [1]
    
    if show_route: # zorder = [2]
        lats, lons = zip(*waypoints)
        ax.plot(lons, lats, 'w-', transform=ccrs.PlateCarree(), linewidth=2.5, zorder=21)
        ax.plot(lons, lats, 'k', transform=ccrs.PlateCarree(), linewidth=1.0, linestyle='--', alpha=0.6, zorder=22)
        
        start_coords = config["waypoints"][0]
        end_coords = config["waypoints"][-1]
        ax.scatter(*start_coords[::-1], color='green', s=100, transform=ccrs.PlateCarree(), zorder=23)
        for waypoint in config["waypoints"][1:-1]:
            ax.scatter(*waypoint[::-1], color='blue', s=100, transform=ccrs.PlateCarree(), zorder=23)
        ax.scatter(*end_coords[::-1], color='red', s=100, transform=ccrs.PlateCarree(), zorder=23)

    if show_qc: # zorder = [3]
        (y_index, x_index), (lat_index, lon_index) = calculate_gridpoint(model_data, qc_latitude, qc_longitude)
        qc_lon = model_data['lon'].isel(x=x_index, y=y_index).values
        qc_lat = model_data['lat'].isel(x=x_index, y=y_index).values
        ax.scatter(qc_lon, qc_lat, color='red', s=100, transform=ccrs.PlateCarree(), zorder=35)

    add_bathymetry(ax, model_data, isobath_levels=[-100, -1000])

    ax.add_feature(cfeature.GSHHSFeature(scale='full'), edgecolor="black", facecolor="tan", zorder=90) # zorder = [9]
    ax.add_feature(cfeature.LAND, edgecolor="black", facecolor="tan", zorder=90) # zorder = [9]
    ax.add_feature(cfeature.BORDERS, edgecolor="black", zorder=90) # zorder = [9]
    ax.add_feature(cfeature.RIVERS, edgecolor="steelblue", zorder=90) # zorder = [9]
    ax.add_feature(cfeature.LAKES, edgecolor="black", facecolor="steelblue", zorder=90) # zorder = [9]
    
    patches = [
        mpatches.Patch(color='yellow', label=f'{mag1} - {mag2} m/s'),
        mpatches.Patch(color='orange', label=f'{mag2} - {mag3} m/s'),
        mpatches.Patch(color='firebrick', label=f'>= {mag3} m/s')]
    legend = ax.legend(handles=patches, loc='upper right', title='Current Magnitude', facecolor='grey', edgecolor='black', framealpha=0.75,)
    legend.set_zorder(100) # zorder = [10]
    legend.get_title().set_color('white')
    for text in legend.get_texts():
        text.set_color('white')

    title_text = f"{config['glider_name']} Mission - Threshold Zones"
    ax.set_title(title_text, pad=25)
    subtitle_text = f"Depth Range: {config['max_depth']}m"
    fig.text(0.5, 0.915, subtitle_text, ha='center', va='center', fontsize=10)
    suptitle_text = f"Generated by the Glider Guidance System (GGS) on {datetime.utcnow().strftime('%m-%d-%Y')} at {datetime.utcnow().strftime('%H:%M')} UTC"
    fig.suptitle(suptitle_text, fontsize='smaller', x=0.5, y=0.01, ha='center', va='bottom', color='gray')
    
    fig_filename = f"{config['glider_name']}_{config['max_depth']}m_threshold_zones.png"
    fig_path = os.path.join(directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')

    plt.close(fig)
    