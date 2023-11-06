# =========================
# X - Imports
# =========================

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import cmocean
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os
import plotly.graph_objects as go
import plotly.io as pio

from sup_functions import DD_to_DDM

# =========================
# GGS Plotting
# =========================

### FUNCTION:
def GGS_plot_route(config, waypoints, rtofs):

    """
    Plot the glider's mission route along with the depth-averaged currents and the advected path.

    Args:
    - rtofs (RTOFS): The RTOFS instance.
    - waypoints (list): A list of coordinates representing the mission route.
    - config (dict): The GGS configuration dictionary.
    - advected_route (list): A list of coordinates representing the advected route.
    """

    longitudes = rtofs.data.lon.values
    latitudes = rtofs.data.lat.values
    u_avg, v_avg = rtofs.rtofs_depth_avg(config)
    magnitude = rtofs.rtofs_current_magnitude(u_avg, v_avg)
    
    print(longitudes.shape)
    print(latitudes.shape)
    print(u_avg.shape)
    print(v_avg.shape)
    print(magnitude.shape)

    u_avg_2d = u_avg[0]
    v_avg_2d = v_avg[0]
    magnitude_2d = magnitude[0]

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    contour = ax.contourf(longitudes, latitudes, magnitude_2d, cmap=cmocean.cm.speed, transform=ccrs.PlateCarree())
    ax.streamplot(longitudes, latitudes, u_avg_2d, v_avg_2d, color='black', transform=ccrs.PlateCarree())
    
    lats, lons = zip(*waypoints)
    ax.plot(lons, lats, 'w-', transform=ccrs.PlateCarree(), linewidth=2.5, zorder=1)
    ax.plot(lons, lats, 'k', transform=ccrs.PlateCarree(), linewidth=1.0, linestyle='--', alpha=0.6, zorder=1)
    
    start_coords = config["waypoints"][0]
    end_coords = config["waypoints"][-1]
    ax.scatter(*start_coords[::-1], color='green', s=100, transform=ccrs.PlateCarree(), zorder=2)
    for waypoint in config["waypoints"][1:-1]:
        ax.scatter(*waypoint[::-1], color='blue', s=100, transform=ccrs.PlateCarree(), zorder=2)
    ax.scatter(*end_coords[::-1], color='red', s=100, transform=ccrs.PlateCarree(), zorder=2)
            
    ax.set_extent([longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()], crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, edgecolor='black', zorder=0)
    ax.add_feature(cfeature.OCEAN, zorder=0, edgecolor='k', facecolor='lightblue')
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    ax.coastlines()
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    cbar_ax = fig.add_axes([box.x0 + box.width + 0.01, box.y0, 0.03, box.height])
    fig.colorbar(contour, cax=cbar_ax, orientation='vertical', label='Depth Averaged Current Magnitude (m/s)')

    ax.set_xticks(np.linspace(longitudes.min(), longitudes.max(), 5), crs=ccrs.PlateCarree())
    ax.set_yticks(np.linspace(latitudes.min(), latitudes.max(), 5), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.set_title(f'{config["glider_name"]} Mission Route Along Depth-Averaged Currents', pad=20)
    
    mission_directory = os.path.join(os.path.expanduser("~"), "Downloads", f"GGS_{config['glider_name']}")
    fig_filename = f"{config['glider_name']}_depth_avg_currents.png"
    fig_path = os.path.join(mission_directory, fig_filename)
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')

### FUNCTION:
def GGS_plot_gauge(config, directory, analysis_results):

    """
    Plot a gauge showing the remaining battery capacity.

    Args:
    - config (dict): The GGS configuration dictionary.
    - directory (str): The directory path where the configuration is saved.
    - analysis_results (list): List of dictionaries containing data for each leg of the route.
    """
    
    total_distance, total_time_seconds, total_time_hours, total_battery_drain = route_analysis_output(config, directory, analysis_results)
    battery_remaining = config["battery_capacity"] - total_battery_drain

    fig_gauge = go.Figure(go.Indicator(
        mode = "gauge+number+delta",
        value = battery_remaining,
        domain = {'x': [0, 1], 'y': [0, 1]},
        title = {'text': "Battery Remaining (Amperage)", 'font': {'size': 24}},
        delta = {'reference': config["battery_capacity"], 'increasing': {'color': "RebeccaPurple"}},
        gauge = {
            'axis': {'range': [0, config["battery_capacity"]], 'tickwidth': 1, 'tickcolor': "darkblue"},
            'bar': {'color': "yellow"},
            'bgcolor': "white",
            'borderwidth': 2,
            'bordercolor': "gray",
            'steps': [
                {'range': [0, 0.5*config["battery_capacity"]], 'color': 'lightgray'},
                {'range': [0.5*config["battery_capacity"], config["battery_capacity"]], 'color': 'gray'}]
        }))

    fig_gauge.update_layout(paper_bgcolor = "lavender", font = {'color': "darkblue", 'family': "Arial"})

    gauge_filename = f"{config['glider_name']}_battery_gauge.png"
    gauge_path = os.path.join(directory, gauge_filename)
    
    pio.write_image(fig_gauge, gauge_path, format='png')

# =========================
# Plotting Subfunctions
# =========================

### FUNCTION:
def DD_to_DMS(deg, coord_type='longitude'):

    '''
    Convert a decimal degree position to a degree-minute-second string.

    Args:
    - deg (float): A decimal degree position.
    - coord_type (str): A string indicating whether the position is a longitude or latitude.

    Returns:
    - str: A degree-minute-second string.
    '''

    degrees = int(deg)
    minutes = int((deg - degrees) * 60)
    seconds = (deg - degrees - minutes/60) * 3600.00
    
    direction = ''
    if degrees > 0:
        if coord_type == 'longitude':
            direction = 'E'
        elif coord_type == 'latitude':
            direction = 'N'
    elif degrees < 0:
        if coord_type == 'longitude':
            direction = 'W'
        elif coord_type == 'latitude':
            direction = 'S'
    
    return f"{abs(degrees)}°{abs(minutes)}'{abs(seconds):.2f}\"{direction}"

### FUNCTION:
def DD_to_DM(deg, coord_type='longitude'):
    '''
    Convert a decimal degree position to a degree-minute string.

    Args:
    - deg (float): A decimal degree position.
    - coord_type (str): A string indicating whether the position is a longitude or latitude.

    Returns:
    - str: A degree-minute string.
    '''

    degrees = int(deg)
    minutes = abs(int((deg - degrees) * 60))

    direction = ''
    if degrees > 0:
        if coord_type == 'longitude':
            direction = 'E'
        elif coord_type == 'latitude':
            direction = 'N'
    elif degrees < 0:
        if coord_type == 'longitude':
            direction = 'W'
        elif coord_type == 'latitude':
            direction = 'S'

    return f"{abs(degrees)}°{minutes}'{direction}"

### FUNCTION:
def calculate_ticks(ax, extent_lon, extent_lat):
    """
    Set the ticks for a map based on the extent without changing the map's bounds, ensuring at least three ticks.

    Args:
    - ax (matplotlib.axes): The axes to set the ticks for.
    - extent_lon (array-like): The longitude extent of the map.
    - extent_lat (array-like): The latitude extent of the map.

    Returns:
    - None
    """

    ### FUNCTION:
    def get_tick_interval(extent):

        '''
        Calculate the tick interval for a given extent.

        Args:
        - extent (array-like): The extent of the map.

        Returns:
        - float: The tick interval.
        '''

        min_extent, max_extent = np.min(extent), np.max(extent)
        range_extent = max_extent - min_extent

        intervals = [0.25, 0.5, 1, 2, 5, 10, 15, 30]
        for interval in intervals:
            num_ticks = int(range_extent / interval) + 1
            if 3 <= num_ticks <= 6:
                return interval
        return intervals[-1]

    lon_interval = get_tick_interval(extent_lon)
    lat_interval = get_tick_interval(extent_lat)

    lon_start = np.ceil(np.min(extent_lon) / lon_interval) * lon_interval
    lon_end = np.floor(np.max(extent_lon) / lon_interval) * lon_interval
    lat_start = np.ceil(np.min(extent_lat) / lat_interval) * lat_interval
    lat_end = np.floor(np.max(extent_lat) / lat_interval) * lat_interval

    lon_ticks = np.arange(lon_start, lon_end + lon_interval/2, lon_interval)
    lat_ticks = np.arange(lat_start, lat_end + lat_interval/2, lat_interval)

    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda val, pos: DD_to_DM(val, 'longitude')))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda val, pos: DD_to_DM(val, 'latitude')))
