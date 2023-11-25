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

# =========================
# QUALITY CONTROL CHECKS
# =========================

### FUNCTION:
def qc_selection(rtofs_qc, num_points=10):

    '''
    Select random latitude and longitude coordinate points to sample from the given 'rtofs_qc' dataset.

    Args:
    - rtofs_qc (xarray.Dataset): RTOFS dataset
    - num_points (int): number of points to select
    
    Returns:
    - selected_coords (dict): dictionary of selected coordinates
    '''

    total_points = len(rtofs_qc.y) * len(rtofs_qc.x)
    num_points = min(num_points, total_points)

    y_indices = np.random.choice(rtofs_qc.y, size=num_points, replace=False)
    x_indices = np.random.choice(rtofs_qc.x, size=num_points, replace=False)
    
    selected_coords = {'y': y_indices, 'x': x_indices}
    
    return selected_coords

### FUNCTION:
def qc_plots(directory, rtofs_qc, selected_coords, u_avg, v_avg):
    
    '''
    Plot the quality control checks for the given 'rtofs_qc' dataset.

    Args:
    - rtofs_qc (xarray.Dataset): RTOFS dataset
    - selected_coords (dict): dictionary of selected coordinates
    - u_avg (xarray.DataArray): depth-averaged u currents
    - v_avg (xarray.DataArray): depth-averaged v currents

    Returns:
    - None
    '''
    
    for i, (y_index, x_index) in enumerate(zip(selected_coords['y'], selected_coords['x'])):
        lat = rtofs_qc.lat.sel(y=y_index, x=x_index).values
        lon = rtofs_qc.lon.sel(y=y_index, x=x_index).values

        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

        for j, variable in enumerate(['u', 'v']):
            coord_data = rtofs_qc.sel(y=y_index, x=x_index)
            layer_thickness = coord_data['depth'].diff(dim='depth', label='upper')
            n1_layer_thickness = layer_thickness.isel(depth=0)
            layer_thickness = xr.concat([n1_layer_thickness, layer_thickness], dim='depth')
            layer_thickness['depth'].values[0] = 0
            avg_value = u_avg.sel(y=y_index, x=x_index) if variable == 'u' else v_avg.sel(y=y_index, x=x_index)

            cmap = cmo.speed
            norm = mcolors.Normalize(vmin=coord_data[variable].min(), vmax=coord_data[variable].max())

            trendline_x = []
            trendline_y = []

            for depth in range(len(coord_data['depth'])):
                magnitude = coord_data[variable].isel(depth=depth).values
                depth_value = coord_data['depth'].isel(depth=depth).values

                color = cmap(norm(magnitude))
                axs[j].bar(magnitude, height=layer_thickness.isel(depth=depth), bottom=depth_value, color=color, edgecolor='black')

                trendline_x.append(magnitude)
                trendline_y.append(depth_value + layer_thickness.isel(depth=depth).values / 2)

            axs[j].plot(trendline_x, trendline_y, color='red', marker='o', linestyle='-', linewidth=2)

            axs[j].axvline(x=avg_value, color='blue', label='Depth-Averaged Value')
            axs[j].set_xlabel(f'{variable}-Magnitude (m/s)')
            axs[j].set_ylabel('Depth (m)')
            axs[j].set_title(f'Avg {variable.upper()} Mag: {avg_value.values:.3f}')
            axs[j].invert_yaxis()

            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            fig.colorbar(sm, ax=axs[j], orientation='vertical', label=f'{variable.upper()} Magnitude (m/s)')

        plt.suptitle(f'Quality Control Plots at: Lat: {lat:.3f}, Lon: {lon:.3f}', fontsize=14)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust the layout

        plt.savefig(os.path.join(directory, f'GGS_QualityControl_Lat{lat:.3f}_Lon{lon:.3f}_{i + 1}.png'))
        plt.close(fig)

# =========================
# X - MAIN
# =========================
# qc_coordinates = qc_selection(rtofs_data.rtofs_qc)
# qc_plots(directory, rtofs_qc, qc_coordinates, u_avg, v_avg)
# =========================