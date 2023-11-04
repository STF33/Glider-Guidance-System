# Imports

# from sup_models import rtofs, RTOFS
# from sup_functions import depth_averaged_current_manual
import cool_maps.plot as cplt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

# Code

extent = [-90, -82, 17, 24] # yucatan

url = "C:/Users/salfr/Downloads/rtofs_glo_3dz_f006_6hrly_hvr_US_east.nc" # pull netcdf from boardwalk or noaa
ds = xr.open_dataset(url)#.set_coords(['lon', 'lat'])
ds.attrs['model'] = 'RTOFS'

ds = ds.sel(Depth=slice(0,1000))

def depth_averaged_current_manual(ds):
    u_currents = ds['u']
    v_currents = ds['v']
    depths = ds['Depth']

    layer_thicknesses = np.diff(depths)
    layer_thicknesses = np.insert(layer_thicknesses, 0, layer_thicknesses[0]) # Add the first thickness to the start

    depth_weighted_u_currents = u_currents * layer_thicknesses[:, np.newaxis, np.newaxis]
    depth_weighted_v_currents = v_currents * layer_thicknesses[:, np.newaxis, np.newaxis]

    u_avg = depth_weighted_u_currents.sum(dim='Depth') / np.sum(layer_thicknesses)
    v_avg = depth_weighted_v_currents.sum(dim='Depth') / np.sum(layer_thicknesses)

    return u_avg, v_avg

u_avg, v_avg = depth_averaged_current_manual(ds)

ds['u_avg'] = (('MT', 'Y', 'X'), u_avg.data)
ds['v_avg'] = (('MT', 'Y', 'X'), v_avg.data)

import cartopy.crs as ccrs
import cmocean
from oceans.ocfis import uv2spdir

ds = ds.squeeze()

_, speed = uv2spdir(ds.u_avg, ds.v_avg)

fig, ax = cplt.create(extent, bathymetry=False)
# fig, ax = plt.subplots(figsize=(10, 10))
ax.quiver(ds.Longitude.data, ds.Latitude.data, ds.u_avg.data, ds.v_avg.data, color='black', transform=ccrs.PlateCarree())

h2 = ax.contourf(
    ds.Longitude.data,
    ds.Latitude.data,
    speed,
    transform=ccrs.PlateCarree(),
    levels=np.arange(0, 0.6, .1),
    cmap=cmocean.cm.speed,
    extend='max'
)

cb = fig.colorbar(h2, ax=ax, orientation="vertical", shrink=.95, aspect=80)#, shrink=0.7, aspect=20*0.7)
cb.ax.tick_params(labelsize=12)
cb.set_label(f'Speed (m/s)', fontsize=12, fontweight="bold")
ax.set_title(f'Depth-Average Currents (RTOFS) - 2023-10-28T00:00:00Z', fontweight='bold')
# plt.suptitle(f'RTOFS Binary (a b file) Comparisons - {ctime.strftime("%Y-%m-%dT%H:%MZ")} - {np.floor(np.abs(z0))}m', fontweight='bold')
# plt.savefig(vname, dpi=300, bbox_inches='tight', pad_inches=0.1)
# plt.close()
plt.show()