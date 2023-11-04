import datetime as dt

import ioos_model_comparisons.configs as conf
import numpy as np
import pandas as pd
from ioos_model_comparisons.calc import lon180to360, lon360to180
from ioos_model_comparisons.models import gofs, rtofs, cmems, amseas
from ioos_model_comparisons.platforms import (get_active_gliders, 
                                  get_argo_floats_by_time,
                                  get_bathymetry)
from ioos_model_comparisons.plotting_ugos import (surface_current_fronts, 
                                      surface_current_fronts_single,
                                      surface_current_15knot_all
                                    #   plot_model_region_comparison_streamplot
                                      )
from ioos_model_comparisons.plotting import plot_model_region_comparison_streamplot
from ioos_model_comparisons.regions import region_config
import matplotlib
import time

startTime = time.time() # Start time to see how long the script took
matplotlib.use('agg')

# Set path to save plots
path_save = (conf.path_plots / "maps")

# Which models should we plot?
plot_rtofs = True
plot_gofs = True
plot_cmems = True
plot_amseas = True
plot_cnaps = False

# initialize keyword arguments for map plots
kwargs = dict()
kwargs['transform'] = conf.projection
kwargs['dpi'] = conf.dpi
kwargs['overwrite'] = False

# For debug purposes. Comment this out when commiting to repo.
conf.days = 2

# Get today and yesterday dates
target_time = dt.time(12, 0, 0)
today = dt.date.today()
date_start = dt.datetime.combine(today - dt.timedelta(days=conf.days), target_time)
date_end = dt.datetime.combine(today + dt.timedelta(days=1), target_time)

# Formatter for time
tstr = '%Y-%m-%d %H:%M:%S'

# Create dates that we want to plot
date_list = pd.date_range(date_start, date_end, freq='24H')

# This is the initial time to start the search for argo/gliders
search_start = date_list[0] - dt.timedelta(hours=conf.search_hours)

region = region_config('gom') #gom, loop_current, yucatan
extent = region['extent']
print(f'Region: {region["name"]}, Extent: {extent}')
kwargs['path_save'] = path_save / region['folder']
if conf.argo:
    argo_data = get_argo_floats_by_time(extent, search_start, date_list[-1])
else:
    argo_data = pd.DataFrame()

if conf.gliders:
    glider_data = get_active_gliders(extent, search_start, date_list[-1], parallel=False)
else:
    glider_data = pd.DataFrame()

if conf.bathy:
    bathy_data = get_bathymetry(extent)

if plot_rtofs:
    # Load RTOFS DataSet
    rds = rtofs() 

    # Save rtofs lon and lat as variables to speed up indexing calculation
    grid_lons = rds.lon.values[0,:]
    grid_lats = rds.lat.values[:,0]
    grid_x = rds.x.values
    grid_y = rds.y.values

if plot_gofs:
    # Load GOFS DataSet
    gds = gofs(rename=True)

if plot_cmems:
    # Load Copernicus
    cds = cmems(rename=True)

if plot_amseas:
    # Load AMSEAS
    am = amseas(rename=True)

# Load TOPS
# import xarray as xr
# fname = '/Users/mikesmith/Downloads/GOM22 Fronts/tops_IAS16_20220904_20220904.nc'
# tops = xr.open_dataset(fname).rename({"sst": "temperature", "sss": "salinity", 'uvel': 'u', 'vvel': 'v'})
# tops.attrs['model'] = 'TOPS'

# Loop through times
for ctime in date_list:

    # Deal with time related variables
    # ctime = pd.Timestamp(2022, 9, 4, 12, 0, 0)
    search_window_t0 = (ctime - dt.timedelta(hours=conf.search_hours)).strftime(tstr)
    search_window_t1 = ctime.strftime(tstr) 

    print(f"Checking if {ctime} exists for each model.")
    try:
        rdt = rds.sel(time=ctime, depth=0)
        print(f"RTOFS: True")
        rdt_flag = True
    except (KeyError, NameError) as error:
        print(f"RTOFS: False - {error}")
        rdt_flag = False

    try:
        gdt = gds.sel(time=ctime, depth=0)
        print(f"GOFS: True")
        gdt_flag = True
    except (KeyError, NameError) as error:
        print(f"GOFS: False - {error}")
        gdt_flag = False

    try:
        cdt = cds.sel(time=ctime, depth=0, method='nearest')
        print(f"CMEMS: True")
        cdt_flag = True
    except (KeyError, NameError) as error:
        print(f"CMEMS: False - {error}")
        cdt_flag = False

    try:
        amt = am.sel(time=ctime, depth=0)
        print(f"AMSEAS: True")
        amt_flag = True
    except (KeyError, NameError) as error:
        print(f"AMSEAS: False - {error}")
        amt_flag = False

    # try:
    #     topst = tops.sel(time=ctime)
    #     print(f"TOPS: True")
    #     tops_flag = True
    # except KeyError as error:
    #     print(f"TOPS: False - {error}")
    #     tops_flag = False
    print("\n")

    if 'eez' in region:
        kwargs["eez"] = region["eez"]

    if 'figure' in region:
        if 'legend' in region['figure']:
            kwargs['cols'] = region['figure']['legend']['columns']

        if 'figsize' in region['figure']:
            kwargs['figsize'] = region['figure']['figsize']

    try:
        kwargs['bathy'] = bathy_data.sel(
            longitude=slice(extent[0] - 1, extent[1] + 1),
            latitude=slice(extent[2] - 1, extent[3] + 1)
        )
    except NameError:
        pass
            
    extended = np.add(extent, [-1, 1, -1, 1]).tolist()
    # Find x, y indexes of the area we want to subset
    lons_ind = np.interp(extended[:2], grid_lons, grid_x)
    lats_ind = np.interp(extended[2:], grid_lats, grid_y)
    
    # convert from 360 to 180 lon
    lon360 = lon180to360(extended[:2]) 

    # Use np.floor on the 1st index and np.ceil on the 2nd index of each slice 
    # in order to widen the area of the extent slightly.
    extent_ind = [
        np.floor(lons_ind[0]).astype(int),
        np.ceil(lons_ind[1]).astype(int),
        np.floor(lats_ind[0]).astype(int),
        np.ceil(lats_ind[1]).astype(int)
        ]

    
    if rdt_flag:
        # Subset each model to the proper extent
        # Use .isel selector on x/y since we know indexes that we want to slice
        rds_sub = rdt.isel(
            x=slice(extent_ind[0], extent_ind[1]), 
            y=slice(extent_ind[2], extent_ind[3])
            ).set_coords(['u', 'v'])

    if gdt_flag:
        # subset dataset to the proper extents for each region
        gds_sub = gdt.sel(
            lon=slice(lon360[0], lon360[1]),
            lat=slice(extended[2], extended[3])
        ).set_coords(['u', 'v'])
        gds_sub['lon'] = lon360to180(gds_sub['lon']) # Convert from 0,360 lon to -180,180

    if cdt_flag:
        cds_sub = cdt.sel(
            lon=slice(extended[0], extended[1]),
            lat=slice(extended[2], extended[3])
        ).set_coords(['u', 'v'])

    if amt_flag:
        am_sub = amt.sel(
            lon=slice(lon360[0], lon360[1]),
            lat=slice(extended[2], extended[3])
        ).set_coords(['u', 'v'])

    # tops_sub = topst.sel(
    #     lon=slice(extended[0], extended[1]),
    #     lat=slice(extended[2], extended[3])
    # ).set_coords(['u', 'v'])

    # Check if any asset data exists and subset to appropriate region and time

    # Was any argo data downloaded?
    if not argo_data.empty:
        argo_lon = argo_data['lon']
        argo_lat = argo_data['lat']
        argo_region = argo_data[
            (extended[0] <= argo_lon) & (argo_lon <= extended[1]) & (extended[2] <= argo_lat) & (argo_lat <= extended[3])
        ]
        argo_region.sort_index(inplace=True)
        idx = pd.IndexSlice
        kwargs['argo'] = argo_region.loc[idx[:, search_window_t0:search_window_t1], :]

    # Was any glider data downloaded?
    if not glider_data.empty:
        glider_lon = glider_data['lon']
        glider_lat = glider_data['lat']
        glider_region = glider_data[
            (extended[0] <= glider_lon) & (glider_lon <= extended[1]) & (extended[2] <= glider_lat) & (glider_lat <= extended[3])
            ]
        glider_region = glider_region[
            (search_window_t0 <= glider_region.index.get_level_values('time'))
            &
            (glider_region.index.get_level_values('time') <= search_window_t1)
            ]
        kwargs['gliders'] = glider_region


    # try:
        # if rdt_flag and gdt_flag and cdt_flag and tops_flag:
        #     surface_current_fronts(
        #         rds_sub, 
        #         gds_sub,
        #         cds_sub,
        #         tops_sub,
        #         region,ummus
        #         **kwargs
        #         )
    try:
        surface_current_fronts_single(rds_sub, region, **kwargs)
    except Exception as e:
        print(f"Failed to process RTOFS at {ctime}")
        print(f"Error: {e}")

    try:
        surface_current_fronts_single(gds_sub, region, **kwargs)
    except Exception as e:
        print(f"Failed to process GOFS at {ctime}")
        print(f"Error: {e}")

    try:
        surface_current_fronts_single(cds_sub, region, **kwargs)
    except Exception as e:
        print(f"Failed to process CMEMS at {ctime}")
        print(f"Error: {e}")
    try:
        surface_current_fronts_single(am_sub, region, **kwargs)
    except Exception as e:
        print(f"Failed to process AMSEAS at {ctime}")
        print(f"Error: {e}")
        
        # surface_current_15knot_all(rds_sub, gds_sub, cds_sub, am_sub, region, **kwargs)

    # url = 'http://3.236.148.88:8080/thredds/dodsC/fmrc/useast_coawst_roms/COAWST-ROMS_SWAN_Forecast_Model_Run_Collection_best.ncd'
    # import xarray as xr
    # # 
    # ds = xr.open_dataset(url)
    # # ds.attrs['model'] = 'CNAPS'
    # tds = ds.sel(time=ctime).isel(eta_rho=slice(100, 300), xi_rho=slice(0, 220))
    # from oceans.ocfis import uv2spdir
    # ang, spd = uv2spdir(tds['u_eastward'], tds['v_northward'])
    # tds['speed'] = (('s_rho', 'eta_rho', 'xi_rho'), spd)
    # tds = tds.drop('ocean_time').squeeze()

    # tds = tds.rename({
    #     "lat_rho": "lat",
    #     "lon_rho": "lon", 
    #     "u_eastward": "u",
    #     "v_northward": "v"
    #     }
    #                  )
    # surface_current_15knot_all(rds_sub, gds_sub, cds_sub, am_sub, tds, region, **kwargs)

    # surface_current_fronts_single(tds, region, **kwargs)
    # # surface_current_fronts_single(tops_sub, region, **kwargs)

    # # except Exception as e:
    #     # print(f"Failed to process RTOFS vs GOFS vs CMEMS vs AMSEAS at {ctime}")
    #     # print(f"Error: {e}")
 