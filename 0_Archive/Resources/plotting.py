import os
# import pickle
import warnings
from itertools import cycle
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import matplotlib.colors
# import matplotlib.lines as mlines
import matplotlib.pyplot as plt
# import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from cartopy.io.shapereader import Reader
from cartopy.mpl.geoaxes import GeoAxesSubplot
from cool_maps.calc import categorical_cmap
from cool_maps.plot import (add_bathymetry, 
                            add_features, 
                            add_ticks,
                            create, 
                            save_fig, 
                            load_fig)
from matplotlib.colors import TwoSlopeNorm
# from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition, inset_axes
from oceans.ocfis import spdir2uv, uv2spdir
from scipy.io import loadmat
from shapely.geometry.polygon import LinearRing

import ioos_model_comparisons.configs as conf
from ioos_model_comparisons.calc import dd2dms
import scipy.ndimage as ndimage
import matplotlib.lines as mlines

# Suppresing warnings for a "pretty output."
warnings.simplefilter("ignore")

proj = dict(
    map=ccrs.Mercator(), # the projection that you want the map to be in
    data=ccrs.PlateCarree() # the projection that the data is. 
    )


def export_fig(path, fname, script=None, dpi=150):
    """
    Helper function to save a figure with some nice formatting.
    Include script to print the script that created the plot for future ref.

    Args:
        path (str or Path): Full file name including path
        script (str, optional): Print name of script on plot. Defaults to None.
        dpi (int, optional): Dots per inch. Defaults to 150.
    """
    
    if isinstance(path, str):
        path = Path(path)
    
    os.makedirs(path, exist_ok=True)
    
    if script:
        import datetime as dt
        now = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
        plt.figtext(.98, 0.20, f"{script} {now}",  fontsize=10, rotation=90)
        
    plt.savefig(path / fname, dpi=dpi, bbox_inches='tight', pad_inches=0.1)


def cmaps(variable):
    if variable == 'salinity':
        cmap = cmocean.cm.haline
    elif variable == 'temperature':
        cmap = cmocean.cm.thermal
    elif variable == 'sea_surface_height':
        cmap = cmocean.cm.balance
    return cmap


def map_add_argo(ax, df, transform=proj['data']):
    tdf = df.reset_index()
    most_recent = tdf.loc[tdf.groupby('argo')['time'].idxmax()]

    if most_recent.shape[0] > 50:
        custom_cmap = matplotlib.colors.ListedColormap('red', N=most_recent.shape[0])
        marker = cycle(['o'])
    else:
        custom_cmap = categorical_cmap(10, 5, cmap="tab10")
        marker = cycle(['o', 'h', 'p'])

    n = 0
    for float in most_recent.itertuples():
        ax.plot(float.lon, float.lat, 
                marker=next(marker), linestyle="None",
                markersize=7, markeredgecolor='black', 
                color=custom_cmap.colors[n],
                label=float.argo,
                transform=transform,
                zorder=10000)
        # map_add_legend(ax)
        n = n + 1


def map_add_all_argo(ax, df, transform=proj['data'], markersize=7):
    grouped = df.groupby(['lon', 'lat'])
    for i, x in grouped:
        ax.plot(i[0], i[1], marker='o', markersize=markersize, markeredgecolor='black', color='green', transform=transform)


def map_add_currents(ax, ds, coarsen=None, ptype="quiver",
                    scale=90, headwidth=2.75, headlength=2.75, headaxislength=2.5,
                    density=2, linewidth=.75, color='black',
                    transform=proj['data']):
    """
    Add currents to map

    Args:
        ax (ax): matplotlib.Axes
        ds (xarray.DataSet): xarray 
        coarsen (_type_, optional): Amount to downsample by. Defaults to None.
        ptype (str, optional): Plot type: "quiver" or "streamplot". Defaults to "quiver".
        scale (int, optional): _description_. Defaults to 90.
        headwidth (float, optional): _description_. Defaults to 2.75.
        headlength (float, optional): _description_. Defaults to 2.75.
        headaxislength (float, optional): _description_. Defaults to 2.5.
        transform (_type_, optional): _description_. Defaults to ccrs.PlateCarree().
        density (int, optional): _description_. Defaults to 3.
        linewidth (float, optional): Line width for streamplot. Defaults to .75.
        color (str, optional): Line color for streamplot. Defaults to 'black'.

    Returns:
        _type_: _description_
    """
    angle, speed = uv2spdir(ds['u'], ds['v'])  # convert u/v to angle and speed
    
    if ptype == "quiver":
        if coarsen:
            try:
                ds = ds.coarsen(lon=coarsen, boundary='pad').mean().coarsen(lat=coarsen, boundary='pad').mean()
                mesh = True
            except ValueError:
                ds = ds.coarsen(x=coarsen, boundary='pad').mean().coarsen(y=coarsen, boundary='pad').mean()
                mesh = False

        u, v = spdir2uv(  # convert angle and speed back to u/v, normalizing the arrow sizes
            np.ones_like(speed),
            angle,
            deg=True
        )

        qargs = {}
        qargs['scale'] = scale
        qargs['headwidth'] = headwidth
        qargs['headlength'] = headlength
        qargs['headaxislength'] = headaxislength
        qargs['transform'] = transform

        if mesh:
            lons, lats = np.meshgrid(ds['lon'], ds['lat'])
            q = ax.quiver(lons, lats, u, v, **qargs)
        else:
            q = ax.quiver(
                ds.lon.squeeze().data,
                ds.lat.squeeze().data, 
                u.squeeze(), 
                v.squeeze(), 
                **qargs)
    elif ptype == "streamplot":
        lons = ds.lon.squeeze().data
        lats = ds.lat.squeeze().data
        u = ds.u.squeeze().data
        v = ds.v.squeeze().data
        
        sargs = {}
        sargs["transform"] = transform
        sargs["density"] = density
        sargs["linewidth"] = linewidth
        if color:
            sargs["color"] = color
        else:
            sargs["color"] = speed
            sargs["cmap"] = cmocean.cm.speed
        q = ax.streamplot(lons, lats, u, v, **sargs)
    return q


def map_add_eez(ax, zorder=1, color='white', linewidth=0.75, linestyle='-.'):
    shape_feature = cfeature.ShapelyFeature(
        Reader(conf.eez_path).geometries(), 
        proj['data'],
        linestyle=linestyle,
        linewidth=linewidth,
        edgecolor=color, 
        facecolor='none'
        )
    h = ax.add_feature(shape_feature, zorder=zorder)
    return h


def map_add_gliders(ax, df, transform=proj['data'], color='white'):
    for g, new_df in df.groupby(level=0):
        q = new_df.iloc[-1]
        ax.plot(new_df['lon'], new_df['lat'], color=color,
                linewidth=1.5, transform=transform, zorder=10000)
        ax.plot(q['lon'], q['lat'], marker='^', markeredgecolor='black',
                markersize=8.5, label=g, transform=transform, zorder=10000)
        # map_add_legend(ax)


def map_add_inset(ax, x=.8, y=.3, size=.5, extent=None, zoom_extent=None):
    """_summary_

    Args:
        ax (_type_): _description_
        x (float, optional): inset x location relative to main plot (ax) in normalized units. Defaults to .8.
        y (float, optional): inset y location relative to main plot (ax) in normalized units. Defaults to .3.
        size (float, optional): _description_. Defaults to 0.5.

    Returns:
        _type_: _description_
    """
    import cartopy
    # Inset Axis
    # axin = plt.axes([0, 0, 1, 1], projection=ccrs.Mercator())
    # position = [x - size / 2, y - size / 2, size, size]
    # ip = InsetPosition(ax, position)
    # axin.set_axes_locator(ip)
    axins = inset_axes(ax, width="40%", height="40%", loc="lower left", 
                       axes_class=cartopy.mpl.geoaxes.GeoAxes, 
                       axes_kwargs=dict(map_projection=ccrs.Mercator()))
    axins.set_extent(extent)

    if zoom_extent:
        lonmin, lonmax, latmin, latmax = zoom_extent

        nvert = 100
        lons = np.r_[np.linspace(lonmin, lonmin, nvert),
                    np.linspace(lonmin, lonmax, nvert),
                    np.linspace(lonmax, lonmax, nvert)].tolist()
        lats = np.r_[np.linspace(latmin, latmax, nvert),
                    np.linspace(latmax, latmax, nvert),
                    np.linspace(latmax, latmin, nvert)].tolist()
        
        ring = LinearRing(list(zip(lons, lats)))
        axins.add_geometries([ring], ccrs.PlateCarree(), facecolor='none', edgecolor='red', linewidth=1, zorder=1000)
    return axins


def map_add_legend(ax):
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


def map_add_transects(ax, transects, transform=proj['data']):
    ax.plot(transects['lon'], transects['lat'], 'r-', transform=transform)


def add_colorbar(ax, h, location="bottom", constrain=True):

    if isinstance(ax, list):
        # Multiple axes are being passed
        multi = True
    elif isinstance(ax, GeoAxesSubplot):
        multi = False

    # We want to minimize input parameters to functions. We need the fig
    # so we can use the .get_figure() method on the axes 
    if multi:
        fig = ax[0].get_figure()
    else:
        fig = ax.get_figure()

    if constrain:
        # Constrain the colorbar to the size of the plot
        if location == "bottom" or location == "top":
            if multi:
                # length = len(ax)
                widths = 250 + sum([a.bbox.width for a in ax])
                ax_ratio = 0.047*(widths/ax[0].bbox.height)
            else:
                ax_ratio = 0.047*(ax.bbox.width/ax.bbox.height)
        elif location == "left" or location == "right":
            if multi:
                ax_ratio = 0.047*(sum([a.bbox.height for a in ax])/ax[0].bbox.width)
            else:
                ax_ratio = 0.047*(ax.bbox.height/ax.bbox.width)
    else:
        ax_ratio = 0.15 # This is the default in matplotlib.
        
    # Add colorbar to axes
    cb = fig.colorbar(h, ax=ax, location=location, fraction=ax_ratio)
    return cb


def plot_model_region(ds, region,
                      bathy=None,
                      argo=None,
                      gliders=None,
                      currents=dict(bool=False),
                      transform=dict(
                          map=proj['map'],
                          data=proj['data']
                          ),
                      legend=True,
                      model='rtofs',
                      path_save=os.getcwd(),
                      dpi=150,
                      t0=None):
    """

    :param lon: longitude
    :param lat: latitude
    :param variable: data variable you want to plot
    :param kwargs:
    :return:
    """
    region_name = region["name"]
    extent = region["extent"]
    time = pd.to_datetime(ds.time.values)

    # Create subdirectory for region
    region_file_str = '_'.join(region_name.lower().split(' '))
    path_save_region = path_save / 'regions' / region_file_str

    if not isinstance(gliders, pd.DataFrame):
        gliders = pd.DataFrame()
    
    if not isinstance(argo, pd.DataFrame):
        argo = pd.DataFrame()

    # Iterate through the region dictionary. This dict contains information
    # on what variables and depths to plot. 
    for key, values in region["variables"].items():
        # Create subdirectory for variable under region directory
        var_str = ' '.join(key.split('_')).title()

        # Iterate through values of the key
        for item in values:
            depth = item['depth']

            # Select variable and depth to plot
            # print(ds[k].name)
            try:
                da = ds[key].sel(depth=depth)
            except KeyError:
                da = ds[key]

            # Create subdirectory for depth under variable subdirectory
            save_dir_final = path_save_region / f"{key}_{depth}m" / time.strftime('%Y/%m')
            os.makedirs(save_dir_final, exist_ok=True)

            # Create a string for the title of the plot
            title_time = time.strftime("%Y-%m-%d %H:%M:%S")
            title = f"{model.upper()} - {var_str} ({depth} m) - {title_time}\n"

            # if not gliders.empty or not argo.empty:
            #         title += f'Assets ({str(t0)} to {str(time)})'

            # Create a file name to save the plot as
            sname = f'{model}-{key}-{time.strftime("%Y-%m-%dT%H%M%SZ")}'
            save_file = save_dir_final / f"{sname}.png"

            # Create a map figure and serialize it if one doesn't already exist
            region_name = "_".join(region["name"].split(' ')).lower()
            path_maps = path_save / "mapfigs"
            os.makedirs(path_maps, exist_ok=True)
            sfig = (path_maps / f"{region_name}_fig.pkl")

            # if not sfig.exists():
                # Create an empty projection within set extent
            fig, ax = create(extent, proj=transform['map'])

            # Add bathymetry
            if bathy:
                add_bathymetry(ax,
                                bathy.longitude.values, 
                                bathy.latitude.values, 
                                bathy.elevation.values,
                                levels=(-1000, -100),
                                zorder=1.5)

            #     save_fig(fig, path_maps, f"{region_name}_fig.pkl")       
            # else:
            #     fig = load_fig(sfig)
            #     ax = fig.axes[0]
                               
            rargs = {}
            rargs['argo'] = argo
            rargs['gliders'] = gliders
            rargs['transform'] = transform['data']  
            plot_regional_assets(ax, **rargs)

            cargs = {}
            cargs['vmin'] = item['limits'][0]
            cargs['vmax'] = item['limits'][1]
            cargs['transform'] = transform['data']
            cargs['cmap'] = cmaps(da.name)
            cargs['levels'] = np.arange(cargs['vmin'], cargs['vmax'], item['limits'][2])
            cargs['extend'] = 'both'

            try:
                cargs.pop('vmin'), cargs.pop('vmax')
            except KeyError:
                pass
            
            # If the xarray DataArray contains data, let's contour the data.
            if da is not None:
                h = ax.contourf(da['lon'], da['lat'], da.squeeze(), **cargs)

                # Create the colorbar
                axins = inset_axes(ax,  # here using axis of the lowest plot
                    width="2.5%",  # width = 5% of parent_bbox width
                    height="100%",  # height : 340% good for a (4x4) Grid
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0., 1, 1),
                    bbox_transform=ax.transAxes,
                    borderpad=0
                    )
                cb = plt.colorbar(h, cax=axins)
                cb.ax.tick_params(labelsize=12)
                cb.set_label(f'{da.name.title()} ({da.units})', fontsize=13)

            ax.set_title(title, fontsize=16, fontweight='bold')

            if legend:
                h, l = ax.get_legend_handles_labels()  # get labels and handles from ax1

                if (len(h) > 0) & (len(l) > 0):
                    # Shrink current axis's height by 10% on the bottom
                    box = ax.get_position()
                    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                                    box.width, box.height * 0.9])

                    # Put a legend below current axis
                    ax.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, -0.05),
                            fancybox=True, shadow=True, ncol=5)
                    legstr = f'Glider/Argo Search Window: {str(t0)} to {str(time)}'
                    plt.figtext(0.5, -0.07, legstr, ha="center", fontsize=10, fontweight='bold')

            fig.savefig(save_file, dpi=dpi, bbox_inches='tight', pad_inches=0.1)
            
            # Add currents
            if currents['bool']:
                quiver_dir = save_dir_final / "currents"
                os.makedirs(quiver_dir, exist_ok=True)
                
                save_file_q = quiver_dir / f"{sname}.png"
                coarsen = currents['coarsen']
                map_add_currents(ax, da, coarsen=coarsen[model], **currents['kwargs'])
                fig.savefig(save_file_q, dpi=dpi, bbox_inches='tight', pad_inches=0.1)
            
            plt.close()

def plot_model_region_currents(region,
                      transform=dict(
                          map=proj['map'],
                          data=proj['data']
                          ),
                      legend=True,
                      path_save=Path(os.getcwd()),
                      dpi=150,
                      t0=None,
                      ):
    """

    :param lon: longitude
    :param lat: latitude
    :param variable: data variable you want to plot
    :param kwargs:
    :return:
    """
    region_name = region["name"]
    extent = region["extent"]

    # Create subdirectory for region
    region_file_str = '_'.join(region_name.lower().split(' '))
    path_save_region = path_save / 'regions' / region_file_str
    
    # Create subdirectory for depth under variable subdirectory
    save_dir_final = path_save_region 
    os.makedirs(save_dir_final, exist_ok=True)

    # Create a string for the title of the plot
 
    # Create a file name to save the plot as
    save_file = save_dir_final / f"passengers.png"

    fig, ax = create(extent, proj=transform['map'], bathymetry=False, figsize=(16,9))
                               
    url = 'https://encdirect.noaa.gov/arcgis/services/encdirect/enc_overview/MapServer/WMSServer?request=GetCapabilities&service=WMS'
    url2 = 'https://encdirect.noaa.gov/arcgis/services/encdirect/enc_general/MapServer/WMSServer?request=GetCapabilities&service=WMS'
    # url = 'https://encdirect.noaa.gov/arcgis/services/encdirect/enc_coastal/MapServer/WMSServer?request=GetCapabilities&service=WMS'

    layers2 = [
        '1',
        '2',
        '3',
        '4',
        '5',
        #  '6',
        '7',
        '8',
        # '9',
        '10',
        '11',
        '12',
        '13',
        # '14',
        '15',
        '16',
        '17',
        '18',
        '19',
        # '20',
        '21',
        # '22',
        '23',
        '24',
        '25',
        '26',
        '27',
        '28',
        '29',
        '30',
        '31',
        # '32',
        '33',
        '34',
        # '35',
        '36',
        '37',
        # '38',
        '39',
        # '40',
        '41',
        '42',
        '43',
        # '44',
        '45',
        '46',
        '47',
        '48',
        # '49',
        '50',
        '51',
        '52',
        # '53',
        '54',
        # '55',
        '56',
        # '57',
        '58',
        '59',
        # '60',
        '61',
        # '62',
        # '63',
        '64',
        # '65',
        '66',
        '67',
        # '68',
        '69',
        # '70',
        '71',
        '72',
        '73',
        # '74',
        '75',
        # '76',
        '77',
        '78',
        # '79',
        '80',
        '81',
        # '82',
        '83',
        '84',
        '85',
        '86',
        # '87',
        '88',
        # '89',
        '90',
        '91',
        '92',
        '93',
        '94',
        # '95',
        '96',
        '97',
        '98',
        '99',
        '100',
        # '101',
        '102',
        '103',
        '104',
        '105',
        '106',
        '107',
        '108',
        # '109',
        '110',
        # '111',
        '112',
        '113',
        '114',
        '115',
        '116',
        '117',
        '118',
        '119',
        '120',
        '121',
        '122',
        '123',
        '124']

    layers = [
        '1',
        '2',
        '3',
        '4',
        '5',
        # '6',
        '7',
        '8',
        # '9',
        '10',
        '11',
        '12',
        '13',
        '14',
        '15',
        # '16',
        '17',
        # '18',
        '19',
        # '20',
        '21',
        # '22',
        '23',
        '24',
        # '25',
        '26',
        # '27',
        '28',
        '29',
        # '30',
        '31',
        # '32',
        '33',
        # '34',
        '35',
        '36',
        # '37',
        '38',
        '39',
        # '40',
        '41',
        '42',
        # '43',
        '44',
        # '45',
        '46',
        # '47',
        '48',
        '49',
        # '50',
        '51',
        '52',
        # '53',
        '54',
        # '55',
        '56',
        # '57',
        '58',
        '59',
        # '60',
        '61',
        '62',
        '63',
        # '64',
        '65',
        # '66',
        '67',
        '68',
        '69',
        '70',
        '71',
        # '72',
        '73',
        '74',
        '75',
        # '76',
        '77',
        # '78',
        '79',
        '80',
        '81',
        '82',
        '83',
        '84',
        # '85',
        '86',
        '87',
        '88',
        '89',
        '90',
        '91',
        '92',
        '93',
        '94',
        '95',
        '96',]


    
    # ax.add_wms(wms=url2, layers=layers2)
    
    # map_add_eez(ax, color='red')

        # # Create the colorbar
        # axins = inset_axes(ax,  # here using axis of the lowest plot
        #     width="2.5%",  # width = 5% of parent_bbox width
        #     height="100%",  # height : 340% good for a (4x4) Grid
        #     loc='lower left',
        #     bbox_to_anchor=(1.05, 0., 1, 1),
        #     bbox_transform=ax.transAxes,
        #     borderpad=0
        #     )
        # cb = plt.colorbar(h, cax=axins)
        # cb.ax.tick_params(labelsize=12)
        # cb.set_label(f'{da.name.title()} ({da.units})', fontsize=13)

    # ax.set_title(title, fontsize=16, fontweight='bold')

    if legend:
        h, l = ax.get_legend_handles_labels()  # get labels and handles from ax1

        if (len(h) > 0) & (len(l) > 0):
            # Shrink current axis's height by 10% on the bottom
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + box.height * 0.1,
                            box.width, box.height * 0.9])

            # Put a legend below current axis
            ax.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, -0.05),
                    fancybox=True, shadow=True, ncol=5)
            legstr = f'Glider/Argo Search Window: {str(t0)} to {str(time)}'
            plt.figtext(0.5, -0.07, legstr, ha="center", fontsize=10, fontweight='bold')

    # fig.savefig(save_file, dpi=dpi, bbox_inches='tight', pad_inches=0.1)
    
    # Add currents
    # coarsen = currents['coarsen']
    # map_add_currents(ax, tds, coarsen=coarsen[model], **currents['kwargs'])

    # field_lon = [-86.6021, -85.5144,-84.9706, -84.9102, -85.1849, -85.5144, -86.6076]
    # field_lat = [21.2004, 22.5968, 21.8544, 21.0056, 20.1880, 19.7544, 21.1901]
    # ax.fill(field_lon, field_lat, color='coral', edgecolor='black', transform=ccrs.PlateCarree(), alpha=0.6, zorder=3000)

    style = "Simple, tail_width=0.5, head_width=8, head_length=8"

    import matplotlib.patches as patches

    # Triangle
    kw = dict(arrowstyle=style, color="red", linewidth=4, transform=ccrs.PlateCarree(), zorder=3000, label='Fixed')
    a1 = patches.FancyArrowPatch((-87, 19.25), (-84.5, 21.5), **kw)
    a2 = patches.FancyArrowPatch((-84.5, 21.5), (-85.25, 19.25), **kw)
    a3 = patches.FancyArrowPatch((-85.25, 19.25), (-87, 19.25), **kw)

    plt.gca().add_patch(a1)
    plt.gca().add_patch(a2)
    plt.gca().add_patch(a3)

    # # Triangle 2
    # kw = dict(arrowstyle=style, color="green", linewidth=4, transform=ccrs.PlateCarree(), zorder=3000)
    # a4 = patches.FancyArrowPatch((-86.25, 19.5), (-82.75, 21), **kw)
    # a5 = patches.FancyArrowPatch((-82.6, 21), (-83.5, 18.5), **kw)
    # a6 = patches.FancyArrowPatch((-83.5, 18.25), (-86.25, 19), **kw)
    
    # plt.gca().add_patch(a4)
    # plt.gca().add_patch(a5)
    # plt.gca().add_patch(a6)

    # Curvy area
    kw = dict(arrowstyle=style, color="magenta", linewidth=4, transform=ccrs.PlateCarree(), zorder=3000, label='Expected')
    a4 = patches.FancyArrowPatch((-86.5, 20.5), (-86, 21.25), **kw)
    # a2 = patches.FancyArrowPatch((-85.75, 21.5), (-88, 25), connectionstyle="angle3, angleA=90, angleB=20", **kw)
    a5 = patches.FancyArrowPatch((-85.75, 21.5), (-88, 25), connectionstyle="angle3, angleA=70, angleB=15", **kw)
    a6 = patches.FancyArrowPatch((-85.75, 21.5), (-84, 23.5), connectionstyle="angle3, angleA=75, angleB=-20", **kw)
    a7 = patches.FancyArrowPatch((-85.75, 21.5), (-84, 24.5), connectionstyle="angle3, angleA=75, angleB=-40", **kw)
    a8 = patches.FancyArrowPatch((-84, 24), (-82.25, 24), **kw)
    # a6 = patches.FancyArrowPatch((-82, 24), (-81.25, 23.5), **kw)
    # a7 = patches.FancyArrowPatch((-80.75, 24), (-80, 25), **kw)
    # a3 = patches.FancyArrowPatch((-0.4, -0.6), (0.4, -0.6), connectionstyle="arc3,rad=.5", **kw)

    # for a in [a1, a2, a3]:
        # plt.gca().add_patch(a)
    plt.gca().add_patch(a4)
    plt.gca().add_patch(a5)
    plt.gca().add_patch(a6)
    plt.gca().add_patch(a7)
    plt.gca().add_patch(a8)
    # plt.gca().add_patch(a6)
    # plt.gca().add_patch(a7)

    # Area 1
    field_lon = [-86.5, -86.5, -87.3, -87.3, -86.5]
    field_lat = [19, 20, 20, 19, 19]
    # ax.fill(field_lon, field_lat, color='coral', edgecolor='black', transform=ccrs.PlateCarree(), alpha=0.6, zorder=3000)
    l1 = ax.plot(field_lon, field_lat, color='maroon', linewidth=4, zorder=3000, transform=ccrs.PlateCarree(), label='Fixed')
    
    # Area 2
    field_lon = [-86.75, -86, -85.5, -86.25, -86.75]
    field_lat = [20.75, 20.75, 21.5, 21.5, 20.75]
    # # ax.fill(field_lon, field_lat, color='coral', edgecolor='black', transform=ccrs.PlateCarree(), alpha=0.6, zorder=3000)
    ax.plot(field_lon, field_lat, color='maroon', linewidth=4, zorder=3000, transform=ccrs.PlateCarree())

    # # Florida
    ax.text(-82, 27.5, 'Florida', fontsize=14, transform=ccrs.PlateCarree(), zorder=3000)
    
    # # Cuba
    ax.text(-81, 22.4, 'Cuba', fontsize=14, transform=ccrs.PlateCarree(), zorder=3000)
    
    # # Mexico
    ax.text(-89.5, 19, 'Mexico', fontsize=14, transform=ccrs.PlateCarree(), zorder=3000)
    
    # # Yucatan/Quintana Roo
    ax.text(-89.9, 21, 'Yucatan/Quintana Roo', fontsize=11, transform=ccrs.PlateCarree(), zorder=3000)


    # # Add loop current contour from WHO Group
    # fname = '/Users/mikesmith/Downloads/GOM front/2023-01-31_fronts.mat'
    # data = loadmat(fname)
 
    # fronts = []
    # for item in data['BZ_all'][0]:
    #     loop_y = item['y'].T
    #     loop_x = item['x'].T

    #     hf = ax.plot(loop_x, loop_y,
    #                  linestyle=item['LineStyle'][0],
    #                  color='black',
    #                  linewidth=1, 
    #                  transform=ccrs.PlateCarree(), 
    #                  zorder=120
    #                  )
    #     fronts.append(hf)

    #     # Add arrows
    #     start_lon = item['bx'].T
    #     start_lat = item['by'].T
    #     end_lon = item['tx'].T
    #     end_lat = item['ty'].T

    #     for count, _ in enumerate(start_lon):
    #         ax.arrow(
    #             start_lon[count][0],
    #             start_lat[count][0],
    #             end_lon[count][0]-start_lon[count][0],
    #             end_lat[count][0]-start_lat[count][0],
    #             linewidth=0, 
    #             head_width=0.1,
    #             shape='full', 
    #             fc='black', 
    #             ec='black',
    #             transform=ccrs.PlateCarree(),
    #             zorder=130,
    #             )
    # fronts.reverse()
    url = 'https://gis.charttools.noaa.gov/arcgis/rest/services/MarineChart_Services/NOAACharts/MapServer/WMTS'
    w = ax.add_wmts(wmts=url, layer_name='MarineChart_Services_NOAACharts')
    # ax.legend([a1, l1[0], a6], ['Fixed', 'Fixed', 'Expected'], loc='upper left', title='Glider Tracks')
    legend_h = []
    import matplotlib.lines as mlines
    legend_h.append(mlines.Line2D([], [], linestyle='-', color='red', linewidth=6))
    legend_h.append(mlines.Line2D([], [], linestyle='-', color='maroon', linewidth=6))
    legend_h.append(mlines.Line2D([], [], linestyle='-', color='magenta', linewidth=6))

    ax.legend(legend_h, ['Fixed', 'Fixed', 'Expected'], loc='upper left', title='Glider Tracks', title_fontproperties={'weight':'bold'})
    plt.figtext(0.25, 0.01, 'Fixed tracks will be followed as shown.\nExpected tracks will approximately follow tracks as shown but may be affected by currents in unpredictable way.', ha="left", fontsize=9, fontstyle='italic')

    fig.savefig(save_file, dpi=150, bbox_inches='tight', pad_inches=0.1)

    plt.close()
 
def remove_quiver_handles(ax):
    for art in ax.get_children():
        if isinstance(art, matplotlib.patches.FancyArrowPatch):
            art.remove()      


def plot_model_region_comparison(ds1, ds2, region,
                                       bathy=None,
                                       argo=None,
                                       gliders=None,
                                       currents=None,
                                       eez=False,
                                       cols=6,
                                       transform=dict(map=proj['map'], 
                                                      data=proj['data']
                                                      ),
                                       path_save=os.getcwd(),
                                       figsize=(14,8),
                                       dpi=150,
                                       colorbar=True,
                                       overwrite=False
                                       ):
    
    # Convert ds.time value to a normal datetime
    time = pd.to_datetime(ds1.time.data)
    extent = region['extent']

    # Formatter for time
    tstr_title = time.strftime('%Y-%m-%d %H:%M:%S')
    tstr_folder = time.strftime('%Y-%m-%dT%H%M%SZ')
    year = time.strftime("%Y")
    month = time.strftime("%m")

    # Create subdirectory for region
    # region_file_str = ('_').join(region_name.lower().split(' '))
    # path_save_region = path_save / region['folder']
    
    # # Create a map figure and serialize it if one doesn't already exist
    # region_name = "_".join(region["name"].split(' ')).lower()
    # mdir = path_save / "mapfigs"
    # os.makedirs(mdir, exist_ok=True)
    # sfig = mdir / f"{region_name}_fig.pkl"

    # if not sfig.exists():
    # Create figure

    grid = """
    RG
    LL
    """

    fig, _ = plt.subplot_mosaic(
        grid,
        figsize=figsize,
        layout="constrained",
        subplot_kw={
            'projection': transform['map']
            },
        gridspec_kw={
            # set the height ratios between the rows
            "height_ratios": [4, 1],
            # set the width ratios between the columns
            # # "width_ratios": [1],
            },
        )
    axs = fig.axes
    ax1 = axs[0] # Model 1
    ax2 = axs[1] # Model 2
    ax3 = axs[2] # Legend for argo/gliders

    # Set map extent
    ax1.set_extent(extent)
    ax2.set_extent(extent)
          
    # Make the map pretty
    add_features(ax1)# zorder=0)
    add_features(ax2)# zorder=0)

    # Add bathymetry lines
    if bathy:
        try:
            add_bathymetry(ax1,
                            bathy.longitude.values, 
                            bathy.latitude.values, 
                            bathy.elevation.values,
                            levels=(-1000, -100),
                            zorder=1.5)
            add_bathymetry(ax2,
                            bathy.longitude.values, 
                            bathy.latitude.values, 
                            bathy.elevation.values,
                            levels=(-1000, -100),
                            zorder=1.5)
        except ValueError:
            print("Bathymetry deeper than specified levels.")

    # Add ticks
    add_ticks(ax1, extent, label_left=True)
    add_ticks(ax2, extent, label_left=False, label_right=True)


    #             with open(sfig, 'wb') as file:
    #                 pickle.dump(fig, file)
    # else:
    #     with open(sfig, "rb") as file:
    #         fig = pickle.load(file)
    #         axs = fig.axes
    #         ax1 = axs[0]
    #         ax2 = axs[1]
    #         ax3 = axs[2]

    # Plot gliders and argo floats
    rargs = {}
    rargs['argo'] = argo
    rargs['gliders'] = gliders
    rargs['transform'] = transform['data']  
    plot_regional_assets(ax1, **rargs)
    plot_regional_assets(ax2, **rargs)

    # Label the subplots
    ax1.set_title(ds1.model, fontsize=16, fontweight="bold")
    ax2.set_title(ds2.model, fontsize=16, fontweight="bold")
    txt = plt.suptitle("", fontsize=22, fontweight="bold")
    
    # Deal with the third axes
    h, l = ax1.get_legend_handles_labels()  # get labels and handles from ax1
    if (len(h) > 0) & (len(l) > 0):
        
        # Add handles to legend
        legend = ax3.legend(h, l, ncol=cols, loc='center', fontsize=8)

        # Add title to legend
        t0 = []
        if isinstance(argo, pd.DataFrame):
            if not argo.empty:
                t0.append(argo.index.min()[1])

        if isinstance(gliders, pd.DataFrame):
            if not gliders.empty:
                t0.append(gliders.index.min()[1])

        if len(t0) > 0:
            t0 = min(t0).strftime('%Y-%m-%d %H:00:00')
        else:
            t0 = None
        legstr = f'Glider/Argo Search Window: {t0} to {str(time)}'
        ax3.set_title(legstr, loc="center", fontsize=9, fontweight="bold", style='italic')
        legend._legend_box.sep = 1
        # plt.figtext(0.5, 0.001, legstr, ha="center", fontsize=10, fontweight='bold')
    ax3.set_axis_off()

    # Iterate through the variables to be plotted for each region. 
    # This dict contains information on what variables and depths to plot. 
    for k, v in region["variables"].items():
        # Create subdirectory for variable under region directory
        var_str = ' '.join(k.split('_')).title()

        for item in v:
            print(f"Plotting {k} @ {item['depth']}")
            depth = item['depth']
            rsub = ds1[k].sel(depth=depth)
            gsub = ds2[k].sel(depth=depth, method='nearest')
            
            # Create subdirectory for depth under variable subdirectory
            save_dir_final = path_save / f"{k}_{depth}m" / time.strftime('%Y/%m')
            os.makedirs(save_dir_final, exist_ok=True)

            # Create a file name to save the plot as
            # sname = f'{ds1.model}_vs_{ds2.model}_{k}-{time.strftime("%Y-%m-%dT%H%M%SZ")}'
            sname = f'{"-".join(region["folder"].split("_"))}_{time.strftime("%Y-%m-%dT%H%M%SZ")}_{k}-{depth}m_{ds1.model.lower()}-vs-{ds2.model.lower()}'
            save_file = save_dir_final / f"{sname}.png"

            if save_file.is_file():
                if not overwrite:
                    print(f"{save_file} exists. Overwrite: False. Skipping.")
                    continue
                else:
                    print(f"{save_file} exists. Overwrite: True. Replotting.")
                    
            # Add the super title (title for both subplots)
            txt.set_text(f"{var_str} ({depth} m) - {tstr_title}\n")

            # Create dictionary for variable argument inputs for contourf
            vargs = {}
            vargs['transform'] = transform['data']
            vargs['transform_first'] = True
            vargs['cmap'] = cmaps(ds1[k].name)
            vargs['extend'] = "both"

            if 'limits' in item:
                vargs['vmin'] = item['limits'][0]
                vargs['vmax'] = item['limits'][1]
                vargs['levels'] = np.arange(vargs['vmin'], vargs['vmax'], item['limits'][2])
        
            # Filled contour for each model variable
            if (rsub['lon'].ndim == 1) & rsub['lat'].ndim == 1:
                rlons, rlats = np.meshgrid(rsub['lon'], rsub['lat'])
            else:
                rlons = rsub['lon']
                rlats = rsub['lat']
            h1 = ax1.contourf(rlons, rlats, rsub.squeeze(), **vargs)

            # Check if ndims are 1, transform_first requires 2d array
            if (gsub['lon'].ndim == 1) & gsub['lat'].ndim == 1:
                glons, glats = np.meshgrid(gsub['lon'], gsub['lat'])
            else:
                glons = gsub['lon']
                glats = gsub['lat']
            h2 = ax2.contourf(glons, glats, gsub.squeeze(), **vargs)

            if colorbar:
                cb = fig.colorbar(h1, ax=axs[:2], orientation="horizontal", shrink=.95, aspect=80)#, shrink=0.7, aspect=20*0.7)
                cb.ax.tick_params(labelsize=12)
                cb.set_label(f'{k.title()} ({rsub.units})', fontsize=12, fontweight="bold")

            # Add EEZ
            if eez:
                eez1 = map_add_eez(ax1, zorder=1)
                eez2 = map_add_eez(ax2, zorder=1)

            # Save the figure. Using fig to savefig allows us to delete any
            # figure handles so that we can reuse the figure.
            # export_fig(save_dir_final, sname, dpi=dpi)
            fig.savefig(save_dir_final / sname, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

            # Add currents as overlays over variable plots
            if currents['bool']:
                quiver_dir = save_dir_final / "currents_overlay"
                os.makedirs(quiver_dir, exist_ok=True)
                
                coarsen = currents['coarsen']
                q1 = map_add_currents(ax1, rsub, coarsen=coarsen['rtofs'], **currents['kwargs'])
                q2 = map_add_currents(ax2, gsub, coarsen=coarsen["gofs"], **currents['kwargs'])

                if eez:
                    eez1._kwargs['edgecolor']= 'white'                
                    eez2._kwargs['edgecolor']= 'white'
                
                fig.savefig(quiver_dir / sname, dpi=dpi, bbox_inches='tight', pad_inches=0.1)
                # export_fig(quiver_dir, f"{sname}.png", dpi=dpi)

                # Remove quiver handles from each axes
                q1.lines.remove(), q2.lines.remove()
                remove_quiver_handles(ax1), remove_quiver_handles(ax2)

            # Remove handles so we can reuse figure
            # Delete contour handles 
            [x.remove() for x in h1.collections] # axes 1
            [x.remove() for x in h2.collections] # axes 2

            if colorbar: 
                cb.remove()
                
            if eez:
                eez1.remove()
                eez2.remove()            

            plt.close()


def plot_regional_assets(ax, argo=None, gliders=None, 
                         transects=None,
                         transform=None):
    if argo is None:
        argo = pd.DataFrame()

    if gliders is None:
        gliders = pd.DataFrame()

    if transects is None:
        transects = pd.DataFrame()

    if not argo.empty:
        map_add_argo(ax, argo, transform)

    if not gliders.empty:
        map_add_gliders(ax, gliders, transform)

    if not transects.empty:
        map_add_transects(ax, transects, transform)


def transect(fig, ax, x, y, z, c, cmap=None, levels=None, isobath=None, flip_y=None):
    cmap = cmap or 'parula'
    levels = levels or dict(deep=np.arange(0, 28), shallow=np.arange(14, 28))
    flip_y = flip_y or True
    isobath = isobath or None
    levels = levels or [26]

    if not isinstance(isobath, list):
        isobath = [isobath]

    offset = TwoSlopeNorm(vcenter=0)

    ax1 = ax.contourf(x, y, c, cmap=cmap, levels=levels['deep'], extend='both', norm=offset)

    if isobath:
        for line in isobath:
            ax.contour(x, y, c, [line], colors='k')  # add contour at 26m

    if flip_y:
        ax.set_ylim(z, 0)

    cb = fig.colorbar(ax1, ax=ax, orientation='vertical')
    cb.set_label('m/s', fontweight='bold')
    return ax


def plot_transect(x, y, c, xlabel, cmap, title=None, save_file=None, flip_y=None, levels=None, isobath=None):
    title = title or 'Transect Plot'
    save_file = save_file or 'transect.png'

    # Initiate transect plot
    fig, ax = plt.subplots(figsize=(12, 6))

    ax = transect(fig, ax, x, y, 1000, c, cmap, levels, isobath, flip_y)

    # Add titles and labels
    plt.suptitle(title, size=16, fontweight='bold')
    ax.set_ylabel('Depth (m)', size=12, fontweight='bold')
    ax.set_xlabel(xlabel, size=12, fontweight='bold')

    plt.savefig(save_file, bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.close()


def plot_transects(x, y, c, xlabel, cmap=None, title=None, save_file=None, flip_y=None, levels=None, isobath=None):
    title = title or 'Transect Plot'
    save_file = save_file or 'transect.png'

    fig, (ax1, ax2) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(11, 8)
        )

    # 1000m subplot
    ax1 = transect(fig, ax1, x, y, 1000, c, cmap, levels, isobath, flip_y)

    # 300m subplot
    ax2 = transect(fig, ax2, x, y, 100, c, cmap, levels, isobath, flip_y)

    # Add titles and labels
    plt.suptitle(title, size=16, fontweight='bold')
    ax1.set_ylabel('Depth (m)', size=12, fontweight='bold')
    ax2.set_ylabel('Depth (m)', size=12, fontweight='bold')
    ax2.set_xlabel(xlabel, size=12, fontweight='bold')

    plt.tight_layout()

    plt.savefig(save_file, bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.close()


def plot_model_region_comparison_streamplot(ds1, ds2, region,
                                                bathy=None,
                                                argo=None,
                                                gliders=None,
                                                currents=None,
                                                eez=False,
                                                cols=6,
                                                transform=dict(map=proj['map'], 
                                                                data=proj['data']
                                                                ),
                                                path_save=os.getcwd(),
                                                figsize=(14,8),
                                                dpi=150,
                                                colorbar=True,
                                                overwrite=False
                                                ):

    time = pd.to_datetime(ds1.time.data)
    extent = region['extent']
    cdict = region['currents']
    
    grid = """
    RG
    LL
    """
    # Add loop current contour from WHO Group
    # fname = '/Users/mikesmith/Downloads/GOM22 Fronts/2022-09-04_fronts.mat'
    # data = loadmat(fname)

    # Iterate through the variables to be plotted for each region. 
    # This dict contains information on what variables and depths to plot. 
    for depth in cdict['depths']:
        print(f"Plotting currents @ {depth}m")
        ds1_depth = ds1.sel(depth=depth)

        try:
            ds2_depth = ds2.sel(depth=depth, method='nearest')
        except:
            ds2_depth = ds2
        
        # Plot currents with magnitude and direction
        quiver_dir = path_save / f"currents_{depth}m" / time.strftime('%Y/%m')
        os.makedirs(quiver_dir, exist_ok=True)

        # Generate descriptive filename
        sname = f'{region["folder"]}_{time.strftime("%Y-%m-%dT%H%M%SZ")}_currents-{depth}m_{ds1.model.lower()}-vs-{ds2.model.lower()}'
        save_file_q = quiver_dir / f"{sname}.png"

        # Check if filename already exists
        if save_file_q.is_file():
            if not overwrite:
                print(f"{sname} exists. Overwrite: False. Skipping.")
                continue
            else:
                print(f"{sname} exists. Overwrite: True. Replotting.")

        # Convert u and v radial velocities to magnitude
        _, mag_r = uv2spdir(ds1_depth['u'], ds1_depth['v'])
        _, mag_g = uv2spdir(ds2_depth['u'], ds2_depth['v'])

        # Initialize qargs dictionary for input into contour plot of magnitude
        qargs = {}
        qargs['transform'] = transform['data']
        qargs['cmap'] = cmocean.cm.speed
        qargs['extend'] = "max"

        if 'limits' in cdict:
            lims = cdict['limits']
            qargs['levels'] = np.arange(lims[0], lims[1]+lims[2], lims[2])

        # Initialize figure
        fig, _ = plt.subplot_mosaic(
            grid,
            figsize=figsize,
            layout="constrained",
            subplot_kw={
                'projection': proj['map']
                },
            gridspec_kw={
                # set the height ratios between the rows
                "height_ratios": [4, 1],
                # set the width ratios between the columns
                # # "width_ratios": [1],
                },
            )
        axs = fig.axes
        ax1 = axs[0] # rtofs
        ax2 = axs[1] # gofs
        ax3 = axs[2] # legend for argo/gliders

        # Set map extents
        ax1.set_extent(extent)
        ax2.set_extent(extent)

        # Make the map pretty  
        add_features(ax1)# zorder=0)
        add_features(ax2)# zorder=0)
        if bathy:       
            try:
                add_bathymetry(ax1,
                            bathy.longitude.values, 
                            bathy.latitude.values, 
                            bathy.elevation.values,
                            levels=(-1000, -100),
                            zorder=1.5
                            )

                add_bathymetry(ax2,
                            bathy.longitude.values, 
                            bathy.latitude.values, 
                            bathy.elevation.values,
                            levels=(-1000, -100),
                            zorder=1.5
                            )
            except ValueError:
                print("Bathymetry deeper than specified levels.")
        add_ticks(ax1, extent)
        add_ticks(ax2, extent, label_left=False, label_right=True)

        # Plot gliders and argo floats
        rargs = {}
        rargs['argo'] = argo
        rargs['gliders'] = gliders
        rargs['transform'] = transform['data']  
        plot_regional_assets(ax1, **rargs)
        plot_regional_assets(ax2, **rargs)

        # Label the subplots
        ax1.set_title(ds1.model.upper(), fontsize=16, fontweight="bold")
        ax2.set_title(ds2.model.upper(), fontsize=16, fontweight="bold")

        # # Add EddyWatch Fronts
        # for ax in [ax1, ax2]:
        #     fronts = []
        #     for item in data['BZ_all'][0]:
        #         loop_y = item['y'].T
        #         loop_x = item['x'].T

        #         hf = ax.plot(loop_x, loop_y,
        #                     linestyle=item['LineStyle'][0],
        #                     color='black',
        #                     linewidth=3, 
        #                     transform=ccrs.PlateCarree(), 
        #                     zorder=120
        #                     )
        #         fronts.append(hf)

        #         # Add arrows
        #         start_lon = item['bx'].T
        #         start_lat = item['by'].T
        #         end_lon = item['tx'].T
        #         end_lat = item['ty'].T

        #         for count, _ in enumerate(start_lon):
        #             ax.arrow(
        #                 start_lon[count][0],
        #                 start_lat[count][0],
        #                 end_lon[count][0]-start_lon[count][0],
        #                 end_lat[count][0]-start_lat[count][0],
        #                 linewidth=0, 
        #                 head_width=0.2,
        #                 shape='full', 
        #                 fc='black', 
        #                 ec='black',
        #                 transform=ccrs.PlateCarree(),
        #                 zorder=130,
        #                 )
        # fronts.reverse()

        # Deal with the third axes
        h, l = ax2.get_legend_handles_labels()  # get labels and handles from ax1

        h_n = h #+ fronts
        l_n = l #+ ["EddyWatch - 1.5 knots (Active)", "EddyWatch - Eddy Remnants"]
        
        if (len(h) > 0) & (len(l) > 0):
            ax3.legend(h_n, l_n, ncol=cols, loc='center', fontsize=8)

            # Add title to legend
            t0 = []
            if isinstance(argo, pd.DataFrame):
                if not argo.empty:
                    t0.append(argo.index.min()[1])

            if isinstance(gliders, pd.DataFrame):
                if not gliders.empty:
                    t0.append(gliders.index.min()[1])

            if len(t0) > 0:
                t0 = min(t0).strftime('%Y-%m-%d %H:00:00')
            else:
                t0 = None
            legstr = f'Glider/Argo Search Window: {str(t0)} to {str(time)}'
            ax3.set_title(legstr, loc="center", fontsize=9, fontweight="bold", style="italic")
            # plt.figtext(0.5, 0.001, legstr, ha="center", fontsize=10, fontweight='bold')
        ax3.set_axis_off()

        # fig.tight_layout()
        # fig.subplots_adjust(top=0.88, wspace=.001)
        
        # Filled contour for each model variable
        m1 = ax1.contourf(ds1_depth["lon"], ds1_depth["lat"], mag_r, **qargs)
        m2 = ax2.contourf(ds2_depth["lon"], ds2_depth["lat"], mag_g, **qargs)

        # Set coarsening configs to a variable
        if 'coarsen' in cdict:
            coarsen = region['currents']['coarsen']
        else:
            coarsen['rtofs'] = 1
            coarsen['gofs'] = 1

        # Add streamlines
        s1 = map_add_currents(ax1, ds1_depth, coarsen=coarsen["rtofs"], **currents["kwargs"])
        s2 = map_add_currents(ax2, ds2_depth, coarsen=coarsen["gofs"], **currents["kwargs"])

        # Add EEZ
        if eez:
            eez1 = map_add_eez(ax1, zorder=1, color='red', linewidth=2, linestyle='-')
            eez2 = map_add_eez(ax2, zorder=1, color='red', linewidth=2, linestyle='-')

        if colorbar:
            cb = fig.colorbar(m1, ax=axs[:2], orientation="horizontal", shrink=.95, aspect=80)#, shrink=0.7, aspect=20*0.7)
            # cb = add_colorbar(axs[:2], m1, location="bottom")
            cb.ax.tick_params(labelsize=12)
            cb.set_label(f'Magnitude (m/s)', fontsize=12, fontweight="bold")

        # Create a string for the title of the plot
        title_time = time.strftime("%Y-%m-%d %H:%M:%S")
        title = f"Currents ({depth} m) - {title_time}\n"
        plt.suptitle(title, fontsize=22, fontweight="bold")

        # subplot 1
        # if depth == 0:
        # ax1.contour(ds1_depth['lon'], ds1_depth['lat'], mag_r, [0.771667], colors='k', linewidths=3, transform=ccrs.PlateCarree(), zorder=100)
        # ax2.contour(ds2_depth['lon'], ds2_depth['lat'], mag_g, [0.771667], colors='k', linewidths=3, transform=ccrs.PlateCarree(), zorder=100)    

        # Save figure
        fig.savefig(save_file_q, dpi=dpi, bbox_inches='tight', pad_inches=0.1)
        if colorbar:
            cb.remove()
            
        # Delete contour handles and remove colorbar axes to use figure
        s1.lines.remove(), s2.lines.remove()
        remove_quiver_handles(ax1), remove_quiver_handles(ax2)
        [x.remove() for x in m1.collections]
        [x.remove() for x in m2.collections]


def salinity_max(ds, extent, region_name,
                 limits=None,
                 bathy=None,
                 argo=None,
                 gliders=None,
                 eez=False,
                 cols=6,
                 path_save=os.getcwd(), 
                 transform=dict(map=ccrs.Mercator(), 
                                data=ccrs.PlateCarree()
                                ), 
                 figsize=(14,8),
                 dpi=150,
                 overwrite=False):
    
    # Convert ds.time value to a normal datetime
    time = pd.to_datetime(ds.time.values)
    
    # Formatter for time
    tstr_title = time.strftime('%Y-%m-%d %H:%M:%S')
    tstr_folder = time.strftime('%Y-%m-%dT%H%M%SZ')
    year = time.strftime("%Y")
    month = time.strftime("%m")

    # Generate filename
    fname = f"{path_save.name}_{tstr_folder}_salinity_max_{ds.model.lower()}.png"

    # Append salinity_max, year, and month to path_save
    path_save = path_save / 'salinity_max' / year / month

    save_file = path_save / fname
    
    if save_file.is_file():
        if not overwrite:
            print(f"{save_file} exists. Overwrite: False. Skipping.")
            return
        else:
            print(f"{save_file} exists. Overwrite: True. Replotting.")
        
    # Make sure path_save exists
    os.makedirs(path_save, exist_ok=True)
    
    print(f"Plotting Salinity Max of {region_name} at {tstr_title}")
    
    # Get the maximum salinity over the dimension 'depth'
    smax = ds['salinity'].max("depth")

    # Get the depth that the maximum salinity occurs at.
    # We find the indexes of the salinity maximum using the .argmax() method.
    # We use the .isel() method to select the depths that the salinity max occured.
    smax_depth = ds['salinity'].idxmax("depth")

    # Initialize figure    
    fig, _ = plt.subplot_mosaic(
        """
        RG
        LL
        """,
        figsize=figsize,
        layout="constrained",
        subplot_kw={
            'projection': transform['map']
            },
        gridspec_kw={
            # set the height ratios between the rows
            "height_ratios": [4, 1],
            # set the width ratios between the columns
            # # "width_ratios": [1],
            },
        )
    axs = fig.axes
    ax1 = axs[0] # rtofs
    ax2 = axs[1] # gofs
    ax3 = axs[2] # legend for argo/gliders

    # Plot gliders and argo floats
    rargs = {}
    rargs['argo'] = argo
    rargs['gliders'] = gliders
    rargs['transform'] = transform['data']  
    
    for ax in [ax1, ax2]:
        ax.set_extent(extent)

        # Add features to both map axes. Land, water, coastlines, etc.
        add_features(ax)
        
        # Add bathymetry lines
        if bathy:
            add_bathymetry(ax,
                           bathy.longitude.values, 
                           bathy.latitude.values, 
                           bathy.elevation.values,
                           levels=(-1000, -100),
                           zorder=1.5
                           )
        # Add eez lines
        if eez:
            map_add_eez(ax, color='white', zorder=10)

        plot_regional_assets(ax, **rargs)

    # Add ticks
    add_ticks(ax1, extent, label_left=True)
    add_ticks(ax2, extent, label_left=False, label_right=True)

    # Deal with the third axes
    h, l = ax1.get_legend_handles_labels()  # get labels and handles from ax1
    if (len(h) > 0) & (len(l) > 0):
        
        # Add handles to legend
        ax3.legend(h, l, ncol=cols, loc='center', fontsize=8)

        # Add title to legend
        t0 = []
        if isinstance(argo, pd.DataFrame):
            if not argo.empty:
                t0.append(argo.index.min()[1])

        if isinstance(gliders, pd.DataFrame):
            if not gliders.empty:
                t0.append(gliders.index.min()[1])

        if len(t0) > 0:
            t0 = min(t0).strftime('%Y-%m-%d %H:00:00')
        else:
            t0 = None
        legstr = f'Glider/Argo Search Window: {t0} to {str(time)}'
        ax3.set_title(legstr, loc="center", fontsize=9, fontweight="bold", style='italic')
    ax3.set_axis_off()

    # Calculate contours
    if limits:
        salt_min = limits[0]
        salt_max = limits[1]
        salt_stride = limits[2]
    else:
        # Calculate the colorbar limits automatically
        percentiles = [.25, .75]
        salinity_quantile = smax.quantile(percentiles)
        salt_min = np.floor(salinity_quantile[0])
        salt_max = np.ceil(salinity_quantile[1])
        salt_stride = .1

    # Calculate salinity contours
    levels = np.arange(
        salt_min,
        salt_max+salt_stride,
        salt_stride
        )

    # Calculate depth contours
    depths = np.arange(40, 200, 20)

    # Salinity Max Plot
    h1 = ax1.contourf(smax['lon'], smax['lat'], smax,
                        levels=levels, 
                        extend="both",
                        cmap=cmocean.cm.haline,
                        transform=transform['data'])

    # Salinity Max Depth Plot
    h2 = ax2.contourf(smax_depth['lon'], smax_depth['lat'], smax_depth, 
                        levels=depths,
                        extend="both",
                        cmap=cmocean.cm.deep,
                        transform=transform['data'])
    
    # Add colorbar to first axes
    cb = add_colorbar(ax1, h1)
    cb.ax.tick_params(labelsize=10)
    cb.set_label(f'Salinity', fontsize=11, fontweight="bold")
    
    # Add colorbar to second axes
    cb = add_colorbar(ax2, h2)
    cb.ax.tick_params(labelsize=10)
    cb.set_label(f'Depth (m)', fontsize=11, fontweight="bold")

    # Set title for each axes
    ax1.set_title("Salinity Maximum", fontsize=16, fontweight='bold')
    ax2.set_title("Depth of Salinity Maximum", fontsize=16, fontweight='bold')
    fig.suptitle(f"{ds.model.upper()} - {time.strftime(tstr_title)}", fontweight="bold", fontsize=20)

    # # Hide x labels and tick labels for top plots and y ticks for right plots.
    # for axs in ax.flat:
    #     axs.label_inner()

    # fig.tight_layout()
    # fig.subplots_adjust(top=0.95)
    
    export_fig(path_save, fname, dpi=dpi)
    plt.close()

def plot_storm_track(ax, storm, zorder=90, proj=proj['data']):
    cone = storm['cone']
    storm_lon = storm['lon']
    storm_lat = storm['lat']
    
    #Plot cone
    cone_2d = cone['cone']
    cone_2d = ndimage.gaussian_filter(cone_2d, sigma=0.5, order=0)

    # Cone shading
    ax.contourf(cone['lon2d'], cone['lat2d'], cone_2d, [0.9,1.1], 
                colors=['#ffffff', '#ffffff'],
                alpha=.25,
                zorder=zorder+1,
                transform=proj)

    # Cone outline
    ax.contour(cone['lon2d'], cone['lat2d'], cone_2d, [0.9], 
                linewidths=1.5,
                colors=['k'],
                zorder=zorder+1,
                transform=proj)

    #Plot forecast center line & account for dateline crossing
    ax.plot(cone['center_lon'], cone['center_lat'],
            color='k',
            linewidth=2.0,
            zorder=zorder+2,
            transform=proj
            )

    #Plot previous track
    ax.plot(storm_lon, storm_lat, color='red',
            linestyle='dotted',
            linewidth=1,
            zorder=zorder+1,
            transform=proj)

def salinity_max_comparison(ds1, ds2, extent, region_name,
                            limits=None,
                            bathy=None,
                            argo=None,
                            gliders=None,
                            eez=False,
                            cols=6,
                            path_save=os.getcwd(), 
                            transform=dict(map=ccrs.Mercator(), 
                                           data=ccrs.PlateCarree()), 
                            figsize=(14,8),
                            dpi=150,
                            overwrite=False,
                            storms=None):
    
    # Convert ds.time value to a normal datetime
    time = pd.to_datetime(ds1.time.values)
    
    # Formatter for time
    tstr_title = time.strftime('%Y-%m-%d %H:%M:%S')
    tstr_folder = time.strftime('%Y-%m-%dT%H%M%SZ')
    year = time.strftime("%Y")
    month = time.strftime("%m")

    # Generate filename
    fname_max = f"{path_save.name}_{tstr_folder}_salinity_max_{ds1.model.lower()}-{ds2.model.lower()}.png"
    fname_depth = f"{path_save.name}_{tstr_folder}_salinity_max_{ds1.model.lower()}-{ds2.model.lower()}-depth.png"

    # Append salinity_max, year, and month to path_save
    path_save = path_save / 'salinity_max' / year / month

    save_file = path_save / fname_max
    
    if save_file.is_file():
        if not overwrite:
            print(f"{save_file} exists. Overwrite: False. Skipping.")
            return
        else:
            print(f"{save_file} exists. Overwrite: True. Replotting.")
        
    # Make sure path_save exists
    os.makedirs(path_save, exist_ok=True)
    
    print(f"Plotting Salinity Max of {region_name} at {tstr_title}")
    
    # Get the maximum salinity over the dimension 'depth'
    smax1 = ds1['salinity'].max("depth")
    smax2 = ds2['salinity'].max("depth")

    # Initialize figure    
    fig, _ = plt.subplot_mosaic(
        """
        RG
        LL
        """,
        figsize=figsize,
        layout="constrained",
        subplot_kw={
            'projection': transform['map']
            },
        gridspec_kw={
            # set the height ratios between the rows
            "height_ratios": [4, 1],
            # set the width ratios between the columns
            # # "width_ratios": [1],
            },
        )
    axs = fig.axes
    ax1 = axs[0] # rtofs
    ax2 = axs[1] # gofs
    ax3 = axs[2] # legend for argo/gliders

    ax1.set_extent(extent)
    ax2.set_extent(extent)

    # Argo/Glider Data Dicts
    rargs = {}
    rargs['argo'] = argo
    rargs['gliders'] = gliders
    rargs['transform'] = transform['data']
    

    for ax in [ax1, ax2]:
        # Make the map pretty  
        add_features(ax)# zorder=0)
        
        # Add features to both map axes. Land, water, coastlines, etc.
        add_features(ax)

        # Add bathymetry lines
        # if bathy:
        #     add_bathymetry(ax,
        #                    bathy.longitude.values, 
        #                    bathy.latitude.values, 
        #                    bathy.elevation.values,
        #                    levels=(-1000, -100),
        #                    zorder=1.5
        #                    )
            
        # Add eez lines
        if eez:
            map_add_eez(ax, color='white', zorder=10, linewidth=1)

        # Plot gliders and argo floats 
        plot_regional_assets(ax, **rargs)

        if storms:
            for s in storms.keys():
                storms['track']
                plot_storm_track(ax, lon, lat, cone)
                
    # Add ticks
    add_ticks(ax1, extent, label_left=True)
    add_ticks(ax2, extent, label_left=False, label_right=True)

    # Deal with the third axes
    h, l = ax1.get_legend_handles_labels()  # get labels and handles from ax1
    if (len(h) > 0) & (len(l) > 0):
        
        # Add handles to legend
        ax3.legend(h, l, ncol=cols, loc='center', fontsize=8)

        # Add title to legend
        t0 = []
        if isinstance(argo, pd.DataFrame):
            if not argo.empty:
                t0.append(argo.index.min()[1])

        if isinstance(gliders, pd.DataFrame):
            if not gliders.empty:
                t0.append(gliders.index.min()[1])

        if len(t0) > 0:
            t0 = min(t0).strftime('%Y-%m-%d %H:00:00')
        else:
            t0 = None
        legstr = f'Glider/Argo Search Window: {t0} to {str(time)}'
        ax3.set_title(legstr, loc="center", fontsize=9, fontweight="bold", style='italic')
    ax3.set_axis_off()

    # Calculate contours
    if limits:
        salt_min = limits[0]
        salt_max = limits[1]
        salt_stride = limits[2]
    else:
        # Calculate the colorbar limits automatically
        percentiles = [.25, .75]
        salinity_quantile = smax1.quantile(percentiles)
        salt_min = np.floor(salinity_quantile[0])
        salt_max = np.ceil(salinity_quantile[1])
        salt_stride = .1

    # Calculate salinity contours
    levels = np.arange(
        salt_min,
        salt_max+salt_stride,
        salt_stride
        )

    # Salinity Max Plot
    h1 = ax1.contourf(smax1['lon'], smax1['lat'], smax1,
                      levels=levels, 
                      extend="both",
                      cmap=cmocean.cm.haline,
                      transform=transform['data'])

    h2 = ax2.contourf(smax2['lon'], smax2['lat'], smax2,
                      levels=levels, 
                      extend="both",
                      cmap=cmocean.cm.haline,
                      transform=transform['data'])

    # Add colorbar to first axes
    cb = fig.colorbar(h1, ax=axs[:2], orientation="horizontal", shrink=.95, aspect=80)#, shrink=0.7, aspect=20*0.7)
    cb.ax.tick_params(labelsize=10)
    cb.set_label('Salinity', fontsize=11, fontweight="bold")
    
    # Set title for each axes
    ax1.set_title(f"{ds1.model.upper()}", fontsize=16, fontweight='bold')
    ax2.set_title(f"{ds2.model.upper()}", fontsize=16, fontweight='bold')
    fig.suptitle(f"Salinity Maximum - {time.strftime(tstr_title)}", fontweight="bold", fontsize=20)
    
    export_fig(path_save, fname_max, dpi=dpi)

    # Remove handles so we can reuse figure
    # Delete contour handles 
    [x.remove() for x in h1.collections] # axes 1
    [x.remove() for x in h2.collections] # axes 2

    cb.remove()

    # Get the depth that the maximum salinity occurs at.
    # We find the indexes of the salinity maximum using the .argmax() method.
    # We use the .isel() method to select the depths that the salinity max occured.
    smax1_depth = ds1['salinity'].idxmax("depth")
    smax2_depth = ds2['salinity'].idxmax("depth")

    # Calculate depth contours
    depths = np.arange(40, 200, 20)
    
    # Salinity Max Depth Plot
    h1 = ax1.contourf(smax1_depth['lon'], smax1_depth['lat'], smax1_depth, 
                      levels=depths,
                      extend="both",
                      cmap=cmocean.cm.deep,
                      transform=transform['data'])
    h2 = ax2.contourf(smax2_depth['lon'], smax2_depth['lat'], smax2_depth, 
                      levels=depths,
                      extend="both",
                      cmap=cmocean.cm.deep,
                      transform=transform['data'])

    # Add colorbar to first axes
    cb = fig.colorbar(h1, ax=axs[:2], orientation="horizontal", shrink=.95, aspect=80)#, shrink=0.7, aspect=20*0.7)
    cb.ax.tick_params(labelsize=10)
    cb.set_label('Depth (m)', fontsize=11, fontweight="bold")
    
    # Set title for each axes
    ax1.set_title(f"{ds1.model.upper()}", fontsize=16, fontweight='bold')
    ax2.set_title(f"{ds2.model.upper()}", fontsize=16, fontweight='bold')
    fig.suptitle(f"Depth of Salinity Maximum - {time.strftime(tstr_title)}", fontweight="bold", fontsize=20)
    
    export_fig(path_save, fname_depth, dpi=dpi)
    
    plt.close()


def plot_ohc(ds1, ds2, extent, region_name,
             limits=None,
             bathy=None,
             argo=None,
             gliders=None,
             eez=False,
             cols=6,
             path_save=os.getcwd(), 
             transform=dict(map=ccrs.Mercator(), 
                            data=ccrs.PlateCarree()
                            ), 
             figsize=(14,8),
             dpi=150,
             overwrite=False
             ):

    # Convert ds.time value to a normal datetime
    time = pd.to_datetime(ds1.time.values)
    
    # Formatter for time
    tstr_title = time.strftime('%Y-%m-%d %H:%M:%S')
    tstr_folder = time.strftime('%Y-%m-%dT%H%M%SZ')
    year = time.strftime("%Y")
    month = time.strftime("%m")

    # Generate filename
    fname = f"{path_save.name}_{tstr_folder}_heat_content_{ds1.model.lower()}-{ds2.model.lower()}.png"

    # Append salinity_max, year, and month to path_save
    path_save = path_save / 'ocean_heat_content' / year / month

    save_file = path_save / fname
    
    if save_file.is_file():
        if not overwrite:
            print(f"{save_file} exists. Overwrite: False. Skipping.")
            return
        else:
            print(f"{save_file} exists. Overwrite: True. Replotting.")    
    
    # Make sure path_save exists
    os.makedirs(path_save, exist_ok=True)
    
    print(f"Plotting Ocean Heat Content of {region_name} at {tstr_title}")
 
    # Initialize figure    
    fig, _ = plt.subplot_mosaic(
        """
        RG
        LL
        """,
        figsize=figsize,
        layout="constrained",
        subplot_kw={
            'projection': transform['map']
            },
        gridspec_kw={
            # set the height ratios between the rows
            "height_ratios": [4, 1],
            # set the width ratios between the columns
            # # "width_ratios": [1],
            },
        )
    axs = fig.axes
    ax1 = axs[0] # rtofs
    ax2 = axs[1] # gofs
    ax3 = axs[2] # legend for argo/gliders

    # Setup keyword arguments dictionary for plot_regional_assets
    rargs = {}
    rargs['argo'] = argo
    rargs['gliders'] = gliders
    rargs['transform'] = transform['data']  
    
    for ax in [ax1, ax2]:
        create(extent, ax=ax1, ticks=False)
        create(extent, ax=ax2, ticks=False)
        
        # Make the map pretty  
        # add_features(ax)# zorder=0)

        # Add features to both map axes. Land, water, coastlines, etc.
        # add_features(ax)

        # Add bathymetry lines
        if bathy:       
            add_bathymetry(ax,
                           bathy.longitude.values, 
                           bathy.latitude.values, 
                           bathy.elevation.values,
                           levels=(-1000, -100),
                           zorder=1.5
                           )
        # Add eez lines
        if eez:
            map_add_eez(ax, color='black', zorder=10)

        # Plot gliders and argo floats
        plot_regional_assets(ax, **rargs)

    # Add ticks
    add_ticks(ax1, extent, label_left=True)
    add_ticks(ax2, extent, label_left=False, label_right=True)

    # Calculate contours
    if limits:
        ohc_min = limits[0]
        ohc_max = limits[1]
        ohc_stride = limits[2]
    else:
        # Calculate the colorbar limits automatically
        percentiles = [.25, .75]
        quantile = ds1['ohc'].quantile(percentiles)
        ohc_min = np.floor(quantile[0])
        ohc_max = np.ceil(quantile[1])
        ohc_stride = 10

    # Calculate salinity contours
    levels = np.arange(
        ohc_min,
        ohc_max+ohc_stride,
        ohc_stride
        )

    cmap = cmocean.cm.thermal

    # Ocean Heat Content Plot
    h1 = ax1.contourf(ds1['lon'], ds1['lat'], ds1['ohc'],
                      levels=levels, 
                      extend="max",
                      cmap=cmap,
                      transform=transform['data'])
    h3 = ax1.contour(ds1['lon'], 
                     ds1['lat'], 
                     ds1['ohc'], 
                     [60],
                     linestyles='-',
                     colors=['silver'],
                     linewidths=1,
                     alpha=1,
                     transform=ccrs.PlateCarree(),
                     zorder=101)


    # Ocean Heat Content Plot
    h2 = ax2.contourf(ds2['lon'], ds2['lat'], ds2['ohc'], 
                      levels=levels,
                      extend="max",
                      cmap=cmap,
                      transform=transform['data']
                      )
    ax2.contour(ds2['lon'], 
                ds2['lat'], 
                ds2['ohc'], 
                [60],
                linestyles='-',
                colors=['silver'],
                linewidths=1,
                alpha=1,
                transform=ccrs.PlateCarree(),
                zorder=101)
    h0 = []
    l0 = []
    
    h0.append(mlines.Line2D([], [], linestyle='-', color='silver', alpha=1, linewidth=1))
    l0.append('60 kJ cm-2')
    # h0.append(mlines.Line2D([], [], linestyle='-', color='white', alpha=1, linewidth=1))
    # l0.append('Past 5 days')
    leg1 = ax1.legend(h0, l0, loc='upper left', fontsize=9)
    leg2 = ax2.legend(h0, l0, loc='upper left', fontsize=9)

    # Deal with the third axes
    h, l = ax1.get_legend_handles_labels()  # get labels and handles from ax1
    if (len(h) > 0) & (len(l) > 0):
        # h.append(mlines.Line2D([], [], linestyle='-', color='silver', alpha=1, linewidth=1))
        # l.append('60 kJ cm-2')
        # h.append(mlines.Line2D([], [], linestyle='-', color='white', alpha=1, linewidth=1))
        # l.append('Tracks - Past 5 days')
        # h.append()
        
        # Add handles to legend
        ax3.legend(h, l, ncol=cols, loc='center', fontsize=8)

        # Add title to legend
        t0 = []
        if isinstance(argo, pd.DataFrame):
            if not argo.empty:
                t0.append(argo.index.min()[1])

        if isinstance(gliders, pd.DataFrame):
            if not gliders.empty:
                t0.append(gliders.index.min()[1])

        if len(t0) > 0:
            t0 = min(t0).strftime('%Y-%m-%d %H:00:00')
        else:
            t0 = None
        legstr = f'Glider/Argo Search Window: {t0} to {str(time)}'
        ax3.set_title(legstr, loc="center", fontsize=9, fontweight="bold", style='italic')
    ax3.set_axis_off()
    # plt.gca().add_artist(leg1)
    # plt.gca().add_artist(leg2)

    leg1.set_zorder(10001)
    leg2.set_zorder(10001)    
    
    # Add colorbar to first axes
    cb = fig.colorbar(h1, ax=axs[:2], orientation="horizontal", shrink=.95, aspect=80)#, shrink=0.7, aspect=20*0.7)
    cb.ax.tick_params(labelsize=12)
    cb.set_label('kJ/cm^2', fontsize=12, fontweight="bold")

    # Set title for each axes
    ax1.set_title(f"{ds1.model.upper()}", fontsize=16, fontweight='bold')
    ax2.set_title(f"{ds2.model.upper()}", fontsize=16, fontweight='bold')
    fig.suptitle(f"Ocean Heat Content - {time.strftime(tstr_title)}", fontweight="bold", fontsize=20)

    # from ioos_model_comparisons.plotting_hurricanes import plot_storms
    # from tropycal import realtime
    # import datetime as dt
    
    # realtime_obj = realtime.Realtime()

    # realtime_obj.list_active_storms(basin='north_atlantic')

    # # realtime_obj.plot_summary(domain={'w':-100,'e':-10,'s':4,'n':60})

    # #Get realtime forecasts
    # forecasts = []
    # for key in realtime_obj.storms:
    #     if realtime_obj[key].invest == False:
    #         try:
    #             forecasts.append(realtime_obj.get_storm(key).get_forecast_realtime(True))
    #         except:
    #             forecasts.append({})
    #     else:
    #         forecasts.append({})
    # forecasts = [entry if 'init' in entry.keys() and (dt.utcnow() - entry['init']).total_seconds() / 3600.0 <= 12 else {} for entry in forecasts]
    # storms = [realtime_obj.get_storm(key) for key in realtime_obj.storms]

    # plot_storms(ax, storms, forecasts, zorder=80)


    # # Hide x labels and tick labels for top plots and y ticks for right plots.
    # for axs in ax.flat:
    #     axs.label_inner()

    # fig.tight_layout()
    # fig.subplots_adjust(top=0.95)
    
    export_fig(path_save, fname, dpi=dpi)
    plt.close()


def plot_model_region_ts_vs_speed(ds,
                                  region,
                                  depths,
                                  bathy=None,
                                  argo=None,
                                  gliders=None,
                                  cols=6,
                                  transform=dict(map=proj['map'], 
                                                 data=proj['data']
                                                 ),
                                  path_save=os.getcwd(),
                                  figsize=(14,8),
                                  dpi=150,
                                  overwrite=False
                                  ):
    
    # Convert ds.time value to a normal datetime
    time = pd.to_datetime(ds.time.data)
    extent = region['extent']

    # Formatter for time
    tstr_title = time.strftime('%Y-%m-%d %H:%M:%S')

    # Create figure
    grid = """
    RG
    LL
    """

    fig, _ = plt.subplot_mosaic(
        grid,
        figsize=figsize,
        layout="constrained",
        subplot_kw={
            'projection': transform['map']
            },
        gridspec_kw={
            # set the height ratios between the rows
            "height_ratios": [4, 1],
            # set the width ratios between the columns
            # # "width_ratios": [1],
            },
        )
    axs = fig.axes
    ax1 = axs[0] # Model 1
    ax2 = axs[1] # Model 2
    ax3 = axs[2] # Legend for argo/gliders
    
    # Make the map pretty  
    add_features(ax1)# zorder=0)
    add_features(ax2)# zorder=0)  

    # Add bathymetry lines
    if bathy:
        add_bathymetry(ax1,
                           bathy.longitude.values, 
                           bathy.latitude.values, 
                           bathy.elevation.values,
                           levels=(-1000, -100),
                           zorder=1.5)
        add_bathymetry(ax2,
                           bathy.longitude.values, 
                           bathy.latitude.values, 
                           bathy.elevation.values,
                           levels=(-1000, -100),
                           zorder=1.5)

    # Add ticks
    add_ticks(ax1, extent, label_left=True)
    add_ticks(ax2, extent, label_left=False, label_right=True)

    # Plot gliders and argo floats
    rargs = {}
    rargs['argo'] = argo
    rargs['gliders'] = gliders
    rargs['transform'] = transform['data']  
    plot_regional_assets(ax1, **rargs)
    plot_regional_assets(ax2, **rargs)

    # Label the subplots
    # ax1.set_title(ds1.model, fontsize=16, fontweight="bold")
    # ax2.set_title(ds2.model, fontsize=16, fontweight="bold")
    txt = plt.suptitle("", fontsize=22, fontweight="bold")
    
    # Deal with the third axes
    h, l = ax1.get_legend_handles_labels()  # get labels and handles from ax1
    if (len(h) > 0) & (len(l) > 0):
        
        # Add handles to legend
        legend = ax3.legend(h, l, ncol=cols, loc='center', fontsize=8)

        # Add title to legend
        t0 = []
        if isinstance(argo, pd.DataFrame):
            if not argo.empty:
                t0.append(argo.index.min()[1])

        if isinstance(gliders, pd.DataFrame):
            if not gliders.empty:
                t0.append(gliders.index.min()[1])

        if len(t0) > 0:
            t0 = min(t0).strftime('%Y-%m-%d %H:00:00')
        else:
            t0 = None
        legstr = f'Glider/Argo Search Window: {t0} to {str(time)}'
        ax3.set_title(legstr, loc="center", fontsize=9, fontweight="bold", style='italic')
        legend._legend_box.sep = 1
        # plt.figtext(0.5, 0.001, legstr, ha="center", fontsize=10, fontweight='bold')
    ax3.set_axis_off()

    # Iterate through the variables to be plotted for each region. 
    # This dict contains information on what variables and depths to plot. 
    for k, v in region["variables"].items():  
        var_str = ' '.join(k.split('_')).title()      
        for item in v:
            if item['depth'] in depths:
                print(f"Plotting {k} @ {item['depth']}")
                depth = item['depth']
                tds = ds[k].sel(depth=depth)
                _, speed = uv2spdir(tds['u'], tds['v'])
                # ds_mag = ds[k].sel(depth=depth, method='nearest')
                
                # Create subdirectory for depth under variable subdirectory
                save_dir_final = path_save / f"{k}_{depth}m" / time.strftime('%Y/%m')
                os.makedirs(save_dir_final, exist_ok=True)

                # Create a file name to save the plot as
                sname = f'{"-".join(region["folder"].split("_"))}_{time.strftime("%Y-%m-%dT%H%M%SZ")}_{k}-{depth}m_{ds.model.lower()}'
                save_file = save_dir_final / f"{sname}.png"

                if save_file.is_file():
                    if not overwrite:
                        print(f"{save_file} exists. Overwrite: False. Skipping.")
                        continue
                    else:
                        print(f"{save_file} exists. Overwrite: True. Replotting.")
                        
                # Add the super title (title for both subplots)
                txt.set_text(f"{var_str} ({depth} m) - {tstr_title}\n")
            
                # Filled contour for each model variable
                # Check if ndims are 1, transform_first requires 2d array
                if (tds['lon'].ndim == 1) & tds['lat'].ndim == 1:
                    rlons, rlats = np.meshgrid(tds['lon'], tds['lat'])
                else:
                    rlons = tds['lon']
                    rlats = tds['lat']

                # Plot first subplot
                # Create dictionary for variable argument inputs for contourf
                vargs = {}
                vargs['transform'] = transform['data']
                vargs['transform_first'] = True
                vargs['cmap'] = cmaps(ds[k].name)
                vargs['extend'] = "both"

                if 'limits' in item:
                    vargs['vmin'] = item['limits'][0]
                    vargs['vmax'] = item['limits'][1]
                    vargs['levels'] = np.arange(vargs['vmin'], vargs['vmax'], item['limits'][2])
                h1 = ax1.contourf(rlons, rlats, tds.squeeze(), **vargs)
                cb1 = fig.colorbar(h1, ax=ax1, orientation="horizontal", shrink=.95, aspect=80)#, shrink=0.7, aspect=20*0.7)
                cb1.ax.tick_params(labelsize=12)
                cb1.set_label(f'{k.title()} ({tds.units})', fontsize=12, fontweight="bold")

                # Plot second subplot
                margs = {}
                margs['transform'] = transform['data']
                margs['transform_first'] = True
                margs['cmap'] = cmocean.cm.speed
                margs['extend'] = "both"
                margs['levels'] = np.arange(0, .6, .05)
                h2 = ax2.contourf(rlons, rlats, speed, **margs)
                cb2 = fig.colorbar(h2, ax=ax2, orientation="horizontal", shrink=.95, aspect=80)#, shrink=0.7, aspect=20*0.7)
                cb2.ax.tick_params(labelsize=12)
                cb2.set_label(f'Speed (m/s)', fontsize=12, fontweight="bold")

                # Plot streamlines
                sargs = {}
                sargs["transform"] = transform['data']
                sargs["density"] = 3
                sargs["linewidth"] = .75
                sargs["color"] = 'red'
                sargs['zorder'] = 100
                    
                q = ax2.streamplot(rlons.values, rlats.values, 
                                   tds['u'].values, tds['v'].values,
                                   **sargs)
                fig.savefig(save_dir_final / sname, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

                # Remove quiver handles from each axes
                q.lines.remove()
                remove_quiver_handles(ax2)

                # Remove handles so we can reuse figure
                # Delete contour handles 
                [x.remove() for x in h1.collections] # axes 1
                [x.remove() for x in h2.collections] # axes 2

                cb1.remove()
                cb2.remove()    

                plt.close()
