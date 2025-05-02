import xarray as xr
import geopandas as gpd
import os
import pandas as pd
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

def load_ucdb_city(root, city):
    """
    Load and filter a city shapefile from the Urban Centre Database (UCDB).

    Parameters:
    root (str): The root directory where the shapefile is located.
    city (str): The name of the city to load.

    Returns:
    gpd.GeoDataFrame: A GeoDataFrame containing the filtered city shapefile.
    """
    ucdb_info = gpd.read_file(root + '/GHS_FUA_UCD/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg')
    ucdb_city = ucdb_info.query(f'UC_NM_MN =="{city}"').to_crs(crs='EPSG:4326')
    if city == 'London':
        ucdb_city = ucdb_city[ucdb_city['CTR_MN_NM'] == 'United Kingdom']
    if city == 'Santiago':
        ucdb_city = ucdb_city[ucdb_city['CTR_MN_NM'] == 'Chile']
    if city == 'Barcelona':
        ucdb_city = ucdb_city[ucdb_city['CTR_MN_NM'] == 'Spain']
    if city == 'Dhaka':
        ucdb_city = ucdb_city[ucdb_city['CTR_MN_NM'] == 'Bangladesh']
    if city == 'Naples':
        ucdb_city = ucdb_city[ucdb_city['CTR_MN_NM'] == 'Italy']
    return ucdb_city


def traverseDir(root):
    for (dirpath, dirnames, filenames) in os.walk(root):
        for file in filenames:
            if file.endswith(('.nc')):
                yield os.path.join(dirpath, file)

def fix_360_longitudes(
    dataset: xr.Dataset, lonname: str = "lon"
) -> xr.Dataset:
    """
    Fix longitude values.

    Function to transform datasets where longitudes are in (0, 360) to (-180, 180).

    Parameters
    ----------
    dataset (xarray.Dataset): data stored by dimensions
    lonname (str): name of the longitude dimension

    Returns
    -------
    dataset (xarray.Dataset): data with the new longitudes
    """
    lon = dataset[lonname]
    if lon.max().values > 180 and lon.min().values >= 0:
        dataset[lonname] = dataset[lonname].where(lon <= 180, other=lon - 360)
    return dataset

def kelvin2degC(ds, variable):
    """
    Convert a variable's units from Kelvin to Celsius in an xarray Dataset.

    Parameters:
    ds (xr.Dataset): The dataset containing the variable to convert.
    variable (str): The name of the variable to convert.

    Returns:
    xr.Dataset: The dataset with the variable converted to Celsius.
    """
    if ds[variable].attrs.get('units') == 'K':
        ds[variable] = ds[variable] -273.15
        ds[variable].attrs['units'] = 'degC'
        
    return ds


def fix_360_longitudes(
    dataset: xr.Dataset, lonname: str = "lon"
) -> xr.Dataset:
    """
    Fix longitude values.

    Function to transform datasets where longitudes are in (0, 360) to (-180, 180).

    Parameters
    ----------
    dataset (xarray.Dataset): data stored by dimensions
    lonname (str): name of the longitude dimension

    Returns
    -------
    dataset (xarray.Dataset): data with the new longitudes
    """
    lon = dataset[lonname]
    if lon.max().values > 180 and lon.min().values >= 0:
        dataset[lonname] = dataset[lonname].where(lon <= 180, other=lon - 360)
    return dataset

def plot_urban_polygon(ds, ax):
    '''
    '''
    # Assume the mask is in the 'urmask' variable
    mask = ds['urmask']
    lon2d = mask.lon.values
    lat2d = mask.lat.values
    
    # Create lists to store polygons for urban areas (mask == 1) and non-urban areas (mask == 0)
    urban_polygons = []
    non_urban_polygons = []
    
    if lon2d.ndim == 1:
        dist_lat = abs(lat2d[1] - lat2d[0]) / 2
        dist_lon = abs(lon2d[1] - lon2d[0]) / 2
        for lon in range(len(lon2d)):
            for lat in range(len(lat2d)):
                square = Polygon([
                    (round(lon2d[lon] - dist_lon, 3), round(lat2d[lat] - dist_lat, 3)),  # bottom-left corner
                    (round(lon2d[lon] + dist_lon, 3), round(lat2d[lat] - dist_lat, 3)),  # bottom-right corner
                    (round(lon2d[lon] + dist_lon, 3), round(lat2d[lat] + dist_lat, 3)),  # top-right corner
                    (round(lon2d[lon] - dist_lon, 3), round(lat2d[lat] + dist_lat, 3)),  # top-left corner
                ])

                # Add the polygon to the corresponding list
                if mask[lat, lon] == 1:
                    urban_polygons.append(square)
                elif mask[lat, lon] == 0:
                    non_urban_polygons.append(square)
    else:
        dist_lat = abs(lat2d[1, 0] - lat2d[0, 0])/2
        dist_lon = abs(lon2d[0, 1] - lon2d[0, 0])/2

        # Iterate through the mask and generate polygons for urban (1) and non-urban (0) cells
        for lat in range(mask.shape[0] - 1):  # Avoid the last index to prevent out-of-bounds errors
            for lon in range(mask.shape[1] - 1):
                # Create a polygon using the 2D lat/lon coordinates of the cell corners
                if pd.isnull(mask.lon[lat, lon]): # If cell contains nans continue
                    continue
                square = Polygon([
                    (mask.lon[lat, lon] - dist_lon, mask.lat[lat, lon] - dist_lat),          # bottom-left corner
                    (mask.lon[lat, lon + 1] - dist_lon, mask.lat[lat, lon + 1] - dist_lat),  # bottom-right corner
                    (mask.lon[lat + 1, lon + 1] - dist_lon, mask.lat[lat + 1, lon + 1] - dist_lat),  # top-right corner
                    (mask.lon[lat + 1, lon] - dist_lon, mask.lat[lat + 1, lon] - dist_lat),  # top-left corner
                ])
                # Add the polygon to the corresponding list
                if mask[lat, lon] == 1:
                    urban_polygons.append(square)
                elif mask[lat, lon] == 0:
                    non_urban_polygons.append(square)
                    
    # Unite all adjacent polygons for urban (mask == 1) and non-urban (mask == 0)
    unified_urban_polygon = unary_union(urban_polygons)
    unified_non_urban_polygon = unary_union(non_urban_polygons)
    # Create GeoDataFrames for the urban and non-urban polygons
    # CRS 'EPSG:4326' specifies the WGS 84 coordinate system, which is widely used for global GPS coordinates (lat/lon)
    gdf_urban = gpd.GeoDataFrame(geometry=[unified_urban_polygon])
    gdf_non_urban = gpd.GeoDataFrame(geometry=[unified_non_urban_polygon])
    # Plot the boundary of the unified non-urban polygon (in blue)
    gdf_non_urban.boundary.plot(ax=ax,color='b', zorder=1, linewidth=2)
    # Plot the boundary of the unified urban polygon (in red) on top of the non-urban
    gdf_urban.boundary.plot(ax=ax,  color='red', zorder=100, linewidth=2)

    return(gdf_urban, gdf_non_urban)

def plot_urban_borders(ds, ax, alpha = 1, linewidth = 2):
    """
    Plot the borders of urban areas on a map.

    Parameters:
    ds (xr.Dataset): The dataset containing longitude, latitude, and urban area data.
    ax (matplotlib.axes._subplots.AxesSubplot): The matplotlib axes on which to plot.

    """
    lon2d = ds.lon.values
    lat2d = ds.lat.values
    dist_lat = lat2d[1, 0] - lat2d[0, 0]
    dist_lon = lon2d[0, 1] - lon2d[0, 0]
    dist_latlon = lat2d[0 ,1] - lat2d[0, 0]
    dist_lonlat = lon2d[1, 0] - lon2d[0, 0]
    # Overlay the cell borders and handle NaNs
    for i in range(len(ds.lat)-1):
        for j in range(len(ds.lon)-1):
            lons = [lon2d[i, j], lon2d[i, j+1], lon2d[i+1, j+1], lon2d[i+1, j], lon2d[i, j]]
            lats = [lat2d[i, j], lat2d[i, j+1], lat2d[i+1, j+1], lat2d[i+1, j], lat2d[i, j]]

            lons = lons - abs(lon2d[i, j] - lon2d[i, j+1])/2              
            lats = lats - abs(lat2d[i, j] - lat2d[i+1, j])/2
            
            data_cell = ds['urmask'].values[i, j]
            
            if data_cell == 1:
                ax.plot(lons, lats, color='grey', zorder = 100, linewidth = linewidth, alpha = alpha)
            elif data_cell == 0:
                ax.plot(lons, lats, color='green', zorder = 1, linewidth = linewidth, alpha = alpha)

     # Plot the rightmost column
    for i in range(len(ds.lat) - 1):
        lons = [lon2d[i, -1], lon2d[i + 1, -1], lon2d[i + 1, -1] + dist_lon, lon2d[i, -1] + dist_lon, lon2d[i, -1]]
        lats = [lat2d[i, -1], lat2d[i + 1, -1], lat2d[i + 1, -1] + dist_latlon, lat2d[i, -1] + dist_latlon, lat2d[i, -1]]

        lons = lons - abs(lon2d[i, -1] - lon2d[i, -1])/2 - dist_lon/2            
        lats = lats - abs(lat2d[i, -1] - lat2d[i+1, -1])/2 
        
        data_cell = ds['urmask'].values[i, -1]
    
        if data_cell == 1:
            ax.plot(lons, lats, color='grey', zorder=100, linewidth = linewidth, alpha = alpha)
        elif data_cell == 0:
            ax.plot(lons, lats, color='green', zorder=1, linewidth = linewidth, alpha = alpha)
    
    # Plot the topmost row
    for j in range(len(ds.lon) - 1):
        lons = [lon2d[-1, j], lon2d[-1, j + 1], lon2d[-1, j + 1]  + dist_lonlat, lon2d[-1, j]  + dist_lonlat, lon2d[-1, j]]
        lats = [lat2d[-1, j], lat2d[-1, j + 1], lat2d[-1, j + 1] + dist_lat, lat2d[-1, j] + dist_lat, lat2d[-1, j]]

        
        lons = lons - abs(lon2d[-1, j] - lon2d[-1, j+1])/2           
        lats = lats - abs(lat2d[-1, j] - lat2d[-1, j])/2 - dist_lat/2  
        
        data_cell = ds['urmask'].values[-1, j]
    
        if data_cell == 1:
            ax.plot(lons, lats, color='red', zorder=100, linewidth=2)
        elif data_cell == 0:
            ax.plot(lons, lats, color='b', zorder=1, linewidth=2)
    
    # Plot the bottom right corner
    lons = [
        lon2d[-1, -1],
        lon2d[-1, -1] + dist_lon,
        lon2d[-1, -1] + dist_lon  + dist_lonlat,
        lon2d[-1, -1]  + dist_lonlat,
        lon2d[-1, -1]
    ]
    lats = [
        lat2d[-1, -1],
        lat2d[-1, -1] + dist_latlon,
        lat2d[-1, -1] + dist_lat + dist_latlon,
        lat2d[-1, -1] + dist_lat,
        lat2d[-1, -1]
    ]
        
    data_cell = ds['urmask'].values[-1, -1]

    lons = lons - abs(lon2d[-1, -1] - lon2d[-1, -1])/2 - dist_lon/2
    lats = lats - abs(lat2d[-1, -1] - lat2d[-1, -1])/2 - dist_lat/2  
    
    if data_cell == 1:
        ax.plot(lons, lats, color='red', zorder=100, linewidth=2)
    elif data_cell == 0:
        ax.plot(lons, lats, color='b', zorder=1, linewidth=2)


RCM_DICT = {
    'EUR-11': 
    { 
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-6',
    },
    'EUR-22': 
    {
        'REMO': 'GERICS_REMO2015',
    },
    'WAS-22': {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
    'EAS-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'KNU_RegCM4-0',
    },
    'CAM-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
    'SAM-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
    'NAM-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-6',
    },
    'AUS-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
    'AFR-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
    'SEA-22':
    {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
}