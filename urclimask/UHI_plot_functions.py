import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import os
import pandas as pd
import xarray as xr
from icecream import ic
from itertools import product
from urclimask.urban_areas import plot_urban_polygon


var_map = {
    'tasmin': 'TMIN',
    'tasmax': 'TMAX'
}

def plot_climatology(ds = None, 
                     variable = None, 
                     urban_vicinity = None, 
                     ucdb_city = [], 
                     valid_stations = [], 
                     time_series = [], 
                     city = None,
                     ax = None, 
                     alpha_urb_borders = 1, 
                     linewidth_urb_borders = 2, 
                     vmax = None
                    ):
    """
    Plot the climatological data.

    Parameters:
        ds (xr.Dataset): Dataset containing the climatological data.
        ucdb_city (gpd.GeoDataFrame): GeoDataFrame of the city boundaries.
        urban_vicinity (object): Object representing urban vicinity.
        obs (pd.DataFrame, optional): DataFrame containing observational data (default is None).

    Returns:
        matplotlib.figure.Figure: The generated figure.
    """
    proj = ccrs.PlateCarree()
    if ax == None: 
        fig, ax = plt.subplots(subplot_kw={'projection': proj}, figsize=(12, 6))
    else:
        fig = None
                         
    # calculate climatology and anomaly from the model
    ds_var_period_mean = ds.mean('time').compute()
    rural_mean = ds_var_period_mean[variable].where(
        urban_vicinity['urmask'] == 0).mean().compute()
    ds_anomaly = ds_var_period_mean[variable] - rural_mean
    if variable in ['huss', 'pr', 'sfcWind']:
        ds_anomaly = (ds_anomaly/ds_var_period_mean[variable])*100
    
    ds_anomaly.attrs['units'] = ds[variable].attrs.get('units', 'unknown')


    # Compute the maximum absolute value
    if vmax:
        max_abs_value = vmax
    else:
        max_abs_value = abs(ds_anomaly).max().item()
                         
    if not isinstance(valid_stations, list):
        # Filter stations that match the selected city
        city_stations = valid_stations[valid_stations['city'] == city]  # Filter by the selected city
        
        # Assign a single color to all selected stations
        station_color = 'g'  # Uniform color for all stations in the city
        
        # Plot all stations that match the selected city
        ax.scatter(city_stations.lon, 
                   city_stations.lat, 
                   c = station_color,  # Uniform color for all stations
                   marker='o', 
                   s = 40, 
                   edgecolors = 'gray', 
                   zorder = 10000)
    
    im1 = ax.pcolormesh(ds_anomaly.lon, ds_anomaly.lat, ds_anomaly.values,
                    cmap='bwr', alpha = 0.7,
                     vmin = -max_abs_value, vmax = max_abs_value)
    if fig:
        cbar = fig.colorbar(im1, ax = ax)
        # Dynamically get the unit from the dataset
        unit = ds_anomaly.attrs.get('units', 'unknown')  # Default to 'unknown' if 'units' is missing
        if variable in ['huss', 'pr', 'sfcWind']:
            cbar.set_label(f"% of {unit}", rotation = 90, fontsize = 14)
        else:
            cbar.set_label(f"{unit}", rotation = 90, fontsize = 14)
            
        if variable == 'tasmin':
            ax.set_title(f"Minimum temperature anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Minimum temperature anomaly ({unit})", fontsize=18)
        elif variable == 'tasmax':
            ax.set_title(f"Maximum temperature anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Maximum temperature anomaly ({unit})", fontsize=18)
        elif variable == 'huss':
            ax.set_title(f"Atmospheric moisture relative anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Atmospheric moisture anomaly (% relative to {unit})", fontsize=18)
        elif variable == 'hurs':
            ax.set_title(f"Relative atmospheric humidity anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Relative atmospheric humidity ({unit})", fontsize=18)
        elif variable == 'sfcWind':
            ax.set_title(f"Wind speed anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Wind speed anomaly (% relative to {unit})", fontsize=18)
        elif variable == 'pr':
            ax.set_title(f"Precipitation anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Precipitation anomaly (% relative to {unit})", fontsize=18)
    #if not isinstance(valid_stations, list):#change
        #ucdb_city.plot(ax=ax, facecolor="none", transform=proj, edgecolor="Green", linewidth=2, zorder = 1000)
    
    ax.coastlines()
    
    

    # Overlay the cell borders and handle NaNs
    #URBAN.plot_urban_borders(urban_vicinity, ax, 
    #                         alpha_urb_borders, 
    #                         linewidth_urb_borders)
    plot_urban_polygon(urban_vicinity, ax)
    
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    
    return fig

def plot_annual_cycle(variable = None, ds = None, urban_vicinity = None, 
                     time_series = [], valid_stations = [], 
                     data_squares = True, percentiles = [], 
                     var_map = var_map, ucdb_city = None, 
                     city = None, ax = None, 
                     vmax = None, vmin = None):
    '''
    Plot time series data with optional urban area overlay and additional time series overlay.

    Parameters:
    ds (xarray.Dataset): Dataset containing the variable data.
    variable (str): Name of the variable of interest.
    urban_vicinity (xarray.Dataset): Dataset containing information about urban areas.
    time_series (list of pandas.DataFrame, optional): List of time series dataframes to overlay on the plot.
    data_squares (bool, optional): Flag indicating whether to plot individual data squares for urban and rural areas.

    Returns:
    matplotlib.figure.Figure: The plotted figure.
    '''
    urban_area_legend = False
    not_urban_area_legend = False

    ds = ds[variable]
    
    is_rural = urban_vicinity['urmask'] == 0
    is_urban = urban_vicinity['urmask'] == 1
    rural_mean = (ds
        .where(is_rural)
        .groupby('time.month')
        .mean(dim = [ds.cf['Y'].name, ds.cf['X'].name, 'time'])
        .compute()
    )
                   
    ds_period_mean = ds.groupby('time.month').mean('time')                  
    ds_annomaly = ds_period_mean - rural_mean
    if variable in ['huss', 'pr', 'sfcWind']:
        ds_annomaly = (ds_annomaly / ds_period_mean) * 100

    rural_anomaly = ds_annomaly.where(is_rural)
    urban_anomaly = ds_annomaly.where(is_urban)

    urban_mean = urban_anomaly.mean(dim = [ds.cf['Y'].name, ds.cf['X'].name]).compute()   
        

    # Plot mean annual cycle (urban and rural)
    if ax == None: 
        fig, ax = plt.subplots(figsize=(15, 7))
    else:
        fig = None
        
    (urban_mean).plot(ax=ax,  color = '#A52A2A', linestyle='-', 
                                     linewidth = 4, label='Urban mean')
                         
    (rural_mean-rural_mean).plot(ax=ax,  color = '#8A8D28', linestyle='-', 
                                 linewidth = 4, label='Vicinity mean')
                         
    # Plot individual data squares for urban and rural areas if requested
    if data_squares:
        # Fill within percentiles
        axis = [rural_anomaly.get_axis_num(rural_anomaly.cf['X'].name),
                rural_anomaly.get_axis_num(rural_anomaly.cf['Y'].name)]
        colors = ['#8A8D28', '#A52A2A']
        for index, anom in enumerate([ rural_anomaly, urban_anomaly]):
            for percentile in percentiles:
                lower_percentile = np.nanpercentile(anom, percentile, axis=axis)
                upper_percentile = np.nanpercentile(anom, 100-percentile, axis=axis)
                ax.fill_between(
                    rural_anomaly['month'],
                    lower_percentile, upper_percentile,
                    color=colors[index], alpha=0.1, linewidth=1, linestyle = '--',
                )
    
                # Plot the lower percentile line
                ax.plot(
                    rural_anomaly['month'], lower_percentile,
                    color=colors[index], alpha=0.5, linewidth=1, linestyle='--', label=f'Lower Percentile')
                # Plot the upper percentile line
                ax.plot(
                    rural_anomaly['month'], upper_percentile,
                    color=colors[index], alpha=0.5, linewidth=1, linestyle='--', label=f'Upper Percentile')
            for i, j in product(anom.cf['X'].values, anom.cf['Y'].values):
                anom_val = anom.sel({ds.cf['X'].name:i,
                                     ds.cf['Y'].name:j})
                if not np.isnan(anom_val[0]):
                    anom_val.plot(ax=ax, color=colors[index], linewidth=0.1, alpha = 0.1)
                         
    #Plot the observation if requested
    if not isinstance(valid_stations, list):
        codes_ins_city = valid_stations.code[valid_stations['inside_city'] == True]
        codes_out_city = valid_stations.code[valid_stations['inside_city'] == False]
        time_series_mon = time_series.groupby(time_series.index.month).mean()
        time_series_mon_mean = pd.DataFrame(index = time_series_mon.index)
        time_series_mon_mean['rural_mean'] = time_series_mon[codes_out_city].mean(axis = 1).values
        time_series_mon_mean['urban_mean'] = time_series_mon[codes_ins_city].mean(axis = 1).values
        time_series_anomaly = time_series_mon.sub(time_series_mon_mean['rural_mean'], axis = 0)
        time_series_mon_mean_anom = time_series_mon_mean.sub(time_series_mon_mean['rural_mean'], axis = 0)
     
        time_series_anomaly[codes_ins_city].plot(ax = ax, marker='o', color = 'k', 
                                                 linestyle='--', linewidth = 2)
        time_series_anomaly[codes_out_city].plot(ax = ax, marker='o', color = 'g', 
                                                 linestyle='--', linewidth = 2)
    
        time_series_mon_mean_anom['urban_mean'].plot(ax = ax, color='k', linestyle='-', 
                                                     linewidth = 4, label='Urban obs. mean', 
                                                     zorder = 2000) 
        time_series_mon_mean_anom['rural_mean'].plot(ax = ax, color='g', linestyle='-', 
                                                     linewidth = 4, label='Vicinity obs. mean', 
                                                     zorder = 2000)
            
        
    ax.set_xticks(np.arange(1, 13))
    ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 
                        'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], fontsize = 18)
    ax.tick_params(axis='y', labelsize=18)
    if vmax is not None and vmin is not None:
        ax.set_ylim(vmin, vmax)

    else:
        unit = ds.attrs.get('units', 'unknown')  # Default to 'unknown' if 'units' is missing
    if fig:
        # Add legend to the plot
        ax.legend(fontsize = 14, loc='center left', bbox_to_anchor=(0, -0.2), prop={'size': 14})
        
        # Customize the plot
        #ax.set_xlabel('Month', fontsize = 18)
        # Dynamically get the unit from the dataset

    
        if variable == 'tasmin':
            ax.set_title(f"Minimum temperature anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Minimum temperature anomaly ({unit})", fontsize=18)
        elif variable == 'tasmax':
            ax.set_title(f"Maximum temperature anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Maximum temperature anomaly ({unit})", fontsize=18)
        elif variable == 'huss':
            ax.set_title(f"Atmospheric moisture relative anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Atmospheric moisture anomaly (% relative to {unit})", fontsize=18)
        elif variable == 'hurs':
            ax.set_title(f"Relative atmospheric humidity anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Relative atmospheric humidity ({unit})", fontsize=18)
        elif variable == 'sfcWind':
            ax.set_title(f"Wind speed anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Wind speed anomaly (% relative to {unit})", fontsize=18)
        elif variable == 'pr':
            ax.set_title(f"Precipitation anomaly for {city}", fontsize=18)
            ax.set_ylabel(f"Precipitation anomaly (% relative to {unit})", fontsize=18)
        

        return fig
    else:
        ax.set_title(f"", fontsize=18)
        ax.set_ylabel(f"", fontsize=18)
        
def plot_daily_cycle(variable, ds_var=None, urban_vicinity=None, 
             time_series=[], valid_stations=[], 
             data_squares=True, percentiles=[], 
             var_map=var_map, ucdb_city=None, 
             city=None, ax=None, vmax=None, 
             vmin=None, period= 'Annual', annomaly = True):
    '''
    Plot daily cycle data with optional urban area overlay and additional time series overlay.

    Parameters:
    ds_var (xarray.Dataset): Dataset containing the variable data.
    variable (str): Name of the variable of interest.
    urban_vicinity (xarray.Dataset): Dataset containing information about urban areas.
    time_series (list of pandas.DataFrame, optional): List of time series dataframes to overlay on the plot.
    data_squares (bool, optional): Flag indicating whether to plot individual data squares for urban and rural areas.

    Returns:
    matplotlib.figure.Figure: The plotted figure.
    '''
    urban_area_legend = False
    not_urban_area_legend = False
    
    if period == 'jja':
        ds_var = ds_var.sel(time=ds_var['time'].dt.month.isin([6,7,8]))
    if period == 'djf':
        ds_var = ds_var.sel(time=ds_var['time'].dt.month.isin([1,2,12]))
    is_rural = urban_vicinity['urmask'] == 0
    is_urban = urban_vicinity['urmask'] == 1
    rural_mean = (ds_var[variable]
        .where(is_rural)
        .groupby('time.hour')  # Group by daily hours
        .mean(dim=[ds_var.cf['Y'].name, ds_var.cf['X'].name, 'time'])
        .compute()
    )

    ds_var_period_mean = ds_var.groupby('time.hour').mean('time')
    # Group by daily hours
    if annomaly == True:
        ds_anomaly = ds_var_period_mean[variable] - rural_mean
        if variable in ['huss', 'pr', 'sfcWind']:
            ds_anomaly = (ds_anomaly / ds_var_period_mean[variable]) * 100
    else:
        ds_anomaly = ds_var_period_mean[variable]

    rural_anomaly = ds_anomaly.where(is_rural)
    urban_anomaly = ds_anomaly.where(is_urban)

    urban_mean = urban_anomaly.mean(dim=[ds_var.cf['Y'].name, ds_var.cf['X'].name]).compute()

    # Plot daily cycle (urban and rural)
    if ax is None: 
        fig, ax = plt.subplots(figsize=(15, 7))
    else:
        fig = None

    urban_mean.plot(ax=ax, color='#A52A2A', linestyle='-', 
                    linewidth=4, label='Urban mean')
    if annomaly == True:
        (rural_mean - rural_mean).plot(ax=ax, color='#8A8D28', linestyle='-', 
                                       linewidth=4, label='Vicinity mean')
    else:
        rural_mean.plot(ax=ax, color='#8A8D28', linestyle='-', 
                                       linewidth=4, label='Vicinity mean')
    # Plot individual data squares for urban and rural areas if requested
    if data_squares:
        axis = [rural_anomaly.get_axis_num(rural_anomaly.cf['X'].name),
                rural_anomaly.get_axis_num(rural_anomaly.cf['Y'].name)]
        colors = ['#8A8D28', '#A52A2A']
        for index, anom in enumerate([rural_anomaly, urban_anomaly]):
            for percentile in percentiles:
                lower_percentile = np.nanpercentile(anom, percentile, axis=axis)
                upper_percentile = np.nanpercentile(anom, 100 - percentile, axis=axis)
                ax.fill_between(
                    rural_anomaly['hour'],
                    lower_percentile, upper_percentile,
                    color=colors[index], alpha=0.1, linewidth=1, linestyle='--',
                )
    
                ax.plot(
                    rural_anomaly['hour'], lower_percentile,
                    color=colors[index], alpha=0.7, linewidth=1, linestyle='--', label=f'Lower Percentile')
    
                ax.plot(
                    rural_anomaly['hour'], upper_percentile,
                    color=colors[index], alpha=0.7, linewidth=1, linestyle='--', label=f'Upper Percentile')
            for i, j in product(anom.cf['X'].values, anom.cf['Y'].values):
                anom_val = anom.sel({ds_var.cf['X'].name: i, ds_var.cf['Y'].name: j})
                if not np.isnan(anom_val[0]):
                    anom_val.plot(ax=ax, color=colors[index], linewidth=0.1, alpha=0.1)

    # Plot observations if requested
    if not isinstance(valid_stations, list):
        codes_city = valid_stations.code[valid_stations['city'] == city]
        time_series_hour = time_series.groupby(time_series.index.hour).mean()
        time_series_hour_mean = pd.DataFrame(index=time_series_hour.index)
        time_series_hour_mean['urban_mean'] = time_series_hour[codes_city].mean(axis=1).values
        time_series_anomaly = time_series_hour.sub(time_series_hour_mean['urban_mean'], axis=0)

        time_series_anomaly[codes_city].plot(ax=ax, marker='o', color='k', 
                                             linestyle='--', linewidth=2)

        time_series_hour_mean['urban_mean'].plot(ax=ax, color='k', linestyle='-', 
                                                 linewidth=4, label='Urban obs. mean', 
                                                 zorder=2000)

    ax.set_xticks(np.arange(0, 24))
    ax.set_xticklabels([f'{h}' for h in range(24)], fontsize=12)
    ax.tick_params(axis='y', labelsize=18)
    if vmax is not None and vmin is not None:
        ax.set_ylim(vmin, vmax)
        
    unit = ds_var.attrs.get('units', 'unknown')

    if fig:
        ax.legend(fontsize=14, loc='center left', bbox_to_anchor=(0, -0.2), prop={'size': 14})

        if variable == 'tas':
            ax.set_title(f"Mean temperature daily cycle for {city} ({period})", fontsize=18)
            ax.set_ylabel(f"Mean temperature anomaly ({unit})", fontsize=18)
        elif variable == 'huss':
            ax.set_title(f"Atmospheric moisture daily cycle ({period})", fontsize=18)
            ax.set_ylabel(f"Atmospheric moisture anomaly (% relative to {unit})", fontsize=18)
        elif variable == 'hurs':
            ax.set_title(f"Relative atmospheric humidity daily cycle ({period})", fontsize=18)
            ax.set_ylabel(f"Relative atmospheric humidity ({unit})", fontsize=18)
        elif variable == 'sfcWind':
            ax.set_title(f"Wind speed daily cycle for ({period})", fontsize=18)
            ax.set_ylabel(f"Wind speed anomaly (% relative to {unit})", fontsize=18)
        elif variable == 'pr':
            ax.set_title(f"Precipitation daily cycle for ({period})", fontsize=18)
            ax.set_ylabel(f"Precipitation anomaly (% relative to {unit})", fontsize=18)

        return fig
    else:
        ax.set_title(f"", fontsize=18)
        ax.set_ylabel(f"", fontsize=18)
