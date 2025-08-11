import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import os
import pandas as pd
import xarray as xr
from icecream import ic
from itertools import product
from urclimask.urban_areas import plot_urban_polygon

class UrbanIsland:
    def __init__(
        self,
        *,
        ds = None, 
        urban_vicinity = None, 
        rcm = None,
        anomaly = 'abs',
        period= 'Annual',
        obs_attributes = pd.DataFrame(), 
        obs_timeseries = pd.DataFrame(), 
        ):
        """
        Parameters:
        - ds (xarray.Dataset): Dataset containing the variable of interest.
        - urban_vicinity (xarray.DataArray or ndarray): Urban mask identifying urban grid cells.
        - anomaly (str): Type of anomaly to compute; either 'rel' (relative) or 'abs' (absolute).
        - period (str): Time period to filter ('Annual', 'jja', 'djf').
        - obs_attributes (list): Metadata for the observation points (e.g., locations).
        - obs_timeseries (list): Time series data for each observation point.
        """       
        self.ds = ds
        self.urban_vicinity = urban_vicinity
        self.rcm = rcm
        self.anomaly = anomaly
        self.period = period
        self.obs_attr = obs_attributes
        self.obs_time = obs_timeseries

    def compute_spatial_climatology(self):
        """
        Computes the spatial climatology and corresponding anomaly from both the model 
        dataset and the observations, if available.
    
        Outputs:
        - self.ds_spatial_climatology: Spatial anomaly computed from the model.
        - self.obs_spatial_climatology: Spatial anomaly computed from the observations (if available).
        """
        # calculate climatology and anomaly from the model
        ds_var_period_mean = self.ds.mean('time').compute()
        rural_mean = ds_var_period_mean.where(
            self.urban_vicinity['urmask'] == 0).mean().compute()
        ds_anomaly = ds_var_period_mean - rural_mean        
        ds_anomaly.attrs['units'] = self.ds.attrs.get('units', 'unknown')

        # calculate climatology from the observations
        obs_anomaly = None
        if not self.obs_attr.empty:
            obs_period_mean = self.obs_time.mean()
            obs_codes = self.obs_attr.loc[self.obs_attr["inside_city"] == False, "code"].astype(str)
            
            obs_rur_mean = obs_period_mean[obs_codes].mean( skipna = True)
            obs_anomaly = (obs_period_mean - obs_rur_mean)

        if self.anomaly == 'rel':
            ds_anomaly = (ds_anomaly/ds_var_period_mean)*100
            ds_anomaly.attrs['units'] = "%"
            if obs_anomaly is not None:
                obs_anomaly = (obs_anomaly/obs_period_mean)*100

        self.ds_spatial_climatology = ds_anomaly.compute()
        if obs_anomaly is not None:
            self.obs_spatial_climatology = obs_anomaly
        
        
    def compute_annual_cycle(self):
        """
        Computes the hourly climatology (annual cycle) and corresponding anomaly 
        from the model dataset and the observation, if available.

        Output:
        - self.ds_annual_cycle: Monthly anomaly for each month of the year from the model.
        - self.obs_annual_cycle: Monthly anomaly for each month of the year from the observation (if available).
        """
        is_rural = self.urban_vicinity['urmask'] == 0
        is_urban = self.urban_vicinity['urmask'] == 1
        rural_mean = (self.ds
            .where(is_rural)
            .groupby('time.month')
            .mean(dim = [self.ds.cf['Y'].name, self.ds.cf['X'].name, 'time'])
            .compute()
        )       
        ds_period_mean = self.ds.groupby('time.month').mean('time')                  
        ds_anomaly = (ds_period_mean - rural_mean).compute()
        ds_anomaly.attrs['units'] = self.ds.attrs.get('units', 'unknown')

        # calculate climatology from the observations
        obs_anomaly = None
        if not self.obs_attr.empty:
            codes_ins_city = self.obs_attr.code[self.obs_attr['inside_city'] == True]
            codes_out_city = self.obs_attr.code[self.obs_attr['inside_city'] == False]
            obs_month = self.obs_time.groupby(self.obs_time.index.month).mean()
            obs_month_mean = pd.DataFrame(index = obs_month.index)
            obs_month_mean['rural_mean'] = obs_month[codes_out_city].mean(axis = 1).values
            obs_month_mean['urban_mean'] = obs_month[codes_ins_city].mean(axis = 1).values
            obs_anomaly = obs_month_mean.sub(obs_month_mean['rural_mean'], axis = 0)
            raw_anomaly = obs_month.sub(obs_month_mean['rural_mean'], axis=0)
            for code in obs_month.columns:
                obs_anomaly[code] = raw_anomaly[code]

        if self.anomaly == 'rel':
            ds_anomaly = (ds_anomaly / ds_period_mean) * 100
            ds_anomaly.attrs['units'] = "%"
            if obs_anomaly is not None:
                obs_anomaly = obs_anomaly.div(obs_month_mean['rural_mean'], axis=0) * 100

        self.ds_annual_cycle = ds_anomaly
        self.obs_annual_cycle = obs_anomaly
    
    
    def compute_daily_cycle(self):
        """
        Computes the hourly climatology (daily cycle) and corresponding anomaly for the specified period
        from the model dataset and the observations, if available.
    
        Output:
        - self.ds_daily_cycle: Hourly anomaly for each hour of the day from the model.
        - self.obs_daily_cycle: Hourly anomaly for each hour of the day from the observations (if available).
        """
        ds_var = self.ds
        if self.period == 'jja':
            ds_var = ds_var.sel(time=ds_var['time'].dt.month.isin([6, 7, 8]))
        elif self.period == 'djf':
            ds_var = ds_var.sel(time=ds_var['time'].dt.month.isin([12, 1, 2]))
    
        is_rural = self.urban_vicinity['urmask'] == 0
        is_urban = self.urban_vicinity['urmask'] == 1
        rounded_time = ds_var.time.dt.round('h')
        ds_var = ds_var.assign_coords(hour=rounded_time.dt.hour)
        rural_mean = (
            ds_var
            .where(is_rural)
            .groupby('time.hour')
            .mean(dim=[ds_var.cf['Y'].name, ds_var.cf['X'].name, 'time'])
            .compute()
        )
    
        ds_period_mean = ds_var.groupby('time.hour').mean('time')
    
        ds_anomaly = (ds_period_mean - rural_mean).compute()
        ds_anomaly.attrs['units'] = ds_var.attrs.get('units', 'unknown')
    
        obs_anomaly = None
    
        if not self.obs_attr.empty:    
            codes_ins_city = self.obs_attr.code[self.obs_attr['inside_city'] == True]
            codes_out_city = self.obs_attr.code[self.obs_attr['inside_city'] == False]
    
            obs_hour = self.obs_time.groupby(self.obs_time.index.hour).mean()
    
            obs_hour_mean = pd.DataFrame(index=obs_hour.index)
            obs_hour_mean['rural_mean'] = obs_hour[codes_out_city].mean(axis=1).values
            obs_hour_mean['urban_mean'] = obs_hour[codes_ins_city].mean(axis=1).values
    
            obs_anomaly = obs_hour_mean.sub(obs_hour_mean['rural_mean'], axis=0)
            raw_anomaly = obs_hour.sub(obs_hour_mean['rural_mean'], axis=0)
            for code in obs_hour.columns:
                obs_anomaly[code] = raw_anomaly[code]
    
        if self.anomaly == 'rel':
            ds_anomaly = (ds_anomaly / ds_period_mean) * 100
            ds_anomaly.attrs['units'] = "%"
            if obs_anomaly is not None:
                obs_anomaly = obs_anomaly.div(obs_hour_mean['rural_mean'], axis=0) * 100
    
        self.ds_daily_cycle = ds_anomaly
        self.obs_daily_cycle = obs_anomaly

    
    def plot_UI_map(
        self,
        *,
        ds_anomaly = None,
        obs_anomaly = None,
        city_name = None,
        ucdb_city = None, 
        alpha_urb_borders = 1, 
        linewidth_urb_borders = 2, 
        vmax = None,
        ):
        """
        Plots the Urban Island effect using model and optional observation anomalies.
    
        Parameters:
        - ds_anomaly (xarray.DataArray, optional): Spatial anomaly from the model. If not provided, it will be computed.
        - obs_anomaly (numpy.ndarray, optional): Spatial anomaly from the observations. If not provided, it will be computed.
        - city_name (str, optional): Name of the city to include in the plot title.
        - ucdb_city: GeoDataFrame containing the city's geometry.
        - alpha_urb_borders (float): Transparency level of the urban borders (0 to 1).
        - linewidth_urb_borders (float): Line width of the urban borders.
        - vmax (float, optional): Maximum absolute value for the color scale. If not provided, it is computed from the data.
    
        Outputs:
        - fig (matplotlib.figure.Figure): The resulting figure object.
        """
        
        if ds_anomaly is None:
            if not hasattr(self, 'ds_spatial_climatology'):
                self.compute_spatial_climatology()
            ds_anomaly = self.ds_spatial_climatology
            
        if not self.obs_attr.empty:
            if obs_anomaly is None:
                if not hasattr(self, 'obs_spatial_climatology'):
                    self.compute_spatial_climatology()
                obs_anomaly = self.obs_spatial_climatology
        
        proj = ccrs.PlateCarree()
        fig, ax = plt.subplots(subplot_kw={'projection': proj}, figsize=(12, 6))
    
        # Compute the maximum absolute value
        if vmax:
            max_abs_value = vmax
        else:
            max_abs_value = abs(ds_anomaly).max().item()
        
        im1 = ax.pcolormesh(ds_anomaly.lon, ds_anomaly.lat, 
                            ds_anomaly.values,
                            cmap='bwr', alpha = 0.7,
                            vmin = -max_abs_value, vmax = max_abs_value)

        if obs_anomaly is not None:
            sc = ax.scatter(
                self.obs_attr.lon, 
                self.obs_attr.lat,
                c=obs_anomaly.values, 
                cmap='bwr',
                vmin=-max_abs_value, vmax=max_abs_value,
                marker='o', 
                s=40, 
                edgecolors='gray', 
                zorder=10000
            )

        if ucdb_city is not None:
            ucdb_city.plot(ax=ax, facecolor="none", transform=proj, edgecolor="#ff66ff", linewidth=2, zorder=1000)

        cbar = fig.colorbar(im1, ax = ax)
        unit = ds_anomaly.attrs.get('units', 'unknown') 
        
        cbar.set_label(f"{unit}", rotation = 90, fontsize = 14)
        
        if self.rcm:
            ax.set_title(f"Urban Island for {city_name} (variable: {ds_anomaly.name}, rcm: {self.rcm})", fontsize=18)
        else:
            ax.set_title(f"Urban Island for {city_name} (variable: {ds_anomaly.name})", fontsize=18)     

        ax.coastlines()
        plot_urban_polygon(self.urban_vicinity, ax)
        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        
        return fig

    def plot_UI_annual_cycle(
        self,
        *,
        ds_anomaly = None,
        obs_anomaly = None,
        percentiles = [5],
        gridcell_series = True, 
        city_name = None, 
        ax = None,            
        vmax = None, 
        vmin = None):
        '''
        Plots the annual cycle of the Urban Island effect using model and optional observation data.
        
        Parameters:
        - ds_anomaly (xarray.DataArray, optional): Monthly anomaly from the model. If not provided, it will be computed.
        - obs_anomaly (pandas.DataFrame, optional): Monthly anomaly from observations. If not provided, it will be computed.
        - percentiles (list of int): Percentiles to shade around the mean values (e.g., [5] for 5th–95th range).
        - gridcell_series (bool): Whether to plot individual grid cell time series as transparent lines.
        - city_name (str, optional): Name of the city to include in the plot title.
        - ax (matplotlib.axes.Axes, optional): Existing axis to plot on. If None, a new figure and axis are created.
        - vmax (float, optional): Upper limit for the y-axis.
        - vmin (float, optional): Lower limit for the y-axis.
        
        Outputs:
        - fig (matplotlib.figure.Figure): The resulting figure object.
        '''

        if ds_anomaly is None:
            if not hasattr(self, 'ds_annual_cycle'):
                self.compute_annual_cycle()
            ds_anomaly = self.ds_annual_cycle
        
        is_rural = self.urban_vicinity['urmask'] == 0
        is_urban = self.urban_vicinity['urmask'] == 1

        rural_anomaly = ds_anomaly.where(is_rural)
        urban_anomaly = ds_anomaly.where(is_urban)
    
        urban_mean = urban_anomaly.mean(dim = [ds_anomaly.cf['Y'].name, ds_anomaly.cf['X'].name]).compute()
        rural_mean = rural_anomaly.mean(dim = [ds_anomaly.cf['Y'].name, ds_anomaly.cf['X'].name]).compute()

        urban_area_legend = False
        not_urban_area_legend = False
    
        fig, ax = plt.subplots(figsize=(15, 7))

        (urban_mean).plot(ax=ax,  color = '#A52A2A', linestyle='-', 
                                         linewidth = 4, label='Urban mean')
                             
        (rural_mean-rural_mean).plot(ax=ax,  color = '#8A8D28', linestyle='-', 
                                     linewidth = 4, label='Vicinity mean')
                             
        if percentiles:
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
                        color=colors[index], alpha=0.5, linewidth=1, linestyle='--',  label=f'{percentile} to {100-percentile} Percentile')
                    # Plot the upper percentile line
                    ax.plot(
                        rural_anomaly['month'], upper_percentile,
                        color=colors[index], alpha=0.5, linewidth=1, linestyle='--')
                if gridcell_series:
                    for i, j in product(anom.cf['X'].values, anom.cf['Y'].values):
                        anom_val = anom.sel({ds_anomaly.cf['X'].name:i,
                                             ds_anomaly.cf['Y'].name:j})
                        if not np.isnan(anom_val[0]):
                            anom_val.plot(ax=ax, color=colors[index], linewidth=0.1, alpha = 0.1)
                             
        #Plot the observation if requested
        if not self.obs_attr.empty:
            if not hasattr(self, 'obs_annual_cycle'):
                self.compute_annual_cycle() 
            obs_anomaly = self.obs_annual_cycle
            codes_ins_city = self.obs_attr.code[self.obs_attr['inside_city'] == True]
            codes_out_city = self.obs_attr.code[self.obs_attr['inside_city'] == False]
            obs_anomaly[codes_ins_city].plot(ax = ax, marker='o', color = '#A52A2A', 
                                                     linestyle='--', linewidth = 1)
            obs_anomaly[codes_out_city].plot(ax = ax, marker='o', color = 'g', 
                                                     linestyle='--', linewidth = 1)
            obs_anomaly['urban_mean'].plot(ax = ax, color= 'k', linestyle='--', marker='o',
                                                         linewidth = 4, label='Urban mean (obs.)', 
                                                         zorder = 2000) 

            obs_anomaly['rural_mean'].plot(ax = ax, color='g', linestyle='--', marker='o',
                                                         linewidth = 4, label='Vicinity mean (obs.)', 
                                                         zorder = 2000)
                
            
        ax.set_xticks(np.arange(1, 13))
        ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 
                            'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], fontsize = 18)
        ax.tick_params(axis='y', labelsize=18)
        if vmax is not None and vmin is not None:
            ax.set_ylim(vmin, vmax)

        # Add legend to the plot
        #ax.legend(fontsize = 14, loc='center left', bbox_to_anchor=(0, -0.2), prop={'size': 14})

        if self.rcm:
            ax.set_title(f"Urban Island for {city_name} (variable: {ds_anomaly.name}, rcm: {self.rcm})", fontsize=18)
        else:
            ax.set_title(f"Urban Island for {city_name} (variable: {ds_anomaly.name})", fontsize=18)    

        unit = ds_anomaly.attrs.get('units', 'unknown')         
        ax.set_ylabel(f"{unit}")
      
        return fig

    def plot_daily_cycle(
        self,
        *,
        ds_anomaly = None,
        obs_anomaly = None,
        percentiles = [5],                      
        gridcell_series = True, 
        city_name = None, 
        ax = None,            
        vmax = None, 
        vmin = None):
        '''
        Plots the daily cycle of a variable, with optional anomalies, urban/rural differentiation, and observational overlays.
        
        Parameters:
        - ds_anomaly (xarray.DataArray, optional): Daily anomaly from the model. If not provided, it will be computed.
        - obs_anomaly (pandas.DataFrame, optional): Daily anomaly from observations. If not provided, it will be computed.
        - percentiles (list of int): Percentiles to shade around the mean values (e.g., [5] for 5th–95th range).
        - gridcell_series (bool): Whether to plot individual grid cell time series as transparent lines.
        - city_name (str, optional): Name of the city to include in the plot title.
        - ax (matplotlib.axes.Axes, optional): Existing axis to plot on. If None, a new figure and axis are created.
        - vmax (float, optional): Upper limit for the y-axis.
        - vmin (float, optional): Lower limit for the y-axis.
        
        Outputs:
        - fig (matplotlib.figure.Figure): The resulting figure object.
        '''
        if ds_anomaly is None:
            if not hasattr(self, 'ds_daily_cycle'):
                self.compute_daily_cycle()
            ds_anomaly = self.ds_daily_cycle
        
        is_rural = self.urban_vicinity['urmask'] == 0
        is_urban = self.urban_vicinity['urmask'] == 1

        rural_anomaly = ds_anomaly.where(is_rural)
        urban_anomaly = ds_anomaly.where(is_urban)
    
        urban_mean = urban_anomaly.mean(dim = [ds_anomaly.cf['Y'].name, ds_anomaly.cf['X'].name]).compute()
        rural_mean = rural_anomaly.mean(dim = [ds_anomaly.cf['Y'].name, ds_anomaly.cf['X'].name]).compute()

        urban_area_legend = False
        not_urban_area_legend = False
    
        fig, ax = plt.subplots(figsize=(15, 7))

        (urban_mean).plot(ax=ax,  color = '#A52A2A', linestyle='-', 
                                         linewidth = 4, label='Urban mean')
                             
        (rural_mean-rural_mean).plot(ax=ax,  color = '#8A8D28', linestyle='-', 
                                     linewidth = 4, label='Vicinity mean')

        
        # Plot individual data squares for urban and rural areas if requested
        if percentiles:
            # Fill within percentiles
            axis = [rural_anomaly.get_axis_num(rural_anomaly.cf['X'].name),
                    rural_anomaly.get_axis_num(rural_anomaly.cf['Y'].name)]
            colors = ['#8A8D28', '#A52A2A']
            for index, anom in enumerate([ rural_anomaly, urban_anomaly]):
                for percentile in percentiles:
                    lower_percentile = np.nanpercentile(anom, percentile, axis=axis)
                    upper_percentile = np.nanpercentile(anom, 100-percentile, axis=axis)
                    ax.fill_between(
                        rural_anomaly['hour'],
                        lower_percentile, upper_percentile,
                        color=colors[index], alpha=0.1, linewidth=1, linestyle = '--',
                    )
        
                    # Plot the lower percentile line
                    ax.plot(
                        rural_anomaly['hour'], lower_percentile,
                        color=colors[index], alpha=0.5, linewidth=1, linestyle='--', label=f'{percentile} to {100-percentile} Percentile')
                    # Plot the upper percentile line
                    ax.plot(
                        rural_anomaly['hour'], upper_percentile,
                        color=colors[index], alpha=0.5, linewidth=1, linestyle='--')
                for i, j in product(anom.cf['X'].values, anom.cf['Y'].values):
                    anom_val = anom.sel({ds_anomaly.cf['X'].name:i,
                                         ds_anomaly.cf['Y'].name:j})
                    if not np.isnan(anom_val[0]):
                        anom_val.plot(ax=ax, color=colors[index], linewidth=0.1, alpha = 0.1)

        # Plot observations if requested
        if not self.obs_attr.empty:
            if not hasattr(self, 'obs_daily_cycle'):
                self.compute_daily_cycle() 
            obs_anomaly = self.obs_daily_cycle
            codes_ins_city = self.obs_attr.code[self.obs_attr['inside_city'] == True]
            codes_out_city = self.obs_attr.code[self.obs_attr['inside_city'] == False]
            obs_anomaly[codes_ins_city].plot(ax = ax, marker='o', color = 'k', 
                                                     linestyle='--', linewidth = 2)
            obs_anomaly[codes_out_city].plot(ax = ax, marker='o', color = 'g', 
                                                     linestyle='--', linewidth = 2)
            obs_anomaly['urban_mean'].plot(ax = ax, color= 'k', linestyle='--', marker='o',
                                                         linewidth = 4, label='Urban mean (obs.)', 
                                                         zorder = 2000) 

            obs_anomaly['rural_mean'].plot(ax = ax, color='g', linestyle='--', marker='o',
                                                         linewidth = 4, label='Vicinity mean (obs.)', 
                                                         zorder = 2000)

        
        ax.set_xticks(np.arange(0, 24))
        ax.set_xticklabels([f'{h}' for h in range(24)], fontsize=12)
        ax.tick_params(axis='y', labelsize=18)
        if vmax is not None and vmin is not None:
            ax.set_ylim(vmin, vmax)

        # Add legend to the plot
        ax.legend(fontsize = 14, loc='center left', bbox_to_anchor=(0, -0.2), prop={'size': 14})

        ax.set_title(f"Urban Island for {city_name} (variable: {ds_anomaly.name})", fontsize=18)

        unit = ds_anomaly.attrs.get('units', 'unknown')         
        ax.set_ylabel(f"{unit}")
      
        return fig