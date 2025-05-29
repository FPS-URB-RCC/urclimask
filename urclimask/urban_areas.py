import cartopy.crs as ccrs
import dask
import pandas as pd
import geopandas as gpd
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr
from icecream import ic
from itertools import product
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.morphology import dilation, square, remove_small_objects
from  skimage import measure, morphology


from shapely.geometry import box
from shapely.geometry import Polygon

from urclimask.utils import plot_urban_polygon

class UrbanVicinity:
    def __init__(
        self,
        *,
        urban_th : float = 0.1,
        urban_sur_th : float = 0.1,
        orog_diff : float = 100.0,
        sftlf_th : float = 70,
        scale : float = 2.0, 
        min_city_size : int = 0, 
        lon_city : float | None = None,
        lat_city : float | None = None,
        lon_lim : float = 1.0,
        lat_lim : float = 1.0,
        model : str | None = None,
        domain : str | None = None,
        urban_var : str | None = None,
    ):
        """
        Hyperparameters requered for urban/rural area selection
    
        Parameters
        ----------
        urban_th : float
            Urban fracción threshold. Cells with urban fraccion values above this threshold are considered urban cells
        urban_sur_th : float
Urban surrounding threshold. Cells with urban fraction values below this threshold might be considered rural surrounding cells.
        orog_diff : float
Altitude difference (m) respects the maximum and minimum elevation of the urban cells.
        sftlf_th : float 
            Minimum fraction of land required to include a cell in the analysis
        scale : float 
            Ratio between rural surrounding  and urban grid boxes
        min_city_size : int
            Remove urban nuclei smaller than the specified size.
        lon_city : float
            Longitude of the city cente        
        lat_city : float
            Latitude of the city cente
        lon_lim : float
            Longitude limit of the study area respect to the city center. Cells outside this area are excluded from the analysis.
        lat_lim : float
            Latitude limit of the study area respect to the city center.
            Cells outside this area are excluded from the analysis.
        model : str
            GCM/RCM model name
        domain : str
            Mode domain (if applicable)
        urban_var : str
            Name of the urban static field (if applicable)
        """        
        self.urban_th = urban_th
        self.urban_sur_th = urban_sur_th
        self.orog_diff = orog_diff
        self.sftlf_th = sftlf_th
        self.scale = scale
        self.min_city_size = min_city_size
        self.lon_city = lon_city
        self.lat_city = lat_city
        self.lon_lim = lon_lim
        self.lat_lim = lat_lim
        self.model = model
        self.domain = domain
        self.urban_var = urban_var

    def crop_area_city(
        self, 
        *,
        ds : xr.DataArray | None = None,
        res : int | None = None
        ) -> xr.DataArray:
        """
        Select area around a central city point.
    
        Parameters
        ----------
        ds : xarray.DataArray 
            xarray with longitud and latitud.
        res : xarray.DataArray
            Domain resolution (e.g. 11/22).
            
        Returns
        -------
        ds : xarray.DataArray
            Cropped xarray.
        """
        # number of cells around the city
        dlon = int(111*self.lon_lim/res)
        dlat = int(111*self.lat_lim/res)
        # select point close the city
        dist = (ds['lon']-self.lon_city)**2 + (ds['lat']-self.lat_city)**2
        [ilat], [ilon] = np.where(dist == np.min(dist))
        if ds.lon.ndim == 2:
        # crop area
            ds = ds.isel(**{
            ds.cf['Y'].name: slice(ilat-dlat,ilat+dlat),
            ds.cf['X'].name : slice(ilon-dlon,ilon+dlon)
            })   
        else:
         		# Crop the area for the city using the domain resolution
            # Define trimming limits
            lat_min = self.lat_city - self.lat_lim
            lat_max = self.lat_city + self.lat_lim
            lon_min = self.lon_city - self.lon_lim
            lon_max = self.lon_city + self.lon_lim
            
            # Crop the dataset
            ds = ds.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
        return ds


    def remove_small_city(self,
                          *, 
                          mask: xr.DataArray) -> xr.DataArray:
        """
        Remove small urban regions from an input binary mask, retaining only the largest 
        or the region closest to the predefined city center if necessary.
    
        Parameters
        ----------
        mask : xarray.DataArray
            A binary mask (values 0 and 1) where 1 indicates urban areas.
        
        Returns
        -------
        xr.DataArray
            A cleaned binary mask where small objects have been removed, 
            preserving only the main urban region.
        """
        # Label connected regions in the mask
        labeled_mask = measure.label(mask.values)
    
        # Remove small objects based on a minimum city size threshold
        cleaned_mask = morphology.remove_small_objects(labeled_mask, min_size=self.min_city_size)
    
        # If all objects were removed, select the region closest to the city center
        if np.max(cleaned_mask) == 0:
            # Calculate pixel indices closest to the city center coordinates
            y_center = np.abs(mask.cf['lat'].values - self.lat_city).argmin()
            x_center = np.abs(mask.cf['lon'] - self.lon_city).argmin()
            center_pixel = np.array([y_center, x_center])
    
            # Compute distances from each region centroid to the city center
            distances = []
            for region in measure.regionprops(labeled_mask):
                region_center = np.array(region.centroid)  # (row, col) format
                distance = np.linalg.norm(region_center - center_pixel)
                distances.append((region.label, distance))
    
            # Select the label of the closest region
            closest_label = min(distances, key=lambda x: x[1])[0]
    
            # Create a mask keeping only the closest region
            cleaned_mask = (labeled_mask == closest_label)
        else:
            # If objects remain, set the cleaned_mask to 1 (urban) and 0 (non-urban)
            cleaned_mask = (cleaned_mask > 0)
    
        # Return the result as an xarray.DataArray with the original coordinates
        return xr.DataArray(
            cleaned_mask.astype(int),
            coords=mask.coords,
            dims=mask.dims)

    def define_masks(
        self, 
        *,
        ds_sfturf : xr.DataArray | None = None, 
        ds_orog : xr.DataArray | None = None, 
        ds_sftlf: xr.DataArray | None = None,
    )-> xr.DataArray:
        """
        Define masks for urban fraction, orography and land-sea mask.
    
        Parameters
        ----------
        ds_sfturf : xarray.DataArray 
            Urban fraction (0-1)
        ds_orog : xarray.DataArray
            Orogrhapy (m)
        ds_sftlf : xarray.DataArray
            Land-sea percentaje (%)
            
        Returns
        -------
        sfturf_mask : xarray.DataArray 
            Binary mask indicating urban areas with 1 and 0 for the rest.
        sfturf_sur_mask : xarray.DataArray 
            Binary mask indicating surroundings of urban areas affected by the urban effect with 1 and 0 for the rest.
        orog_mask : xarray.DataArray
            Binary mask indicating of orography with .
        sftlf_mask : xarray.DataArray
            Binary mask indicating sea areas.
        """
        # sfturf
        sfturf_mask = ds_sfturf[self.urban_var] > self.urban_th
        # Remove small objects
        sfturf_mask_rem_small = UrbanVicinity.remove_small_city(self,mask = sfturf_mask.astype(bool))
        sfturf_mask.data = sfturf_mask_rem_small
        deleted_small = ~sfturf_mask_rem_small*(ds_sfturf[self.urban_var] > self.urban_th)
        # Calculate surrounding mask and delete small objects from it
        sfturf_sur_mask_1 = ds_sfturf[self.urban_var] <= self.urban_th
        sfturf_sur_mask_2 = ds_sfturf[self.urban_var] > self.urban_sur_th
        sfturf_sur_mask_th = sfturf_sur_mask_1*sfturf_sur_mask_2
        sfturf_sur_mask = xr.where(deleted_small, True, sfturf_sur_mask_th)
        # orog
        urban_elev_max = ds_orog["orog"].where(sfturf_mask).max().item()
        urban_elev_min = ds_orog["orog"].where(sfturf_mask).min().item()
        orog_mask1 = ds_orog["orog"] < (self.orog_diff + urban_elev_max)
        orog_mask2 = ds_orog["orog"] > (urban_elev_min - self.orog_diff)
        orog_mask = orog_mask1 & orog_mask2
        
        #sftlf
        sftlf_mask = ds_sftlf["sftlf"] > self.sftlf_th   
        
        # Apply orog and sftlf thresholds to the urban_mask
        sfturf_mask = sfturf_mask*sftlf_mask
        sfturf_sur_mask = sfturf_sur_mask*sftlf_mask

        self.urban_elev_min = urban_elev_min 
        self.urban_elev_max = urban_elev_max 
    
        return sfturf_mask, sfturf_sur_mask, orog_mask, sftlf_mask


    def select_urban_vicinity(
        self,
        *,
        sfturf_mask : xr.DataArray | None = None,
        orog_mask : xr.DataArray | None = None,
        sftlf_mask : xr.DataArray | None = None,
        sfturf_sur_mask : xr.DataArray | None = None,
        scale: int | None = None
    ) -> xr.DataArray:
        """
        Funtion to select a number of non-urban cells based on surrounding urban areas using a dilation operation and excluding large water bodies, mountains and small urban nuclei.
    
        Parameters
        ----------
        sfturf_mask : xarray.DataArray 
            Binary mask indicating urban areas with 1 and 0 for the rest.
        orog_mask : xarray.DataArray
            Binary mask indicating of orography with .
        sftlf_mask : xarray.DataArray
            Binary mask indicating sea areas.
        sfturf_sur_mask : xarray.DataArray 
            Binary mask indicating surroundings of urban areas affected by the urban effect with 1 and 0 for the rest.
        scale : int 
            Urban-rural ratio of grid boxes.
    
        Returns
        -------
        xarray.DataArray 
            Mask of urban (1) surrounding non-urban cells (0) and the rest (NaN).
        """
        def delete_surrounding_intersect(dilated_data, sfturf_sur_mask):
            """
            Delete surroundings intersecting with dilated data
            """
            # Delete surroundings which intersect with dilated data
            dilated_data_surr = dilated_data * sfturf_sur_mask.astype(int)
            dilated_data_surr_opposite = xr.where(dilated_data_surr == 0, 1, 
                                                  xr.where(dilated_data_surr == 1, 0, dilated_data_surr))
            dilated_data = dilated_data * dilated_data_surr_opposite
            return dilated_data
        
        kernel1 = np.array([[0, 1, 0],
                            [1, 1, 1],
                            [0, 1, 0]])        
        kernel2 = np.array([[1, 1, 1],
                            [1, 1, 1],
                            [1, 1, 1]])
        
        if scale is None:
            scale = self.scale
        
        data_array = xr.DataArray(sfturf_mask).astype(int)

        urban_cells = np.sum(sfturf_mask).values
        non_urban_cells = 0
        counter = 0
        while non_urban_cells <= urban_cells * scale:
            # Dilation (Try with kernel 1)
            dilated_data = xr.apply_ufunc(dilation, 
                                          data_array if counter == 0 else dilated_data, 
                                          kwargs={'footprint': kernel1})
            # Delete fixed variables
            dilated_data = (dilated_data * orog_mask * sftlf_mask).astype(int)
            
            if np.sum(dilated_data) - urban_cells == non_urban_cells:
                #Try with kernel2
                dilated_data = xr.apply_ufunc(dilation, 
                                              data_array if counter == 0 else dilated_data, 
                                              kwargs={'footprint': kernel2})
                # Delete fixed variables
                dilated_data = (dilated_data * orog_mask * sftlf_mask).astype(int)
                
                if np.sum(dilated_data) - urban_cells  == non_urban_cells:
                    print(f"Warning: No more non-urban cells can be found in iteration number {counter}")
                    break
                    
            # Number of surrounding cells intersecting dilated data
            dilated_data_surr_cells = np.sum(dilated_data * sfturf_sur_mask.astype(int))
            non_urban_cells = (np.sum(dilated_data) - urban_cells).values - dilated_data_surr_cells
            counter += 1

        # Delete surrounding intersectig with dilated data
        dilated_data = delete_surrounding_intersect(dilated_data, sfturf_sur_mask)  
        # Assing rural cells (1), vicinity (0) and the rest (nan) 
        non_urban_mask = xr.DataArray(dilated_data.where(~sfturf_mask).fillna(0))
        urban_area = sfturf_mask.astype(int).where(sfturf_mask.astype(int) == 1, np.nan)
        urban_area = urban_area.where(non_urban_mask.astype(int) == 0, 0)
        urban_area = urban_area.to_dataset(name='urmask')
        # Add attributes
        urban_area = UrbanVicinity.netcdf_attrs(self, urban_area)
        
        return urban_area
            
    

    def netcdf_attrs(self, ds):        
        """
        Add metadata to urban area file.
    
        Parameters
        ----------
        urban_area : xarray.Dataset 
            Binary mask indicating urban areas (1), non-urban (vicinity) areas (0) and NaN for the rest.
        """
        # add attribtes
        ds['urmask'].attrs['long_name'] = 'Urban vs. vicinity. 1 corresponds to urban areas and 0 to the surrounding areas'
        
        attrs_list = ["urban_th", "urban_sur_th", "orog_diff", "sftlf_th", "sftlf_th", "scale", 
                      "min_city_size", "lon_city", "lat_city", "lon_lim", "lat_lim", "model", "domain"]
        
        for attr in attrs_list:
            if getattr(self, attr):
                ds['urmask'].attrs[attr] = getattr(self, attr)
            
        return ds
    def plot_static_variables(self,
                          *,
                          ds_sfturf, 
                          ds_orog, 
                          ds_sftlf,
                          sfturf_mask, 
                          orog_mask, 
                          sftlf_mask, 
                          urban_areas=None,
                          ucdb_city=None,
                          composite=False,
                          urban_th_plot=10,
                          ax=None):
        """
        Plot static variables: urban fraction, orography, and land-sea mask.
    
        Parameters
        ----------
        composite : bool, optional
            If True, plots all layers in a single composite axis. If False, creates 6 subplots.
        """
        colors = ['#278908', '#faf998', '#66473b']
        custom_cmap = LinearSegmentedColormap.from_list("custom_terrain", colors)
        proj = ccrs.PlateCarree()
    
        if composite:
            if ax is None:
                fig, ax = plt.subplots(subplot_kw={'projection': proj}, figsize=(16, 10))
            else:
                fig = ax.figure
    
            orog_data = ds_orog['orog'].where(orog_mask == 1, np.nan)
            im_orog = ax.pcolormesh(ds_orog.lon, ds_orog.lat, orog_data, cmap=custom_cmap,
                                    vmin=0, vmax=275, zorder=0, alpha=0.9)
            
            urban_data = ds_sfturf[self.urban_var].where(ds_sfturf[self.urban_var] > urban_th_plot, np.nan)
            im_urb = ax.pcolormesh(ds_sfturf.lon, ds_sfturf.lat, urban_data, cmap='binary',
                                   vmin=0, vmax=100, zorder=1)
    
            sftlf_inverse = np.where(sftlf_mask == 1, np.nan, ds_sftlf["sftlf"])
            im_landsea = ax.pcolormesh(ds_sftlf.lon, ds_sftlf.lat, sftlf_inverse, cmap='Blues_r',
                                       vmin=0, vmax=100, zorder=2, alpha=1)
    
            if ucdb_city is not None:
                ucdb_city.plot(ax=ax, facecolor="none", transform=proj, edgecolor="#ff66ff", linewidth=2, zorder=1000)
    
            if urban_areas:
                plot_urban_polygon(urban_areas, ax)
    
            ax.coastlines(zorder=1000)
            elev_lim_min = self.urban_elev_min - self.orog_diff
            elev_lim_max = self.urban_elev_max + self.orog_diff
            ax.set_title(
                f"Orography + Urban + Land-sea\n"
                f"Orography: {elev_lim_min:.0f}–{elev_lim_max:.0f} m | "
                f"Urban sfturf > {self.urban_th}, surr. ≤ {self.urban_sur_th} | "
                f"Land-sea ≤ {self.sftlf_th}%",
                fontsize=12
            )
    
            return fig
    
        else:
            fig, axes = plt.subplots(2, 3, subplot_kw={'projection': proj}, figsize=(20, 10))
    
            im1 = axes[0, 0].pcolormesh(ds_sfturf.lon, ds_sfturf.lat,
                                        ds_sfturf[self.urban_var],
                                        cmap='binary',
                                        vmin=np.nanmin(ds_sfturf[self.urban_var]),
                                        vmax=np.nanmax(ds_sfturf[self.urban_var]))
            fig.colorbar(im1, ax=axes[0, 0])
            axes[0, 0].set_title('Urban Fraction')
            axes[0, 0].coastlines()
    
            im2 = axes[0, 1].pcolormesh(ds_orog.lon, ds_orog.lat,
                                        ds_orog['orog'],
                                        cmap=custom_cmap,
                                        vmin=np.nanmin(ds_orog['orog']),
                                        vmax=np.nanmax(ds_orog['orog']))
            fig.colorbar(im2, ax=axes[0, 1])
            axes[0, 1].set_title('Orography')
            axes[0, 1].coastlines()
    
            im3 = axes[0, 2].pcolormesh(ds_sftlf.lon, ds_sftlf.lat,
                                        ds_sftlf["sftlf"],
                                        cmap='winter',
                                        vmin=np.nanmin(ds_sftlf["sftlf"]),
                                        vmax=np.nanmax(ds_sftlf["sftlf"]))
            fig.colorbar(im3, ax=axes[0, 2])
            axes[0, 2].set_title('Land-sea')
            axes[0, 2].coastlines()
    
            im4 = axes[1, 0].pcolormesh(ds_sfturf.lon, ds_sfturf.lat,
                                        ds_sfturf[self.urban_var].where(sfturf_mask == 1, np.nan),
                                        cmap='binary',
                                        vmin=np.nanmin(ds_sfturf[self.urban_var]),
                                        vmax=np.nanmax(ds_sfturf[self.urban_var]))
            fig.colorbar(im4, ax=axes[1, 0])
            axes[1, 0].set_title(f"Urban Fraction (urb_th > {self.urban_th})" if not urban_areas else
                                 f"Urban Fraction\n(urb_th = {self.urban_th}, urb_sur_th = {self.urban_sur_th}\nscale = {self.scale}, max_city = {self.min_city_size})")
            axes[1, 0].coastlines()
    
            im5 = axes[1, 1].pcolormesh(ds_orog.lon, ds_orog.lat,
                                        ds_orog['orog'].where(orog_mask == 1, np.nan),
                                        cmap=custom_cmap,
                                        vmin=np.nanmin(ds_orog['orog']),
                                        vmax=np.nanmax(ds_orog['orog']))
            fig.colorbar(im5, ax=axes[1, 1])
            elev_lim_min = self.urban_elev_min - self.orog_diff
            elev_lim_max = self.urban_elev_max + self.orog_diff
            axes[1, 1].set_title(f'Orography\n(orog_diff = {self.orog_diff};\n{elev_lim_min:.0f} m < orog < {elev_lim_max:.0f} m)')
            axes[1, 1].coastlines()
    
            im6 = axes[1, 2].pcolormesh(ds_sftlf.lon, ds_sftlf.lat,
                                        ds_sftlf["sftlf"].where(sftlf_mask == 1, np.nan),
                                        cmap='winter',
                                        vmin=np.nanmin(ds_sftlf["sftlf"]),
                                        vmax=np.nanmax(ds_sftlf["sftlf"]))
            fig.colorbar(im6, ax=axes[1, 2])
            axes[1, 2].set_title(f'Land-sea\n(sftlf > {self.sftlf_th}%)')
            axes[1, 2].coastlines()
    
            if urban_areas:
                for k in range(3):
                    plot_urban_polygon(urban_areas, axes[1, k])
    
            plt.subplots_adjust(wspace=0.1, hspace=0.1)
            return fig


    def create_urban_dataset(self,ucdb_city, ds):
        '''
        Creates a dataset for urban areas where each cell value represents the percentage of the cell area
        that is within the city limits.
    
        Parameters:
        ucdb_city: GeoDataFrame containing the city's geometry.
        ds: Dataset from which we are extracting the coordinates.
    
        Returns:
            ds_urban: An xarray dataset with cells representing the percentage of urban area coverage.
        '''
    
        # Create a grid of latitude and longitude edges to define cell boundaries
        lon_grid, lat_grid = np.meshgrid(ds['lon'].values, ds['lat'].values)
    
        # Combine all geometries into a single geometry (in case there are multiple city polygons)
        city_geometry = ucdb_city.geometry.unary_union
    
        # Create an array to store the percentage of each cell covered by the city
        urban_data = np.zeros((ds['lat'].size, ds['lon'].size))
    
        # Calculate percentage for each cell
        for i in range(ds['lat'].size):
            for j in range(ds['lon'].size):
                # Define the boundaries of the cell as a polygon
                cell_poly = box(
                    lon_grid[i, j] - 0.5 * np.diff(ds['lon'])[0],  # left
                    lat_grid[i, j] - 0.5 * np.diff(ds['lat'])[0],  # bottom
                    lon_grid[i, j] + 0.5 * np.diff(ds['lon'])[0],  # right
                    lat_grid[i, j] + 0.5 * np.diff(ds['lat'])[0]   # top
                )
    
                # Calculate the area of intersection with the city geometry
                intersection_area = cell_poly.intersection(city_geometry).area
                cell_area = cell_poly.area
    
                # Calculate the percentage of the cell covered by the city
                urban_data[i, j] = (intersection_area / cell_area) * 100
    
        # Create the final xarray dataset containing the urban percentage information
        ds_urban = xr.Dataset(
            {
                'sfturf': (['lat', 'lon'], urban_data)
            },
            coords={
                'lat': ds['lat'],
                'lon': ds['lon'],
            }
        )
        
        return ds_urban
