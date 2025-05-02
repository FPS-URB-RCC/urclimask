import pandas as pd
import geopandas as gpd
import xarray as xr
import numpy as np
from shapely.geometry import Point

var_map = {
    'tasmin': 'TMIN',
    'tasmax': 'TMAX'
}

def load_ghcnd_stations(lon, lat, radious = 0.5):
    '''
    Load GHCND stations near a specific location.

    Parameters:
    lon (float): Longitude of the selected city.
    lat (float): Latitude of the selected city.
    radius (float): Maximum distance allowed.

    Returns:
    gpd.GeoDataFrame: Geospatial DataFrame of nearby GHCND stations.
    '''
    ghcnd_stations_url = 'https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/doc/ghcnd-stations.txt'
    ghcnd_stations_url = '/lustre/gmeteo/WORK/diezsj/research/cordex-fps-urb-rcc/old/CORDEX-CORE_WG_delete/ghcnd-stations.txt'
    ghcnd_stations_column_names = ['code', 'lat', 'lon', 'elev', 'name', 'net', 'numcode']
    ghcnd_stations_column_widths = [   11,     9,    10,      7,     34,     4,       10 ]
    df = pd.read_fwf(ghcnd_stations_url, header = 0, widths = ghcnd_stations_column_widths, names = ghcnd_stations_column_names)
    ghcnd_stations=gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs = 'EPSG:4326')
    rval = ghcnd_stations.assign(dist = ghcnd_stations.distance(Point(lon, lat)))
    rval.sort_values(by = 'dist', inplace = True)
    rval = rval[rval.dist < radious].to_crs(epsg=3857)
    return rval

def get_ghcnd_df(code):
    '''
    Load GHCND data for a specific station.

    Parameters:
    code (str): The station code.

    Returns:
    pd.DataFrame: DataFrame containing the GHCND data for the specified station.
    '''
    baseurl = '/lustre/gmeteo/WORK/WWW/chus/ghcnd/data'
    try:
        # Attempt to load the compressed file from the original location
        rval = pd.read_csv(f'{baseurl}/{code[0]}/{code}.csv.gz',
                           compression='gzip',
                           index_col='DATE',
                           parse_dates=True,
                           low_memory=False  # Avoid warnings for mixed data types in some columns
                           )
    except:
       try:
            # Load the combined_temp_data.csv file
            file_path = 'PARIS_surface_weather_data/combined_temp_data.csv'
            combined_data = pd.read_csv(file_path)

            # Filter the data by the station code
            rval = combined_data[combined_data['code'] == int(code)]
            rval = rval.set_index('DATE')

            # Raise an error if no data is found for the station code
       except Exception as e:
            # Handle any errors during the second attempt
            print(f"Error loading data for {code}: {e}")
            rval = pd.DataFrame()  # Return an empty DataFrame if everything fails

    return rval

def get_valid_timeseries(city, stations, ds_var, variable = 'tasmin', valid_threshold=0.8, idate='1979-01-01', fdate='2014-12-31',var_map=var_map, divide=10.0):
    '''
    Retrieves valid time series data for a specific variable from GHCND stations for a given city.

    Parameters:
    city (str): The name of the city for which the data is to be retrieved.
    stations (GeoDataFrame): A GeoDataFrame containing station metadata.
    var (str): The variable of interest (default is 'PRCP' for precipitation).
    valid_threshold (float): The threshold proportion of valid records required (default is 0.8).
    idate (str): The start date for the period of interest (default is '1979-01-01').
    fdate (str): The end date for the period of interest (default is '2014-12-31').
    var_map (dict): A dictionary mapping the variable names from the input to the dataset variable names.
    
    Returns:
    tuple: A tuple containing:
        - GeoDataFrame: The subset of stations with valid data.
        - pd.DataFrame: A DataFrame of valid time series data.
        - xr.Dataset: The subset of the dataset containing the selected period.
    '''
    var = var_map.get(variable, None)
    period = slice(idate, fdate)
    ds_var_period=ds_var.sel(time=period)
    ndays = (pd.to_datetime(fdate)-pd.to_datetime(idate)).days
    valid_codes, valid_time_series = [], []
    for stn_code in stations.code:
        stn_data = get_ghcnd_df(stn_code)
        if stn_data.empty:
            continue
        availvars = available_vars(stn_data)
        if var in availvars:
            valid_records = stn_data[var].loc[period].notna().sum()/ndays
            if valid_records > valid_threshold:
                print(f'{city} -- {stn_data.NAME[0]} - {var} has {100*valid_records:.1f}% valid records in {idate} to {fdate}')
                valid_codes.append(stn_code)
                valid_time_series.append({'data':stn_data[var].loc[period]/divide,'code':stn_code})

    #convert list in a dataframe
    if valid_time_series:
        for n_s, serie in enumerate(valid_time_series):
            if n_s == 0:
                freq = xr.infer_freq(serie['data'].index)
                df_time_series_obs = pd.DataFrame(
                        index = pd.date_range(period.start, period.stop, 
                        freq = freq)
                    )
            serie['data'].index = pd.to_datetime(serie['data'].index)

            serie['data'] = serie['data'].reindex(df_time_series_obs.index)
            df_time_series_obs[serie['code']] = serie['data']

    else:
        df_time_series_obs = valid_time_series
        
    return(stations[stations.code.isin(valid_codes)], df_time_series_obs, ds_var_period)

def available_vars(station):
    """
    Determines which variables are available in the station's dataset.

    Parameters:
    station (DataFrame): The DataFrame containing the station's data.

    Returns:
    set: A set of available variables that intersect with the known set of variables.
    """
    return(set(station.columns).intersection({'PRCP', 'TAVG', 'TMAX', 'TMIN', 'SNWD'}))


def inside_city(valid_obs, ucdb_city):
    """
    Add a column to the dataframe with the atributes of the series 
    including in they are inside or outside the city.
    """

    valid_obs['inside_city'] = np.nan
    for index, obs in valid_obs.iterrows():
        # Create a point with latitude and longuitude
        point = Point(obs['lon'], obs['lat'])
        is_inside = ucdb_city.contains(point)
        # Clasify in urban and vicinity
        valid_obs.loc[valid_obs['code'].str.contains(obs.code), 'inside_city'] = is_inside.values[0]

    n_series_inside = (valid_obs['inside_city'].values == True).sum()
    n_series_ouside = (valid_obs['inside_city'].values == False).sum()
    
    print(f"There are {n_series_inside} series inside the city and {n_series_ouside} outside")

    if (n_series_inside == 0) or (n_series_ouside == 0):
        print(f"The number of inside/outside observations is 0 therefore valid_obs is empty")
        valid_obs = []

    return valid_obs


