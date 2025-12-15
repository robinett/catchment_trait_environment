import random
import netCDF4 as nc
import numpy as np
import sys
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import os
import copy

class choose_unreg:
    def __init__(self):
        random.seed(123)
        np.random.seed(190)
    def get_basin_info(self,fname):
        s_info = nc.Dataset(fname)
        print(s_info)
        sys.exit()
        return s_info
    def get_tile_info(self,fname):
        tile_info = pd.read_csv(fname)
        tile_info = tile_info.set_index('tile_id')
        return tile_info
    def get_boundary_info(self,fname):
        boundary_info = gpd.read_file(fname)
        boundary_info = boundary_info.to_crs('EPSG:4326')
        return boundary_info
    def get_station_names(self,s_info):
        all_names = np.array(s_info['Station_ID']).astype(np.float64)
        station_names  = np.zeros(np.shape(all_names)[0])
        for s,stat in enumerate(all_names):
            stat_str = ''
            for num in stat:
                stat_str = stat_str + str(int(num))
            station_names[s] = int(stat_str)
        return station_names
    def select_random_stations(self,station_names,num_stations):
        selected = np.random.choice(station_names,num_stations,replace=False)
        print(selected)
        return selected
    def confirm_nc(self,station_names,s_info,tile_info,boundary_info,
                   states_shp_fname,plots_dir):
        # print our variable information
        #print(s_info.variables)
        # intialize our dictionaries and arrays that will hold data
        intersection_gdfs = {}
        intersecting_tiles = np.zeros(0)
        total_tiles = 0
        # get the data that we need from the nc4
        all_station_tiles = np.array(s_info['SMAPID'])
        all_station_tiles_perc = np.array(s_info['SMAPFrac'])
        all_station_num_tiles = np.array(s_info['N_cells_basin'])
        # get where we have station information and boundary information
        bound_sites = np.array(boundary_info['BasinID'])
        bound_sites_int = np.zeros(len(bound_sites),dtype='int')
        for b,bo in enumerate(bound_sites):
            try:
                b_int = int(bo)
            except:
                b_int = 0
            bound_sites_int[b] = b_int
        station_with_bound = np.zeros(0)
        for s in station_names:
            if s in bound_sites_int:
                station_with_bound = np.append(station_with_bound,s)
        print(station_with_bound)
        print(boundary_info)
        # get the information that we need for each of the selected basins
        # for each huc in the subset
        for s,sel in enumerate(station_with_bound):
            plt.figure()
            # plot the united states outline
            # get the shapefile
            states = gpd.read_file(states_shp_fname)
            states = states.to_crs('EPSG:4326')
            # get rid of non-conus states since not considering
            non_conus = ['HI','VI','MP','GU','AK','AS','PR']
            states_conus = states
            for n in non_conus:
                states_conus = states_conus[states_conus.STUSPS != n]
            # add this to the plot
            states_conus.boundary.plot()
            station_idx = np.where(station_names == sel)
            this_station_tiles = all_station_tiles[station_idx,:]
            tiles_idx = np.where(this_station_tiles != -2147483647)
            this_station_tiles = this_station_tiles[tiles_idx]
            this_station_tiles_perc = all_station_tiles_perc[station_idx,:]
            perc_idx = np.where(this_station_tiles_perc != 9.9692100e+36)
            this_station_tiles_perc = this_station_tiles_perc[perc_idx]
            this_station_num_tiles = all_station_num_tiles[station_idx]
            print(this_station_tiles)
            print(this_station_tiles_perc)
            print(this_station_num_tiles)
            tiles_unique = np.zeros(0)
            tiles_perc = np.zeros(0)
            for t,ti in enumerate(this_station_tiles):
                this_tile_idx = np.where(this_station_tiles == ti)
                this_tile_perc = this_station_tiles_perc[this_tile_idx]
                this_tot_perc = np.sum(this_tile_perc)
                if ti not in tiles_unique:
                    tiles_unique = np.append(tiles_unique,ti)
                    tiles_perc = np.append(tiles_perc,this_tot_perc)
            print(tiles_unique)
            print(tiles_perc)
            print(tile_info)
            # get the polygon for the basin
            sel_str = str(sel)
            sel_str = sel_str[:-2]
            sel_str = sel_str.zfill(8)
            this_basin_idx = np.where(boundary_info['BasinID'] == str(sel_str))
            this_basin_poly = list(
                boundary_info['geometry'].iloc[this_basin_idx]
            )
            this_basin_poly = this_basin_poly[0]
            # create the polyon for the tiles of interest
            polygons = []
            min_lon = np.inf
            min_lat = np.inf
            max_lon = -np.inf
            max_lat = -np.inf
            for t,ti in enumerate(tiles_unique):
                this_lon_vals = [
                    tile_info['min_lon'].loc[ti],
                    tile_info['max_lon'].loc[ti],
                    tile_info['max_lon'].loc[ti],
                    tile_info['min_lon'].loc[ti]
                ]
                this_lat_vals = [
                    tile_info['max_lat'].loc[ti],
                    tile_info['max_lat'].loc[ti],
                    tile_info['min_lat'].loc[ti],
                    tile_info['min_lat'].loc[ti]
                ]
                this_min_lon = np.min(this_lon_vals)
                this_max_lon = np.max(this_lon_vals)
                this_min_lat = np.min(this_lat_vals)
                this_max_lat = np.max(this_lat_vals)
                if this_min_lon < min_lon:
                    min_lon = copy.deepcopy(this_min_lon)
                if this_min_lat < min_lat:
                    min_lat = copy.deepcopy(this_min_lat)
                if this_max_lon > max_lon:
                    max_lon = copy.deepcopy(this_max_lon)
                if this_max_lat > max_lat:
                    max_lat = copy.deepcopy(this_max_lat)
                # make this lat/lon into a polygon
                this_polygon = Polygon(
                    zip(this_lon_vals,this_lat_vals)
                )
                # keep track for safekeeping
                polygons.append(this_polygon)
            min_lon = min_lon - .5
            min_lat = min_lat - .5
            max_lon = max_lon + .5
            max_lat = max_lat + .5
            df = pd.DataFrame({
                'tile':tiles_unique,
            })
            tile_gdf = gpd.GeoDataFrame(
                df,geometry=polygons,crs='EPSG:4326'
            )
            tile_gdf = tile_gdf.set_index('tile')
            # plot the polygon, tiles of interest with the percents in the
            # tiles to judge whether this information is accurate
            for t,ti in enumerate(tiles_unique):
                this_geom = tile_gdf['geometry'].loc[ti]
                print(this_geom)
                plt.plot(*this_geom.exterior.xy)
                # get the centroid and plot the percentager here
                centroid = this_geom.centroid
                print(centroid)
                print(tiles_perc[t])
                plt.text(
                    centroid.x,
                    centroid.y,
                    '{:.2f}'.format(
                        tiles_perc[t]
                    ),
                    fontsize=2
                )
            plt.plot(*this_basin_poly.exterior.xy)
            save_name = os.path.join(
                plots_dir,
                'basin_and_tiles_9km_{}.png'.format(
                    sel
                )
            )
            plt.ylim(bottom=min_lat,top=max_lat)
            plt.xlim(left=min_lon,right=max_lon)
            plt.savefig(save_name,dpi=350,bbox_inches='tight')
            plt.close()
