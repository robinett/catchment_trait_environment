import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
import sys
import random
import pandas as pd
from shapely.geometry import Polygon
import pickle
import os
import pandas as pd
import copy

class choose:
    def __init__(self):
        # set the random seed for repredocible results
        np.random.seed(123)
    def get_catch_tiles(self,tile_fname,camels_boundaries_fname,max_tiles):
        # get the catchment tile info
        tile_info = pd.read_csv(tile_fname)
        tile_info = tile_info.set_index('tile_id')
        tile_idx = tile_info.index
        print(tile_info)
        # get the camels boundaries
        camels_bounds = gpd.read_file(camels_boundaries_fname)
        camels_bounds = camels_bounds.set_index('hru_id')
        camels_sheds = np.array(camels_bounds.index)
        camels_sheds_random = copy.deepcopy(camels_sheds)
        np.random.shuffle(camels_sheds_random)
        # create a gpd geodataframe so can check for intersection
        polygons = []
        # set up dictionary that will store intersection informat
        # each key in the dictionary is the smaller, subset huc (huc6 or huc8
        # dependeing on choice
        # the first list in list of lists is intersecting catchment pixels for
        # that huc
        # the second list in the list of lists is the percent of that catchment
        # pixel within that huc
        # the third list in the list of lists is the area (in XXXX) of the
        # initersection between Catchment pixel and huc
        intersection_info = {}
        # set up numpy array that will keep track of all tiles, just so we know
        # how many total we are running
        all_intersecting = np.zeros(0)
        # loop over all tiles
        for t,ti in enumerate(tile_idx):
            # get lat and lon for the four points that will define this polygon
            # top left, top right, bottom right, bottom left
            if ti%1000 == 0:
                print('creating polygon for tile {}'.format(ti))
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
            # make this lat/lon into a polygon
            this_polygon = Polygon(
                zip(this_lon_vals,this_lat_vals)
            )
            # keep track for safekeeping
            polygons.append(this_polygon)
        # create a geopandas geodf from these polygons
        df = pd.DataFrame({
            'tile':tile_idx
        })
        tile_gdf = gpd.GeoDataFrame(
            df,geometry=polygons,crs='EPSG:4326'
        )
        tile_gdf = tile_gdf.set_index('tile')
        print(tile_gdf)
        # convert to coordinates that  can be used for m2
        tile_gdf_m = tile_gdf.to_crs({'proj':'cea'})
        print(tile_gdf_m)
        # we're going to want to save these areas for later use
        camels_areas = np.zeros(0)
        chosen_camels = []
        for h,huc in enumerate(camels_sheds_random):
            print('check interesection for huc {}'.format(huc))
            # create intersection info for this
            if h < len(camels_sheds_random) - 1:
                this_camel_gdf = camels_bounds.loc[[
                    huc,camels_sheds_random[h+1]
                ]]
            else:
                this_camel_gdf = camels_bounds.loc[[
                    huc,camels_sheds_random[0]
                ]]
            this_camel_gdf_m = this_camel_gdf.to_crs({'proj':'cea'})
            this_camel_m = this_camel_gdf_m['geometry'].loc[huc]
            this_camel_area_km = this_camel_m.area/(10**6)
            intersects = list(tile_gdf_m.intersects(this_camel_m))
            for i,it in enumerate(intersects):
                if it and this_camel_area_km > 1000:
                    # get the polygon for the intersecting tile
                    this_tile = tile_idx[i]
                    this_tile_m = tile_gdf_m['geometry'].loc[tile_idx[i]]
                    # get the overlapping area in km**2
                    tile_area_km = this_tile_m.area/(10**6)
                    intersection_area_km = (
                        this_camel_m.intersection(this_tile_m).area/(10**6)
                    )
                    perc_in_camel = intersection_area_km/tile_area_km
                    if huc not in intersection_info.keys():
                        intersection_info[huc] = [[],[],[]]
                        chosen_camels.append(huc)
                        camels_areas = np.append(
                            camels_areas,this_camel_area_km
                        )
                    intersection_info[huc][0].append(this_tile)
                    intersection_info[huc][1].append(perc_in_camel)
                    intersection_info[huc][2].append(intersection_area_km)
                    if this_tile not in all_intersecting:
                        all_intersecting = np.append(
                            all_intersecting,this_tile
                        )
            print('size of this camel (km): {}'.format(this_camel_area_km))
            print(
                'number of camels iterated so far: {}'.format(
                    h+1
                )
            )
            print(
                'number of camels selected so far: {}'.format(
                    len(intersection_info)
                )
            )
            print(
                'number of tiles we need to run so far: {}'.format(
                    len(all_intersecting)
                )
            )
        print(intersection_info)
        # sort the list of catchment pixels in ascending order
        all_intersecting = np.sort(all_intersecting)
        really_chosen_tiles_shapes = tile_gdf.loc[all_intersecting]
        # make sure this list in unique--don't run a tile twice if it
        # intersects multiple watershed
        # should have taken care of this in the newest version of the code but
        # we will keep it just to be sure
        all_intersecting = np.unique(all_intersecting)
        # get the list of camels and sort in ascenting order
        camels_final = pd.DataFrame(columns=['area'],index=chosen_camels)
        camels_final['area'] = camels_areas
        camels_final.index.name = 'camel'
        return [
            intersection_info,all_intersecting,tile_gdf,
            camels_final,camels_bounds,
            really_chosen_tiles_shapes
        ]
    def save_dictionary(self,dictionary,save_name):
        print('saving dictionary to {}'.format(save_name))
        with open(save_name,'wb') as f:
            pickle.dump(dictionary,f)
    def save_array_as_csv(self,to_save,save_name):
        as_df = pd.DataFrame(to_save)
        as_df.to_csv(save_name,index=False,header=False)
    def save_df_as_csv(self,to_save,save_name):
        to_save.to_csv(save_name)
    def create_include_file(self,all_intersecting,include_fname):
        '''
        create a file in the style of the 'include' files needed for running
        Catchment-CN only at specified tiles
        '''
        # remove the file if it already exists
        if os.path.isfile(include_fname):
            os.remove(include_fname)
        # loop over all tiles to run
        for t,ti in enumerate(all_intersecting):
            # turn our itle number into a string
            ti_str = str(int(ti))
            # format it into a len(12) string with empty characters leading
            empty_str = ' '*12
            len_ti = len(ti_str)
            first_part = empty_str[:-len_ti]
            formatted_str = first_part + ti_str
            # write this to the include file
            with open(include_fname,'a') as f:
                if t != (len(all_intersecting) - 1):
                    f.write(formatted_str + '\n')
                else:
                    f.write(formatted_str)
    def save_gdf(self,gdf,fname):
        gdf.to_file(fname)
    def load_dictionary(self,fname):
        with open(fname,'rb') as f:
            d = pickle.load(f)
        return d
    def load_gdf(self,fname,index_name='False'):
        gdf = gpd.read_file(fname)
        if index_name != 'False':
            gdf = gdf.set_index(index_name)
        return gdf
    def plot_hucs_tiles(self,intersection_info,tiles_shp,camels_shp,
                        states,plots_dir):
        # plot the united states outline
        # get the shapefile
        #print('converting states shp')
        #states = states.to_crs('EPSG:4326')
        # get rid of non-conus states since not considering
        non_conus = ['HI','VI','MP','GU','AK','AS','PR']
        states_conus = states
        for n in non_conus:
            states_conus = states_conus[states_conus.STUSPS != n]
        # open a figure
        plt.figure()
        # add this to the plot
        print('plotting states shp')
        states_conus.boundary.plot()
        # for each huc in the subset
        tiles_plotted = np.zeros(0)
        for key in intersection_info.keys():
            # for each tile that intersects this huc
            for ti in intersection_info[key][0]:
                if ti not in tiles_plotted:
                    # get the tile geometry
                    this_tile = tiles_shp['geometry'].loc[ti]
                    # plot the tile geometry
                    plt.plot(*this_tile.exterior.xy)
                    tiles_plotted = np.append(
                        tiles_plotted,ti
                    )
            # get geometry of this huc
            #print(camels_shp.index)
            this_geom = camels_shp['geometry'].loc[key]
            if this_geom.geom_type == 'Polygon':
                # plot this geometry
                plt.plot(*this_geom.exterior.xy)
            elif this_geom.geom_type == 'MultiPolygon':
                for this_this_geom in this_geom.geoms:
                    plt.plot(*this_this_geom.exterior.xy)
            else:
                raise IOError('Shape is not a polygon.')
            # get the centroid and plot the huc number there
            #centroid = this_geom.centroid
            #print(centroid)
            #print(idx)
            #plt.text(centroid.x,centroid.y,idx)
        # save and close
        save_name = os.path.join(
            plots_dir,
            'selected_camels_and_tiles_conus.png'
        )
        plt.savefig(save_name,dpi=350,bbox_inches='tight')
        plt.close()
        # now let's zoom in to make a plot for every watershed to check our
        # method
        for key in intersection_info.keys():
            print('plotting zoomed in for watershed {}'.format(key))
            plt.figure()
            for t,ti in enumerate(intersection_info[key][0]):
                this_tile = tiles_shp['geometry'].loc[ti]
                this_tile_perc = intersection_info[key][1][t]
                # plot the tile geometry
                plt.plot(*this_tile.exterior.xy)
                centroid = this_tile.centroid
                this_text = '{:.2f}'.format(this_tile_perc)
                plt.text(centroid.x,centroid.y,this_text)
            # get geometry of this huc
            this_geom = camels_shp['geometry'].loc[key]
            if this_geom.geom_type == 'Polygon':
                # plot this geometry
                plt.plot(*this_geom.exterior.xy)
            elif this_geom.geom_type == 'MultiPolygon':
                for this_this_geom in this_geom.geoms:
                    plt.plot(*this_this_geom.exterior.xy)
            else:
                raise IOError('Shape is not a polygon.')
            save_name = os.path.join(
                plots_dir,
                'camel_and_tile_{}.png'.format(
                    key
                )
            )
            plt.savefig(save_name,dpi=350,bbox_inches='tight')
            plt.close()




