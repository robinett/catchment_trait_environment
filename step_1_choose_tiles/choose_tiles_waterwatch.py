import matplotlib
matplotlib.use('Agg')
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

class choose:
    def __init__(self):
        # seed random for reproducable results
        random.seed(123)
    def choose_subset(self,fname,huc_2s,out_dir,gdf_fname):
        '''
        function that randomly chooses one smaller huc watershed for each of
        the larger huc2 watersheds listed
        '''
        # we will need this later on
        self.huc_2s = huc_2s
        # which variables do we car about keeping?
        vars_keep = ['areasqkm','states','huc6','name','geometry']
        # loop over all the huc2s proveded
        chosen_hucs = pd.DataFrame(columns=vars_keep)
        for h,huc in enumerate(huc_2s):
            print('subsetting for HUC2 {}'.format(huc))
            # format the name correctly
            this_fname = fname.format(huc)
            # open the shapefile
            this_f  = gpd.read_file(this_fname)
            # get the index
            this_ind = this_f.index
            # pick a random smaller huc within the huc2
            this_huc_idx = random.randint(this_ind[0],this_ind[-1])
            # get the information for that smaller huc
            this_huc = this_f[vars_keep].loc[this_huc_idx]
            # add that to a list of all smaller hucs
            chosen_hucs.loc[huc] = this_huc
        # get the hucs that have been subsetted and save them
        subset_hucs = chosen_hucs['huc6']
        subset_hucs.to_csv(
            os.path.join(
                out_dir,'subset_hucs.csv'
            ),
            index=False
        )
        # put this into geodataframe for later use
        chosen_hucs = gpd.GeoDataFrame(chosen_hucs)
        # tell geopandas that this is currently in degree units
        chosen_hucs = chosen_hucs.set_crs('epsg:4326')
        # let's save the chosen hucs
        chosen_hucs.to_file(gdf_fname,driver='GeoJSON')
        self.chosen_hucs = chosen_hucs
        return chosen_hucs
    def get_catch_tiles(self,tile_fname,huc_subset):
        # get the catchment tile info
        tile_info = pd.read_csv(tile_fname)
        tile_info = tile_info.set_index('tile_id')
        tile_idx = tile_info.index
        print(tile_info)
        # for area calculations, we will need huc areas in km2, so lets do that
        # conversion here
        huc_subset_m = huc_subset.to_crs({'proj':'cea'})
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
        # get an array that is the smaller, subset hucs. we will use this as
        # the key for our intersection info dictionary
        subset_hucs = np.array(huc_subset['huc6'])
        for h,huc in enumerate(subset_hucs):
            intersection_info[huc] = [[],[],[]]
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
        for t,ti in enumerate(tile_idx):
            if ti%1000 == 0:
                print('checking intersection for tile {}'.format(ti))
            # get the polygon that we are talking about
            this_polygon_m = tile_gdf_m.loc[ti]['geometry']
            #print(this_polygon)
            #print(huc_subset_km)
            # check if this intersects with any of the huc boundaries
            intersects = list(huc_subset_m.intersects(this_polygon_m))
            for i,it in enumerate(intersects):
                if it:
                    #get the polygon for the intersecting sub-huc
                    this_huc_m = huc_subset_m['geometry'].loc[self.huc_2s[i]]
                    # get the overlapping area of these polygons
                    # all in km2
                    pixel_area_km = this_polygon_m.area/(10**6)
                    intersection_area_km = (
                        this_huc_m.intersection(this_polygon_m).area/(10**6)
                    )
                    perc_in_huc = intersection_area_km/pixel_area_km
                    intersection_info[subset_hucs[i]][0].append(ti)
                    intersection_info[subset_hucs[i]][1].append(perc_in_huc)
                    intersection_info[subset_hucs[i]][2].append(intersection_area_km)
                    all_intersecting = np.append(all_intersecting,ti)
        print(intersection_info)
        sys.exit()
        # sort the list of catchment pixels in ascending order
        all_intersecting = np.sort(all_intersecting)
        # make sure this list in unique--don't run a tile twice if it
        # intersects multiple watershed
        all_intersecting = np.unique(all_intersecting)
        return [intersection_info,all_intersecting,tile_gdf]
    def plot_hucs_tiles(self,subset,intersection_info,tile_shapes,
                        states_shp,tile_plot_fname):
        # open a figure
        plt.figure()
        # plot the united states outline
        # get the shapefile
        states = gpd.read_file(states_shp)
        states = states.to_crs('EPSG:4326')
        # get rid of non-conus states since not considering
        non_conus = ['HI','VI','MP','GU','AK','AS','PR']
        states_conus = states
        for n in non_conus:
            states_conus = states_conus[states_conus.STUSPS != n]
        # add this to the plot
        states_conus.boundary.plot()
        # for each huc in the subset
        for key in intersection_info.keys():
            # for each tile that intersects this huc
            for ti in intersection_info[key][0]:
                print('plotting tile {}'.format(ti))
                # get the tile geometry
                this_tile = tile_shapes['geometry'].loc[ti]
                # plot the tile geometry
                plt.plot(*this_tile.exterior.xy)
        # for each huc in the subset
        for idx in subset.index:
            # get geometry of this huc
            this_geom = subset['geometry'].loc[idx]
            # plot this geometry
            plt.plot(*this_geom.exterior.xy)
            # get the centroid and plot the huc number there
            centroid = this_geom.centroid
            print(centroid)
            print(idx)
            plt.text(centroid.x,centroid.y,idx)
        # save and close
        plt.savefig(tile_plot_fname,dpi=350,bbox_inches='tight')
        plt.close()
    def save_dictionary(self,dictionary,save_name):
        print('saving dictionary to {}'.format(save_name))
        with open(save_name,'wb') as f:
            pickle.dump(dictionary,f)
    def save_csv(self,to_save,save_name):
        as_df = pd.DataFrame(to_save)
        as_df.to_csv(save_name,index=False,header=False)
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



