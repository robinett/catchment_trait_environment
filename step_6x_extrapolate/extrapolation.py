import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import copy

class extrapolate:
    def get_conus_pix(self,conus_gdf,pixel_fname,out_dir,save_name,
                      plots_dir,plot_name,save_or_load='save'):
        if save_or_load == 'save':
            conus_geo = conus_gdf['geometry']
            print(conus_geo)
            pixel_info = pd.read_csv(pixel_fname)
            pixel_info = pixel_info.set_index('tile_id')
            all_pixels = list(pixel_info.index)
            print(pixel_info)
            # things we will use to hold
            conus_pixels = np.zeros(0)
            conus_pixels_polygons = []
            for p,pi in enumerate(all_pixels):
                # get lat and lon for the four points that will define this polygon
                # top left, top right, bottom right, bottom left
                if pi%1000 == 0:
                    print('creating polygon for pixel {}'.format(pi))
                this_lon_vals = [
                    pixel_info['min_lon'].loc[pi],
                    pixel_info['max_lon'].loc[pi],
                    pixel_info['max_lon'].loc[pi],
                    pixel_info['min_lon'].loc[pi]
                ]
                this_lat_vals = [
                    pixel_info['max_lat'].loc[pi],
                    pixel_info['max_lat'].loc[pi],
                    pixel_info['min_lat'].loc[pi],
                    pixel_info['min_lat'].loc[pi]
                ]
                # make this lat/lon into a polygon
                this_polygon = Polygon(
                    zip(this_lon_vals,this_lat_vals)
                )
                # check its intersecpion area with conus
                this_intersection_area = this_polygon.intersection(
                    conus_geo
                ).area
                this_pixel_area = this_polygon.area
                this_intersection_area = this_intersection_area.iloc[0]
                intersection_perc = this_intersection_area/this_pixel_area
                # if >0 add it to a list
                if intersection_perc >= 0.75:
                    conus_pixels = np.append(conus_pixels,pi)
                    conus_pixels_polygons.append(this_polygon)
            # put these polygons into geodataframe
            df = pd.DataFrame({
                'pixel':conus_pixels
            })
            conus_tile_gdf = gpd.GeoDataFrame(
                df,geometry=conus_pixels_polygons,crs='EPSG:4326'
            )
            conus_tile_gdf = conus_tile_gdf.set_index('pixel')
            print(conus_tile_gdf)
            conus_tile_gdf.to_file(
                os.path.join(
                    out_dir,
                    save_name
                )
            )
        elif save_or_load == 'load':
            conus_tile_gdf = gpd.read_file(
                os.path.join(
                    out_dir,
                    save_name
                )
            )
            conus_tile_gdf = conus_tile_gdf.set_index('pixel')
            conus_pixels = list(conus_tile_gdf.index)
        gpds_to_plot = [
            conus_gdf,
            conus_tile_gdf
        ]
        self.plot_gpds(gpds_to_plot,plots_dir,plot_name)
        return [conus_tile_gdf,conus_pixels]
    def plot_gpds(self,gpds,plots_dir,save_name):
        plt.figure()
        for g,gpd in enumerate(gpds):
            this_rows = list(gpd.index)
            for i,idx in enumerate(this_rows):
                this_geo = gpd['geometry'].loc[idx]
                if this_geo.geom_type == 'Polygon':
                    plt.plot(*this_geo.exterior.xy)
                elif this_geo.geom_type == 'MultiPolygon':
                    for this_this_geo in this_geo.geoms:
                        plt.plot(*this_this_geo.exterior.xy)
                else:
                    raise Exception('This shape is not a Polygon')
        plt.savefig(
            os.path.join(
                plots_dir,
                save_name
            ),
            dpi=350
        )
        plt.close()
    def get_conus_gdf(self,us_shp_fname,plots_dir):
        us_gdf = gpd.read_file(us_shp_fname)
        us_geom = us_gdf['geometry'].iloc[0]
        largest_area = 0
        for this_geo in us_geom.geoms:
            this_area = this_geo.area
            if this_area > largest_area:
                largest_area = copy.deepcopy(this_area)
                conus_geom = copy.deepcopy(this_geo)
        df = pd.DataFrame({
            'num':[0]
        })
        conus_gdf = gpd.GeoDataFrame(
            df,geometry=[conus_geom],crs='EPSG:4326'
        )
        conus_gdf = conus_gdf.set_index('num')
        #self.plot_gpds([conus_gdf],plots_dir,'conus_outline.png')
        return conus_gdf














