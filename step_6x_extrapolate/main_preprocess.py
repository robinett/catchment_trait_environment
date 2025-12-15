import os
import sys
from extrapolation import extrapolate
import numpy as np
import pandas as pd
import datetime
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/pso/step_5.1_analyze_outputs/funcs'
)
from general_functions import gen_funcs
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large'
)
from choose_tiles_camels import choose
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/pso/step_2x_env_covariates'
)
from get_env_covariates import get_covariates

def main():
    # where are we?
    base_dir = '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate'
    # where do we put outputs?
    out_dir = os.path.join(
        base_dir,
        'outputs'
    )
    # where do we put plots?
    plots_dir = os.path.join(
        base_dir,
        'plots'
    )
    # where is the data
    data_dir = os.path.join(
        base_dir,
        'data'
    )
    # where are the pixels that we already ran located?
    previous_pixels_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/' +
        'step_1x_choose_tiles_large/' +
        'outputs/intersecting_catch_tiles.csv'
    )
    # where is the conus outline located?
    us_outline_fname = os.path.join(
        data_dir,
        'cb_2023_us_nation_20m'
    )
    # what is the info for all pixels?
    all_pixel_info_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/' +
        'step_1_choose_tiles/data/tile_coord.csv'
    )
    # where are we going to save the conus pixels gdf?
    conus_pixels_gdf_fname = 'conus_pixels.gdf'
    # where are we going to save the extrap pixels gdf?
    extrap_pix_gdf_fname = os.path.join(
        out_dir,'extrap_pixels.gpkg'
    )
    # where are we going to save the plot of all conus pixels?
    conus_pixels_map_fname = 'conus_pixels.png'
    # where do we save the list of all the pixels in conus?
    conus_pix_save_fname = os.path.join(
        out_dir,'conus_pix.csv'
    )
    # where is the preicp info located?
    precip_fname = os.path.join(
        data_dir,
        'gldas_precip_kg_m2_s.csv'
    )
    # what are the dates over which to average gldas?
    precip_start = np.datetime64('1985-01-01')
    precip_end = np.datetime64('2014-12-31')
    # where is the canopy info located?
    canopy_fname = os.path.join(
        data_dir,
        'canopy_height.csv'
    )
    # where should we save precip info?
    precip_out_fname = os.path.join(
        out_dir,
        'gldas_avg_precip.nc4'
    )
    # where should we save precip info normalized?
    precip_norm_out_fname = os.path.join(
        out_dir,
        'gldas_avg_precip_normalized.nc4'
    )
    # where should we save canopy info?
    canopy_out_fname = os.path.join(
        out_dir,
        'canopy_height.nc4'
    )
    # where should we save canopy info normalized?
    canopy_norm_out_fname = os.path.join(
        out_dir,
        'canopy_height_normalized.nc4'
    )
    # get all pixels that are conus
    e = extrapolate()
    # get just conus as opposed to all of US
    conus_gdf = e.get_conus_gdf(us_outline_fname,plots_dir)
    # get pixels that intersect with CONUS
    conus_tile_gdf,conus_pix = e.get_conus_pix(
        conus_gdf,all_pixel_info_fname,out_dir,
        conus_pixels_gdf_fname,plots_dir,
        conus_pixels_map_fname,save_or_load='load'
    )
    # let's save the conus pixels so we can use them later
    c = choose()
    c.save_array_as_csv(conus_pix,conus_pix_save_fname)
    # get rid of pixels that we have already run
    g = gen_funcs()
    print('number of conus pix:')
    print(len(conus_pix))
    previous_pix = g.get_pixels(previous_pixels_fname)
    remaining_pix = [
        p for p in conus_pix if p not in previous_pix
    ]
    print('number of pix after removing previous:')
    print(len(remaining_pix))
    # let's create a gdf of this and plot just to be safe
    rem_conus_tile_gdf = conus_tile_gdf.loc[remaining_pix]
    c.save_gdf(
        rem_conus_tile_gdf,
        extrap_pix_gdf_fname
    )
    e.plot_gpds(
        [conus_gdf,rem_conus_tile_gdf],
        plots_dir,
        'remaining_conus_pix.png'
    )
    # create an include file from this
    c.create_include_file(
        remaining_pix,
        os.path.join(
            out_dir,
            'include_extrapolation'
        )
    )
    # save the tiles that we are going to run as a .csv as well
    extrapolation_pix_fname = os.path.join(
        out_dir,
        'extrapolation_pix.csv'
    )
    c.save_array_as_csv(
        remaining_pix,
        extrapolation_pix_fname
    )
    # yeehaw. Now for our EF extrapolation, we're going to need evironmental
    # covariates
    g_c = get_covariates(tile_overide=np.array(remaining_pix,dtype=int))
    g_c.get_gldas_info(
        precip_start,
        precip_end,
        precip_fname,
        precip_out_fname,
        precip_norm_out_fname
    )
    g_c.get_canopy_info(
        canopy_fname,
        canopy_out_fname,
        canopy_norm_out_fname
    )











if __name__ == '__main__':
    main()
