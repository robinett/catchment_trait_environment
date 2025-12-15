import os
import datetime
from compare_lai import compare
import numpy as np

def main():
    # what is the name of the directory where we are located?
    base_dir = '/shared/pso/other_analyses/lai_compare'
    # where is the data located?
    data_dir = os.path.join(
        base_dir,
        'data'
    )
    # where should we save outputs from these scripts?
    out_dir = os.path.join(
        base_dir,
        'outputs'
    )
    # where do we save plots?
    plots_dir = os.path.join(
        base_dir,
        'plots'
    )
    # what are the start and end dates for comparison?
    start = datetime.date(2002,7,4)
    end = datetime.date(2015,12,31)
    # what is the filename hwere the lai values are stored?
    lai_fname = os.path.join(
        data_dir,
        'cn_lai.nc4'
    )
    # what catchment experiment are we going to use for our copmarison?
    # this must be an experiment that is already saved as a .pkl file
    catch_exp_dir = (
        '/shared/pso/step_5_analyze_outputs/saved_timeseries/'+
        'med-1996-2015.pkl'
    )
    # what is the file that has the tiles that we are running?
    tile_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/intersecting_catch_tiles.csv'
    )
    # what is the file that has the tile information?
    tile_info_fname = (
        '/shared/pso/step_1_choose_tiles/data/tile_coord.csv'
    )
    # where the lai fnames so that we can get dates?
    lai_fnames ='/shared/pso/other_analyses/lai_compare/data/lai_fnames.txt'
    # where are the cn lat vals?
    cn_lats = '/shared/pso/other_analyses/lai_compare/data/catch_cn_lats.txt'
    # where are the cn lon vals?
    cn_lons = '/shared/pso/other_analyses/lai_compare/data/catch_cn_lons.txt'
    # lets get an instance of the class and get the lai vals
    comp = compare()
    dates = comp.get_dates(lai_fnames)
    modis_lai_vals = comp.get_lai_truth(lai_fname)
    all_lons,all_lats = comp.get_lai_lons_lats(cn_lons,cn_lats)
    lons_mesh,lats_mesh = np.meshgrid(all_lons,all_lats)
    lons_flat = np.ndarray.flatten(lons_mesh)
    lats_flat = np.ndarray.flatten(lats_mesh)
    #comp.plot_modis_lai(lons_flat,lats_flat,modis_lai_vals,dates,plots_dir)
    catch_lai_vals = comp.get_lai_catch(catch_exp_dir,start,end)
    tiles,tiles_i,tiles_j,tile_info = comp.get_tiles_and_indices(
        tile_fname,tile_info_fname
    )
    comp_df,tile_lats,tile_lons = comp.get_comp_df(
        modis_lai_vals,catch_lai_vals,all_lons,all_lats,tile_info,tiles,dates
    )
    #comp.plot_comparison_timeseries(comp_df,tiles,plots_dir)
    tile_lons_flat = np.ndarray.flatten(tile_lons)
    tile_lats_flat = np.ndarray.flatten(tile_lats)
    comp.plot_comparison_map(
        comp_df,tiles,tile_lons_flat,tile_lats_flat,plots_dir
    )

if __name__ == '__main__':
    main()
