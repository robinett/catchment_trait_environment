import os
import sys
from choose_tiles_waterwatch import choose

def main():
    # set main directory
    main_dir = (
        '/shared/pso/step_1_choose_tiles'
    )
    # set directory where data is stored
    data_dir = os.path.join(
        main_dir,'data'
    )
    # set directory to save outputs
    out_dir = os.path.join(
        main_dir,'outputs'
    )
    # set directory to save plots
    plot_dir = os.path.join(
        main_dir,'plots'
    )
    # file names for watershed boundaries
    # {} is filled in by HUC2 number for that data
    # can change number in file name for different huc sizes
    hucs_fname = os.path.join(data_dir,'HUC{}/Shape/WBDHU6.shp')
    # HUC-2s to choose a random representative from
    huc_2s = [
        1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18
    ]
    # file for catchment tile information
    tile_fname = os.path.join(data_dir,'tile_coord.csv')
    # file name for the final output dictionary
    dictionary_name = os.path.join(out_dir,'intersection_info.pkl')
    # file name for output of all Catchment tiles that need to be run
    all_tile_fname = os.path.join(out_dir,'intersecting_catch_tiles.csv')
    # name for the shapefile that contains state information for backgroun of plot
    state_shp = os.path.join(data_dir,'state_shp/tl_2022_us_state.shp')
    # name for where to save plot of tiles currently in the include file
    tile_plot_fname = os.path.join(plot_dir,'tiles_to_run.png')
    # name of where to save the formatted include file for pso run
    include_fname = os.path.join(out_dir,'include')
    # name of where to save geodataframe with all the chosen HUC6 information
    gdf_fname = os.path.join(out_dir,'subset_huc6s.geojson')
    # initiate an instance of class choose
    ch = choose()
    # get the subset of hucs we will use from each huc2
    subset = ch.choose_subset(hucs_fname,huc_2s,out_dir,gdf_fname)
    # get this catchment tile intersection infor from these
    intersection_info,all_intersecting,tile_shapes = (
        ch.get_catch_tiles(tile_fname,subset)
    )
    print('number of total tiles needed to run:')
    print(len(all_intersecting))
    # plot the hucs and tiles
    ch.plot_hucs_tiles(
        subset,intersection_info,tile_shapes,state_shp,tile_plot_fname
    )
    # save the catchment tile info to be used in future steps
    ch.save_dictionary(intersection_info,dictionary_name)
    # save the list of tiles to be run by Catchment in a later step
    ch.save_csv(all_intersecting,all_tile_fname)
    ch.create_include_file(all_intersecting,include_fname)
if __name__ == '__main__':
    main()
