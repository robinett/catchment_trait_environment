import os
from choose_tiles_unregulated import choose_unreg
import numpy as np
import sys

def main():
    # where is this script located?
    base_dir = '/shared/pso/step_1_choose_tiles'
    # where is the data located?
    data_dir = os.path.join(
        base_dir,
        'data'
    )
    # where is the station infor file located?
    s_info_fname = os.path.join(
        data_dir,
        'Riverflow_Station_Information_NA-M09.nc'
    )
    # where is the info on the Catchment tiles?
    tile_info_fname = os.path.join(
        data_dir,
        'tile_coord_M09.csv'
    )
    # where is the shapefile for the watershed boundaries for our potential
    # stations of interest?
    boundaries_fname = os.path.join(
        data_dir,
        'SW-Trends-watershed-bounds-shapefile'
    )
    # where is the shapefile for the states
    states_shp_fname = os.path.join(
        data_dir,
        'state_shp'
    )
    # where are we putting plots?
    plots_dir = os.path.join(
        base_dir,
        'plots'
    )
    # how many random stations do we want to use as our truth?
    num_stations = 30
    # do we want to run the analysis to confirm the data from the .nc file?
    confirm_nc_info = True

    # initiate a class
    c_unreg = choose_unreg()
    s_info = c_unreg.get_basin_info(s_info_fname)
    tile_info = c_unreg.get_tile_info(tile_info_fname)
    boundary_info = c_unreg.get_boundary_info(boundaries_fname)
    station_names = c_unreg.get_station_names(s_info)
    if confirm_nc_info:
        boundary_stations = np.array(boundary_info['BasinID'])
        c_unreg.confirm_nc(
            station_names,s_info,tile_info,boundary_info,states_shp_fname,
            plots_dir
        )
    selected = c_unreg.select_random_stations(station_names,num_stations)
    #c_ungagued.plot_selected_tiles(stations,intersecting_tiles)

if __name__ == '__main__':
    main()
