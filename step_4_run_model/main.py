from confirm_tiles import confirm

def main():
    # where are we?
    base_dir = '/shared/pso/step_4_run_model'
    # what are the tiles that Catchment is going to run?
    run_tiles_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/' +
        'intersecting_catch_tiles.csv'
    )
    # what are the tiles that we have point forcing for currently?
    # just give an example file
    point_force_example_fname = (
        '/shared/point_forcing_huc6/' +
        '1980_point_forcing_data.nc4'
    )
    # get an instance and check tiles
    co = confirm()
    run_tiles,force_tiles = co.get_tiles(
        run_tiles_fname,point_force_example_fname
    )
    co.check_tiles(run_tiles,force_tiles)

if __name__ == '__main__':
    main()
