import os
import datetime
from process_oco2 import process_oco2

def main():
    base_dir = (
        '/discover/nobackup/trobinet/from_aws/'
        'pso/step_7_reviews'
    )
    oco2_dir = (
        '/discover/nobackup/trobinet/from_aws/pso/'
        'step_7_reviews/data/csif_raw/data'
    )
    # outputs/plots
    output_dir = os.path.join(base_dir, 'outputs')
    os.makedirs(output_dir, exist_ok=True)

    # catchment tiles & geometry (same as GLEAM)
    catch_tiles_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/'
        'step_1x_choose_tiles_large/outputs/intersecting_catch_tiles.csv'
    )
    catch_tile_info_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/'
        'step_1x_choose_tiles_large/data/tile_coord.csv'
    )

    # date span to materialize (continuous)
    start_date = datetime.date(2001, 1, 1)
    end_date   = datetime.date(2014, 12, 31)

    # file patterns and names
    # example file to read lat/lon
    example_oco2_fname = os.path.join(
        oco2_dir, 'OCO2.SIF.all.daily.2001.nc'
    )
    # yearly files to open during processing
    oco2_fname_fmt = os.path.join(
        oco2_dir, 'OCO2.SIF.all.daily.{:04d}.nc'
    )

    # weights cache
    weights_fname = os.path.join(
        output_dir, 'oco2_weights_allpix.pkl'
    )
    # output truth CSV (selected/extrap tiles)
    truth_fname_selected_tiles = os.path.join(
        output_dir,
        'sif_truth_oco2_all_daily_'
        '{start}_{end}_extrapolation_tiles.csv'
    )

    # toggles
    load_weights = False
    save_weights = True
    save_selected_truth = True

    pro = process_oco2(base_dir)
    pro.get_tiles(catch_tiles_fname)

    pro.get_weights(
        catch_tile_info_fname=catch_tile_info_fname,
        load_weights=load_weights,
        save_weights=save_weights,
        weights_fname=weights_fname,
        example_oco2_fname=example_oco2_fname
    )

    pro.create_truth_oco2(
        start=start_date,
        end=end_date,
        oco2_fname_fmt=oco2_fname_fmt,
        save_selected_truth=save_selected_truth,
        truth_fname_selected_tiles=truth_fname_selected_tiles
    )

if __name__ == '__main__':
    main()

