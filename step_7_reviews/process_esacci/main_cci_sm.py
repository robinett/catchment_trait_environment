import os
import datetime
from process_cci_sm import process_cci_sm

def main():
    base_dir = (
        '/discover/nobackup/trobinet/from_aws/'
        'pso/step_7_reviews'
    )
    # root of ESA CCI ACTIVE v09.1 files (year subdirs inside)
    cci_root_dir = (
        '/discover/nobackup/trobinet/from_aws/pso/'
        'step_7_reviews/data/sm_raw/data'
    )
    # outputs/plots
    output_dir = os.path.join(base_dir, 'outputs')
    os.makedirs(output_dir, exist_ok=True)

    catch_tiles_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/'
        'step_1x_choose_tiles_large/outputs/intersecting_catch_tiles.csv'
    )
    catch_tile_info_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/'
        'step_1x_choose_tiles_large/data/tile_coord.csv'
    )

    # date span to materialize (continuous daily index)
    start_date = datetime.date(1991, 10, 1)
    end_date   = datetime.date(2014, 12, 31)

    example_cci_sm_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/'
        'step_7_reviews/data/sm_raw/data/1991/'
        'ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-19911001000000-fv09.1.nc'
    )

    # weights cache
    weights_fname = os.path.join(
        output_dir, 'cci_sm_weights_allpix.pkl'
    )
    # output truth CSV (selected/extrap tiles)
    truth_fname_selected_tiles = os.path.join(
        output_dir,
        'sm_truth_cci_active_pct_'
        '{start}_{end}_extrapolation_tiles.csv'
    )

    load_weights = True
    save_weights = True
    save_selected_truth = True

    pro = process_cci_sm(base_dir)
    pro.get_tiles(catch_tiles_fname)

    pro.get_weights(
        catch_tile_info_fname=catch_tile_info_fname,
        load_weights=load_weights,
        save_weights=save_weights,
        weights_fname=weights_fname,
        example_cci_sm_fname=example_cci_sm_fname
    )

    pro.create_truth_cci_sm(
        start=start_date,
        end=end_date,
        cci_root_dir=cci_root_dir,
        save_selected_truth=save_selected_truth,
        truth_fname_selected_tiles=truth_fname_selected_tiles
    )

if __name__ == '__main__':
    main()

