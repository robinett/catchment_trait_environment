import os
import sys
import datetime
from process_camels import camels

def main():
    # where are we now?
    base_dir = '/shared/pso/step_3.1.1x_process_camels'
    # where is the data located?
    data_dir = '/shared/pso/step_3.1.1_process_camels/data'
    # where will we save outputs?
    out_dir = os.path.join(
        base_dir,
        'outputs'
    )
    # what is the filename for the camels data?
    camels_strm_dir = os.path.join(
        data_dir,
        'basin_dataset_public_v1p2',
        'usgs_streamflow'
    )
    # what is the filename for the intersecting tiles?
    tiles_fname = (
        '/shared/pso/step_1x_choose_tiles_large/outputs/intersecting_catch_tiles.csv'
    )
    # what is the filename for the camels basins we have chosen?
    chosen_camels_fname = (
        '/shared/pso/step_1x_choose_tiles_large/outputs/chosen_camels.csv'
    )
    # what are the start and end dates?
    start = datetime.date(1995,1,1)
    end = datetime.date(2014,12,31)
    # let's start the analysis
    c = camels(start,end)
    chosen_camels_df = c.get_camels(chosen_camels_fname)
    camels_truth = c.get_streamflow(
        camels_strm_dir,chosen_camels_df
    )
    c.save_truth(camels_truth,out_dir)


if __name__ == '__main__':
    main()
