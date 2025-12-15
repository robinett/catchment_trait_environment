from process_streamflow import process_streamflow
import os
import datetime

def main():
    # set the base dir
    base_dir = '/shared/pso/step_3.1_process_streamflow'
    # set the dir where the data is stored
    data_dir = os.path.join(base_dir,'data')
    out_dir = os.path.join(base_dir,'outputs')
    # set the location where the streamflow raw data is stored
    strm_data_dir = os.path.join(
        data_dir,'mv01d_row_data.txt'
    )
    strm_flag_dir = os.path.join(
        data_dir,'mv01d_row_flag.txt'
    )
    tile_dir = (
        '/shared/pso/step_1_choose_tiles/outputs/intersecting_catch_tiles.csv'
    )
    huc_dir = (
        '/shared/pso/step_1_choose_tiles/outputs/subset_hucs.csv'
    )
    save_strm_data_dir = os.path.join(
        out_dir,'streamflow_data_mm_per_month_{start}_{end}.csv'
    )
    save_strm_flags_dir = os.path.join(
        out_dir,'streamflow_data_flags_{start}_{end}.csv'
    )
    # start and end date that we want for the data (data is montly resolution)
    start = datetime.date(2002,1,1) # inclusive
    end = datetime.date(2002,12,31) # inclusive
    # create an instance of the class
    pro = process_streamflow()
    # get tiles of interest
    tiles = pro.get_tiles(tile_dir)
    # get the subsetted hucs of interest
    hucs = pro.get_hucs(huc_dir)
    # get the streamflow data
    strm_data = pro.get_streamflow(strm_data_dir)
    # get the streamflow flags
    strm_flags = pro.get_streamflow(strm_flag_dir)
    # get data at only times and hucs of interest
    strm_data_trim = pro.trim_data(strm_data,start,end,hucs)
    strm_flags_trim = pro.trim_data(strm_flags,start,end,hucs)
    # save the data to csvs
    pro.save_data(strm_data_trim,save_strm_data_dir,start,end)
    pro.save_data(strm_flags_trim,save_strm_flags_dir,start,end)

if __name__ == '__main__':
    main()
