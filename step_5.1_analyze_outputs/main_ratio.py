import sys
sys.path.append('/shared/pso/step_5.1_analyze_outputs/funcs')
sys.path.append('/shared/pso/step_5.1_analyze_outputs/tasks')
import os
import datetime
from get_timeseries import get_timeseries
from rainfall_runoff_ratio import rainfall_runoff

def main():
    # what is the base dir?
    base_dir = '/shared/pso/step_5.1_analyze_outputs'
    # where should we save plots?
    plots_dir = os.path.join(
        base_dir,'plots'
    )
    # where should we save computed outputs?
    out_dir = os.path.join(
        base_dir,'outputs'
    )
    save_dir = (
        '/shared/pso/step_5_analyze_outputs/saved_timeseries'
    )
    # what is the start date?
    start = datetime.date(1992,1,1)
    # what is the end date?
    end = datetime.date(2014,12,31)
    # where are the analysis pixels located?
    pixels_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/intersecting_catch_tiles.csv'
    )
    # where is the intersection info located?
    intersection_info_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/intersection_info.pkl'
    )
    # where is the streamflow truth located?
    streamflow_truth_fname = (
        '/shared/pso/step_3.1.1_process_camels/outputs/' +
        'camels_truth_yearly_1995-01-01_2014-12-31_mm_day.csv'
    )
    # what do we want to save the ratio as?
    save_ratio_fname = os.path.join(
        out_dir,
        'streamflow_over_rainfall.csv'
    )
    # what are the names and info of the timeseries that we are going to load?
    timeseries_info = {
        'test_med_default_2006':{
            'dir':(
                '/lustre/catchment/exps/' +
                'GEOSldas_CN45_med_default_pft_g1_1992_2014_camels/0'
            ),
            'load_or_save':'load',
            'default_type':'default',
            'read_met_forcing':True,
            'timeseries_dir':os.path.join(
                save_dir,'med-default-pft-g1-1992-2014'
            )
        }
    }
    # load the timeseries
    g = get_timeseries()
    to_load = list(timeseries_info.keys())
    for l,load in enumerate(to_load):
        # get the raw timeseries at the pixel scales
        this_pixel_raw_timeseries = g.get_catchcn(
            timeseries_info[load]['dir'],
            timeseries_info[load]['load_or_save'],
            timeseries_info[load]['default_type'],
            timeseries_info[load]['read_met_forcing'],
            timeseries_info[load]['timeseries_dir'],
            start,
            end
        )
        timeseries_info[load]['pixel_raw_timeseries'] = this_pixel_raw_timeseries
        # get the raw timeseries at the watershed scale
        this_wat_raw_timeseries = g.get_watershed(
            timeseries_info[load]['pixel_raw_timeseries'],
            timeseries_info[load]['read_met_forcing'],
            intersection_info_fname
        )
        timeseries_info[load]['wat_raw_timeseries'] = this_wat_raw_timeseries
    # let's get the streamflow truth timeseries
    strm_obs = g.get_streamflow_obs(streamflow_truth_fname)
    # let's use the functions to plot the ratio of streamflow to precipitation
    # as a scatter and a map
    rainfall_runoff(
        timeseries_info[to_load[0]]['wat_raw_timeseries']['rainfsnowf_yr'],
        strm_obs,
        plots_dir,
        save_ratio_fname
    )

if __name__ == '__main__':
    main()
