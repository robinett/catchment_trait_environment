import sys
sys.path.append('/shared/pso/step_5.1_analyze_outputs/funcs')
sys.path.append('/shared/pso/step_5.1_analyze_outputs/tasks')
import os
import datetime
from get_timeseries import get_timeseries
from rainfall_runoff_ratio import rainfall_runoff
from plot_timeseries import timeseries
from get_averages_and_error import averages_and_error
from plot_other import plot_other

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
    # what are the start and end dates for when we should compute error?
    start_err = datetime.date(1995,1,1)
    end_err = datetime.date(2014,12,31)
    # where are the analysis pixels located?
    pixels_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/intersecting_catch_tiles.csv'
    )
    # where is the intersection info located?
    intersection_info_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/intersection_info.pkl'
    )
    # where is the le truth located?
    le_truth_fname = (
        '/shared/pso/step_3_process_gleam/outputs/' +
        'le_truth_gleam_38a_watts_per_m2_1995-01-01_2014-12-31_' +
        'camels_tiles.csv'
    )
    # where is the streamflow truth located?
    streamflow_truth_fname = (
        '/shared/pso/step_3.1.1_process_camels/outputs/' +
        'camels_truth_yearly_1995-01-01_2014-12-31_mm_day.csv'
    )
    # what are the names and info of the timeseries that we are going to load?
    timeseries_info = {
        'med-default-pft-g1-1992-2014':{
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
    # let's get the le truth timeseries
    g = get_timeseries()
    le_obs = g.get_le_obs(le_truth_fname)
    # let's get the streamflow truth timeseries
    strm_obs = g.get_streamflow_obs(streamflow_truth_fname)
    # load the timeseries
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
        # get the pixel-averaged and error information for pixel-scale outputs
        get_metrics = averages_and_error()
        #pixel_avgs = get_metrics.var_avgs(this_pixel_raw_timeseries)
        #timeseries_info[load]['pixel_avgs'] = pixel_avgs
        # get the watershed-averaged and error information for watershed-scale
        # ouptts
        #wat_avgs = get_metrics.var_avgs(this_wat_raw_timeseries)
        #timeseries_info[load]['wat_avgs'] = wat_avgs
        #this_le_err_df = get_metrics.var_error(
        #    this_pixel_raw_timeseries['le'],le_obs,start_err,end_err
        #)
        #timeseries_info[load]['pixel_le_errors'] = this_le_err_df
        #this_wat_err_df = get_metrics.var_error(
        #    this_wat_raw_timeseries['strm_yr'],strm_obs,start_err,end_err
        #)
        #timeseries_info[load]['wat_strm_errors'] = this_wat_err_df
        timestep_le_err = get_metrics.get_timestep_error(
            this_pixel_raw_timeseries['le'],le_obs,start_err,end_err
        )
        timestep_strm_err = get_metrics.get_timestep_error(
            this_wat_raw_timeseries['strm_yr'],strm_obs,start_err,end_err
        )
        # cool we have the information we need
        # let's make the error distribution historgrams
        p = plot_other()
        #cols = list(timestep_le_err['err'].columns)
        #bin_width = 1
        #for c,col in enumerate(cols):
        #    this_err = timestep_le_err['err'][col]
        #    this_abs_err = timestep_le_err['abs_err'][col]
        #    print('making histogram for {}'.format(col))
        #    p.histogram(
        #        this_err,plots_dir,'err_le_hist_{}'.format(col),
        #        x_label='LE error (W/m2)',y_label='num of pixels',
        #        bin_width=bin_width
        #    )
        #    p.histogram(
        #        this_abs_err,plots_dir,'abs_err_le_hist_{}'.format(col),
        #        x_label='absolute LE error (W/m2)',y_label='num of days',
        #        bin_width=bin_width
        #    )
        cols = list(timestep_strm_err['err'].columns)
        num_bins = 7
        for c,col in enumerate(cols):
            this_err = timestep_strm_err['err'][col]
            this_abs_err = timestep_strm_err['abs_err'][col]
            print('making histogram for {}'.format(col))
            p.histogram(
                this_err,plots_dir,'err_strm_hist_{}'.format(col),
                x_label='streamflow error (mm/day)',y_label='num of years',
                num_bins=num_bins
            )
            p.histogram(
                this_abs_err,plots_dir,'abs_err_strm_hist_{}'.format(col),
                x_label='absolute streamflow error (mm/day)',y_label='num of years',
                num_bins=num_bins
            )
    #cool, we have the information that we need
    # let's make a histogram of our MAEs for both ET and streamflow

if __name__ == '__main__':
    main()
