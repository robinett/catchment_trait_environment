from get_timeseries import get_timeseries
import datetime
import os
from analyze_timeseries_pixel_scale import analyze_pix
from plot import plot

def main():
    # where is this script?
    base_dir = '/shared/pso/step_5_analyze_outputs'
    # where will we be saving plots from this script?
    plots_dir = os.path.join(
        base_dir,'plots'
    )
    # what experiment are we loading for comparison?
    default_name = 'med-2006-2006'
    # what are we going to call the bonetti timeseries?
    bonetti_name = 'bonetti-2006-2006'
    # where is the default experiment located?
    default_exp_dir = (
        '/lustre/catchment/exps/GEOSldas_CN45_med_default_2006/0'
    )
    # where is the bonettit timeseries located?
    bonetti_exp_dir = (
        '/lustre/catchment/exps/GEOSldas_CN45_bonetti_ksat/0'
    )
    # where is the fluxcom timeseries (for comparison)
    fluxcom_dir = (
        '/shared/pso/step_3_process_fluxcom/outputs/'+
        'le_truth_fluxcom_rs_meteo_ensemble_watts_per_m2_2001-01-01_'+
        '2009-12-31_selected_tiles.csv'
    )
    # where are the names of the current pixels located
    pixel_dir = (
        '/shared/pso/step_1_choose_tiles/outputs/intersecting_catch_tiles.csv'
    )
    # information for creating/saving/loading timeseries. Default experiment
    # information first, then bonetti information
    # should we save the timeseries?
    save_timeseries = [False,False]
    # should we load the bonetti timeseries?
    load_timeseries = [True,True]
    # is this experiment in the default format
    default_format = [True,True]
    # does this experiment have met forcing output associated with it?
    read_met_forcing = [False,False]
    # put directories and names together for loading purposes
    all_names = [default_name,bonetti_name]
    all_dir = [default_exp_dir,bonetti_exp_dir]
    all_timeseries = {}
    # where are timeseries saved?
    timeseries_dir = os.path.join(
        base_dir,'saved_timeseries'
    )
    # set the start and end dates (inclusive)
    start = datetime.date(2006,1,1)
    end = datetime.date(2006,12,31)
    # let's get the bonetti and default timeseries
    get = get_timeseries(start,end)
    pixels = get.get_pixels(pixel_dir)
    for e,exp in enumerate(all_names):
        all_timeseries[exp] = get.get_catch_timeseries(
            all_dir[e],all_names[e],save_timeseries[e],load_timeseries[e],
            default_format[e],read_met_forcing[e],timeseries_dir
        )
    # let's get the fluxcom timeseries just for comparison
    fluxcom_timeseries = get.get_fluxcom_timeseries(fluxcom_dir)
    # now let's see how these two timeseries differ
    a_pix = analyze_pix()
    # these are both general experiments, not pso experiments, so we need to
    # tell this to the function
    experiment_type = 'general'
    default_df = a_pix.get_rmse_dict(
        all_names[0],all_timeseries[all_names[0]],fluxcom_timeseries,
        experiment_type
    )
    bonetti_df = a_pix.get_rmse_dict(
        all_names[1],all_timeseries[all_names[1]],fluxcom_timeseries,
        experiment_type
    )
    # let's make some plots to see how this changed ET and runoff
    a_pix.plot_general_comparison_maps(
        all_names,default_df,bonetti_df,plots_dir
    )

if __name__ == '__main__':
    main()
