import sys
sys.path.append('/shared/pso/step_5.1_analyze_outputs/funcs')
sys.path.append('/shared/pso/step_5.1_analyze_outputs/tasks')
import os
import datetime
import pandas as pd
import numpy as np
from get_timeseries import get_timeseries
from rainfall_runoff_ratio import rainfall_runoff
from plot_timeseries import timeseries
from get_averages_and_error import averages_and_error
from plot_other import plot_other
np.random.seed(123)

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
    # what are the start and end dates for our error metrics?
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
    # where is the streamflow truth located?
    streamflow_truth_fname = (
        '/shared/pso/step_3.1.1_process_camels/outputs/' +
        'camels_truth_yearly_1995-01-01_2014-12-31_mm_day.csv'
    )
    # where is le truth?
    le_obs_fname = (
        '/shared/pso/step_3_process_gleam/outputs/le_truth_gleam_38a_watts' +
        '_per_m2_1995-01-01_2014-12-31_camels_tiles.csv'
    )
    # what are the names and info of the timeseries that we are going to load?
    timeseries_info = {}
    num_particles = 30
    for p in range(num_particles):
        this_name = (
            'g1-a0-a1-et-str-spin19921994-test19952014-init-{}'.format(p)
        )
        timeseries_info[this_name] = {
            'dir':(
                '/lustre/catchment/exps/' +
                'GEOSldas_CN45_pso_g1_a0_a1_et_strm_camels_spin19921994_' +
                'test19952014_oneit/{}'.format(p)
            ),
            'load_or_save':'load',
            'default_type':'default',
            'read_met_forcing':False,
            'timeseries_dir':os.path.join(
                save_dir,this_name
            )
        }
    # load the timeseries
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
        # get some le error statistics
        get_metrics = averages_and_error()
    # let's get the streamflow truth timeseries
    strm_obs = g.get_streamflow_obs(streamflow_truth_fname)
    # get the le truth timeseries
    le_obs = g.get_le_obs(le_obs_fname)
    # let's plot the all the different initial particles
    #start_plot = datetime.date(2011,1,1)
    #end_plot = datetime.date(2012,12,31)
    #to_plot = []
    #small_preds = []
    #names = []
    #for p in range(num_particles):
    #    to_plot.append(
    #        timeseries_info[to_load[p]]['pixel_raw_timeseries']
    #    )
    #    small_preds.append(0.3)
    #    names.append(to_load[p])
    #to_plot.append({'le':le_obs})
    #small_preds.append(1)
    #obs_idx = num_particles
    #names.append('GLEAM')
    #p_t = timeseries()
    #p_t.plot_one_var(
    #    to_plot,names,'le','W/m2','le_init_spread',plots_dir,
    #    small_preds=small_preds,obs_idx=obs_idx,
    #    start=start_plot,end=end_plot
    #)
    g_error = averages_and_error()
    err_bias_df = pd.DataFrame(
        columns = ['ubMAE','MAE','additive','multiplicative']
    )
    pixels = np.array(
        timeseries_info[to_load[0]]['pixel_raw_timeseries']['le'].columns
    )
    num_to_choose = 200
    pixels_chosen = np.random.choice(pixels,size=num_to_choose)
    plot_o = plot_other()
    avg_pix_add_bias = np.zeros(len(pixels_chosen))
    avg_pix_abs_add_bias = np.zeros(len(pixels_chosen))
    avg_pix_mae = np.zeros(len(pixels_chosen))
    avg_pix_ubmae = np.zeros(len(pixels_chosen))
    avg_pix_diff_mae_ubmae = np.zeros(len(pixels_chosen))
    avg_pix_abs_diff_mae_ubmae = np.zeros(len(pixels_chosen))
    for p,pix in enumerate(pixels_chosen):
        print(
            'analysis for pix {}; this is {} out of {} pixels'.format(
                pix,p,num_to_choose
            )
        )
        all_add_bias = np.zeros(len(to_load))
        all_mult_bias = np.zeros(len(to_load))
        all_mae = np.zeros(len(to_load))
        all_ubmae = np.zeros(len(to_load))
        for l,load in enumerate(to_load):
            add_bias,mult_bias = g_error.get_add_and_mult_bias(
                timeseries_info[load]['pixel_raw_timeseries']['le'][pix],
                le_obs[pix],
                plots_dir,
                load+'_{}'.format(pix),
                start=start_err,
                end=end_err,
                plot=True
            )
            all_add_bias[l] = add_bias
            all_mult_bias[l] = mult_bias
            le_error = g_error.var_error(
                timeseries_info[load]['pixel_raw_timeseries']['le'],
                le_obs,
                start_err,
                end_err
            )
            all_mae[l] = le_error[pix].loc['mae']
            all_ubmae[l] = le_error[pix].loc['ubmae']
        diff_mae_ubmae = all_mae - all_ubmae
        abs_diff_mae_ubmae = np.abs(diff_mae_ubmae)
        abs_add_bias = np.abs(all_add_bias)
        save_name = 'mae_vs_add_bias_{}.png'.format(pix)
        plot_o.scatter(
            all_mae,all_add_bias,plots_dir,save_name,'MAE (W/m2)',
            'additive bias (W/m2)', one_to_one_line=True
        )
        save_name = 'mae_vs_abs_add_bias_{}.png'.format(pix)
        plot_o.scatter(
            all_mae,abs_add_bias,plots_dir,save_name,'MAE (W/m2)',
            '|additive bias| (W/m2)', one_to_one_line=True
        )
        save_name = 'add_bias_vs_mae_minus_ubmae_{}.png'.format(pix)
        plot_o.scatter(
            all_add_bias,diff_mae_ubmae,plots_dir,save_name,
            'additive bias (W/m2)',
            'MAE - ubMAE (W/m2)',one_to_one_line=True
        )
        save_name = 'abs_add_bias_vs_mae_minus_ubmae_{}.png'.format(pix)
        plot_o.scatter(
            abs_add_bias,diff_mae_ubmae,plots_dir,save_name,
            '|additive bias| (W/m2)',
            'MAE - ubMAE (W/m2)',one_to_one_line=True
        )
        save_name = 'mae_vs_ubmae_{}.png'.format(pix)
        plot_o.scatter(
            all_mae,all_ubmae,plots_dir,save_name,'MAE (W/m2)',
            'ubMAE (W/m2)',one_to_one_line=True
        )
        avg_add_bias = np.mean(all_add_bias)
        avg_pix_add_bias[p] = avg_add_bias
        avg_pix_abs_add_bias[p] = np.abs(avg_add_bias)
        avg_mae = np.mean(all_mae)
        avg_pix_mae[p] = avg_mae
        avg_ubmae = np.mean(all_ubmae)
        avg_pix_ubmae[p] = avg_ubmae
        avg_diff_mae_ubmae = np.mean(diff_mae_ubmae)
        avg_pix_diff_mae_ubmae[p] = avg_diff_mae_ubmae
        avg_pix_abs_diff_mae_ubmae[p] = np.abs(avg_diff_mae_ubmae)
    save_name = 'mae_vs_ubmae_all.png'
    plot_o.scatter(
        avg_pix_mae,avg_pix_ubmae,plots_dir,save_name,
        'average LE MAE across 30 runs (W/m2)',
        'average LE ubMAE across 30 runs (W/m2)',
        one_to_one_line=True
    )
    save_name = 'mae_vs_add_bias_all.png'
    plot_o.scatter(
        avg_pix_mae,avg_pix_add_bias,plots_dir,save_name,
        'average LE MAE across 30 runs (W/m2)',
        'average LE add. bias across 30 runs (W/m2)',
        one_to_one_line=True
    )
    save_name = 'mae_vs_abs_add_bias_all.png'
    plot_o.scatter(
        avg_pix_mae,avg_pix_abs_add_bias,plots_dir,save_name,
        'average LE MAE across 30 runs (W/m2)',
        '|average LE add. bias| across 30 runs (W/m2)',
        one_to_one_line=True
    )
    save_name = 'diff_mae_ubmae_vs_add_bias_all.png'
    plot_o.scatter(
        avg_pix_diff_mae_ubmae,avg_pix_add_bias,plots_dir,save_name,
        'average (LE MAE - LE ubMAE) across 30 runs (W/m2)',
        'average LE add. bias across 30 runs (W/m2)',
        one_to_one_line=True
    )
    save_name = 'abs_diff_mae_ubmae_vs_abs_add_bias_all.png'
    plot_o.scatter(
        avg_pix_abs_diff_mae_ubmae,avg_pix_abs_add_bias,plots_dir,save_name,
        '|average (LE MAE - LE ubMAE)| across 30 runs (W/m2)',
        '|average LE add. bias| across 30 runs (W/m2)',
        one_to_one_line=True
    )

if __name__ == '__main__':
    main()
