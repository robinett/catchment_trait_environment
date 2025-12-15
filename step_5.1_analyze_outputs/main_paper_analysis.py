import sys
sys.path.append('/shared/pso/step_5.1_analyze_outputs/funcs')
sys.path.append('/shared/pso/step_5.1_analyze_outputs/tasks')
import os
import datetime
import numpy as np
import pickle as pickle
import pandas as pd
import netCDF4 as nc
import statsmodels.formula.api as sm
import copy
import matplotlib as mpl
from get_timeseries import get_timeseries
from plot_timeseries import timeseries
from get_averages_and_error import averages_and_error
from plot_watersheds import plot_wat
from plot_other import plot_other
from general_functions import gen_funcs
from plot_pixels import plot_pixels
from plot_timeseries import timeseries

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
        '/shared/pso/step_1x_choose_tiles_large/outputs/intersecting_catch_tiles.csv'
    )
    # where is the intersection info located?
    intersection_info_fname = (
        '/shared/pso/step_1x_choose_tiles_large/outputs/intersection_info.pkl'
    )
    # where is the le truth located?
    le_truth_fname = (
        '/shared/pso/step_3x_process_gleam/outputs/' +
        'le_truth_gleam_38a_watts_per_m2_1995-01-01_2014-12-31_' +
        'camels_tiles.csv'
    )
    # where is the streamflow truth located?
    streamflow_truth_fname = (
        '/shared/pso/step_3.1.1x_process_camels/outputs/' +
        'camels_truth_yearly_1995-01-01_2014-12-31_mm_day.csv'
    )
    # where is the tile info?
    tile_pft_info_fname = (
        '/shared/pso/step_5_analyze_outputs/outputs/pft_distribution.csv'
    )
    # where is the precip info for calc g1?
    precip_fname = (
        '/shared/pso/step_2x_env_covariates/outputs/' +
        'gldas_avg_precip_normalized.nc4'
    )
    # where is the canopy height info for g1 calc?
    canopy_fname = (
        '/shared/pso/step_2x_env_covariates/outputs/' +
        'canopy_height_normalized.nc4'
    )
    # give me an example restart file (so we can get wetness at wilting point)
    restart_fname = (
        '/lustre/catchment/exps/GEOSldas_CN45_pso_g1_a0_a1_et_strm_' +
        'camels_spin19921994_test19952014_mae/0/output/SMAP_EASEv2_M36/' +
        'rs/ens0000/Y2015/M01/0.catchcnclm45_internal_rst.20150101_0000'
    )
    # what are the names and info of the timeseries that we are going to load?
    timeseries_info = {
        'med-default-pft-g1-1992-2014':{
            'dir':(
                'nan'
            ),
            'load_or_save':'load',
            'default_type':'nan',
            'read_met_forcing':True,
            'timeseries_dir':os.path.join(
                save_dir,
                'med-default-pft-g1-1992-2014'
            ),
            'optimization_type':'nan',
            'iteration_number':'nan'
        },
        'g1-ai-et-strm-spin19921994-test19952014-mae-num8':{
            'dir':(
                '/shared/pso_outputs/g1_ai_et_strm_camels_spin19921994_test19952914_mae/num_8/29'
            ),
            'load_or_save':'load',
            'default_type':'pso_output',
            'read_met_forcing':False,
            'timeseries_dir':os.path.join(
                save_dir,
                'g1-ai-et-strm-spin19921994-test19952014-mae-num8'
            ),
            'optimization_type':'pft',
            'iteration_number':8
        },
        'g1-a0-a1-et-strm-spin19921994-test19952014-mae-num6':{
            'dir':(
                '/shared2/pso_outputs/g1_a0_a1_et_strm_camels_spin19921994_test19952914_mae/num_6/17'
            ),
            'load_or_save':'load',
            'default_type':'pso_output',
            'read_met_forcing':False,
            'timeseries_dir':os.path.join(
                save_dir,
                'g1-a0-a1-et-strm-spin19921994-test19952014-mae-num6'
            ),
            'optimization_type':'ef',
            'iteration_number':6
        }
    }
    # let's get the le truth timeseries
    gen = gen_funcs()
    g = get_timeseries()
    le_obs = g.get_le_obs(le_truth_fname)
    # let's get the streamflow truth timeseries
    strm_obs = g.get_streamflow_obs(streamflow_truth_fname)
    strm_obs_mean = strm_obs.mean()
    # and why not pixels while we are at it
    pixels = gen.get_pixels(pixels_fname)
    # load the timeseries
    exps = list(timeseries_info.keys())
    timeseries_info = g.get_all_timeseries_info(
        timeseries_info,start,end,pixels_fname,
        intersection_info_fname,precip_fname,canopy_fname,
        tile_pft_info_fname
    )
    g_a = averages_and_error()
    timeseries_info = g_a.get_all_averages_and_error(
        timeseries_info,le_obs,strm_obs,start_err,end_err
    )
    # get ef j for le and strm
    le_j_ef = (
        timeseries_info[exps[2]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
    strm_j_ef = (
        timeseries_info[exps[2]]['wat_strm_errors'].loc['mae']/
        strm_obs.mean()
    )
    # get  pft j for le and strm
    le_j_pft = (
        timeseries_info[exps[1]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
    overall_le_j_pft = np.nanmean(le_j_pft)
    strm_j_pft = (
        timeseries_info[exps[1]]['wat_strm_errors'].loc['mae']/
        strm_obs.mean()
    )
    overall_strm_j_pft = np.nanmean(strm_j_pft)
    # get differences between j for EF and PFT
    watersheds = np.array(strm_j_pft.index)
    le_j_diff_ef_pft = le_j_ef - le_j_pft
    strm_j_diff_ef_pft = strm_j_ef - strm_j_pft

    ########################################
    ####  INITIATE CLASSES FOR PLOTTING ####
    ########################################
    p_p = plot_pixels()
    #p_w = plot_wat()
    p_o = plot_other()

    ###############################################################
    ####  MAKING THE PLOTS FOR SUB-SECTION 2 OF PROPOSED PAPER ####
    ###############################################################
    # standard maps of the difference between the two methods
    #p_p.plot_map(
    #    'le_j_diff_ef_pft',pixels,le_j_diff_ef_pft*-1,
    #    np.nanmean(le_j_diff_ef_pft*-1),plots_dir,
    #    vmin=-0.2,vmax=0.2,cmap='PiYG'
    #)
    #p_w.plot_map(
    #    'strm_j_diff_ef_pft',watersheds,np.array(strm_j_diff_ef_pft)*-1,
    #    np.nanmean(np.array(strm_j_diff_ef_pft)*-1),plots_dir,
    #    cmap='PiYG',vmin=-.4,vmax=.4
    #)
    #gen.add_to_gdf_and_save(
    #    '/shared/pso/step_1x_choose_tiles_large/outputs/' +
    #    'really_chosen_tiles.geojson',
    #    os.path.join(
    #        out_dir,
    #        'le_j_diff_ef_pft.gdf'
    #    ),
    #    le_j_diff_ef_pft*-1,
    #    vmin=-0.4,
    #    vmax=0.4,
    #    cmap='PiYG'
    #)
    #gen.add_to_gdf_and_save(
    #    '/shared/pso/step_1x_choose_tiles_large/outputs/' +
    #    'chosen_camels.geojson',
    #    os.path.join(
    #        out_dir,
    #        'strm_j_diff_ef_pft.gdf'
    #    ),
    #    strm_j_diff_ef_pft*-1,
    #    vmin=-0.4,
    #    vmax=0.4,
    #    cmap='PiYG',
    #    subselection=watersheds
    #)
    # is there a correaltion between basin-wide le and improvement and
    # basin-wide streamflow improvement?
    # first we need GLEAM at watershed scale
    intersection_info = gen.get_intersection_info(
        intersection_info_fname
    )
    le_obs_wat = gen.pix_df_to_wat_df(
        le_obs,intersection_info
    )
    le_obs_wat_yr = gen.df_to_yearly(le_obs_wat)
    # the ef j at the watershed scale
    et_preds_wat_ef = timeseries_info[exps[2]]['wat_raw_timeseries']['le']
    le_preds_wat_ef = et_preds_wat_ef*28.94
    le_mae_wat_ef = g_a.mae(
        le_preds_wat_ef,
        le_obs_wat
    )
    le_j_wat_ef = le_mae_wat_ef/le_obs_wat.mean()
    le_preds_wat_yr_ef = gen.df_to_yearly(le_preds_wat_ef)
    le_mae_wat_yr_ef = g_a.mae(
        le_preds_wat_yr_ef,
        le_obs_wat_yr
    )
    le_j_wat_yr_ef = le_mae_wat_yr_ef/le_obs_wat_yr.mean()
    # the pft j at the watershed scale
    et_preds_wat_pft = timeseries_info[exps[1]]['wat_raw_timeseries']['le']
    le_preds_wat_pft = et_preds_wat_pft*28.94
    le_mae_wat_pft = g_a.mae(
        le_preds_wat_pft,
        le_obs_wat
    )
    le_j_wat_pft = le_mae_wat_pft/le_obs_wat.mean()
    le_preds_wat_yr_pft = gen.df_to_yearly(le_preds_wat_pft)
    le_mae_wat_yr_pft = g_a.mae(
        le_preds_wat_yr_pft,
        le_obs_wat_yr
    )
    le_j_wat_yr_pft = le_mae_wat_yr_pft/le_obs_wat_yr.mean()
    # difference between EF and PFT for LE at watershed scale
    le_j_wat_diff_ef_pft = le_j_wat_ef - le_j_wat_pft
    le_j_wat_yr_diff_ef_pft = le_j_wat_yr_ef - le_j_wat_yr_pft
    # for later, the stream default bias
    strm_bias_default = timeseries_info[exps[0]]['wat_strm_errors'].loc['bias']
    # scatter plot of the two
    #p_o.scatter(
    #    np.array(strm_j_diff_ef_pft),
    #    np.array(le_j_wat_diff_ef_pft),
    #    plots_dir,
    #    'strm_j_diff_ef_pft_vs_le_j_wat_diff_ef_pft',
    #    'Difference between Jstrm for PFT and EF',
    #    'Difference between JET at wat for PFT and EF',
    #    quadrant_lines=True,
    #    xlim=[-0.9,0.9],
    #    ylim=[-0.9,0.9],
    #    dot_size=2
    #)
    # what about if we assign streamflow J to each pixel instead?
    # also need to assign streamflow bias
    strm_j_pix_ef_all = [ [] for i in range(len(pixels))]
    strm_j_pix_pft_all = [ [] for i in range(len(pixels))]
    strm_bias_pix_default_all = [ [] for i in range(len(pixels))]
    strm_j_pix_weights = [ [] for i in range(len(pixels))]
    for w,wat in enumerate(watersheds):
        this_pixels = intersection_info[wat][0]
        this_perces = intersection_info[wat][1]
        this_strm_j_ef = strm_j_ef[wat]
        this_strm_j_pft = strm_j_pft[wat]
        this_strm_bias_default = strm_bias_default[wat]
        for p,pix in enumerate(this_pixels):
            this_pix_idx = np.where(pixels == pix)[0][0]
            this_perc = this_perces[p]
            strm_j_pix_ef_all[this_pix_idx].append(this_strm_j_ef)
            strm_j_pix_pft_all[this_pix_idx].append(this_strm_j_pft)
            strm_bias_pix_default_all[this_pix_idx].append(this_strm_bias_default)
            strm_j_pix_weights[this_pix_idx].append(this_perc)
    strm_j_pix_ef = np.zeros(len(pixels))
    strm_j_pix_pft = np.zeros(len(pixels))
    strm_bias_pix_default = np.zeros(len(pixels))
    for p,pix in enumerate(pixels):
        this_mean_ef = np.average(
            strm_j_pix_ef_all[p],
            weights=strm_j_pix_weights[p]
        )
        strm_j_pix_ef[p] = this_mean_ef
        this_mean_pft = np.average(
            strm_j_pix_pft_all[p],
            weights=strm_j_pix_weights[p]
        )
        strm_j_pix_pft[p] = this_mean_pft
        this_mean_bias = np.average(
            strm_bias_pix_default_all[p],
            weights=strm_j_pix_weights[p]
        )
        strm_bias_pix_default[p] = this_mean_bias
    strm_j_pix_diff_ef_pft = strm_j_pix_ef - strm_j_pix_pft
    #p_o.scatter(
    #    strm_j_pix_diff_ef_pft,
    #    np.array(le_j_diff_ef_pft),
    #    plots_dir,
    #    'strm_j_pix_diff_ef_pft_vs_le_j_diff_ef_pft',
    #    'Difference between Jstrm at pix for PFT and EF',
    #    'Difference between JET for PFT and EF',
    #    quadrant_lines=True,
    #    xlim=[-0.9,0.9],
    #    ylim=[-0.9,0.9]
    #)
    # just the pixels where there is large improvement
    improve_pix_idx = np.where(
        strm_j_pix_diff_ef_pft <= -0.3
    )
    strm_j_pix_diff_ef_pft_improve = (
        strm_j_pix_diff_ef_pft[improve_pix_idx]
    )
    improve_pix = pixels[improve_pix_idx]
    #p_p.plot_map(
    #    'strm_j_pix_diff_ef_pft_improve_pix',
    #    improve_pix,
    #    strm_j_pix_diff_ef_pft_improve,
    #    np.nanmean(strm_j_pix_diff_ef_pft_improve),
    #    plots_dir
    #)
    # what % of basins fall intoG the 4 categories of improvement?
    num_watersheds = len(watersheds)
    strm_j_diff_ef_pft_max = np.nanmax(strm_j_diff_ef_pft)
    big_improve_idx = np.where(
        np.array(strm_j_diff_ef_pft) < -strm_j_diff_ef_pft_max
    )
    big_improve_wat = watersheds[big_improve_idx]
    #print('# of basins < -{:.2f}:'.format(strm_j_diff_ef_pft_max))
    #print(len(big_improve_wat))
    #small_improve_idx = np.where(
    #    (np.array(strm_j_diff_ef_pft) >= -strm_j_diff_ef_pft_max) &
    #    (np.array(strm_j_diff_ef_pft) < 0)
    #)
    #small_improve_wat = watersheds[small_improve_idx]
    #print('# of basins >= -{:.2f} & < 0:'.format(strm_j_diff_ef_pft_max))
    #print(len(small_improve_wat))
    #small_degrade_idx = np.where(
    #    (np.array(strm_j_diff_ef_pft) <= strm_j_diff_ef_pft_max) &
    #    (np.array(strm_j_diff_ef_pft) > 0)
    #)
    #small_degrade_wat = watersheds[small_degrade_idx]
    #print('# of basins <= {:.2f} & > 0:'.format(strm_j_diff_ef_pft_max))
    #print(len(small_degrade_wat))
    #big_degrade_idx = np.where(
    #    np.array(strm_j_diff_ef_pft) > strm_j_diff_ef_pft_max
    #)
    #big_degrade_wat = watersheds[big_degrade_idx]
    #print('# of basins > {:.2f}:'.format(strm_j_diff_ef_pft_max))
    #print(len(big_degrade_wat))
    ## what % of pixels fall into the 4 categories of improvement?
    strm_j_pix_diff_ef_pft_max = np.nanmax(strm_j_pix_diff_ef_pft)
    num_pixels = len(pixels)
    big_improve_idx = np.where(
        strm_j_pix_diff_ef_pft < -strm_j_pix_diff_ef_pft_max
    )
    big_improve_pix = pixels[big_improve_idx]
    #print('# of pixels < -{:.2f}:'.format(strm_j_pix_diff_ef_pft_max))
    #print(len(big_improve_pix))
    #small_improve_idx = np.where(
    #    (strm_j_pix_diff_ef_pft >= -strm_j_pix_diff_ef_pft_max) &
    #    (strm_j_pix_diff_ef_pft < 0)
    #)
    #small_improve_pix = pixels[small_improve_idx]
    #print('# of pixels >= -{:.2f} & < 0:'.format(strm_j_pix_diff_ef_pft_max))
    #print(len(small_improve_pix))
    #small_degrade_idx = np.where(
    #    (strm_j_pix_diff_ef_pft <= strm_j_pix_diff_ef_pft_max) &
    #    (strm_j_pix_diff_ef_pft > 0)
    #)
    #small_degrade_pix = pixels[small_degrade_idx]
    #print('# of pixels <= {:.2f} & > 0:'.format(strm_j_pix_diff_ef_pft_max))
    #print(len(small_degrade_pix))
    #big_degrade_idx = np.where(
    #    strm_j_pix_diff_ef_pft > strm_j_pix_diff_ef_pft_max
    #)
    #big_degrade_pix = pixels[big_degrade_idx]
    #print('# of pixels > {:.2f}:'.format(strm_j_pix_diff_ef_pft_max))
    #print(len(big_degrade_pix))
    # what would the jdiff be if we excluded the outliars?
    strm_j_diff_ef_pft_no_outliar_idx = np.where(
        np.array(strm_j_diff_ef_pft) > -strm_j_diff_ef_pft_max
    )
    strm_j_diff_ef_pft_no_outliar = np.array(strm_j_diff_ef_pft)[
        strm_j_diff_ef_pft_no_outliar_idx
    ]
    print('Jstrm EF - Jstrm PFT:')
    print(np.nanmean(strm_j_diff_ef_pft))
    print('Jstrm no_outliar EF - Jstrm no_outliar PFT:')
    print(np.nanmean(strm_j_diff_ef_pft_no_outliar))
    strm_j_pix_diff_ef_pft_no_outliar_idx = np.where(
        strm_j_pix_diff_ef_pft > -strm_j_pix_diff_ef_pft_max
    )
    strm_j_pix_diff_ef_pft_no_outliar = strm_j_pix_diff_ef_pft[
        strm_j_pix_diff_ef_pft_no_outliar_idx
    ]
    print('Jstrm,pix EF - Jstrm,pix PFT:')
    print(np.nanmean(strm_j_pix_diff_ef_pft))
    perc_strm_j_pix_diff_ef_pft = (
        np.nanmean(strm_j_pix_diff_ef_pft)/overall_strm_j_pft*100
    )
    print(
        'This is an overall change of {:.2f}%'.format(
            perc_strm_j_pix_diff_ef_pft
        )
    )
    #print('Jstrm,pix no_outliar EF - Jstrm,pix no_outliar PFT:')
    #print(np.nanmean(strm_j_pix_diff_ef_pft_no_outliar))
    #perc_strm_j_pix_diff_ef_pft_no_outliar = (
    #    np.nanmean(strm_j_pix_diff_ef_pft_no_outliar)/overall_strm_j_pft*100
    #)
    #print(
    #    'This is an overall change of {:.2f}%'.format(
    #        perc_strm_j_pix_diff_ef_pft_no_outliar
    #    )
    #)
    le_j_diff_ef_pft_no_outliar = np.array(le_j_diff_ef_pft)[
        strm_j_pix_diff_ef_pft_no_outliar_idx
    ]
    perc_le_j_diff_ef_pft_no_outliar = (
        np.nanmean(le_j_diff_ef_pft_no_outliar)/overall_le_j_pft*100
    )
    print('Jle EF - Jle PFT:')
    print(np.nanmean(le_j_diff_ef_pft))
    perc_le_j_diff_ef_pft = (
        np.nanmean(le_j_diff_ef_pft)/overall_le_j_pft*100
    )
    print(
        'This is an overall change of {:.2f}%'.format(
            perc_le_j_diff_ef_pft
        )
    )
    #print('Jle no_outliar EF - Jle no_outliar PFT:')
    #print(np.nanmean(le_j_diff_ef_pft_no_outliar))
    #perc_le_j_diff_ef_pft_no_outliar = (
    #    np.nanmean(le_j_diff_ef_pft_no_outliar)/overall_strm_j_pft*100
    #)
    #print(
    #    'This is an overall change of {:.2f}%'.format(
    #        perc_le_j_diff_ef_pft_no_outliar
    #    )
    #)
    print('The outliar watersheds are:')
    print(big_improve_wat)
    # plot the timeseries for these watersheds
    p_t = timeseries()
    strm_obs_yr = gen.df_to_yearly(strm_obs)
    to_plot = [
        timeseries_info[exps[0]]['wat_raw_timeseries'],
        timeseries_info[exps[1]]['wat_raw_timeseries'],
        timeseries_info[exps[2]]['wat_raw_timeseries'],
        {'strm_yr':strm_obs_yr}
    ]
    names = [
        'default',
        'PFT converged',
        'EF converged',
        'CAMELS'
    ]
    sizes = np.repeat(0.5,len(names))
    colors = ['b','m','g','k']
    start_plot = datetime.date(1992,1,1)
    end_plot = datetime.date(2014,12,31)
    #p_t.plot_one_var(
    #    to_plot,names,'strm_yr','mm/day',
    #    'converged_timeseries',plots_dir,
    #    locations=big_improve_wat,
    #    start=start_plot,end=end_plot,small_preds=sizes,
    #    colors=colors
    #)
    # alright, why do these pixels/basins matter
    # let's first compare g1 values
    all_pix_colors = ['k' for i in range(len(pixels))]
    not_big_improve_pix = np.zeros(0)
    for p,pix in enumerate(pixels):
        if pix in big_improve_pix:
            all_pix_colors[p] = '#42bcf5'
        else:
            not_big_improve_pix = np.append(not_big_improve_pix,pix)
    all_wat_colors = ['k' for i in range(len(watersheds))]
    not_big_improve_wat = np.zeros(0)
    for w,wat in enumerate(watersheds):
        if wat in big_improve_wat:
            all_wat_colors[w] = '#42bcf5'
        else:
            not_big_improve_wat = np.append(not_big_improve_wat,wat)
    g1_pft_df = timeseries_info[exps[1]]['g1_map']
    g1_pft = g1_pft_df.loc['g1']
    g1_ef_df = timeseries_info[exps[2]]['g1_map']
    g1_ef = g1_ef_df.loc['g1']
    g1_diff_ef_pft_df = g1_ef_df - g1_pft_df
    g1_diff_ef_pft = g1_diff_ef_pft_df.loc['g1']
    # convert g1 difference to the watershed scale
    g1_wat_pft_df = gen.pix_df_to_wat_df(
        g1_pft_df,intersection_info
    )
    g1_wat_pft = g1_wat_pft_df.loc['g1']
    g1_wat_ef_df = gen.pix_df_to_wat_df(
        g1_ef_df,intersection_info
    )
    g1_wat_ef = g1_wat_ef_df.loc['g1']
    g1_wat_diff_ef_pft_df = g1_wat_ef_df - g1_wat_pft_df
    g1_wat_diff_ef_pft = g1_wat_ef - g1_wat_pft
    # how do these g1 values compare at the watershed scale and at the pixel
    # scale?
    #p_o.scatter(
    #    g1_pft,
    #    g1_ef,
    #    plots_dir,
    #    'g1_pft_vs_g1_ef',
    #    'g1 PFT (sqrt(kPa))',
    #    'g1 EF (sqrt(kPa))',
    #    color=all_pix_colors,
    #    xlim=[0,12],
    #    ylim=[0,12]
    #)
    #p_o.scatter(
    #    g1_wat_pft,
    #    g1_wat_ef,
    #    plots_dir,
    #    'g1_wat_pft_vs_g1_wat_ef',
    #    'g1 PFT at wat (sqrt(kPa))',
    #    'g1 EF at wat (sqrt(kPa))',
    #    color=all_wat_colors,
    #    xlim=[0,8.5],
    #    ylim=[0,8.5]
    #)
    # does combination of g1 differences and poor performance by PFTs explain
    # the 6 important basins?
    #p_o.scatter(
    #    g1_wat_diff_ef_pft,
    #    strm_j_pft,
    #    plots_dir,
    #    'g1_wat_diff_ef_pft_vs_strm_j_pft',
    #    'g1 EF - g1 PFT at wat (sqrt(kPa))',
    #    'strm_j_pft',
    #    color=all_wat_colors
    #)
    # is the overall response of streamflow to change in g1 consistant at these
    # bains?
    # fit a line to the non-special basins
    # plot the scatter of g1 diff vs. streamflow diff with that line
    

    strm_yr_pft = timeseries_info[exps[1]]['wat_raw_timeseries']['strm_yr']
    strm_avg_pft = strm_yr_pft.mean()
    strm_yr_ef = timeseries_info[exps[2]]['wat_raw_timeseries']['strm_yr']
    strm_avg_ef = strm_yr_ef.mean()
    strm_avg_diff_ef_pft = strm_avg_ef - strm_avg_pft
    strm_mae_pft = g_a.mae(strm_yr_ef,strm_obs)
    strm_mae_norm_pft = strm_mae_pft/strm_obs.mean()
    avg_strm_obs = strm_obs.mean()
    strm_avg_diff_norm_ef_pft = strm_avg_diff_ef_pft/avg_strm_obs
    g1_wat_diff_ef_pft_not_big = g1_wat_diff_ef_pft[not_big_improve_wat]
    strm_avg_diff_ef_pft_not_big = strm_avg_diff_ef_pft[not_big_improve_wat]
    #data = pd.DataFrame({
    #    'g1_diff':g1_wat_diff_ef_pft_not_big,
    #    'strm_diff':strm_avg_diff_ef_pft_not_big
    #})
    #model = sm.ols(
    #    'strm_diff ~ g1_diff',
    #    data
    #).fit()
    #coefs = model._results.params
    #b = coefs[0]
    #m = coefs[1]
    #p_o.scatter(
    #    g1_wat_diff_ef_pft,
    #    strm_avg_diff_ef_pft,
    #    plots_dir,
    #    'g1_wat_diff_ef_pft_vs_strm_avg_diff_ef_pft',
    #    'g1 EF - g1 PFT at wat (sqrt(kPa))',
    #    'avg(strm_ef) - avg(strm_pft)',
    #    quadrant_lines=True,
    #    color=all_wat_colors,
    #    line_coefficients = [b,m]
    #)
    strm_avg_diff_norm_ef_pft_not_big = strm_avg_diff_norm_ef_pft[not_big_improve_wat]
    #data = pd.DataFrame({
    #    'g1_diff':g1_wat_diff_ef_pft_not_big,
    #    'strm_diff_norm':strm_avg_diff_norm_ef_pft_not_big
    #})
    #model = sm.ols(
    #    'strm_diff_norm ~ g1_diff',
    #    data
    #).fit()
    #coefs = model._results.params
    #b = coefs[0]
    #m = coefs[1]
    #p_o.scatter(
    #    g1_wat_diff_ef_pft,
    #    strm_avg_diff_norm_ef_pft,
    #    plots_dir,
    #    'g1_wat_diff_ef_pft_vs_strm_avg_diff_norm_ef_pft',
    #    'g1 EF - g1 PFT at wat (sqrt(kPa))',
    #    '(avg(strm_ef) - avg(strm_pft))/avg(strm_obs)',
    #    quadrant_lines=True,
    #    color=all_wat_colors,
    #    line_coefficients = [b,m]
    #)
    # okay, let's define sensitivity
    # first we only want to do this for basins wehre the g1 change >1 to avoid
    # our value from really blowing up at small g1 changes
    g1_wat_diff_abs_ef_pft = g1_wat_diff_ef_pft.abs()
    big_g1_diff_idx = g1_wat_diff_ef_pft[
        g1_wat_diff_abs_ef_pft > 0.5
    ].index
    big_g1_diff_idx_num = np.where(g1_wat_diff_abs_ef_pft > 0.5)[0]
    big_g1_diff_idx_int = [int(i) for i in big_g1_diff_idx_num]
    # get these actual g1 values and define our colors
    g1_wat_diff_big_ef_pft = g1_wat_diff_ef_pft[big_g1_diff_idx]
    all_wat_colors_big = np.array(all_wat_colors)[big_g1_diff_idx_int]
    # define our metric at these basins
    # first we need the g1 differences
    g1_wat_big_diff_ef_pft = g1_wat_diff_ef_pft[big_g1_diff_idx]
    strm_avg_diff_norm_big_ef_pft = strm_avg_diff_norm_ef_pft[
        big_g1_diff_idx
    ]
    #now compute our metric
    strm_sens_to_g1 = strm_avg_diff_norm_big_ef_pft/g1_wat_big_diff_ef_pft
    # what if we use a different sensitivity metric
    #dist_from_exp_sens = (
    #    (
    #        -m*g1_wat_diff_ef_pft +
    #        strm_avg_diff_norm_ef_pft -
    #        b
    #    )/(
    #        np.sqrt((m)**2)
    #    )
    #)
    #print(dist_from_exp_sens)
    #p_o.scatter(
    #    g1_wat_diff_ef_pft,
    #    strm_avg_diff_norm_ef_pft,
    #    plots_dir,
    #    'g1_wat_diff_ef_pft_vs_strm_avg_diff_norm_ef_pft',
    #    'g1 EF - g1 PFT at wat (sqrt(kPa))',
    #    '(avg(strm_ef) - avg(strm_pft))/avg(strm_obs)',
    #    quadrant_lines=True,
    #    color=np.array(dist_from_exp_sens),
    #    line_coefficients = [b,m]
    #)
    # let's see if we can exlpain this metric
    #start with mean streamflow
    strm_obs_mean_big = strm_obs_mean[big_g1_diff_idx]
    #p_o.scatter(
    #    strm_obs_mean_big,
    #    strm_sens_to_g1,
    #    plots_dir,
    #    'strm_obs_mean_big_vs_strm_sens_to_g1',
    #    'avg(strm_obs) where abs_g1_diff > 1',
    #    'streamflow sensitivity',
    #    color=all_wat_colors_big,
    #    x_log=True
    #)
    #p_o.scatter(
    #    strm_obs_mean,
    #    dist_from_exp_sens,
    #    plots_dir,
    #    'strm_obs_mean_vs_dist_from_exp_sens',
    #    'avg(strm_obs) where abs_g1_diff > 1',
    #    'distance from expected sensitivity',
    #    color=all_wat_colors
    #)
    # we also want to test if it matter the degree to which the basin is
    # homogenous
    # create a dataframe where each column is a basin and each row is the % of
    # each of the 5 pfts in that basin
    tile_pft_info = pd.read_csv(tile_pft_info_fname)
    tile_pft_info = tile_pft_info.set_index('tile')
    simple_map = {
        'Needleleaf evergreen temperate tree':'needleleaf_trees',
        'Needleleaf evergreen boreal tree':'needleleaf_trees',
        'Needleleaf deciduous boreal tree':'broadleaf_trees',
        'Broadleaf evergreen tropical tree':'broaleaf_trees',
        'Broadleaf evergreen temperate tree':'broadleaf_trees',
        'Broadleaf deciduous tropical tree':'broadleaf_trees',
        'Broadleaf deciduous temperate tree':'broadleaf_trees',
        'Broadleaf deciduous boreal tree':'broadleaf_trees',
        'Broadleaf evergreen temperate shrub':'shrub',
        'Broadleaf deciduous boreal shrub':'shrub',
        'Broadleaf deciduous temperate shrub':'shrub',
        'Broadleaf deciduous temperate shrub[moisture +deciduous]': (
            'shrub'
        ),
        'Broadleaf deciduous temperate shrub[moisture stress only]': (
            'shrub'
        ),
        'Arctic c3 grass':'c3_grass',
        'Cool c3 grass':'c3_grass',
        'Cool c3 grass [moisture + deciduous]':'c3_grass',
        'Cool c3 grass [moisture stress only]':'c3_grass',
        'Warm c4 grass [moisture + deciduous]':'c4_grass',
        'Warm c4 grass [moisture stress only]':'c4_grass',
        'Warm c4 grass':'c4_grass',
        'Crop':'crop',
        'Crop [moisture + deciduous]':'crop',
        'Crop [moisture stress only]':'crop',
        '(Spring temperate cereal)':'crop',
        '(Irrigated corn)':'crop',
        '(Soybean)':'crop',
        '(Corn)':'crop'
    }
    # make a df where each column is a tile and the row is the pft %
    tile_pft_info_df = pd.DataFrame(columns=pixels)
    tile_pft_info_df['pfts'] = [
        'needleleaf_trees','broadleaf_trees','shrub','c3_grass','c4_grass',
        'crop'
    ]
    tile_pft_info_df = tile_pft_info_df.set_index('pfts')
    for p,pix in enumerate(pixels):
        for pft in range(4):
            p = pft + 1
            this_pft = tile_pft_info['pft_{}_name'.format(p)].loc[pix]
            this_pft_simple = simple_map[this_pft]
            this_perc = tile_pft_info['pft_{}_perc'.format(p)].loc[pix]
            if pd.isnull(tile_pft_info_df[pix].loc[this_pft_simple]) == True:
                tile_pft_info_df[pix].loc[this_pft_simple] = this_perc
            else:
                tile_pft_info_df[pix].loc[this_pft_simple] += this_perc
    # now we can convert this to the watershed space
    tile_pft_info_wat_df = gen.pix_df_to_wat_df(
        tile_pft_info_df,intersection_info
    )
    tile_pft_info_wat_df = tile_pft_info_wat_df/tile_pft_info_wat_df.sum()
    tile_pft_info_wat_df = tile_pft_info_wat_df*100
    # get the % of the largest PFT in each pixel
    largest_pft_perc = tile_pft_info_wat_df.max()
    largest_pft_perc_big = largest_pft_perc[big_g1_diff_idx]
    #p_o.scatter(
    #    largest_pft_perc,
    #    dist_from_exp_sens,
    #    plots_dir,
    #    'largest_pft_perc_vs_dist_from_exp_sens',
    #    '% of basin covered by most abundant PFT',
    #    'distance from expected sensitivity',
    #    color=all_wat_colors
    #)
    #p_o.scatter(
    #    largest_pft_perc_big,
    #    strm_sens_to_g1,
    #    plots_dir,
    #    'largest_pft_perc_big_vs_strm_sens_to_g1',
    #    '% of basin covered by most abundant PFT',
    #    'strm sensitivity to g1 change',
    #    color=all_wat_colors_big
    #)
    # we also want to test if it matters how homogenous g1 is within a basin
    # calculate the std of the g1 values in a given basin
    # let's get a weighted standard deviation though so that we value the
    # percent of each pixel in the basin
    # we can do this for pft, ef, and then the difference in the std
    g1_wat_ef_std = np.zeros(len(watersheds))
    g1_wat_pft_std = np.zeros(len(watersheds))
    for w,wat in enumerate(watersheds):
        # get the g1 values that contribute
        this_pixels = intersection_info[wat][0]
        this_perc = intersection_info[wat][1]
        # for pfts
        this_g1_pft = g1_pft[this_pixels]
        this_g1_pft_avg = np.average(this_g1_pft,weights=this_perc)
        this_g1_pft_var = np.average(
            (this_g1_pft - this_g1_pft_avg)**2,
            weights=this_perc
        )
        this_g1_pft_std = np.sqrt(this_g1_pft_var)
        g1_wat_pft_std[w] = this_g1_pft_std
        # for ef
        this_g1_ef = g1_ef[this_pixels]
        this_g1_ef_avg = np.average(this_g1_ef,weights=this_perc)
        this_g1_ef_var = np.average(
            (this_g1_ef - this_g1_ef_avg)**2,
            weights=this_perc
        )
        this_g1_ef_std = np.sqrt(this_g1_ef_var)
        g1_wat_ef_std[w] = this_g1_ef_std
    g1_wat_pft_df.loc['std'] = g1_wat_pft_std
    g1_wat_ef_df.loc['std'] = g1_wat_ef_std
    g1_wat_pft_df.loc['std_norm'] = g1_wat_pft_std/g1_wat_pft
    g1_wat_ef_df.loc['std_norm'] = g1_wat_ef_std/g1_wat_ef
    g1_wat_pft_std_big = g1_wat_pft_df[big_g1_diff_idx].loc['std']
    g1_wat_ef_std_big = g1_wat_ef_df[big_g1_diff_idx].loc['std']
    g1_wat_pft_std_norm_big = g1_wat_pft_df[big_g1_diff_idx].loc['std_norm']
    g1_wat_ef_std_norm_big = g1_wat_ef_df[big_g1_diff_idx].loc['std_norm']
    #p_o.scatter(
    #    g1_wat_pft_std_norm_big,
    #    strm_sens_to_g1,
    #    plots_dir,
    #    'g1_wat_pft_std_norm_big_vs_strm_sens_to_g1',
    #    'std(g1 PFT)/avg(g1 PFT)',
    #    'strm sensitivity to g1',
    #    color=all_wat_colors_big
    #)
    #p_o.scatter(
    #    g1_wat_ef_std,
    #    dist_from_exp_sens,
    #    plots_dir,
    #    'g1_wat_ef_std_vs_dist_from_exp_sens',
    #    'std(g1 EF)',
    #    'distance from expected sensitivity',
    #    color=all_wat_colors
    #)
    #p_o.scatter(
    #    g1_wat_ef_std_norm_big,
    #    strm_sens_to_g1,
    #    plots_dir,
    #    'g1_wat_ef_std_norm_big_vs_strm_sens_to_g1',
    #    'std(g1 EF)/avg(g1 EF)',
    #    'strm sensitivity to g1',
    #    color=all_wat_colors_big
    #)
    # let's check with precip as well
    precip = np.array(
        nc.Dataset(
            precip_fname
        )['vals']
    )
    precip_df = pd.DataFrame(columns=pixels)
    precip_df.loc['precip'] = precip
    precip_wat_df = gen.pix_df_to_wat_df(
        precip_df,intersection_info
    )
    precip_wat = precip_wat_df.loc['precip']
    precip_wat_big = precip_wat[big_g1_diff_idx]
    # load canopy height while we are at it
    canopy = np.array(
        nc.Dataset(
            canopy_fname
        )['vals']
    )
    canopy_df = pd.DataFrame(columns=pixels)
    canopy_df.loc['canopy'] = canopy
    canopy_wat_df = gen.pix_df_to_wat_df(
        canopy_df,intersection_info
    )
    canopy_wat = canopy_wat_df.loc['canopy']
    #p_o.scatter(
    #    precip_wat,
    #    dist_from_exp_sens,
    #    plots_dir,
    #    'precip_wat_vs_dist_from_exp_sens',
    #    'mean annual precipitation (mm/day)',
    #    'distance from expected sensitivity',
    #    color=all_wat_colors
    #)
    #p_o.scatter(
    #    precip_wat_big,
    #    strm_sens_to_g1,
    #    plots_dir,
    #    'precip_wat_vs_strm_sens_to_g1',
    #    'mean annual precipitation (mm/day)',
    #    'strm sensitivity to g1',
    #    color=all_wat_colors_big
    #)
    # get the basins where each PFT is > x%
    pft_thresh_perc = 50
    pfts = list(tile_pft_info_wat_df.index)
    pfts_plus = copy.deepcopy(pfts)
    pfts_plus.append('no_majority')
    pft_wats = {}
    pft_pix = {}
    pft_strm_j_pft_mean = [0 for k in range(len(pfts_plus))]
    pft_strm_j_ef_mean = [0 for k in range(len(pfts_plus))]
    pft_strm_j_diff_ef_pft_mean = [0 for k in range(len(pfts_plus))]
    pft_strm_j_pft_std = [0 for k in range(len(pfts_plus))]
    pft_strm_j_ef_std = [0 for k in range(len(pfts_plus))]
    pft_strm_j_diff_ef_pft_std = [0 for k in range(len(pfts_plus))]
    pft_num_wats = [0 for k in range(len(pfts_plus))]
    wats_used = []
    for p,pft in enumerate(pfts_plus):
        if pft == 'no_majority':
            this_pft_wats = [
                w for w in watersheds if w not in wats_used
            ]
        else:
            this_pft_wats = tile_pft_info_wat_df.loc[pft][
                tile_pft_info_wat_df.loc[pft] > pft_thresh_perc
            ].index
        this_pft_wats = list(this_pft_wats)
        wats_used.extend(this_pft_wats)
        this_num_wats = len(this_pft_wats)
        if this_num_wats == 0:
            pft_wats[pft] = 'none'
            continue
        this_pft_wats_df = pd.DataFrame(columns=this_pft_wats)
        this_pft_wats_df.loc['strm_j_pft'] = strm_j_pft[this_pft_wats]
        this_pft_wats_df.loc['strm_j_ef'] = strm_j_ef[this_pft_wats]
        this_pft_wats_df.loc['strm_j_diff_ef_pft'] = (
            strm_j_diff_ef_pft[this_pft_wats]
        )
        pft_wats[pft] = this_pft_wats_df
        pft_strm_j_pft_mean[p] = (
            this_pft_wats_df.loc['strm_j_pft'].mean()
        )
        pft_strm_j_ef_mean[p] = (
            this_pft_wats_df.loc['strm_j_ef'].mean()
        )
        pft_strm_j_diff_ef_pft_mean[p] = (
            this_pft_wats_df.loc['strm_j_diff_ef_pft'].mean()
        )
        pft_strm_j_pft_std[p] = (
            this_pft_wats_df.loc['strm_j_pft'].std()
        )
        pft_strm_j_ef_std[p] = (
            this_pft_wats_df.loc['strm_j_ef'].std()
        )
        pft_strm_j_diff_ef_pft_std[p] = (
            this_pft_wats_df.loc['strm_j_diff_ef_pft'].std()
        )
        pft_num_wats[p] = this_num_wats
    pft_strm_summary = pd.DataFrame(columns=pfts_plus)
    pft_strm_summary.loc['strm_j_pft_mean'] = pft_strm_j_pft_mean
    pft_strm_summary.loc['strm_j_ef_mean'] = pft_strm_j_ef_mean
    pft_strm_summary.loc['strm_j_diff_ef_pft_mean'] = (
        pft_strm_j_diff_ef_pft_mean
    )
    pft_strm_summary.loc['strm_j_pft_std'] = pft_strm_j_pft_std
    pft_strm_summary.loc['strm_j_ef_std'] = pft_strm_j_ef_std
    pft_strm_summary.loc['strm_j_diff_ef_pft_std'] = (
        pft_strm_j_diff_ef_pft_std
    )
    pft_strm_summary.loc['num_wats'] = pft_num_wats
    # let's check the same thing for le
    pft_le_j_pft_mean = [0 for k in range(len(pfts_plus))]
    pft_le_j_ef_mean = [0 for k in range(len(pfts_plus))]
    pft_le_j_diff_ef_pft_mean = [0 for k in range(len(pfts_plus))]
    pft_le_j_pft_std = [0 for k in range(len(pfts_plus))]
    pft_le_j_ef_std = [0 for k in range(len(pfts_plus))]
    pft_le_j_diff_ef_pft_std = [0 for k in range(len(pfts_plus))]
    pft_num_pix = [0 for k in range(len(pfts_plus))]
    pft_g1_std = [0 for k in range(len(pfts_plus))]
    ef_g1_std = [0 for k in range(len(pfts_plus))]
    pft_precip_std = [0 for k in range(len(pfts_plus))]
    pft_canopy_std = [0 for k in range(len(pfts_plus))]
    pft_g1_median = [0 for k in range(len(pfts_plus))]
    ef_g1_median = [0 for k in range(len(pfts_plus))]
    pix_used = []
    for p,pft in enumerate(pfts_plus):
        if pft == 'no_majority':
            this_pft_pix = [
                pix for pix in pixels if pix not in pix_used
            ]
        else:
            this_pft_pix = tile_pft_info_df.loc[pft][
                tile_pft_info_df.loc[pft] > pft_thresh_perc
            ].index
        this_pft_pix = list(this_pft_pix)
        pix_used.extend(this_pft_pix)
        this_num_pix = len(this_pft_pix)
        if this_num_pix == 0:
            pft_pix[pft] = 'none'
            continue
        this_pft_pix_df = pd.DataFrame(columns=this_pft_pix)
        this_pft_pix_df.loc['le_j_pft'] = le_j_pft[this_pft_pix]
        this_pft_pix_df.loc['le_j_ef'] = le_j_ef[this_pft_pix]
        this_pft_pix_df.loc['le_j_diff_ef_pft'] = (
            le_j_diff_ef_pft[this_pft_pix]
        )
        this_pft_pix_df.loc['pft_g1'] = g1_pft[this_pft_pix]
        this_pft_pix_df.loc['ef_g1'] = g1_ef[this_pft_pix]
        this_pft_pix_df.loc['precip'] = precip_df[this_pft_pix].loc['precip']
        this_pft_pix_df.loc['canopy'] = canopy_df[this_pft_pix].loc['canopy']
        pft_pix[pft] = this_pft_pix_df
        pft_le_j_pft_mean[p] = (
            this_pft_pix_df.loc['le_j_pft'].mean()
        )
        pft_le_j_ef_mean[p] = (
            this_pft_pix_df.loc['le_j_ef'].mean()
        )
        pft_le_j_diff_ef_pft_mean[p] = (
            this_pft_pix_df.loc['le_j_diff_ef_pft'].mean()
        )
        pft_le_j_pft_std[p] = (
            this_pft_pix_df.loc['le_j_pft'].std()
        )
        pft_le_j_ef_std[p] = (
            this_pft_pix_df.loc['le_j_ef'].std()
        )
        pft_le_j_diff_ef_pft_std[p] = (
            this_pft_pix_df.loc['le_j_diff_ef_pft'].std()
        )
        pft_num_pix[p] = this_num_pix
        pft_g1_std[p] = (
            this_pft_pix_df.loc['pft_g1'].std()
        )
        ef_g1_std[p] = (
            this_pft_pix_df.loc['ef_g1'].std()
        )
        pft_precip_std[p] = (
            this_pft_pix_df.loc['precip'].std()
        )
        pft_canopy_std[p] = (
            this_pft_pix_df.loc['canopy'].std()
        )
        pft_g1_median[p] = (
            this_pft_pix_df.loc['pft_g1'].median()
        )
        ef_g1_median[p] = (
            this_pft_pix_df.loc['ef_g1'].median()
        )
    pft_le_summary = pd.DataFrame(columns=pfts_plus)
    pft_le_summary.loc['le_j_pft_mean'] = pft_le_j_pft_mean
    pft_le_summary.loc['le_j_ef_mean'] = pft_le_j_ef_mean
    pft_le_summary.loc['le_j_diff_ef_pft_mean'] = (
        pft_le_j_diff_ef_pft_mean
    )
    pft_le_summary.loc['le_j_pft_std'] = pft_le_j_pft_std
    pft_le_summary.loc['le_j_ef_std'] = pft_le_j_ef_std
    pft_le_summary.loc['le_j_diff_ef_pft_std'] = (
        pft_le_j_diff_ef_pft_std
    )
    pft_le_summary.loc['num_pix'] = pft_num_pix
    pft_le_summary.loc['pft_g1_std'] = pft_g1_std
    pft_le_summary.loc['ef_g1_std'] = ef_g1_std
    pft_le_summary.loc['precip_std'] = pft_precip_std
    pft_le_summary.loc['canopy_std'] = pft_canopy_std
    pft_le_summary.loc['pft_g1_median'] = pft_g1_median
    pft_le_summary.loc['ef_g1_median'] = ef_g1_median
    # bar plot of all of those values
    vals = np.array(
        [
            pft_strm_summary.loc['strm_j_pft_mean'],
            pft_strm_summary.loc['strm_j_ef_mean']
        ]
    )
    categories = list(pft_strm_summary.columns)
    annotation = [0 for k in range(len(pfts_plus))]
    for p,pft in enumerate(pfts_plus):
        annotation[p] = 'n={}'.format(
            pft_strm_summary[pft].loc['num_wats']
        )
    p_o.multibar_plot(
        vals,categories,['PFT','EF'],plots_dir,
        'bar_strm_j_pft_mean_vs_strm_j_ef_mean',
        'PFTs','Jstrm',annotation=annotation
    )
    # plot the Jstrm for each of those basins
    vals = np.array(
        [
            pft_strm_summary.loc['strm_j_pft_std'],
            pft_strm_summary.loc['strm_j_ef_std']
        ]
    )
    categories = list(pft_strm_summary.columns)
    annotation = [0 for k in range(len(pfts_plus))]
    for p,pft in enumerate(pfts_plus):
        annotation[p] = 'n={}'.format(
            np.round(pft_strm_summary[pft].loc['num_wats'])
        )
    p_o.multibar_plot(
        vals,categories,['PFT','EF'],plots_dir,
        'bar_strm_j_pft_std_vs_strm_j_ef_std',
        'PFTs','std(Jstrm)',annotation=annotation
    )
    # plot the Jstrm for each of those basins
    vals = np.array(
        [
            pft_le_summary.loc['pft_g1_std'],
            pft_le_summary.loc['ef_g1_std']
        ]
    )
    categories = list(pft_strm_summary.columns)
    annotation = [0 for k in range(len(pfts_plus))]
    for p,pft in enumerate(pfts_plus):
        annotation[p] = 'n={}'.format(
            np.round(pft_le_summary[pft].loc['num_pix'])
        )
    p_o.multibar_plot(
        vals,categories,['PFT','EF'],plots_dir,
        'bar_pft_g1_std_vs_ef_g1_std',
        'PFTs','std(g1)',annotation=annotation
    )
    # for the environmental predictors
    vals = np.array(
        [
            pft_le_summary.loc['precip_std'],
            pft_le_summary.loc['canopy_std']
        ]
    )
    categories = list(pft_strm_summary.columns)
    annotation = [0 for k in range(len(pfts_plus))]
    for p,pft in enumerate(pfts_plus):
        annotation[p] = 'n={}'.format(
            np.round(pft_le_summary[pft].loc['num_pix'])
        )
    p_o.multibar_plot(
        vals,categories,['precip','canopy'],plots_dir,
        'bar_precip_std_vs_canopy_std',
        'PFTs','std(norm(env)',annotation=annotation
    )
    # for the medians
    vals = np.array(
        [
            pft_le_summary.loc['pft_g1_median'],
            pft_le_summary.loc['ef_g1_median']
        ]
    )
    categories = list(pft_strm_summary.columns)
    annotation = [0 for k in range(len(pfts_plus))]
    for p,pft in enumerate(pfts_plus):
        annotation[p] = 'n={}'.format(
            np.round(pft_le_summary[pft].loc['num_pix'])
        )
    p_o.multibar_plot(
        vals,categories,['PFT','EF'],plots_dir,
        'bar_pft_g1_median_vs_ef_g1_median',
        'PFTs','median(g1)',annotation=annotation
    )
    # histogram of Jstrm for EF and for PFT
    bins = np.linspace(0,2,21)
    p_o.double_histogram(
        strm_j_pft,strm_j_ef,'PFT','EF',bins,
        plots_dir,'hist_strm_j_pft_vs_strm_j_ef',
        x_label='Jstrm',y_label='Count'
    )
    # let's make the g1 map for PFTs and EF
    #p_p.plot_map(
    #    'g1_pft_pixel_map',
    #    pixels,
    #    g1_pft,
    #    np.nanmean(g1_pft),
    #    plots_dir,
    #    vmin=0.5,
    #    vmax=6
    #)
    #p_p.plot_map(
    #    'g1_ef_pixel_map',
    #    pixels,
    #    g1_ef,
    #    np.nanmean(g1_ef),
    #    plots_dir,
    #    vmin=0.5,
    #    vmax=6
    #)
    # we need map of the basins assigned to different land cover types
    wat_pft_cats = []
    wats_reordered = []
    for p,pft in enumerate(list(pft_wats.keys())):
        this_wat = pft_wats[pft]
        if type(this_wat) == str:
            continue
        wats_reordered.extend(this_wat.columns)
        this_pft_rep = [pft for n in range(len(this_wat.columns))]
        wat_pft_cats.extend(this_pft_rep)
    #p_w = plot_wat()
    #p_w.plot_map_categorical(
    #    'pft_classification',
    #    wats_reordered,
    #    wat_pft_cats,
    #    plots_dir
    #)
    print(len(g1_wat_ef))
    p_o.scatter(
        g1_wat_ef,
        strm_j_diff_ef_pft,
        plots_dir,
        'g1_ef_vs_strm_j_diff_ef_pft',
        'g1 from EF',
        'Jstrm EF - Jstrm PFT',
        color=all_wat_colors
    )
    p_o.scatter(
        g1_wat_diff_ef_pft,
        strm_j_diff_ef_pft,
        plots_dir,
        'g1_wat_diff_ef_pft_vs_strm_j_diff_ef_pft',
        'EF g1 - PFT g1',
        'Jstrm EF - Jstrm PFT',
        color=all_wat_colors
    )
    p_o.scatter(
        strm_j_pft,
        strm_j_diff_ef_pft,
        plots_dir,
        'strm_j_pft_vs_strm_j_diff_ef_pft',
        'Jstrm PFT',
        'Jstrm EF - Jstrm PFT',
        color=all_wat_colors
    )

    sys.exit()
















    strm_yr_pft = timeseries_info[exps[1]]['wat_raw_timeseries']['strm_yr']
    strm_avg_pft = strm_yr_pft.mean()
    strm_yr_ef = timeseries_info[exps[2]]['wat_raw_timeseries']['strm_yr']
    strm_avg_ef = strm_yr_ef.mean()
    strm_avg_diff_ef_pft = strm_avg_ef - strm_avg_pft
    strm_mae_pft = g_a.mae(strm_yr_ef,strm_obs)
    strm_mae_norm_pft = strm_mae_pft/strm_obs.mean()
    avg_strm_obs = strm_obs.mean()
    strm_avg_diff_norm_ef_pft = strm_avg_diff_ef_pft/avg_strm_obs
    p_o.scatter(
        g1_wat_diff_ef_pft_df,
        strm_avg_diff_norm_ef_pft,
        plots_dir,
        'g1_wat_diff_ef_pft_vs_strm_avg_diff_norm_ef_pft',
        'g1 EF - g1 PFT (sqrt(kPa))',
        '(strm_yr EF - strm_yr PFT)/avg(strm_obs) (-)',
        quadrant_lines=True,
        color=all_wat_colors
    )
    p_o.scatter(
        g1_wat_diff_ef_pft_df,
        avg_strm_obs,
        plots_dir,
        'g1_wat_diff_ef_pft_vs_strm_obs_avg',
        'g1 EF - g1 PFT (sqrt(kPa))',
        'avg(strm_obs) (mm/day)',
        quadrant_lines=True,
        color=all_wat_colors
    )
    p_o.scatter(
        g1_wat_diff_ef_pft_df,
        strm_mae_pft,
        plots_dir,
        'g1_wat_diff_ef_pft_vs_strm_mae_pft',
        'g1 EF - g1 PFT (sqrt(kPa))',
        'strm MAE for PFT (mm/day)',
        quadrant_lines=True,
        color=all_wat_colors
    )
    p_o.scatter(
        g1_wat_diff_ef_pft_df,
        strm_mae_norm_pft,
        plots_dir,
        'g1_wat_diff_ef_pft_vs_strm_mae_norm_pft',
        'g1 EF - g1 PFT (sqrt(kPa))',
        'normalized strm MAE for PFT (mm/day)',
        quadrant_lines=True,
        color=all_wat_colors
    )
    # alright; now let's just isolate the basins that fit our criteria for
    # potentially mattering a lot for EF and PFT Jstrm differences. that is
    # pixels that have g1 difference between two methods of > 2 and streamflow
    # < 1.5 mm/day
    g1_wat_diff_ef_pft = np.array(g1_wat_diff_ef_pft_df.loc['g1'])
    avg_strm_obs = np.array(avg_strm_obs)
    strm_j_diff_ef_pft = np.array(strm_j_diff_ef_pft)
    matters_idx = np.where(
        (g1_wat_diff_ef_pft >= 2) &
        (avg_strm_obs <= 1.5)
    )
    matters_wat = watersheds[matters_idx]
    matters_colors = ['#42bcf5' for i in range(len(matters_wat))]
    for w,wat in enumerate(matters_wat):
        if wat in big_improve_wat:
            matters_colors[w] = '#a65db7'
    g1_wat_diff_ef_pft_matters = g1_wat_diff_ef_pft[matters_idx]
    strm_j_diff_ef_pft_matters = strm_j_diff_ef_pft[matters_idx]
    p_o.scatter(
        g1_wat_diff_ef_pft_matters,
        strm_j_diff_ef_pft_matters,
        plots_dir,
        'g1_wat_diff_ef_pft_matters_vs_strm_j_diff_ef_pft_matters',
        'g1 EF - g1 PFT (sqrt(kPa))',
        'Jstrm EF - Jstrm PFT',
        quadrant_lines=True,
        color=matters_colors
    )
    # expand the timeseries plotting to be for all of these
    #p_t.plot_one_var(
    #    to_plot,names,'strm_yr','mm/day',
    #    'converged_timeseries',plots_dir,
    #    locations=matters_wat,
    #    start=start_plot,end=end_plot,small_preds=sizes,
    #    colors=colors
    #)
    # is it a function of how much g1 difference actaully effects streamflow
    # differences in these basins?
    strm_avg_norm_diff_ef_pft_matters = np.zeros(len(matters_wat))
    pft_bias_matters = np.zeros(len(matters_wat))
    ef_bias_matters = np.zeros(len(matters_wat))
    pft_bias_norm_matters = np.zeros(len(matters_wat))
    ef_bias_norm_matters = np.zeros(len(matters_wat))
    strm_yr_pft = timeseries_info[exps[1]]['wat_raw_timeseries']['strm_yr']
    strm_yr_ef = timeseries_info[exps[2]]['wat_raw_timeseries']['strm_yr']
    for w,wat in enumerate(matters_wat):
        this_ef = strm_yr_ef[wat]
        this_pft = strm_yr_pft[wat]
        this_diff = this_ef - this_pft
        this_diff_avg = this_diff.mean()
        this_strm_obs = strm_obs[wat]
        this_strm_obs_avg = this_strm_obs.mean()
        this_diff_avg_norm = this_diff_avg/this_strm_obs_avg
        strm_avg_norm_diff_ef_pft_matters[w] = this_diff_avg_norm
        this_pft_bias = g_a.bias(this_pft,this_strm_obs)
        this_ef_bias = g_a.bias(this_ef,this_strm_obs)
        pft_bias_matters[w] = this_pft_bias
        ef_bias_matters[w] = this_ef_bias
        this_pft_bias_norm = this_pft_bias/this_strm_obs_avg
        this_ef_bias_norm = this_ef_bias/this_strm_obs_avg
        pft_bias_norm_matters[w] = this_pft_bias_norm
        ef_bias_norm_matters[w] = this_ef_bias_norm
    p_o.scatter(
        g1_wat_diff_ef_pft_matters,
        strm_avg_norm_diff_ef_pft_matters,
        plots_dir,
        'g1_wat_diff_ef_pft_matters_vs_strm_avg_norm_diff_ef_pft_matters',
        'EF g1 - PFT g1 (sqrt(kPa))',
        'avg(EF_strm - PFT_strm)/avg(strm_obs)',
        color=matters_colors
    )
    p_o.scatter(
        pft_bias_matters,
        strm_j_diff_ef_pft_matters,
        plots_dir,
        'pft_bias_matters_vs_strm_j_diff_ef_pft_matters',
        'PFT strm bias (mm/day)',
        'Jstrm EF - Jstrm PFT',
        quadrant_lines=True,
        color=matters_colors
    )
    p_o.scatter(
        ef_bias_matters,
        strm_j_diff_ef_pft_matters,
        plots_dir,
        'ef_bias_matters_vs_strm_j_diff_ef_pft_matters',
        'PFT strm bias (mm/day)',
        'Jstrm EF - Jstrm PFT',
        quadrant_lines=True,
        color=matters_colors
    )
    p_o.scatter(
        pft_bias_matters,
        ef_bias_matters,
        plots_dir,
        'pft_bias_matters_vs_ef_bias_matters',
        'PFT strm bias (mm/day)',
        'EF strm bias (mm/day)',
        quadrant_lines=True,
        color=matters_colors
    )
    p_o.scatter(
        pft_bias_norm_matters,
        ef_bias_norm_matters,
        plots_dir,
        'pft_bias_norm_matters_vs_ef_bias_norm_matters',
        'PFT strm bias normalized',
        'EF strm bias normalized',
        quadrant_lines=True,
        color=matters_colors,
        one_to_one_line=True
    )
    # for completeness, let's take a look at the g1 maps
    #p_p.plot_map(
    #    'g1_pft_trimmed',
    #    pixels,
    #    g1_pft,
    #    np.nanmean(g1_pft),
    #    plots_dir,
    #    vmin=0.5,
    #    vmax=4
    #)
    #p_p.plot_map(
    #    'g1_ef_trimmed',
    #    pixels,
    #    g1_ef,
    #    np.nanmean(g1_ef),
    #    plots_dir,
    #    vmin=0.5,
    #    vmax=4
    #)
    #p_p.plot_map(
    #    'g1_pft_full',
    #    pixels,
    #    g1_pft,
    #    np.nanmean(g1_pft),
    #    plots_dir,
    #    vmin=0.5,
    #    vmax=12
    #)
    #p_p.plot_map(
    #    'g1_ef_full',
    #    pixels,
    #    g1_ef,
    #    np.nanmean(g1_ef),
    #    plots_dir,
    #    vmin=0.5,
    #    vmax=12
    #)
    #p_p.plot_map(
    #    'g1_diff_ef_pft',
    #    pixels,
    #    g1_diff_ef_pft,
    #    np.nanmean(g1_diff_ef_pft),
    #    plots_dir,
    #    vmin=-5,
    #    vmax=5,
    #    cmap='PiYG'
    #)
    p_o.scatter(
        g1_wat_diff_ef_pft,
        strm_avg_diff_norm_ef_pft,
        plots_dir,
        'g1_wat_diff_ef_pft_vs_strm_avg_diff_norm_ef_pft_best_fit',
        'g1 EF - g1 PFT (sqrt(kPa))',
        '(strm_yr EF - strm_yr PFT)/avg(strm_obs) (-)',
        quadrant_lines=True,
        best_fit_line=True
    )
    le_pft = timeseries_info[exps[1]]['pixel_raw_timeseries']['le']
    le_ef = timeseries_info[exps[2]]['pixel_raw_timeseries']['le']
    le_diff_ef_pft = le_ef - le_pft
    le_avg_diff_ef_pft = le_diff_ef_pft.mean()
    le_obs_avg = le_obs.mean()
    le_avg_diff_norm_ef_pft = le_avg_diff_ef_pft/le_obs_avg
    p_o.scatter(
        g1_diff_ef_pft,
        le_avg_diff_norm_ef_pft,
        plots_dir,
        'g1_diff_ef_pft_vs_le_avg_diff_norm_ef_pft_best_fit',
        'g1 EF - g1 PFT (sqrt(kPa))',
        '(LE EF - LE PFT)/avg(le_obs) (-)',
        quadrant_lines=True,
        best_fit_line=True
    )
    le_wat_pft = gen.pix_df_to_wat_df(
        le_pft,intersection_info
    )
    le_wat_ef = gen.pix_df_to_wat_df(
        le_ef,intersection_info
    )
    le_wat_diff_ef_pft = le_wat_ef - le_wat_pft
    le_wat_avg_diff_ef_pft = le_wat_diff_ef_pft.mean()
    le_wat_avg_diff_norm_ef_pft = le_wat_avg_diff_ef_pft/le_obs_wat.mean()
    p_o.scatter(
        strm_avg_diff_norm_ef_pft,
        le_wat_avg_diff_norm_ef_pft,
        plots_dir,
        'strm_avg_diff_norm_ef_pft_vs_le_wat_avg_diff_norm_ef_pft',
        '(EF LE - PFT LE)/avg(le_obs) for each watershed',
        '(EF strm - PFT strm)/avg(strm_obs)',
        quadrant_lines=True,
        xlim=[-1.25,1.25],
        ylim=[-1.25,1.25],
        dot_size=2
    )
    #p_w = plot_wat()
    #p_w.plot_map(
    #    'le_j_wat_diff_ef_pft',
    #    watersheds,
    #    np.array(le_j_wat_diff_ef_pft)*-1,
    #    np.nanmean(le_j_wat_diff_ef_pft*-1),
    #    plots_dir,
    #    cmap='PiYG',
    #    vmin=-0.2,
    #    vmax=0.2
    #)
    #p_w.plot_map(
    #    'le_j_wat_yr_diff_ef_pft',
    #    watersheds,
    #    np.array(le_j_wat_yr_diff_ef_pft)*-1,
    #    np.nanmean(le_j_wat_yr_diff_ef_pft*-1),
    #    plots_dir,
    #    cmap='PiYG',
    #    vmin=-0.2,
    #    vmax=0.2
    #)
    # let's make the maps of the g1 values where each of our 5 pfts are over
    # 50% of than land cover in that pixel
    tile_pft_info = pd.read_csv(tile_pft_info_fname)
    tile_pft_info = tile_pft_info.set_index('tile')
    simple_map = {
        'Needleleaf evergreen temperate tree':'needleleaf_trees',
        'Needleleaf evergreen boreal tree':'needleleaf_trees',
        'Needleleaf deciduous boreal tree':'broadleaf_trees',
        'Broadleaf evergreen tropical tree':'broaleaf_trees',
        'Broadleaf evergreen temperate tree':'broadleaf_trees',
        'Broadleaf deciduous tropical tree':'broadleaf_trees',
        'Broadleaf deciduous temperate tree':'broadleaf_trees',
        'Broadleaf deciduous boreal tree':'broadleaf_trees',
        'Broadleaf evergreen temperate shrub':'shrub',
        'Broadleaf deciduous boreal shrub':'shrub',
        'Broadleaf deciduous temperate shrub':'shrub',
        'Broadleaf deciduous temperate shrub[moisture +deciduous]': (
            'shrub'
        ),
        'Broadleaf deciduous temperate shrub[moisture stress only]': (
            'shrub'
        ),
        'Arctic c3 grass':'c3_grass',
        'Cool c3 grass':'c3_grass',
        'Cool c3 grass [moisture + deciduous]':'c3_grass',
        'Cool c3 grass [moisture stress only]':'c3_grass',
        'Warm c4 grass [moisture + deciduous]':'c4_grass',
        'Warm c4 grass [moisture stress only]':'c4_grass',
        'Warm c4 grass':'c4_grass',
        'Crop':'crop',
        'Crop [moisture + deciduous]':'crop',
        'Crop [moisture stress only]':'crop',
        '(Spring temperate cereal)':'crop',
        '(Irrigated corn)':'crop',
        '(Soybean)':'crop',
        '(Corn)':'crop'
    }
    # let's see what the breakdown of these pixels is!
    needleleaf_perc = np.zeros(len(pixels))
    needleleaf_pix = np.zeros(0)
    needleleaf_g1_pft = np.zeros(0)
    needleleaf_g1_ef = np.zeros(0)
    broadleaf_perc = np.zeros(len(pixels))
    broadleaf_pix = np.zeros(0)
    broadleaf_g1_pft = np.zeros(0)
    broadleaf_g1_ef = np.zeros(0)
    c3_perc = np.zeros(len(pixels))
    c3_pix = np.zeros(0)
    c3_g1_pft = np.zeros(0)
    c3_g1_ef = np.zeros(0)
    c4_perc = np.zeros(len(pixels))
    c4_pix = np.zeros(0)
    c4_g1_pft = np.zeros(0)
    c4_g1_ef = np.zeros(0)
    crop_perc = np.zeros(len(pixels))
    crop_pix = np.zeros(0)
    crop_g1_pft = np.zeros(0)
    crop_g1_ef = np.zeros(0)
    g1_pft_np = np.array(g1_pft.loc['g1'])
    g1_ef_np = np.array(g1_ef.loc['g1'])
    largest_pft_perc_pix = np.zeros(len(pixels))
    threshold_perc = 50
    for p,pix in enumerate(pixels):
        #print(pix)
        this_need_perc = 0
        this_broad_perc = 0
        this_c3_perc = 0
        this_c4_perc = 0
        this_crop_perc = 0
        for pft in range(4):
            this_pft = tile_pft_info[
                'pft_{}_name'.format(
                    pft+1
                )
            ].loc[pix]
            this_pft = simple_map[this_pft]
            this_perc = tile_pft_info[
                'pft_{}_perc'.format(
                    pft+1
                )
            ].loc[pix]
            if this_pft == 'needleleaf_trees':
                this_need_perc += this_perc
            elif this_pft == 'broadleaf_trees':
                this_broad_perc += this_perc
            elif this_pft == 'c3_grass':
                this_c3_perc += this_perc
            elif this_pft == 'c4_grass':
                this_c4_perc += this_perc
            elif this_pft == 'crop':
                this_crop_perc += this_perc
            #print('pft {} name: {}'.format(pft+1,this_pft))
            #print('pft {} perc: {}'.format(pft+1,this_perc))
        needleleaf_perc[p] = this_need_perc
        broadleaf_perc[p] = this_broad_perc
        c3_perc[p] = this_c3_perc
        c4_perc[p] = this_c4_perc
        crop_perc[p] = this_crop_perc
        this_largest_perc = np.max(
            [
                this_need_perc,
                this_broad_perc,
                this_c3_perc,
                this_c4_perc,
                this_crop_perc
            ]
        )
        largest_pft_perc_pix[p] = this_largest_perc
        # lets get the % coverage of the most dominant pft in each pixel
        # okay, now lets assign pixels based off of whether they are majority
        # for any of the land cover types
        if needleleaf_perc[p] >= threshold_perc:
            needleleaf_pix = np.append(needleleaf_pix,pix)
            needleleaf_g1_pft = np.append(needleleaf_g1_pft,g1_pft_np[p])
            needleleaf_g1_ef = np.append(needleleaf_g1_ef,g1_ef_np[p])
        if broadleaf_perc[p] >= threshold_perc:
            broadleaf_pix = np.append(broadleaf_pix,pix)
            broadleaf_g1_pft = np.append(broadleaf_g1_pft,g1_pft_np[p])
            broadleaf_g1_ef = np.append(broadleaf_g1_ef,g1_ef_np[p])
        if c3_perc[p] >= threshold_perc:
            c3_pix = np.append(c3_pix,pix)
            c3_g1_pft = np.append(c3_g1_pft,g1_pft_np[p])
            c3_g1_ef = np.append(c3_g1_ef,g1_ef_np[p])
        if c4_perc[p] >= threshold_perc:
            c4_pix = np.append(c4_pix,pix)
            c4_g1_pft = np.append(c4_g1_pft,g1_pft_np[p])
            c4_g1_ef = np.append(c4_g1_ef,g1_ef_np[p])
        if crop_perc[p] >= threshold_perc:
            crop_pix = np.append(crop_pix,pix)
            crop_g1_pft = np.append(crop_g1_pft,g1_pft_np[p])
            crop_g1_ef = np.append(crop_g1_ef,g1_ef_np[p])
    # alright, let's plot the g1 values for each of the landcover types where
    # they are majority
    p_p.plot_map(
        'g1_pft_needleleaf_over_{}'.format(threshold_perc),
        needleleaf_pix,
        needleleaf_g1_pft,
        np.nanmean(needleleaf_g1_pft),
        plots_dir,
        vmin=0.5,
        vmax=10
    )
    p_p.plot_map(
        'g1_ef_needleleaf_over_{}'.format(threshold_perc),
        needleleaf_pix,
        needleleaf_g1_ef,
        np.nanmean(needleleaf_g1_ef),
        plots_dir,
        vmin=0.5,
        vmax=10
    )
    needleleaf_pix_matters_c = []
    for p,pix in enumerate(needleleaf_pix):
        if pix in big_improve_pix:
            needleleaf_pix_matters_c.append('#a65db7')
        else:
            needleleaf_pix_matters_c.append('#42bcf5')
    p_o.scatter(
        needleleaf_g1_pft,
        needleleaf_g1_ef,
        plots_dir,
        'needleleaf_g1_pft_over_{perc}_vs_needleleaf_g1_ef_over_{perc}'.format(
            perc=threshold_perc
        ),
        'g1 for PFT, pixels where needleleaf > {perc}%'.format(
            perc=threshold_perc
        ),
        'g1 for EF, pixels where needleleaf > {perc}%'.format(
            perc=threshold_perc
        ),
        color=needleleaf_pix_matters_c
    )
    p_p.plot_map(
        'g1_pft_broadleaf_over_{}'.format(threshold_perc),
        broadleaf_pix,
        broadleaf_g1_pft,
        np.nanmean(broadleaf_g1_pft),
        plots_dir,
        vmin=0.5,
        vmax=6.5
    )
    p_p.plot_map(
        'g1_ef_broadleaf_over_{}'.format(threshold_perc),
        broadleaf_pix,
        broadleaf_g1_ef,
        np.nanmean(broadleaf_g1_ef),
        plots_dir,
        vmin=0.5,
        vmax=6.5
    )
    broadleaf_pix_matters_c = []
    for p,pix in enumerate(broadleaf_pix):
        if pix in big_improve_pix:
            broadleaf_pix_matters_c.append('#a65db7')
        else:
            broadleaf_pix_matters_c.append('#42bcf5')
    p_o.scatter(
        broadleaf_g1_pft,
        broadleaf_g1_ef,
        plots_dir,
        'broadleaf_g1_pft_over_{perc}_vs_broadleaf_g1_ef_over_{perc}'.format(
            perc=threshold_perc
        ),
        'g1 for PFT, pixels where broadleaf > {perc}%'.format(
            perc=threshold_perc
        ),
        'g1 for EF, pixels where broadleaf > {perc}%'.format(
            perc=threshold_perc
        ),
        color=broadleaf_pix_matters_c
    )
    p_p.plot_map(
        'g1_pft_c3_over_{}'.format(threshold_perc),
        c3_pix,
        c3_g1_pft,
        np.nanmean(c3_g1_pft),
        plots_dir,
        vmin=0.5,
        vmax=3.75
    )
    p_p.plot_map(
        'g1_ef_c3_over_{}'.format(threshold_perc),
        c3_pix,
        c3_g1_ef,
        np.nanmean(c3_g1_ef),
        plots_dir,
        vmin=0.5,
        vmax=3.75
    )
    p_p.plot_map(
        'g1_pft_c4_over_{}'.format(threshold_perc),
        c4_pix,
        c4_g1_pft,
        np.nanmean(c4_g1_pft),
        plots_dir,
        vmin=1,
        vmax=3
    )
    p_p.plot_map(
        'g1_ef_c4_over_{}'.format(threshold_perc),
        c4_pix,
        c4_g1_ef,
        np.nanmean(c4_g1_ef),
        plots_dir,
        vmin=1,
        vmax=3
    )
    p_p.plot_map(
        'g1_pft_crop_over_{}'.format(threshold_perc),
        crop_pix,
        crop_g1_pft,
        np.nanmean(crop_g1_pft),
        plots_dir,
        vmin=0.5,
        vmax=4.5
    )
    p_p.plot_map(
        'g1_ef_crop_over_{}'.format(threshold_perc),
        crop_pix,
        crop_g1_ef,
        np.nanmean(crop_g1_ef),
        plots_dir,
        vmin=0.5,
        vmax=4.5
    )
    crop_pix_matters_c = []
    for p,pix in enumerate(pixels):
        if pix in big_improve_pix:
            crop_pix_matters_c.append('#a65db7')
        else:
            crop_pix_matters_c.append('#42bcf5')
    #p_o.scatter(
    #    crop_g1_pft,
    #    crop_g1_ef,
    #    plots_dir,
    #    'crop_g1_pft_over_{perc}_vs_crop_g1_ef_over_{perc}'.format(
    #        perc=threshold_perc
    #    ),
    #    'g1 for PFT, pixels where crop > {perc}%'.format(
    #        perc=threshold_perc
    #    ),
    #    'g1 for EF, pixels where crop > {perc}%'.format(
    #        perc=threshold_perc
    #    ),
    #    color=crop_pix_matters_c
    #)
    # one last check--how does combined objective value do as pixel becomes
    # more homogenous?
    comb_j_pix_ef = strm_j_pix_ef + le_j_ef
    comb_j_pix_pft = strm_j_pix_pft + le_j_pft
    comb_j_pix_diff_ef_pft = comb_j_pix_ef - comb_j_pix_pft
    p_o.scatter(
        largest_pft_perc_pix,
        comb_j_pix_diff_ef_pft,
        plots_dir,
        'largest_pft_perc_pix_vs_comb_j_pix_diff_ef_pft',
        'PFT % for largest PFT in that pixel',
        'Jtot for that pixel',
        quadrant_lines=True
    )
    bias_wat_pft = timeseries_info[exps[1]]['wat_strm_errors'].loc['bias']
    bias_wat_pft = np.array(bias_wat_pft)
    bias_wat_norm_pft = bias_wat_pft/np.array(strm_obs.mean())
    strm_obs_mean = np.array(strm_obs.mean())
    g1_diff_idx = np.where(np.array(g1_wat_diff_ef_pft) > 2)
    strm_obs_mean_g1_diff = strm_obs_mean[g1_diff_idx]
    bias_wat_pft_g1_diff = bias_wat_pft[g1_diff_idx]
    all_pix_colors_g1_diff = np.array(all_wat_colors)[g1_diff_idx]
    p_o.scatter(
        strm_obs_mean_g1_diff,
        bias_wat_pft_g1_diff,
        plots_dir,
        'strm_obs_mean_vs_bias_wat_pft',
        'avg(strm_obs) (mm/day)',
        'bias in streamflow pred. (mm/day)',
        color=all_pix_colors_g1_diff,
        quadrant_lines=True,
        x_log=True
    )
    p_o.scatter(
        g1_wat_diff_ef_pft,
        bias_wat_norm_pft,
        plots_dir,
        'g1_wat_diff_ef_pft_vs_bias_wat_norm_pft',
        'g1 EF - g1 PFT (sqrt(kPa))',
        'bias in streamflow pred (mm/day)',
        color=all_wat_colors,
        dot_size=2
    )
    p_o.scatter(
        g1_wat_diff_ef_pft,
        strm_j_pft,
        plots_dir,
        'g1_wat_diff_ef_pft_vs_strm_j_pft',
        'g1 EF - g1 PFT (sqrt(kPa))',
        'Jstrm for PFT',
        color=all_wat_colors
    )
    p_o.scatter(
        strm_j_pft,
        strm_j_ef,
        plots_dir,
        'strm_j_pft_vs_strm_j_ef',
        'Jstrm for PFT',
        'Jstrm for EF',
        color=all_wat_colors,
        one_to_one_line=True
    )
    canopy_height = np.array(
        nc.Dataset(
            canopy_fname
        )['vals']
    )
    canopy_height_df = pd.DataFrame(columns=pixels)
    canopy_height_df.loc['canopy_height'] = canopy_height
    precip = np.array(
        nc.Dataset(
            precip_fname
        )['vals']
    )
    precip_df = pd.DataFrame(columns=pixels)
    precip_df.loc['precip'] = precip
    canopy_height_wat_df = gen.pix_df_to_wat_df(
        canopy_height_df,intersection_info
    )
    precip_wat_df = gen.pix_df_to_wat_df(
        precip_df,intersection_info
    )
    canopy_height_wat = np.array(canopy_height_wat_df.loc['canopy_height'])
    precip_wat = np.array(precip_wat_df.loc['precip'])
    p_o.scatter(
        canopy_height_wat,
        precip_wat,
        plots_dir,
        'canopy_height_wat_vs_precip_wat',
        'Canopy height (normalized)',
        'Precipitation (normalized)',
        color=all_wat_colors,
        quadrant_lines=True
    )
    g1_wat_ef = gen.pix_df_to_wat_df(
        g1_ef,intersection_info
    )
    g1_wat_pft = gen.pix_df_to_wat_df(
        g1_pft,intersection_info
    )
    p_o.scatter(
        g1_wat_pft,
        g1_wat_ef,
        plots_dir,
        'g1_wat_pft_vs_g1_wat_ef',
        'g1 from PFT at wat scale',
        'g1 from EF at wat scale',
        color=all_wat_colors,
        xlim=[0,8.3],
        ylim=[0,8.3]
    )
    le_wat_yr_avg_diff_ef_pft = (
        le_preds_wat_yr_ef - le_preds_wat_yr_pft
    ).mean()
    print(le_preds_wat_yr_ef)
    print(le_preds_wat_yr_pft)
    print(le_wat_yr_avg_diff_ef_pft)
    et_wat_yr_avg_diff_ef_pft = le_wat_yr_avg_diff_ef_pft/28.94
    le_wat_yr_avg_diff_norm_ef_pft = (
        le_wat_yr_avg_diff_ef_pft/le_obs_wat_yr.mean()
    )
    p_o.scatter(
        strm_avg_diff_ef_pft,
        et_wat_yr_avg_diff_ef_pft,
        plots_dir,
        'strm_avg_diff_ef_pft_vs_le_wat_yr_avg_diff_ef_pft',
        'EF strm - PFT strm',
        'EF LE - PFT LE at yearly for each watershed',
        quadrant_lines=True,
        xlim=[-0.5,1.3],
        ylim=[-0.5,1.3],
        dot_size=2
    )
    p_o.scatter(
        strm_avg_diff_norm_ef_pft,
        le_wat_yr_avg_diff_norm_ef_pft,
        plots_dir,
        'strm_avg_diff_norm_ef_pft_vs_le_wat_yr_avg_diff_norm_ef_pft',
        '(EF strm - PFT strm)/avg(strm_obs)',
        '(EF LE - PFT LE)/avg(le_obs) at yearly for each watershed',
        quadrant_lines=True,
        xlim=[-1.25,1.25],
        ylim=[-1.25,1.25],
        dot_size=2
    )
    p_o.scatter(
        g1_wat_diff_ef_pft,
        strm_avg_diff_ef_pft,
        plots_dir,
        'g1_wat_diff_ef_pft_vs_strm_avg_diff_ef_pft',
        'g1 EF - g1 PFT (sqrt(kPa))',
        'EF strm - PFT strm',
        quadrant_lines=True,
        dot_size=2,
        color=all_wat_colors
    )
    p_o.scatter(
        g1_wat_diff_ef_pft,
        strm_avg_diff_norm_ef_pft,
        plots_dir,
        'g1_wat_diff_ef_pft_vs_strm_avg_diff_norm_ef_pft',
        'g1 EF - g1 PFT (sqrt(kPa))',
        '(EF strm - PFT strm)/avg(strm_obs)',
        quadrant_lines=True,
        dot_size=2,
        color=all_wat_colors
    )








if __name__ == '__main__':
    main()
