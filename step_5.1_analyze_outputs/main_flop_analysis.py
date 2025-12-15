import sys
sys.path.append('/shared/pso/step_5.1_analyze_outputs/funcs')
sys.path.append('/shared/pso/step_5.1_analyze_outputs/tasks')
import os
import datetime
import numpy as np
import pickle as pickle
import math
import pandas as pd
import matplotlib as mpl
import netCDF4 as nc
from get_timeseries import get_timeseries
from rainfall_runoff_ratio import rainfall_runoff
from plot_timeseries import timeseries
from get_averages_and_error import averages_and_error
from plot_watersheds import plot_wat
from plot_other import plot_other
from general_functions import gen_funcs
from plot_pixels import plot_pixels
from get_drought import drought

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
        'g1-ai-et-strm-spin19921994-test19952014-mae-num1':{
            'dir':(
                '/shared/pso_outputs/g1_ai_et_strm_camels_spin19921994_' +
                'test19952914_mae/num_1/7'
            ),
            'load_or_save':'load',
            'default_type':'pso_output',
            'read_met_forcing':False,
            'timeseries_dir':os.path.join(
                save_dir,
                'g1-ai-et-strm-spin19921994-test19952014-mae-num1'
            ),
            'optimization_type':'pft',
            'iteration_number':1
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
        'g1-a0-a1-et-strm-spin19921994-test19952014-mae-num1':{
            'dir':(
                '/shared2/pso_outputs/g1_a0_a1_et_strm_camels_' +
                'spin19921994_test19952914_mae/num_1/29'
            ),
            'load_or_save':'load',
            'default_type':'pso_output',
            'read_met_forcing':False,
            'timeseries_dir':os.path.join(
                save_dir,
                'g1-a0-a1-et-strm-spin19921994-test19952014-mae-num1'
            ),
            'optimization_type':'ef',
            'iteration_number':1
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
    g = get_timeseries()
    le_obs = g.get_le_obs(le_truth_fname)
    # let's get the streamflow truth timeseries
    strm_obs = g.get_streamflow_obs(streamflow_truth_fname)
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
    # cool we have the information we need
    # let's calcaute the 1/3 of years that are wet/dry at each pixel
    d = drought()
    gen = gen_funcs()
    precip = timeseries_info[exps[0]]['pixel_raw_timeseries']['rainfsnowf']
    precip_mm = precip*86400
    precip_mm_week = gen.df_to_weekly(precip_mm)
    precip_mm_week_short = precip_mm_week[precip_mm_week.index.year >= 1995]
    precip_mm_mon = gen.df_to_monthly(precip_mm)
    precip_mm_mon = precip_mm_mon[
        (precip_mm_mon.index.month != 10) &
        (precip_mm_mon.index.month != 11) &
        (precip_mm_mon.index.month != 12) &
        (precip_mm_mon.index.month != 1) &
        (precip_mm_mon.index.month != 2) &
        (precip_mm_mon.index.month != 2)
    ]
    precip_mm_mon_short = precip_mm_mon[precip_mm_mon.index.year >= 1995]
    precip_mm_mon_short = gen.df_to_yearly(precip_mm_mon_short)
    precip_mm_yr = gen.df_to_yearly(precip_mm)
    precip_mm_yr_short = precip_mm_yr[precip_mm_yr.index.year >= 1995]
    rainfall = timeseries_info[exps[0]]['pixel_raw_timeseries']['rainf']
    rainf_mm = rainfall*86400
    # let's get the wettest and driest third of years by precip for shading
    # purposes
    pixels = list(precip.columns)
    num_pixels = len(precip.columns)
    num_years = len(list(precip_mm_mon_short.index))
    dry_frac = math.floor(2)
    total_colors = dry_frac
    wet_frac = num_years - dry_frac
    precip_mm_mon_ranked = d.rank(precip_mm_mon)
    precip_mm_mon_short_ranked = d.rank(precip_mm_mon_short)
    years = list(precip_mm_mon_short_ranked.index.year)
    for p,pix in enumerate(pixels):
        precip_mm_mon_short_ranked_d_idx = np.where(
            precip_mm_mon_short_ranked[pix] <= dry_frac
        )
        precip_mm_mon_short_ranked_w_idx = np.where(
            precip_mm_mon_short_ranked[pix] > wet_frac
        )
        precip_mm_mon_short_ranked[pix].iloc[precip_mm_mon_short_ranked_d_idx] = 'd'
        precip_mm_mon_short_ranked[pix].iloc[precip_mm_mon_short_ranked_w_idx] = 'w'
    color_years_start_stop = [
        [
            [
                0 for i in range(2)
            ] for j in range(total_colors)
        ] for k in range(num_pixels)
    ]
    back_colors = [
        [
            0 for j in range(total_colors)
        ] for k in range(num_pixels)
    ]
    for p,pix in enumerate(pixels):
        this_rank = precip_mm_mon_short_ranked[pix]
        this_color_year = 0
        for i,idx in enumerate(years):
            if this_rank.iloc[i] == 'd':
                color_start = datetime.date(idx,4,1)
                color_end = datetime.date(idx,9,30)
                color_years_start_stop[p][this_color_year][:] = [
                    color_start,color_end
                ]
                back_colors[p][this_color_year] = 'r'
                this_color_year += 1
    # get g1 maps for EF and PFT, respectively
    pft_g1_map = timeseries_info[exps[2]]['g1_map']
    ef_g1_map = timeseries_info[exps[4]]['g1_map']
    # calculate the difference between the g1 maps
    g1_diff = ef_g1_map - pft_g1_map
    abs_g1_diff = np.abs(g1_diff)
    # get average precip at each pixel
    avg_precip = precip_mm.mean()
    # standard deviation of precipitation at each pixel
    std_precip = precip_mm_yr.std()
    std_precip_norm = std_precip/avg_precip
    # get average snowfall at each pixel
    avg_rainf = rainf_mm.mean()
    avg_snowf = avg_precip - avg_rainf
    avg_snowf_perc = avg_snowf/avg_precip
    # get average temperature at each pixel
    temp = timeseries_info[exps[0]]['pixel_raw_timeseries']['tair']
    avg_temp = temp.mean()
    # get average Februrary temp at each pixel
    feb_temp = temp[temp.index.month == 2]
    avg_feb_temp = feb_temp.mean()
    # calcualte max ET at each site
    pixel_max_et = le_obs.max()
    # get wetness at wilting point
    restart = nc.Dataset(restart_fname)
    wilting_point = np.array(restart['WPWET'])
    # Jet of each pixel
    avg_le_obs = le_obs.mean()
    le_mae_pft = timeseries_info[exps[2]]['pixel_le_errors'].loc['mae']
    j_le_pft = le_mae_pft/avg_le_obs
    le_mae_ef = timeseries_info[exps[4]]['pixel_le_errors'].loc['mae']
    j_le_ef = le_mae_ef/avg_le_obs
    j_le_diff = j_le_ef - j_le_pft
    j_le_diff_abs = j_le_diff.abs()
    # get the percent difference of each day between EF and PFT
    pft_le = timeseries_info[exps[2]]['pixel_raw_timeseries']['le']
    ef_le = timeseries_info[exps[4]]['pixel_raw_timeseries']['le']
    ef_pft_le_diff = ef_le - pft_le
    perc_diff = ef_pft_le_diff/pixel_max_et
    abs_perc_diff = perc_diff.abs()
    # array for percentage that one is greater
    perc_greater = np.zeros(len(pixels))
    perc_greater_pixels = np.zeros(0)
    num_days_req = 0.2
    method_diff_req = 0.1
    # loop over each site
    for p,pi in enumerate(pixels):
        this_diff = perc_diff[pi]
        this_abs_diff = abs_perc_diff[pi]
        # caculate percent of days where difference is > 10%
        num_days = len(this_abs_diff)
        greater_10 = this_abs_diff.loc[this_abs_diff > method_diff_req]
        num_greater_10 = len(greater_10)
        perc_greater_10 = num_greater_10/num_days
        # for if more than 10% of pixels have greater than 10% difference,
        # calculate how often one method is larger in those days with at least
        # 10% difference
        if perc_greater_10 > num_days_req:
            perc_greater_pixels = np.append(
                perc_greater_pixels,pi
            )
            big_vals = this_diff.loc[this_abs_diff > method_diff_req]
            num_big_vals = len(big_vals)
            pos_big_vals = big_vals[big_vals > 0]
            perc_pos_big_vals = len(pos_big_vals)/num_big_vals
            if perc_pos_big_vals < 0.5:
                perc_pos_big_vals = 1 - perc_pos_big_vals
            perc_greater[p] = perc_pos_big_vals
        else:
            perc_greater[p] = np.nan
    # plot g1 diff vs precip, colored by % of days
    perc_greater_norm = perc_greater/np.nanmax(perc_greater)
    cmap = mpl.cm.get_cmap('rainbow')
    colors = [
        [
            0 for r in range(4)
        ] for c in range(len(perc_greater))
    ]
    for p,per in enumerate(perc_greater_norm):
        if np.isnan(per):
            colors[p][:] = mpl.colors.to_rgba('silver')
        else:
            this_c = cmap(per)
            colors[p][:] = this_c
    # make all of the scatters
    p_o = plot_other()
    abs_g1_diff = abs_g1_diff.loc['g1']
    # abs g1 diff versus average precip
    p_o.scatter(
        abs_g1_diff,avg_precip,plots_dir,'abs_g1_diff_vs_avg_precip',
        'abs(g1 diff) EF vs. PFT','average daily precip (mm)',
        color=perc_greater
    )
    # abs g1 diff versus std of yearly precip
    p_o.scatter(
        abs_g1_diff,std_precip_norm,plots_dir,'abs_g1_diff_vs_std_precip',
        'abs(g1 diff) EF vs. PFT',
        'standard deviation of yearly precip normalized by average precip',
        color=perc_greater
    )
    # abs g1 diff versus average snowfall
    p_o.scatter(
        abs_g1_diff,avg_snowf,plots_dir,'abs_g1_diff_vs_avg_snowf',
        'abs(g1 diff) EF vs. PFT','average daily snowfall (mm)',
        color=perc_greater
    )
    # abs g1 diff versus average temp
    p_o.scatter(
        abs_g1_diff,avg_temp,plots_dir,'abs_g1_diff_vs_avg_temp',
        'abs(g1 diff) EF vs. PFT','average daily temperature (K)',
        color=perc_greater
    )
    # abs g1 diff versus average feb temp
    p_o.scatter(
        abs_g1_diff,avg_feb_temp,plots_dir,'abs_g1_diff_vs_avg_feb_temp',
        'abs(g1 diff) EF vs. PFT','average daily Feb. temp (K)',
        color=perc_greater
    )
    # abs g1 diff versus average % of precip as snow
    p_o.scatter(
        abs_g1_diff,avg_snowf_perc,plots_dir,'abs_g1_diff_vs_avg_snowf_perc',
        'abs(g1 diff) EF vs. PFT','average % of precip as snow (-)',
        color=perc_greater
    )
    # abs g1 diff versus wetness at wilting point
    p_o.scatter(
        abs_g1_diff,wilting_point,plots_dir,'abs_g1_diff_vs_wilting_point',
        'abs(g1 diff) EF vs. PFT','soil wetness at wilting point (-)',
        color=perc_greater
    )
    # abs g1 diff versus wetness at wilting point
    p_o.scatter(
        avg_precip,wilting_point,plots_dir,'avg_precip_vs_wilting_point',
        'average daily precip (mm)','soil wetness at wilting point (-)',
        color=perc_greater
    )
    # abs g1 diff versus std of yearly precip
    p_o.scatter(
        std_precip_norm,wilting_point,plots_dir,'std_precip_vs_wilting_point',
        'standard deviation of yearly precip normalized by average precip',
        'soil wetness at wilting point (-)',
        color=perc_greater
    )
    # abs g1 diff versus std of yearly precip
    p_o.scatter(
        abs_g1_diff,j_le_diff_abs,plots_dir,'abs_g1_diff_vs_abs_j_diff',
        'absolute value of g1 diff (sqrt(kPa))',
        'absolute value of J_LE diff (-)',
        color=perc_greater
    )
    # abs g1 diff versus std of yearly precip
    p_o.scatter(
        perc_greater,j_le_diff_abs,plots_dir,'perc_greater_vs_abs_j_le_diff',
        'perc_greater',
        'absolute value of J_LE diff (-)'
    )
    # map of the high difference pixels
    p_p = plot_pixels()
    nan_idx = np.where(np.isnan(perc_greater) == False)
    non_nan_perc_greater = perc_greater[nan_idx]
    p_p.plot_map(
        'percent_similar_difference',perc_greater_pixels,
        non_nan_perc_greater,np.nanmean(perc_greater),plots_dir
    )
    # calculate difference during wet and dry years
    wet_pft_j = np.zeros(len(pixels))
    wet_ef_j = np.zeros(len(pixels))
    dry_pft_j = np.zeros(len(pixels))
    dry_ef_j = np.zeros(len(pixels))
    non_dry_pft_j = np.zeros(len(pixels))
    non_dry_ef_j = np.zeros(len(pixels))
    for p,pix in enumerate(pixels):
        print(pix)
        this_pft_le = (
            timeseries_info[exps[2]]['pixel_raw_timeseries']['le'][pix]
        )
        this_pft_le = this_pft_le[
            (this_pft_le.index.month != 1) &
            (this_pft_le.index.month != 2) &
            (this_pft_le.index.month != 3) &
            (this_pft_le.index.month != 10) &
            (this_pft_le.index.month != 11) &
            (this_pft_le.index.month != 12)
        ]
        this_ef_le = (
            timeseries_info[exps[4]]['pixel_raw_timeseries']['le'][pix]
        )
        this_ef_le = this_ef_le[
            (this_ef_le.index.month != 1) &
            (this_ef_le.index.month != 2) &
            (this_ef_le.index.month != 3) &
            (this_ef_le.index.month != 10) &
            (this_ef_le.index.month != 11) &
            (this_ef_le.index.month != 12)
        ]
        this_obs_le = le_obs[pix]
        this_obs_le = this_obs_le[
            (this_obs_le.index.month != 1) &
            (this_obs_le.index.month != 2) &
            (this_obs_le.index.month != 3) &
            (this_obs_le.index.month != 10) &
            (this_obs_le.index.month != 11) &
            (this_obs_le.index.month != 12)
        ]
        wet_years = precip_mm_mon_short_ranked[
            precip_mm_mon_short_ranked[pix] == 'w'
        ].index.to_list()
        wet_pft_le = np.zeros(len(wet_years)*200)
        wet_ef_le = np.zeros(len(wet_years)*200)
        wet_obs_le = np.zeros(len(wet_years)*200)
        total_num = 0
        total_num_end = 0
        for y,yr in enumerate(wet_years):
            this_wet_pft_le = np.array(
                this_pft_le[
                    #(this_pft_le.index.month == yr.month) &
                    (this_pft_le.index.year == yr.year)
                ]
            )
            num_this_year = len(this_wet_pft_le)
            total_num_end += num_this_year
            wet_pft_le[total_num:total_num_end] = this_wet_pft_le
            this_wet_ef_le = np.array(
                this_ef_le[
                    #(this_ef_le.index.month == yr.month) &
                    (this_ef_le.index.year == yr.year)
                ]
            )
            wet_ef_le[total_num:total_num_end] = this_wet_ef_le
            this_wet_obs_le = np.array(
                this_obs_le[
                    #(this_obs_le.index.month == yr.month) &
                    (this_obs_le.index.year == yr.year)
                ]
            )
            wet_obs_le[total_num:total_num_end] = this_wet_obs_le
            total_num += num_this_year
        len_arrays = len(wet_years)*200-1
        wet_pft_le = np.delete(
            wet_pft_le,
            [total_num,len_arrays]
        )
        wet_ef_le = np.delete(
            wet_ef_le,
            [total_num,len_arrays]
        )
        wet_obs_le = np.delete(
            wet_obs_le,
            [total_num,len_arrays]
        )
        wet_pft_le = pd.Series(wet_pft_le)
        wet_ef_le = pd.Series(wet_ef_le)
        wet_obs_le = pd.Series(wet_obs_le)
        this_obs_avg = wet_obs_le.mean()
        this_wet_pft_mae = g_a.mae(wet_pft_le,wet_obs_le)
        this_wet_pft_j = this_wet_pft_mae/this_obs_avg
        wet_pft_j[p] = this_wet_pft_j
        this_wet_ef_mae = g_a.mae(wet_ef_le,wet_obs_le)
        this_wet_ef_j = this_wet_ef_mae/this_obs_avg
        wet_ef_j[p] = this_wet_ef_j
        #now for dry
        dry_years = precip_mm_mon_short_ranked[
            precip_mm_mon_short_ranked[pix] == 'd'
        ].index.to_list()
        dry_pft_le = np.zeros(len(dry_years)*200)
        dry_ef_le = np.zeros(len(dry_years)*200)
        dry_obs_le = np.zeros(len(dry_years)*200)
        total_num = 0
        total_num_end = 0
        for y,yr in enumerate(dry_years):
            this_dry_pft_le = np.array(
                this_pft_le[
                    #(this_pft_le.index.month == yr.month) &
                    (this_pft_le.index.year == yr.year)
                ]
            )
            num_this_year = len(this_dry_pft_le)
            total_num_end += num_this_year
            dry_pft_le[total_num:total_num_end] = this_dry_pft_le
            this_dry_ef_le = np.array(
                this_ef_le[
                    #(this_ef_le.index.month == yr.month) &
                    (this_ef_le.index.year == yr.year)
                ]
            )
            dry_ef_le[total_num:total_num_end] = this_dry_ef_le
            this_dry_obs_le = np.array(
                this_obs_le[
                    #(this_obs_le.index.month == yr.month) &
                    (this_obs_le.index.year == yr.year)
                ]
            )
            dry_obs_le[total_num:total_num_end] = this_dry_obs_le
            total_num += num_this_year
        len_arrays = len(dry_years)*200-1
        dry_pft_le = np.delete(
            dry_pft_le,
            [total_num,len_arrays]
        )
        dry_ef_le = np.delete(
            dry_ef_le,
            [total_num,len_arrays]
        )
        dry_obs_le = np.delete(
            dry_obs_le,
            [total_num,len_arrays]
        )
        dry_pft_le = pd.Series(dry_pft_le)
        dry_ef_le = pd.Series(dry_ef_le)
        dry_obs_le = pd.Series(dry_obs_le)
        this_obs_avg = dry_obs_le.mean()
        this_dry_pft_mae = g_a.mae(dry_pft_le,dry_obs_le)
        this_dry_pft_j = this_dry_pft_mae/this_obs_avg
        dry_pft_j[p] = this_dry_pft_j
        this_dry_ef_mae = g_a.mae(dry_ef_le,dry_obs_le)
        this_dry_ef_j = this_dry_ef_mae/this_obs_avg
        dry_ef_j[p] = this_dry_ef_j
        #now for non-dry
        non_dry_years = precip_mm_mon_short_ranked[
            precip_mm_mon_short_ranked[pix] != 'd'
        ].index.to_list()
        non_dry_pft_le = np.zeros(len(non_dry_years)*200)
        non_dry_ef_le = np.zeros(len(non_dry_years)*200)
        non_dry_obs_le = np.zeros(len(non_dry_years)*200)
        total_num = 0
        total_num_end = 0
        for y,yr in enumerate(non_dry_years):
            this_non_dry_pft_le = np.array(
                this_pft_le[
                    #(this_pft_le.index.month == yr.month) &
                    (this_pft_le.index.year == yr.year)
                ]
            )
            num_this_year = len(this_non_dry_pft_le)
            total_num_end += num_this_year
            non_dry_pft_le[total_num:total_num_end] = this_non_dry_pft_le
            this_non_dry_ef_le = np.array(
                this_ef_le[
                    #(this_ef_le.index.month == yr.month) &
                    (this_ef_le.index.year == yr.year)
                ]
            )
            non_dry_ef_le[total_num:total_num_end] = this_non_dry_ef_le
            this_non_dry_obs_le = np.array(
                this_obs_le[
                    #(this_obs_le.index.month == yr.month) &
                    (this_obs_le.index.year == yr.year)
                ]
            )
            non_dry_obs_le[total_num:total_num_end] = this_non_dry_obs_le
            total_num += num_this_year
        len_arrays = len(non_dry_years)*200-1
        non_dry_pft_le = np.delete(
            non_dry_pft_le,
            [total_num,len_arrays]
        )
        non_dry_ef_le = np.delete(
            non_dry_ef_le,
            [total_num,len_arrays]
        )
        non_dry_obs_le = np.delete(
            non_dry_obs_le,
            [total_num,len_arrays]
        )
        non_dry_pft_le = pd.Series(non_dry_pft_le)
        non_dry_ef_le = pd.Series(non_dry_ef_le)
        non_dry_obs_le = pd.Series(non_dry_obs_le)
        this_obs_avg = non_dry_obs_le.mean()
        this_non_dry_pft_mae = g_a.mae(non_dry_pft_le,non_dry_obs_le)
        this_non_dry_pft_j = this_non_dry_pft_mae/this_obs_avg
        non_dry_pft_j[p] = this_non_dry_pft_j
        this_non_dry_ef_mae = g_a.mae(non_dry_ef_le,non_dry_obs_le)
        this_non_dry_ef_j = this_non_dry_ef_mae/this_obs_avg
        non_dry_ef_j[p] = this_non_dry_ef_j
    pft_j_wet_dry_diff = wet_pft_j - dry_pft_j
    pft_j_wet_dry_abs_diff = np.abs(pft_j_wet_dry_diff)
    ef_j_wet_dry_diff = wet_ef_j - dry_ef_j
    ef_j_wet_dry_abs_diff = np.abs(ef_j_wet_dry_diff)
    pft_j_nondry_dry_diff = non_dry_pft_j - dry_pft_j
    ef_j_nondry_dry_diff = non_dry_ef_j - dry_ef_j
    diff_ef_j_nondry_dry_pft_j_nondry_dry = (
        ef_j_nondry_dry_diff - pft_j_nondry_dry_diff
    )
    diff_ef_pft_j_dry = dry_ef_j - dry_pft_j
    diff_ef_pft_j_nondry = non_dry_ef_j - non_dry_pft_j
    all_pft_j = timeseries_info[exps[2]]['pixel_le_errors'].loc['mae']
    all_pft_j = all_pft_j/le_obs.mean()
    all_ef_j = timeseries_info[exps[4]]['pixel_le_errors'].loc['mae']
    all_ef_j = all_ef_j/le_obs.mean()
    all_default_j = timeseries_info[exps[0]]['pixel_le_errors'].loc['mae']
    all_default_j = all_default_j/le_obs.mean()
    all_pft_j_strm = timeseries_info[exps[2]]['wat_strm_errors'].loc['mae']
    all_pft_j_strm = all_pft_j_strm/strm_obs.mean()
    all_ef_j_strm = timeseries_info[exps[4]]['wat_strm_errors'].loc['mae']
    all_ef_j_strm = all_ef_j_strm/strm_obs.mean()
    all_default_j_strm = timeseries_info[exps[0]]['wat_strm_errors'].loc['mae']
    all_default_j_strm = all_default_j_strm/strm_obs.mean()
    diff_ef_pft_j = all_ef_j - all_pft_j
    diff_ef_default_j = all_ef_j - all_default_j
    diff_ef_pft_j_strm = all_ef_j_strm - all_pft_j_strm
    diff_ef_default_j_strm = all_ef_j_strm - all_default_j_strm
    ef_j_all_dry_diff = all_ef_j - dry_ef_j
    pft_j_all_dry_diff = all_pft_j - dry_pft_j
    diff_ef_j_all_dry_pft_j_all_dry = (
        ef_j_all_dry_diff - pft_j_all_dry_diff
    )
    # scatter plots of wet vs. dry differences as a function of flip flop
    # metric
    num_each = len(all_pft_j)
    num_dots = len(all_pft_j)*2
    c_pft_ef = [0 for j in range(num_dots)]
    for c in range(num_each):
        c_pft_ef[c] = 'r'
    for c in range(num_each):
        c_pft_ef[c+num_each] = 'b'
    all_pft_plus_ef_j_et = np.concatenate(
        (all_pft_j,all_ef_j)
    )
    dry_pft_plus_ef_j_et = np.concatenate(
        (dry_pft_j,dry_ef_j)
    )
    non_dry_pft_plus_ef_j_et = np.concatenate(
        (non_dry_pft_j,non_dry_ef_j)
    )
    p_o.scatter(
        all_pft_plus_ef_j_et,
        dry_pft_plus_ef_j_et,
        plots_dir,
        'all_j_vs_dry_j',
        'j_le_all',
        'j_le_dry',
        color=c_pft_ef,
        one_to_one_line=True,
        dot_size=0.1
    )
    p_o.scatter(
        diff_ef_pft_j_dry,
        diff_ef_pft_j_nondry,
        plots_dir,
        'diff_ef_pft_j_dry_vs_diff_ef_pft_j_nondry',
        'diff_ef_pft_j_dry',
        'diff_ef_pft_j_nondry',
        quadrant_lines=True,
        one_to_one_line=True,
        dot_size=0.2
    )
    p_o.scatter(
        wet_pft_j,dry_pft_j,plots_dir,
        'pft_wet_j_vs_dry_j',
        'j_le_pft_wet',
        'j_le_pft_dry',
        color=perc_greater,
        one_to_one_line=True
    )
    p_o.scatter(
        wet_ef_j,dry_ef_j,plots_dir,
        'ef_wet_j_vs_dry_j',
        'j_le_ef_wet',
        'j_le_ef_dry',
        color=perc_greater,
        one_to_one_line=True
    )
    p_o.scatter(
        dry_pft_j,non_dry_pft_j,plots_dir,
        'pft_dry_j_vs_non_dry_j',
        'j_le_pft_dry',
        'j_le_pft_non_dry',
        color=perc_greater,
        one_to_one_line=True
    )
    p_o.scatter(
        pft_j_all_dry_diff,ef_j_all_dry_diff,plots_dir,
        'pft_j_all_dry_diff_vs_ef_j_all_dry_diff',
        'pft_j_all_dry_diff',
        'ef_j_all_dry_diff',
        one_to_one_line=True,
        quadrant_lines=True,
        dot_size = 0.5
    )
    p_o.scatter(
        dry_ef_j,non_dry_ef_j,plots_dir,
        'ef_dry_j_vs_non_dry_j',
        'j_le_ef_dry',
        'j_le_ef_non_dry',
        color=perc_greater,
        one_to_one_line=True
    )
    # let's plot wilting point vs difference between EF_J and PFT_J, colored by
    # the g1 difference
    idx = np.where(np.abs(diff_ef_j_nondry_dry_pft_j_nondry_dry) > 1)
    print(idx)
    print(np.array(pixels)[idx])
    print('default strm obj:')
    print(np.nanmean(all_default_j_strm))
    print('default et obj:')
    print(np.nanmean(all_default_j))
    print('pft strm obj:')
    print(np.nanmean(all_pft_j))
    print('pft et obj:')
    print(np.nanmean(all_pft_j_strm))
    p_o.scatter(
        wilting_point,diff_ef_pft_j_dry,plots_dir,
        'wilting_point_vs_diff_dry_ef_j_dry_pft_j',
        'soil wetness at wilting point (-)',
        'ef_et_j_dry - pft_et_j_dry',
        color=np.array(std_precip_norm),
        ylim=[-1.25,1.25]
    )
    p_o.scatter(
        abs_g1_diff,diff_ef_pft_j_dry,plots_dir,
        'abs_diff_g1_vs_diff_dry_ef_j_dry_pft_j',
        'abs(diff g1)',
        'ef_et_j_dry - pft_et_j_dry',
        ylim=[-1.25,1.25]
    )
    p_o.scatter(
        wilting_point,diff_ef_pft_j,
        plots_dir,
        'wilting_point_vs_diff_ef_j_all_pft_j_all',
        'soil wetness at wilting point (-)',
        'ef_j_nondry_minus_dry - pft_j_nondry_minus_dry',
        color=np.array(std_precip_norm),
        ylim=[-1.25,1.25]
    )
    p_p.plot_map(
        'diff_ef_default_j_et',pixels,diff_ef_default_j*-1,
        np.nanmean(diff_ef_default_j),plots_dir,
        cmap='PiYG',vmin=-0.25,vmax=0.25
    )
    gen.add_to_gdf_and_save(
        '/shared/pso/step_1x_choose_tiles_large/outputs/really_chosen_tiles.geojson',
        './outputs/diff_ef_default_j_et.gdf',diff_ef_default_j*-1,vmin=-0.25,vmax=0.25,
        cmap='PiYG'
    )
    p_p.plot_map(
        'diff_ef_pft_j_et',pixels,diff_ef_pft_j*-1,
        np.nanmean(diff_ef_pft_j),plots_dir,
        cmap='PiYG',vmin=-0.2,vmax=0.2
    )
    gen.add_to_gdf_and_save(
        '/shared/pso/step_1x_choose_tiles_large/outputs/really_chosen_tiles.geojson',
        './outputs/diff_ef_pft_j_et.gdf',diff_ef_pft_j*-1,vmin=-0.2,vmax=0.2,
        cmap='PiYG'
    )
    p_p.plot_map(
        'diff_ef_pft_j_et_washed',pixels,diff_ef_pft_j*-1,
        np.nanmean(diff_ef_pft_j),plots_dir,
        cmap='PiYG',vmin=-0.4,vmax=0.4
    )
    gen.add_to_gdf_and_save(
        '/shared/pso/step_1x_choose_tiles_large/outputs/really_chosen_tiles.geojson',
        './outputs/diff_ef_pft_j_et_washed.gdf',diff_ef_pft_j*-1,vmin=-0.4,vmax=0.4,
        cmap='PiYG'
    )
    p_p.plot_map(
        'diff_ef_pft_j_dry_et',pixels,diff_ef_pft_j_dry*-1,
        np.nanmean(diff_ef_pft_j_dry),plots_dir,
        cmap='PiYG',vmin=-0.4,vmax=0.4
    )
    gen.add_to_gdf_and_save(
        '/shared/pso/step_1x_choose_tiles_large/outputs/really_chosen_tiles.geojson',
        './outputs/diff_ef_pft_j_et_dry.gdf',diff_ef_pft_j_dry*-1,vmin=-0.4,vmax=0.4,
        cmap='PiYG'
    )
    p_p.plot_map(
        'diff_ef_pft_j_nondry_et',pixels,diff_ef_pft_j_nondry*-1,
        np.nanmean(diff_ef_pft_j_dry*-1),plots_dir,
        cmap='PiYG',vmin=-0.2,vmax=0.2
    )
    le_obs_yr = gen.df_to_yearly(le_obs)
    le_ef_yr = gen.df_to_yearly(
        timeseries_info[exps[4]]['pixel_raw_timeseries']['le']
    )
    le_ef_yr_mae = g_a.var_error(
         le_ef_yr,le_obs_yr,start_err,end_err
    ).loc['mae']
    le_ef_yr_j = le_ef_yr_mae/le_obs_yr.mean()
    le_pft_yr = gen.df_to_yearly(
        timeseries_info[exps[2]]['pixel_raw_timeseries']['le']
    )
    le_pft_yr_mae = g_a.var_error(
         le_pft_yr,le_obs_yr,start_err,end_err
    ).loc['mae']
    le_pft_yr_j = le_pft_yr_mae/le_obs_yr.mean()
    le_ef_pft_yr_j = le_ef_yr_j - le_pft_yr_j
    p_p.plot_map(
        'diff_le_mae_ef_pft_j',pixels,le_ef_pft_yr_j,
        le_ef_pft_yr_j.mean(),plots_dir,vmin=-0.5,vmax=0.5,
        cmap='bwr'
    )
    
    sys.exit()
    p_w = plot_wat()
    hucs = list(strm_obs.columns)
    #p_w.plot_map(
    #    'diff_ef_pft_j_strm_cbarfull',hucs,np.array(diff_ef_pft_j_strm),
    #    np.mean(diff_ef_pft_j_strm),plots_dir,
    #    cmap='PiYG',vmin=-0.9,vmax=0.9
    #)
    p_w.plot_map(
        'diff_ef_pft_j_strm_cbarclipped',hucs,np.array(diff_ef_pft_j_strm)*-1,
        np.mean(diff_ef_pft_j_strm),plots_dir,
        cmap='PiYG',vmin=-0.4,vmax=0.4
    )
    gen.add_to_gdf_and_save(
        '/shared/pso/step_1x_choose_tiles_large/outputs/chosen_camels.geojson',
        './outputs/diff_ef_pft_j_strm.gdf',np.array(diff_ef_pft_j_strm)*-1,vmin=-0.4,vmax=0.4,
        cmap='PiYG',subselection=hucs
    )
    p_w.plot_map(
        'diff_ef_default_j_strm',hucs,np.array(diff_ef_default_j_strm)*-1,
        np.mean(diff_ef_default_j_strm),plots_dir,
        cmap='PiYG',vmin=-0.65,vmax=0.65
    )
    gen.add_to_gdf_and_save(
        '/shared/pso/step_1x_choose_tiles_large/outputs/chosen_camels.geojson',
        './outputs/diff_ef_default_j_strm.gdf',np.array(diff_ef_default_j_strm)*-1,vmin=-0.4,vmax=0.4,
        cmap='PiYG',subselection=hucs
    )
    sys.exit()
    large_idx = np.where(np.abs(diff_ef_pft_dry) > 0.2)
    print(np.array(pixels)[large_idx])
    print(diff_ef_pft_dry[large_idx])
    to_plot_one = [
        timeseries_info[exps[2]]['pixel_raw_timeseries'],
        timeseries_info[exps[4]]['pixel_raw_timeseries'],
        {'le':le_obs}
    ]
    names_one = [
        'PFT itertaion 5',
        'EF iteration 6',
        'GLEAM'
    ]
    sizes = np.repeat(0.5,len(names_one))
    colors = ['r','g','k']
    to_plot_two = {'rainfsnowf':precip_mm_week}
    names_two = 'MERRA2 weekly precip'
    start_plot = datetime.date(1992,1,1)
    end_plot = datetime.date(2014,12,31)
    p_t = timeseries()
    p_t.plot_two_var(
        to_plot_one,names_one,'le','W/m2',
        to_plot_two,names_two,'rainfsnowf','mm/day',
        'timeseries_ef_and_pft_drought',
        plots_dir,start=start_plot,end=end_plot,
        small_preds=sizes,colors=colors,
        figsize=(60,4),locations=pixels,
        size_two=5,yearly_ticks=True,
        back_times=color_years_start_stop,
        back_colors=back_colors
    )
    sys.exit()
    # lets isolate our corner of interest
    corner_idx = np.where(
        (std_precip_norm > 0.175) &
        (wilting_point > 0.23)
    )[0]
    #corner_idx = np.where(
    #    (std_precip_norm > 0.175) &
    #    (wilting_point > 0.23) &
    #    (perc_greater > 0.95)
    #)[0]
    pixels = np.array(pixels)
    corner_pixels = pixels[corner_idx]
    corner_perc = perc_greater[corner_idx]
    std_precip_norm = np.array(std_precip_norm)
    corner_std_precip_norm = std_precip_norm[corner_idx]
    corner_wilt = wilting_point[corner_idx]
    avg_precip = np.array(avg_precip)
    corner_precip = avg_precip[corner_idx]
    avg_feb_temp = np.array(avg_feb_temp)
    corner_feb_temp = avg_feb_temp[corner_idx]
    pft_g1 = np.array(pft_g1_map.loc['g1'])
    corner_pft_g1 = pft_g1[corner_idx]
    ef_g1 = np.array(ef_g1_map.loc['g1'])
    corner_ef_g1 = ef_g1[corner_idx]
    not_nan_idx = np.where(np.isnan(corner_perc) == False)
    corner_pixels = corner_pixels[not_nan_idx]
    corner_perc = corner_perc[not_nan_idx]
    corner_std_precip_norm = corner_std_precip_norm[not_nan_idx]
    corner_wilt = corner_wilt[not_nan_idx]
    corner_precip = corner_precip[not_nan_idx]
    corner_feb_temp = corner_feb_temp[not_nan_idx]
    corner_pft_g1 = corner_pft_g1[not_nan_idx]
    corner_ef_g1 = corner_ef_g1[not_nan_idx]
    print('corner pixels:')
    print(corner_pixels)
    print('corner percentages:')
    print(corner_perc)
    print('corner norm std precip:')
    print(corner_std_precip_norm)
    print('corner wilting point:')
    print(corner_wilt)
    print('corner avg precip:')
    print(corner_precip)
    print('corner pft g1:')
    print(corner_pft_g1)
    print('corner ef g1:')
    print(corner_ef_g1)
    # abs g1 diff versus std of yearly precip
    p_o.scatter(
        corner_precip,corner_perc,plots_dir,'corner_precip_vs_corner_perc',
        'average yearly precip (mm)',
        '% of days where one param method is greater',
        best_fit_line=True
    )
    # abs g1 diff versus std of yearly precip
    p_o.scatter(
        corner_feb_temp,corner_perc,plots_dir,'corner_feb_temp_vs_corner_perc',
        'average Feb. temp (K)',
        '% of days where one param method is greater',
        best_fit_line=True
    )
    # let's test whether there are big J differences for ET during wet and dry
    # years where perc_greater < 0.8

    # Let's plot the full LE timeseries, shaded by wet/dry
    pixels = corner_pixels
    plot_timeseries = False
    # lets plot pft, ef, and obs on the same timeseries
    # for le
    to_plot_one = [
        timeseries_info[exps[2]]['pixel_raw_timeseries'],
        timeseries_info[exps[4]]['pixel_raw_timeseries'],
        {'le':le_obs}
    ]
    names_one = [
        'PFT itertaion 5',
        'EF iteration 6',
        'GLEAM'
    ]
    sizes = np.repeat(0.5,len(names_one))
    colors = ['r','g','k']
    to_plot_two = {'rainfsnowf':precip_mm_week}
    names_two = 'MERRA2 weekly precip'
    start_plot = datetime.date(1992,1,1)
    end_plot = datetime.date(2014,12,31)
    if plot_timeseries:
        p_t.plot_two_var(
            to_plot_one,names_one,'le','W/m2',
            to_plot_two,names_two,'rainfsnowf','mm/day',
            'timeseries_ef_and_pft_drought',
            plots_dir,start=start_plot,end=end_plot,
            small_preds=sizes,colors=colors,
            figsize=(60,4),locations=pixels,
            size_two=5,yearly_ticks=True
        )
    # what about for soil moisture
    to_plot_one = [
        timeseries_info[exps[2]]['pixel_raw_timeseries'],
        timeseries_info[exps[4]]['pixel_raw_timeseries']
    ]
    names_one = [
        'PFT itertaion 5',
        'EF iteration 6'
    ]
    sizes = np.repeat(0.5,len(names_one))
    colors = ['r','g']
    if plot_timeseries:
        p_t.plot_two_var(
            to_plot_one,names_one,'root_sm','-',
            to_plot_two,names_two,'rainfsnowf','mm/day',
            'timeseries_ef_and_pft_drought',
            plots_dir,start=start_plot,end=end_plot,
            small_preds=sizes,colors=colors,
            figsize=(60,4),locations=pixels,
            size_two=5,yearly_ticks=True
        )
    # and finally for beta
    to_plot_one = [
        timeseries_info[exps[2]]['pixel_raw_timeseries'],
        timeseries_info[exps[4]]['pixel_raw_timeseries']
    ]
    names_one = [
        'PFT itertaion 5',
        'EF iteration 6'
    ]
    sizes = np.repeat(0.5,len(names_one))
    colors = ['r','g']
    if plot_timeseries:
        p_t.plot_two_var(
            to_plot_one,names_one,'beta','-',
            to_plot_two,names_two,'rainfsnowf','mm/day',
            'timeseries_ef_and_pft_drought',
            plots_dir,start=start_plot,end=end_plot,
            small_preds=sizes,colors=colors,
            figsize=(60,4),locations=pixels,
            size_two=5,yearly_ticks=True
        )
    # and finally for runoff
    to_plot_one = [
        timeseries_info[exps[2]]['pixel_raw_timeseries'],
        timeseries_info[exps[4]]['pixel_raw_timeseries']
    ]
    names_one = [
        'PFT itertaion 5',
        'EF iteration 6'
    ]
    sizes = np.repeat(0.5,len(names_one))
    colors = ['r','g']
    if plot_timeseries:
        p_t.plot_two_var(
            to_plot_one,names_one,'runoff','-',
            to_plot_two,names_two,'rainfsnowf','mm/day',
            'timeseries_ef_and_pft_drought',
            plots_dir,start=start_plot,end=end_plot,
            small_preds=sizes,colors=colors,
            figsize=(60,4),locations=pixels,
            size_two=5,yearly_ticks=True
        )






if __name__ == '__main__':
    main()
