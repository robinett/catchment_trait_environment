import sys
sys.path.append('/discover/nobackup/trobinet/from_aws/pso/step_5.1_analyze_outputs/funcs')
sys.path.append('/discover/nobackup/trobinet/from_aws/pso/step_5.1_analyze_outputs/tasks')
import os
import datetime
import numpy as np
import pickle as pickle
import pandas as pd
import netCDF4 as nc
import statsmodels.formula.api as sm
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from get_timeseries import get_timeseries
from plot_timeseries import timeseries
from get_averages_and_error import averages_and_error
from plot_watersheds import plot_wat
from plot_other import plot_other
from general_functions import gen_funcs
from plot_pixels import plot_pixels
from plot_timeseries import timeseries
<<<<<<< HEAD
from scipy.stats import wilcoxon
import xarray as xr
import numpy as np

def paired_permutation_p(A, B, n_perm=100000, alternative="less", weights=None):
    A, B = np.asarray(A), np.asarray(B)
    m = np.isfinite(A) & np.isfinite(B)
    A, B = A[m], B[m]
    w = np.ones_like(A) if weights is None else np.asarray(weights)[m]
    w = w / w.sum()

    D = A - B
    obs = np.sum(w * D)  # weighted mean difference
    signs = (np.random.rand(n_perm, D.size) < 0.5).astype(int)*2 - 1
    perm = (signs * D).dot(w)

    if alternative == "less":   # A < B is better (errors)
        p = (np.sum(perm <= obs) + 1) / (n_perm + 1)
    elif alternative == "greater":
        p = (np.sum(perm >= obs) + 1) / (n_perm + 1)
    else:  # two-sided
        p = (np.sum(np.abs(perm) >= abs(obs)) + 1) / (n_perm + 1)
    return obs, p

=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)

def main():
    # what is the base dir?
    base_dir = '/discover/nobackup/trobinet/from_aws/pso/step_5.1_analyze_outputs'
    # where should we save plots?
    plots_dir = os.path.join(
        base_dir,'plots'
    )
    # where should we save computed outputs?
    out_dir = os.path.join(
        base_dir,'outputs'
    )
    save_dir = (
        '/discover/nobackup/trobinet/from_aws/pso/step_5_analyze_outputs/saved_timeseries'
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
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/intersecting_catch_tiles.csv'
    )
    # where is the intersection info located?
    intersection_info_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/intersection_info.pkl'
    )
    # where is the le truth located?
    le_truth_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_3x_process_gleam/outputs/' +
        'le_truth_gleam_38a_watts_per_m2_1995-01-01_2014-12-31_' +
        'camels_tiles.csv'
    )
    # where is the streamflow truth located?
    streamflow_truth_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_3.1.1x_process_camels/outputs/' +
        'camels_truth_yearly_1995-01-01_2014-12-31_mm_day.csv'
    )
<<<<<<< HEAD
    sm_truth_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_7_reviews/outputs/' +
        'sm_truth_cci_active_pct_1991-10-01_2014-12-31.csv'
    )
    sif_truth_fname= (
        '/discover/nobackup/trobinet/from_aws/pso/step_7_reviews/outputs/' +
        'sif_truth_oco2_all_daily_2001-01-01_2014-12-31.csv'
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # where is the tile info?
    tile_pft_info_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_5_analyze_outputs/outputs/pft_distribution.csv'
    )
    # where is the precip info for calc g1?
    precip_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_2x_env_covariates/outputs/' +
        'gldas_avg_precip_normalized.nc4'
    )
    # where is the canopy height info for g1 calc?
    canopy_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_2x_env_covariates/outputs/' +
        'canopy_height_normalized.nc4'
    )
<<<<<<< HEAD
    canopy_nonorm_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_2x_env_covariates/outputs/' +
        'canopy_height.nc4'
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # give me an example restart file (so we can get wetness at wilting point)
    restart_fname = (
        '/lustre/catchment/exps/GEOSldas_CN45_pso_g1_a0_a1_et_strm_' +
        'camels_spin19921994_test19952014_mae/0/output/SMAP_EASEv2_M36/' +
        'rs/ens0000/Y2015/M01/0.catchcnclm45_internal_rst.20150101_0000'
    )
    # we also need to know the positions that we used for pft and ef respectively
    positions_pft_fname = (
        '/discover/nobackup/trobinet/GEOSldas_pso_final_pft/' +
        'GEOSldas/positions.csv'
    )
    positions_ef_fname = (
        '/discover/nobackup/trobinet/GEOSldas_pso_final_ef/' +
        'GEOSldas/positions.csv'
    )
    # what are the names and info of the timeseries that we are going to load?
    timeseries_info = {
        'med-default-pft-g1-1992-2014':{
            'dir':(
<<<<<<< HEAD
                '/discover/nobackup/trobinet/from_aws/exps/GEOSldas_CN45_med_default_pft_g1_1992_2014_camels/0'
            ),
            'load_or_save':'load',
            'default_type':'default',
=======
                'nan'
            ),
            'load_or_save':'load',
            'default_type':'nan',
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
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
                '/discover/nobackup/trobinet/from_aws/pso_outputs/g1_ai_et_strm_camels_spin19921994_test19952914_mae/num_8/29'
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
                '/discover/nobackup/trobinet/from_aws/pso_outputs/g1_a0_a1_et_strm_camels_spin19921994_test19952914_mae/num_6/17'
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
    # lets also include the default Catchment-CN4.5 g1 values
    default_g1s = {
        'needleleaf_trees':2.3,
        'broadleaf_trees':4.25,
        'shrub':4.7,
        'c3_grass':5.3,
        'c4_grass':1.6,
        'crop':5.79
    }
    default_g1s_coefs = {
        'a0_needleaf_trees':[[2.3]],
        'a0_broadleaf_trees':[[4.25]],
        'a0_shrub':[[4.7]],
        'a0_c3_grass':[[5.3]],
        'a0_c4_grass':[[1.6]],
        'a0_crop':[[5.79]]
    }
<<<<<<< HEAD
    # just plot the canopy height and precip to get started
    p_p = plot_pixels()
    gen = gen_funcs()
    pixels = gen.get_pixels(pixels_fname)
    precip = np.array(
        nc.Dataset(
            precip_fname
        )['vals']
    )
    precip_df = pd.DataFrame(columns=pixels)
    precip_df.loc['precip'] = precip
    canopy = np.array(
        nc.Dataset(
            canopy_fname
        )['vals']
    )
    canopy_df = pd.DataFrame(columns=pixels)
    canopy_df.loc['canopy'] = canopy
    p_p.plot_map(
        'precip',
        pixels,
        precip,
        np.nanmean(precip),
        plots_dir
    )
    p_p.plot_map(
        'canopy',
        pixels,
        canopy,
        np.nanmean(canopy),
        plots_dir
    )
    canopy_nonorm = np.array(
        nc.Dataset(
            canopy_nonorm_fname
        )['vals']
    )
    canopy_nonorm_df = pd.DataFrame(columns=pixels)
    canopy_nonorm_df.loc['canopy'] = canopy
    p_p.plot_map(
        'canopy',
        pixels,
        canopy_nonorm,
        np.nanmean(canopy_nonorm),
        plots_dir
    )
    sys.exit()

    # let's get the le truth timeseries
=======
    # let's get the le truth timeseries
    gen = gen_funcs()
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    g = get_timeseries()
    le_obs = g.get_le_obs(le_truth_fname)
    # let's get the streamflow truth timeseries
    strm_obs = g.get_streamflow_obs(streamflow_truth_fname)
    strm_obs_mean = strm_obs.mean()
<<<<<<< HEAD
    # sm truth
    sm_obs = pd.read_csv(sm_truth_fname)
    sm_obs['time'] = pd.to_datetime(sm_obs['time'])
    sm_obs = sm_obs.set_index('time')
    # make all the column names strings
    sm_obs.columns = [int(col) for col in sm_obs.columns]
    # sif truth
    sif_obs = pd.read_csv(sif_truth_fname)
    sif_obs['time'] = pd.to_datetime(sif_obs['time'])
    sif_obs = sif_obs.set_index('time')
    # make all the column names strings
    sif_obs.columns = [int(col) for col in sif_obs.columns]
    # and why not pixels while we are at it
=======
    # and why not pixels while we are at it
    pixels = gen.get_pixels(pixels_fname)
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # and the intersection info
    intersection_info = gen.get_intersection_info(
        intersection_info_fname
    )
    # load the timeseries
    exps = list(timeseries_info.keys())
    timeseries_info = g.get_all_timeseries_info(
        timeseries_info,start,end,pixels_fname,
        intersection_info_fname,precip_fname,canopy_fname,
        tile_pft_info_fname
    )
<<<<<<< HEAD
    timeseries_info_drought = g.get_all_timeseries_info(
        timeseries_info,start,end,pixels_fname,
        intersection_info_fname,precip_fname,canopy_fname,
        tile_pft_info_fname
    )
    g_a = averages_and_error()
    print('getting overall errors')
    timeseries_info = g_a.get_all_averages_and_error(
        timeseries_info,le_obs,strm_obs,sm_obs,sif_obs,start_err,end_err
    )
    timeseries_info = copy.deepcopy(timeseries_info)
    g_a = averages_and_error()
    print('getting drought errors')
    timeseries_info_drought = g_a.get_all_averages_and_error(
        timeseries_info_drought,le_obs,strm_obs,sm_obs,sif_obs,start_err,end_err,
        during_drought=True
    )
    print(timeseries_info[exps[2]]['pixel_le_errors'])
    #print(timeseries_info_drought[exps[2]]['pixel_le_errors'])
=======
    g_a = averages_and_error()
    timeseries_info = g_a.get_all_averages_and_error(
        timeseries_info,le_obs,strm_obs,start_err,end_err
    )
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # get ef j for le and strm
    le_j_ef = (
        timeseries_info[exps[2]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
    strm_j_ef = (
        timeseries_info[exps[2]]['wat_strm_errors'].loc['mae']/
        strm_obs.mean()
    )
<<<<<<< HEAD
    le_j_ef_drought = (
        timeseries_info_drought[exps[2]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
    strm_j_ef_drought = (
        timeseries_info_drought[exps[2]]['wat_strm_errors'].loc['mae']/
        strm_obs.mean()
    )
    
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
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
<<<<<<< HEAD
    le_j_pft_drought = (
        timeseries_info_drought[exps[1]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
    overall_le_j_pft_drought = np.nanmean(le_j_pft_drought)
    strm_j_pft_drought = (
        timeseries_info_drought[exps[1]]['wat_strm_errors'].loc['mae']/
        strm_obs.mean()
    )
    overall_strm_j_pft_drought = np.nanmean(strm_j_pft_drought)
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # also get jLE at the watershed scale
    le_j_pft_df = pd.DataFrame(columns=pixels)
    le_j_pft_df.loc['le_j_pft'] = le_j_pft
    le_j_wat_pft_df = gen.pix_df_to_wat_df(
        le_j_pft_df,
        intersection_info
    )
    le_j_wat_pft = le_j_wat_pft_df.loc['le_j_pft']
<<<<<<< HEAD
    le_j_ef_df = pd.DataFrame(columns=pixels)
    le_j_ef_df.loc['le_j_ef'] = le_j_ef
    le_j_wat_ef_df = gen.pix_df_to_wat_df(
        le_j_ef_df,
        intersection_info
    )
    le_j_wat_ef = le_j_wat_ef_df.loc['le_j_ef']
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # get default j for le and strm
    le_j_def = (
        timeseries_info[exps[0]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
    strm_j_def = (
        timeseries_info[exps[0]]['wat_strm_errors'].loc['mae']/
        strm_obs.mean()
    )
<<<<<<< HEAD
    le_j_def_drought = (
        timeseries_info_drought[exps[0]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
    strm_j_def_drought = (
        timeseries_info_drought[exps[0]]['wat_strm_errors'].loc['mae']/
        strm_obs.mean()
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # get differences between j for EF and PFT
    watersheds = np.array(strm_j_pft.index)
    le_j_diff_ef_pft = le_j_ef - le_j_pft
    strm_j_diff_ef_pft = strm_j_ef - strm_j_pft
    le_j_diff_ef_def = le_j_ef - le_j_def
    strm_j_diff_ef_def = strm_j_ef - strm_j_def
<<<<<<< HEAD
    le_j_diff_pft_def = le_j_pft - le_j_def
    strm_j_diff_pft_def = strm_j_pft - strm_j_def
    le_j_diff_ef_pft_drought = le_j_ef_drought - le_j_pft_drought
    strm_j_diff_ef_pft_drought = strm_j_ef_drought - strm_j_pft_drought
    le_j_diff_ef_def_drought = le_j_ef_drought - le_j_def_drought
    strm_j_diff_ef_def_drought = strm_j_ef_drought - strm_j_def_drought
    le_j_diff_pft_def_drought = le_j_pft_drought - le_j_def_drought
    strm_j_diff_pft_def_drought = strm_j_pft_drought - strm_j_def_drought
    # also get signficance
    mask = np.isfinite(le_j_ef) & np.isfinite(le_j_pft)
    _, le_j_diff_ef_pft_p = paired_permutation_p(
        le_j_ef[mask],
        le_j_pft[mask],
        alternative='greater'
    )
    mask = np.isfinite(strm_j_ef) & np.isfinite(strm_j_pft)
    _, strm_j_diff_ef_pft_p = paired_permutation_p(
        strm_j_ef[mask],
        strm_j_pft[mask],
        alternative='less'
    )
    mask = np.isfinite(le_j_ef) & np.isfinite(le_j_def)
    _, le_j_diff_ef_def_p = paired_permutation_p(
        le_j_ef[mask],
        le_j_def[mask],
        alternative='less'
    )
    mask = np.isfinite(strm_j_ef) & np.isfinite(strm_j_def)
    _, strm_j_diff_ef_def_p = paired_permutation_p(
        strm_j_ef[mask],
        strm_j_def[mask],
        alternative='less'
    )
    mask = np.isfinite(le_j_ef_drought) & np.isfinite(le_j_pft_drought)
    _, le_j_diff_ef_pft_p_drought = paired_permutation_p(
        le_j_ef_drought[mask],
        le_j_pft_drought[mask],
        alternative='greater'
    )
    mask = np.isfinite(strm_j_ef_drought) & np.isfinite(strm_j_pft_drought)
    _, strm_j_diff_ef_pft_p_drought = paired_permutation_p(
        strm_j_ef_drought[mask],
        strm_j_pft_drought[mask],
        alternative='less'
    )
    mask = np.isfinite(le_j_ef_drought) & np.isfinite(le_j_def_drought)
    _, le_j_diff_ef_def_p_drought = paired_permutation_p(
        le_j_ef_drought[mask],
        le_j_def_drought[mask],
        alternative='greater'
    )
    mask = np.isfinite(strm_j_ef_drought) & np.isfinite(strm_j_def_drought)
    _, strm_j_diff_ef_def_p_drought = paired_permutation_p(
        strm_j_ef_drought[mask],
        strm_j_def_drought[mask],
        alternative='less'
    )
    mask = np.isfinite(le_j_pft_drought) & np.isfinite(le_j_def_drought)
    _, le_j_diff_pft_def_p_drought = paired_permutation_p(
        le_j_pft_drought[mask],
        le_j_def_drought[mask],
        alternative='greater'
    )
    mask = np.isfinite(strm_j_pft_drought) & np.isfinite(strm_j_def_drought)
    _, strm_j_diff_pft_def_p_drought = paired_permutation_p(
        strm_j_pft_drought[mask],
        strm_j_def_drought[mask],
        alternative='less'
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # let's convert le j diff to the watershed scale
    le_j_diff_ef_pft_df = pd.DataFrame(columns=pixels)
    le_j_diff_ef_pft_df.loc['le_j_diff'] = le_j_diff_ef_pft
    le_j_wat_diff_ef_pft_df = gen.pix_df_to_wat_df(
        le_j_diff_ef_pft_df,
        intersection_info
    )
    le_j_wat_diff_ef_pft = le_j_wat_diff_ef_pft_df.loc['le_j_diff']
    # also for compared to default
    le_j_diff_ef_def_df = pd.DataFrame(columns=pixels)
    le_j_diff_ef_def_df.loc['le_j_diff'] = le_j_diff_ef_def
    le_j_wat_diff_ef_def_df = gen.pix_df_to_wat_df(
        le_j_diff_ef_def_df,
        intersection_info
    )
    le_j_wat_diff_ef_def = le_j_wat_diff_ef_def_df.loc['le_j_diff']
    # print the % improvement over PFTs for LE and strm
    le_perc_imp = le_j_diff_ef_pft.mean()/le_j_pft.mean()*-100
    strm_perc_imp = strm_j_diff_ef_pft.mean()/strm_j_pft.mean()*-100
    le_perc_imp_def = le_j_diff_ef_def.mean()/le_j_def.mean()*-100
    strm_perc_imp_def = strm_j_diff_ef_def.mean()/strm_j_def.mean()*-100
<<<<<<< HEAD
    le_perc_imp_drought = le_j_diff_ef_pft_drought.mean()/le_j_pft_drought.mean()*-100
    strm_perc_imp_drought = strm_j_diff_ef_pft_drought.mean()/strm_j_pft_drought.mean()*-100
    le_perc_imp_def_drought = le_j_diff_ef_def_drought.mean()/le_j_def_drought.mean()*-100
    strm_perc_imp_def_drought = strm_j_diff_ef_def_drought.mean()/strm_j_def_drought.mean()*-100
    le_perc_imp_pft_def = le_j_diff_pft_def.mean()/le_j_def.mean()*-100
    strm_perc_imp_pft_def = strm_j_diff_pft_def.mean()/strm_j_def.mean()*-100
    le_perc_imp_pft_def_drought = le_j_diff_pft_def_drought.mean()/le_j_def_drought.mean()*-100
    strm_perc_imp_pft_def_drought = strm_j_diff_pft_def_drought.mean()/strm_j_def_drought.mean()*-100
    sm_corr_ef_pft = (
        timeseries_info[exps[2]]['pixel_sm_errors'].loc['corr'] -
        timeseries_info[exps[1]]['pixel_sm_errors'].loc['corr']
    )
    sm_corr_ef_def = (
        timeseries_info[exps[2]]['pixel_sm_errors'].loc['corr'] -
        timeseries_info[exps[0]]['pixel_sm_errors'].loc['corr']
    )
    sif_corr_ef_pft = (
        timeseries_info[exps[2]]['pixel_sif_errors'].loc['corr'] -
        timeseries_info[exps[1]]['pixel_sif_errors'].loc['corr']
    )
    sif_corr_ef_def = (
        timeseries_info[exps[2]]['pixel_sif_errors'].loc['corr'] -
        timeseries_info[exps[0]]['pixel_sif_errors'].loc['corr']
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    print(
        'EF improved LE by {}%'.format(le_perc_imp)
    )
    print(
<<<<<<< HEAD
        'p value of {}'.format(le_j_diff_ef_pft_p)
    )
    print(
        'EF improved strm by {}%'.format(strm_perc_imp)
    )
    print(
        'p value of {}'.format(strm_j_diff_ef_pft_p)
    )
    print(
        'EF improved LE for default by {}%'.format(le_perc_imp_def)
    )
    print(
        'p value of {}'.format(le_j_diff_ef_def_p)
    )
    print(
        'EF improved strm for default by {}%'.format(strm_perc_imp_def)
    )
    print(
        'p value of {}'.format(strm_j_diff_ef_def_p)
    )
    print(
        'PFT improved LE over default by {}%'.format(le_perc_imp_pft_def)
    )
    print(
        'PFT improved strm over default by {}%'.format(strm_perc_imp_pft_def)
    )
    print(
        'EF improved LE during drought by {}%'.format(le_perc_imp_drought)
    )
    print(
        'p value of {}'.format(le_j_diff_ef_pft_p_drought)
    )
    print(
        'EF improved strm during drought by {}%'.format(strm_perc_imp_drought)
    )
    print(
        'p value of {}'.format(strm_j_diff_ef_pft_p_drought)
    )
    print(
        'EF improved LE for default during drought by {}%'.format(le_perc_imp_def_drought)
    )
    print(
        'p value of {}'.format(le_j_diff_ef_def_p_drought)
    )
    print(
        'EF improved strm for default during drought by {}%'.format(strm_perc_imp_def_drought)
    )
    print(
        'p value of {}'.format(strm_j_diff_ef_def_p_drought)
    )
    print(
        'PFT improved LE over default by {} for drought%'.format(le_perc_imp_pft_def_drought)
    )
    print(
        'p value of {}'.format(le_j_diff_pft_def_p_drought)
    )
    print(
        'PFT improved strm over default by {}% for drought'.format(strm_perc_imp_pft_def_drought)
    )
    print(
        'p value of {}'.format(strm_j_diff_pft_def_p_drought)
    )
    print(
        'EF corr for SM: {}'.format(
            timeseries_info[exps[2]]['pixel_sm_errors'].loc['corr'].mean()
        )
    )
    print(
        'PFT corr for SM: {}'.format(
            timeseries_info[exps[1]]['pixel_sm_errors'].loc['corr'].mean()
        )
    )
    print(
        'Default corr for SM: {}'.format(
            timeseries_info[exps[0]]['pixel_sm_errors'].loc['corr'].mean()
        )
    )
    print(
        'EF corr for SIF: {}'.format(
            timeseries_info[exps[2]]['pixel_sif_errors'].loc['corr'].mean()
        )
    )
    print(
        'PFT corr for SIF: {}'.format(
            timeseries_info[exps[1]]['pixel_sif_errors'].loc['corr'].mean()
        )
    )
    print(
        'Default corr for SIF: {}'.format(
            timeseries_info[exps[0]]['pixel_sif_errors'].loc['corr'].mean()
        )
    )
    # get the pixel/basin-wise percentages
    le_perc = le_j_diff_ef_pft/le_j_pft*-100
    strm_perc = strm_j_diff_ef_pft/strm_j_pft*-100
    le_perc_def = le_j_diff_ef_def/le_j_def*-100
    strm_perc_def = strm_j_diff_ef_def/strm_j_def*-100
    le_perc_drought = le_j_diff_ef_pft_drought/le_j_pft_drought*-100
    strm_perc_drought = strm_j_diff_ef_pft_drought/strm_j_pft_drought*-100
    le_perc_def_drought = le_j_diff_ef_def_drought/le_j_def_drought*-100
    strm_perc_def_drought = strm_j_diff_ef_def_drought/strm_j_def_drought*-100
    le_perc_pft_def = le_j_diff_pft_def/le_j_def*-100
    strm_perc_pft_def = strm_j_diff_pft_def/strm_j_def*-100
    le_perc_pft_def_drought = le_j_diff_pft_def_drought/le_j_def_drought*-100
    strm_perc_pft_def_drought = strm_j_diff_pft_def_drought/strm_j_def_drought*-100
    # get just the overalls
    strm_perc_imp = strm_j_diff_ef_pft.mean()/strm_j_pft.mean()*-100
    le_perc_imp_def = le_j_diff_ef_def.mean()/le_j_def.mean()*-100
    strm_perc_imp_def = strm_j_diff_ef_def.mean()/strm_j_def.mean()*-100
    le_perc_imp_drought = le_j_diff_ef_pft_drought.mean()/le_j_pft_drought.mean()*-100
    strm_perc_imp_drought = strm_j_diff_ef_pft_drought.mean()/strm_j_pft_drought.mean()*-100
    le_perc_imp_def_drought = le_j_diff_ef_def_drought.mean()/le_j_def_drought.mean()*-100
    strm_perc_imp_def_drought = strm_j_diff_ef_def_drought.mean()/strm_j_def_drought.mean()*-100
    le_perc_imp_pft_def = le_j_diff_pft_def.mean()/le_j_def.mean()*-100
    strm_perc_imp_pft_def = strm_j_diff_pft_def.mean()/strm_j_def.mean()*-100
    le_perc_imp_pft_def_drought = le_j_diff_pft_def_drought.mean()/le_j_def_drought.mean()*-100
    strm_perc_imp_pft_def_drought = strm_j_diff_pft_def_drought.mean()/strm_j_def_drought.mean()*-100
=======
        'EF improved strm by {}%'.format(strm_perc_imp)
    )
    print(
        'EF improved LE for default by {}%'.format(le_perc_imp_def)
    )
    print(
        'EF improved strm for default by {}%'.format(strm_perc_imp_def)
    )
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    tot_j_def = le_j_def.mean() + strm_j_def.mean()
    tot_j_pft = le_j_pft.mean() + strm_j_pft.mean()
    tot_j_ef = le_j_ef.mean() + strm_j_ef.mean()
    tot_perc_imp = (tot_j_ef - tot_j_pft)/tot_j_pft*-100
    tot_perc_imp_def = (tot_j_ef - tot_j_def)/tot_j_def*-100
<<<<<<< HEAD
    tot_j_def_drought = le_j_def_drought.mean() + strm_j_def_drought.mean()
    tot_j_pft_drought = le_j_pft_drought.mean() + strm_j_pft_drought.mean()
    tot_j_ef_drought = le_j_ef_drought.mean() + strm_j_ef_drought.mean()
    tot_perc_imp_drought = (tot_j_ef_drought - tot_j_pft_drought)/tot_j_pft_drought*-100
    tot_perc_imp_def_drought = (tot_j_ef_drought - tot_j_def_drought)/tot_j_def_drought*-100
    tot_perc_imp_pft_def_drought = (tot_j_pft_drought - tot_j_def_drought)/tot_j_def_drought*-100
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    print(
        'EF improved TOT by {}%'.format(tot_perc_imp)
    )
    print(
        'EF improved TOT for def by {}%'.format(tot_perc_imp_def)
    )
<<<<<<< HEAD
    print(
        'EF improved TOT during drought by {}%'.format(tot_perc_imp_drought)
    )
    print(
        'EF improved TOT for def during drought by {}%'.format(tot_perc_imp_def_drought)
    )
    print(
        'PFT improved TOT for def during drought by {}%'.format(tot_perc_imp_pft_def_drought)
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # let's make figure 2a: PFT-based g1 vs. EF-based g1
    # first let's get the g1 information
    g1_pft_df = timeseries_info[exps[1]]['g1_map']
    g1_pft = g1_pft_df.loc['g1']
    g1_ef_df = timeseries_info[exps[2]]['g1_map']
    g1_ef = g1_ef_df.loc['g1']
    fake_obj = {}
    fake_obj['all'] = [[1]]
    g1_default_df = g.get_g1_map(
        default_g1s_coefs,
        fake_obj,
        precip_fname,
        canopy_fname,
        tile_pft_info_fname,
        'pft',
        0
    )
    g1_default = g1_default_df.loc['g1']
    g1_diff_ef_pft_df = g1_ef_df - g1_pft_df
    g1_diff_ef_pft = g1_diff_ef_pft_df.loc['g1']
<<<<<<< HEAD
    g1_diff_ef_default_df = g1_ef_df - g1_default_df
    g1_diff_ef_default = g1_diff_ef_default_df.loc['g1']
    g1_diff_pft_default_df = g1_pft_df - g1_default_df
    g1_diff_pft_default = g1_diff_pft_default_df.loc['g1']
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
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
<<<<<<< HEAD
    # this is a neutral point, so let's make it 
=======
    # this is a neutral point, so let's make it purple
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    all_purple = ['#8441B3' for i in range(len(pixels))]
    all_purple_wat = ['#8441B3' for i in range(len(watersheds))]
    # okay let's make the plot
    p_o = plot_other()
    p_o.scatter(
        g1_pft,
        g1_ef,
        plots_dir,
        'g1_pft_vs_g1_ef',
        'g1 PFT (sqrt(kPa))',
        'g1 EF (sqrt(kPa))',
        color=all_purple,
        xlim=[0,12],
        ylim=[0,12]
    )
    # now we need to save the gdf's so that they can be put into qgis to make
    # some pretty maps
<<<<<<< HEAD
=======
    p_p = plot_pixels()
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    p_p.plot_map(
        'g1_default',
        pixels,
        g1_default,
        np.nanmean(g1_default),
        plots_dir,
        vmin=0.5,
        vmax=7
    )
    p_p.plot_map(
        'g1_pft',
        pixels,
        g1_pft,
        np.nanmean(g1_pft),
        plots_dir,
        vmin=0.5,
        vmax=7
    )
    p_p.plot_map(
        'g1_ef',
        pixels,
        g1_ef,
        np.nanmean(g1_ef),
        plots_dir,
        vmin=0.5,
        vmax=7
    )
<<<<<<< HEAD
    p_p.plot_map(
        'g1_diff_ef_default',
        pixels,
        g1_diff_ef_default,
        np.nanmean(g1_diff_ef_default),
        plots_dir,
        cmap='PiYG',
        vmin=-5,
        vmax=5
    )
    p_p.plot_map(
        'g1_diff_pft_default',
        pixels,
        g1_diff_pft_default,
        np.nanmean(g1_diff_pft_default),
        plots_dir,
        cmap='PiYG',
        vmin=-5,
        vmax=5
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'g1_default.gpkg'
        ),
        g1_default,
        vmin=0.5,
        vmax=7,
        cmap='rainbow'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'g1_ef.gpkg'
        ),
        g1_ef,
        vmin=0.5,
        vmax=7,
        cmap='rainbow'
    )   
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'g1_pft.gpkg'
        ),
        g1_pft,
        vmin=0.5,
        vmax=7,
        cmap='rainbow'
    )
<<<<<<< HEAD
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'g1_diff_pft_default.gpkg'
        ),
        g1_diff_pft_default,
        vmin=-5,
        vmax=5,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'g1_diff_ef_default.gpkg'
        ),
        g1_diff_ef_default,
        vmin=-5,
        vmax=5,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'sm_corr_ef_pft.gpkg'
        ),
        sm_corr_ef_pft,
        vmin=-0.15,
        vmax=0.15,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'sm_corr_ef_def.gpkg'
        ),
        sm_corr_ef_def,
        vmin=-0.15,
        vmax=0.15,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'sif_corr_ef_pft.gpkg'
        ),
        sif_corr_ef_pft,
        vmin=-0.15,
        vmax=0.15,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'sif_corr_ef_def.gpkg'
        ),
        sif_corr_ef_def,
        vmin=-0.15,
        vmax=0.15,
        cmap='PiYG'
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # make the pdf for the different g1s
    all_g1_vals = [
        g1_default,g1_pft,g1_ef
    ]
    all_g1_labels = [
        'default','pft','ef'
    ]
    all_g1_colors = [
        '#000000',
        '#dc70ab',
        '#86c148'
    ]
    kernel_widths = [
        1,1,1
    ]
<<<<<<< HEAD
    line_styles = [
        '-','-','-'
    ]
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    line_widths = [1,1,1]
    fills = [False,False,False]
    p_o.plot_pdf(
        all_g1_vals,
        all_g1_labels,
        all_g1_colors,
        kernel_widths,
        line_widths,
        fills,
<<<<<<< HEAD
        line_styles,
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
        os.path.join(
            plots_dir,
            'g1_pdfs.svg'
        ),
        xlim=[0,8],
        min_val=[0.5]
    )
<<<<<<< HEAD
    # get changes in average le
    le_default = timeseries_info[exps[0]]['pixel_avgs'].loc['le']
    et_default = le_default/28.94
    le_pft = timeseries_info[exps[1]]['pixel_avgs'].loc['le']
    et_pft = le_pft/28.94
    le_ef = timeseries_info[exps[2]]['pixel_avgs'].loc['le']
    et_ef = le_ef/28.94
    et_diff_pft_default = et_pft - et_default
    et_diff_ef_default = et_ef - et_default
    p_p.plot_map(
        'et_diff_pft_default',
        pixels,
        et_diff_pft_default,
        np.nanmean(et_diff_pft_default),
        plots_dir,
        cmap='PiYG',
        vmin=-0.5,
        vmax=0.5
    )
    p_p.plot_map(
        'et_diff_ef_default',
        pixels,
        et_diff_ef_default,
        np.nanmean(et_diff_ef_default),
        plots_dir,
        cmap='PiYG',
        vmin=-0.5,
        vmax=0.5
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'et_diff_pft_default.gpkg'
        ),
        et_diff_pft_default,
        vmin=-0.6,
        vmax=0.6,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'et_diff_ef_default.gpkg'
        ),
        et_diff_ef_default,
        vmin=-0.6,
        vmax=0.6,
        cmap='PiYG'
    )
    # get changes in average strm
    strm_default = timeseries_info[exps[0]]['wat_avgs'].loc['strm']
    strm_pft = timeseries_info[exps[1]]['wat_avgs'].loc['strm']
    strm_ef = timeseries_info[exps[2]]['wat_avgs'].loc['strm']
    strm_diff_pft_default = strm_pft - strm_default
    strm_diff_ef_default = strm_ef - strm_default
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'chosen_camels.geojson',
        os.path.join(
            out_dir,
            'strm_diff_pft_default.gpkg'
        ),
        strm_diff_pft_default,
        vmin=-0.5,
        vmax=0.5,
        cmap='PiYG',
        subselection=watersheds
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'chosen_camels.geojson',
        os.path.join(
            out_dir,
            'strm_diff_ef_default.gpkg'
        ),
        strm_diff_ef_default,
        vmin=-0.5,
        vmax=0.5,
        cmap='PiYG',
        subselection=watersheds
    )
    #p_w = plot_wat()
    #p_w.plot_map(
    #    'strm_diff_pft_default',
    #    watersheds,
    #    strm_diff_pft_default,
    #    np.nanmean(strm_diff_pft_default),
    #    plots_dir,
    #    cmap='PiYG',
    #)
    #p_w.plot_map(
    #    'strm_diff_ef_default',
    #    watersheds,
    #    strm_diff_ef_default,
    #    np.nanmean(strm_diff_ef_default),
    #    plots_dir,
    #    cmap='PiYG'
    #)
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # alright perfect, now we need to make the figure 3 histograms
    # we need to get the information by PFT
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
    # we're also going to need the environmental covariates
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
    # get the basins where each PFT is > x%
    pft_thresh_perc = 50
    pfts = list(tile_pft_info_wat_df.index)
    pfts_plus = copy.deepcopy(pfts)
    #pfts_plus.append('no_majority')
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
    # okay great, now we can make the bar plots
    # first for the standard deviation of g1 in each of the pfts
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
    # and for the standard deviation of normalized environmental predictors for
    # each of the pfts
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
    # let's also make the g1 distribution on the training pixels
    # for each of the PFTs
    # okay let's make the maps of Jstrm and JET difference
    # we also want to do this for each pixel that contains any % of that PFT
    print(tile_pft_info_df)
    pfts = list(tile_pft_info_df.index)
    print(pfts)
    ef_columns = copy.deepcopy(pfts)
    ef_columns += ['b0']
    ef_columns += ['b1']
    ef_columns += ['b2']
    print(ef_columns)
    positions_pft = pd.read_csv(positions_pft_fname,names=pfts,header=None)
    positions_ef = pd.read_csv(positions_ef_fname,names=ef_columns,header=None)
    print(positions_pft)
    print(positions_ef)
    g1_by_pft = {
        'needleleaf_trees':pd.DataFrame(
            {'type':['pft','ef','perc']}
        ).set_index('type'),
        'broadleaf_trees':pd.DataFrame(
            {'type':['pft','ef','perc']}
        ).set_index('type'),
        'shrub':pd.DataFrame(
            {'type':['pft','ef','perc']}
        ).set_index('type'),
        'c3_grass':pd.DataFrame(
            {'type':['pft','ef','perc']}
        ).set_index('type'),
        'c4_grass':pd.DataFrame(
            {'type':['pft','ef','perc']}
        ).set_index('type'),
        'crop':pd.DataFrame(
            {'type':['pft','ef','perc']}
        ).set_index('type'),
    }
    # for each pixel
    for p,pix in enumerate(pixels):
        # for each pft
        for p,pft in enumerate(tile_pft_info_df.index.to_list()):
            # if not nan, calculate ef-based g1 and pft-based g1 for
            # that pft at that pixel
            this_val = tile_pft_info_df[pix].loc[pft]
            if pd.isnull(this_val) == False:
                # add to the pft-based series or g1-based series
                this_pft_g1 = positions_pft[pft]
                this_ef_g1 = (
                    positions_ef[pft]*(
                        positions_ef['b0'] +
                        positions_ef['b1']*precip_df[pix].loc['precip'] +
                        positions_ef['b2']*canopy_df[pix].loc['canopy']
                    )
                )
                this_pft_g1 = list(this_pft_g1)[0]
                this_ef_g1  = list(this_ef_g1)[0]
                if this_pft_g1 < 0.5:
                    this_pft_g1 = 0.5
                if this_ef_g1 < 0.5:
                    this_ef_g1 = 0.5
                comb_g1 = [this_pft_g1,this_ef_g1,this_val]
                g1_by_pft[pft][pix] = comb_g1
    # let's make a histogram of the distribution for each of the PFTs
    # first we need the Lin and default values
    lin_g1s = {
        'needleleaf_trees':2.35,
        'broadleaf_trees':3.97,
        'shrub':4.22,
        'c4_grass':1.62,
        'c3_grass':4.5,
        'crop':5.79
    }
    print(g1_by_pft)
    # get consolidated g1 values
    ef_medians = np.zeros(0)
    pft_vals = np.zeros(0)
    default_vals = np.zeros(0)
    lin_vals = np.zeros(0)
    for p,pft in enumerate(pfts):
        default_vals = np.append(default_vals,default_g1s[pft])
        lin_vals = np.append(lin_vals,lin_g1s[pft])
        this_pft_g1 = g1_by_pft[pft].loc['pft'].dropna()
        this_pft_g1_mean = this_pft_g1.mean()
        pft_vals = np.append(pft_vals,this_pft_g1_mean)
        this_ef_g1 = g1_by_pft[pft].loc['ef'].dropna()
        this_ef_g1_med = this_ef_g1.median()
        ef_medians = np.append(ef_medians,this_ef_g1_med)
    lin_pft_diff = pft_vals - lin_vals
    lin_pft_diff = np.abs(lin_pft_diff)
    lin_pft_diff_avg = np.nanmean(lin_pft_diff)
    lin_std = np.std(lin_vals)
    lin_ef_diff = ef_medians - lin_vals
    lin_ef_diff = np.abs(lin_ef_diff)
    lin_ef_diff_avg = np.nanmean(lin_ef_diff)
    print('lin_pft_diff_avg: {}'.format(lin_pft_diff_avg))
    print('lin_ef_diff_avg: {}'.format(lin_ef_diff_avg))
    print('lin_std: {}'.format(lin_std))
    p_o = plot_other()
    bin_sizes = [0.5,0.5,0.1,0.1,0.2,0.2]
    num_bins = 20
    for p,pft in enumerate(pfts):
        this_default_g1 = default_g1s[pft]
        this_lin_g1 = lin_g1s[pft]
        this_ef_g1 = g1_by_pft[pft].loc['ef'].dropna()
        this_pft_g1 = g1_by_pft[pft].loc['pft'].dropna()
        this_perc = g1_by_pft[pft].loc['perc'].dropna()
        this_ef_g1_avg = this_ef_g1.mean()
        this_ef_g1_median = this_ef_g1.median()
        this_pft_g1_avg = this_pft_g1.mean()
        this_avg_perc = this_perc.mean()
        tot_pix = len(this_ef_g1)
        print('total pix for {} is {} at {}%'.format(pft,tot_pix,this_avg_perc))
        save_name = os.path.join(
            plots_dir,
            'opt_{}_g1_pdf.svg'
        )
        p_o.plot_pdf(
            [this_ef_g1],
            ['ef_g1'],
            ['#87c56d'],
            [0.5],
            [3],
            [True],
<<<<<<< HEAD
            ['-'],
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
            save_name.format(pft),
            vert_lines=[
                this_ef_g1_median,this_pft_g1_avg,
                this_default_g1,this_lin_g1
            ],
            vert_line_labels=[
                'ef_median','pft_avg',
                'default','lin'
            ],
            vert_line_colors=[
                '#358414','#8e5e96',
                '#080a0b','#ef3f29'
            ],
            xlim=[0,8],
            min_val=[0.5]
        )

    # let's mae PiYG into Pink Grey Green
    old_cmap = plt.get_cmap('PiYG')
    colors = old_cmap(np.linspace(0,1,256))
    # weighted interpolation to bring in grey
    n_colors = len(colors)
    grey = [0.9,0.9,0.9,1.0]
    mid_idx = n_colors // 2
    # modify pink to grey
    for i in range(mid_idx):
        weight = (i/mid_idx)**1
        colors[i] = (1 - weight)*np.array(colors[i]) + weight*np.array(grey)
    # modify grey to green
    for i in range(mid_idx,n_colors):
        weight = ((i - mid_idx)/mid_idx)**1
        colors[i] = (1 - weight)*np.array(grey) + weight*np.array(colors[i])
    new_colors = np.vstack((
        np.linspace(colors[0],grey,128),
        np.linspace(grey,colors[-1],128)
    ))
    PiYG_grey = LinearSegmentedColormap.from_list('PiYG_grey',colors)

    p_p.plot_map(
<<<<<<< HEAD
        'le_perc',pixels,le_perc,
        np.nanmean(le_perc),plots_dir,
        vmin=-50,vmax=50,cmap='PiYG'
    )
    p_p.plot_map(
        'le_perc_drought',pixels,le_perc_drought,
        np.nanmean(le_perc_drought),plots_dir,
        vmin=-50,vmax=50,cmap='PiYG'
    )
=======
        'le_j_diff_ef_pft',pixels,le_j_diff_ef_pft*-1,
        np.nanmean(le_j_diff_ef_pft*-1),plots_dir,
        vmin=-0.5,vmax=0.5,cmap='PiYG'
    )
    #p_w = plot_wat()
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    #p_w.plot_map(
    #    'strm_j_diff_ef_pft',watersheds,np.array(strm_j_diff_ef_pft)*-1,
    #    np.nanmean(np.array(strm_j_diff_ef_pft)*-1),plots_dir,
    #    cmap='PiYG',vmin=-.4,vmax=.4
    #)
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
<<<<<<< HEAD
            'le_j_diff_ef_pft_perc.gpkg'
        ),
        le_perc,
        vmin=-50,
        vmax=50,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'chosen_camels.geojson',
        os.path.join(
            out_dir,
            'strm_j_diff_ef_pft_perc.gpkg'
        ),
        strm_perc,
        vmin=50,
        vmax=50,
        cmap='PiYG',
        subselection=watersheds
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'le_j_diff_ef_def_perc.gpkg'
        ),
        le_perc_def,
        vmin=-50,
        vmax=50,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'chosen_camels.geojson',
        os.path.join(
            out_dir,
            'strm_j_diff_ef_def_perc.gpkg'
        ),
        strm_perc_def,
        vmin=-50,
        vmax=50,
        cmap='PiYG',
        subselection=watersheds
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'le_j_diff_ef_pft_drought.gpkg'
        ),
        le_j_diff_ef_pft_drought*-1,
=======
            'le_j_diff_ef_pft.gpkg'
        ),
        le_j_diff_ef_pft*-1,
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
        vmin=-0.5,
        vmax=0.5,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
<<<<<<< HEAD
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'le_j_diff_ef_def_drought.gpkg'
        ),
        le_j_diff_ef_def_drought*-1,
        vmin=-0.50,
        vmax=0.50,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'chosen_camels.geojson',
        os.path.join(
            out_dir,
            'strm_j_diff_ef_pft_drought.gpkg'
        ),
        strm_j_diff_ef_pft_drought*-1,
=======
        'chosen_camels.geojson',
        os.path.join(
            out_dir,
            'strm_j_diff_ef_pft.gpkg'
        ),
        strm_j_diff_ef_pft*-1,
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
        vmin=-0.5,
        vmax=0.5,
        cmap='PiYG',
        subselection=watersheds
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
<<<<<<< HEAD
        'chosen_camels.geojson',
        os.path.join(
            out_dir,
            'strm_j_diff_ef_def_drought.gpkg'
        ),
        strm_j_diff_ef_def_drought*-1,
        vmin=-0.50,
        vmax=0.50,
=======
        'really_chosen_tiles.geojson',
        os.path.join(
            out_dir,
            'le_j_diff_ef_def.gpkg'
        ),
        le_j_diff_ef_def*-1,
        vmin=-0.5,
        vmax=0.5,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'chosen_camels.geojson',
        os.path.join(
            out_dir,
            'strm_j_diff_ef_def.gpkg'
        ),
        strm_j_diff_ef_def*-1,
        vmin=-0.5,
        vmax=0.5,
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
        cmap='PiYG',
        subselection=watersheds
    )
    # let's make the scatter plot of this
    p_o.scatter(
        strm_j_diff_ef_pft*-1,
        le_j_wat_diff_ef_pft*-1,
        plots_dir,
        'strm_j_diff_ef_pft_vs_le_j_wat_diff_ef_pft',
        'Jstrm PFT - Jstrm EF',
        'JET wat PFT - JET wat EF',
        xlim=[-0.9,0.9],
<<<<<<< HEAD
        ylim=[-0.4,0.4],
        color=all_purple_wat,
        quadrant_lines=True,
        one_to_one_line=True
=======
        ylim=[-0.9,0.9],
        color=all_purple_wat,
        quadrant_lines=True
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    )
    p_o.scatter(
        strm_j_diff_ef_def*-1,
        le_j_wat_diff_ef_def*-1,
        plots_dir,
        'strm_j_diff_ef_def_vs_le_j_wat_diff_ef_def',
        'Jstrm PFT - Jstrm EF',
        'JET wat PFT - JET wat EF',
        xlim=[-0.9,0.9],
<<<<<<< HEAD
        ylim=[-0.4,0.4],
        color=all_purple_wat,
        quadrant_lines=True,
        one_to_one_line=True
=======
        ylim=[-0.9,0.9],
        color=all_purple_wat,
        quadrant_lines=True
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    )
    # let's get the Jstrm and JET for the 6 highest performing basins
    sorted_vals = strm_j_diff_ef_pft.sort_values().head(6)
    sorted_idx = sorted_vals.index
    strm_j_diff_ef_pft_big = strm_j_diff_ef_pft.loc[sorted_idx]
    strm_j_pft_big = strm_j_pft.loc[sorted_idx]
    le_j_wat_pft_big = le_j_wat_pft.loc[sorted_idx]
<<<<<<< HEAD
    strm_j_ef_big = strm_j_ef.loc[sorted_idx]
    le_j_wat_ef_big = le_j_wat_ef.loc[sorted_idx]
    le_j_wat_diff_ef_pft_big = le_j_wat_diff_ef_pft.loc[sorted_idx]
    strm_perc_imp_big = strm_j_diff_ef_pft_big.mean()/strm_j_pft_big.mean()*-100
    le_perc_imp_big = le_j_wat_diff_ef_pft_big.mean()/le_j_wat_pft_big.mean()*-100
    mask = np.isfinite(le_j_wat_ef_big) & np.isfinite(le_j_wat_pft_big)
    _, le_j_diff_ef_pft_big_p = paired_permutation_p(
        le_j_wat_ef_big[mask],
        le_j_wat_pft_big[mask],
        alternative='less'
    )
    mask = np.isfinite(strm_j_ef_big) & np.isfinite(strm_j_pft_big)
    _, strm_j_diff_ef_pft_big_p = paired_permutation_p(
        strm_j_ef_big[mask],
        strm_j_pft_big[mask],
        alternative='less'
    )
=======
    le_j_wat_diff_ef_pft_big = le_j_wat_diff_ef_pft.loc[sorted_idx]
    strm_perc_imp_big = strm_j_diff_ef_pft_big.mean()/strm_j_pft_big.mean()*-100
    le_perc_imp_big = le_j_wat_diff_ef_pft_big.mean()/le_j_wat_pft_big.mean()*-100
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    print(
        'Jstrm improvement in big basins is {}%'.format(strm_perc_imp_big)
    )
    print(
<<<<<<< HEAD
        'p value of {}'.format(strm_j_diff_ef_pft_big_p)
    )
    print(
        'JLE improvement in big basins is {}%'.format(le_perc_imp_big)
    )
    print(
        'p value of {}'.format(le_j_diff_ef_pft_big_p)
    )
=======
        'JLE improvement in big basins is {}%'.format(le_perc_imp_big)
    )
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # and then for the rest of the basins
    strm_j_diff_ef_pft_small = strm_j_diff_ef_pft.drop(sorted_idx)
    strm_j_pft_small = strm_j_pft.drop(sorted_idx)
    le_j_wat_diff_ef_pft_small = le_j_wat_diff_ef_pft.drop(sorted_idx)
    le_j_wat_pft_small = le_j_wat_pft.drop(sorted_idx)
<<<<<<< HEAD
    strm_j_ef_small = strm_j_ef.drop(sorted_idx)
    le_j_wat_ef_small = le_j_wat_ef.drop(sorted_idx)
    strm_perc_imp_small = strm_j_diff_ef_pft_small.mean()/strm_j_pft_small.mean()*-100
    le_perc_imp_small = le_j_wat_diff_ef_pft_small.mean()/le_j_wat_pft_small.mean()*-100
    mask = np.isfinite(le_j_wat_ef_small) & np.isfinite(le_j_wat_pft_small)
    _, le_j_diff_ef_pft_small_p = paired_permutation_p(
        le_j_wat_ef_small[mask],
        le_j_wat_pft_small[mask],
        alternative='greater'
    )
    mask = np.isfinite(strm_j_ef_small) & np.isfinite(strm_j_pft_small)
    _, strm_j_diff_ef_pft_small_p = paired_permutation_p(
        strm_j_ef_small[mask],
        strm_j_pft_small[mask],
        alternative='greater'
    )
=======
    strm_perc_imp_small = strm_j_diff_ef_pft_small.mean()/strm_j_pft_small.mean()*-100
    le_perc_imp_small = le_j_wat_diff_ef_pft_small.mean()/le_j_wat_pft_small.mean()*-100
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    print(
        'Jstrm imp in small basisns is {}%'.format(strm_perc_imp_small)
    )
    print(
<<<<<<< HEAD
        'p value of {}'.format(strm_j_diff_ef_pft_small_p)
    )
    print(
        'JLE imp in small basins is {}%'.format(le_perc_imp_small)
    )
    print(
        'p value of {}'.format(le_j_diff_ef_pft_small_p)
    )
=======
        'JLE imp in small basins is {}%'.format(le_perc_imp_small)
    )
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # okay, now for Figure 5. Let's make the plot of the majority PFT in each
    # different watershed
    wat_pft_cats = []
    wats_reordered = []
    for p,pft in enumerate(list(pft_wats.keys())):
        this_wat = pft_wats[pft]
        if type(this_wat) == str:
            continue
        wats_reordered.extend(this_wat.columns)
        this_pft_rep = [pft for n in range(len(this_wat.columns))]
        wat_pft_cats.extend(this_pft_rep)
    # let's get this back in the right order
    wat_pft_cats_series = pd.Series(
        data=wat_pft_cats,
        index=wats_reordered
    )
    wat_pft_cats_df = pd.DataFrame(
        columns=watersheds
    )
    wat_pft_cats_df.loc['cats'] = wat_pft_cats_series
    wat_pft_cats_df = wat_pft_cats_df.fillna('no_majority')
    wat_pft_cats = wat_pft_cats_df.loc['cats']
    cats_colors_dict = {
        'needleleaf_trees':'#86c36c',
        'broadleaf_trees':'#cb9ac2',
        'c3_grass':'#824c8b',
        'no_majority':'#73c5eb'
    }
    gen.add_cat_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large/outputs/' +
        'chosen_camels.geojson',
        os.path.join(
            out_dir,
            'wat_pft_cats.gpkg'
        ),
        wat_pft_cats,
        cats_colors_dict,
        subselection=watersheds
    )
    # now we need to make the plots that show why these watersheds are special
    # get these 6 watersheds
    num_watersheds = len(watersheds)
    strm_j_diff_ef_pft_max = np.nanmax(strm_j_diff_ef_pft)
    big_improve_idx = np.where(
        np.array(strm_j_diff_ef_pft) < -strm_j_diff_ef_pft_max
    )
    big_improve_wat = watersheds[big_improve_idx]
    # assign them a different color
<<<<<<< HEAD
    all_wat_colors = ['#73c5eb' for i in range(len(watersheds))]
    not_big_improve_wat = np.zeros(0)
    for w,wat in enumerate(watersheds):
        if wat in big_improve_wat:
            all_wat_colors[w] = 'k'
=======
    all_wat_colors = ['k' for i in range(len(watersheds))]
    not_big_improve_wat = np.zeros(0)
    for w,wat in enumerate(watersheds):
        if wat in big_improve_wat:
            all_wat_colors[w] = '#73c5eb'
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    p_o.scatter(
        g1_wat_diff_ef_pft,
        strm_j_pft,
        plots_dir,
        'g1_diff_ef_pft_vs_strm_j_pft',
        'g1 EF - g1 PFT',
        'Jstrm PFT',
        color=all_wat_colors
    )
    # plot PDF of basin-average g1 for EF-optimized
    big_improve_g1s_ef = g1_wat_ef[big_improve_wat]
    p_o.plot_pdf(
        [g1_wat_ef],
        ['g1_watershed'],
        ['#86c148'],
        [1],
        [1],
        [False],
        os.path.join(
            plots_dir,
            'g1_wat_ef_pdf.svg'
        ),
        vert_lines=big_improve_g1s_ef,
        vert_line_labels=[1,2,3,4,5,6],
        vert_line_colors=['k','k','k','k','k','k'],
        min_val=[0.5],
        xlim=[0,8]
    )
    # plot PDF of basin-average g1 for PFT-optimized
    big_improve_g1s_pft = g1_wat_pft[big_improve_wat]
    p_o.plot_pdf(
        [g1_wat_pft],
        ['g1_watershed'],
        ['#86c148'],
        [1],
        [1],
        [False],
        os.path.join(
            plots_dir,
            'g1_wat_pft_pdf.svg'
        ),
        vert_lines=big_improve_g1s_pft,
        vert_line_labels=[1,2,3,4,5,6],
        vert_line_colors=['k','k','k','k','k','k'],
        min_val=[0.5],
        xlim=[0,8]
    )
    # plot pdfs for basin-average precip and canopy
    big_improve_precip = precip_wat[big_improve_wat]
    p_o.plot_pdf(
        [precip_wat],
        ['precip_watershed'],
        ['#86c148'],
        [1],
        [1],
        [False],
        os.path.join(
            plots_dir,
            'precip_wat_pdf.svg'
        ),
        vert_lines=big_improve_precip,
        vert_line_labels=[1,2,3,4,5,6],
        vert_line_colors=['k','k','k','k','k','k']
    )
    big_improve_canopy = canopy_wat[big_improve_wat]
    p_o.plot_pdf(
        [canopy_wat],
        ['canopy_watershed'],
        ['#86c148'],
        [1],
        [1],
        [False],
        os.path.join(
            plots_dir,
            'canopy_wat_pdf.svg'
        ),
        vert_lines=big_improve_canopy,
        vert_line_labels=[1,2,3,4,5,6],
        vert_line_colors=['k','k','k','k','k','k']
    )
    # for the six basins of interest, plot the %need and 
    # %canopy for each of those basins
    print(tile_pft_info_wat_df)
    tile_pft_info_wat_big_df = tile_pft_info_wat_df[
        big_improve_wat
    ]
    need_vals = np.array(
        tile_pft_info_wat_big_df.loc['needleleaf_trees']
    )
    broad_vals = np.array(
        tile_pft_info_wat_big_df.loc['broadleaf_trees']
    )
    vals = np.array([need_vals,broad_vals])
    p_o.multibar_plot(
        vals,
        big_improve_wat,
        ['need','broad'],
        plots_dir,
        'need_broad_big_improve_perc.svg',
        'big improve watersheds',
        '% PFT coverage'
    )















if __name__ == '__main__':
    main()
