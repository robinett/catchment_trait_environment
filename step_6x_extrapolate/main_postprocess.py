import sys
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/pso/step_5.1_analyze_outputs/funcs'
)
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large'
)
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/pso/step_2x_env_covariates'
)
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/pso/step_5_analyze_outputs'
)
import os
import datetime
import numpy as np
import pickle as pickle
import pandas as pd
import netCDF4 as nc
import statsmodels.formula.api as sm
import copy
import matplotlib as mpl
import geopandas as gpd
from tqdm import tqdm
from get_env_covariates import get_covariates
from extrapolation import extrapolate
from general_functions import gen_funcs
from choose_tiles_camels import choose
from get_timeseries import get_timeseries
from get_averages_and_error import averages_and_error
from plot_watersheds import plot_wat
from plot_other import plot_other
from plot_pixels import plot_pixels
from analyze_timeseries_pixel_scale import analyze_pix

<<<<<<< HEAD
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
    # where are we?
    base_dir = '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate'
    # where do we put outputs?
    out_dir = os.path.join(
        base_dir,
        'outputs'
    )
    # where do we put plots?
    plots_dir = os.path.join(
        base_dir,
        'plots'
    )
    # where is the data
    data_dir = os.path.join(
        base_dir,
        'data'
    )
    # where are all saved timeseries kept?
    save_dir = (
        '/discover/nobackup/trobinet/from_aws/pso/' +
        'step_5_analyze_outputs/saved_timeseries'
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
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/' +
        'outputs/extrapolation_pix.csv'
    )
    # where is the le truth located?
    le_truth_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_3x_process_gleam/outputs/' +
        'le_truth_gleam_38a_watts_per_m2_1995-01-01_2014-12-31_' +
        'extrapolation_tiles.csv'
    )
    # where is the tile info?
    raw_pft_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_5_analyze_outputs/' +
        'data/CLM4.5_veg_typs_fracs'
    )
    processed_pft_fname = os.path.join(
        out_dir,
        'processed_pft_info_extrap.csv'
    )
    # where is the tile info?
    tile_pft_info_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_5_analyze_outputs/' +
        'outputs/pft_distribution.csv'
    )
    # where is the precip info for calc g1?
    precip_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/' +
        'outputs/gldas_avg_precip_normalized.nc4'
    )
    # where is the canopy height info for g1 calc?
    canopy_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/' +
        'outputs/canopy_height_normalized.nc4'
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
    # pass a fake intersection info fname because we aren't going to get 
    # it anyways
    intersection_info_fname = 'fake'
    # what are the names and info of the timeseries that we are going to load?
    timeseries_info = {
        'g1-default-extrapolate-1992-2014':{
            'dir':(
                '/discover/nobackup/trobinet/exps/GEOSldas_CN45_med_extrapolate_1992_2014/0'
            ),
            'load_or_save':'load',
            'default_type':'default',
            'read_met_forcing':False,
            'timeseries_dir':os.path.join(
                save_dir,
                'g1-med-extrapolate-1992-2014'
            ),
            'optimization_type':'nan',
            'iteration_number':'nan'
        },
        'g1-ai-extrapolate-1992-2014':{
            'dir':(
                '/discover/nobackup/trobinet/exps/GEOSldas_CN45_ai_et_extrapolate_1992_2014/0'
            ),
            'load_or_save':'load',
            'default_type':'default',
            'read_met_forcing':False,
            'timeseries_dir':os.path.join(
                save_dir,
                'g1-ai-extrapolate-1992-2014'
            ),
            'optimization_type':'nan',
            'iteration_number':'nan'
        },
        'g1-a0-a1-extrapolate-1992-2014':{
            'dir':(
                '/discover/nobackup/trobinet/exps/GEOSldas_CN45_a0_a1_et_extrapolate_1992_2014/0'
            ),
            'load_or_save':'load',
            'default_type':'default',
            'read_met_forcing':False,
            'timeseries_dir':os.path.join(
                save_dir,
                'g1-a0-a1-extrapolate-1992-2014'
            ),
            'optimization_type':'nan',
            'iteration_number':'nan'
        }
    }
    gen = gen_funcs()
    g = get_timeseries()
    # let's start by getting our g1 values
    # get our pixels
    pixels = gen.get_pixels(pixels_fname)
    # get g1 map for pfts
    # first we need the pft breakdown at our pixels of interest
    a_pix = analyze_pix()
    pft_info = a_pix.get_pft_info(raw_pft_fname,processed_pft_fname,pixels)
    # start with positions
    pos_pft = pd.read_csv(positions_pft_fname,header=None)
    pos_pft = np.array(pos_pft)[0]
    # now put it in the weird format so that our g1 decoder can accept it
    position_pft_names = [
        'a0_needleaf_trees',
        'a0_broadleaf_trees',
        'a0_shrub',
        'a0_c3_grass',
        'a0_c4_grass',
        'a0_crop'
    ]
    pos_pft_final = {}
    for p,pft in enumerate(position_pft_names):
        this_pos = pos_pft[p]
        this_pos = [this_pos]
        this_pos = [this_pos]
        pos_pft_final[pft] = this_pos
    # now we need a fake objective value to pass in this weird format
    fake_obj = {}
    fake_obj['all'] = [[1]]
    g1_pft_df = g.get_g1_map(
        pos_pft_final,
        fake_obj,
        precip_fname,
        canopy_fname,
        processed_pft_fname,
        'pft',0
    )
    g1_pft = g1_pft_df.loc['g1']
    # now let's do this for ef
    pos_ef = pd.read_csv(positions_ef_fname,header=None)
    pos_ef = np.array(pos_ef)[0]
    # now put it in the weird format so that our g1 decoder can accept it
    position_ef_names = [
        'b_needleleaf_trees',
        'b_broadleaf_trees',
        'b_shrub',
        'b_c3_grass',
        'b_c4_grass',
        'b_crop',
        'a0_intercept',
        'a1_precip_coef',
        'a2_canopy_coef'
    ]
    pos_ef_final = {}
    for p,pft in enumerate(position_ef_names):
        this_pos = pos_ef[p]
        this_pos = [this_pos]
        this_pos = [this_pos]
        pos_ef_final[pft] = this_pos
    g1_ef_df = g.get_g1_map(
        pos_ef_final,
        fake_obj,
        precip_fname,
        canopy_fname,
        processed_pft_fname,
        'ef',0
    )
    g1_ef = g1_ef_df.loc['g1']
    g1_diff_ef_pft = g1_ef - g1_pft
    # for default g1
    default_g1s_coefs = {
        'a0_needleaf_trees':[[2.3]],
        'a0_broadleaf_trees':[[4.25]],
        'a0_shrub':[[4.7]],
        'a0_c3_grass':[[5.3]],
        'a0_c4_grass':[[1.6]],
        'a0_crop':[[5.79]]
    }
    fake_obj = {}
    fake_obj['all'] = [[1]]
    g1_default_df = g.get_g1_map(
        default_g1s_coefs,
        fake_obj,
        precip_fname,
        canopy_fname,
        processed_pft_fname,
        'pft',
        0
    )
    g1_default = g1_default_df.loc['g1']
<<<<<<< HEAD
    g1_diff_pft_default = g1_pft - g1_default
    g1_diff_ef_default = g1_ef - g1_default
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # make plots of each and their difference
    #p_p = plot_pixels()
    #p_p.plot_map(
    #    'g1_pft_extrap',
    #    pixels,
    #    g1_pft,
    #    np.nanmean(g1_pft),
    #    plots_dir,
    #    vmin=0.5,
    #    vmax=12
    #)
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf',
        os.path.join(
            out_dir,
            'g1_pft_extrap.gpkg'
        ),
        g1_pft,
        vmin=0,
        vmax=7,
        cmap='rainbow'
    )
    #p_p.plot_map(
    #    'g1_ef_extrap',
    #    pixels,
    #    g1_ef,
    #    np.nanmean(g1_ef),
    #    plots_dir,
    #    vmin=0.5,
    #    vmax=12
    #)
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf',
        os.path.join(
            out_dir,
            'g1_ef_extrap.gpkg'
        ),
        g1_ef,
        vmin=0,
        vmax=7,
        cmap='rainbow'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf',
        os.path.join(
            out_dir,
            'g1_default_extrap.gpkg'
        ),
        g1_default,
        vmin=0,
        vmax=7,
        cmap='rainbow'
    )
    # make the pdf for the different g1s
<<<<<<< HEAD
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf',
        os.path.join(
            out_dir,
            'g1_diff_pft_default.gpkg'
        ),
        g1_diff_pft_default,
        vmin=-5,
        vmax=5,
        cmap='PiYG'
    )
    # make the pdf for the different g1s
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf',
        os.path.join(
            out_dir,
            'g1_diff_ef_default.gpkg'
        ),
        g1_diff_ef_default,
        vmin=-5,
        vmax=5,
        cmap='PiYG'
    )
    # make the pdf for the different g1s
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    p_o = plot_other()
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
    line_widths = [1,1,1]
    fills = [False,False,False]
<<<<<<< HEAD
    line_styles = [
        '-','-','-'
    ]
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
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
    #p_p.plot_map(
    #    'g1_diff_ef_pft_extrap',
    #    pixels,
    #    g1_diff_ef_pft,
    #    np.nanmean(g1_diff_ef_pft),
    #    plots_dir,
    #    vmin=-8,
    #    vmax=8,
    #    cmap='PiYG'
    #)
    #p_o = plot_other()
    #p_o.scatter(
    #    g1_pft,
    #    g1_ef,
    #    plots_dir,
    #    'g1_ef_extrap_vs_g1_pft_extrap',
    #    'g1 from PFT for extrap',
    #    'g1 from EF for extrap',
    #    xlim=[0,12],
    #    ylim=[0,12]
    #)
    # and now move on to getting the le values
    # load the timeseries
    exps = list(timeseries_info.keys())
    timeseries_info = g.get_all_timeseries_info(
        timeseries_info,start,end,pixels_fname,
        intersection_info_fname,precip_fname,canopy_fname,
        processed_pft_fname,get_watershed=False
    )
<<<<<<< HEAD
    timeseries_info_drought = copy.deepcopy(timeseries_info)
    timeseries_info_drought = g.get_all_timeseries_info(
        timeseries_info_drought,start,end,pixels_fname,
        intersection_info_fname,precip_fname,canopy_fname,
        tile_pft_info_fname,get_watershed=False
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # get the le obs
    #print('getting le obs')
    le_obs = g.get_le_obs(le_truth_fname)
    ## for now let's just make some fake data
    #old_le_obs = g.get_le_obs(le_truth_fname)
    #print('making fake obs')
    #le_obs = pd.DataFrame(columns=pixels)
    #le_obs['time'] = old_le_obs.index
    #le_obs = le_obs.set_index('time')
    #for t,ti in enumerate(le_obs.index.tolist()):
    #    this_fake_obs = np.random.rand(len(pixels))
    #    le_obs.loc[ti] = this_fake_obs*100
    g_a = averages_and_error()
    for e,exp in enumerate(exps):
        print('getting error for experiment {}'.format(e))
        this_le_err_df = g_a.var_error(
            timeseries_info[exp]['pixel_raw_timeseries']['le'],
<<<<<<< HEAD
            le_obs,start_err,end_err,var_type='daily'
        )
        timeseries_info[exp]['pixel_le_errors'] = this_le_err_df
        this_le_err_df_drought = g_a.var_error(
            timeseries_info_drought[exp]['pixel_raw_timeseries']['le'],
            le_obs,start_err,end_err,during_drought=True,var_type='daily'
        )
        timeseries_info_drought[exp]['pixel_le_errors'] = this_le_err_df_drought
=======
            le_obs,start_err,end_err
        )
        timeseries_info[exp]['pixel_le_errors'] = this_le_err_df
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    print('getting j')
    # get ef j for le and strm
    le_j_ef = (
        timeseries_info[exps[2]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
<<<<<<< HEAD
    le_j_ef_drought = (
        timeseries_info_drought[exps[2]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # get  pft j for le and strm
    le_j_pft = (
        timeseries_info[exps[1]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
<<<<<<< HEAD
    le_j_pft_drought = (
        timeseries_info_drought[exps[1]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    le_j_default = (
        timeseries_info[exps[0]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
<<<<<<< HEAD
    le_j_default_drought = (
        timeseries_info_drought[exps[0]]['pixel_le_errors'].loc['mae']/
        le_obs.mean()
    )
    overall_le_j_pft = np.nanmean(le_j_pft)
    # get differences between j for EF and PFT
    le_j_diff_ef_pft = le_j_ef - le_j_pft
    le_j_diff_ef_pft_drought = le_j_ef_drought - le_j_pft_drought
    le_j_diff_ef_default_drought = le_j_pft - le_j_default
    le_j_ef_mean = le_j_ef.mean()
    le_j_ef_mean_drought = le_j_ef_drought.mean()
    le_j_pft_mean = le_j_pft.mean()
    le_j_pft_mean_drought = le_j_pft_drought.mean()
    le_j_default_mean_drought = le_j_default_drought.mean()
    perc_le_j_diff_ef_pft = (
        (le_j_ef_mean - le_j_pft_mean)/le_j_pft_mean
    )
    perc_le_j_diff_ef_pft_drought = (
        (le_j_ef_mean_drought - le_j_pft_mean_drought)/le_j_pft_mean_drought
    )
    perc_le_j_diff_pft_default_drought = (
        (le_j_pft_mean_drought - le_j_default_mean_drought)/le_j_default_mean_drought
    )
    perc_le_j_diff_ef_pft = perc_le_j_diff_ef_pft*-100
    perc_le_j_diff_ef_pft_drought = perc_le_j_diff_ef_pft_drought*-100
    perc_le_j_diff_pft_default_drought = perc_le_j_diff_pft_default_drought*-100
    # get the significanc3
    mask = np.isfinite(le_j_ef) & np.isfinite(le_j_pft)
    _, le_j_diff_ef_pft_p = paired_permutation_p(
        le_j_ef[mask],
        le_j_pft[mask],
        alternative='greater'
    )
    mask = np.isfinite(le_j_ef_drought) & np.isfinite(le_j_pft_drought)
    _, le_j_diff_ef_pft_p_drought = paired_permutation_p(
        le_j_ef_drought[mask],
        le_j_pft_drought[mask],
        alternative='greater'
    )
    mask = np.isfinite(le_j_pft_drought) & np.isfinite(le_j_default_drought)
    _, le_j_diff_pft_default_p_drought = paired_permutation_p(
        le_j_pft_drought[mask],
        le_j_default_drought[mask],
        alternative='less'
    )
    print(
        'ef is better than pft for extrap by {}%'.format(perc_le_j_diff_ef_pft)
    )
    print(
        'p value of {}'.format(le_j_diff_ef_pft_p)
    )
    print(
        'ef is better than pft for extrap during drought by {}%'.format(perc_le_j_diff_ef_pft_drought)
    )
    print(
        'p value of {}'.format(le_j_diff_ef_pft_p_drought)
    )
    print(
        'pft is better than default for extrap during drought by {}%'.format(perc_le_j_diff_pft_default_drought)
    )
    print(
        'p value of {}'.format(le_j_diff_pft_default_p_drought)
    )
=======
    overall_le_j_pft = np.nanmean(le_j_pft)
    # get differences between j for EF and PFT
    le_j_diff_ef_pft = le_j_ef - le_j_pft
    le_j_ef_mean = le_j_ef.mean()
    le_j_pft_mean = le_j_pft.mean()
    perc_le_j_diff_ef_pft = (
        (le_j_ef_mean - le_j_pft_mean)/le_j_pft_mean
    )
    perc_le_j_diff_ef_pft = perc_le_j_diff_ef_pft*-100
    print(
        'ef is better than pft for extrap by {}%'.format(perc_le_j_diff_ef_pft)
    )
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # get differences between j for EF and default
    le_j_diff_ef_default = le_j_ef - le_j_default
    le_j_ef_mean = le_j_ef.mean()
    le_j_default_mean = le_j_default.mean()
    perc_le_j_diff_ef_default = (
        (le_j_ef_mean - le_j_default_mean)/le_j_default_mean
    )
    perc_le_j_diff_ef_default = perc_le_j_diff_ef_default*-100
<<<<<<< HEAD
    mask = np.isfinite(le_j_ef) & np.isfinite(le_j_default)
    _, le_j_diff_ef_default_p = paired_permutation_p(
        le_j_ef[mask],
        le_j_default[mask],
        alternative='less'
    )
    print(
        'ef is better than default for extrap by {}%'.format(perc_le_j_diff_ef_default)
    )
    print(
        'p value of {}'.format(le_j_diff_ef_default_p)
    )
    le_j_diff_ef_default_drought = le_j_ef_drought - le_j_default_drought
    le_j_ef_mean_drought = le_j_ef_drought.mean()
    le_j_default_mean_drought = le_j_default_drought.mean()
    perc_le_j_diff_ef_default_drought = (
        (le_j_ef_mean_drought - le_j_default_mean_drought)/le_j_default_mean_drought
    )
    perc_le_j_diff_ef_default_drought = perc_le_j_diff_ef_default_drought*-100
    mask = np.isfinite(le_j_ef_drought) & np.isfinite(le_j_default_drought)
    _, le_j_diff_ef_default_p_drought = paired_permutation_p(
        le_j_ef_drought[mask],
        le_j_default_drought[mask],
        alternative='greater'
    )
    print(
        'ef is better than default for extrap during drought by {}%'.format(perc_le_j_diff_ef_default_drought)
    )
    print(
        'p value of {}'.format(le_j_diff_ef_default_p_drought)
    )
=======
    print(
        'ef is better than default for extrap by {}%'.format(perc_le_j_diff_ef_default)
    )
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    #p_p = plot_pixels()
    #p_p.plot_map(
    #    'le_j_pft_extrap',
    #    pixels,
    #    le_j_pft,
    #    np.nanmean(le_j_pft),
    #    plots_dir,
    #    vmin=0,
    #    vmax=1.7
    #)
    #p_p.plot_map(
    #    'le_j_ef_extrap',
    #    pixels,
    #    le_j_ef,
    #    np.nanmean(le_j_ef),
    #    plots_dir,
    #    vmin=0,
    #    vmax=1.7
    #)
    #p_p.plot_map(
    #    'le_j_diff_ef_pft_extrap',
    #    pixels,
    #    le_j_diff_ef_pft,
    #    np.nanmean(le_j_diff_ef_pft),
    #    plots_dir,
    #    vmin=-0.5,
    #    vmax=0.5,
    #    cmap='PiYG'
    #)
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf',
        os.path.join(
            out_dir,
            'le_j_diff_ef_default_extrap.gpkg'
        ),
        le_j_diff_ef_default*-1,
        vmin=-0.5,
        vmax=0.5,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf',
        os.path.join(
            out_dir,
            'le_j_diff_ef_pft_extrap.gpkg'
        ),
        le_j_diff_ef_pft*-1,
        vmin=-0.5,
        vmax=0.5,
        cmap='PiYG'
    )
<<<<<<< HEAD
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf',
        os.path.join(
            out_dir,
            'le_j_diff_ef_default_extrap_drought.gpkg'
        ),
        le_j_diff_ef_default_drought*-1,
        vmin=-0.5,
        vmax=0.5,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf',
        os.path.join(
            out_dir,
            'le_j_diff_ef_pft_extrap_drought.gpkg'
        ),
        le_j_diff_ef_pft_drought*-1,
        vmin=-0.5,
        vmax=0.5,
        cmap='PiYG'
    )
    print(len(le_j_diff_ef_pft))
    # get changes in average le
    #print(timeseries_info[exps[0]])
    le_default_avg = timeseries_info[exps[0]]['pixel_raw_timeseries']['le']
    #print(le_default_avg)
    le_default_avg = le_default_avg.mean(axis=0)
    le_pft_avg = timeseries_info[exps[1]]['pixel_raw_timeseries']['le']
    #print(le_pft_avg)
    le_pft_avg = le_pft_avg.mean(axis=0)
    le_ef_avg = timeseries_info[exps[2]]['pixel_raw_timeseries']['le']
    #print(le_ef_avg)
    le_ef_avg = le_ef_avg.mean(axis=0)
    et_default = le_default_avg/28.94
    et_pft = le_pft_avg/28.94
    et_ef = le_ef_avg/28.94
    et_diff_pft_default = et_pft - et_default
    et_diff_ef_default = et_ef - et_default
    print(len(et_diff_pft_default))
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf',
        os.path.join(
            out_dir,
            'et_diff_pft_default_extrap.gpkg'
        ),
        et_diff_pft_default,
        vmin=-0.6,
        vmax=0.6,
        cmap='PiYG'
    )
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf',
        os.path.join(
            out_dir,
            'et_diff_ef_default_extrap.gpkg'
        ),
        et_diff_ef_default,
        vmin=-0.6,
        vmax=0.6,
        cmap='PiYG'
    )
    sys.exit()
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # let's check by land cover just to see what happens
    # lets start by getting our pft info
    # make a df where each column is a tile and the row is the pft %
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
        '(Corn)':'crop',
        '(Irrigated spring temperate cereal)':'crop'
    }
    # make a df where each column is a tile and the row is the pft %
    tile_pft_info = pd.read_csv(processed_pft_fname)
    tile_pft_info = tile_pft_info.set_index('tile')
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
    # we're also going to need the environmental covariates
    #precip = np.array(
    #    nc.Dataset(
    #        precip_fname
    #    )['vals']
    #)
    #precip_df = pd.DataFrame(columns=pixels)
    #precip_df.loc['precip'] = precip
    #canopy = np.array(
    #    nc.Dataset(
    #        canopy_fname
    #    )['vals']
    #)
    #canopy_df = pd.DataFrame(columns=pixels)
    #canopy_df.loc['canopy'] = canopy
    # now let's eliminate all pixels that are at least 25% croplands
    # report the difference between the PFT-optimized and EF-optimized
    # an plot the map of these pixels
    # identify pixels where croplands > 25%
    #crop_perc = tile_pft_info_df.loc['crop'].to_numpy()
    #print(crop_perc)
    #gen.add_to_gdf_and_save(
    #    '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
    #    'extrap_pixels.gdf',
    #    os.path.join(
    #        out_dir,
    #        'crop_perc_extrap.gpkg'
    #    ),
    #    crop_perc,
    #    vmin=0,
    #    vmax=100,
    #    cmap='rainbow'
    #)
    no_crop_df = tile_pft_info_df.loc['crop'][
        (tile_pft_info_df.loc['crop'] < 10) |
        (tile_pft_info_df.loc['crop'].isna())
    ]
    no_crop_pix = no_crop_df.index.tolist()
    # re-calculate and print the difference with these pixels removed
    le_j_default_no_crop = le_j_default[no_crop_pix]
    le_j_pft_no_crop = le_j_pft[no_crop_pix]
    le_j_ef_no_crop = le_j_ef[no_crop_pix]
    le_j_default_no_crop_mean = le_j_default_no_crop.mean()
    le_j_pft_no_crop_mean = le_j_pft_no_crop.mean()
    le_j_ef_no_crop_mean = le_j_ef_no_crop.mean()
    perc_le_j_diff_ef_pft_no_crop = (
        (le_j_ef_no_crop_mean - le_j_pft_no_crop_mean)/le_j_pft_no_crop_mean
    )
    perc_le_j_diff_ef_pft_no_crop = perc_le_j_diff_ef_pft_no_crop*-100
    print(
        'ef is better than pft without crop for extrap by {}%'.format(
            perc_le_j_diff_ef_pft_no_crop
        )
    )
    perc_le_j_diff_ef_default_no_crop = (
        (le_j_ef_no_crop_mean - le_j_default_no_crop_mean)/le_j_default_no_crop_mean
    )
    perc_le_j_diff_ef_default_no_crop = perc_le_j_diff_ef_default_no_crop*-100
    print(
        'ef is better than default without crop for extrap by {}%'.format(
            perc_le_j_diff_ef_default_no_crop
        )
    )
    # remove those pixels from extrap_pixels.gdf
    extrap_pix_gdf = gpd.read_file(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels.gdf'
    )
    print(extrap_pix_gdf)
    extrap_pix_no_crop_gdf = extrap_pix_gdf.set_index('pixel')
    extrap_pix_no_crop_gdf = extrap_pix_no_crop_gdf.loc[no_crop_pix]
    extrap_pix_no_crop_gdf.to_file(
        (
            '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/' +
            'outputs/extrap_pixels_no_crop.gpkg'
        ),
        driver='GPKG'
    )
    print(extrap_pix_no_crop_gdf)
    # remove those pixels from le_j_diff_ef_pft
    le_j_diff_ef_pft_no_crop = le_j_diff_ef_pft[no_crop_pix]
    # save le_j_diff_ef_pft as gpkg
    gen.add_to_gdf_and_save(
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/outputs/' +
        'extrap_pixels_no_crop.gpkg',
        os.path.join(
            out_dir,
            'le_j_diff_ef_pft_extrap_no_crop.gpkg'
        ),
        le_j_diff_ef_pft_no_crop*-1,
        vmin=-0.5,
        vmax=0.5,
        cmap='PiYG'
    )
    sys.exit()


    ## okay now do the actual analysis
    #pft_pix = {}
    #pft_thresh_perc = 50
    #pfts = list(tile_pft_info_df.index)
    #pfts_plus = copy.deepcopy(pfts)
    #pft_le_j_pft_mean = [0 for k in range(len(pfts_plus))]
    #pft_le_j_ef_mean = [0 for k in range(len(pfts_plus))]
    #pft_le_j_diff_ef_pft_mean = [0 for k in range(len(pfts_plus))]
    #pft_le_j_pft_std = [0 for k in range(len(pfts_plus))]
    #pft_le_j_ef_std = [0 for k in range(len(pfts_plus))]
    #pft_le_j_diff_ef_pft_std = [0 for k in range(len(pfts_plus))]
    #pft_num_pix = [0 for k in range(len(pfts_plus))]
    #pft_g1_std = [0 for k in range(len(pfts_plus))]
    #ef_g1_std = [0 for k in range(len(pfts_plus))]
    #pft_precip_std = [0 for k in range(len(pfts_plus))]
    #pft_canopy_std = [0 for k in range(len(pfts_plus))]
    #pft_g1_median = [0 for k in range(len(pfts_plus))]
    #ef_g1_median = [0 for k in range(len(pfts_plus))]
    #pix_used = []
    #for p,pft in enumerate(pfts_plus):
    #    if pft == 'no_majority':
    #        this_pft_pix = [
    #            pix for pix in pixels if pix not in pix_used
    #        ]
    #    else:
    #        this_pft_pix = tile_pft_info_df.loc[pft][
    #            tile_pft_info_df.loc[pft] > pft_thresh_perc
    #        ].index
    #    this_pft_pix = list(this_pft_pix)
    #    pix_used.extend(this_pft_pix)
    #    this_num_pix = len(this_pft_pix)
    #    if this_num_pix == 0:
    #        pft_pix[pft] = 'none'
    #        continue
    #    this_pft_pix_df = pd.DataFrame(columns=this_pft_pix)
    #    this_pft_pix_df.loc['le_j_pft'] = le_j_pft[this_pft_pix]
    #    this_pft_pix_df.loc['le_j_ef'] = le_j_ef[this_pft_pix]
    #    this_pft_pix_df.loc['le_j_diff_ef_pft'] = (
    #        le_j_diff_ef_pft[this_pft_pix]
    #    )
    #    this_pft_pix_df.loc['pft_g1'] = g1_pft[this_pft_pix]
    #    this_pft_pix_df.loc['ef_g1'] = g1_ef[this_pft_pix]
    #    this_pft_pix_df.loc['precip'] = precip_df[this_pft_pix].loc['precip']
    #    this_pft_pix_df.loc['canopy'] = canopy_df[this_pft_pix].loc['canopy']
    #    pft_pix[pft] = this_pft_pix_df
    #    pft_le_j_pft_mean[p] = (
    #        this_pft_pix_df.loc['le_j_pft'].mean()
    #    )
    #    pft_le_j_ef_mean[p] = (
    #        this_pft_pix_df.loc['le_j_ef'].mean()
    #    )
    #    pft_le_j_diff_ef_pft_mean[p] = (
    #        this_pft_pix_df.loc['le_j_diff_ef_pft'].mean()
    #    )
    #    pft_le_j_pft_std[p] = (
    #        this_pft_pix_df.loc['le_j_pft'].std()
    #    )
    #    pft_le_j_ef_std[p] = (
    #        this_pft_pix_df.loc['le_j_ef'].std()
    #    )
    #    pft_le_j_diff_ef_pft_std[p] = (
    #        this_pft_pix_df.loc['le_j_diff_ef_pft'].std()
    #    )
    #    pft_num_pix[p] = this_num_pix
    #    pft_g1_std[p] = (
    #        this_pft_pix_df.loc['pft_g1'].std()
    #    )
    #    ef_g1_std[p] = (
    #        this_pft_pix_df.loc['ef_g1'].std()
    #    )
    #    pft_precip_std[p] = (
    #        this_pft_pix_df.loc['precip'].std()
    #    )
    #    pft_canopy_std[p] = (
    #        this_pft_pix_df.loc['canopy'].std()
    #    )
    #    pft_g1_median[p] = (
    #        this_pft_pix_df.loc['pft_g1'].median()
    #    )
    #    ef_g1_median[p] = (
    #        this_pft_pix_df.loc['ef_g1'].median()
    #    )
    #pft_le_summary = pd.DataFrame(columns=pfts_plus)
    #pft_le_summary.loc['le_j_pft_mean'] = pft_le_j_pft_mean
    #pft_le_summary.loc['le_j_ef_mean'] = pft_le_j_ef_mean
    #pft_le_summary.loc['le_j_diff_ef_pft_mean'] = (
    #    pft_le_j_diff_ef_pft_mean
    #)
    #pft_le_summary.loc['le_j_pft_std'] = pft_le_j_pft_std
    #pft_le_summary.loc['le_j_ef_std'] = pft_le_j_ef_std
    #pft_le_summary.loc['le_j_diff_ef_pft_std'] = (
    #    pft_le_j_diff_ef_pft_std
    #)
    #pft_le_summary.loc['num_pix'] = pft_num_pix
    #pft_le_summary.loc['pft_g1_std'] = pft_g1_std
    #pft_le_summary.loc['ef_g1_std'] = ef_g1_std
    #pft_le_summary.loc['precip_std'] = pft_precip_std
    #pft_le_summary.loc['canopy_std'] = pft_canopy_std
    #pft_le_summary.loc['pft_g1_median'] = pft_g1_median
    #pft_le_summary.loc['ef_g1_median'] = ef_g1_median
    ## okay great, now we can make the bar plots
    ## first for the standard deviation of g1 in each of the pfts
    #vals = np.array(
    #    [
    #        pft_le_summary.loc['pft_g1_std'],
    #        pft_le_summary.loc['ef_g1_std']
    #    ]
    #)
    #categories = list(pft_le_summary.columns)
    #annotation = [0 for k in range(len(pfts_plus))]
    #for p,pft in enumerate(pfts_plus):
    #    annotation[p] = 'n={}'.format(
    #        np.round(pft_le_summary[pft].loc['num_pix'])
    #    )
    #p_o.multibar_plot(
    #    vals,categories,['PFT','EF'],plots_dir,
    #    'bar_pft_g1_std_extrap_vs_ef_g1_std_extrap',
    #    'PFTs','std(g1)',annotation=annotation
    #)
    ## and for the standard deviation of normalized environmental predictors for
    ## each of the pfts
    #vals = np.array(
    #    [
    #        pft_le_summary.loc['precip_std'],
    #        pft_le_summary.loc['canopy_std']
    #    ]
    #)
    #categories = list(pft_le_summary.columns)
    #annotation = [0 for k in range(len(pfts_plus))]
    #for p,pft in enumerate(pfts_plus):
    #    annotation[p] = 'n={}'.format(
    #        np.round(pft_le_summary[pft].loc['num_pix'])
    #    )
    #p_o.multibar_plot(
    #    vals,categories,['precip','canopy'],plots_dir,
    #    'bar_precip_std_extrap_vs_canopy_std_extrap',
    #    'PFTs','std(norm(env)',annotation=annotation
    #)
    ## and finally just for J_LE for both
    #vals = np.array(
    #    [
    #        pft_le_summary.loc['le_j_pft_mean'],
    #        pft_le_summary.loc['le_j_ef_mean']
    #    ]
    #)
    #categories = list(pft_le_summary.columns)
    #annotation = [0 for k in range(len(pfts_plus))]
    #for p,pft in enumerate(pfts_plus):
    #    annotation[p] = 'n={}'.format(
    #        np.round(pft_le_summary[pft].loc['num_pix'])
    #    )
    #p_o.multibar_plot(
    #    vals,categories,['PFT','EF'],plots_dir,
    #    'bar_le_j_pft_mean_extrap_vs_le_j_pft_mean_extrap',
    #    'PFTs','J_LE',annotation=annotation
    #)
    ## let's just make a map of the majority pixels too
    #pix_pft_cats = []
    #pix_reordered = []
    #counter = 0
    #for p,pft in enumerate(list(pft_pix.keys())):
    #    this_pix = pft_pix[pft]
    #    pix_reordered.extend(this_pix.columns)
    #    this_pft_rep = [counter for n in range(len(this_pix.columns))]
    #    pix_pft_cats.extend(this_pft_rep)
    #    counter += 1
    #    print(pft)
    #    print(counter)
    #p_p.plot_map(
    #    'pixel_pft_cats',
    #    pix_reordered,
    #    pix_pft_cats,
    #    np.nanmean(pix_pft_cats),
    #    plots_dir
    #)
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
    print(g1_by_pft)
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
    # let's make a pdf of the distribution for each pft

    # first we need the Lin and default values
    lin_g1s = {
        'needleleaf_trees':2.35,
        'broadleaf_trees':3.97,
        'shrub':4.22,
        'c4_grass':1.62,
        'c3_grass':4.5,
        'crop':5.79
    }
    default_g1s = {
        'needleleaf_trees':2.3,
        'broadleaf_trees':4.25,
        'shrub':4.7,
        'c3_grass':5.3,
        'c4_grass':1.6,
        'crop':5.79
    }
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
            'ext_{}_g1_pdf.svg'
        )
        p_o.plot_pdf(
            [this_ef_g1],
            ['ef_g1'],
            ['#87c56d'],
            [0.5],
            [3],
            [True],
            save_name.format(pft),
            fill_alphas=[0.5],
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
    # let's make a histogram of the distribution for each of the PFTs
    p_o = plot_other()
    bin_sizes = [0.5,0.5,0.1,0.1,0.2,0.2]
    num_bins = 20
    for p,pft in enumerate(pfts):
        this_ef_g1 = g1_by_pft[pft].loc['ef'].dropna()
        this_pft_g1 = g1_by_pft[pft].loc['pft'].dropna()
        this_perc = g1_by_pft[pft].loc['perc'].dropna()
        this_ef_g1_avg = this_ef_g1.mean()
        this_ef_g1_median = this_ef_g1.median()
        this_pft_g1_avg = this_pft_g1.mean()
        this_avg_perc = this_perc.mean()
        tot_pix = len(this_ef_g1)
        print('total pix for {} is {} at {}%'.format(pft,tot_pix,this_avg_perc))
        p_o.histogram(
            this_ef_g1,
            plots_dir,
            'g1_ef_hist_{}.png'.format(pft),
            x_label='g1 (sqrt(kPa))',
            y_label='num pixels',
            #bin_width=bin_sizes[p],
            num_bins=num_bins,
            vert_lines=[this_ef_g1_median,this_pft_g1_avg],
            vert_line_labels=['ef_median','pft_avg'],
            vert_line_colors=['#86c36c','#cb9ac2']
        )



    sys.exit()





if __name__ == '__main__':
    main()



