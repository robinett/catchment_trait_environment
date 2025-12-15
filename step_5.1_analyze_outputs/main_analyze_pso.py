import sys
sys.path.append('/shared/pso/step_5.1_analyze_outputs/funcs')
sys.path.append('/shared/pso/step_5.1_analyze_outputs/tasks')
import os
import datetime
import numpy as np
import pickle as pickle
from get_timeseries import get_timeseries
from rainfall_runoff_ratio import rainfall_runoff
from plot_timeseries import timeseries
from get_averages_and_error import averages_and_error
from plot_watersheds import plot_wat
from plot_other import plot_other
from general_functions import gen_funcs
from plot_pixels import plot_pixels

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
            'optimization_type':'nan'
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
            'optimization_type':'pft'
        },
        'g1-ai-et-strm-spin19921994-test19952014-mae-num5':{
            'dir':(
                '/shared/pso_outputs/g1_ai_et_strm_camels_spin19921994_test19952914_mae/num_5/10'
            ),
            'load_or_save':'load',
            'default_type':'pso_output',
            'read_met_forcing':False,
            'timeseries_dir':os.path.join(
                save_dir,
                'g1-ai-et-strm-spin19921994-test19952014-mae-num5'
            ),
            'optimization_type':'pft'
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
            'optimization_type':'ef'
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
            'optimization_type':'ef'
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
        intersection_info_fname
    )
    g_a = averages_and_error()
    timeseries_info = g_a.get_all_averages_and_error(
        timeseries_info,le_obs,strm_obs,start_err,end_err
    )
    # cool we have the information we need

    ##########################################
    #  The section where we plot timeseries  #
    ##########################################
    #p_t = timeseries()
    ## le timeseries for PFT
    #to_plot = [
    #    timeseries_info[to_load[0]]['pixel_raw_timeseries'],
    #    timeseries_info[to_load[1]]['pixel_raw_timeseries'],
    #    timeseries_info[to_load[2]]['pixel_raw_timeseries'],
    #    {'le':le_obs}
    #]
    #names = [
    #    'default',
    #    'iteration 1',
    #    'iteration 5',
    #    'GLEAM'
    #]
    #sizes = np.repeat(0.5,len(names))
    #colors = ['r','y','g','k']
    #start_plot = datetime.date(2010,1,1)
    #end_plot = datetime.date(2011,12,31)
    #p_t.plot_one_var(
    #    to_plot,names,'le','W/m2','timeseries_pft_mae_it1_it5',
    #    plots_dir,start=start_plot,end=end_plot,
    #    small_preds=sizes,colors=colors
    #)
    ## streamflow timeseries for PFT
    #to_plot = [
    #    timeseries_info[to_load[0]]['wat_raw_timeseries'],
    #    timeseries_info[to_load[1]]['wat_raw_timeseries'],
    #    timeseries_info[to_load[2]]['wat_raw_timeseries'],
    #    {'strm_yr':strm_obs}
    #]
    #names = [
    #    'default',
    #    'iteration 1',
    #    'iteration 5',
    #    'CAMELS'
    #]
    #sizes = np.repeat(0.5,len(names))
    #p_t.plot_one_var(
    #    to_plot,names,'strm_yr','mm/day','timeseries_pft_mae_it1_it5',
    #    plots_dir,start=start,end=end,
    #    small_preds=sizes,colors=colors
    #)
    ## le timeseries for EF
    #to_plot = [
    #    timeseries_info[to_load[0]]['pixel_raw_timeseries'],
    #    timeseries_info[to_load[3]]['pixel_raw_timeseries'],
    #    timeseries_info[to_load[4]]['pixel_raw_timeseries'],
    #    {'le':le_obs}
    #]
    #names = [
    #    'default',
    #    'iteration 1',
    #    'iteration 6',
    #    'GLEAM'
    #]
    #sizes = np.repeat(0.5,len(names))
    #start_plot = datetime.date(2010,1,1)
    #end_plot = datetime.date(2011,12,31)
    #p_t.plot_one_var(
    #    to_plot,names,'le','W/m2','timeseries_ef_mae_it1_it6',
    #    plots_dir,start=start_plot,end=end_plot,
    #    small_preds=sizes,colors=colors
    #)
    ## streamflow timeseries for EF
    #to_plot = [
    #    timeseries_info[to_load[0]]['wat_raw_timeseries'],
    #    timeseries_info[to_load[3]]['wat_raw_timeseries'],
    #    timeseries_info[to_load[4]]['wat_raw_timeseries'],
    #    {'strm_yr':strm_obs}
    #]
    #names = [
    #    'default',
    #    'iteration 1',
    #    'iteration 6',
    #    'CAMELS'
    #]
    #sizes = np.repeat(0.5,len(names))
    #p_t.plot_one_var(
    #    to_plot,names,'strm_yr','mm/day','timeseries_ef_mae_it1_it6',
    #    plots_dir,start=start,end=end,
    #    small_preds=sizes,colors=colors
    #)
    ## lets plot pft, ef, and obs on the same timeseries
    ## for le
    #to_plot = [
    #    timeseries_info[to_load[2]]['pixel_raw_timeseries'],
    #    timeseries_info[to_load[4]]['pixel_raw_timeseries'],
    #    {'le':le_obs}
    #]
    #names = [
    #    'PFT itertaion 5',
    #    'EF iteration 6',
    #    'GLEAM'
    #]
    #sizes = np.repeat(0.5,len(names))
    #colors = ['r','g','k']
    #start_plot = datetime.date(2010,1,1)
    #end_plot = datetime.date(2011,12,31)
    #p_t.plot_one_var(
    #    to_plot,names,'le','W/m2','timeseries_ef_and_pft',
    #    plots_dir,start=start_plot,end=end_plot,
    #    small_preds=sizes,colors=colors
    #)
    ## for strm
    #to_plot = [
    #    timeseries_info[to_load[2]]['wat_raw_timeseries'],
    #    timeseries_info[to_load[4]]['wat_raw_timeseries'],
    #    {'strm_yr':strm_obs}
    #]
    #names = [
    #    'PFT iteration 5',
    #    'EF iteration 6',
    #    'CAMELS'
    #]
    #sizes = np.repeat(0.5,len(names))
    #p_t.plot_one_var(
    #    to_plot,names,'strm_yr','mm/day','timeseries_ef_and_pft',
    #    plots_dir,start=start,end=end,
    #    small_preds=sizes,colors=colors
    #)

    ##############################################
    #  The section where we plot watershed maps  #
    ##############################################
    # initiate the class
    #p_w = plot_wat()
    # get the average stream obs, which we will need for a ton of these
    # calculations
    avg_strm = strm_obs.mean()
    ## get the hucs, which we will need to specify in all maps
    #hucs = list(strm_obs.columns)
    #hucs = [
    #    int(h) for h in hucs
    #]
    ## let's plot the change in J value from iteration 1 to 5 for streamflow for
    ## PFT
    j_strm_1_pft = (
        timeseries_info[exps[1]]['wat_strm_errors'].loc['mae']/
        avg_strm
    )
    j_strm_5_pft = (
        timeseries_info[exps[2]]['wat_strm_errors'].loc['mae']/
        avg_strm
    )
    change_j_strm_pft = j_strm_5_pft - j_strm_1_pft
    change_j_strm_pft = np.array(change_j_strm_pft)
    avg_change_j_strm_pft = np.mean(change_j_strm_pft)
    #p_w.plot_map(
    #    'change_j_strm_pft',hucs,change_j_strm_pft,
    #    avg_change_j_strm_pft,plots_dir,vmin=-0.45,vmax=0.45,
    #    cmap='bwr'
    #)
    # let's plot the change in J value from iteration 1 to 5 for streamflow for
    # EF
    j_strm_1_ef = (
        timeseries_info[exps[3]]['wat_strm_errors'].loc['mae']/
        avg_strm
    )
    j_strm_6_ef = (
        timeseries_info[exps[4]]['wat_strm_errors'].loc['mae']/
        avg_strm
    )
    change_j_strm_ef = j_strm_6_ef - j_strm_1_ef
    change_j_strm_ef = np.array(change_j_strm_ef)
    avg_change_j_strm_ef = np.mean(change_j_strm_ef)
    #p_w.plot_map(
    #    'change_j_strm_ef',hucs,change_j_strm_ef,
    #    avg_change_j_strm_ef,plots_dir,vmin=-0.25,vmax=0.25,
    #    cmap='bwr'
    #)
    # let's plot the difference in objective values between PFT and EF
    # optimizations
    diff_j_strm_ef_pft = j_strm_6_ef - j_strm_5_pft
    diff_j_strm_ef_pft = np.array(diff_j_strm_ef_pft)
    avg_diff_j_strm_ef_pft = np.nanmean(diff_j_strm_ef_pft)
    print(np.nanmin(diff_j_strm_ef_pft))
    print(np.nanmax(diff_j_strm_ef_pft))
    print(avg_diff_j_strm_ef_pft)
    #p_w.plot_map(
    #    'diff_j_strm_ef_pft',hucs,diff_j_strm_ef_pft,
    #    avg_diff_j_strm_ef_pft,plots_dir,vmin=-.7,vmax=.7,
    #    cmap='bwr'
    #)
    ## let's plot the % change compared to PFT at each pixel
    #perc_diff_j_strm_ef_pft = diff_j_strm_ef_pft/j_strm_5_pft
    #perc_diff_j_strm_ef_pft = np.array(perc_diff_j_strm_ef_pft)
    #avg_perc_diff_j_strm_ef_pft = np.nanmean(perc_diff_j_strm_ef_pft)
    #p_w.plot_map(
    #    'perc_diff_j_strm_ef_pft',hucs,perc_diff_j_strm_ef_pft,
    #    avg_perc_diff_j_strm_ef_pft,plots_dir,
    #    cmap='bwr')

    ##########################################
    #  The section where we plot pixel maps  #
    ##########################################
    p_p = plot_pixels()
    # get the average le, which we will need for a ton of the calculations
    avg_et = le_obs.mean()
    # get the pixels, which we will need for a ton of the calculations
    gen = gen_funcs()
    pixels = gen.get_pixels(pixels_fname)
    # let's plot the change in J value from iteration 1 to 5 for streamflow for
    # PFT
    j_et_1_pft = (
        timeseries_info[exps[1]]['pixel_le_errors'].loc['mae']/
        avg_et
    )
    j_et_5_pft = (
        timeseries_info[exps[2]]['pixel_le_errors'].loc['mae']/
        avg_et
    )
    change_j_et_pft = j_et_5_pft - j_et_1_pft
    change_j_et_pft = np.array(change_j_et_pft)
    avg_change_j_et_pft = np.mean(change_j_et_pft)
    #p_p.plot_map(
    #    'change_j_et_pft',pixels,change_j_et_pft,
    #    avg_change_j_et_pft,plots_dir,vmin=-0.08,vmax=0.08,
    #    cmap='bwr'
    #)
    # let's plot the change in J value from iteration 1 to 5 for streamflow for
    # EF
    j_et_1_ef = (
        timeseries_info[exps[3]]['pixel_le_errors'].loc['mae']/
        avg_et
    )
    j_et_6_ef = (
        timeseries_info[exps[4]]['pixel_le_errors'].loc['mae']/
        avg_et
    )
    change_j_et_ef = j_et_6_ef - j_et_1_ef
    change_j_et_ef = np.array(change_j_et_ef)
    avg_change_j_et_ef = np.mean(change_j_et_ef)
    #p_p.plot_map(
    #    'change_j_et_ef',pixels,change_j_et_ef,
    #    avg_change_j_et_ef,plots_dir,vmin=-0.12,vmax=0.12,
    #    cmap='bwr'
    #)
    # let's plot the difference in objective values between PFT and EF
    # optimizations
    diff_j_et_ef_pft = j_et_6_ef - j_et_5_pft
    diff_j_et_ef_pft = np.array(diff_j_et_ef_pft)
    avg_diff_j_et_ef_pft = np.nanmean(diff_j_et_ef_pft)
    #p_p.plot_map(
    #    'diff_j_et_ef_pft',pixels,diff_j_et_ef_pft,
    #    avg_diff_j_et_ef_pft,plots_dir,vmin=-.18,vmax=.18,
    #    cmap='bwr'
    #)
    # let's plot the % difference in objective values between PFT and EF
    perc_diff_j_et_ef_pft = diff_j_et_ef_pft/j_et_5_pft
    perc_diff_j_et_ef_pft = np.array(perc_diff_j_et_ef_pft)
    print(np.nanmax(perc_diff_j_et_ef_pft))
    print(np.nanmin(perc_diff_j_et_ef_pft))
    avg_perc_diff_j_et_ef_pft = np.nanmean(perc_diff_j_et_ef_pft)
    print(avg_perc_diff_j_et_ef_pft)
    p_p.plot_map(
        'perc_diff_j_et_ef_pft',pixels,perc_diff_j_et_ef_pft,
        avg_perc_diff_j_et_ef_pft,plots_dir,
        cmap='bwr'
    )

    #######################################################
    #  The section where we plot movements per iteration  #
    #######################################################
    #p_o = plot_other()
    #pft_obj = timeseries_info[to_load[2]]['obj_vals']['all']
    #pft_obj_strm = timeseries_info[to_load[2]]['obj_vals']['strm']
    #pft_obj_et = timeseries_info[to_load[2]]['obj_vals']['et']
    #p_o.iteration_plot(
    #    pft_obj,'objective value',plots_dir,
    #    'pft_obj_val_iteration.png'
    #)
    #p_o.iteration_plot(
    #    pft_obj_strm,'strm objective value',plots_dir,
    #    'pft_strm_obj_val_iteration.png'
    #)
    #p_o.iteration_plot(
    #    pft_obj_et,'objective value',plots_dir,
    #    'pft_et_obj_val_iteration.png'
    #)
    ## make the plot of change in objective value for EF
    #ef_obj = timeseries_info[to_load[4]]['obj_vals']['all']
    #ef_obj_strm = timeseries_info[to_load[4]]['obj_vals']['strm']
    #ef_obj_et = timeseries_info[to_load[4]]['obj_vals']['et']
    #pft_a0_needleleaf = timeseries_info[to_load[2]]['positions']['a0_needleaf_trees']
    #pft_a0_crop = timeseries_info[to_load[2]]['positions']['a0_crop']
    #ef_a1 = timeseries_info[to_load[4]]['positions']['a1_precip_coef']
    #ef_a2 = timeseries_info[to_load[4]]['positions']['a2_canopy_coef']
    #p_o.iteration_plot(
    #    ef_obj,'objective value',plots_dir,
    #    'ef_obj_val_iteration.png'
    #)
    #p_o.iteration_plot(
    #    ef_obj_strm,'strm objective value',plots_dir,
    #    'ef_strm_obj_val_iteration.png'
    #)
    #p_o.iteration_plot(
    #    ef_obj_et,'objective value',plots_dir,
    #    'ef_et_obj_val_iteration.png'
    #)
    #p_o.iteration_plot(
    #    pft_a0_needleleaf,'a0 needleleaf trees',plots_dir,
    #    'pft_a0_needleleaf_iteration.png'
    #)
    #p_o.iteration_plot(
    #    pft_a0_crop,'a0 crop',plots_dir,
    #    'pft_a0_crop_iteration.png'
    #)
    #p_o.iteration_plot(
    #    ef_a1,'a1 precip coef',plots_dir,
    #    'ef_a1_iteration.png'
    #)
    #p_o.iteration_plot(
    #    ef_a2,'a2 precip coef',plots_dir,
    #    'ef_a2_iteration.png'
    #)

    ##################################################################
    #  Checking the importance of small vs. large watersheds in PSO  #
    ##################################################################
    ## let's check how big the change in j_strm_pft is versus the size of the
    ## watershed
    #abs_change_j_strm_pft = np.abs(change_j_strm_pft)
    #p_o.scatter(
    #    avg_strm,abs_change_j_strm_pft,plots_dir,'change_j_strm_pft_vs_avg_strm',
    #    'average basin streamflow (mm/day)','change J strm PFT',
    #    xlim=[-.5,10],ylim=[-.02,.37]
    #)
    #abs_change_j_strm_ef = np.abs(change_j_strm_ef)
    #p_o.scatter(
    #    avg_strm,abs_change_j_strm_ef,plots_dir,'change_j_strm_ef_vs_avg_strm',
    #    'average basin streamflow (mm/day)','change J strm EF',
    #    xlim=[-.5,10],ylim=[-.02,.37]
    #)


if __name__ == '__main__':
    main()
