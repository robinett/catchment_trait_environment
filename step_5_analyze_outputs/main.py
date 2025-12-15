from get_timeseries import get_timeseries
from analyze_timeseries_pixel_scale import analyze_pix
from analyze_timeseries_watershed_scale import analyze_watershed
from analyze_pso import analyze_pso
from compare_experiments import compare_experiments
from drought_analysis import drought
import datetime
import os
import numpy as np
import sys

def main():
    # highest level directory for this step
    base_dir = '/shared/pso/step_5_analyze_outputs/'
    # lets get this information at all of the pixels from step 1. List the step
    # 1 pixels directory here
    step_1_pixels_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/intersecting_catch_tiles.csv'
    )
    # where to save the catchment timeseries generated
    timeseries_dir = os.path.join(
        base_dir,'saved_timeseries'
    )
    # where to save plots
    plots_dir = os.path.join(
        base_dir,'plots'
    )
    # directory of miscalaneous data to be used
    data_dir = os.path.join(
        base_dir,'data'
    )
    # where to save outputs
    out_dir = os.path.join(
        base_dir,'outputs'
    )
    # where raw PFT information exists
    raw_pft_dir = os.path.join(
        data_dir,'CLM4.5_veg_typs_fracs'
    )
    # where analyzed pft info exists
    analyzed_pft_dir = os.path.join(
        out_dir,'pft_distribution.csv'
    )
    # the directory for comparison run
    default_experiment = (
        '/lustre/catchment/exps/GEOSldas_CN45_med_default_pft_g1_2000_2014_camels/0'
    )
    # the directory for the first iteration of the PSO run
    first_iteration_experiment = (
        '/shared/pso_outputs/g1_a0_a1_et_camels_spin2006_test2007_v2/num_1/4'
    )
    # the directory for the final iteration of the PSO run
    final_iteration_experiment = (
        '/shared/pso_outputs/g1_a0_a1_et_camels_spin2006_test2007_v2/num_17/9'
    )
    # set the start date and end date for our period of interest
    # start date and ned date are inclusive
    start = datetime.date(2000,1,1)
    end = datetime.date(2014,12,31)
    # start and end dates for when to compute eror statistics
    start_error = datetime.date(2003,1,1)
    end_error = datetime.date(2014,12,31)
    # let's set the experiments that we need to load
    # the name, dir, load_bools, and save_bools below are all lists where each
    # element of that list corresponds to:
        # 0: the default, medlyn experiment
        # 1: the first iteration of the PSO run
        # 2: the final iteration of the PSO run
    # what are the experiment names to use
    exp_names = [
        'med-mbb-default-pft-g1-2000-2014',
        'med-mbb-default-pft-g1-2000-2014',
        'med-mbb-default-pft-g1-2000-2014'
    ]
    # location of the all_info.pkl file containing pso information
    all_info_fname = (
        '/shared/pso_outputs/g1_a0_a1_et_camels_spin2006_test2007_v2/all_info.pkl'
    )
    # what are the directories where these experiments exist
    exp_dirs = [
        default_experiment,first_iteration_experiment,final_iteration_experiment
    ]
    # Are these default experiments? A defualt experiment is one that is still
    # in the place it was created (not copied to a different directory during
    # PSO optimization)?
    is_default_experiment = [True,False,False]
    # do these experiments have met forcing data associated with them?
    read_met_forcing = [False,False,False]
    # if raw model outputs have already been processed, 
    # should processed data be loaded instead?
    load_bools =  [False,True,True]
    # if generating a new timeseries, should it be saved?
    save_bools = [True,False,False]
    # default experiment must be filled in as a real directory the other
    # two must be fake and the experiments must be pre-saved as .pkl files
    default_experiment_compare = default_experiment
    # the directory for the first iteration of the PSO run
    first_iteration_experiment_compare = 'fake'
    # the directory for the final iteration of the PSO run
    final_iteration_experiment_compare = 'fake'
    exp_dirs_compare = [
        default_experiment_compare,first_iteration_experiment_compare,
        final_iteration_experiment_compare
    ]
    exp_names_compare = [
        'med-mbb-default-pft-g1',
        'g1-ai-et-camels-v3-spin2006-test2007-first',
        'g1-ai-et-camels-v3-spin2006-test2007-final'
    ]
    save_bools_compare = [False,False,False]
    load_bools_compare = [True,True,True]
    is_default_experiment_compare = [True,False,False]
    read_met_forcing_compare = [False,False,False]
    # what is the all info for this comparison experiment?
    comparison_all_info_fname = (
        '/shared/pso_outputs/g1_ai_et_camels_spin2006_test2007_v3/all_info.pkl'
    )
    # is this a pso experiment? obviously, yes it is
    experiment_type = 'pso'
    # where is the fluxcom data located?
    fluxcom_dir = (
        '/shared/pso/step_3_process_fluxcom/outputs/' +
        'le_truth_fluxcom_rs_CRUNCEP_ensemble_watts_per_m2_' +
        '2006-01-01_2007-12-31_camels_tiles.csv'
    )
    # where is the streamflow data located?
    strm_dir = (
        '/shared/pso/step_3.1.1_process_camels/outputs/' +
        'camels_truth_2006-01-01_2007-12-31_mm_day.csv'
    )
    # where is the dictionary containing information about which tiles
    # intersected which hucs?
    intersection_dict_dir = (
        '/shared/pso/step_1_choose_tiles/outputs/intersection_info.pkl'
    )
    # what is the base dir for covariate maps?t
    env_covariate_base_dir = (
        '/shared/pso/step_2_env_covariates/outputs'
    )
    # where is precip covariate map located?
    precip_map_fname = os.path.join(
        env_covariate_base_dir,
        'averaged_precip.nc4'
    )
    # where is lai covariate map located?
    lai_map_fname = os.path.join(
        env_covariate_base_dir,
        'averaged_lai.nc4'
    )
    # where is sand fraction covariate map located?
    sand_map_fname = os.path.join(
        env_covariate_base_dir,
        'averaged_sand.nc4'
    )
    # where is the k0 covariate map located?
    k0_map_fname = os.path.join(
        env_covariate_base_dir,
        'averaged_ksat.nc4'
    )
    # where are the error comparison maps located?
    error_comparisons_dir = (
        '/shared/pso/other_analyses/errors_vs_covariates/outputs'
    )
    # where is the gldas precip located?
    gldas_precip_dir = os.path.join(
        error_comparisons_dir,
        'gldas_precip_kg_m2_s.csv'
    )
    # where is the gldas temp loated?
    gldas_temp_dir = os.path.join(
        error_comparisons_dir,
        'gldas_air_temp_K.csv'
    )
    # where is the gldas potential ET located?
    gldas_pet_dir = os.path.join(
        error_comparisons_dir,
        'gldas_potential_et_W_m2.csv'
    )
    # where is the canopy height located?
    canopy_height_dir = os.path.join(
        error_comparisons_dir,
        'canopy_height.csv'
    )
    # where is the geojson file with the information for the huc6 located?
    geojson_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/chosen_camels.geojson'
    )
    # where is the directory with shapefile with information on the US states?
    states_shp_fname = (
        '/shared/pso/step_1_choose_tiles/data/state_shp/'
    )
    # how long a of a rolling window do we want for spei calc (in days)
    roll_window = 60
    # should we save the spei calculations?
    save_spei = True
    # should we load a previous spei calcualtion?
    load_spei = False
    # what is the filename for the load/save location for spei information?
    spei_fname = os.path.join(
        out_dir,
        'spei_df.csv'
    )
    is_drought_fname = os.path.join(
        out_dir,
        'is_drought_df.csv'
    )
    # where is the geodataframe with tile information stored?
    tile_gdf_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/really_chosen_tiles.geojson'
    )
    # what should we save pixel-scale error metric gdf as?
    pixel_error_gdf_fname = os.path.join(
        out_dir,'pixel_error.gdf'
    )
    # what is the name to save the pixel-scale error dict?
    # what analysis sections do we want to run?
    run_pix_analysis = True
    plot_pix_analysis = False
    run_watershed_analysis = True
    plot_watershed_analysis = False
    # must run both pix analysis and watershed analysis for pso analysis to
    # work correctly
    run_pso_analysis = True
    compare_multiple_experiments = False
    analyze_during_drought = False
    plot_spei = False
    # should the all the experiments be plotted?
    skinny_plot = False

    # initiate an instance of the class
    t = get_timeseries(start,end)
    # get the pixels to work with from step 1
    pixels = t.get_pixels(step_1_pixels_fname)
    # print the variables that exist for your info, if desired
    # exits after doing this
    t.print_catch_variables(default_experiment)
    # put the out_dicts for each experiment into a larger dictionary where the
    # key for this dictionary is the name of the experiment
    catch_timeseries = {}
    for e in range(len(exp_names)):
        print('getting timeseries for {}'.format(exp_names[e]))
        this_out_dict = t.get_catch_timeseries(
            exp_dirs[e],exp_names[e],save_bools[e],load_bools[e],
            is_default_experiment[e],read_met_forcing[e],timeseries_dir
        )
        catch_timeseries[exp_names[e]] = this_out_dict
    # get the fluxcom timeseries
    print('getting fluxcom timeseries')
    fluxcom_timeseries = t.get_fluxcom_timeseries(fluxcom_dir)
    # get the streamflow timeseries
    print('getting streamflow timeseries')
    streamflow_timeseries = t.get_strm_timeseries(strm_dir)
    # initiate a class of analyze for looking at daily/pixels scale
    if run_pix_analysis:
        print('running pixel scale analysis')
        a_pix = analyze_pix()
        # save the pft info for later analysis
        pft_info = a_pix.get_pft_info(raw_pft_dir,analyzed_pft_dir,pixels)
        # get the rmse dicts for pixel/daily scale
        print('getting rmse dicts for daily pixel scale')
        default_df_pix,pso_df_pix,fluxcom_df = a_pix.get_rmse_dict(
            exp_names,catch_timeseries,fluxcom_timeseries,experiment_type,
            start_error,end_error
        )
        if plot_pix_analysis:
            # plot LE timeseries and their changes
            #a_pix.plot_le_timeseries(
            #    exp_names,catch_timeseries,fluxcom_timeseries,plots_dir
            #)
            # plot map of rmse and changes for pixel/daily scale
            #print('plotting maps of rmse and average changes at pixel scale')
            a_pix.plot_pso_maps(
                exp_names,default_df_pix,pso_df_pix,plots_dir,skinny_plot
            )
            # plot error and average LE as it progresses through PSO
            # do this for pixel/daily scale
            # do this by PFT
            print('plotting scatter of pso progression by pft')
            a_pix.plot_pso_progression(
                exp_names,default_df_pix,pso_df_pix,fluxcom_df,plots_dir,pft_info
            )
        # let's just plot some general comparisons
        default_df_general = a_pix.get_rmse_dict(
            exp_names[0],catch_timeseries[exp_names[0]],fluxcom_timeseries,
            'general',start_error,end_error
        )
        pso_df_general = a_pix.get_rmse_dict(
            exp_names[2],catch_timeseries[exp_names[2]],fluxcom_timeseries,
            'general',start_error,end_error
        )
        if plot_pix_analysis:
            a_pix.plot_general_comparison_maps(
                [exp_names[0],exp_names[2]],default_df_general,pso_df_general,
                plots_dir,tile_gdf_fname,pixel_error_gdf_fname
            )

    # initiate a class of analyze for watershed/monthly scale
    if run_watershed_analysis:
        print('analyzing outputs at the watershed scale')
        a_water = analyze_watershed()
        # get the intersection info so that we know which pixels correspond to
        # which watersheds
        intersection_info = a_water.get_intersection(intersection_dict_dir)
        # get the rmse dicts for watershed/monthly scale
        print('getting watershed scale model outputs')
        wat_scale_outs = a_water.get_model_preds_watershed(
            start,end,exp_names,catch_timeseries,streamflow_timeseries,
            fluxcom_timeseries,intersection_info
        )
        print('plotting streamflow timeseries')
        #if plot_watershed_analysis:
        #    a_water.plot_streamflow_timeseries(wat_scale_outs,plots_dir,exp_names)
        all_tiles = a_water.get_tiles(step_1_pixels_fname)
        print('getting watershed scale rmse dicts')
        returned_rmse_dfs = a_water.get_rmse_dict(
            wat_scale_outs,lai_map_fname,intersection_info,all_tiles,
            start_error,end_error
        )
        # extract the different rmse dfs that were returned
        # information for the return rmse list:
            # default_rmse_df = returned_rmse[0]
            # pso_init_rmse_df = returned_rmse[1]
            # pso_final_rmse_df = returned_rmse[2]
            # waterwatch_df = returned_rmse[3]
        # plot maps for watershed-scale outputs
        #print('plotting maps for watershed scale outputs')
        if plot_watershed_analysis:
            a_water.plot_water_metrics(
                wat_scale_outs,plots_dir
            )
            #a_water.plot_maps(
            #    returned_rmse_dfs,plots_dir,geojson_fname,exp_names,
            #    states_shp_fname,plot_watershed_analysis,skinny_plot
            #)
    if run_pso_analysis:
        if not run_pix_analysis or not run_watershed_analysis:
            print('you must run both pixel analysis and watershed analysis')
            print('in order to run pso analysis! Change this. Exiting.')
            sys.exit()
        # initiate a class looking at what the pso did
        pso = analyze_pso()
        pso_parameter_names = pso.get_parameter_names()
        pso_all_info = pso.get_all_info(all_info_fname)
        if compare_multiple_experiments:
            comparison_all_info = pso.get_all_info(comparison_all_info_fname)
        pso.plot_parameter_convergence(
            pso_all_info,pso_parameter_names,plots_dir,exp_names
        )
        if not compare_multiple_experiments:
            pso.plot_parameter_maps(
                pso_all_info,pixels,precip_map_fname,lai_map_fname,
                sand_map_fname,k0_map_fname,
                pft_info,default_df_pix,plots_dir,exp_names
            )
        else:
            pso.plot_parameter_maps(
                pso_all_info,pixels,precip_map_fname,lai_map_fname,
                sand_map_fname,k0_map_fname,
                pft_info,default_df_pix,plots_dir,exp_names,
                comparison_all_info
            )
    if compare_multiple_experiments:
        print('comparing to the second experiment')
        # thing to fill in for this compare multiple experiments section
        # the directory for comparison run
        catch_timeseries_compare = {}
        for e in range(len(exp_names)):
            print('getting timeseries for {}'.format(exp_names[e]))
            this_out_dict = t.get_catch_timeseries(
                exp_dirs_compare[e],exp_names_compare[e],save_bools_compare[e],
                load_bools_compare[e],is_default_experiment_compare[e],
                read_met_forcing_compare[e],timeseries_dir
            )
            catch_timeseries_compare[exp_names_compare[e]] = this_out_dict
        a_pix = analyze_pix()
        # save the pft info for later analysis
        pft_info = a_pix.get_pft_info(raw_pft_dir,analyzed_pft_dir,pixels)
        # get the rmse dicts for pixel/daily scale
        print('getting rmse dicts for daily pixel scale')
        default_df_pix_compare,pso_df_pix_compare,fluxcom_df_compare = a_pix.get_rmse_dict(
            exp_names_compare,catch_timeseries_compare,fluxcom_timeseries,
            experiment_type,start_error,end_error
        )
        a_water = analyze_watershed()
        # get the intersection info so that we know which pixels correspond to
        # which watersheds
        intersection_info = a_water.get_intersection(intersection_dict_dir)
        # get the rmse dicts for watershed/monthly scale
        print('getting watershed scale model outputs')
        wat_scale_outs_compare = a_water.get_model_preds_watershed(
            start,end,exp_names_compare,catch_timeseries_compare,
            streamflow_timeseries,fluxcom_timeseries,intersection_info
        )
        all_tiles = a_water.get_tiles(step_1_pixels_fname)
        print('getting watershed scale rmse dicts')
        returned_rmse_dfs_compare = a_water.get_rmse_dict(
            wat_scale_outs_compare,lai_map_fname,intersection_info,all_tiles,
            start_error,end_error
        )
        a_water.plot_maps(
            returned_rmse_dfs_compare,plots_dir,geojson_fname,exp_names_compare,
            states_shp_fname,False,skinny_plot
        )
        comp = compare_experiments()
        comp.plot_diff(
            exp_names,exp_names_compare,pso_df_pix,pso_df_pix_compare,
            returned_rmse_dfs[2],returned_rmse_dfs_compare[2],
            geojson_fname,states_shp_fname,plots_dir
        )
    if analyze_during_drought:
        d = drought()
        spei_df,is_drought_df = d.find_drought(
            fluxcom_timeseries,catch_timeseries[exp_names[0]],roll_window,
            plots_dir,plot_spei,save_spei,load_spei,
            spei_fname,is_drought_fname
        )
        d.trim_to_drought(catch_timeseries[0],is_drought_df)
    print('done with all analysis!')

if __name__ == '__main__':
    main()
