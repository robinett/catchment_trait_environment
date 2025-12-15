import os
from check_water_balance import water_balance
import datetime
import sys

def main():
    # what is the base dir
    base_dir = '/shared/pso/other_analyses/water_balance'
    # where is the data located
    data_dir = os.path.join(
        base_dir,'data'
    )
    # where should we direct outputs?
    out_dir = os.path.join(
        base_dir,'outputs'
    )
    # where should we direct plots?
    plots_dir = os.path.join(
        base_dir,'plots'
    )
    # where should we direct and look for saved timeseries?
    saved_time_dir = (
        '/shared/pso/step_5_analyze_outputs/saved_timeseries'
    )
    # when should we start and ned this analysis?
    start = datetime.date(2006,1,1)
    end = datetime.date(2006,12,31)
    # what pixels did this run have?
    pixels_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/intersecting_catch_tiles.csv'
    )
    # the first experiment: one run with the default everything:
        # run with point forcing
        # run with restart from 1980
        # run on Discover
    default_exp_dir = (
        '/shared/pso/other_analyses/water_balance/data/' +
        'GEOSldas_CN45_res_pf_g1_ksat'
    )
    # what is the name of this experiment for saving and loading/
    default_exp_name = 'res-pf-g1-ksat'
    # what is the directory for the experiment run:
        # with point forcing
        # with restart from 2006
        # on AWS
    pf_exp_dir = (
        '/lustre/catchment/exps/GEOSldas_CN45_med_2006_2006_v3/0'
    )
    # what is the name of this experiment?
    pf_exp_name = 'med-2006-2006-pf-v3'
    # where is the point forcing information located?
    point_forcing_dir = (
        '/shared/point_forcing_huc6/{}_point_forcing_data.nc4'
    )
    # for the following four variables, first element is for deafult and second
    # element is for experiment using point forcing
    # should we load the timeseries?
    load_timeseries = [True,True]
    # should we save the timeseries?
    save_timeseries = [False,False]
    # is the directory that for a default catchment output?
    is_default_experiment = [False,True]
    # we obviously need to read met forcing for this exercise
    read_met_forcing = [True,True]
    # let's put the experiment directories and names in one list
    exp_dirs = [default_exp_dir,pf_exp_dir]
    exp_names = [default_exp_name,pf_exp_name]
    # when comparing to single point driver met forcing files, do these need to
    # be loaded from .nc4? If true, they will be. if false, there should be a
    # .npy file corresponding to these in ./outputs
    single_point_from_nc4 = True
    # where is the intersection fname?
    intersection_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/intersection_info.pkl'
    )
    fluxcom_dir = (
        '/shared/pso/step_3_process_fluxcom/outputs/' +
        'le_truth_fluxcom_rs_meteo_ensemble_watts_per_m2_2001-01-01_' +
        '2009-12-31_selected_tiles.csv'
    )
    # where is the streamflow data located?
    strm_dir = (
        '/shared/pso/step_3.1_process_streamflow/outputs/'+
        'streamflow_data_mm_per_month.csv'
    )
    huc6_shape_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/subset_huc6s.geojson'
    )
    states_shp_fname = (
        '/shared/pso/step_1_choose_tiles/data/state_shp'
    )

    # initiate a class
    ch = water_balance(start,end)
    outs = ch.get_catch_met(
        pixels_fname,exp_dirs,exp_names,save_timeseries,load_timeseries,
        is_default_experiment,read_met_forcing,saved_time_dir,
        intersection_fname,fluxcom_dir,strm_dir
    )
    default_watershed_timeseries = outs[0]
    pf_watershed_timeseries = outs[1]
    waterwatch_timeseries = outs[2]
    default_timeseries_pixel = outs[3]
    pf_timeseries_pixel = outs[4]
    intersection_info = outs[5]
    # get the water balance timeseries
    outs = ch.get_water_balance_timeseries(
        pf_watershed_timeseries,waterwatch_timeseries
    )
    # this is all observation based--has nothing to do with Catch-CN output yet
    # outs in order: merra2_precip, fluxcom_le, waterwatch_strm, incoming,
    # outgoing, diff
    #ch.plot_water_balance(outs,plots_dir,huc6_shape_fname,states_shp_fname)
    ## let's figure out what's going on when we convert these units
    ## first for our default case
    #ch.plot_unit_conversion(
    #    default_timeseries_pixel,plots_dir,intersection_info,
    #    waterwatch_timeseries,exp_names[0]
    #)
    # and next for the pf case
    ch.plot_unit_conversion(
        pf_timeseries_pixel,plots_dir,intersection_info,
        waterwatch_timeseries,exp_names[1]
    )
    ## let's compare what was going on in our met_timeseries_pixel met forcing
    ## with the met forcing from our single point driver output file
    #ch.point_forcing_vs_default(
    #    default_timeseries_pixel,pf_timeseries_pixel,point_forcing_dir,
    #    single_point_from_nc4,out_dir,plots_dir
    #)

if __name__ == '__main__':
    main()
