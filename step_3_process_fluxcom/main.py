import os
from process_fluxcom import process
import datetime

def main():
    ##### INFORMATION FOR USER TO EDIT #####
    # general decisions to be made about inputs
    # highest directory where data and script is stored
    base_dir = '/shared/pso/step_3_process_fluxcom/'
    # directory where the data exists
    data_dir = os.path.join(
        base_dir,'data'
    )
    # directory where to store outputs
    output_dir = os.path.join(
        base_dir,'outputs'
    )
    # directory to store plots made during process
    plots_dir = os.path.join(
        base_dir,'plots'
    )
    # base dir for catchment files
    experiment_dir = (
        '/shared/medComplex/archived_experiments/'+
        'GEOSldas_CN45_default_1980_2021_segFIXED_v2/output/SMAP_EASEv2_M36/'+
        'cat/ens0000/'
    )
    # directory for file that holds catchment pixel information
    # location of the catchment tiles to run
    catch_tiles_fname = os.path.join(
        '/shared/pso/step_1_choose_tiles/outputs/intersecting_catch_tiles.csv'
    )
    # location of the cathment tile info file
    catch_tile_info_fname = os.path.join(
        '/shared/pso/step_1_choose_tiles/data/tile_coord.csv'
    )
    # dates to start and end, inclusive
    start_date = datetime.date(2003,1,1)
    end_date = datetime.date(2014,12,31)
    # example fname for fluxcom. if changed, shoudl be in same format where date can be
    #   changed via fname.format()
    example_fluxcom_fname = os.path.join(
        data_dir,'LE.RS_METEO.EBC-ALL.MLM-ALL.METEO-CRUNCEP_v8.720_360.'+
        'daily.{}.nc'
    )
    # a general fname for a Catchment experiment, if we want to compare fluxcom to catchment
    experiment_fname = os.path.join(
        experiment_dir,'Y{y}/M{m}/GEOSldas_CN45_default_1980_2021_segFIXED_v2.'+
        'tavg24_1d_lnd_Nt.monthly.{y}{m}.nc'
    )
    #name to save spatial plots. Pixel number changed via .format() entry
    spatial_plot_name = os.path.join(
        plots_dir,'spatial_weights_{}.png'
    )
    # name to save pot of fluxcom average for summer day
    fluxcom_plot_name = os.path.join(
        plots_dir,'le_fluxcom_watts_per_m2_doy_180.png'
    )
    # directory and fname where to save calculated weights
    weights_fname = os.path.join(output_dir,'fluxcom_weights_allpix.pkl')
    # where to save truth to. start and end should be editable via .format()
    truth_fname_all_tiles = os.path.join(
        output_dir,'le_truth_fluxcom_rs_CRUNCEP_ensemble_watts_per_m2_'+
        '{start}_{end}_all_tiles.csv'
    )
    truth_fname_selected_tiles = os.path.join(
        output_dir,'le_truth_fluxcom_rs_CRUNCEP_ensemble_watts_per_m2_'+
        '{start}_{end}_camels_tiles.csv'
    )
    # Basic functioning of code:
        # will get weights for all tiles
        # will create truth datasets both for all tiles and for only selected tiles
    # choose what to recreate, load, and be saved
    load_weights = False
    save_weights = False
    load_truth = False
    save_all_truth = True
    save_selected_truth = True

    # choose what sections to actually run
    inspect_fluxcom = False
    get_weights = True
    pixel_plots = False
    create_truth = True
    ########################################

    # see if catchment experiment data is montly
    if 'monthly' in experiment_fname:
        is_monthly = True
    else:
        is_monthly = False

    pro = process(base_dir)
    pro.get_tiles(catch_tiles_fname)
    if inspect_fluxcom:
        pro.inspect_fluxcom(example_fluxcom_fname)
    if get_weights:
        pro.get_weights(
            catch_tile_info_fname,
            load_weights,
            save_weights,
            weights_fname,
            pixel_plots,
            spatial_plot_name,
            fluxcom_plot_name,
            example_fluxcom_fname
        )
    if create_truth:
        pro.create_truth(
            start_date,
            end_date,
            example_fluxcom_fname,
            save_all_truth,
            save_selected_truth,
            truth_fname_all_tiles,
            truth_fname_selected_tiles
        )

if __name__ == '__main__':
    main()
