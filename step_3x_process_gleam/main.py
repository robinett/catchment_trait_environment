import os
from process_gleam import process
import datetime

def main():
    ##### INFORMATION FOR USER TO EDIT #####
    # general decisions to be made about inputs
    # highest directory where data and script is stored
    base_dir = '/discover/nobackup/trobinet/from_aws/pso/step_3x_process_gleam/'
    # directory where the data exists
    data_dir = '/discover/nobackup/trobinet/from_aws/pso/step_3_process_gleam/data'
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
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/' +
        'outputs/extrapolation_pix.csv'
    )
    # location of the cathment tile info file
    catch_tile_info_fname = os.path.join(
        '/discover/nobackup/trobinet/from_aws/pso/' +
        'step_1_choose_tiles/data/tile_coord.csv'
    )
    # dates to start and end, inclusive
    start_date = datetime.date(1995,1,1)
    end_date = datetime.date(2014,12,31)
    # start and end dates for the gleam files that we have
    gleam_start = datetime.date(1995,1,1)
    gleam_end = datetime.date(2014,12,31)
    # example fname for fluxcom. if changed, shoudl be in same format where date can be
    #   changed via fname.format()
    example_fluxcom_fname = os.path.join(
        data_dir,
        'gleam_38a',
        'E_{}_GLEAM_v3.8a.nc'
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
        plots_dir,'le_gleam_watts_per_m2_doy_180.png'
    )
    # directory and fname where to save calculated weights
    weights_fname = os.path.join(output_dir,'fluxcom_weights_allpix.pkl')
    # where to save truth to. start and end should be editable via .format()
    truth_fname_all_tiles = os.path.join(
        output_dir,'le_truth_gleam_38a_watts_per_m2_'+
        '{start}_{end}_all_tiles.csv'
    )
    truth_fname_selected_tiles = os.path.join(
        output_dir,'le_truth_gleam_38a_watts_per_m2_'+
        '{start}_{end}_extrapolation_tiles.csv'
    )
    # Basic functioning of code:
        # will get weights for all tiles
        # will create truth datasets both for all tiles and for only selected tiles
    # choose what to recreate, load, and be saved
    load_weights = True
    save_weights = False
    load_truth = False
    save_all_truth = False
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
            truth_fname_selected_tiles,
            gleam_start,
            gleam_end
        )

if __name__ == '__main__':
    main()
