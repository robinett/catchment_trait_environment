import os
from get_datasets import get
import datetime

def main():
    # where are we?
    base_dir = '/shared/pso/other_analyses/errors_vs_covariates'
    # where is the data?
    data_dir = os.path.join(
        base_dir,'data'
    )
    # where should we store outputs?
    out_dir = os.path.join(
        base_dir,'outputs'
    )
    # where should we make plots?
    plots_dir = os.path.join(
        base_dir,'plots'
    )
    # what are the catchment-CN tiles of interest?
    tiles_fname = (
        '/shared/pso/step_6x_extrapolate/outputs/conus_pix.csv'
    )
    # what is the fname for the file that contains the information on this?
    tile_info_fname = (
        '/shared/pso/step_1_choose_tiles/data/tile_coord.csv'
    )
    # let's start with getting the gldas information
    # general gldas info
    gldas_ex_fname = os.path.join(
        data_dir,
        'gldas',
        'GLDAS_NOAH025_M.A{year}{month:02d}.020.nc4'
    )
    # what are the start and end dates for gldas?
    start_gldas = datetime.date(1985,1,1)
    end_gldas = datetime.date(2014,12,31)
    # do we need to create weights? If not, they should be able to be laoded
    # from the filename provided here
    gldas_create_weights = False
    gldas_weights_fname = os.path.join(
        out_dir,'gldas_weights.csv'
    )
    # plot the weights to make sure the code is working correctly? This only
    # works if create_weights = True
    gldas_plot_weights = False
    gldas_weights_plot_fname = os.path.join(
        plots_dir,'gldas_weights_plot_{}.png'
    )
    # what is the generic way to save the output .csv files?
    gldas_df_fname = os.path.join(
        out_dir,'gldas_{var}.csv'
    )
    gldas_map_fname = os.path.join(
        plots_dir,'gldas_{var}_map.png'
    )
    # information for gldas temperature dataframe
    temp_var_name = 'Tair_f_inst'
    temp_common_name = 'air_temp_K'
    # information for the gldas pet dataframe
    pet_var_name = 'PotEvap_tavg'
    pet_common_name = 'potential_et_W_m2'
    # information for the precip datafram
    precip_var_name = 'Rainf_f_tavg'
    precip_common_name = 'precip_kg_m2_s'
    # information for the gldas latent heat df
    lh_var_name = 'Qle_tavg'
    lh_common_name = 'latent_heat_W_m2'
    # put the gldas info into lists so that it can be passed to the function
    gldas_tech_var_names = [
        temp_var_name,pet_var_name,precip_var_name,
        lh_var_name
    ]
    gldas_common_var_names = [
        temp_common_name,pet_common_name,precip_common_name,
        lh_common_name
    ]
    # where is the canopy height information located?
    canopy_height_fname = os.path.join(
        data_dir,
        'canopy_height',
        'gedi_canopy_height_2019.nc'
    )
    canopy_df_fname = os.path.join(
        out_dir,'canopy_height.csv'
    )
    canopy_map_fname = os.path.join(
        plots_dir,'canopy_height.png'
    )
    canopy_tech_var_name = 'Band1'
    canopy_common_var_name = 'canopy_height_m'
    # do we need to create weights? If not, they should be able to be laoded
    # from the filename provided here
    canopy_create_weights = False
    canopy_weights_fname = os.path.join(
        out_dir,'canopy_weights.csv'
    )
    # plot the weights to make sure the code is working correctly? This only
    # works if create_weights = True
    canopy_plot_weights = False
    canopy_weights_plot_fname = os.path.join(
        plots_dir,'canopy_weights_plot_{}.png'
    )
    # what we need to create
    create_gldas = True
    create_canopy_height = True
    # initiate a class
    g = get(base_dir)
    # get the tiles of interest
    tiles = g.get_tiles(tiles_fname)
    if create_gldas:
        g.average_gldas(
            start_gldas,end_gldas,gldas_ex_fname,gldas_df_fname,
            gldas_tech_var_names,gldas_common_var_names,
            out_dir,plots_dir,tile_info_fname,
            gldas_create_weights,
            gldas_plot_weights,gldas_weights_fname,gldas_weights_plot_fname,
            gldas_map_fname
        )
    if create_canopy_height:
        g.get_canopy_height(
            canopy_height_fname,canopy_df_fname,canopy_tech_var_name,
            canopy_common_var_name,out_dir,plots_dir,tile_info_fname,
            canopy_create_weights,canopy_plot_weights,
            canopy_weights_fname,canopy_weights_plot_fname,
            canopy_map_fname
        )


    ## what is the example fname for catchment-CN data?
    ## this will be used for precip and temp
    #catch_ex_fname = os.path.join(
    #    data_dir,
    #    'catchment_cn_forcing',
    #    'GEOSldas_CN45_default_1980_2021_segFIXED_v2.' +
    #    'tavg24_1d_lnd_Nt.monthly.{year}{month:02d}.nc4'
    #)
    ## what is the example fname for GLDAS data?
    ## this will be used for PET
    #if create_precip:
    #    rg.average_catchcn(
    #        start,end,catch_ex_fname,save_precip_fname,catch_precip_var_name
    #    )
    #else:
    #    precip = rg.load_env(save_precip_fname)

if __name__ == '__main__':
    main()
