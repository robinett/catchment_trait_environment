import os
from get_env_covariates import get_covariates
import datetime
import numpy as np

def main():
    # what is the base name of this directory?
    main_dir = (
        '/shared/pso/step_2_env_covariates'
    )
    # directory where data is stored
    data_dir = os.path.join(main_dir,'data')
    # directory where catchment-cn forcing is stored
    catch_dir = os.path.join(data_dir,'catchment_cn_forcing')
    # directory where lai data is stored
    lai_dir = os.path.join(data_dir,'lai')
    # directory to save outputs
    out_dir = os.path.join(main_dir,'outputs')
    # what is the .csv file that lists the tiles at which we will be running?
    tile_fname = (
        '/shared/pso/step_1_choose_tiles/outputs/' +
        'intersecting_catch_tiles.csv'
    )
    # what is the file name that defines the catchment-cn soil information?
    # FYI, the description for the columns of this dataset can be found at:
    # /discover/nobackup/ltakacs/bcs/Icarus-NLv3/Icarus-NLv3_EASE/
    # SMAP_EASEv2_M36/clsm/README
    soil_fname = os.path.join(
        data_dir,'soil_param.dat'
    )
    save_ksat = os.path.join(
        out_dir,'averaged_ksat.nc4'
    )
    save_sand = os.path.join(
        out_dir,'averaged_sand.nc4'
    )
    save_ksat_normal = os.path.join(
        out_dir,'averaged_ksat_normalized.nc4'
    )
    save_sand = os.path.join(
        out_dir,'averaged_sand_normalized.nc4'
    )
    # for canopy height, what is the name of the file?
    canopy_dir = os.path.join(
        data_dir,'canopy_height.csv'
    )
    # and what should we save canopy as?
    save_canopy = os.path.join(
        out_dir,'canopy_height.nc4'
    )
    save_canopy_normal = os.path.join(
        out_dir,'canopy_height_normalized.nc4'
    )
    # what is gldas precip coming in as?
    gldas_precip_dir = os.path.join(
        data_dir,'gldas_precip_kg_m2_s.csv'
    )
    # and what should we save gldas precip as
    gldas_precip_save = os.path.join(
        out_dir,'gldas_avg_precip.nc4'
    )
    gldas_precip_save_normal = os.path.join(
        out_dir,'gldas_avg_precip_normalized.nc4'
    )
    # what is the period over which we want to define average behavior?
    start = np.datetime64('1985-01-01')
    end = np.datetime64('2014-12-31')
    # what is the generic filename for the catchment-cn forcing files?
    generic_forcing = os.path.join(
        catch_dir,
        (
            'GEOSldas_CN45_default_1980_2021_segFIXED_v2.' +
            'tavg24_1d_lnd_Nt.monthly.{year}{month:02d}.nc4'
        )
    )
    # where should we save the averaged precipitation data?
    save_precip = os.path.join(
        out_dir,'averaged_precip.nc4'
    )
    save_precip_norm = os.path.join(
        out_dir,'averaged_precip_norm.nc4'
    )
    # where should we save the averaged lai data?
    save_lai = os.path.join(
        out_dir,'averaged_lai.nc4'
    )

    # create an instance of get_covariates
    cov = get_covariates()
    # get the tiles at which we are going to run
    tiles = cov.get_tiles(tile_fname)
    # get the soil information and save it as nc4
    cov.get_soil_info(
        soil_fname,save_ksat,save_ksat_normal,save_sand,save_ksat_normal
    )
    # get and save the MAP values from Catchment-CN for EF of g1
    #cov.get_env_info(
    #    start,end,generic_forcing,save_precip,save_norm_precip,save_lai
    #)
    # get and save the gldas precip
    cov.get_gldas_info(
        start,end,gldas_precip_dir,gldas_precip_save,gldas_precip_save_normal
    )
    # get and save the canopy height vals
    cov.get_canopy_info(canopy_dir,save_canopy,save_canopy_normal)

if __name__ == '__main__':
    main()
