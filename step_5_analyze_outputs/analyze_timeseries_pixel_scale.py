import sys
sys.path.insert(0,'/shared/pso/step_1_choose_tiles')
import pandas as pd
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime
import geopandas as gpd
from choose_tiles_camels import choose

class analyze_pix:
    def __init__(self):
        pass
    def get_pft_info(self,original_fname,save_fname,pixels):
        # get the pft data as displayed in Catchment-CN
        orig_data = np.genfromtxt(original_fname)
        # link the pft codes in this dataset to the names of these pft codes
        pft_codes = {
            1  : 'Needleleaf evergreen temperate tree',
            2  : 'Needleleaf evergreen boreal tree',
            3  : 'Needleleaf deciduous boreal tree',
            4  : 'Broadleaf evergreen tropical tree',
            5  : 'Broadleaf evergreen temperate tree',
            6  : 'Broadleaf deciduous tropical tree',
            7  : 'Broadleaf deciduous temperate tree',
            8  : 'Broadleaf deciduous boreal tree',
            9  : 'Broadleaf evergreen temperate shrub',
            10 : 'Broadleaf deciduous temperate shrub',
            11 : 'Broadleaf deciduous temperate shrub[moisture stress only]',
            12 : 'Broadleaf deciduous boreal shrub',
            13 : 'Arctic c3 grass',
            14 : 'Cool c3 grass',
            15 : 'Cool c3 grass [moisture stress only]',
            16 : 'Warm c4 grass',
            17 : 'Warm c4 grass [moisture stress only]',
            18 : 'Crop',
            19 : 'Crop [moisture stress only]',
            20 : '(Corn)',
            21 : '(Irrigated corn)',
            22 : '(Spring temperate cereal)',
            23 : '(Irrigated spring temperate cereal)',
            24 : '(winter temperate cereal)',
            25 : '(Irrigated winter temperate cereal)',
            26 : '(Soybean)',
            27 : '(Irrigated Soybean)',
        }
        # link the pft codes in the original dataset to the simple
        # classification names that we use
        pft_codes_simple = {
            1  : 'Forest',
            2  : 'Forest',
            3  : 'Forest',
            4  : 'Forest',
            5  : 'Forest',
            6  : 'Forest',
            7  : 'Forest',
            8  : 'Forest',
            9  : 'Shrublands',
            10 : 'Shrublands',
            11 : 'Shrublands',
            12 : 'Shrublands',
            13 : 'Grasslands',
            14 : 'Grasslands',
            15 : 'Grasslands',
            16 : 'Grasslands',
            17 : 'Grasslands',
            18 : 'Croplands',
            19 : 'Croplands',
            20 : 'Croplands',
            21 : 'Croplands',
            22 : 'Croplands',
            23 : 'Croplands',
            24 : 'Croplands',
            25 : 'Croplands',
            26 : 'Croplands',
            27 : 'Croplands'
        }
        # extract the data of interest from this original file
        # first get the index of the pixels that we are intersted in
        pixels_idx = pixels - 1
        tiles = orig_data[pixels_idx,0]
        # make sure this matches the pixels that we shoud be looking for
        if (np.all(tiles==pixels)) != True:
            print('There is an error in getting the correct pixels from' +
                  'the pft information. Please check error and continue.')
            sys.exit()
        # get the code corresponding to the pfts at each pixel
        veg_1_code = orig_data[pixels_idx,2]
        veg_2_code = orig_data[pixels_idx,3]
        veg_3_code = orig_data[pixels_idx,4]
        veg_4_code = orig_data[pixels_idx,5]
        # get the precent of each pft at this pixel
        veg_1_perc = orig_data[pixels_idx,6]
        veg_2_perc = orig_data[pixels_idx,7]
        veg_3_perc = orig_data[pixels_idx,8]
        veg_4_perc = orig_data[pixels_idx,9]
        # get the specific pft name at each pixel
        veg_1_name = [pft_codes[pft] for pft in veg_1_code]
        veg_2_name = [pft_codes[pft] for pft in veg_2_code]
        veg_3_name = [pft_codes[pft] for pft in veg_3_code]
        veg_4_name = [pft_codes[pft] for pft in veg_4_code]
        # get the simple pft name at each pixel
        veg_1_simple = [pft_codes_simple[pft] for pft in veg_1_code]
        veg_2_simple = [pft_codes_simple[pft] for pft in veg_2_code]
        veg_3_simple = [pft_codes_simple[pft] for pft in veg_3_code]
        veg_4_simple = [pft_codes_simple[pft] for pft in veg_4_code]
        # put all of this information into a dataframe
        pft_df = pd.DataFrame({
            'tile':tiles,
            'pft_1_perc':veg_1_perc,
            'pft_2_perc':veg_2_perc,
            'pft_3_perc':veg_3_perc,
            'pft_4_perc':veg_4_perc,
            'pft_1_name':veg_1_name,
            'pft_2_name':veg_2_name,
            'pft_3_name':veg_3_name,
            'pft_4_name':veg_4_name,
            'pft_1_simple':veg_1_simple,
            'pft_2_simple':veg_2_simple,
            'pft_3_simple':veg_3_simple,
            'pft_4_simple':veg_4_simple
        })
        # get the simple name of the primary pft for each pixels
        all_max_perc_pft = np.zeros(len(pixels),dtype='U25')
        for p,pix in enumerate(pixels):
            all_perc = [
                veg_1_perc[p],
                veg_2_perc[p],
                veg_3_perc[p],
                veg_4_perc[p]
            ]
            max_idx = np.where(
                all_perc==np.amax(all_perc)
            )
            max_idx = max_idx[0][0]
            max_pft_num = max_idx + 1
            max_perc_pft = pft_df['pft_{}_simple'.format(max_pft_num)].iloc[p]
            all_max_perc_pft[p] = max_perc_pft
        # add the name of each primary pft to the df
        pft_df['primary_pft'] = all_max_perc_pft
        # and finally let's make the tiles the index
        pft_df = pft_df.set_index('tile')
        pft_df.to_csv(save_fname)
        return pft_df
    def get_rmse_dict(self,exp_names,exps_dict,fluxcom_timeseries,
                      experiment_type,start_error,end_error,out_dir,
                      save_default_error):
        '''
        Function to generate a dictionary containing all of the RMSE
        information that we could want. For this dictionary, each key will be
        an experiment name that points to a dataframe. in this df, each column
        will be a different pixel (plus one column that is 'all') and the rows
        are [average le, rmse, % change rmse, rmse after first iteration, %
        change rmse after first iteration]
        Options for experiment_type are: 'pso','general'
        '''
        # get the pixels for this analysis
        pixels = list(fluxcom_timeseries.columns)
        # get rid of time from this list of pixels
        pixels = pixels[1:]
        this_cols = list(pixels)
        # add a column that is the metric for all columns
        this_cols.append('all')
        # set up the DataFrame
        default_df = pd.DataFrame(columns=this_cols)
        # get the information that we want for the default experiment
        # first extract the information for the default experiment
        if experiment_type == 'general':
            default_dict = exps_dict
        elif experiment_type == 'pso':
            default_dict = exps_dict[exp_names[0]]
        else:
            print(
                'ERROR! Must specify experiment type as either PSO'+
                'general'
            )
            sys.exit()
        # extract the variables that we care about
        default_le = default_dict['le']
        default_ave_sm = default_dict['ave_sm']
        default_root_sm = default_dict['root_sm']
        default_surf_sm = default_dict['surf_sm']
        default_infil = default_dict['infil']
        default_runoff = default_dict['runoff']
        default_baseflow = default_dict['baseflow']
        times = list(default_le.index)
        times = pd.to_datetime(times)
        # let's see where is our start and end idx
        for t,ti in enumerate(times):
            if ti.date() == start_error:
                start_idx = t
            if ti.date() == end_error:
                end_idx = t
        default_le = default_le.iloc[start_idx:end_idx+1]
        default_ave_sm = default_ave_sm.iloc[start_idx:end_idx+1]
        default_root_sm = default_root_sm.iloc[start_idx:end_idx+1]
        default_surf_sm = default_surf_sm.iloc[start_idx:end_idx+1]
        default_infil = default_infil.iloc[start_idx:end_idx+1]
        default_runoff = default_runoff.iloc[start_idx:end_idx+1]
        default_baseflow = default_baseflow.iloc[start_idx:end_idx+1]
        # and take the averages of these variables that we care about
        ave_le = list(default_le.mean())
        ave_ave_sm = list(default_ave_sm.mean())
        ave_root_sm = list(default_root_sm.mean())
        ave_surf_sm = list(default_surf_sm.mean())
        ave_infil = list(default_infil.mean())
        ave_runoff = list(default_runoff.mean())
        ave_baseflow = list(default_baseflow.mean())
        # add the average of all for the 'all' column
        ave_le.append(np.mean(ave_le))
        default_df.loc['ave_le'] = ave_le
        ave_ave_sm.append(np.mean(ave_ave_sm))
        default_df.loc['ave_ave_sm'] = ave_ave_sm
        ave_root_sm.append(np.mean(ave_root_sm))
        default_df.loc['ave_root_sm'] = ave_root_sm
        ave_surf_sm.append(np.mean(ave_surf_sm))
        default_df.loc['ave_surf_sm'] = ave_surf_sm
        ave_infil.append(np.mean(ave_infil))
        default_df.loc['ave_infil'] = ave_infil
        ave_runoff.append(np.mean(ave_runoff))
        default_df.loc['ave_runoff'] = ave_runoff
        ave_baseflow.append(np.mean(ave_baseflow))
        default_df.loc['ave_baseflow'] = ave_baseflow
        # now let's calculate rmse as compared to fluxcom
        fluxcom_timeseries = fluxcom_timeseries.set_index('time')
        fluxcom_times_str = list(fluxcom_timeseries.index)
        fluxcom_times = [
            datetime.datetime.strptime(
                d,'%Y-%m-%d'
            ) for d in fluxcom_times_str
        ]
        # let's see where is our start and end idx
        for t,ti in enumerate(fluxcom_times):
            if ti.date() == start_error:
                fluxcom_start_idx = t
            if ti.date() == end_error:
                fluxcom_end_idx = t
        fluxcom_timeseries = fluxcom_timeseries.iloc[
            fluxcom_start_idx:fluxcom_end_idx+1
        ]
        # to hold all rmse
        default_all_le_rmse = np.zeros(len(this_cols))
        default_all_le_r2 = np.zeros(len(this_cols))
        default_all_le_corr = np.zeros(len(this_cols))
        default_all_le_ubrmse = np.zeros(len(this_cols))
        default_all_le_ubrmse_norm = np.zeros(len(this_cols))
        for p,pix in enumerate(pixels):
            # get the default and fluxcom le
            this_default_le = np.array(default_le[pix])
            this_fluxcom = np.array(fluxcom_timeseries[pix])
            this_avg_obs = np.nanmean(this_fluxcom)
            # get the locations where we want to calculate
            this_calc_idx = np.where(
                (this_default_le > 0) & (np.isnan(this_fluxcom) == False)
            )
            # if these criteria can't be met, assign nan. otherwise calcluate
            # the rmse
            if len(this_calc_idx[0]) > 0:
                # calculate rmse
                this_rmse = np.sqrt(
                    (
                        (
                            this_default_le[this_calc_idx] - 
                            this_fluxcom[this_calc_idx]
                        )**2
                    ).mean()
                )
                # calculate r2
                this_avg_y = np.mean(this_fluxcom[this_calc_idx])
                this_r2 = 1 - (
                    (
                        np.sum(
                            np.square(
                                this_fluxcom[this_calc_idx] -
                                this_default_le[this_calc_idx]
                            )
                        )
                    )/(
                        np.sum(
                            np.square(
                                this_fluxcom[this_calc_idx] -
                                this_avg_y
                            )
                        )
                    )
                )
                # calculate correlation coefficient
                this_avg_x = np.mean(this_default_le[this_calc_idx])
                numerator = np.sum(
                    (this_default_le[this_calc_idx] - this_avg_x)*
                    (this_fluxcom[this_calc_idx] - this_avg_y)
                )
                denominator_1 = np.sqrt(
                    np.sum(
                        np.square(
                            this_default_le[this_calc_idx] - this_avg_x
                        )
                    )
                )
                denominator_2 = np.sqrt(
                    np.sum(
                        np.square(
                            this_fluxcom[this_calc_idx] - this_avg_y
                        )
                    )
                )
                denominator = denominator_1*denominator_2
                this_le_corr = numerator/denominator
                # calculate ubrmse
                this_ubrmse = np.sqrt(
                    (
                        (
                            (
                                this_default_le[this_calc_idx] -
                                this_avg_x
                            ) - (
                                this_fluxcom[this_calc_idx] -
                                this_avg_y
                            )
                        )**2
                    ).mean()
                )
                # get pixel-normalized ubrmse
                this_ubrmse_norm = this_ubrmse/this_avg_obs
            else:
                this_rmse = np.nan
                this_r2 = np.nan
                this_le_corr = np.nan
                this_ubrmse = np.nan
            default_all_le_rmse[p] = this_rmse
            default_all_le_r2[p] = this_r2
            default_all_le_corr[p] = this_le_corr
            default_all_le_ubrmse[p] = this_ubrmse
            default_all_le_ubrmse_norm[p] = this_ubrmse_norm
        # now let's get the rmse for all pixels for the 'all' column
        fluxcom_np = np.array(fluxcom_timeseries)
        default_le_np = np.array(default_le)
        calc_idx = np.where(
            (np.isnan(fluxcom_np) == False) & (default_le_np > 0)
        )
        # calculate the rmse across all pixels and all times
        rmse_all = np.sqrt(
            (
                (
                    default_le_np[calc_idx] -
                    fluxcom_np[calc_idx]
                )**2
            ).mean()
        )
        # calcualte the r2 across all pixels and all times
        all_avg_y = np.average(fluxcom_np[calc_idx])
        le_r2_all = 1 - (
            (
                np.sum(
                    np.square(
                        fluxcom_np[calc_idx] -
                        default_le_np[calc_idx]
                    )
                )
            )/(
                np.sum(
                    np.square(
                        fluxcom_np[calc_idx] -
                        all_avg_y
                    )
                )
            )
        )
        # calculate the correlation accross all pixels and all times
        all_avg_x = np.mean(default_le_np[calc_idx])
        numerator = np.sum(
            (default_le_np[calc_idx] - all_avg_x)*
            (fluxcom_np[calc_idx] - all_avg_y)
        )
        denominator_1 = np.sqrt(
            np.sum(
                np.square(
                    default_le_np[calc_idx] - all_avg_x
                )
            )
        )
        denominator_2 = np.sqrt(
            np.sum(
                np.square(
                    fluxcom_np[calc_idx] - all_avg_y
                )
            )
        )
        denominator = denominator_1*denominator_2
        le_corr_all = numerator/denominator
        # calculate the ubrmse across all pixels and all times
        ubrmse_all = np.sqrt(
            (
                (
                    (
                        default_le_np[calc_idx] -
                        all_avg_x
                    ) - (
                        fluxcom_np[calc_idx] -
                        all_avg_y
                    )
                )**2
            ).mean()
        )
        # calculate the normalized ubrmse everywhere
        ubrmse_norm_all = np.nanmean(default_all_le_ubrmse_norm[:-1])
        # add these to their arrays
        default_all_le_rmse[-1] = rmse_all
        default_all_le_r2[-1] = le_r2_all
        default_all_le_corr[-1] = le_corr_all
        default_all_le_ubrmse[-1] = ubrmse_all
        default_all_le_ubrmse_norm[-1] = ubrmse_norm_all
        # add this to the default dataframe
        default_df.loc['le_rmse'] = default_all_le_rmse
        default_df.loc['le_r2'] = default_all_le_r2
        default_df.loc['le_corr'] = default_all_le_corr
        default_df.loc['le_ubrmse'] = default_all_le_ubrmse
        default_df.loc['le_ubrmse_norm'] = default_all_le_ubrmse_norm
        if save_default_error != 'none':
            error_no_all = default_all_le_ubrmse[:-1]
            error_no_all.tofile(
                os.path.join(
                    out_dir,
                    save_default_error
                ),
                sep = ','
            )
        # for this, there is no percent change in rmse, because this is the
        # baseline off which percent change will be caluclated. will fill with
        # nans
        default_perc_change = np.repeat(np.nan,len(this_cols))
        default_df.loc['perc_change_le_rmse'] = default_perc_change
        # also no iterations on this because it is default, so no rmse after
        # first iteration and no percent change after the first iteration
        default_rmse_first_it = np.repeat(np.nan,len(this_cols))
        default_first_it_change = np.repeat(np.nan,len(this_cols))
        default_df.loc['le_rmse_first_iteration'] = default_rmse_first_it
        default_df.loc['perc_change_le_rmse_first_it'] = default_first_it_change
        default_df.loc['le_r2_first_iteration'] = default_rmse_first_it
        default_df.loc['perc_change_le_r2_first_it'] = default_first_it_change
        # now let's make one for the fluxcom observations, but just include the
        # average LE
        fluxcom_le_means = list(fluxcom_timeseries.mean())
        fluxcom_le_means.append(np.nanmean(fluxcom_le_means))
        fluxcom_df = pd.DataFrame(columns=this_cols)
        fluxcom_df.loc['ave_le'] = fluxcom_le_means
        default_df.loc['le_obs'] = fluxcom_le_means
        # add lat and lon for plotting purposes later
        # get lats and lons
        all_lats = default_dict['lat']
        all_lons = default_dict['lon']
        lats = list(all_lats.iloc[0])
        lons = list(all_lons.iloc[0])
        lats.append(np.nan)
        lons.append(np.nan)
        # add lats and lons to dfs
        default_df.loc['lon'] = lons
        default_df.loc['lat'] = lats
        fluxcom_df.loc['lon'] = lons
        fluxcom_df.loc['lats'] = lats
        # let's get average le difference
        default_df.loc['ave_le_diff'] = (
            default_df.loc['ave_le'] - fluxcom_df.loc['ave_le']
        )
        # if we're just going for general experiment type return that here and
        # move on
        if experiment_type == 'general':
            return default_df

        # now let's do this for the comparison experiment
        pso_df = pd.DataFrame(columns=this_cols)
        # get the information that we want for the pso experiment
        # first extract the information for the pso experiment
        pso_dict = exps_dict[exp_names[2]]
        pso_le = pso_dict['le']
        pso_le = pso_le.iloc[start_idx:end_idx+1]
        # let's get the average le
        ave_le = list(pso_le.mean())
        ave_le.append(np.mean(ave_le))
        pso_df.loc['ave_le'] = ave_le
        pso_df.loc['le_obs'] = fluxcom_le_means
        # to hold all rmse
        pso_all_le_rmse = np.zeros(len(this_cols))
        pso_all_le_r2 = np.zeros(len(this_cols))
        pso_all_le_corr = np.zeros(len(this_cols))
        pso_all_le_ubrmse = np.zeros(len(this_cols))
        pso_all_le_ubrmse_norm = np.zeros(len(this_cols))
        for p,pix in enumerate(pixels):
            this_pso_le = np.array(pso_le[pix])
            this_fluxcom = np.array(fluxcom_timeseries[pix])
            this_avg_obs = np.nanmean(this_fluxcom)
            # get where we want to calculate
            this_calc_idx = np.where(
                (this_pso_le > 0) & (np.isnan(this_fluxcom) == False)
            )
            # if these criteria can't be met, assign nan. otherwise calcluate
            # the rmse
            if len(this_calc_idx[0]) > 0:
                this_rmse = np.sqrt(
                    (
                        (
                            this_pso_le[this_calc_idx] - 
                            this_fluxcom[this_calc_idx]
                        )**2
                    ).mean()
                )
                this_avg_y = np.mean(this_fluxcom[this_calc_idx])
                this_r2 = 1 - (
                    (
                        np.sum(
                            np.square(
                                this_fluxcom[this_calc_idx] -
                                this_pso_le[this_calc_idx]
                            )
                        )
                    )/(
                        np.sum(
                            np.square(
                                this_fluxcom[this_calc_idx] -
                                this_avg_y
                            )
                        )
                    )
                )
                # calculate correlation coefficient
                this_avg_x = np.mean(this_pso_le[this_calc_idx])
                numerator = np.sum(
                    (this_pso_le[this_calc_idx] - this_avg_x)*
                    (this_fluxcom[this_calc_idx] - this_avg_y)
                )
                denominator_1 = np.sqrt(
                    np.sum(
                        np.square(
                            this_pso_le[this_calc_idx] - this_avg_x
                        )
                    )
                )
                denominator_2 = np.sqrt(
                    np.sum(
                        np.square(
                            this_fluxcom[this_calc_idx] - this_avg_y
                        )
                    )
                )
                denominator = denominator_1*denominator_2
                this_le_corr = numerator/denominator
                # calculate unbiased rmse
                this_ubrmse = np.sqrt(
                    (
                        (
                            (
                                this_pso_le[this_calc_idx] -
                                this_avg_x
                            ) - (
                                this_fluxcom[this_calc_idx] -
                                this_avg_y
                            )
                        )**2
                    ).mean()
                )
                # calculate the normalized unbiased rmse
                this_ubrmse_norm = this_ubrmse/this_avg_obs
            else:
                this_rmse = np.nan
                this_r2 = np.nan
                this_le_corr = np.nan
                this_ubrmse = np.nan
            pso_all_le_rmse[p] = this_rmse
            pso_all_le_r2[p] = this_r2
            pso_all_le_corr[p] = this_le_corr
            pso_all_le_ubrmse[p] = this_ubrmse
            pso_all_le_ubrmse_norm[p] = this_ubrmse_norm
        # now let's get the rmse for all pixels for the 'all' column
        fluxcom_np = np.array(fluxcom_timeseries)
        pso_le_np = np.array(pso_le)
        calc_idx = np.where(
            (np.isnan(fluxcom_np) == False) & (pso_le_np > 0)
        )
        rmse_all = np.sqrt(
            (
                (
                    pso_le_np[calc_idx] - 
                    fluxcom_np[calc_idx]
                )**2
            ).mean()
        )
        all_avg_y = np.average(fluxcom_np[calc_idx])
        le_r2_all = 1 - (
            (
                np.sum(
                    np.square(
                        fluxcom_np[calc_idx] -
                        pso_le_np[calc_idx]
                    )
                )
            )/(
                np.sum(
                    np.square(
                        fluxcom_np[calc_idx] -
                        all_avg_y
                    )
                )
            )
        )
        # calculate the correlation accross all pixels and all times
        all_avg_x = np.mean(pso_le_np[calc_idx])
        numerator = np.sum(
            (pso_le_np[calc_idx] - all_avg_x)*
            (fluxcom_np[calc_idx] - all_avg_y)
        )
        denominator_1 = np.sqrt(
            np.sum(
                np.square(
                    pso_le_np[calc_idx] - all_avg_x
                )
            )
        )
        denominator_2 = np.sqrt(
            np.sum(
                np.square(
                    fluxcom_np[calc_idx] - all_avg_y
                )
            )
        )
        denominator = denominator_1*denominator_2
        le_corr_all = numerator/denominator
        # calculate the ubrmse across all pixels and all times
        ubrmse_all = np.sqrt(
            (
                (
                    (
                        pso_le_np[calc_idx] -
                        all_avg_x
                    ) - (
                        fluxcom_np[calc_idx] -
                        all_avg_y
                    )
                )**2
            ).mean()
        )
        # calculate the normalized ubrmse as the average
        ubrmse_norm_all = np.mean(pso_all_le_ubrmse_norm[:-1])
        # add these to their respective arrays
        pso_all_le_rmse[-1] = rmse_all
        pso_all_le_r2[-1] = le_r2_all
        pso_all_le_corr[-1] = le_corr_all
        pso_all_le_ubrmse[-1] = ubrmse_all
        pso_all_le_ubrmse_norm[-1] = ubrmse_norm_all
        # add this to the pso dataframe
        pso_df.loc['le_rmse'] = pso_all_le_rmse
        pso_df.loc['le_r2'] = pso_all_le_r2
        pso_df.loc['le_corr'] = pso_all_le_corr
        pso_df.loc['le_ubrmse'] = pso_all_le_ubrmse
        pso_df.loc['le_ubrmse_norm'] = pso_all_le_ubrmse_norm
        # calculate what percent change in rmse this is compared to the first
        # iteration
        change_le_rmse = pso_all_le_rmse - default_all_le_rmse
        perc_change_le_rmse = (
            (pso_all_le_rmse - default_all_le_rmse)/default_all_le_rmse
        )
        pso_df.loc['change_le_rmse'] = change_le_rmse
        pso_df.loc['perc_change_le_rmse'] = perc_change_le_rmse
        change_le_r2 = pso_all_le_r2 - default_all_le_r2
        change_le_corr = pso_all_le_corr - default_all_le_corr
        change_le_ubrmse = pso_all_le_ubrmse - default_all_le_ubrmse
        perc_change_le_ubrmse = (
            change_le_ubrmse/default_all_le_ubrmse
        )
        change_le_ubrsme_norm = (
            pso_all_le_ubrmse_norm - default_all_le_ubrmse_norm
        )
        pso_df.loc['change_le_r2'] = change_le_r2
        pso_df.loc['change_le_corr'] = change_le_corr
        pso_df.loc['change_le_ubrmse'] = change_le_ubrmse
        pso_df.loc['perc_change_le_ubrmse'] = perc_change_le_ubrmse
        pso_df.loc['change_le_ubrmse_norm'] = (
            change_le_ubrsme_norm
        )
        # now caluclate the rmse after the first iteration
        # get information
        first_it_dict = exps_dict[exp_names[1]]
        first_it_le = first_it_dict['le']
        first_it_le = first_it_le.iloc[start_idx:end_idx+1]
        # let's get the average le for this first iteration
        ave_le = list(first_it_le.mean())
        ave_le.append(np.mean(ave_le))
        pso_df.loc['ave_le_first_it'] = ave_le
        # calculate rmse
        first_it_all_le_rmse = np.zeros(len(this_cols))
        for p,pix in enumerate(pixels):
            this_first_it_le = np.array(first_it_le[pix])
            this_fluxcom = np.array(fluxcom_timeseries[pix])
            # get the locations where we want to calculate
            this_calc_idx = np.where(
                (this_first_it_le > 0) & (np.isnan(this_fluxcom) == False)
            )
            # if these criteria can't be met, assign nan. otherwise calcluate
            # the rmse
            if len(this_calc_idx[0]) > 0:
                this_rmse = np.sqrt(
                    (
                        (
                            this_first_it_le[this_calc_idx] - 
                            this_fluxcom[this_calc_idx]
                        )**2
                    ).mean()
                )
            else:
                this_rmse = np.nan
            first_it_all_le_rmse[p] = this_rmse
        # now let's get the rmse for all pixels for the 'all' column
        fluxcom_np = np.array(fluxcom_timeseries)
        first_it_le_np = np.array(first_it_le)
        calc_idx = np.where(
            (np.isnan(fluxcom_np) == False) & (first_it_le_np > 0)
        )
        rmse_all = np.sqrt(
            (
                (
                    first_it_le_np[calc_idx] - 
                    fluxcom_np[calc_idx]
                )**2
            ).mean()
        )
        first_it_all_le_rmse[-1] = rmse_all
        pso_df.loc['le_rmse_first_iteration'] = first_it_all_le_rmse
        # get the percent change after the first iteration
        perc_change_le_rmse_it_1 = (
            (first_it_all_le_rmse - default_all_le_rmse)/default_all_le_rmse
        )
        pso_df.loc['perc_change_le_rmse_first_it'] = perc_change_le_rmse_it_1

        # add lon lats to these new dfs
        pso_df.loc['lon'] = lons
        pso_df.loc['lat'] = lats
        # add the average le difference to these new dataframes
        pso_df.loc['ave_le_diff'] = (
            pso_df.loc['ave_le'] - fluxcom_df.loc['ave_le']
        )
        pso_df.loc['ave_le_diff_first_it'] = (
            pso_df.loc['ave_le_first_it'] - fluxcom_df.loc['ave_le']
        )
        return [default_df,pso_df,fluxcom_df]
    def plot_le_timeseries(self,exp_names,catch_timeseries,fluxcom_timeseries,
                           plots_dir):
        #print(catch_timeseries[exp_names[0]])
        #sys.exit()
        default_timeseries_le = catch_timeseries[exp_names[0]]['le']
        first_it_timeseries_le = catch_timeseries[exp_names[1]]['le']
        final_it_timeseries_le = catch_timeseries[exp_names[2]]['le']
        #default_timeseries_temp = catch_timeseries[exp_names[0]]['tair']
        fluxcom_timeseries = fluxcom_timeseries.set_index('time')
        #tdefault_timeseries_q = catch_timeseries[exp_names[0]]['qair']
        times = list(default_timeseries_le.index)
        pixels = list(default_timeseries_le.columns)
        times_str = [str(t) for t in times]
        dates_dtm = [
            datetime.datetime.strptime(d[:10],'%Y-%m-%d') for d in times_str
        ]
        this_exp = exp_names[2]
        day_diff = (
            len(default_timeseries_le[pixels[0]]) -
            len(fluxcom_timeseries[pixels[0]])
        )
        if day_diff != 0:
            fluxcom_extended = np.zeros((
                len(default_timeseries_le[pixels[0]]),
                len(default_timeseries_le.columns)
            ))
            for p,pix in enumerate(pixels):
                fluxcom_extended[:,p] = np.append(
                    np.zeros(day_diff) + np.nan,
                    fluxcom_timeseries[pix]
                )
        for p,pix in enumerate(pixels):
            print('plotting et timeseries for pixel {}'.format(pix))
            this_default = default_timeseries_le[pix]
            this_first_it = first_it_timeseries_le[pix]
            this_final_it = final_it_timeseries_le[pix]
            if day_diff == 0:
                this_fluxcom = fluxcom_timeseries[pix]
            else:
                this_fluxcom = fluxcom_extended[:,p]
            this_times = np.array(fluxcom_timeseries.index)
            this_savename = os.path.join(
                plots_dir,'le_timeseries_{}_{}.png'.format(
                    this_exp,pix
                )
            )
            plt.figure()
            plt.plot(
                dates_dtm,this_default,label='CN4.5 EF It. 1',
                c='r',linewidth=.8
            )
            plt.plot(
                dates_dtm,this_fluxcom,label='GLEAM',
                c='k',linewidth=.8
            )
            plt.plot(
                dates_dtm,this_first_it,label='CN4.5 EF It. 4',
                c='y',linewidth=.8
            )
            plt.plot(
                dates_dtm,this_final_it,label='CN4.5 EF It. 8',
                c='g',linewidth=.8
            )
            plt.legend()
            plt.xlabel('Date')
            plt.ylabel('LE (W/m2)')
            plt.savefig(this_savename,dpi=300,bbox_inches='tight')
            plt.close()
            start_2010 = datetime.date(2010,1,1)
            start_2010_str = datetime.datetime.strftime(start_2010,'%Y-%m-%d')
            end_2012 = datetime.date(2012,12,31)
            end_2012_str = datetime.datetime.strftime(end_2012,'%Y-%m-%d')
            start_idx = np.where(this_times == start_2010_str)[0][0]
            end_idx = np.where(this_times == end_2012_str)[0][0]
            start_idx_catch = start_idx - day_diff
            end_idx_catch = end_idx - day_diff
            this_default_2010_2012 = default_timeseries_le[pix].iloc[
                start_idx_catch:end_idx_catch+1
            ]
            this_first_it_2010_2012 = first_it_timeseries_le[pix].iloc[
                start_idx_catch:end_idx_catch+1
            ]
            this_final_it_2010_2012 = final_it_timeseries_le[pix].iloc[
                start_idx_catch:end_idx_catch+1
            ]
            this_fluxcom_2010_2012 = fluxcom_timeseries[pix].iloc[
                start_idx:end_idx+1
            ]
            dates_dtm_2010_2012 = dates_dtm[start_idx:end_idx+1]
            this_savename = os.path.join(
                plots_dir,'le_timeseries_2010_2012_{}_{}'.format(
                    this_exp,pix
                )
            )
            plt.figure()
            plt.plot(
                dates_dtm_2010_2012,this_default_2010_2012,
                label='CN4.5 EF It. 1',
                c='r',linewidth=.5
            )
            plt.plot(
                dates_dtm_2010_2012,this_fluxcom_2010_2012,label='GLEAM',
                c='k',linewidth=.5
            )
            plt.plot(
                dates_dtm_2010_2012,this_first_it_2010_2012,
                label='CN4.5 EF It. 4',
                c='y',linewidth=.5
            )
            plt.plot(
                dates_dtm_2010_2012,this_final_it_2010_2012,
                label='CN4.5 EF It. 8',
                c='g',linewidth=.5
            )
            plt.legend()
            plt.xlabel('Date')
            plt.ylabel('LE (W/m2)')
            plt.savefig(this_savename,dpi=300,bbox_inches='tight')
            plt.close()
        #for p,pix in enumerate(pixels):
        #    print('plotting tair timeseries for pixel {}'.format(pix))
        #    this_default = default_timeseries_temp[pix]
        #    this_savename = os.path.join(
        #        plots_dir,'tair_timeseries_{}_{}.png'.format(this_exp,pix)
        #    )
        #    plt.figure()
        #    plt.plot(
        #        dates_dtm,this_default,
        #        c='r',linewidth=.8
        #    )
        #    plt.xlabel('Date')
        #    plt.ylabel('temperature')
        #    plt.savefig(this_savename,dpi=300,bbox_inches='tight')
        #    plt.close()
        #for p,pix in enumerate(pixels):
        #    print('plotting q timeseries for pixel {}'.format(pix))
        #    this_default = default_timeseries_q[pix]
        #    this_savename = os.path.join(
        #        plots_dir,'q_timeseries_{}_{}.png'.format(this_exp,pix)
        #    )
        #    plt.figure()
        #    plt.plot(
        #        dates_dtm,this_default,
        #        c='r',linewidth=.8
        #    )
        #    plt.xlabel('Date')
        #    plt.ylabel('Q')
        #    plt.savefig(this_savename,dpi=300,bbox_inches='tight')
        #    plt.close()
    def plot_pso_maps(self,exp_names,default_df,pso_df,plots_dir,
                      plot_trim,extent='conus'
                     ):
        # define the lats and lons for the points
        lons = default_df.loc['lon']
        lons = lons.drop(labels=['all'])
        lats = default_df.loc['lat']
        lats - lats.drop(labels=['all'])
        # now get the different values that we want to plot in the different
        # plots
        # default le rmse
        default_le_rmse = default_df.loc['le_rmse']
        avg_default_le_rmse = default_le_rmse['all']
        # default le r2
        default_le_r2 = default_df.loc['le_r2']
        avg_default_le_r2 = default_le_r2['all']
        # default le corr
        default_le_corr = default_df.loc['le_corr']
        avg_default_le_corr = default_le_corr['all']
        # default le ubrmse
        default_le_ubrmse = default_df.loc['le_ubrmse']
        avg_defaut_le_ubrmse = default_le_ubrmse['all']
        # final pso rmse
        pso_le_rmse = pso_df.loc['le_rmse']
        avg_pso_le_rmse = pso_le_rmse['all']
        # final pso le r2
        pso_le_r2 = pso_df.loc['le_r2']
        avg_pso_le_r2 = pso_le_r2['all']
        # final pso le corr
        pso_le_corr = pso_df.loc['le_corr']
        avg_pso_le_corr = pso_le_corr['all']
        # final pso le ubrmse
        pso_le_ubrmse = pso_df.loc['le_ubrmse']
        avg_pso_le_ubrmse = pso_le_ubrmse['all']
        # first iteration pso rmse
        pso_first_it_le_rmse = pso_df.loc['le_rmse_first_iteration']
        avg_pso_first_it_le_rmse = pso_first_it_le_rmse['all']
        # perc change rmse
        perc_change_le_rmse = pso_df.loc['perc_change_le_rmse']
        avg_perc_change_le_rmse = perc_change_le_rmse['all']
        # change in r2
        change_le_r2 = pso_df.loc['change_le_r2']
        avg_change_le_r2 = change_le_r2['all']
        # change in corr
        change_le_corr = pso_df.loc['change_le_corr']
        avg_change_le_corr = change_le_corr['all']
        # perc change in ubrmse
        perc_change_le_ubrmse = pso_df.loc['perc_change_le_ubrmse']
        avg_perc_change_le_ubrmse = perc_change_le_ubrmse['all']
        # perc change rmse first iteration
        perc_change_le_rmse_first_it = pso_df.loc['perc_change_le_rmse_first_it']
        avg_perc_change_le_rmse_first_it = perc_change_le_rmse_first_it['all']
        # difference in le between fluxcom and average
        diff_le_default = default_df.loc['ave_le_diff']
        avg_diff_le_default = diff_le_default['all']
        # difference in le between fluxcom and PSO
        diff_le_pso = pso_df.loc['ave_le_diff']
        avg_diff_le_pso = diff_le_pso['all']
        # change in ubrmse
        change_le_ubrmse = pso_df.loc['change_le_ubrmse']
        avg_change_le_ubrmse = change_le_ubrmse['all']
        # change in ubrmse norm
        change_le_ubrmse_norm = pso_df.loc['change_le_ubrmse_norm']
        avg_change_le_ubrmse_norm = change_le_ubrmse_norm['all']
        # combine all relevant information into lists to be plotted:
            # exps are the name of the experiments
            # vals are the values of the actual poitns
            # avgs are the averages computed above
            # vmins are the mins for the colorbar
            # vmaxs are the max for the colorbars
            # names are the names of the statistic that we are plotting
        if not plot_trim:
            exps = [
                exp_names[0],exp_names[2],exp_names[1],exp_names[2],exp_names[1],
                exp_names[0],exp_names[0],exp_names[0],exp_names[2],exp_names[2],
                exp_names[2],exp_names[2],exp_names[2],exp_names[2],exp_names[2],
                exp_names[2],exp_names[2],exp_names[2]
            ]
            vals = [
                default_le_rmse,pso_le_rmse,pso_first_it_le_rmse,perc_change_le_rmse,
                perc_change_le_rmse_first_it,diff_le_default,default_le_r2,
                default_le_corr,pso_le_r2,pso_le_corr,change_le_r2,change_le_corr,
                default_le_ubrmse,pso_le_ubrmse,perc_change_le_ubrmse,
                diff_le_pso,change_le_ubrmse,change_le_ubrmse_norm
            ]
            avgs = [
                avg_default_le_rmse,avg_pso_le_rmse,avg_pso_first_it_le_rmse,
                avg_perc_change_le_rmse,
                avg_perc_change_le_rmse_first_it,avg_diff_le_default,
                avg_default_le_r2,avg_default_le_corr,avg_pso_le_r2,
                avg_pso_le_corr,avg_change_le_r2,avg_change_le_corr,
                avg_default_le_rmse,avg_pso_le_ubrmse,
                avg_perc_change_le_ubrmse,avg_diff_le_pso,
                avg_change_le_ubrmse,avg_change_le_ubrmse_norm
            ]
            names = [
                'default_le_rmse','pso_le_rmse','pso_first_it_le_rmse',
                'perce_change_le_rmse','perc_change_le_rmse_first_it',
                'diff_le_default','default_le_r2','default_le_corr',
                'pso_le_r2','pso_le_corr','change_le_r2','change_le_corr',
                'default_le_ubrmse','pso_le_ubrmse',
                'perc_change_le_ubrmse','diff_le_pso','change_le_ubrmse',
                'change_le_ubrmse_norm'
            ]
            types = [
                'le_rmse','le_rmse','le_rmse','le_perc_change','le_perc_change',
                'le_diff','le_perc_change','le_perc_change','le_perc_change',
                'le_perc_change','le_perc_change','le_perc_change',
                'le_rmse','le_rmse','le_perc_change','le_diff','diff_le_rmse',
                'diff_le_rmse_norm'
            ]
            cmaps = {
                'le_rmse':'rainbow',
                'le_perc_change':'bwr',
                'le_diff':'bwr',
                'diff_le_rmse':'bwr',
                'diff_le_rmse_norm':'bwr'
            }
            vmins = {
                'le_rmse':0,
                'le_perc_change':-1,
                'le_diff':-40,
                'diff_le_rmse':-20,
                'diff_le_rmse_norm':-.5
            }
            vmaxs = {
                'le_rmse':50,
                'le_perc_change':1,
                'le_diff':40,
                'diff_le_rmse':20,
                'diff_le_rmse_norm':.5
            }
            compare = {
                'le_rmse':'fluxcom',
                'le_perc_change':'fluxcom',
                'le_diff':'fluxcom',
                'diff_le_rmse':'fluxcom',
                'diff_le_rmse_norm':'fluxcom'
            }
        if plot_trim:
            exps = [
                exp_names[2],
                exp_names[2],
                exp_names[2],
                exp_names[2]
            ]
            vals = [
                perc_change_le_rmse,
                change_le_r2,
                change_le_corr,
                perc_change_le_ubrmse
            ]
            avgs = [
                avg_perc_change_le_rmse,
                avg_change_le_r2,
                avg_change_le_corr,
                avg_perc_change_le_ubrmse
            ]
            names = [
                'perc_change_le_rmse',
                'change_le_r2',
                'change_le_corr',
                'perc_change_le_ubrmse'
            ]
            types = names
            cmaps = {
                'perc_change_le_rmse':'bwr',
                'change_le_r2':'bwr',
                'change_le_corr':'bwr',
                'perc_change_le_ubrmse':'bwr'
            }
            vmins = {
                'perc_change_le_rmse':-.75,
                'change_le_r2':-.75,
                'change_le_corr':-.4,
                'perc_change_le_ubrmse':-.75
            }
            vmaxs = {
                'perc_change_le_rmse':.75,
                'change_le_r2':.75,
                'change_le_corr':.4,
                'perc_change_le_ubrmse':.75
            }


        for p in range(len(vals)):
            # let's first plot the rmse of default experiment versus fluxcom
            # create the figure and set the projection
            fig = plt.figure()
            ax = plt.axes(projection=ccrs.PlateCarree())
            # add coastline and set the limits
            ax.add_feature(cfeature.COASTLINE)
            ax.add_feature(cfeature.STATES,linewidth=0.1)
            ax.set_ylim(bottom=20,top=55)
            ax.set_xlim(left=-127,right=-60)
            # set extent isn't working with this weird bootleg version of cartopy
            # i've had to install on AWS--if it is fixed then could go back to this
            #if extent == 'conus':
            #    ax.set_extent([-127,-60,20,55],crs=ccrs.PlateCarree)
            #elif extent == 'global':
            #    ax.set_global()
            # define lats, lons, values
            lons = default_df.loc['lon']
            lats = default_df.loc['lat']
            val = vals[p]
            # make the scatterplot
            scatter = plt.scatter(
                lons,lats,marker='s',s=1,c=val,transform=ccrs.PlateCarree(),
                cmap=cmaps[types[p]],vmin=vmins[types[p]],vmax=vmaxs[types[p]]
            )
            # add colorbar
            plt.colorbar(scatter)
            # add the average error as text
            ax.text(
                -127+2,20+2,'Average {name}: {val:.2f}'.format(
                    name=names[p],val=avgs[p]
                ),
                bbox=dict(facecolor='white')
            )
            # save
            savename = '{name}_{exp}_iterations.png'.format(
                name=names[p],exp=exps[p]
            )
            savename = os.path.join(
                plots_dir,savename
            )
            plt.savefig(savename,dpi=300,bbox_inches='tight')
            plt.close()
        print('created maps of error and percent change')
    def plot_pso_progression(self,exp_names,default_df,pso_df,fluxcom_df,
                             plots_dir,pft_info):
        # get the unique pfts that we will be plotting for
        all_pfts = pft_info['primary_pft']
        unique_pfts = list(set(all_pfts))
        unique_pfts_idx = np.arange(len(unique_pfts))
        # get the average LE values and average RMSE values for each PFT
        # and get those for initial, first PFT iteration, final PFT iteration
        # also do this for the fluxcom observations
        # first set up the arrays that will hold these values
        default_means = np.zeros(len(unique_pfts))
        step_1_means = np.zeros(len(unique_pfts))
        pso_means = np.zeros(len(unique_pfts))
        fluxcom_means = np.zeros(len(unique_pfts))
        default_le_rmses = np.zeros(len(unique_pfts))
        step_1_le_rmses = np.zeros(len(unique_pfts))
        pso_le_rmses = np.zeros(len(unique_pfts))
        for p,pft in enumerate(unique_pfts):
            this_pft_idx = np.where(pft_info['primary_pft']==pft)
            pft_tiles = list(pft_info.index[this_pft_idx])
            # get the information for default mean le
            default_pft_mean_le_all = default_df[pft_tiles].loc['ave_le']
            default_pft_mean_le = np.mean(default_pft_mean_le_all)
            default_means[p] = default_pft_mean_le
            # get the information for step 1 mean le
            step_1_pft_mean_le_all = pso_df[pft_tiles].loc['ave_le_first_it']
            step_1_pft_mean_le = np.mean(step_1_pft_mean_le_all)
            step_1_means[p] = step_1_pft_mean_le
            # get the information for final pso mean le
            pso_pft_mean_le_all = pso_df[pft_tiles].loc['ave_le']
            pso_pft_mean_le = np.mean(pso_pft_mean_le_all)
            pso_means[p] = pso_pft_mean_le
            # get the information for the fluxcom observations mean le
            fluxcom_pft_mean_le_all = fluxcom_df[pft_tiles].loc['ave_le']
            fluxcom_pft_mean_le = np.mean(fluxcom_pft_mean_le_all)
            fluxcom_means[p] = fluxcom_pft_mean_le
            # get the information for default le rmse
            default_pft_le_rmse_all = default_df[pft_tiles].loc['le_rmse']
            default_pft_le_rmse = np.mean(default_pft_le_rmse_all)
            default_le_rmses[p] = default_pft_le_rmse
            # get the information for step 1 le rmse
            step_1_pft_le_rmse_all = pso_df[pft_tiles].loc[
                'le_rmse_first_iteration'
            ]
            step_1_pft_le_rmse = np.mean(step_1_pft_le_rmse_all)
            step_1_le_rmses[p] = step_1_pft_le_rmse
            # get the information for final pso le rmse
            pso_pft_le_rmse_all = pso_df[pft_tiles].loc['le_rmse']
            pso_pft_le_rmse = np.mean(pso_pft_le_rmse_all)
            pso_le_rmses[p] = pso_pft_le_rmse
        # the names of the two plots that we will make
        this_exp = exp_names[2]
        plot_names = ['ave_le','le_rmse']
        # dictionary that lists the values to plot for each of the plot names
        y_values = {
            'ave_le':[
                default_means,step_1_means,pso_means,fluxcom_means
            ],
            'le_rmse':[
                default_le_rmses,step_1_le_rmses,pso_le_rmses
            ]
        }
        # dictionary that lists the colors for the different y values
        colors = {
            'ave_le':[
                'r','y','g','k'
            ],
            'le_rmse':[
                'r','y','g'
            ]
        }
        # dictionary that lists the legend labels for the different exps
        labels = {
            'ave_le':[
                'Default Catchment-CN','PSO after first iteration',
                'Converged PSO','Fluxcom Observations'
            ],
            'le_rmse':[
                'Default Catchment-CN','PSO after first iteration',
                'Converged PSO'
            ]
        }
        y_labels = {
            'ave_le':'Average LE (W/m2)',
            'le_rmse':'LE RMSE (W/m2)'
        }
        for p,plot in enumerate(plot_names):
            fig,axs = plt.subplots()
            # plot the three different mean LE's
            this_y_vals = y_values[plot]
            this_colors = colors[plot]
            this_labels = labels[plot]
            this_y_label = y_labels[plot]
            for v,val in enumerate(this_y_vals):
                plt.scatter(
                    unique_pfts_idx,val,c=this_colors[v],
                    label=this_labels[v]
                )
            # change the axis names to be the different pfts
            plt.xticks(unique_pfts_idx,unique_pfts)
            # add a legend that says the different experiments
            plt.legend()
            # add x and y labels
            plt.ylabel(this_y_label)
            plt.xlabel('Plant Functional Type')
            # save the figure
            savename = '{plot}_change_pso_{exp}.png'.format(
                plot=plot,exp=this_exp
            )
            savename = os.path.join(
                plots_dir,savename
            )
            plt.savefig(savename)
        print('created scatter plots')
    def plot_general_comparison_maps(self,exp_names,default_df,exp_df,plots_dir,
                                     tile_gdf_fname,pixel_error_gdf_fname,extent='conus'
                                    ):
        # define the lats and lons for the points
        lons = default_df.loc['lon']
        lons = lons.drop(labels=['all'])
        lats = default_df.loc['lat']
        lats - lats.drop(labels=['all'])
        # now get the different values that we want to plot in the different
        # plots
        # default rmse
        default_le_rmse = default_df.loc['le_rmse']
        avg_default_le_rmse = default_le_rmse['all']
        # experiment rmse
        exp_le_rmse = exp_df.loc['le_rmse']
        avg_exp_le_rmse = exp_le_rmse['all']
        # difference between the two le's
        diff_le_runs = exp_df.loc['ave_le'] - default_df.loc['ave_le']
        avg_diff_le_runs = diff_le_runs['all']
        perc_diff_le_runs = diff_le_runs/default_df.loc['ave_le']
        avg_perc_diff_le_runs =  perc_diff_le_runs['all']
        # difference between the two rmse's
        diff_le_rmse_runs = exp_df.loc['le_rmse'] - default_df.loc['le_rmse']
        avg_diff_le_rmse_runs = diff_le_rmse_runs['all']
        perc_diff_le_rmse = diff_le_rmse_runs/default_df.loc['le_rmse']
        avg_perc_diff_le_rmse = perc_diff_le_rmse['all']
        # difference between the two ave_sm
        diff_ave_sm_runs = (
            exp_df.loc['ave_ave_sm'] - default_df.loc['ave_ave_sm']
        )
        avg_diff_ave_sm_runs = diff_ave_sm_runs['all']
        perc_diff_ave_sm_runs = diff_ave_sm_runs/default_df.loc['ave_ave_sm']
        avg_perc_diff_ave_sm_runs = perc_diff_ave_sm_runs['all']
        # difference between the two root_sm
        diff_root_sm_runs = (
            exp_df.loc['ave_root_sm'] - default_df.loc['ave_root_sm']
        )
        avg_diff_root_sm_runs = diff_root_sm_runs['all']
        perc_diff_root_sm_runs = diff_root_sm_runs/default_df.loc['ave_root_sm']
        avg_perc_diff_root_sm_runs = perc_diff_root_sm_runs['all']
        # difference between the two surf_sm
        diff_surf_sm_runs = (
            exp_df.loc['ave_surf_sm'] - default_df.loc['ave_surf_sm']
        )
        avg_diff_surf_sm_runs = diff_surf_sm_runs['all']
        perc_diff_surf_sm_runs = diff_surf_sm_runs/default_df.loc['ave_surf_sm']
        avg_perc_diff_surf_sm_runs = perc_diff_surf_sm_runs['all']
        # difference between the two infil
        diff_infil_runs = (
            exp_df.loc['ave_infil'] - default_df.loc['ave_infil']
        )
        avg_diff_infil_runs = diff_infil_runs['all']
        perc_diff_infil_runs = diff_infil_runs/default_df.loc['ave_infil']
        avg_perc_diff_infil_runs = perc_diff_infil_runs['all']
        # difference between the two runoff
        diff_runoff_runs = (
            exp_df.loc['ave_runoff'] - default_df.loc['ave_runoff']
        )
        avg_diff_runoff_runs = diff_runoff_runs['all']
        perc_diff_runoff_runs = diff_runoff_runs/default_df.loc['ave_runoff']
        avg_perc_diff_runoff_runs = perc_diff_runoff_runs['all']
        # difference between the two baseflow
        diff_baseflow_runs = (
            exp_df.loc['ave_baseflow'] - default_df.loc['ave_baseflow']
        )
        avg_diff_baseflow_runs = diff_baseflow_runs['all']
        perc_diff_baseflow_runs = diff_baseflow_runs/default_df.loc[
            'ave_baseflow'
        ]
        avg_perc_diff_baseflow_runs = perc_diff_baseflow_runs['all']
        # difference between changes in runoff and changes in baseflow
        diff_baseflow_runoff_runs = (
            np.abs(diff_baseflow_runs) - np.abs(diff_runoff_runs)
        )
        avg_diff_baseflow_runoff_runs = diff_baseflow_runoff_runs['all']
        norm_diff_baseflow_runoff_runs = (
            diff_baseflow_runoff_runs/(
                default_df.loc['ave_baseflow'] + default_df.loc['ave_runoff']
            )
        )
        avg_norm_diff_baseflow_runoff_runs = norm_diff_baseflow_runoff_runs[
            'all'
        ]
        # let's convert mm/day for this presentation
        diff_le_rmse_runs_mm_day = diff_le_rmse_runs/28.94
        avg_diff_le_rmse_runs_mm_day = diff_le_rmse_runs_mm_day['all']
        # put values to be plotted into list for plotting
        vals = [
            diff_le_runs,diff_le_rmse_runs,
            diff_le_rmse_runs_mm_day,diff_ave_sm_runs,
            diff_root_sm_runs,diff_surf_sm_runs,diff_infil_runs,diff_runoff_runs,
            diff_baseflow_runs,perc_diff_le_runs,
            perc_diff_le_rmse,perc_diff_ave_sm_runs,
            perc_diff_root_sm_runs,perc_diff_surf_sm_runs,perc_diff_infil_runs,
            perc_diff_runoff_runs,perc_diff_baseflow_runs,
            diff_baseflow_runoff_runs,norm_diff_baseflow_runoff_runs
        ]
        # put the averages that correspond to these values
        avgs = [
            avg_diff_le_runs,avg_diff_le_rmse_runs,
            avg_diff_le_rmse_runs_mm_day,avg_diff_ave_sm_runs,
            avg_diff_root_sm_runs,avg_diff_surf_sm_runs,
            avg_diff_infil_runs,avg_diff_runoff_runs,
            avg_diff_baseflow_runs,avg_perc_diff_le_runs,
            avg_perc_diff_le_rmse,
            avg_perc_diff_ave_sm_runs,avg_perc_diff_root_sm_runs,
            avg_perc_diff_surf_sm_runs,avg_perc_diff_infil_runs,
            avg_perc_diff_runoff_runs,avg_perc_diff_baseflow_runs,
            avg_diff_baseflow_runoff_runs,
            avg_norm_diff_baseflow_runoff_runs
        ]
        # put the name the corresponds to each of these values
        names = [
            'diff_le','diff_le_rmse',
            'diff_le_rmse_mm_day','diff_ave_sm','diff_root_sm','diff_surf_sm',
            'diff_infil',
            'diff_runoff','diff_baseflow','perc_diff_le_runs',
            'perc_diff_le_rmse',
            'perc_diff_ave_sm','perc_diff_root_sm',
            'perc_diff_surf_sm','perc_diff_infil','perc_diff_runoff',
            'perc_diff_baseflow','diff_baseflow_runoff',
            'norm_diff_baseflow_runoff'
        ]
        # vals, avgs, and names need to all be the same length. if this isn't
        # true stop here and inform the user
        len_vals = len(vals)
        len_avgs = len(avgs)
        len_names = len(names)
        if ((len_vals != len_avgs) or
            (len_vals != len_names) or
            (len_avgs != len_names)):
            print('vals, avgs, and names must all be the same length!')
            print('go back and correct this!')
            sys.exit()
        types = [
            'le_diff','le_rmse_diff',
            'le_rmse_diff_mm_day','sm_diff','sm_diff',
            'sm_diff','flow_diff','flow_diff','flow_diff',
            'perc_diff_le','perc_diff_le_rmse',
            'perc_diff_sm','perc_diff_sm','perc_diff_sm',
            'perc_diff_flow','perc_diff_flow','perc_diff_flow',
            'flow_diff','perc_diff_flow'
        ]
        cmaps = {
            'le_diff':'bwr',
            'le_rmse_diff':'bwr',
            'le_rmse_diff_mm_day':'PiYG_r',
            'sm_diff':'bwr',
            'flow_diff':'bwr',
            'perc_diff_le':'bwr',
            'perc_diff_le_rmse':'bwr',
            'perc_diff_sm':'bwr',
            'perc_diff_flow':'bwr'
        }
        vmins = {
            'le_diff':-20,
            'le_rmse_diff':-20,
            'le_rmse_diff_mm_day':-.75,
            'sm_diff':-.5,
            'flow_diff':-0.00005,
            'perc_diff_le':-0.3,
            'perc_diff_le_rmse':-.6,
            'perc_diff_sm':-0.75,
            'perc_diff_flow':-0.75
        }
        vmaxs = {
            'le_diff':20,
            'le_rmse_diff':20,
            'le_rmse_diff_mm_day':.75,
            'sm_diff':.5,
            'flow_diff':0.00005,
            'perc_diff_le':0.3,
            'perc_diff_le_rmse':.6,
            'perc_diff_sm':0.75,
            'perc_diff_flow':0.75
        }
        # lets save all of these things to a geopandas for nice use wiht qgis
        tile_gdf = gpd.read_file(
            tile_gdf_fname
        )
        print('min of gdf:')
        print(np.min(diff_le_rmse_runs_mm_day))
        print('max of gdf:')
        print(np.max(diff_le_rmse_runs_mm_day))
        tile_gdf['err'] = np.array(diff_le_rmse_runs_mm_day)[:-1]
        norm = mpl.colors.Normalize(
            vmin=vmins['le_rmse_diff_mm_day'],vmax=vmaxs['le_rmse_diff_mm_day']
        )
        this_cmap = mpl.cm.get_cmap(cmaps['le_rmse_diff_mm_day'])
        this_vals_norm = norm(np.array(diff_le_rmse_runs_mm_day)[:-1])
        this_colors = this_cmap(this_vals_norm)
        this_colors[:,:-1] = this_colors[:,:-1]*255
        this_colors_str = []
        #for c in range(np.shape(this_colors)[0]):
        #    this_str = '{} {} {}'.format(
        #        round(this_colors[c][0]),
        #        round(this_colors[c][1]),
        #        round(this_colors[c][2])
        #    )
        #    this_colors_str.append(this_str)
        tile_gdf['err_r'] = this_colors[:,0]
        tile_gdf['err_g'] = this_colors[:,1]
        tile_gdf['err_b'] = this_colors[:,2]
        ch = choose()
        ch.save_gdf(tile_gdf,pixel_error_gdf_fname)
        for p in range(len(vals)):
            # let's first plot the rmse of default experiment versus fluxcom
            # create the figure and set the projection
            fig = plt.figure()
            ax = plt.axes(projection=ccrs.PlateCarree())
            # add coastline and set the limits
            ax.add_feature(cfeature.COASTLINE)
            ax.add_feature(cfeature.STATES,linewidth=0.1)
            ax.set_ylim(bottom=20,top=55)
            ax.set_xlim(left=-127,right=-60)
            # set extent isn't working with this weird bootleg version of cartopy
            # i've had to install on AWS--if it is fixed then could go back to this
            #if extent == 'conus':
            #    ax.set_extent([-127,-60,20,55],crs=ccrs.PlateCarree)
            #elif extent == 'global':
            #    ax.set_global()
            # define lats, lons, values
            lons = default_df.loc['lon']
            lats = default_df.loc['lat']
            val = vals[p]
            # make the scatterplot
            scatter = plt.scatter(
                lons,lats,marker='s',s=1,c=val,transform=ccrs.PlateCarree(),
                cmap=cmaps[types[p]],vmin=vmins[types[p]],vmax=vmaxs[types[p]]
            )
            # add colorbar
            plt.colorbar(scatter)
            # add the average error as text
            ax.text(
                -127+2,20+2,'Average {name}: {val:.2f}'.format(
                    name=names[p],val=avgs[p]
                ),
                bbox=dict(facecolor='white')
            )
            # save
            savename = '{name}_{exp}_vs_{comp}_optimization_pixels.png'.format(
                name=names[p],exp=exp_names[0],comp=exp_names[1]
            )
            savename = os.path.join(
                plots_dir,savename
            )
            plt.savefig(savename,dpi=300,bbox_inches='tight')
            plt.close()
