import netCDF4 as nc
import os
import xarray as xr
import numpy as np
import sys
import pandas as pd
import copy
import choose as ch

def main():
    main_dir = (
        '/discover/nobackup/trobinet/from_aws/pso/step_0_choose_predictors'
    )
    data_dir = os.path.join(
        main_dir,
        'data'
    )
    plots_dir = os.path.join(
        main_dir,
        'plots'
    )
    ## open g1 from Yanlan
    #print('getting g1')
    #g1 = xr.open_dataset(
    #    os.path.join(
    #        data_dir,
    #        'retrieved_g1.nc'
    #    )
    #)
    #print('getting canopy height')
    #ch.add_canopy_height(
    #    os.path.join(
    #        data_dir,
    #        'canopy_height.nc'
    #    ),
    #    os.path.join(
    #        data_dir,
    #        'canopy_height_regrid.nc'
    #    ),
    #    g1,
    #    regrid=True
    #)
    #print('getting species richness')
    #ch.add_species_richness(
    #    os.path.join(
    #        data_dir,
    #        'species_richness.csv'
    #    ),
    #    os.path.join(
    #        data_dir,
    #        'species_richness.nc'
    #    ),
    #    g1,
    #    regrid=False
    #)
    #print('getting gldas')
    #ch.add_gldas_vars(
    #    os.path.join(
    #        data_dir,
    #        'gldas',
    #        '*.nc4'
    #    ),
    #    g1,
    #    plot=False
    #)
    ## we need to mask out all ocean values to get accurate normalized
    ## tree percentage
    ## use everywhere that we have gldas values as this mask
    #print('getting tree %')
    #ch.add_tree_perc(
    #    os.path.join(
    #        data_dir,
    #        'modis_veg.tif'
    #    ),
    #    os.path.join(
    #        data_dir,
    #        'modis_veg_perc.nc'
    #    ),
    #    g1,
    #    regrid=False
    #)
    #print(g1)
    ## get the top 10% of g1 confidence according to interquartile range
    ## we only want to do this where all predictors also have values, so
    ## turn all g1 where not all predictors have value to nan
    ## now, at those predictors, find the relevant interquartile range and the
    ## acceptable pixels
    #print('getting confidence in g1')
    #g1['range'] = g1['g1_75'] - g1['g1_25']
    #range_data = g1['range']
    #g1_vals = g1['g1_50']
    #perc_10 = range_data.quantile(0.1)
    #g1['acceptable'] = g1_vals.where(range_data <= perc_10)
    #ch.plot_map_from_xarray(
    #    g1,
    #    'acceptable',
    #    os.path.join(
    #        plots_dir,
    #        'g1_acceptable_range.png'
    #    )
    #)
    ## now get the pixels where g1 is in the acceptable range and
    ## we have all environmental predictors
    #print('selecting regression pixels')
    #regression_vars = [
    #    'h_norm',
    #    'inv_h_norm',
    #    'sr_norm',
    #    'tp_norm',
    #    'Tair_f_inst_mean_norm',
    #    'Rainf_tavg_mean_norm',
    #    'Snowf_tavg_mean_norm',
    #    'Rainf_tavg_Snowf_tavg_mean_norm',
    #    'Rainf_Snowf_over_PotEvap_mean_norm',
    #    'PotEvap_minus_Evap_mean_norm',
    #    'acceptable'
    #]
    #filter_vars = [
    #    'h_norm',
    #    'sr_norm',
    #    'tp_norm',
    #    'Tair_f_inst_mean_norm',
    #    'acceptable'
    #]
    #regression_ds = g1[regression_vars]
    #filter_ds = g1[filter_vars]
    #valid_mask = filter_ds.notnull().to_array().all(dim='variable')
    #computed_mask = valid_mask.compute()
    #filtered_ds = regression_ds.where(computed_mask,drop=True)
    #ch.plot_map_from_xarray(
    #    filtered_ds,
    #    'acceptable',
    #    os.path.join(
    #        plots_dir,
    #        'g1_for_regression.png'
    #    )
    #)
    #print('turning regression variables to dataframe')
    #regression_df = filtered_ds.to_dataframe().reset_index()
    #regression_df = regression_df.dropna()
    #print(regression_df)
    #regression_df.to_csv(
    #    os.path.join(
    #        data_dir,
    #        'regression_vals.csv'
    #    ),
    #    index=False
    #)
    # check for multicolinearity between the predictors
    regression_df = pd.read_csv(
        os.path.join(
            data_dir,
            'regression_vals.csv'
        )
    )
    print(regression_df)
    print(regression_df.columns)
    all_predictors = [
        'h_norm',
        'inv_h_norm',
        'sr_norm',
        'tp_norm',
        'Tair_f_inst_mean_norm',
        'Rainf_tavg_Snowf_tavg_mean_norm',
        'Rainf_Snowf_over_PotEvap_mean_norm',
        'PotEvap_minus_Evap_mean_norm'
    ]
    all_pos_pred_df = regression_df[all_predictors]
    corr_mat = all_pos_pred_df.corr()
    ch.plot_corr_mat(
        corr_mat,
        os.path.join(
            plots_dir,
            'pred_corr_matrix.svg'
        ),
        bold_thresh=0.75
    )
    # okay we need to choose either h or h_inv.
    # let's see which one leads to a better fit
    predictors = [
        'h_norm'
    ]
    y_col = 'acceptable'
    model_out = ch.get_vif_and_fit_model(
        regression_df,
        predictors,
        y_col
    )
    all_predictors = [
        'inv_h_norm'
    ]
    y_col = 'acceptable'
    model_out = ch.get_vif_and_fit_model(
        regression_df,
        predictors,
        y_col
    )
    # using inv_h alone leads to a R2 of 0.147 while using h_norm
    # alone leads to a R2 of 0.219. So we will go with h.
    # alright, with only one h and one aridity
    predictors = [
        'h_norm',
        'sr_norm',
        'tp_norm',
        'Tair_f_inst_mean_norm',
        'Rainf_tavg_Snowf_tavg_mean_norm',
        'PotEvap_minus_Evap_mean_norm'
    ]
    model_out = ch.get_vif_and_fit_model(regression_df,predictors,y_col)
    aics = np.zeros(0)
    aics = np.append(aics,model_out['aic'].iloc[0])
    r2s = np.zeros(0)
    r2s = np.append(r2s,model_out['r2'].iloc[0])
    p_vals = model_out['p-value'].to_frame(name='1')
    # aridity is the least relevant
    predictors = [
        'h_norm',
        'sr_norm',
        'tp_norm',
        'Tair_f_inst_mean_norm',
        'Rainf_tavg_Snowf_tavg_mean_norm'
    ]
    model_out = ch.get_vif_and_fit_model(regression_df,predictors,y_col)
    aics = np.append(aics,model_out['aic'].iloc[0])
    r2s = np.append(r2s,model_out['r2'].iloc[0])
    p_vals['2'] = model_out['p-value']
    # now tair is the least relevant
    predictors = [
        'h_norm',
        'sr_norm',
        'tp_norm',
        'Rainf_tavg_Snowf_tavg_mean_norm'
    ]
    model_out = ch.get_vif_and_fit_model(regression_df,predictors,y_col)
    aics = np.append(aics,model_out['aic'].iloc[0])
    r2s = np.append(r2s,model_out['r2'].iloc[0])
    p_vals['3'] = model_out['p-value']
    # next, get rid of tp
    predictors = [
        'h_norm',
        'sr_norm',
        'Rainf_tavg_Snowf_tavg_mean_norm'
    ]
    model_out = ch.get_vif_and_fit_model(regression_df,predictors,y_col)
    aics = np.append(aics,model_out['aic'].iloc[0])
    r2s = np.append(r2s,model_out['r2'].iloc[0])
    p_vals['4'] = model_out['p-value']
    # then sr
    predictors = [
        'h_norm',
        'Rainf_tavg_Snowf_tavg_mean_norm'
    ]
    model_out = ch.get_vif_and_fit_model(regression_df,predictors,y_col)
    aics = np.append(aics,model_out['aic'].iloc[0])
    r2s = np.append(r2s,model_out['r2'].iloc[0])
    p_vals['5'] = model_out['p-value']
    # and try both one-predictor models
    predictors = [
        'h_norm'
    ]
    model_out = ch.get_vif_and_fit_model(regression_df,predictors,y_col)
    aics = np.append(aics,model_out['aic'].iloc[0])
    r2s = np.append(r2s,model_out['r2'].iloc[0])
    p_vals['6'] = model_out['p-value']
    # with h_norm and rainf snowf
    predictors = [
        'Rainf_tavg_Snowf_tavg_mean_norm'
    ]
    model_out = ch.get_vif_and_fit_model(regression_df,predictors,y_col)
    aics = np.append(aics,model_out['aic'].iloc[0])
    r2s = np.append(r2s,model_out['r2'].iloc[0])
    p_vals['7'] = model_out['p-value']
    # print what we found
    pd.set_option("display.float_format", "{:.2f}".format)
    print(aics)
    print(r2s)
    print(p_vals)
    # get the number of predictors for each corresponding aic value
    num_predictors = np.sum(~p_vals.isna(),axis=0).values
    num_predictors = num_predictors - 1
    print(num_predictors)
    ch.aic_plot(
        aics,
        num_predictors,
        os.path.join(
            plots_dir,
            'aic_plot.svg'
        )
    )
    ch.aic_plot(
        r2s,
        num_predictors,
        os.path.join(
            plots_dir,
            'r2_plot.svg'
        )
    )



if __name__ == '__main__':
    main()
