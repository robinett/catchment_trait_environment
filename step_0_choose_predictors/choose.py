import xarray as xr
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import xesmf as xe
import rioxarray as rio
import sys
import os
import copy
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

def add_canopy_height(
    orig_fname,
    regrid_fname,
    g1_ds,
    regrid=True
):
    if regrid == False:
        h_regrid = xr.open_dataset(
            regrid_fname
        )
    else:
        h = xr.open_dataset(orig_fname)
        target_grid = g1_ds[["lat","lon"]]
        regridder = xe.Regridder(
            h,target_grid,method='bilinear'
        )
        h_regrid = regridder(h)
        h_vals = h_regrid['Band1'].values
        h_vals = np.where(h_vals == 0,np.nan,h_vals)
        h_inv = 1/h_vals
        h_regrid['inv_Band1'] = (('lat','lon'),h_inv)
        h_regrid = norm_ds_var(
            h_regrid,
            'Band1'
        )
        h_regrid = norm_ds_var(
            h_regrid,
            'inv_Band1'
        )
        h_regrid.to_netcdf(regrid_fname)
        plot_map_from_xarray(
            h_regrid,
            'Band1_norm',
            (
                '/discover/nobackup/trobinet/from_aws/pso/' +
                'step_0_choose_predictors/plots/h_norm.png'
            )
        )
        plot_map_from_xarray(
            h_regrid,
            'inv_Band1_norm',
            (
                '/discover/nobackup/trobinet/from_aws/pso/' +
                'step_0_choose_predictors/plots/inv_h_norm.png'
            )
        )
    g1_ds['h_norm'] = h_regrid['Band1_norm']
    g1_ds['inv_h_norm'] = h_regrid['inv_Band1_norm']
    return g1_ds
def add_species_richness(
    orig_fname,
    regrid_fname,
    g1_ds,
    regrid=True
):
    if regrid == False:
        sr = xr.open_dataset(
            regrid_fname
        )
        g1_ds['sr_norm'] = sr['sr_norm']
        return g1_ds
    sr = pd.read_csv(orig_fname)
    # coords and values
    lon = sr['X']
    lat = sr['Y']
    vals = sr['N']
    # target grid to interpolate
    lon_t = g1_ds['lon'].values
    lat_t = g1_ds['lat'].values
    lon_t_mesh,lat_t_mesh = np.meshgrid(lon_t,lat_t)
    # interpolate
    vals_interp = griddata(
        (lon,lat),
        vals,
        (lon_t_mesh,lat_t_mesh),
        method='linear'
    )
    data_array = xr.DataArray(
        vals_interp,
        coords={'lat':lat_t,'lon':lon_t},
        dims=['lat','lon'],
        name='sr'
    )
    sr_ds = xr.Dataset({'sr':data_array})
    sr_ds = norm_ds_var(
        sr_ds,
        'sr'
    )
    sr_ds.to_netcdf(regrid_fname)
    plot_map_from_xarray(
        sr_ds,
        'sr_norm',
        (
            '/discover/nobackup/trobinet/from_aws/pso/' +
            'step_0_choose_predictors/plots/sr_norm.png'
        )
    )
    # add to g1_ds
    g1_ds['sr_norm'] = sr_ds['sr_norm']
    return g1_ds
def add_tree_perc(
    modis_in_fname,
    modis_out_fname,
    g1_ds,
    regrid=True
):
    if regrid == False:
        tp = xr.open_dataset(
            modis_out_fname
        )
        g1_ds['tp_norm'] = tp['tp_w_nan_norm']
        return g1_ds
    tp = rio.open_rasterio(modis_in_fname)
    tree_perc = tp.isel(band=0)
    tree_perc_ds = tree_perc.to_dataset(name='tp')
    tree_perc_ds = tree_perc_ds.rename({'y':'lat','x':'lon'})
    target_grid = g1_ds[["lat","lon"]]
    regridder = xe.Regridder(
        tree_perc_ds,target_grid,method='bilinear'
    )
    tp_regrid = regridder(tree_perc_ds)
    # turn ocean to nan
    tp_regrid['tp_w_nan'] = tp_regrid['tp'].where(
        ~np.isnan(g1_ds['Tair_f_inst_mean_norm'])
    )
    tp_regrid = norm_ds_var(
        tp_regrid,
        'tp_w_nan'
    )
    tp_regrid.to_netcdf(modis_out_fname)
    plot_map_from_xarray(
        tp_regrid,
        'tp_w_nan_norm',
        (
            '/discover/nobackup/trobinet/from_aws/pso/' +
            'step_0_choose_predictors/plots/tp_norm.png'
        )
    )
    g1_ds['tp_norm'] = tp_regrid['tp_w_nan_norm']
def add_gldas_vars(
    gldas_fname,
    g1_ds,
    plot=False
):
    # open gldas file
    gldas = xr.open_mfdataset(
        os.path.join(
            '/discover/nobackup/trobinet/from_aws/pso/step_0_choose_predictors',
            'data',
            'gldas',
            '*.nc4'
        )
    )
    # average relevant vars over 30-year time period
    vars_to_average = [
        'Tair_f_inst',
        'Rainf_tavg',
        'Snowf_tavg',
        'Evap_tavg',
        'PotEvap_tavg'
    ]
    meaned_names = []
    normed_names = []
    for v,var in enumerate(vars_to_average):
        this_var = gldas[var].mean(dim='time')
        this_name_avg = copy.deepcopy(var) + '_mean'
        this_name_norm = copy.deepcopy(this_name_avg) + '_norm'
        gldas[this_name_avg] = this_var
        meaned_names.append(this_name_avg)
        normed_names.append(this_name_norm)
    # now for the vars we need to manually compute
    # first rainf + snowf
    rainf = gldas['Rainf_tavg'].mean(dim='time')
    snowf = gldas['Snowf_tavg'].mean(dim='time')
    rainf_snowf = rainf + snowf
    this_name = 'Rainf_tavg_Snowf_tavg'
    this_name_avg = copy.deepcopy(this_name) + '_mean'
    this_name_norm = copy.deepcopy(this_name_avg) + '_norm'
    gldas[this_name_avg] = rainf_snowf
    meaned_names.append(this_name_avg)
    normed_names.append(this_name_norm)
    # next p_pet
    p_pet = (
        gldas['Rainf_tavg_Snowf_tavg_mean']/
        gldas['PotEvap_tavg_mean']
    )
    this_name = 'Rainf_Snowf_over_PotEvap'
    this_name_avg = copy.deepcopy(this_name) + '_mean'
    this_name_norm = copy.deepcopy(this_name_avg) + '_norm'
    gldas[this_name_avg] = p_pet
    meaned_names.append(this_name_avg)
    normed_names.append(this_name_norm)
    # finally pet_et
    pet_et = (
        gldas['PotEvap_tavg_mean'] -
        gldas['Evap_tavg_mean']
    )
    this_name = 'PotEvap_minus_Evap'
    this_name_avg = copy.deepcopy(this_name) + '_mean'
    this_name_norm = copy.deepcopy(this_name_avg) + '_norm'
    gldas[this_name_avg] = pet_et
    meaned_names.append(this_name_avg)
    normed_names.append(this_name_norm)
    for v,var in enumerate(meaned_names):
        gldas = norm_ds_var(
            gldas,
            var
        )
        if plot:
            plot_map_from_xarray(
                gldas,
                normed_names[v],
                (
                    '/discover/nobackup/trobinet/from_aws/pso/' +
                    'step_0_choose_predictors/plots/{}.png'.format(
                        normed_names[v]
                    )
                ),
                vmin=-4,
                vmax=4
            )
        g1_ds[normed_names[v]] = gldas[normed_names[v]]
def norm_ds_var(ds,var_name):
    this_var_mean = ds[var_name].mean()
    this_var_std = ds[var_name].std()
    this_var_name_new = var_name + '_norm'
    ds[this_var_name_new] = (
        ds[var_name] - this_var_mean
    )/this_var_std
    return ds
def plot_map_from_xarray(
    ds,
    var_name,
    save_name,
    vmin=np.nan,
    vmax=np.nan
):
    projection = ccrs.PlateCarree()
    fig,ax = plt.subplots(subplot_kw={'projection':projection})
    ds[var_name].plot(
        ax=ax,
        transform=ccrs.PlateCarree()
    )
    ax.add_feature(cfeature.COASTLINE)
    plt.savefig(
        save_name,
        dpi=300,
        bbox_inches='tight'
    )
    plt.close()
def plot_corr_mat(
    corr_matrix,
    fname,
    bold_thresh=1.1
):
    annot = np.array([
        [
            format_bold_func(corr_matrix.iloc[i,j],bold_thresh)
            for j in range(corr_matrix.shape[1])
        ]
        for i in range(corr_matrix.shape[0])
    ])

    plt.figure(figsize=(10,8))
    sns.heatmap(
        corr_matrix,
        annot=annot,
        fmt="",
        cmap='coolwarm',
        square=True,
        annot_kws={"fontsize": 10, "va": "center", "ha": "center"},
        vmin=-1,
        vmax=1
    )
    plt.savefig(
        fname,
        dpi=250,
        bbox_inches='tight'
    )
    plt.close()
def format_bold_func(value,thresh):
    if abs(value) > thresh:
        return f"*{value:.2f}*"
    return f"{value:.2f}"
def get_vif_and_fit_model(
    regression_df,
    x_cols,
    y_col
):
    predictor_df = copy.deepcopy(regression_df)
    predictor_df = predictor_df[x_cols]
    X = sm.add_constant(predictor_df)
    #vif_data = pd.DataFrame()
    #vif_data['Feature'] = X.columns
    #vif_data['VIF'] = [
    #    variance_inflation_factor(X.values,i) for i in range(X.shape[1])
    #]
    #print(vif_data)
    # fit regression model and check variable important
    model = sm.OLS(regression_df[y_col],X).fit()
    print(model.summary())
    # Coefficients table as DataFrame
    coefficients = pd.DataFrame({
        "p-value": model.pvalues,
        "aic":model.aic,
        "r2":model.rsquared
    })
    # Customize the display format
    pd.set_option("display.float_format", "{:.100f}".format)
    print(coefficients)
    # return predictor p-values and aic
    return coefficients
def aic_plot(
    aic_vals,
    num_predictors,
    fname
):
    plt.figure()
    nums = np.arange(len(aic_vals))
    plt.plot(nums,aic_vals)
    plt.xticks(ticks=nums,labels=num_predictors)
    plt.xlabel('number of predictors')
    plt.ylabel('AIC')
    plt.savefig(fname,dpi=350,bbox_inches='tight')
    plt.close()




