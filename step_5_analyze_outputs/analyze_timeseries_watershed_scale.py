import pickle
import sys
import pandas as pd
import datetime
import copy
import numpy as np
from dateutil import relativedelta
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import netCDF4 as nc
from scipy.stats import linregress
import geopandas as gpd
from mpl_toolkits.axes_grid1 import make_axes_locatable

class analyze_watershed:
    def __init__(self):
        pass
    def get_tiles(self,tile_fname):
        # get the tiles that we are going to run as defined in step 1
        tiles = pd.read_csv(tile_fname,header=None)
        # turn this into an np array
        tiles = np.array(tiles).astype(int)
        # make a nice np array
        tiles = tiles.T
        tiles = tiles[0]
        # save to self
        self.tiles = tiles
        # index for this tile would start at 0, so subtract 1 to get index
        tiles_idx = tiles - 1
        self.tiles_idx = tiles_idx
        return tiles
    def get_intersection(self,int_dir):
        with open(int_dir,'rb') as f:
            out = pickle.load(f)
        return out
    def get_model_preds_watershed(self,start,end,exp_names,catch_timeseries,
                                  streamflow_timeseries,
                                  fluxcom_timeseries,intersection_info,
                                  start_err,end_err,
                                  get_metf=False):
        '''
        function that gets RMSE dict at the monthly/watershed level
        we will do this for both streamflow as well as LE
        '''
        # one edit to fluxcom timeseries for simplicity
        fluxcom_timeseries = fluxcom_timeseries.set_index('time')
        # get the lits of watersheds as strings to be used as the dict keys for
        # intersection_info
        watersheds_str = list(intersection_info.keys())
        # get the list of watersheds as integers
        watersheds = [int(k) for k in intersection_info.keys()]
        # for each watershed let's get the monthly streamflow and LE:
            # from waterwatch/fluxcom respectively
            # from default
            # from step 1 pso
            # from final step pso
        # get what will be our columns for dfs
        cols = copy.deepcopy(watersheds)
        cols.append('all')
        # get the information that we need and create the arrays that will hold
        # the final watershed information
        num_watersheds = len(watersheds)
        delta = end - start
        num_days = delta.days + 1
        delta_err = end_err - start_err
        num_days_err = delta_err.days + 1
        # placeholders for everything that we need to extract and convert
        default_all_strm = np.zeros((num_days,num_watersheds))
        default_all_le = np.zeros((num_days,num_watersheds))
        default_all_runoff = np.zeros((num_days,num_watersheds))
        default_all_baseflow = np.zeros((num_days,num_watersheds))
        fluxcom_all_le = np.zeros((num_days_err,num_watersheds))
        pso_init_all_strm = np.zeros((num_days,num_watersheds))
        pso_init_all_le = np.zeros((num_days,num_watersheds))
        pso_final_all_strm = np.zeros((num_days,num_watersheds))
        pso_final_all_le = np.zeros((num_days,num_watersheds))
        if get_metf:
            default_all_rainfsnowf = np.zeros((num_days,num_watersheds))
        for w,wat in enumerate(watersheds):
            # get the waterwatch streamflow for this watershed
            this_strm_data = streamflow_timeseries[str(wat)]
            # get the tiles in this watershed
            this_tiles = intersection_info[watersheds_str[w]][0]
            # get the percent of each of these tiles in this watershed
            this_perc = intersection_info[watersheds_str[w]][1]
            # get the fluxcom-defined le for each of these tiles
            fluxcom_le = fluxcom_timeseries[this_tiles]
            # get the default runoff and le for each of these tiles and average
            # runoff to monthly resolution and correct units
            # let's do this for default first
            default_timeseries = catch_timeseries[exp_names[0]]
            default_runoff = default_timeseries['runoff']
            default_baseflow = default_timeseries['baseflow']
            default_strm = default_runoff + default_baseflow
            default_runoff = default_runoff[this_tiles]
            default_baseflow = default_baseflow[this_tiles]
            default_strm = default_strm[this_tiles]
            default_le = default_timeseries['le']
            default_le = default_le[this_tiles]
            #print('default_runoff')
            #print(default_runoff)
            #print('default_baseflow')
            #print(default_baseflow)
            #sys.exit()
            # for default met forcing
            if get_metf:
                # for incoming precip
                default_rainfsnowf = default_timeseries['rainfsnowf']
                default_rainfsnowf = default_rainfsnowf[this_tiles]
            # let's do this for pso init
            pso_init_timeseries = catch_timeseries[exp_names[1]]
            pso_init_runoff = pso_init_timeseries['runoff']
            pso_init_baseflow = pso_init_timeseries['baseflow']
            pso_init_strm = pso_init_runoff + pso_init_baseflow
            pso_init_strm = pso_init_strm[this_tiles]
            pso_init_le = pso_init_timeseries['le']
            pso_init_le = pso_init_le[this_tiles]
            # and dofor pso final
            pso_final_timeseries = catch_timeseries[exp_names[2]]
            pso_final_runoff = pso_final_timeseries['runoff']
            pso_final_baseflow = pso_final_timeseries['baseflow']
            pso_final_strm = pso_final_runoff + pso_final_baseflow
            pso_final_strm = pso_final_strm[this_tiles]
            pso_final_le = pso_final_timeseries['le']
            pso_final_le = pso_final_le[this_tiles]
            # model-generated runoff is currently in kg/m2/s. we need it in mm
            # for the whole month
            # so let's start by changing each of these to totals for that day
            # and then we can add each day to get the monthly total
            default_strm = default_strm*86400 # now in kg/m2/day
            default_runoff = default_runoff*86400 # now in kg/m2/day
            default_baseflow = default_baseflow*86400 # now in kg/m2/day
            pso_init_strm = pso_init_strm*86400 # now in kg/m2/day
            pso_final_strm = pso_final_strm*86400 # now in kg/m2/day
            if get_metf:
                default_rainfsnowf = default_rainfsnowf*86400 # now in kg/m2/day
            ## get the weighted average over the whole watershed and add to
            ## array
            ## for default streamflow
            default_strm_np = np.array(default_strm)
            default_strm_avg = np.average(
                default_strm_np,axis=1,weights=this_perc
            )
            default_all_strm[:,w] = default_strm_avg
            # for default runoff
            default_runoff_np = np.array(default_runoff)
            default_runoff_avg = np.average(
                default_runoff_np,axis=1,weights=this_perc
            )
            default_all_runoff[:,w] = default_runoff_avg
            # for default baseflow
            default_baseflow_np = np.array(default_baseflow)
            default_baseflow_avg = np.average(
                default_baseflow_np,axis=1,weights=this_perc
            )
            default_all_baseflow[:,w] = default_baseflow_avg
            # for default metf
            if get_metf:
                default_rainfsnowf_np = np.array(default_rainfsnowf)
                default_rainfsnowf_avg = np.average(
                    default_rainfsnowf_np,axis=1,weights=this_perc
                )
                default_all_rainfsnowf[:,w] = default_rainfsnowf_avg
            # for default le
            default_le_np = np.array(default_le)
            default_le_avg = np.average(
                default_le_np,axis=1,weights=this_perc
            )
            default_all_le[:,w] = default_le_avg
            # for fluxcom le
            fluxcom_le_np = np.array(fluxcom_le)
            fluxcom_le_avg = np.average(
                fluxcom_le_np,axis=1,weights=this_perc
            )
            fluxcom_all_le[:,w] = fluxcom_le_avg
            # for pso_init streamflow
            pso_init_strm_np = np.array(pso_init_strm)
            pso_init_strm_avg = np.average(
                pso_init_strm_np,axis=1,weights=this_perc
            )
            pso_init_all_strm[:,w] = pso_init_strm_avg
            # for pso_init le
            pso_init_le_np = np.array(pso_init_le)
            pso_init_le_avg = np.average(
                pso_init_le_np,axis=1,weights=this_perc
            )
            pso_init_all_le[:,w] = pso_init_le_avg
            # for pso_final streamflow
            pso_final_strm_np = np.array(pso_final_strm)
            pso_final_strm_avg = np.average(
                pso_final_strm_np,axis=1,weights=this_perc
            )
            pso_final_all_strm[:,w] = pso_final_strm_avg
            # for pso_final le
            pso_final_le_np = np.array(pso_final_le)
            pso_final_le_avg = np.average(
                pso_final_le_np,axis=1,weights=this_perc
            )
            pso_final_all_le[:,w] = pso_final_le_avg
            # for observed precip
        times = list(pso_final_runoff.index)
        times_err = list(fluxcom_le.index)
        # for default stream
        default_strm_df = pd.DataFrame(
            default_all_strm,index=times,columns=watersheds
        )
        default_strm_df.index.name = 'date'
        # for default runoff
        default_runoff_df = pd.DataFrame(
            default_all_runoff,index=times,columns=watersheds
        )
        default_runoff_df.index.name = 'date'
        # for default baseflow
        default_baseflow_df = pd.DataFrame(
            default_all_baseflow,index=times,columns=watersheds
        )
        default_baseflow_df.index.name = 'date'
        # for default le
        default_le_df = pd.DataFrame(
            default_all_le,index=times,columns=watersheds
        )
        default_le_df.index.name = 'date'
        # met rainfsnowf
        if get_metf:
            default_rainfsnowf_df = pd.DataFrame(
                default_all_rainfsnowf,index=times,columns=watersheds
            )
            default_rainfsnowf_df.index.name = 'date'
        # for pso_init stream
        pso_init_strm_df = pd.DataFrame(
            pso_init_all_strm,index=times,columns=watersheds
        )
        pso_init_strm_df.index.name = 'date'
        # for pso_init le
        pso_init_le_df = pd.DataFrame(
            pso_init_all_le,index=times,columns=watersheds
        )
        pso_init_le_df.index.name = 'date'
        # for pso_final stream
        pso_final_strm_df = pd.DataFrame(
            pso_final_all_strm,index=times,columns=watersheds
        )
        pso_final_strm_df.index.name = 'date'
        # for pso_final le
        pso_final_le_df = pd.DataFrame(
            pso_final_all_le,index=times,columns=watersheds
        )
        pso_final_le_df.index.name = 'date'
        # for fluxcom le
        fluxcom_le_df = pd.DataFrame(
            fluxcom_all_le,index=times_err,columns=watersheds
        )
        fluxcom_le_df.index.name = 'date'
        outs = {
            'default_strm':default_strm_df,
            'default_runoff':default_runoff_df,
            'default_baseflow':default_baseflow_df,
            'default_le':default_le_df,
            'fluxcom_le':fluxcom_le_df,
            'waterwatch_strm':streamflow_timeseries,
            'pso_init_strm':pso_init_strm_df,
            'pso_init_le':pso_init_le_df,
            'pso_final_strm':pso_final_strm_df,
            'pso_final_le':pso_final_le_df,
        }
        if get_metf:
            outs['default_rainfsnowf'] = default_rainfsnowf_df
        return outs
    def plot_streamflow_timeseries(self,outs,plots_dir,exp_names):
        default_le = outs['default_le']
        pso_init_le = outs['pso_init_strm']
        pso_final_le = outs['pso_final_le']
        fluxcom_le = outs['fluxcom_le']
        watersheds = list(default_le.columns)
        for w,wat in enumerate(watersheds):
            print('plotting le watershed-scale timeseries for wat {}'.format(
                wat
            ))
            this_default_le = default_le[wat]
            this_fluxcom_le = fluxcom_le[wat]
            this_pso_init_le = pso_init_le[wat]
            this_pso_final_le = pso_final_le[wat]
            dates_dtm = list(default_le.index)
            this_exp = exp_names[0]
            savename = os.path.join(
                plots_dir,'le_timeseries_watershed_scale_{}_{}.png'.format(this_exp,wat)
            )
            plt.figure()
            plt.plot(
                dates_dtm,this_default_le,label='Default Catchment-CN4.5',
                c='r',linewidth=.8
            )
            plt.plot(
                dates_dtm,this_fluxcom_le,label='GLEAM',
                c='k',linewidth=.8
            )
            #plt.plot(
            #    dates_dtm,this_pso_init_le,label='CN4.5 Ball-Berry',
            #    c='y',linewidth=.8
            #)
            #plt.plot(
            #    dates_dtm,this_pso_final_le,label='CN4.5 Medlyn',
            #    c='g',linewidth=.8
            #)
            plt.legend()
            plt.xlabel('Date')
            plt.ylabel('LE (W/m2)')
            plt.savefig(savename,dpi=300,bbox_inches='tight')
            plt.close()
    def get_rmse_dict(self,outs,lai_fname,intersection_info,all_tiles,
                      start_err,end_err,plot_timeseries,exp_names,
                      plots_dir,plot_watershed_error_hist,get_metf,
                      out_dir,out_err_fname):
        # let's start by analyzing stream data
        default_strm = outs['default_strm']
        waterwatch_strm = outs['waterwatch_strm']
        pso_init_strm = outs['pso_init_strm']
        pso_final_strm = outs['pso_final_strm']
        default_le = outs['default_le']
        pso_le = outs['pso_final_le']
        fluxcom_le = outs['fluxcom_le']
        if get_metf:
            obs_precip = outs['default_rainfsnowf']
        # let's trim these to only our analysis days
        start_fmt = start_err.strftime('%Y%m%d')
        end_fmt = end_err.strftime('%Y%m%d')
        default_strm = default_strm.loc[start_fmt:end_fmt]
        waterwatch_strm = waterwatch_strm.loc[start_fmt:end_fmt]
        pso_init_strm = pso_init_strm.loc[start_fmt:end_fmt]
        pso_final_strm = pso_final_strm.loc[start_fmt:end_fmt]
        # let's start by returning avg streamflow, streamflow rmse, and lai 
        # values

        # get the lits of watersheds as strings to be used as the dict keys for
        # intersection_info
        watersheds_str = list(intersection_info.keys())
        # get the list of watersheds as integers
        watersheds = [int(k) for k in intersection_info.keys()]
        # set up the df that we are going to add everything to
        cols = copy.deepcopy(watersheds)
        cols.append('all')
        default_df = pd.DataFrame(columns=cols)
        pso_init_df = pd.DataFrame(columns=cols)
        pso_final_df = pd.DataFrame(columns=cols)
        waterwatch_df = pd.DataFrame(columns=cols)
        # lets start by getting the avg lai for each watershed
        #lai_ds = nc.Dataset(lai_fname)
        #lai_vals = np.array(lai_ds['lai'])
        #avg_lai = np.zeros(len(cols))
        for w,wat in enumerate(watersheds):
            # get the tiles in this watershed
            this_tiles = intersection_info[watersheds_str[w]][0]
            # get the percent of each of these tiles in this watershed
            this_perc = intersection_info[watersheds_str[w]][1]
            tiles_idx = np.zeros(0,dtype=np.int64)
            for t,ti in enumerate(this_tiles):
                this_ti_idx = int(np.where(all_tiles == ti)[0][0])
                tiles_idx = np.append(tiles_idx,this_ti_idx)
            # get the average lai for this watershed
            #this_lai = lai_vals[tiles_idx]
            #this_avg_lai = np.average(this_lai,weights=this_perc)
            #avg_lai[w] = this_avg_lai
        #avg_lai_all = np.mean(avg_lai)
        #avg_lai[-1] = avg_lai_all
        #default_df.loc['lai'] = avg_lai
        # let's get average streamflow
        # for the default senario
        default_strm_avg = np.array(default_strm.mean())
        default_strm_avg = np.append(
            default_strm_avg,default_strm_avg.mean())
        default_df.loc['strm'] = default_strm_avg
        # for the pso_init senario
        pso_init_strm_avg = np.array(pso_init_strm.mean())
        pso_init_strm_avg = np.append(
            pso_init_strm_avg,pso_init_strm_avg.mean())
        pso_init_df.loc['strm'] = pso_init_strm_avg
        # for the pso_final senario
        pso_final_strm_avg = np.array(pso_final_strm.mean())
        pso_final_strm_avg = np.append(
            pso_final_strm_avg,pso_final_strm_avg.mean())
        pso_final_df.loc['strm'] = pso_final_strm_avg
        # finally for waterwatch
        waterwatch_strm_avg = np.array(waterwatch_strm.mean())
        waterwatch_strm_avg = np.append(
            waterwatch_strm_avg,waterwatch_strm_avg.mean())
        waterwatch_df.loc['strm'] = waterwatch_strm_avg
        default_df.loc['strm_obs'] = waterwatch_strm_avg
        pso_final_df.loc['strm_obs'] = waterwatch_strm_avg
        # also for default et
        default_le_avg = np.array(default_le.mean())
        default_le_avg = np.append(
            default_le_avg,default_le_avg.mean())
        default_df.loc['le'] = default_le_avg
        # and for pso et
        pso_le_avg = np.array(pso_le.mean())
        pso_le_avg = np.append(
            pso_le_avg,pso_le_avg.mean())
        pso_final_df.loc['le'] = pso_le_avg
        # and for observed et
        fluxcom_le_avg = np.array(fluxcom_le.mean())
        fluxcom_le_avg = np.append(
            fluxcom_le_avg,fluxcom_le_avg.mean())
        default_df.loc['le_obs'] = fluxcom_le_avg
        pso_final_df.loc['le_obs'] = fluxcom_le_avg
        # and for observed precipitation
        if get_metf:
            obs_precip_avg = np.array(obs_precip.mean())
            obs_precip_avg = np.append(
                obs_precip_avg,obs_precip_avg.mean()
            )
            default_df.loc['precip_obs'] = obs_precip_avg
        # let's get the rmse for the different senarios
        catch_cols = list(default_strm.columns)
        waterwatch_strm_int = waterwatch_strm.set_axis(catch_cols,axis=1)
        # okay, we need to average to yearly now
        year_start = copy.deepcopy(start_err)
        year_end = (
            year_start +
            relativedelta.relativedelta(years=1) -
            relativedelta.relativedelta(days=1)
        )
        model_times = np.array(pso_final_strm.index)
        cols_yr = cols[:-1]
        cols_yr = np.append(cols_yr,'time')
        default_strm_yr = pd.DataFrame(columns=cols_yr)
        default_strm_yr = default_strm_yr.set_index('time')
        default_strm_yr = default_strm_yr.set_axis(catch_cols,axis=1)
        pso_init_strm_yr = pd.DataFrame(columns=cols_yr)
        pso_init_strm_yr = pso_init_strm_yr.set_index('time')
        pso_init_strm_yr = pso_init_strm_yr.set_axis(catch_cols,axis=1)
        pso_final_strm_yr = pd.DataFrame(columns=cols_yr)
        pso_final_strm_yr = pso_final_strm_yr.set_index('time')
        pso_final_strm_yr = pso_final_strm_yr.set_axis(catch_cols,axis=1)
        while year_end <= end_err:
            start_idx = np.where(
                model_times == np.datetime64(year_start)
            )[0][0]
            end_idx = np.where(
                model_times == np.datetime64(year_end)
            )[0][0]
            default_this_year_vals = default_strm.iloc[
                start_idx:end_idx+1
            ]
            pso_init_this_year_vals = pso_init_strm.iloc[
                start_idx:end_idx+1
            ]
            pso_final_this_year_vals = pso_final_strm.iloc[
                start_idx:end_idx+1
            ]
            avg_default_strm = np.array(default_this_year_vals.mean())
            avg_pso_init_strm = np.array(pso_init_this_year_vals.mean())
            avg_pso_final_strm = np.array(pso_final_this_year_vals.mean())
            idx_name = year_start.strftime('%Y%m%d')
            idx_name = np.int64(idx_name)
            default_strm_yr.loc[idx_name] = avg_default_strm
            pso_init_strm_yr.loc[idx_name] = avg_pso_init_strm
            pso_final_strm_yr.loc[idx_name] = avg_pso_final_strm
            year_start += relativedelta.relativedelta(years=1)
            year_end += relativedelta.relativedelta(years=1)
        dates = list(default_strm_yr.index)
        dates_dtm = [
            datetime.datetime.strptime(str(d),'%Y%m%d') for d in dates
        ]
        if plot_timeseries:
            for w,wat in enumerate(watersheds):
                print(
                    'plotting timeseries for watershed {}'.format(
                        wat
                    )
                )
                this_default_strm = default_strm_yr[wat]
                this_pso_init_strm = pso_init_strm_yr[wat]
                this_pso_final_strm = pso_final_strm_yr[wat]
                this_waterwatch_strm = waterwatch_strm_int[wat]
                this_exp = exp_names[2]
                savename = os.path.join(
                    plots_dir,'streamflow_yearly_timeseries_{}_{}.png'.format(
                        this_exp,wat
                    )
                )
                plt.figure()
                plt.plot(
                    dates_dtm,this_default_strm,label='CN4.5 EF It. 1',
                    c='r',marker='o'#linewidth=.8
                )
                plt.plot(
                    dates_dtm,this_waterwatch_strm,label='CAMELS',
                    c='k',marker='o'#linewidth=.8
                )
                plt.plot(
                    dates_dtm,this_pso_init_strm,label='CN4.5 EF It. 4',
                    c='y',marker='o'#,linewidth=.8
                )
                plt.plot(
                    dates_dtm,this_pso_final_strm,label='CN4.5 EF It. 8',
                    c='g',marker='o'#linewidth=.8
                )
                plt.legend()
                plt.xlabel('Date')
                plt.ylabel('mm/day')
                plt.savefig(savename,dpi=300,bbox_inches='tight')
                plt.close()
        # caluclate rmse for default catchment-cn
        #default_diff = np.array(default_strm_yr - waterwatch_strm_int)
        #default_diff_ave = np.mean(default_diff,axis=0)
        #default_diff_ave_avg = np.mean(default_diff_ave)
        #default_diff_ave = np.append(default_diff_ave,default_diff_ave_avg)
        #default_df.loc['strm_diff'] = default_diff_ave
        #default_se = default_diff**2
        #default_mse = np.mean(default_se,axis=0)
        #default_rmse = np.sqrt(default_mse)
        #default_rmse_all = np.mean(default_rmse)
        #default_rmse = np.append(default_rmse,default_rmse_all)
        #default_df.loc['strm_rmse'] = default_rmse
        # calculate r2 for default catchment-CN
        default_diff = np.zeros(len(catch_cols))
        pso_diff = np.zeros(len(catch_cols))
        default_rmse = np.zeros(len(catch_cols))
        pso_rmse = np.zeros(len(catch_cols))
        default_rmse_norm_ind = np.zeros(len(catch_cols))
        pso_rmse_norm_ind = np.zeros(len(catch_cols))
        default_r2 = np.zeros(len(catch_cols))
        pso_r2 = np.zeros(len(catch_cols))
        default_corr = np.zeros(len(catch_cols))
        pso_corr = np.zeros(len(catch_cols))
        default_ubrmse = np.zeros(len(catch_cols))
        pso_ubrmse = np.zeros(len(catch_cols))
        default_nse = np.zeros(len(catch_cols))
        pso_nse = np.zeros(len(catch_cols))
        default_nse_avg = np.zeros(len(catch_cols))
        pso_nse_avg = np.zeros(len(catch_cols))
        default_mae = np.zeros(len(catch_cols))
        pso_mae = np.zeros(len(catch_cols))
        pso_mae_norm_def = np.zeros(len(catch_cols))
        pso_mae_norm_def = np.zeros(len(catch_cols))
        default_diff_le = np.zeros(len(catch_cols))
        pso_diff_le = np.zeros(len(catch_cols))
        default_rmse_le = np.zeros(len(catch_cols))
        pso_rmse_le = np.zeros(len(catch_cols))
        default_r2_le = np.zeros(len(catch_cols))
        pso_r2_le = np.zeros(len(catch_cols))
        default_corr_le = np.zeros(len(catch_cols))
        pso_corr_le = np.zeros(len(catch_cols))
        default_ubrmse_le = np.zeros(len(catch_cols))
        pso_ubrmse_le = np.zeros(len(catch_cols))
        for c,col in enumerate(catch_cols):
            this_default = np.array(default_strm_yr[col])
            this_pso = np.array(pso_final_strm_yr[col])
            this_waterwatch = np.array(waterwatch_strm_int[col])
            not_nan_idx = np.where(
                np.isnan(this_waterwatch) == False
            )
            this_default_le = np.array(default_le[col])
            this_pso_le = np.array(pso_le[col])
            this_fluxcom = np.array(fluxcom_le[col])
            not_nan_idx_le = np.where(
                np.isnan(this_fluxcom) == False
            )
            this_default = this_default[not_nan_idx]
            this_pso = this_pso[not_nan_idx]
            this_waterwatch = this_waterwatch[not_nan_idx]
            this_default_le = this_default_le[not_nan_idx_le]
            this_pso_le = this_pso_le[not_nan_idx_le]
            this_fluxcom = this_fluxcom[not_nan_idx_le]
            # calculate the default difference in average
            element_default_diff = this_default - this_waterwatch
            this_default_diff = np.mean(element_default_diff)
            # calculate the pso different in average
            element_pso_diff = this_pso - this_waterwatch
            this_pso_diff = np.mean(element_pso_diff)
            # calculate the default rmse
            default_se = element_default_diff**2
            default_mse = np.mean(default_se)
            this_default_rmse = np.sqrt(default_mse)
            this_default_rmse_norm_ind = (
                this_default_rmse/waterwatch_strm_avg[c]
            )
            # calculate the pso rmse
            pso_se = element_pso_diff**2
            pso_mse = np.mean(pso_se)
            this_pso_rmse = np.sqrt(pso_mse)
            this_pso_rmse_norm_ind = (
                this_pso_rmse/waterwatch_strm_avg[c]
            )
            # calculate the r2 for default
            this_rss = np.sum(
                np.square(
                    this_waterwatch - this_default
                )
            )
            this_avg_y = np.average(this_waterwatch)
            this_tss = np.sum(
                np.square(
                    this_waterwatch - this_avg_y
                )
            )
            this_default_r2 = 1 - this_rss/this_tss
            # calculate the r2 for pso
            this_rss = np.sum(
                np.square(
                    this_waterwatch - this_pso
                )
            )
            this_avg_y = np.average(this_waterwatch)
            this_tss = np.sum(
                np.square(
                    this_waterwatch - this_avg_y
                )
            )
            this_pso_r2 = 1 - this_rss/this_tss
            # calcaulte the default correlation coefficient
            this_avg_x = np.mean(this_default)
            numerator = np.sum(
                (this_default - this_avg_x)*
                (this_waterwatch - this_avg_y)
            )
            denominator_1 = np.sqrt(
                np.sum(
                    np.square(
                        this_default - this_avg_x
                    )
                )
            )
            denominator_2 = np.sqrt(
                np.sum(
                    np.square(
                        this_waterwatch - this_avg_y
                    )
                )
            )
            denominator = denominator_1*denominator_2
            this_default_corr = numerator/denominator
            # calcaulte the pso correlation coefficient
            this_avg_x = np.mean(this_pso)
            numerator = np.sum(
                (this_pso - this_avg_x)*
                (this_waterwatch - this_avg_y)
            )
            denominator_1 = np.sqrt(
                np.sum(
                    np.square(
                        this_pso - this_avg_x
                    )
                )
            )
            denominator_2 = np.sqrt(
                np.sum(
                    np.square(
                        this_waterwatch - this_avg_y
                    )
                )
            )
            denominator = denominator_1*denominator_2
            this_pso_corr = numerator/denominator
            # calculate default ubrmse
            this_default_ubrmse = np.sqrt(
                (
                    (
                        (
                            this_default -
                            this_avg_x
                        ) - (
                            this_waterwatch -
                            this_avg_y
                        )
                    )**2
                ).mean()
            )
            # calculate pso ubrmse
            this_pso_ubrmse = np.sqrt(
                (
                    (
                        (
                            this_pso -
                            this_avg_x
                        ) - (
                            this_waterwatch -
                            this_avg_y
                        )
                    )**2
                ).mean()
            )
            # calculate the default nse
            nse_numerator = np.sum(
                np.square(
                    this_waterwatch - this_default
                )
            )
            nse_denom = np.sum(
                np.square(
                    this_waterwatch - this_avg_y
                )
            )
            this_default_nse = 1 - nse_numerator/nse_denom
            # calculate the pso nse
            nse_numerator = np.sum(
                np.square(
                    this_waterwatch - this_pso
                )
            )
            nse_denom = np.sum(
                np.square(
                    this_waterwatch - this_avg_y
                )
            )
            this_pso_nse = 1 - nse_numerator/nse_denom
            # calculate the default mae
            this_default_mae = np.mean(
                np.abs(
                    this_waterwatch - this_default
                )
            )
            # calculate the pso mae
            this_pso_mae = np.mean(
                np.abs(
                    this_waterwatch - this_pso
                )
            )
            # let's just look at the error distributions
            this_wat_err = np.abs(
                this_waterwatch - this_default
            )
            # calculate hte pso mae normalized by default error
            this_pso_mae_norm_def = np.mean(
                np.abs(
                    this_waterwatch - this_pso
                )
            )
            this_pso_mae_norm_def = (
                this_pso_mae_norm_def/this_default_mae
            )
            #this_wat_err = this_wat_err - np.mean(this_wat_err)
            #this_wat_err = this_wat_err/np.std(this_wat_err)
            if plot_watershed_error_hist:
                save_name = os.path.join(
                    plots_dir,
                    'wat_error_hist_{}.png'.format(col)
                )
                bin_width = 0.1
                plt.figure()
                plt.hist(
                    this_wat_err,
                    bins=7
                    #bins=np.arange(
                    #    min(this_wat_err),
                    #    max(this_wat_err) + bin_width,
                    #    bin_width
                    #)
                )
                plt.savefig(save_name)
                plt.close()
                this_wat_err_log = np.abs(
                    this_waterwatch - this_default
                )
                this_wat_err_log = np.log(this_wat_err_log)
                #this_wat_err_log = this_wat_err_log - np.mean(
                #    this_wat_err_log
                #)
                #this_wat_err_log = this_wat_err_log/np.std(
                #    this_wat_err_log
                #)
                save_name = os.path.join(
                    plots_dir,
                    'wat_error_hist_{}_log.png'.format(col)
                )
                bin_width_log = 0.05
                plt.figure()
                plt.hist(
                    this_wat_err_log,
                    bins=7
                    #bins=np.arange(
                    #    min(this_wat_err_log),
                    #    max(this_wat_err_log) + bin_width_log,
                    #    bin_width_log
                    #)
                )
                plt.savefig(save_name)
                plt.close()
            # add to running error metrics
            default_diff[c] = this_default_diff
            pso_diff[c] = this_pso_diff
            default_rmse[c] = this_default_rmse
            pso_rmse[c] = this_pso_rmse
            default_r2[c] = this_default_r2
            default_rmse_norm_ind[c] = this_default_rmse_norm_ind
            pso_rmse_norm_ind[c] = this_pso_rmse_norm_ind
            pso_r2[c] = this_pso_r2
            default_corr[c] = this_default_corr
            pso_corr[c] = this_pso_corr
            default_ubrmse[c] = this_default_ubrmse
            pso_ubrmse[c] = this_pso_ubrmse
            default_nse[c] = this_default_nse
            pso_nse[c] = this_pso_nse
            default_mae[c] = this_default_mae
            pso_mae[c] = this_pso_mae
            pso_mae_norm_def[c] = this_pso_mae_norm_def
            # calculate the default difference in average
            element_default_diff_le = this_default_le - this_fluxcom
            this_default_diff_le = np.mean(element_default_diff_le)
            # calculate the pso different in average
            element_pso_diff_le = this_pso_le - this_fluxcom
            this_pso_diff_le = np.mean(element_pso_diff_le)
            # calculate the default rmse
            default_se_le = element_default_diff_le**2
            default_mse_le = np.mean(default_se_le)
            this_default_rmse_le = np.sqrt(default_mse_le)
            # calculate the pso rmse
            pso_se_le = element_pso_diff_le**2
            pso_mse_le = np.mean(pso_se_le)
            this_pso_rmse_le = np.sqrt(pso_mse_le)
            # calculate the r2 for default
            this_rss_le = np.sum(
                np.square(
                    this_fluxcom - this_default_le
                )
            )
            this_avg_y_le = np.average(this_fluxcom)
            this_tss_le = np.sum(
                np.square(
                    this_fluxcom - this_avg_y_le
                )
            )
            this_default_r2_le = 1 - this_rss_le/this_tss_le
            # calculate the r2 for pso
            this_rss_le = np.sum(
                np.square(
                    this_fluxcom - this_pso_le
                )
            )
            this_avg_y_le = np.average(this_fluxcom)
            this_tss_le = np.sum(
                np.square(
                    this_fluxcom - this_avg_y_le
                )
            )
            this_pso_r2_le = 1 - this_rss_le/this_tss_le
            # calcaulte the default correlation coefficient
            this_avg_x_le = np.mean(this_default_le)
            numerator_le = np.sum(
                (this_default_le - this_avg_x_le)*
                (this_fluxcom - this_avg_y_le)
            )
            denominator_1_le = np.sqrt(
                np.sum(
                    np.square(
                        this_default_le - this_avg_x_le
                    )
                )
            )
            denominator_2_le = np.sqrt(
                np.sum(
                    np.square(
                        this_fluxcom - this_avg_y_le
                    )
                )
            )
            denominator_le = denominator_1_le*denominator_2_le
            this_default_corr_le = numerator_le/denominator_le
            # calcaulte the pso correlation coefficient
            this_avg_x_le = np.mean(this_pso_le)
            numerator_le = np.sum(
                (this_pso_le - this_avg_x_le)*
                (this_fluxcom - this_avg_y_le)
            )
            denominator_1_le = np.sqrt(
                np.sum(
                    np.square(
                        this_pso_le - this_avg_x_le
                    )
                )
            )
            denominator_2_le = np.sqrt(
                np.sum(
                    np.square(
                        this_fluxcom - this_avg_y_le
                    )
                )
            )
            denominator_le = denominator_1_le*denominator_2_le
            this_pso_corr_le = numerator_le/denominator_le
            # calculate default ubrmse
            this_default_ubrmse_le = np.sqrt(
                (
                    (
                        (
                            this_default_le -
                            this_avg_x_le
                        ) - (
                            this_fluxcom -
                            this_avg_y_le
                        )
                    )**2
                ).mean()
            )
            # calculate pso ubrmse
            this_pso_ubrmse_le = np.sqrt(
                (
                    (
                        (
                            this_pso_le -
                            this_avg_x_le
                        ) - (
                            this_fluxcom -
                            this_avg_y_le
                        )
                    )**2
                ).mean()
            )
            # add to running error metrics
            default_diff_le[c] = this_default_diff
            pso_diff_le[c] = this_pso_diff_le
            default_rmse_le[c] = this_default_rmse_le
            pso_rmse_le[c] = this_pso_rmse_le
            default_r2_le[c] = this_default_r2_le
            pso_r2_le[c] = this_pso_r2_le
            default_corr_le[c] = this_default_corr_le
            pso_corr_le[c] = this_pso_corr_le
            default_ubrmse_le[c] = this_default_ubrmse_le
            pso_ubrmse_le[c] = this_pso_ubrmse_le
        # we're going to need some np arrays for this
        waterwatch_strm_int_np = np.array(waterwatch_strm_int)
        default_strm_yr_np = np.array(default_strm_yr)
        pso_final_strm_yr_np = np.array(pso_final_strm_yr)
        all_not_nan_idx = np.where(
            np.isnan(waterwatch_strm_int_np) == False
        )
        waterwatch_strm_int_np = waterwatch_strm_int_np[
            all_not_nan_idx
        ]
        default_strm_yr_np = default_strm_yr_np[all_not_nan_idx]
        pso_final_strm_yr_np = pso_final_strm_yr_np[all_not_nan_idx]
        fluxcom_le_np = np.array(fluxcom_le)
        default_le_np = np.array(default_le)
        pso_le_np = np.array(pso_le)
        all_not_nan_idx_le = np.where(
            np.isnan(fluxcom_le_np) == False
        )
        fluxcom_le_np = fluxcom_le_np[
            all_not_nan_idx_le
        ]
        default_le_np = default_le_np[
            all_not_nan_idx_le
        ]
        pso_le_np = pso_le_np[
            all_not_nan_idx_le
        ]
        # calculate the diff for default all
        all_avg_waterwatch = np.average(waterwatch_strm_int_np)
        all_avg_default = np.average(default_strm_yr_np)
        all_default_diff = all_avg_default - all_avg_waterwatch
        default_diff = np.append(default_diff,all_default_diff)
        # calculate the diff for pso all
        all_avg_pso = np.average(pso_final_strm_yr_np)
        all_pso_diff = all_avg_pso - all_avg_waterwatch
        pso_diff = np.append(pso_diff,all_pso_diff)
        # calculate the default rmse for all
        all_default_diff = default_strm_yr_np - waterwatch_strm_int_np
        all_default_se = all_default_diff**2
        all_default_mse = np.mean(all_default_se)
        all_default_rmse = np.sqrt(all_default_mse)
        default_rmse = np.append(default_rmse,all_default_rmse)
        default_rmse_norm = default_rmse/waterwatch_strm_avg[-1]
        all_default_rmse_norm_ind = np.average(default_rmse_norm_ind)
        default_rmse_norm_ind = np.append(
            default_rmse_norm_ind,all_default_rmse_norm_ind
        )
        # calculate the pso rmse for all
        all_pso_diff = pso_final_strm_yr_np - waterwatch_strm_int_np
        all_pso_se = all_pso_diff**2
        all_pso_mse = np.mean(all_pso_se)
        all_pso_rmse = np.sqrt(all_pso_mse)
        pso_rmse = np.append(pso_rmse,all_pso_rmse)
        pso_rmse_norm = pso_rmse/waterwatch_strm_avg[-1]
        all_pso_rmse_norm_ind = np.average(pso_rmse_norm_ind)
        pso_rmse_norm_ind = np.append(
            pso_rmse_norm_ind,all_pso_rmse_norm_ind
        )
        # calculate the r2 for default all
        all_rss = np.sum(
            np.square(
                np.array(waterwatch_strm_int_np - default_strm_yr_np)
            )
        )
        all_avg_y = np.average(waterwatch_strm_int_np)
        all_tss = np.sum(
            np.square(
                np.array(waterwatch_strm_int_np - all_avg_y)
            )
        )
        default_all_r2 = 1 - all_rss/all_tss
        default_r2 = np.append(default_r2,default_all_r2)
        # calculate the r2 for pso all
        all_rss = np.sum(
            np.square(
                np.array(waterwatch_strm_int_np - pso_final_strm_yr_np)
            )
        )
        all_avg_y = np.average(waterwatch_strm_int_np)
        all_tss = np.sum(
            np.square(
                np.array(waterwatch_strm_int_np - all_avg_y)
            )
        )
        pso_all_r2 = 1 - all_rss/all_tss
        pso_r2 = np.append(pso_r2,pso_all_r2)
        # calcaulte the default correlation coefficient
        all_avg_x = np.mean(default_strm_yr_np)
        numerator = np.sum(
            (default_strm_yr_np - all_avg_x)*
            (waterwatch_strm_int_np - all_avg_y)
        )
        denominator_1 = np.sqrt(
            np.sum(
                np.square(
                    default_strm_yr_np - all_avg_x
                )
            )
        )
        denominator_2 = np.sqrt(
            np.sum(
                np.square(
                    waterwatch_strm_int_np - all_avg_y
                )
            )
        )
        denominator = denominator_1*denominator_2
        all_default_corr = numerator/denominator
        default_corr = np.append(default_corr,all_default_corr)
        # calcaulte the pso correlation coefficient
        all_avg_x = np.mean(pso_final_strm_yr_np)
        numerator = np.sum(
            (pso_final_strm_yr_np - all_avg_x)*
            (waterwatch_strm_int_np - all_avg_y)
        )
        denominator_1 = np.sqrt(
            np.sum(
                np.square(
                    pso_final_strm_yr_np - all_avg_x
                )
            )
        )
        denominator_2 = np.sqrt(
            np.sum(
                np.square(
                    waterwatch_strm_int_np - all_avg_y
                )
            )
        )
        denominator = denominator_1*denominator_2
        all_pso_corr = numerator/denominator
        pso_corr = np.append(pso_corr,all_pso_corr)
        # calculate default ubrmse
        all_default_ubrmse = np.sqrt(
            (
                (
                    (
                        default_strm_yr_np -
                        all_avg_x
                    ) - (
                        waterwatch_strm_int_np -
                        all_avg_y
                    )
                )**2
            ).mean()
        )
        default_ubrmse = np.append(default_ubrmse,all_default_ubrmse)
        # calculate pso ubrmse
        all_pso_ubrmse = np.sqrt(
            (
                (
                    (
                        pso_final_strm_yr_np -
                        all_avg_x
                    ) - (
                        waterwatch_strm_int_np -
                        all_avg_y
                    )
                )**2
            ).mean()
        )
        pso_ubrmse = np.append(pso_ubrmse,all_pso_ubrmse)
        # calculate the default nse
        nse_numerator = np.sum(
            np.square(
                waterwatch_strm_int_np - default_strm_yr_np
            )
        )
        nse_denom = np.sum(
            np.square(
                waterwatch_strm_int_np - all_avg_y
            )
        )
        all_default_nse = 1 - nse_numerator/nse_denom
        all_default_nse_ind = np.average(default_nse)
        default_nse = np.append(default_nse,all_default_nse)
        default_nse_ind = copy.deepcopy(default_nse)
        default_nse_ind[-1] = all_default_nse_ind
        # calculate the pso nse
        nse_numerator = np.sum(
            np.square(
                waterwatch_strm_int_np - pso_final_strm_yr_np
            )
        )
        nse_denom = np.sum(
            np.square(
                waterwatch_strm_int_np - all_avg_y
            )
        )
        all_pso_nse = 1 - nse_numerator/nse_denom
        all_pso_nse_ind = np.average(pso_nse)
        pso_nse = np.append(pso_nse,all_pso_nse)
        pso_nse_ind = copy.deepcopy(pso_nse)
        pso_nse_ind[-1] = all_pso_nse_ind
        # calculate the mae for default all
        all_default_mae = np.mean(default_mae)
        default_mae = np.append(default_mae,all_default_mae)
        if out_err_fname != 'none':
            to_save = default_mae[:-1]
            to_save.tofile(
                os.path.join(
                    out_dir,
                    out_err_fname
                ),
                sep = ','
            )
        # calculate the mae for pso all
        all_pso_mae = np.mean(pso_mae)
        pso_mae = np.append(pso_mae,all_pso_mae)
        # calculate the default mae normalized mae
        all_pso_mae_norm_def = np.mean(pso_mae_norm_def)
        pso_mae_norm_def = np.append(pso_mae_norm_def,all_pso_mae_norm_def)
        # calculate the diff for default all
        all_avg_fluxcom_le = np.average(fluxcom_le_np)
        all_avg_default_le = np.average(default_le_np)
        all_default_diff_le = all_avg_default_le - all_avg_fluxcom_le
        default_diff_le = np.append(default_diff_le,all_default_diff_le)
        # calculate the diff for pso all
        all_avg_pso_le = np.average(pso_le_np)
        all_pso_diff_le = all_avg_pso_le - all_avg_fluxcom_le
        pso_diff_le = np.append(pso_diff_le,all_pso_diff_le)
        # calculate the default rmse for all
        all_default_diff_le = default_le_np - fluxcom_le_np
        all_default_se_le = all_default_diff_le**2
        all_default_mse_le = np.mean(all_default_se_le)
        all_default_rmse_le = np.sqrt(all_default_mse_le)
        default_rmse_le = np.append(default_rmse_le,all_default_rmse_le)
        # calculate the pso rmse for all
        all_pso_diff_le = pso_le_np - fluxcom_le_np
        all_pso_se_le = all_pso_diff_le**2
        all_pso_mse_le = np.mean(all_pso_se_le)
        all_pso_rmse_le = np.sqrt(all_pso_mse_le)
        pso_rmse_le = np.append(pso_rmse_le,all_pso_rmse_le)
        # calculate the r2 for default all
        all_rss_le = np.sum(
            np.square(
                np.array(fluxcom_le_np - default_le_np)
            )
        )
        all_avg_y_le = np.average(fluxcom_le_np)
        all_tss_le = np.sum(
            np.square(
                np.array(fluxcom_le_np - all_avg_y_le)
            )
        )
        default_all_r2_le = 1 - all_rss_le/all_tss_le
        default_r2_le = np.append(default_r2_le,default_all_r2_le)
        # calculate the r2 for pso all
        all_rss_le = np.sum(
            np.square(
                np.array(fluxcom_le_np - pso_le_np)
            )
        )
        all_avg_y_le = np.average(fluxcom_le_np)
        all_tss_le = np.sum(
            np.square(
                np.array(fluxcom_le_np - all_avg_y_le)
            )
        )
        pso_all_r2_le = 1 - all_rss_le/all_tss_le
        pso_r2_le = np.append(pso_r2_le,pso_all_r2_le)
        # calcaulte the default correlation coefficient
        all_avg_x_le = np.mean(default_le_np)
        numerator_le = np.sum(
            (default_le_np - all_avg_x_le)*
            (fluxcom_le_np - all_avg_y_le)
        )
        denominator_1_le = np.sqrt(
            np.sum(
                np.square(
                    default_le_np - all_avg_x_le
                )
            )
        )
        denominator_2_le = np.sqrt(
            np.sum(
                np.square(
                    fluxcom_le_np - all_avg_y_le
                )
            )
        )
        denominator_le = denominator_1_le*denominator_2_le
        all_default_corr_le = numerator_le/denominator_le
        default_corr_le = np.append(default_corr_le,all_default_corr_le)
        # calcaulte the pso correlation coefficient
        all_avg_x_le = np.mean(pso_le_np)
        numerator_le = np.sum(
            (pso_le_np - all_avg_x_le)*
            (fluxcom_le_np - all_avg_y_le)
        )
        denominator_1_le = np.sqrt(
            np.sum(
                np.square(
                    pso_le_np - all_avg_x_le
                )
            )
        )
        denominator_2_le = np.sqrt(
            np.sum(
                np.square(
                    fluxcom_le_np - all_avg_y_le
                )
            )
        )
        denominator_le = denominator_1_le*denominator_2_le
        all_pso_corr_le = numerator_le/denominator_le
        pso_corr_le = np.append(pso_corr_le,all_pso_corr_le)
        # calculate default ubrmse
        all_default_ubrmse_le = np.sqrt(
            (
                (
                    (
                        default_le_np -
                        all_avg_x_le
                    ) - (
                        fluxcom_le_np -
                        all_avg_y_le
                    )
                )**2
            ).mean()
        )
        default_ubrmse_le = np.append(default_ubrmse_le,all_default_ubrmse_le)
        # calculate pso ubrmse
        all_pso_ubrmse_le = np.sqrt(
            (
                (
                    (
                        pso_le_np -
                        all_avg_x_le
                    ) - (
                        fluxcom_le_np -
                        all_avg_y_le
                    )
                )**2
            ).mean()
        )
        pso_ubrmse_le = np.append(pso_ubrmse_le,all_pso_ubrmse_le)
        # calculate the default normalized rmse
        default_rmse_normal = default_rmse/all_avg_waterwatch
        # calculate the pso normalized rmse
        pso_rmse_normal = pso_rmse/all_avg_waterwatch
        # lets add all these new metrics to the dataframes
        default_df.loc['strm_avg_diff'] = default_diff
        default_df.loc['strm_rmse'] = default_rmse
        default_df.loc['strm_rmse_norm'] = default_rmse_norm
        default_df.loc['strm_rmse_norm_ind'] = default_rmse_norm_ind
        default_df.loc['strm_r2'] = default_r2
        default_df.loc['strm_corr'] = default_corr
        default_df.loc['strm_ubrmse'] = default_ubrmse
        default_df.loc['strm_nse'] = default_nse
        default_df.loc['strm_nse_ind'] = default_nse_ind
        default_df.loc['strm_mae'] = default_mae
        default_df.loc['le_avg_diff'] = default_diff_le
        default_df.loc['le_rmse'] = default_rmse_le
        default_df.loc['le_r2'] = default_r2_le
        default_df.loc['le_corr'] = default_corr_le
        default_df.loc['le_ubrmse'] = default_ubrmse_le
        default_df.loc['strm_rmse_norm'] = default_rmse_normal
        pso_final_df.loc['strm_avg_diff'] = pso_diff
        pso_final_df.loc['strm_rmse'] = pso_rmse
        pso_final_df.loc['strm_rmse_norm'] = pso_rmse_norm
        pso_final_df.loc['strm_rmse_norm_ind'] = pso_rmse_norm_ind
        pso_final_df.loc['strm_r2'] = pso_r2
        pso_final_df.loc['strm_corr'] = pso_corr
        pso_final_df.loc['strm_ubrmse'] = pso_ubrmse
        pso_final_df.loc['strm_nse'] = pso_nse
        pso_final_df.loc['strm_nse_ind'] = pso_nse_ind
        pso_final_df.loc['strm_mae'] = pso_mae
        pso_final_df.loc['strm_mae_norm_def'] = pso_mae_norm_def
        pso_final_df.loc['le_avg_diff'] = pso_diff_le
        pso_final_df.loc['le_rmse'] = pso_rmse_le
        pso_final_df.loc['le_r2'] = pso_r2_le
        pso_final_df.loc['le_corr'] = pso_corr_le
        pso_final_df.loc['le_ubrmse'] = pso_ubrmse_le
        pso_final_df.loc['strm_rmse_norm'] = pso_rmse_normal
        # caluclate rmse for pso_init catchment-cn
        pso_init_diff = np.array(pso_init_strm_yr - waterwatch_strm_int)
        pso_init_se = pso_init_diff**2
        pso_init_mse = np.mean(pso_init_se,axis=0)
        pso_init_rmse = np.sqrt(pso_init_mse)
        pso_init_rmse_all = np.mean(pso_init_rmse)
        pso_init_rmse = np.append(pso_init_rmse,pso_init_rmse_all)
        pso_init_df.loc['strm_rmse'] = pso_init_rmse
        # put these in their respective dataframes
        to_return = [
            default_df,pso_init_df,pso_final_df,waterwatch_df
        ]
        return to_return
    def plot_water_metrics(self,water_dict,plots_dir):
        # let's first make a plot of average lai versus streamflow
        # under-prediction to see if increasing ksat via lai might help with
        # catchment-cn underprediction problem
        # get the information we need
        default_strm = water_dict['default_strm']
        camels_strm = water_dict['waterwatch_strm']
        default_le = water_dict['default_le']
        fluxcom_le = water_dict['fluxcom_le']
        # get the dates information neecessary for dividing into 2006 and 2007
        dates = np.array(default_strm.index)
        watersheds = np.array(default_strm.columns)
        dates_2006_idx = np.where(dates < 20070000)[0]
        dates_2007_idx = np.where(dates >= 20070000)[0]
        # get the averages for each watershed over the desired years for
        # streamflow
        default_strm_2006 = default_strm.iloc[dates_2006_idx]
        default_strm_2006_avg = np.array(default_strm_2006.mean())
        default_strm_2007 = default_strm.iloc[dates_2007_idx]
        default_strm_2007_avg = np.array(default_strm_2007.mean())
        camels_strm_2006 = camels_strm.iloc[dates_2006_idx]
        camels_strm_2006_avg = np.array(camels_strm_2006.mean())
        camels_strm_2007 = camels_strm.iloc[dates_2007_idx]
        camels_strm_2007_avg = np.array(camels_strm_2007.mean())
        default_strm_2006_2007_avg = np.array(default_strm.mean())
        camels_strm_2006_2007_avg = np.array(camels_strm.mean())
        # get the averages over each watershed over the desired years for le
        default_le_2006 = default_le.iloc[dates_2006_idx]
        default_le_2006_avg = np.array(default_le_2006.mean())
        default_le_2007 = default_le.iloc[dates_2007_idx]
        default_le_2007_avg = np.array(default_le_2007.mean())
        default_le_2006_2007_avg = np.array(default_le.mean())
        fluxcom_le_2006 = fluxcom_le.iloc[dates_2006_idx]
        fluxcom_le_2006_avg = np.array(fluxcom_le_2006.mean())
        fluxcom_le_2007 = fluxcom_le.iloc[dates_2007_idx]
        fluxcom_le_2007_avg = np.array(fluxcom_le_2007.mean())
        fluxcom_le_2006_2007_avg = np.array(fluxcom_le.mean())
        # get the difference in avg le/strm between 2006 and 2007
        default_diff_strm_2006_2007 = (
            np.sum(
                np.absolute(
                    (
                        default_strm_2006_avg - default_strm_2007_avg
                    )
                    /default_strm_2006_avg
                )
                *100
            )
            /len(default_strm_2006_avg)
        )
        camels_diff_strm_2006_2007 = (
            np.sum(
                np.absolute(
                    (
                        camels_strm_2006_avg - camels_strm_2007_avg
                    )
                    /camels_strm_2006_avg
                )
                *100
            )
            /len(camels_strm_2006_avg)
        )
        default_diff_le_2006_2007 = (
            default_le_2006_avg - default_le_2007_avg
        )
        fluxcom_diff_le_2006_2007 = (
            np.sum(
                np.absolute(
                    (
                        fluxcom_le_2006_avg - fluxcom_le_2007_avg
                    )
                    /fluxcom_le_2006_avg
                )
                *100
            )
            /len(fluxcom_le_2006_avg)
        )
        # what do we want to save these as?
        savename_strm_2006 = os.path.join(
            plots_dir,'camels_strm_vs_default_strm_2006.png'
        )
        savename_strm_2007 = os.path.join(
            plots_dir,'camels_strm_vs_default_strm_2007.png'
        )
        savename_strm_2006_2007 = os.path.join(
            plots_dir,'camels_strm_vs_default_strm_2006_2007.png'
        )
        savename_le_2006 = os.path.join(
            plots_dir,'fluxcom_le_vs_default_le_2006.png'
        )
        savename_le_2007 = os.path.join(
            plots_dir,'fluxcom_le_vs_default_le_2007.png'
        )
        savename_le_2006_2007 = os.path.join(
            plots_dir,'fluxcom_le_vs_default_le_2006_2007.png'
        )
        savename_default_diff_strm = os.path.join(
            plots_dir,'default_avg_strm_2006_vs_default_avg_strm_2007.png'
        )
        savename_camels_diff_strm = os.path.join(
            plots_dir,'camels_avg_strm_2006_vs_camels_avg_strm_2007.png'
        )
        savename_default_diff_le = os.path.join(
            plots_dir,'default_avg_le_2006_vs_default_avg_le_2007.png'
        )
        savename_fluxcom_diff_le = os.path.join(
            plots_dir,'fluxcom_avg_le_2006_vs_fluxcom_avg_le_2007.png'
        )
        # make the plot for 2006 streamflow
        plt.figure()
        plt.scatter(camels_strm_2006_avg,default_strm_2006_avg,s=2)
        #plt.xlim([0,40])
        #plt.ylim([0,40])
        max_val_x = np.max(camels_strm_2006_avg)
        max_val_y = np.max(default_strm_2006_avg)
        max_overall = np.max([max_val_x,max_val_y])
        plt.plot([0,max_overall],[0,max_overall],label='1:1 line')
        plt.ylabel('Default Catchment-CN4.5 (mm/day)')
        plt.xlabel('CAMELS (mm/day)')
        plt.legend()
        plt.savefig(savename_strm_2006)
        plt.close()
        # make the plot for 2007 streamflow
        plt.figure()
        plt.scatter(camels_strm_2007_avg,default_strm_2007_avg,s=2)
        #plt.xlim([0,40])
        #plt.ylim([0,40])
        max_val_x = np.max(camels_strm_2007_avg)
        max_val_y = np.max(default_strm_2007_avg)
        max_overall = np.max([max_val_x,max_val_y])
        plt.plot([0,max_overall],[0,max_overall],label='1:1 line')
        plt.ylabel('Default Catchment-CN4.5 (mm/day)')
        plt.xlabel('CAMELS (mm/day)')
        plt.legend()
        plt.savefig(savename_strm_2007)
        plt.close()
        # make the plot for 2006 and 2007 streamflow
        plt.figure()
        plt.scatter(camels_strm_2006_2007_avg,default_strm_2006_2007_avg,s=2)
        #plt.xlim([0,40])
        #plt.ylim([0,40])
        max_val_x = np.max(camels_strm_2006_2007_avg)
        max_val_y = np.max(default_strm_2006_2007_avg)
        max_overall = np.max([max_val_x,max_val_y])
        plt.plot([0,max_overall],[0,max_overall],label='1:1 line')
        plt.ylabel('Default Catchment-CN4.5 (mm/day)')
        plt.xlabel('CAMELS (mm/day)')
        plt.legend()
        plt.savefig(savename_strm_2006_2007)
        plt.close()
        # make the plot for 2006 le
        plt.figure()
        plt.scatter(fluxcom_le_2006_avg,default_le_2006_avg,s=2)
        #plt.xlim([0,40])
        #plt.ylim([0,40])
        max_val_x = np.max(fluxcom_le_2006_avg)
        max_val_y = np.max(default_le_2006_avg)
        max_overall = np.max([max_val_x,max_val_y])
        plt.plot([0,max_overall],[0,max_overall],label='1:1 line')
        plt.ylabel('Default Catchment-CN4.5 (W/m2)')
        plt.xlabel('Fluxcom (W/m2)')
        plt.legend()
        plt.savefig(savename_le_2006)
        plt.close()
        # make the plot for 2007 le
        plt.figure()
        plt.scatter(fluxcom_le_2007_avg,default_le_2007_avg,s=2)
        #plt.xlim([0,40])
        #plt.ylim([0,40])
        max_val_x = np.max(fluxcom_le_2007_avg)
        max_val_y = np.max(default_le_2007_avg)
        max_overall = np.max([max_val_x,max_val_y])
        plt.plot([0,max_overall],[0,max_overall],label='1:1 line')
        plt.ylabel('Default Catchment-CN4.5 (W/m2)')
        plt.xlabel('Fluxcom (W/m2)')
        plt.legend()
        plt.savefig(savename_le_2007)
        plt.close()
        # make the plot for 2006 and 2007 le
        plt.figure()
        plt.scatter(fluxcom_le_2006_2007_avg,default_le_2006_2007_avg,s=2)
        #plt.xlim([0,40])
        #plt.ylim([0,40])
        max_val_x = np.max(fluxcom_le_2006_2007_avg)
        max_val_y = np.max(default_le_2006_2007_avg)
        max_overall = np.max([max_val_x,max_val_y])
        plt.plot([0,max_overall],[0,max_overall],label='1:1 line')
        plt.ylabel('Default Catchment-CN4.5 (W/m2)')
        plt.xlabel('Fluxcom (W/m2)')
        plt.legend()
        plt.savefig(savename_le_2006_2007)
        plt.close()
        # make the plot for 2006 vs 2007 camels strm
        plt.figure()
        plt.scatter(camels_strm_2006_avg,camels_strm_2007_avg,s=2)
        #plt.xlim([0,40])
        #plt.ylim([0,40])
        max_val_x = np.max(camels_strm_2006_avg)
        max_val_y = np.max(camels_strm_2007_avg)
        max_overall = np.max([max_val_x,max_val_y])
        plt.plot([0,max_overall],[0,max_overall],label='1:1 line')
        plt.ylabel('CAMELS avgerage 2006 (mm/day)')
        plt.xlabel('CAMELS average 2007 (mm/day)')
        plt.text(
            5,2,'Mean percent difference: {:.2f}'.format(
                camels_diff_strm_2006_2007
            )
        )
        plt.legend()
        plt.savefig(savename_camels_diff_strm)
        plt.close()
        # make the plot for 2006 vs 2007 camels strm
        plt.figure()
        plt.scatter(fluxcom_le_2006_avg,fluxcom_le_2007_avg,s=2)
        #plt.xlim([0,40])
        #plt.ylim([0,40])
        max_val_x = np.max(fluxcom_le_2006_avg)
        max_val_y = np.max(fluxcom_le_2007_avg)
        max_overall = np.max([max_val_x,max_val_y])
        plt.plot([0,max_overall],[0,max_overall],label='1:1 line')
        plt.ylabel('FluxCom avgerage 2006 (W/m2)')
        plt.xlabel('FluxCom average 2007 (W/m2)')
        plt.text(
            40,10,'Mean percent difference: {:.2f}'.format(
                fluxcom_diff_le_2006_2007
            )
        )
        plt.legend()
        plt.savefig(savename_fluxcom_diff_le)
        plt.close()
    def plot_maps(self,rmse_dfs,plots_dir,geojson_fname,exp_names,states_shp,
                  make_plots,plot_trim,get_metf,chosen_huc6s_fname):
        # extract the different error dfs
        default_df = rmse_dfs[0]
        pso_init_df = rmse_dfs[1]
        pso_final_df = rmse_dfs[2]
        waterwatch_df = rmse_dfs[3]
        # what do we want to plot? All for streamflow:
        # overall rmse for default
        default_strm_rmse = np.array(default_df.loc['strm_rmse'])
        default_strm_rmse_avg = default_strm_rmse[-1]
        default_strm_rmse = default_strm_rmse[:-1]
        # normalized rmse for default
        default_strm_rmse_norm = np.array(default_df.loc['strm_rmse_norm'])
        default_strm_rmse_norm_avg = default_strm_rmse_norm[-1]
        default_strm_rmse_norm = default_strm_rmse_norm[:-1]
        # noramlized and averaged rmse for default
        default_strm_rmse_norm_ind = np.array(
            default_df.loc['strm_rmse_norm_ind']
        )
        default_strm_rmse_norm_ind_avg = default_strm_rmse_norm_ind[-1]
        default_strm_rmse_norm_ind = default_strm_rmse_norm_ind[:-1]
        # r2 for the default
        default_strm_r2 = np.array(default_df.loc['strm_r2'])
        default_strm_r2_avg = default_strm_r2[-1]
        default_strm_r2 = default_strm_r2[:-1]
        # corr for the default
        default_strm_corr = np.array(default_df.loc['strm_corr'])
        default_strm_corr_avg = default_strm_corr[-1]
        default_strm_corr = default_strm_corr[:-1]
        # ubrmse for the default
        default_strm_ubrmse = np.array(default_df.loc['strm_ubrmse'])
        default_strm_ubrmse_avg = default_strm_ubrmse[-1]
        default_strm_ubrmse = default_strm_ubrmse[:-1]
        # nse for the default
        default_strm_nse = np.array(default_df.loc['strm_nse'])
        default_strm_nse_avg = default_strm_nse[-1]
        default_strm_nse = default_strm_nse[:-1]
        # nse ind for the default
        default_strm_nse_ind = np.array(default_df.loc['strm_nse_ind'])
        default_strm_nse_ind_avg = default_strm_nse_ind[-1]
        default_strm_nse_ind = default_strm_nse_ind[:-1]
        # mae for the default
        default_strm_mae = np.array(default_df.loc['strm_mae'])
        default_strm_mae_avg = default_strm_mae[-1]
        default_strm_mae = default_strm_mae[:-1]
        # avg prediction difference for default
        default_strm_avg_diff = np.array(default_df.loc['strm_avg_diff'])
        default_strm_avg_diff_avg = default_strm_avg_diff[-1]
        default_strm_avg_diff = default_strm_avg_diff[:-1]
        # overall rmse for firstloc
        init_strm_rmse = np.array(pso_init_df.loc['strm_rmse'])
        init_strm_rmse_avg = init_strm_rmse[-1]
        init_strm_rmse = init_strm_rmse[:-1]
        # overall rmse for final it
        final_strm_rmse = np.array(pso_final_df.loc['strm_rmse'])
        final_strm_rmse_avg = final_strm_rmse[-1]
        final_strm_rmse = final_strm_rmse[:-1]
        # normalized rmse for final
        final_strm_rmse_norm = np.array(pso_final_df.loc['strm_rmse_norm'])
        final_strm_rmse_norm_avg = final_strm_rmse_norm[-1]
        final_strm_rmse_norm = final_strm_rmse_norm[:-1]
        # noramlized and averaged rmse for final
        final_strm_rmse_norm_ind = np.array(
            pso_final_df.loc['strm_rmse_norm_ind']
        )
        final_strm_rmse_norm_ind_avg = final_strm_rmse_norm_ind[-1]
        final_strm_rmse_norm_ind = final_strm_rmse_norm_ind[:-1]
        # r2 for the final
        final_strm_r2 = np.array(pso_final_df.loc['strm_r2'])
        final_strm_r2_avg = final_strm_r2[-1]
        final_strm_r2 = final_strm_r2[:-1]
        # corr for the final
        final_strm_corr = np.array(pso_final_df.loc['strm_corr'])
        final_strm_corr_avg = final_strm_corr[-1]
        final_strm_corr = final_strm_corr[:-1]
        # ubrmse for the final
        final_strm_ubrmse = np.array(pso_final_df.loc['strm_ubrmse'])
        final_strm_ubrmse_avg = final_strm_ubrmse[-1]
        final_strm_ubrmse = final_strm_ubrmse[:-1]
        # nse for the final
        final_strm_nse = np.array(pso_final_df.loc['strm_nse'])
        final_strm_nse_avg = final_strm_nse[-1]
        final_strm_nse = final_strm_nse[:-1]
        # nse ind for the final
        final_strm_nse_ind = np.array(pso_final_df.loc['strm_nse_ind'])
        final_strm_nse_ind_avg = final_strm_nse_ind[-1]
        final_strm_nse_ind = final_strm_nse_ind[:-1]
        # avg. streamflow for default
        default_strm = np.array(default_df.loc['strm'])
        default_strm_avg = default_strm[-1]
        default_strm = default_strm[:-1]
        # avg. streamflow for first
        init_strm = np.array(pso_init_df.loc['strm'])
        init_strm_avg = init_strm[-1]
        init_strm = init_strm[:-1]
        # avg. streamflow for final
        final_strm = np.array(pso_final_df.loc['strm'])
        final_strm_avg = final_strm[-1]
        final_strm = final_strm[:-1]
        # mae normalized by defualt for final
        final_mae_norm_def = np.array(pso_final_df.loc['strm_mae_norm_def'])
        final_mae_norm_def_avg = final_mae_norm_def[-1]
        final_mae_norm_def = final_mae_norm_def[:-1]
        # changes in rmse for default to first
        diff_rmse_strm_init = init_strm_rmse - default_strm_rmse
        diff_rmse_strm_init_avg = init_strm_rmse_avg - default_strm_rmse_avg
        # changes in rmse for default to final
        diff_rmse_strm_final = final_strm_rmse - default_strm_rmse
        diff_rmse_strm_final_avg = final_strm_rmse_avg - default_strm_rmse_avg
        all_diff_rmse_strm_final = np.append(
            diff_rmse_strm_final,diff_rmse_strm_final_avg
        )
        pso_final_df.loc['change_strm_rmse'] = all_diff_rmse_strm_final
        # change in normalized rmse for default to final
        diff_rmse_norm_strm_final = (
            final_strm_rmse_norm - default_strm_rmse_norm
        )
        diff_rmse_norm_strm_final_avg = (
            final_strm_rmse_norm_avg - default_strm_rmse_norm_avg
        )
        all_diff_rmse_norm_strm_final = np.append(
            diff_rmse_norm_strm_final,diff_rmse_norm_strm_final_avg
        )
        pso_final_df.loc['change_strm_rmse_norm'] = (
            all_diff_rmse_norm_strm_final
        )
        # change in individaully normalized rmse default to final
        diff_rmse_norm_ind_strm_final = (
            final_strm_rmse_norm_ind - default_strm_rmse_norm_ind
        )
        diff_rmse_norm_ind_strm_final_avg = (
            final_strm_rmse_norm_ind_avg - default_strm_rmse_norm_ind_avg
        )
        all_diff_rmse_norm_ind_strm_final = np.append(
            diff_rmse_norm_ind_strm_final,diff_rmse_norm_ind_strm_final_avg
        )
        pso_final_df.loc['change_strm_rmse_norm_ind'] = (
            all_diff_rmse_norm_ind_strm_final
        )
        # perc. changes in rmse for default to first
        perc_diff_rmse_strm_init = (
            (init_strm_rmse - default_strm_rmse)/default_strm_rmse
        )
        perc_diff_rmse_strm_init_avg = (
            (init_strm_rmse_avg - default_strm_rmse_avg)/default_strm_rmse_avg
        )
        # perc. changes in rmse for default to final
        perc_diff_rmse_strm_final = (
            (final_strm_rmse - default_strm_rmse)/default_strm_rmse
        )
        perc_diff_rmse_strm_final_avg = (
            (final_strm_rmse_avg - default_strm_rmse_avg)/default_strm_rmse_avg
        )
        all_perc_diff_rmse_strm_final = np.append(
            perc_diff_rmse_strm_final,perc_diff_rmse_strm_final_avg
        )
        pso_final_df.loc['perc_change_strm_rmse'] = all_perc_diff_rmse_strm_final
        # change r2 default to final
        diff_r2_strm_final = final_strm_r2 - default_strm_r2
        diff_r2_strm_final_avg = final_strm_r2_avg - default_strm_r2_avg
        all_diff_r2_strm_final = np.append(
            diff_r2_strm_final,diff_r2_strm_final_avg
        )
        pso_final_df.loc['change_strm_r2'] = all_diff_r2_strm_final
        # change corr default to final
        diff_corr_strm_final = final_strm_corr - default_strm_corr
        diff_corr_strm_final_avg = final_strm_corr_avg - default_strm_corr_avg
        all_diff_corr_strm_final = np.append(
            diff_corr_strm_final,diff_corr_strm_final_avg
        )
        pso_final_df.loc['change_strm_corr'] = all_diff_corr_strm_final
        # perc change ubrmse default to final
        perc_diff_ubrmse_strm_final = (
            (final_strm_ubrmse - default_strm_ubrmse)/default_strm_ubrmse
        )
        perc_diff_ubrmse_strm_final_avg = (
            (final_strm_ubrmse_avg - default_strm_ubrmse_avg)/default_strm_ubrmse_avg
        )
        all_perc_diff_ubrmse_strm_final = np.append(
            perc_diff_ubrmse_strm_final,perc_diff_ubrmse_strm_final_avg
        )
        pso_final_df.loc['perc_change_strm_ubrmse'] = (
            all_perc_diff_ubrmse_strm_final
        )
        # change nse defautl to final
        diff_nse_strm_final = final_strm_nse - default_strm_nse
        diff_nse_strm_final_avg = final_strm_nse_avg - default_strm_nse_avg
        all_diff_nse_strm_final = np.append(
            diff_nse_strm_final,diff_nse_strm_final_avg
        )
        pso_final_df.loc['change_strm_nse'] = all_diff_nse_strm_final
        # change in individaully averaged nse from default to final
        diff_nse_ind_strm_final = (
            final_strm_nse_ind - default_strm_nse_ind
        )
        diff_nse_ind_strm_final_avg = (
            final_strm_nse_ind_avg - default_strm_nse_ind_avg
        )
        all_diff_nse_ind_strm_final = np.append(
            diff_nse_ind_strm_final,diff_nse_ind_strm_final_avg
        )
        pso_final_df.loc['change_strm_nse_ind'] = (
            all_diff_nse_ind_strm_final
        )
        # changes in avg. strm default to first
        diff_strm_init = init_strm - default_strm
        diff_strm_init_avg = init_strm_avg - default_strm_avg
        # changes in avg. strm deafult to final
        diff_strm_final = final_strm - default_strm
        diff_strm_final_avg = final_strm_avg - default_strm_avg
        all_diff_strm_final = np.append(
            diff_strm_final,diff_strm_final_avg
        )
        pso_final_df.loc['diff_strm'] = all_diff_strm_final
        # perc. changes in avg. strm default to first
        perc_diff_strm_init = (
            (init_strm - default_strm)/default_strm
        )
        perc_diff_strm_init_avg = (
            (init_strm_avg - default_strm_avg)/default_strm_avg
        )
        # perc. chagnes in avg. strm default to final
        perc_diff_strm_final = (
            (final_strm - default_strm)/default_strm
        )
        perc_diff_strm_final_avg = (
            (final_strm_avg - default_strm_avg)/default_strm_avg
        )
        # default normalized rmse
        default_strm_rmse_norm = np.array(default_df.loc['strm_rmse_norm'])
        default_strm_rmse_norm_avg = default_strm_rmse_norm[-1]
        default_strm_rmse_norm = default_strm_rmse_norm[:-1]
        # final normalized rmse
        final_strm_rmse_norm = np.array(pso_final_df.loc['strm_rmse_norm'])
        final_strm_rmse_norm_avg = final_strm_rmse_norm[-1]
        final_strm_rmse_norm = final_strm_rmse_norm[:-1]
        # difference in normalized rmse
        diff_strm_rmse_norm_final = (
            final_strm_rmse_norm - default_strm_rmse_norm
        )
        diff_strm_rmse_norm_final_avg = (
            final_strm_rmse_norm_avg - default_strm_rmse_norm_avg
        )
        all_diff_strm_rmse_norm_final = np.append(
            diff_strm_rmse_norm_final,diff_strm_rmse_norm_final_avg
        )
        #print(diff_rmse_norm_ind_strm_final_avg)
        pso_final_df.loc['diff_strm_rmse_norm'] = all_diff_strm_rmse_norm_final
        # stream observations
        huc6s = np.array(default_df.columns)
        strm_obs = np.array(pso_final_df.loc['strm_obs'])
        strm_obs_avg = strm_obs[-1]
        strm_obs = strm_obs[:-1]
        print('for mae normalized by default mae')
        print(final_mae_norm_def_avg)
        deleting_vals = copy.deepcopy(final_mae_norm_def)
        num_delete = 4
        iteration = 1
        while iteration <= num_delete:
            deleting_idx = np.where(deleting_vals != np.min(deleting_vals))
            deleting_idx_opp = np.where(deleting_vals == np.min(deleting_vals))
            deleting_vals = deleting_vals[deleting_idx]
            deleting_avg = np.mean(deleting_vals)
            print('basin: {}'.format(huc6s[deleting_idx_opp]))
            print('default mae: {}'.format(default_strm_mae[deleting_idx_opp]))
            print('{} deleted: {}'.format(iteration,deleting_avg))
            iteration += 1
        #sys.exit()
        print('for rmse normalized by average streamflow')
        print(final_strm_rmse_norm_ind)
        print(diff_rmse_norm_ind_strm_final_avg)
        small_obs = np.where(strm_obs < 0.0105)
        big_obs = np.array(np.where(strm_obs > 0.0105))[0]
        print(small_obs)
        print(big_obs)
        print(huc6s[small_obs])
        print(strm_obs[small_obs])
        print(final_strm_rmse_norm_ind[small_obs])
        rmse_big_obs_final = final_strm_rmse_norm_ind[big_obs]
        rmse_big_obs_final_avg = np.mean(rmse_big_obs_final)
        rmse_big_obs_default = default_strm_rmse_norm_ind[big_obs]
        rmse_big_obs_default_avg = np.mean(rmse_big_obs_default)
        diff_strm_rmse_norm_ind_final_big = (
            rmse_big_obs_final_avg - rmse_big_obs_default_avg
        )
        #all_diff_strm_rmse_norm_ind_final[-1] = (
        #    diff_strm_rmse_norm_ind_final_big
        #)
        print(diff_strm_rmse_norm_ind_final_big)
        #sys.exit()
        #pso_final_df.loc['diff_strm_rmse_norm'] = all_diff_strm_rmse_norm_final
        # default ratio
        default_strm_ratio = default_strm/strm_obs
        default_strm_ratio_avg = default_strm_avg/strm_obs_avg
        # let's plot streamflow versus area-weighted precip to see how things
        # are looking on the measurement side
        huc6s = gpd.read_file(geojson_fname)
        huc6s_area_df = pd.read_csv(chosen_huc6s_fname)
        hucs = list(default_df.columns)
        print(huc6s)
        print(huc6s_area_df)
        print(default_df.loc['precip_obs'])
        #sys.exit()
        hucs = hucs[:-1]
        if get_metf:
            areas = np.zeros(len(hucs))
            for h,huc in enumerate(hucs):
                this_idx = np.where(
                    huc6s_area_df['camel'] == huc
                )[0][0]
                areas[h] = huc6s_area_df['area'].iloc[this_idx]
            obs_precip = np.array(default_df.loc['precip_obs'])[:-1]
            area_weighted_precip = areas*obs_precip
            obs_strm = np.array(default_df.loc['strm_obs'])[:-1]
            print(default_df.loc['precip_obs'])
            print(obs_strm)
            sys.exit()
            norm_area_weighted_precip = obs_strm/area_weighted_precip
            norm_area_weighted_precip_avg = np.mean(norm_area_weighted_precip)
            x_vals = np.arange(len(obs_precip))
            save_name = os.path.join(
                plots_dir,
                'precip_runoff_ratio_scatter.png'
            )
            plt.figure(figsize=(40,5))
            plt.scatter(x_vals,norm_area_weighted_precip)
            plt.xticks(x_vals,hucs,rotation=90)
            plt.savefig(save_name,bbox_inches='tight')
            plt.close()
        if make_plots:
            # now let's get the shapes that we need for plotting
            if not plot_trim:
                # now let's get everythin in arrays for proper plotting
                names = [
                    'default_strm_rmse','init_strm_rmse','final_strm_rmse',
                    'default_strm','init_strm','final_strm',
                    'diff_rmse_strm_init','diff_rmse_strm_final',
                    'perc_diff_rmse_strm_init','perc_diff_rmse_strm_final',
                    'diff_strm_init','diff_strm_final',
                    'perc_diff_strm_init','perc_diff_strm_final',
                    'default_strm_avg_diff',
                    'default_strm_r2','default_strm_corr','default_strm_ubrmse',
                    'default_strm_nse','pso_strm_r2','pso_strm_corr',
                    'pso_strm_ubrmse','pso_strm_nse','diff_r2_strm_final',
                    'diff_corr_strm_final','perc_diff_ubrmse_strm_final',
                    'diff_nse_strm_final'
                ]
                vals = [
                    default_strm_rmse,init_strm_rmse,final_strm_rmse,
                    default_strm,init_strm,final_strm,
                    diff_rmse_strm_init,diff_rmse_strm_final,
                    perc_diff_rmse_strm_init,perc_diff_rmse_strm_final,
                    diff_strm_init,diff_strm_final,
                    perc_diff_strm_init,perc_diff_strm_final,
                    default_strm_avg_diff,
                    default_strm_r2,default_strm_corr,default_strm_ubrmse,
                    default_strm_nse,final_strm_r2,final_strm_corr,
                    final_strm_ubrmse,final_strm_nse,diff_r2_strm_final,
                    diff_corr_strm_final,perc_diff_ubrmse_strm_final,
                    diff_nse_strm_final

                ]
                avgs = [
                    default_strm_rmse_avg,init_strm_rmse_avg,final_strm_rmse_avg,
                    default_strm_avg,init_strm_avg,final_strm_avg,
                    diff_rmse_strm_init_avg,diff_rmse_strm_final_avg,
                    perc_diff_rmse_strm_init_avg,perc_diff_rmse_strm_final_avg,
                    diff_strm_init_avg,diff_strm_final_avg,
                    perc_diff_strm_init_avg,perc_diff_strm_final_avg,
                    default_strm_avg_diff_avg,
                    default_strm_r2_avg,default_strm_corr_avg,default_strm_ubrmse_avg,
                    default_strm_nse_avg,final_strm_r2_avg,final_strm_corr_avg,
                    final_strm_ubrmse_avg,final_strm_nse_avg,diff_r2_strm_final_avg,
                    diff_corr_strm_final_avg,perc_diff_ubrmse_strm_final_avg,
                    diff_nse_strm_final_avg
                ]
                types = [
                    'strm_rmse','strm_rmse','strm_rmse',
                    'strm','strm','strm',
                    'diff_strm_rmse','diff_strm_rmse',
                    'perc_diff_strm_rmse','perc_diff_strm_rmse',
                    'diff_strm','diff_strm',
                    'perc_diff_strm','perc_diff_strm',
                    'diff_strm',
                    'strm_r2','strm_r2','strm_rmse',
                    'strm_r2','strm_r2','strm_r2',
                    'strm_rmse','strm_r2','perc_diff_strm_rmse',
                    'perc_diff_strm_rmse','perc_diff_strm_rmse',
                    'perc_diff_strm_rmse'
                ]
                cmaps = {
                    'strm_rmse':'rainbow',
                    'strm':'rainbow',
                    'diff_strm_rmse':'bwr',
                    'perc_diff_strm_rmse':'bwr',
                    'diff_strm':'bwr',
                    'perc_diff_strm':'bwr',
                    'strm_r2':'rainbow'
                }
                vmins = {
                    'strm_rmse':0,
                    'strm':0,
                    'diff_strm_rmse':-2,
                    'perc_diff_strm_rmse':-.25,
                    'diff_strm':-2,
                    'perc_diff_strm':-.25,
                    'strm_r2':0
                }
                vmaxs = {
                    'strm_rmse':2,
                    'strm':2,
                    'diff_strm_rmse':2,
                    'perc_diff_strm_rmse':.25,
                    'diff_strm':2,
                    'perc_diff_strm':.25,
                    'strm_r2':1
                }
            if plot_trim:
                names = [
                    #'final_strm_rmse',
                    #'final_strm_rmse_norm_ind',
                    #'diff_rmse_norm_ind_strm_final'
                    #'default_strm'
                    #'diff_strm_final'
                    #'default_strm_obs',
                    #'default_strm_ratio'
                    #'norm_area_weighted_precip'
                    'final_mae_norm_def'
                ]
                vals = [
                    #final_strm_rmse,
                    #final_strm_rmse_norm_ind,
                    #diff_rmse_norm_ind_strm_final
                    #default_strm
                    #diff_strm_final
                    #strm_obs,
                    #default_strm_ratio
                    #norm_area_weighted_precip
                    final_mae_norm_def
                ]
                avgs = [
                    #final_strm_rmse_avg,
                    #final_strm_rmse_norm_ind_avg,
                    #diff_rmse_norm_ind_strm_final_avg
                    #default_strm_avg
                    #diff_strm_final_avg
                    #strm_obs_avg,
                    #default_strm_ratio_avg
                    #norm_area_weighted_precip_avg
                    final_mae_norm_def_avg
                ]
                types = names
                cmaps = {
                    'final_strm_rmse':'rainbow',
                    'final_strm_rmse_norm_ind':'rainbow',
                    'diff_rmse_norm_ind_strm_final':'bwr',
                    'default_strm':'rainbow',
                    'diff_strm_final':'bwr',
                    'default_strm_obs':'rainbow',
                    'default_strm_ratio':'rainbow',
                    'norm_area_weighted_precip':'rainbow',
                    'final_mae_norm_def':'rainbow'
                }
                vmins = {
                    'final_strm_rmse':0,
                    'final_strm_rmse_norm_ind':0,
                    'diff_rmse_norm_ind_strm_final':-5,
                    'default_strm':0,
                    'diff_strm_final':-1,
                    'default_strm_obs':0,
                    'default_strm_ratio':0,
                    'norm_area_weighted_precip':0,
                    'final_mae_norm_def':0
                }
                vmaxs = {
                    'final_strm_rmse':1.5,
                    'final_strm_rmse_norm_ind':1.5,
                    'diff_rmse_norm_ind_strm_final':5,
                    'default_strm':0.1,
                    'diff_strm_final':1,
                    'default_strm_obs':0.1,
                    'default_strm_ratio':30,
                    'norm_area_weighted_precip':1e12,
                    'final_mae_norm_def':3
                }
            print('reading states')
            states = gpd.read_file(states_shp)
            states = states.to_crs('EPSG:4326')
            # get rid of non-conus states since not considering
            non_conus = ['HI','VI','MP','GU','AK','AS','PR']
            states_conus = states
            print('looping non conus')
            for n in non_conus:
                states_conus = states_conus[states_conus.STUSPS != n]
            # get a list of all the hucs
            all_hucs = np.array(default_df.columns)
            all_hucs = all_hucs[:-1]
            for n,name in enumerate(names):
                print(name)
                fig,ax = plt.subplots()
                divider = make_axes_locatable(ax)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                state_ids = list(states_conus['GEOID'])
                for s,sid in enumerate(state_ids):
                    this_geom = states_conus['geometry'].iloc[s]
                    try:
                        xs,ys = this_geom.exterior.xy
                        ax.fill(xs,ys,fc='none',ec='k',linewidth=.1)
                    except:
                        for geom in this_geom.geoms:
                            xs,ys = geom.exterior.xy
                            ax.fill(xs,ys,fc='none',ec='k',linewidth=.1)
                # get our normalize function for getting colors
                norm = mpl.colors.Normalize(
                    vmin=vmins[types[n]],vmax=vmaxs[types[n]]
                )
                this_cmap = mpl.cm.get_cmap(cmaps[types[n]])
                for h,huc in enumerate(all_hucs):
                    idx = np.where(
                        huc6s['hru_id'] == huc
                    )[0][0]
                    this_geom = huc6s['geometry'].iloc[idx]
                    this_val = vals[n][h]
                    this_val_norm = norm(this_val)
                    this_color = this_cmap(this_val_norm)
                    if this_geom.geom_type == 'Polygon':
                        xs,ys = this_geom.exterior.xy
                        ax.fill(xs,ys,fc=this_color,ec='k',linewidth=0.2)
                    elif this_geom.geom_type == 'MultiPolygon':
                        for this_this_geom in this_geom.geoms:
                            xs,ys = this_this_geom.exterior.xy
                            ax.fill(xs,ys,fc=this_color,ec='k',linewidth=0.2)
                    else:
                        raise IOError('Shape is not a polygon')
                ax.text(
                    -127+2,20+4,'Average {name}: {val:.2f}'.format(
                        name=names[n],val=avgs[n]
                    ),
                    bbox=dict(facecolor='white')
                )
                fig.colorbar(
                    mpl.cm.ScalarMappable(norm=norm, cmap=this_cmap),
                    cax=cax, orientation='vertical'
                )
                this_savename = os.path.join(
                    plots_dir,
                    '{name}_{exp}_iteration.png'.format(
                        name=names[n],exp=exp_names[2]
                    )
                )
                plt.savefig(this_savename,dpi=350,bbox_inches='tight')
                plt.close()
