import sys
sys.path.append('/shared/pso/step_5_analyze_outputs/')
from get_timeseries import get_timeseries
from analyze_timeseries_watershed_scale import analyze_watershed
from analyze_timeseries_pixel_scale import analyze_pix
import datetime
from dateutil.relativedelta import *
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import geopandas as gpd
import os
import numpy as np
import pandas as pd
import netCDF4 as nc

class water_balance:
    def __init__(self,start,end):
        self.start = start
        self.end = end
    def get_catch_met(self,pixels_fname,exp_dirs,exp_names,save_timeseries,
                      load_timeseries,is_default_experiment,
                      read_met_forcing,saved_time_dir,intersection_fname,
                      fluxcom_dir,strm_dir):
        print('getting pixel scale products')
        # initiate a class of get_timeseries
        g = get_timeseries(self.start,self.end)
        pixels = g.get_pixels(pixels_fname)
        # we are going to load first for the default case
        default_timeseries_pixel = g.get_catch_timeseries(
            exp_dirs[0],exp_names[0],save_timeseries[0],load_timeseries[0],
            is_default_experiment[0],
            read_met_forcing[0],saved_time_dir
        )
        fluxcom_timeseries = g.get_fluxcom_timeseries(fluxcom_dir)
        streamflow_timeseries = g.get_strm_timeseries(strm_dir)
        a_water = analyze_watershed()
        intersection_info = a_water.get_intersection(intersection_fname)
        # need a fake exp name
        default_exp_names_repeat = [exp_names[0],exp_names[0],exp_names[0]]
        default_timeseries_pixel_repeat = {
            exp_names[0]:default_timeseries_pixel,
            exp_names[0]:default_timeseries_pixel,
            exp_names[0]:default_timeseries_pixel,
        }
        # we want watershed scale met forcing
        get_metf = True
        default_watershed_timeseries = a_water.get_model_preds_watershed(
            self.start,self.end,default_exp_names_repeat,default_timeseries_pixel_repeat,streamflow_timeseries,
            fluxcom_timeseries,intersection_info,get_metf
        )
        # now lets do for the pf case
        pf_timeseries_pixel = g.get_catch_timeseries(
            exp_dirs[1],exp_names[1],save_timeseries[1],load_timeseries[1],
            is_default_experiment[1],
            read_met_forcing[1],saved_time_dir
        )
        pf_exp_names_repeat = [exp_names[1],exp_names[1],exp_names[1]]
        pf_timeseries_pixel_repeat = {
            exp_names[1]:pf_timeseries_pixel,
            exp_names[1]:pf_timeseries_pixel,
            exp_names[1]:pf_timeseries_pixel,
        }
        pf_watershed_timeseries = a_water.get_model_preds_watershed(
            self.start,self.end,pf_exp_names_repeat,pf_timeseries_pixel_repeat,streamflow_timeseries,
            fluxcom_timeseries,intersection_info,get_metf
        )
        outs = [
            default_watershed_timeseries,pf_watershed_timeseries,streamflow_timeseries,
            default_timeseries_pixel,pf_timeseries_pixel,intersection_info
        ]
        return outs
    def get_water_balance_timeseries(self,watershed_timeseries,
                                     waterwatch_timeseries):
        print('getting watershed scale products')
        # get the different elements of the water balance
        merra2_precip = watershed_timeseries['default_rainfsnowf']
        fluxcom_le = watershed_timeseries['fluxcom_le']
        waterwatch_strm = waterwatch_timeseries
        waterwatch_strm.index.names = ['date']
        # we need fluxcom in mm/month, not W/m2
        fluxcom_le = fluxcom_le/28.94 # now in mm/day
        fluxcom_months = list(fluxcom_le.index)
        for d,date in enumerate(fluxcom_months):
            this_date = datetime.datetime.strptime(str(date),'%Y%m')
            next_mon = this_date + relativedelta(months=1)
            diff = next_mon - this_date
            diff = diff.days
            fluxcom_le.loc[date] = fluxcom_le.loc[date]*diff # now in mm/month
        # waterwatch columns need to be ints not strings
        wat_cols = list(waterwatch_strm.columns)
        new_wat_cols = [
            int(col) for col in wat_cols
        ]
        waterwatch_strm = waterwatch_strm.set_axis(new_wat_cols,axis=1)
        incoming = copy.deepcopy(merra2_precip)
        outgoing = fluxcom_le+ waterwatch_strm
        diff = incoming - outgoing
        return [
            merra2_precip,fluxcom_le,waterwatch_strm,incoming,outgoing,diff
        ]
    def plot_water_balance(self,metrics,plots_dir,huc6_shape_fname,
                           states_shp_fname):
        start_str = self.start.strftime('%Y%m')
        end_str = self.end.strftime('%Y%m')
        plot_timeseries = True
        plot_map = True
        merra2_precip = metrics[0]
        fluxcom_le = metrics[1]
        waterwatch_strm = metrics[2]
        incoming = metrics[3]
        outgoing = metrics[4]
        diff = metrics[5]
        if plot_timeseries:
            # let's first plot the three different products for each watershed
            print('plotting product timeseries')
            watersheds = list(merra2_precip.columns)
            months = list(merra2_precip.index)
            months_dt = [
                datetime.datetime.strptime(str(mon),'%Y%m') for mon in months
            ]
            for w,wat in enumerate(watersheds):
                this_savename = os.path.join(
                    plots_dir,'water_balance_products_{}_{}_{}.png'.format(
                        wat,start_str,end_str
                    )
                )
                plt.figure()
                plt.plot(months_dt,merra2_precip[wat],label='MERRA2 precip')
                plt.plot(months_dt,fluxcom_le[wat],label='FluxCom LE')
                plt.plot(
                    months_dt,waterwatch_strm[wat],label='WaterWatch streamflow'
                )
                plt.legend()
                plt.savefig(this_savename,dpi=300,bbox_inches='tight')
                plt.close()
            # and now let's plot incoming, outgoing, diff for each watershed
            print('plotting water balance timeseries')
            for w,wat in enumerate(watersheds):
                this_savename = os.path.join(
                    plots_dir,'water_balance_{}_{}_{}.png'.format(
                        wat,start_str,end_str
                    )
                )
                plt.figure()
                plt.plot(
                    months_dt,incoming[wat],label='incoming',c='g'
                )
                plt.plot(months_dt,outgoing[wat],label='outgoing',c='r')
                plt.plot(
                    months_dt,diff[wat],label='difference',c='k'
                )
                plt.axhline(c='b',linestyle='--')
                plt.legend()
                plt.savefig(this_savename,dpi=300,bbox_inches='tight')
                plt.close()
        # now let's plot the diff as a map for each watershed
        if plot_map:
            # now let's get the shapes that we need for plotting
            huc6s = gpd.read_file(huc6_shape_fname)
            states = gpd.read_file(states_shp_fname)
            states = states.to_crs('EPSG:4326')
            # we need the average difference for each watershed
            avg_diff_wat = list(diff.mean())
            avg_diff_all = np.mean(avg_diff_wat)
            names = ['diff_water_balance']
            vals = [avg_diff_wat]
            avgs = [avg_diff_all]
            types = ['water_balance_diff']
            cmaps = {
                'water_balance_diff':'bwr'
            }
            vmins = {
                'water_balance_diff':-60
            }
            vmaxs = {
                'water_balance_diff':60
            }
            # get rid of non-conus states since not considering
            non_conus = ['HI','VI','MP','GU','AK','AS','PR']
            states_conus = states
            print('looping non conus')
            for n in non_conus:
                states_conus = states_conus[states_conus.STUSPS != n]
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
                        ax.fill(xs,ys,fc='none',ec='k',linewidth=.2)
                    except:
                        for geom in this_geom.geoms:
                            xs,ys = geom.exterior.xy
                            ax.fill(xs,ys,fc='none',ec='k',linewidth=.2)
                # get a list of all the hucs
                all_hucs = list(huc6s['huc6'])
                # get our normalize function for getting colors
                print('normalizing and getting cmap')
                print(types[n])
                norm = mpl.colors.Normalize(
                    vmin=vmins[types[n]],vmax=vmaxs[types[n]]
                )
                this_cmap = mpl.cm.get_cmap(cmaps[types[n]])
                for h,huc in enumerate(all_hucs):
                    this_geom = huc6s['geometry'].iloc[h]
                    xs,ys = this_geom.exterior.xy
                    this_val = vals[n][h]
                    this_val_norm = norm(this_val)
                    this_color = this_cmap(this_val_norm)
                    ax.fill(xs,ys,fc=this_color,ec='k',linewidth=0.5)
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
                    '{}_huc6_map_{}_{}.png'.format(
                        names[n],start_str,end_str
                    )
                )
                plt.savefig(this_savename,dpi=350,bbox_inches='tight')
                plt.close()
    def plot_unit_conversion(self,pixel_timeseries,plots_dir,intersection_info,
                             streamflow_timeseries,exp_name):
        plot_strm_precip_kg_m2_s = True
        plot_incoming_outgoing_kg_m2_s = True
        plot_incoming_outgoing_mm_day = True
        plot_incoming_outgoing_mm_month = True
        plot_incoming_outgoing_mm_month_huc6 = True
        plot_waterwatch_vs_catchment_mm_month_huc6 = True
        print('plotting unit conversion progression')
        # let's get the different things we need
        native_precip = pixel_timeseries['rainfsnowf'] # in kg/m2/s
        precip_kg_m2_s = native_precip
        native_runoff = pixel_timeseries['runoff'] # in kg/m2/s
        runoff_kg_m2_s = native_runoff
        native_baseflow = pixel_timeseries['baseflow'] # in kg/m2/s
        baseflow_kg_m2_s = native_baseflow
        strm_kg_m2_s = runoff_kg_m2_s + baseflow_kg_m2_s
        native_le = pixel_timeseries['le']
        le_w_m2 = native_le
        le_kg_m2_s = le_w_m2/28.94/86400
        incoming_kg_m2_s = precip_kg_m2_s
        outgoing_kg_m2_s = le_kg_m2_s + strm_kg_m2_s
        diff_strm_precip_kg_m2_s = precip_kg_m2_s - strm_kg_m2_s
        diff_incoming_outgoing_kg_m2_s = incoming_kg_m2_s - outgoing_kg_m2_s
        # now let's plot here with the native outputs
        times_daily = list(strm_kg_m2_s.index)
        pixels_all = list(strm_kg_m2_s.columns)
        if plot_strm_precip_kg_m2_s:
            for p,pix in enumerate(pixels_all):
                print('precip vs strm for pixel {}'.format(pix))
                this_fname = (
                    'MERRA2_vs_Catchment_strm_kg_m2_s_{exp}_{pix}' +
                    '_{start}_{end}.png'
                )
                this_savename = os.path.join(
                    plots_dir,
                    (
                        this_fname.format(
                            exp=exp_name,
                            pix=pix,
                            start=self.start.strftime('%Y%m%d'),
                            end=self.end.strftime('%Y%m%d')
                        )
                    )
                )
                this_avg_diff = np.mean(diff_strm_precip_kg_m2_s[pix])
                fig,ax = plt.subplots()
                plt.plot(times_daily,runoff_kg_m2_s[pix],c='y',label='Catchment-CN runoff')
                plt.plot(
                    times_daily,baseflow_kg_m2_s[pix],c='g',label='Catchment-CN baseflow'
                )
                plt.plot(times_daily,precip_kg_m2_s[pix],c='k',label='MERRA2 precip')
                plt.plot(
                    times_daily,strm_kg_m2_s[pix],c='b',label='Catchment-CN streamflow'
                )
                plt.plot(times_daily,diff_strm_precip_kg_m2_s[pix],c='r',label='precip - streamflow')
                plt.legend(loc='upper right')
                plt.xlabel('date')
                plt.ylabel('kg/m2/s')
                plt.text(
                    0.01,0.99,'Average difference: {:.2f}'.format(this_avg_diff),
                    ha='left',va='top',transform=ax.transAxes
                )
                plt.savefig(this_savename)
                plt.close()
        if plot_incoming_outgoing_kg_m2_s:
            this_fname = (
                'MERRA2_incoming_vs_Catchment_outgoing_kg_m2_s_' +
                '{exp}_{pix}_{start}_{end}.png'
            )
            for p,pix in enumerate(pixels_all):
                print('plotting kg/m2/s for pixel {}'.format(pix))
                this_savename = os.path.join(
                    plots_dir,
                    (
                        this_fname.format(
                            exp=exp_name,
                            pix=pix,
                            start=self.start.strftime('%Y%m%d'),
                            end=self.end.strftime('%Y%m%d')
                        )
                    )
                )
                this_avg_diff = np.mean(diff_incoming_outgoing_kg_m2_s[pix])
                fig,ax = plt.subplots()
                plt.plot(times_daily,incoming_kg_m2_s[pix],c='b',label='incoming')
                plt.plot(
                    times_daily,outgoing_kg_m2_s[pix],c='r',label='outgoing'
                )
                plt.plot(
                    times_daily,diff_incoming_outgoing_kg_m2_s[pix],c='k',
                    label='incoming - outgoing'
                )
                plt.legend(loc='upper right')
                plt.xlabel('date')
                plt.ylabel('kg/m2/s')
                plt.text(
                    0.01,0.99,'Average difference: {:.7f}'.format(this_avg_diff),
                    ha='left',va='top',transform=ax.transAxes
                )
                plt.savefig(this_savename)
                plt.close()
        # let's convert to mm/day and check again
        precip_mm_day = precip_kg_m2_s*86400
        strm_mm_day = strm_kg_m2_s*86400
        le_mm_day = le_kg_m2_s*86400
        incoming_mm_day = precip_mm_day
        outgoing_mm_day = strm_mm_day + le_mm_day
        diff_incoming_outgoing_mm_day = incoming_mm_day - outgoing_mm_day
        if plot_incoming_outgoing_mm_day:
            this_fname = (
                'MERRA2_incoming_vs_Catchment_outgoing_mm_day_' +
                '{exp}_{pix}_{start}_{end}.png'
            )
            for p,pix in enumerate(pixels_all):
                print('plotting mm/day for pixel {}'.format(pix))
                this_savename = os.path.join(
                    plots_dir,
                    (
                        this_fname.format(
                            exp=exp_name,
                            pix=pix,
                            start=self.start.strftime('%Y%m%d'),
                            end=self.end.strftime('%Y%m%d')
                        )
                    )
                )
                this_avg_diff = np.mean(diff_incoming_outgoing_mm_day[pix])
                fig,ax = plt.subplots()
                plt.plot(times_daily,incoming_mm_day[pix],c='b',label='incoming')
                plt.plot(
                    times_daily,outgoing_mm_day[pix],c='r',label='outgoing'
                )
                plt.plot(
                    times_daily,diff_incoming_outgoing_mm_day[pix],c='k',
                    label='incoming - outgoing'
                )
                plt.legend(loc='upper right')
                plt.xlabel('date')
                plt.ylabel('mm/day')
                plt.text(
                    0.01,0.99,'Average difference: {:.7f}'.format(this_avg_diff),
                    ha='left',va='top',transform=ax.transAxes
                )
                plt.savefig(this_savename)
                plt.close()
        # let's convert to mm/month and check
        # goign to have to hold our new monthly datasets
        precip_mm_month = pd.DataFrame(columns=pixels_all)
        le_mm_month = pd.DataFrame(columns=pixels_all)
        strm_mm_month = pd.DataFrame(columns=pixels_all)
        curr = copy.deepcopy(self.start)
        first_day = curr.strftime('%Y-%m-%d')
        last_month = curr.strftime('%m')
        while curr <= (self.end + datetime.timedelta(days=1)):
            this_date_fmt = curr.strftime('%Y-%m-%d')
            this_month = curr.strftime('%m')
            if this_month != last_month:
                last_day = curr + datetime.timedelta(days=-1)
                last_day_str = last_day.strftime('%Y-%m-%d')
                this_month_str = last_day.strftime('%Y%m')
                # for precip
                this_month_precip_mm_day = precip_mm_day.loc[
                    first_day:last_day_str
                ]
                this_month_precip_mm_month = this_month_precip_mm_day.sum()
                precip_mm_month.loc[this_month_str] = (
                    this_month_precip_mm_month
                )
                # for le
                this_month_le_mm_day = le_mm_day.loc[
                    first_day:last_day_str
                ]
                this_month_le_mm_month = this_month_le_mm_day.sum()
                le_mm_month.loc[this_month_str] = (
                    this_month_le_mm_month
                )
                # for strm
                this_month_strm_mm_day = strm_mm_day.loc[
                    first_day:last_day_str
                ]
                this_month_strm_mm_month = this_month_strm_mm_day.sum()
                strm_mm_month.loc[this_month_str] = (
                    this_month_strm_mm_month
                )
                # update last month holder
                first_day = curr.strftime('%Y-%m-%d')
                last_month = curr.strftime('%m')
            # update curr and mvoe to next step
            curr += datetime.timedelta(days=1)
        incoming_mm_month = precip_mm_month
        outgoing_mm_month = le_mm_month + strm_mm_month
        diff_incoming_outgoing_mm_month = incoming_mm_month - outgoing_mm_month
        times_monthly = list(incoming_mm_month.index)
        if plot_incoming_outgoing_mm_month:
            this_fname = (
                'MERRA2_incoming_vs_Catchment_outgoing_mm_month_' +
                '{exp}_{pix}_{start}_{end}.png'
            )
            for p,pix in enumerate(pixels_all):
                print('plotting mm/month for pixel {}'.format(pix))
                this_savename = os.path.join(
                    plots_dir,
                    (
                        this_fname.format(
                            exp=exp_name,
                            pix=pix,
                            start=self.start.strftime('%Y%m%d'),
                            end=self.end.strftime('%Y%m%d')
                        )
                    )
                )
                this_avg_diff = np.mean(diff_incoming_outgoing_mm_month[pix])
                fig,ax = plt.subplots()
                plt.plot(times_monthly,incoming_mm_month[pix],c='b',label='incoming')
                plt.plot(
                    times_monthly,outgoing_mm_month[pix],c='r',label='outgoing'
                )
                plt.plot(
                    times_monthly,diff_incoming_outgoing_mm_month[pix],c='k',
                    label='incoming - outgoing'
                )
                plt.legend(loc='upper right')
                plt.xlabel('date')
                plt.ylabel('mm/day')
                plt.text(
                    0.01,0.99,'Average difference: {:.7f}'.format(this_avg_diff),
                    ha='left',va='top',transform=ax.transAxes
                )
                plt.savefig(this_savename)
                plt.close()
        # now we gotta get to the basin scale
        hucs = [
            int(h) for h in list(intersection_info.keys())
        ]
        precip_mm_month_huc6 = pd.DataFrame(columns=hucs,index=times_monthly)
        le_mm_month_huc6 = pd.DataFrame(columns=hucs,index=times_monthly)
        strm_mm_month_huc6 = pd.DataFrame(columns=hucs,index=times_monthly)
        for h,huc in enumerate(hucs):
            this_huc_str = str(huc)
            this_huc_tiles = intersection_info[this_huc_str.zfill(6)][0]
            this_huc_perc = intersection_info[str(this_huc_str.zfill(6))][1]
            # average for precip
            this_tiles_precip = precip_mm_month[this_huc_tiles]
            this_huc_precip = np.average(
                this_tiles_precip,axis=1,weights=this_huc_perc
            )
            precip_mm_month_huc6[huc] = this_huc_precip
            # average for le
            this_tiles_le = le_mm_month[this_huc_tiles]
            this_huc_le = np.average(
                this_tiles_le,axis=1,weights=this_huc_perc
            )
            le_mm_month_huc6[huc] = this_huc_le
            # average for strm
            this_tiles_strm = strm_mm_month[this_huc_tiles]
            this_huc_strm = np.average(
                this_tiles_strm,axis=1,weights=this_huc_perc
            )
            strm_mm_month_huc6[huc] = this_huc_strm
        incoming_mm_month_huc6 = precip_mm_month_huc6
        outgoing_mm_month_huc6 = le_mm_month_huc6 + strm_mm_month_huc6
        diff_incoming_outgoing_mm_month_huc6 = (
            incoming_mm_month_huc6 - outgoing_mm_month_huc6
        )
        if plot_incoming_outgoing_mm_month_huc6:
            this_fname = (
                'MERRA2_incoming_vs_Catchment_outgoing_mm_month_huc6_' +
                '{exp}_{huc}_{start}_{end}.png'
            )
            for h,huc in enumerate(hucs):
                print('plotting mm/month for huc {}'.format(huc))
                this_savename = os.path.join(
                    plots_dir,
                    (
                        this_fname.format(
                            exp=exp_name,
                            huc=huc,
                            start=self.start.strftime('%Y%m%d'),
                            end=self.end.strftime('%Y%m%d')
                        )
                    )
                )
                this_avg_diff = np.mean(
                    diff_incoming_outgoing_mm_month_huc6[huc]
                )
                fig,ax = plt.subplots()
                plt.plot(times_monthly,incoming_mm_month_huc6[huc],c='b',label='incoming')
                plt.plot(
                    times_monthly,outgoing_mm_month_huc6[huc],c='r',label='outgoing'
                )
                plt.plot(
                    times_monthly,diff_incoming_outgoing_mm_month_huc6[huc],c='k',
                    label='incoming - outgoing'
                )
                plt.legend(loc='upper right')
                plt.xlabel('date')
                plt.ylabel('mm/day')
                plt.text(
                    0.01,0.99,'Average difference: {:.7f}'.format(this_avg_diff),
                    ha='left',va='top',transform=ax.transAxes
                )
                plt.savefig(this_savename)
                plt.close()
        streamflow_timeseries.columns = hucs
        streamflow_timeseries.index = times_monthly
        diff_waterwatch_cn_mm_month_hucs = (
            streamflow_timeseries - strm_mm_month_huc6
        )
        if plot_waterwatch_vs_catchment_mm_month_huc6:
            this_fname = (
                'Waterwatch_vs_Catchment_mm_month_huc6_' +
                '{exp}_{huc}_{start}_{end}.png'
            )
            for h,huc in enumerate(hucs):
                print('plotting mm/month vs truth for huc {}'.format(huc))
                this_savename = os.path.join(
                    plots_dir,
                    (
                        this_fname.format(
                            exp=exp_name,
                            huc=huc,
                            start=self.start.strftime('%Y%m%d'),
                            end=self.end.strftime('%Y%m%d')
                        )
                    )
                )
                this_avg_diff = np.mean(
                    diff_waterwatch_cn_mm_month_hucs[huc]
                )
                fig,ax = plt.subplots()
                plt.plot(times_monthly,streamflow_timeseries[huc],c='b',label='WaterWatch')
                plt.plot(
                    times_monthly,strm_mm_month_huc6[huc],c='r',label='Catchment-CN'
                )
                plt.plot(
                    times_monthly,diff_waterwatch_cn_mm_month_hucs[huc],c='k',
                    label='WaterWatch - Catchment-CN'
                )
                plt.legend(loc='upper right')
                plt.xlabel('date')
                plt.ylabel('mm/day')
                plt.text(
                    0.01,0.99,'Average difference: {:.7f}'.format(this_avg_diff),
                    ha='left',va='top',transform=ax.transAxes
                )
                plt.savefig(this_savename)
                plt.close()
    def point_forcing_vs_default(self,default_timeseries_pixel,pf_timeseries_pixel,
                                 point_forcing_dir,
                                 single_point_from_nc4,out_dir,plots_dir):
        # note: we are only going to compare for the first pixel in pf because
        # otherwise this would be an incredible amount of computation and
        # memory required. we can move on from here if we decide that we need
        # to
        # which zero-indexed pixel do we want to compare here?
        #this_pix = 2874
        default_catch_tot_precip = default_timeseries_pixel['rainfsnowf']
        pf_catch_tot_precip = pf_timeseries_pixel['rainfsnowf']
        pixels = np.array(default_catch_tot_precip.columns)
        for p,pix in enumerate(pixels):
            print('comparing pso files with model runs for {}'.format(
                pix)
            )
            this_pix = pix
            default_catch_tot_precip = default_timeseries_pixel['rainfsnowf']
            pf_catch_tot_precip = pf_timeseries_pixel['rainfsnowf']
            pix_idx = np.where(pixels == this_pix)
            pix_idx = pix_idx[0][0]
            # let's first get the point forcing
            start_year = int(self.start.strftime('%Y'))
            end_year = int(self.end.strftime('%Y'))
            if start_year != end_year:
                print('point_forcing_vs_default is only set up for one year')
                print('please choose start and end to be beginning and end of')
                print('the same year')
                print('Exiting')
                sys.exit()
            this_year = copy.deepcopy(start_year)
            while this_year <= end_year:
                this_pf_date_savename = os.path.join(
                    out_dir,'point_forcing_date_{}_{}.npy'.format(
                        this_year,this_pix
                    )
                )
                this_pf_rainf_savename = os.path.join(
                    out_dir,'point_forcing_rainf_{}_{}.npy'.format(
                        this_year,this_pix
                    )
                )
                this_pf_rainfc_savename = os.path.join(
                    out_dir,'point_forcing_rainfc_{}_{}.npy'.format(
                        this_year,this_pix
                    )
                )
                this_pf_snowf_savename = os.path.join(
                    out_dir,'point_forcing_snowf_{}_{}.npy'.format(
                        this_year,this_pix
                    )
                )
                this_pf_tile_savename = os.path.join(
                    out_dir,'point_forcing_tile_{}_{}.npy'.format(
                        this_year,this_pix
                    )
                )
                # do we need to load these from the .nc4? if so, they will be saved
                if single_point_from_nc4:
                    this_fname = point_forcing_dir.format(this_year)
                    print('getting nc.Dataset')
                    this_pf = nc.Dataset(this_fname)
                    print('extracting tile to np array')
                    this_pf_tile = np.array(this_pf['tile'][:,pix_idx])
                    this_pf_tile = [int(t) for t in this_pf_tile]
                    this_pf_tile = np.array(this_pf_tile)
                    # now save this tile info
                    np.save(this_pf_tile_savename,this_pf_tile)
                    print('extracting date to np array')
                    this_pf_date = np.array(this_pf['date_int'][:,pix_idx])
                    this_pf_date = [int(d) for d in this_pf_date]
                    this_pf_date = np.array(this_pf_date)
                    # save the dates
                    np.save(this_pf_date_savename,this_pf_date)
                    print('extracting rainf to np arary')
                    this_pf_rainf = np.array(this_pf['Rainf'][:,pix_idx])
                    np.save(this_pf_rainf_savename,this_pf_rainf)
                    print('extracting rainf_c to np array')
                    this_pf_rainfc = np.array(this_pf['Rainf_C'][:,pix_idx])
                    np.save(this_pf_rainfc_savename,this_pf_rainfc)
                    print('extracting snowf to np array')
                    this_pf_snowf = np.array(this_pf['Snowf'][:,pix_idx])
                    np.save(this_pf_snowf_savename,this_pf_snowf)
                else:
                    this_pf_tile = np.load(this_pf_tile_savename)
                    this_pf_date = np.load(this_pf_date_savename)
                    this_pf_rainf = np.load(this_pf_rainf_savename)
                    this_pf_rainfc = np.load(this_pf_rainfc_savename)
                    this_pf_snowf = np.load(this_pf_snowf_savename)
                this_year += 1
            # we need to make sure that we are looking at the same time in the
            # default and in the forcing
            this_pf_tile_check = this_pf_tile[0]
            if this_pf_tile_check != this_pix:
                print('point forcing tile is not the same as default tile')
                print('check what is going on!')
                print('exiting')
                sys.exit()
            this_pf_date_dt = [
                datetime.datetime.strptime(
                    str(d),'%Y%m%d%H%M%S'
                ) for d in this_pf_date
            ]
            this_years = np.array([d.year for d in this_pf_date_dt])
            this_months = np.array([d.month for d in this_pf_date_dt])
            this_days = np.array([d.day for d in this_pf_date_dt])
            this_pf_tot_precip = this_pf_snowf + this_pf_rainf
            #print(this_pf_date_dt)
            # concatenate to daily and add to df
            dates_day = list(default_catch_tot_precip.index)
            dates_day = [str(d) for d in dates_day]
            dates_day = [d[:10] for d in dates_day]
            dates_day_dt = [
                datetime.datetime.strptime(
                    str(d),'%Y-%m-%d'
                ) for d in dates_day
            ]
            default_catch_tot_precip = np.array(default_catch_tot_precip[this_pix])
            pf_catch_tot_precip = np.array(pf_catch_tot_precip[this_pix])
            pf_tot_precip = np.zeros(len(dates_day_dt))
            for d,dt in enumerate(dates_day_dt):
                this_day_idx = np.where(
                    (this_years == dt.year) &
                    (this_months == dt.month) &
                    (this_days == dt.day)
                )
                this_day_pf_rainf = this_pf_rainf[this_day_idx]
                this_day_pf_snowf = this_pf_snowf[this_day_idx]
                this_day_pf_rainfsnowf = (
                    this_day_pf_rainf + this_day_pf_snowf
                )
                avg_this_day_pf_rainfsnowf = np.average(
                    this_day_pf_rainfsnowf
                )
                pf_tot_precip[d] = avg_this_day_pf_rainfsnowf
            diff_pf_default_catch = (
                pf_tot_precip - default_catch_tot_precip
            )
            diff_pf_pf_catch = (
                pf_tot_precip - pf_catch_tot_precip
            )
            plot_fname = (
                os.path.join(
                    plots_dir,
                    'default_catch_vs_point_forcing_files_{}_{}.png'.format(
                        start_year,this_pix
                    )
                )
            )
            plt.figure()
            plt.plot(dates_day_dt,pf_tot_precip,label='point forcing',c='b')
            plt.plot(dates_day_dt,default_catch_tot_precip,label='default Catch',c='g')
            plt.plot(dates_day_dt,diff_pf_default_catch,label='diff pf default-Catch',c='k')
            
            plt.xlabel('kg/m2/s')
            plt.legend()
            plt.savefig(plot_fname)
            plt.close()

            plot_fname = (
                os.path.join(
                    plots_dir,
                    'pf_catch_vs_point_forcing_files_{}_{}.png'.format(
                        start_year,this_pix
                    )
                )
            )
            
            plt.figure()
            plt.plot(dates_day_dt,pf_tot_precip,label='point forcing',c='b')
            plt.plot(dates_day_dt,pf_catch_tot_precip,label='pf Catch',c='g')
            plt.plot(
                dates_day_dt,diff_pf_pf_catch,label='diff pf pf-Catch',c='brown'
            )
            plt.xlabel('kg/m2/s')
            plt.legend()
            plt.savefig(plot_fname)
            plt.close()
            
            plot_fname = (
                os.path.join(
                    plots_dir,
                    'point_forcing_files_{}_{}.png'.format(
                        start_year,this_pix
                    )
                )
            )
            
            plt.figure()
            plt.plot(dates_day_dt,pf_tot_precip,label='point forcing',c='b')
            plt.xlabel('kg/m2/s')
            plt.legend()
            plt.savefig(plot_fname)
            plt.close()
