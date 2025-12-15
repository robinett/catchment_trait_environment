import netCDF4 as nc
import pickle as pkl
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import datetime

class compare:
    def get_dates(self,lai_fnames):
        fnames = list(np.genfromtxt(lai_fnames,dtype='str'))
        dates = []
        for f,fn in enumerate(fnames):
            this_stng = str(fn[-11:-4])
            dates.append(this_stng)
        return dates
    def get_lai_truth(self,lai_fname):
        lai_ds = nc.Dataset(lai_fname)
        lai_vals = np.array(lai_ds['lai'])
        nans_idx = np.where(lai_vals > 100)
        lai_vals[nans_idx] = np.nan
        lai_vals = lai_vals*0.1
        return lai_vals
    def get_lai_catch(self,catch_dir,start,end):
        with open(catch_dir,'rb') as f:
            catch_dict = pkl.load(f)
        print(catch_dict)
        catch_lai = catch_dict['laiv']
        print(catch_lai)
        #trim to start and end dates
        #catch_dates = list(catch_lai.index)
        #catch_dates_dt = []
        #for d in catch_dates:
        #    ts = (
        #        (d - np.datetime64('1970-01-01T00:00:00Z'))/
        #        np.timedelta64(1,'s')
        #    )
        #    catch_dates_dt.append(datetime.datetime.utcfromtimestamp(ts))
        #start_idx_all = np.where(np.array(catch_dates_dt) < start)[0][0]
        #start_idx = start_idx_all[-1] + 1
        #catch_lai = catch_lai.iloc[start_idx:]
        print(catch_lai)
        return catch_lai
    def get_lai_lons_lats(self,lons_fname,lats_fname):
        lons = np.genfromtxt(lons_fname)
        lats = np.genfromtxt(lats_fname)
        return [lons,lats]
    def get_tiles_and_indices(self,tiles_fname,tile_info_fname):
        # get the tiles that we are going to run as defined in step 1
        tiles = pd.read_csv(tiles_fname,header=None)
        # turn this into an np array
        tiles = np.array(tiles).astype(int)
        # make a nice np array
        tiles = tiles.T
        tiles = tiles[0]
        tile_info = pd.read_csv(tile_info_fname)
        tiles_i = np.zeros(len(tiles))
        tiles_j = np.zeros(len(tiles))
        for t,ti in enumerate(tiles):
            this_idx = np.where(tile_info['tile_id'] == ti)[0][0]
            tiles_i[t] = tile_info['i_indg'].iloc[this_idx]
            tiles_j[t] = tile_info['j_indg'].iloc[this_idx]
        return [tiles,tiles_i,tiles_j,tile_info]
    def plot_modis_lai(self,cn_lons,cn_lats,vals,dates,plots_dir):
        cmap = 'winter'
        vmin = 0
        vmax = 10
        for p in range(np.shape(vals)[2]):
            print('plotting modis lai for {}'.format(dates[p]))
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
            val = vals[:,:,p]
            val = np.ndarray.flatten(val)
            # make the scatterplot
            scatter = plt.scatter(
                cn_lons,cn_lats,marker='s',s=1,c=val,transform=ccrs.PlateCarree(),
                cmap=cmap,vmin=vmin,vmax=vmax
            )
            # add colorbar
            plt.colorbar(scatter)
            # save
            savename = 'lai_conus_nans_{}.png'.format(
                dates[p]
            )
            savename = os.path.join(
                plots_dir,savename
            )
            plt.savefig(savename,dpi=300,bbox_inches='tight')
            plt.close()
    def get_comp_df(self,modis_vals,catch_vals,lons,lats,tile_info,tiles,dates):
        # get cn dates
        cn_dates = catch_vals.index
        comp_df = pd.DataFrame(index=cn_dates)
        # find where cn dates are the same as modis dates
        dates_idx = np.zeros(len(dates),dtype='int')
        tile_lons = np.zeros(0)
        tile_lats = np.zeros(0)
        for d,dt in enumerate(dates):
            dt_dt = datetime.datetime.strptime(dt,'%Y%j')
            dt_str = dt_dt.strftime('%Y-%m-%d')
            this_idx = np.where(cn_dates == dt_str)[0][0]
            dates_idx[d] = int(this_idx)
        for t,ti in enumerate(tiles):
            this_lat = tile_info['com_lat'].iloc[ti-1]
            this_lon = tile_info['com_lon'].iloc[ti-1]
            tile_lats = np.append(tile_lats,this_lat)
            tile_lons = np.append(tile_lons,this_lon)
            lat_idx = np.where(lats == this_lat)[0][0]
            lon_idx = np.where(lons == this_lon)[0][0]
            this_modis = modis_vals[lat_idx,lon_idx,:]
            this_cn = catch_vals[ti]
            this_modis_fill = np.repeat(np.nan,len(this_cn))
            this_modis_fill[dates_idx] = this_modis
            comp_df['{}_modis'.format(ti)] = this_modis_fill
            comp_df['{}_cn'.format(ti)] = this_cn
            comp_df = comp_df.copy()
        return [comp_df,tile_lats,tile_lons]
    def plot_comparison_timeseries(self,comp_df,tiles,plots_dir):
        dates = list(comp_df.index)
        for t,ti in enumerate(tiles):
            print('plotting timeseries for {}'.format(ti))
            this_modis = comp_df['{}_modis'.format(ti)]
            this_cn = comp_df['{}_cn'.format(ti)]
            plt.figure()
            plt.plot(dates,this_cn,label='Catch-CN LAI')
            plt.scatter(dates,this_modis,label='MODIS LAI')
            plt.legend()
            savename = os.path.join(
                plots_dir,
                'lai_timeseries_{}.png'.format(ti)
            )
            plt.savefig(savename)
            plt.close()
    def plot_comparison_map(self,comp_df,tiles,cn_lons,cn_lats,plots_dir):
        modis_avg_vals = np.zeros(len(tiles))
        cn_avg_vals = np.zeros(len(tiles))
        for t,ti in enumerate(tiles):
            this_modis = comp_df['{}_modis'.format(ti)]
            this_cn = comp_df['{}_cn'.format(ti)]
            this_avg_modis = np.nanmean(this_modis)
            this_avg_cn = np.nanmean(this_cn)
            modis_avg_vals[t] = this_avg_modis
            cn_avg_vals[t] = this_avg_cn
        lai_diff = cn_avg_vals - modis_avg_vals
        names = [
            'cn_avg_lai',
            'modis_avg_lai',
            'lai_diff_cn_modis'
        ]
        vals = [
            cn_avg_vals,
            modis_avg_vals,
            lai_diff
        ]
        cmaps = [
            'winter',
            'winter',
            'bwr'
        ]
        vmins = [
            0,
            0,
            -5
        ]
        vmaxs = [
            10,
            10,
            5
        ]
        for n in range(len(names)):
            print('plotting avg lai map for {}'.format(names[n]))
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
            val = vals[n]
            # make the scatterplot
            scatter = plt.scatter(
                cn_lons,cn_lats,marker='s',s=1,c=val,transform=ccrs.PlateCarree(),
                cmap=cmaps[n],vmin=vmins[n],vmax=vmaxs[n]
            )
            # add colorbar
            plt.colorbar(scatter)
            # save
            savename = '{}.png'.format(
                names[n]
            )
            savename = os.path.join(
                plots_dir,savename
            )
            plt.savefig(savename,dpi=300,bbox_inches='tight')
            plt.close()
