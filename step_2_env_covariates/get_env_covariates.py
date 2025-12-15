import pandas as pd
import numpy as np
from dateutil.relativedelta import relativedelta
import netCDF4 as nc
import sys
import datetime
import xarray as xr

class get_covariates:
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
    def get_soil_info(self,fname,save_ksat,save_ksat_normal,save_sand,
                      save_sand_normal):
        '''
        Load soil texture data and returns data at selected tiles
        '''
        dat_cols = [
            'tile','k_sat','sand_perc_root'
        ]
        with open(fname) as fp:
            print('loading soil data')
            counter = 0
            for l in fp:
                counter += 1
            print(counter)
        tile = np.zeros(counter)
        k_sat = np.zeros(counter)
        sand_perc_root = np.zeros(counter)
        with open(fname) as fp:
            for l,line in enumerate(fp):
                line_txt = line.strip()
                line_list = line_txt.split()
                tile[l] = int(line_list[0])
                k_sat[l] = line_list[7]
                sand_perc_root[l] = line_list[15]
        data = np.array([
            tile,k_sat,sand_perc_root
        ])
        data = np.transpose(data)
        soil_df = pd.DataFrame(data,columns=dat_cols)
        soil_df = soil_df.set_index('tile')
        # save k_sat to nc4 to be used by pso
        ksat_array = np.array(soil_df['k_sat'].loc[self.tiles])
        ksat_df = xr.Dataset(
            data_vars=dict(
                k_sat=(['tile_num'],ksat_array)
            ),
            coords=dict(
                tile=(['tile_num'],self.tiles)
            )
        )
        print(ksat_df)
        ksat_df.to_netcdf(save_ksat)
        # lets make a standardized df
        ksat_normal_array = np.array(soil_df['k_sat'].loc[self.tiles])
        ksat_max = np.max(ksat_normal_array)
        ksat_normal_array = ksat_normal_array/ksat_max
        ksat_normal_df = xr.Dataset(
            data_vars=dict(
                k_sat_norm=(['tile_num'],ksat_normal_array)
            ),
            coords=dict(
                tile=(['tile_num'],self.tiles)
            )
        )
        print(ksat_normal_df)
        ksat_normal_df.to_netcdf(save_ksat_normal)
        # save sand fraction to nc4 to be used by pso
        sand_array = np.array(soil_df['sand_perc_root'].loc[self.tiles])
        sand_df = xr.Dataset(
            data_vars=dict(
                sand_perc=(['tile_num'],sand_array)
            ),
            coords=dict(
                tile=(['tile_num'],self.tiles)
            )
        )
        print(sand_df)
        sand_df.to_netcdf(save_sand)
        # normalized array for sand
        sand_normal_array = np.array(soil_df['sand_perc_root'].loc[self.tiles])
        sand_max = np.max(sand_normal_array)
        sand_normal_array = sand_normal_array/sand_max
        sand_normal_df = xr.Dataset(
            data_vars=dict(
                sand_perc_norm=(['tile_num'],sand_normal_array)
            ),
            coords=dict(
                tile=(['tile_num'],self.tiles)
            )
        )
        print(sand_normal_df)
        sand_normal_df.to_netcdf(save_sand_normal)
    def get_env_info(self,start,end,generic_fname,save_precip,
                     save_precip_norm,save_lai):
        # first get the static covariates needed for Ksat EF
        # this is sand fraction and ksat from only textural effects (default)

        # use start and end to get all dates
        all_dates = np.arange(
            start,end,np.timedelta64(1,'M'),dtype='datetime64[M]'
        )
        # get the years and months associated with these dates
        all_years = all_dates.astype('datetime64[Y]').astype(int) + 1970
        all_months = all_dates.astype('datetime64[M]').astype(int) % 12 + 1
        # initialize an array to hold loaded precipitation
        all_precip = np.zeros((len(all_dates),len(self.tiles)))
        all_lai = np.zeros((len(all_dates),len(self.tiles)))
        # loop over each month, since we have a precip file for each month
        for d,date in enumerate(all_dates):
            print('working for env covariates for year {}.'.format(date))
            # get this year and month
            this_year = all_years[d]
            this_month = all_months[d]
            # get next year and month. used for converting units
            next_date = date + np.timedelta64(1,'M')
            next_year = next_date.astype('datetime64[Y]').astype(int) + 1970
            next_month = next_date.astype('datetime64[M]').astype(int) % 12 + 1
            # open the dataset for this month
            this_ds = nc.Dataset(
                generic_fname.format(year=this_year,month=this_month)
            )
            #for var in this_ds.variables:
            #    print(this_ds[var].long_name)
            #sys.exit()
            # turn vals of interest into an np array
            precip = np.array(this_ds['PRECTOTLAND'])[0]
            lai = np.array(this_ds['LAI'])[0]
            # select vals for just tiles of interest
            precip = precip[self.tiles_idx]
            lai = lai[self.tiles_idx]
            # get number of secodns in this month for converting units
            dt = (
                datetime.date(next_year,next_month,1) -
                datetime.date(this_year,this_month,1)
            )
            seconds_month = dt.total_seconds()
            # convert units for precipitation to match the coefficients 
            # found offline
            precip = precip*seconds_month
            # save to all holders
            all_precip[d,:] = precip
            all_lai[d,:] = lai
        # average the vals for each tile
        avg_precip = np.mean(all_precip,axis=0)
        avg_lai = np.mean(all_lai,axis=0)
        # save to an nc4 file so it can be opened by the pso in Catchment
        # first for precip
        precip_df = xr.Dataset(
            data_vars=dict(
                precip=(['tile_num'],avg_precip)
            ),
            coords=dict(
                tile=(['tile_num'],self.tiles)
            )
        )
        print(precip_df)
        precip_df.to_netcdf(save_precip)
        # and precip normalized
        norm_precip = np.copy(avg_precip)
        max_precip = np.max(norm_precip)
        norm_precip = norm_precip/max_precip
        norm_precip_df = xr.Dataset(
            data_vars=dict(
                norm_precip=(['tile_num'],norm_precip)
            ),
            coords=dict(
                tile=(['tile_num'],self.tiles)
            )
        )
        print(norm_precip_df)
        norm_precip_df.to_netcdf(save_norm)
        # and now for lai
        lai_df = xr.Dataset(
            data_vars=dict(
                lai=(['tile_num'],avg_lai)
            ),
            coords=dict(
                tile=(['tile_num'],self.tiles)
            )
        )
        print(lai_df)
        lai_df.to_netcdf(save_lai)
    def get_gldas_info(self,start,end,fname,save_name,save_name_norm):
        this_data = pd.read_csv(fname)
        this_data = this_data.set_index('time')
        tiles_data = np.array(this_data.columns)
        tiles_times = np.array(this_data.index)
        start = start.astype(datetime.date)
        start_str = datetime.datetime.strftime(start,'%Y-%m-%d')
        start_idx = np.where(
            tiles_times == start_str
        )[0][0]
        end = end.astype(datetime.date)
        if end.day != 1:
            end = datetime.date(end.year,end.month,1)
        end_str = datetime.datetime.strftime(end,'%Y-%m-%d')
        end_idx = np.where(
            tiles_times == end_str
        )[0][0]
        this_data = this_data.iloc[start_idx:end_idx+1]
        averaged_data = pd.DataFrame(columns=self.tiles)
        averaged_data.loc['average'] = np.zeros(
            len(self.tiles)
        )
        for t,ti in enumerate(self.tiles):
            this_tile_info = np.array(this_data[str(ti)])
            this_avg = np.average(this_tile_info)
            averaged_data[ti].loc['average'] = this_avg
        var_array = np.array(averaged_data[self.tiles].loc['average'])
        var_ds = xr.Dataset(
            data_vars=dict(
                vals=(['tile_num'],var_array)
            ),
            coords=dict(
                tile=(['tile_num'],self.tiles)
            )
        )
        var_ds.to_netcdf(save_name)
        norm_var_array = np.copy(var_array)
        var_mean = np.mean(norm_var_array)
        var_std = np.std(norm_var_array)
        norm_var_array = norm_var_array - var_mean
        norm_var_array = norm_var_array/var_std
        norm_var_ds = xr.Dataset(
            data_vars=dict(
                vals=(['tile_num'],norm_var_array)
            ),
            coords=dict(
                tile=(['tile_num'],self.tiles)
            )
        )
        print(norm_var_ds)
        norm_var_ds.to_netcdf(save_name_norm)
    def get_canopy_info(self,fname,save_name,save_name_norm):
        canopy_data = pd.read_csv(fname)
        canopy_data = canopy_data.set_index('time')
        tiles_str = [
            str(t) for t in self.tiles
        ]
        val_array = np.array(canopy_data[tiles_str].loc['static'])
        canopy_ds = xr.Dataset(
            data_vars=dict(
                vals=(['tile_num'],val_array)
            ),
            coords=dict(
                tile=(['tile_num'],self.tiles)
            )
        )
        print(canopy_ds)
        canopy_ds.to_netcdf(save_name)
        # normalize and save
        norm_var_array = np.copy(val_array)
        var_mean = np.mean(norm_var_array)
        var_std = np.std(norm_var_array)
        norm_var_array = norm_var_array - var_mean
        norm_var_array = norm_var_array/var_std
        norm_canopy_ds = xr.Dataset(
            data_vars=dict(
                vals=(['tile_num'],norm_var_array)
            ),
            coords=dict(
                tile=(['tile_num'],self.tiles)
            )
        )
        print(norm_canopy_ds)
        norm_canopy_ds.to_netcdf(save_name_norm)



