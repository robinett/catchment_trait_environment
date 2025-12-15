import numpy as np
import os
import netCDF4 as nc
import sys
import numpy as np
import copy
import datetime
import pickle as pkl
import pandas as pd

class get_timeseries:
    def __init__(self,start,end):
        self.start = start
        self.end = end
        self.output_middle_dir = (
            'output/SMAP_EASEv2_M36/cat/ens0000'
        )
    def get_pixels(self,step_1_pixels_fname):
        # get the tiles that we are going to run as defined in step 1
        tiles = pd.read_csv(step_1_pixels_fname,header=None)
        # turn this into an np array
        tiles = np.array(tiles).astype(int)
        # make a nice np array
        tiles = tiles.T
        tiles = tiles[0]
        # save to self as pixels for later use
        self.pixels = tiles
        return tiles
    def print_catch_variables(self,exp_dir):
        '''
        Function that will print variables in the Catchment dataset an then
        exit. Prints all varaibles with full description, then just a list of
        all the variables names.
        '''
        exp_results_dir = os.path.join(
            exp_dir,self.output_middle_dir
        )
        exp_tech_name = os.path.basename(os.path.normpath(exp_dir))
        self.exp_tech_name = exp_tech_name
        curr_str = self.start.strftime("%Y%m%d")
        curr_yr = self.start.strftime("%Y")
        curr_mon = self.start.strftime("%m")
        this_file = os.path.join(
            exp_results_dir,
            'Y'+curr_yr,
            'M'+curr_mon,
            exp_tech_name+'.tavg24_1d_lnd_Nt.'+curr_str+'_1200z.nc4'
        )
        ds_pix = nc.Dataset(this_file)
        print('forcing variables in this catchment dataset:')
        print(ds_pix.variables)
        for var in ds_pix.variables:
            print(var)
        this_file = os.path.join(
            exp_results_dir,
            'Y'+curr_yr,
            'M'+curr_mon,
            exp_tech_name+'.tavg24_1d_lfs_Nt.'+curr_str+'_1200z.nc4'
        )
        ds_pix = nc.Dataset(this_file)
        print('variables in this catchment dataset:')
        print(ds_pix.variables)
        for var in ds_pix.variables:
            print(var)
        sys.exit()
    def get_catch_timeseries(self,exp_dir,exp_name,save_bool,
                             load_bool,is_default_experiment,read_met_forcing,
                             timeseries_dir):
        '''
        Function that creates a dictionary of dictionries of compiled model output. The keys to the
        first layer of the dictionary are the model pixel numbers, and the keys to the second
        layer of the dictionaries are the variables of interest.
        Inputs:
            All described in "main.py"
        Outputs:
            out_dir: as described in the function description above.
        '''
        exp_results_dir = os.path.join(
            exp_dir,'output','SMAP_EASEv2_M36',
            'cat','ens0000'
        )
        # get the true pixel numbers
        # needed no matter if we load or not
        true_tiles = np.zeros(0)
        # load the catchment tile numbers run in the experiment from the
        # ldas_domain.txt file
        exp_tech_name = os.path.basename(os.path.normpath(exp_dir))
        if is_default_experiment:
            ldas_domain_fname = os.path.join(
                exp_dir,'output/SMAP_EASEv2_M36/rc_out',
                exp_tech_name+'.ldas_domain.txt'
            )
        else:
            ldas_domain_fname = os.path.join(
                (
                    '/lustre/catchment/exps/' +
                    'GEOSldas_CN45_med_default_pft_g1_2006_2007_camels/0/output/'+
                    'SMAP_EASEv2_M36/rc_out'
                ),'0.ldas_domain.txt'
            )
        with open(ldas_domain_fname) as fp:
            for line in fp:
                line_txt = line.strip()
                true_tiles = np.append(true_tiles,int(line_txt[0:8]))
        true_tiles = true_tiles[1:]
        # dummy holder that will become fake tiles
        catch_pix_idx = np.zeros(len(true_tiles),dtype='int32')
        # for all step 1 pixels, get idx of those in context of all 
        # tiles from model
        for p,pi in enumerate(self.pixels):
            try:
                idx = np.where(true_tiles == pi)[0][0]
                catch_pix_idx[p] = idx
            except:
                raise Exception('User-specified pixel {} does not exist in model output.'.format(pi))
        #self.pixels = pixels
        self.catch_pix_idx = catch_pix_idx
        # if datset hasn't been saved then generate it
        if load_bool == False:
            # set the number of pixels that we will be working with and number
            # of days over which we are running
            num_pixels = len(self.pixels)
            num_days = (self.end - self.start).days + 1
            # for trns
            all_times_trns = np.zeros((num_days,num_pixels))
            # for ave_sm
            all_times_ave_sm = np.zeros((num_days,num_pixels))
            # for root_sm
            all_times_root_sm = np.zeros((num_days,num_pixels))
            # for surf_sm
            all_times_surf_sm = np.zeros((num_days,num_pixels))
            # for beta
            all_times_beta = np.zeros((num_days,num_pixels))
            # for infil
            all_times_infil = np.zeros((num_days,num_pixels))
            # for le
            all_times_le = np.zeros((num_days,num_pixels))
            # for runoff
            all_times_runoff = np.zeros((num_days,num_pixels))
            # for baseflow
            all_times_baseflow = np.zeros((num_days,num_pixels))
            # for visible lai
            all_times_laiv = np.zeros((num_days,num_pixels))
            # for total lai
            all_times_lait = np.zeros((num_days,num_pixels))
            # for lat
            all_times_lat = np.zeros((num_days,num_pixels))
            # for lon
            all_times_lon = np.zeros((num_days,num_pixels))
            # if instructed to get met forcing information, let's do this as
            # well
            if read_met_forcing:
                # for tair
                all_times_tair = np.zeros((num_days,num_pixels))
                # for qair
                all_times_qair = np.zeros((num_days,num_pixels))
                # for lwdown
                all_times_lwdown = np.zeros((num_days,num_pixels))
                # for swdown
                all_times_swdown = np.zeros((num_days,num_pixels))
                # for psurf
                all_times_psurf = np.zeros((num_days,num_pixels))
                # for rainf
                all_times_rainf = np.zeros((num_days,num_pixels))
                # for rainfsnowf
                all_times_rainfsnowf = np.zeros((num_days,num_pixels))
            # increment by one day for each loop
            delta = datetime.timedelta(days=1)
            # set current as the start date
            curr = self.start
            #set the varaibles that we are going to use to keep
            #track of the different time series
            #model generate varaibles
            time = np.arange(
                self.start,self.end+datetime.timedelta(days=1),
                datetime.timedelta(days=1)
            )
            time = pd.to_datetime(time)
            d = 0
            # for each day
            while (curr <= self.end):
                # get this_time
                this_time = time[d]
                if curr.day == 1:
                    # let the user know where we are at
                    print('working on directory {}'.format(
                        exp_dir)
                    )
                    print('working on year {}'.format(
                        curr.year)
                    )
                    print('working on month {}'.format(
                        curr.month)
                    )
                # get the current info
                curr_str = curr.strftime("%Y%m%d")
                curr_yr = curr.strftime("%Y")
                curr_mon = curr.strftime("%m")
                # open the results and forcing files
                # again, make sure formatted correctly 
                if is_default_experiment:
                    this_file = os.path.join(
                        exp_results_dir,
                        'Y'+curr_yr,
                        'M'+curr_mon,
                        exp_tech_name+'.tavg24_1d_lnd_Nt.'+curr_str+'_1200z.nc4'
                    )
                    # turn the files into a dataset
                    ds = nc.Dataset(this_file,mode='r')
                else:
                    this_file = os.path.join(
                        exp_dir,
                        'Y'+curr_yr,
                        'M'+curr_mon,
                        exp_tech_name+'.tavg24_1d_lnd_Nt.'+curr_str+'_1200z.nc4'
                    )
                    # turn the files into a dataset
                    ds = nc.Dataset(this_file,mode='r')
                # if we want to read met forcing, get that file
                if read_met_forcing and is_default_experiment:
                    this_force = os.path.join(
                        exp_results_dir,
                        'Y'+curr_yr,
                        'M'+curr_mon,
                        exp_tech_name+'.tavg24_1d_lfs_Nt.'+curr_str+'_1200z.nc4'
                    )
                    ds_force = nc.Dataset(this_force,mode='r')
                elif read_met_forcing and not is_default_experiment:
                    this_force = os.path.join(
                        exp_dir,
                        'Y'+curr_yr,
                        'M'+curr_mon,
                        exp_tech_name+'.tavg24_1d_lfs_Nt.'+curr_str+'_1200z.nc4'
                    )
                    ds_force = nc.Dataset(this_force,mode='r')
                # for each file, loop over each pixel for that day and get the
                # information
                all_trns = np.array(ds['EVPTRNS'][:])[0]
                this_trns = all_trns[catch_pix_idx]
                all_times_trns[d,:] = this_trns
                all_ave_sm = np.array(ds['PRMC'][:])[0]
                this_ave_sm = all_ave_sm[catch_pix_idx]
                all_times_ave_sm[d,:] = this_ave_sm
                all_root_sm = np.array(ds['RZMC'][:])[0]
                this_root_sm = all_root_sm[catch_pix_idx]
                all_times_root_sm[d,:] = this_root_sm
                all_surf_sm = np.array(ds['SFMC'][:])[0]
                this_surf_sm = all_surf_sm[catch_pix_idx]
                all_times_surf_sm[d,:] = this_surf_sm
                all_le = np.array(ds['LHLAND'][:])[0]
                this_le = all_le[catch_pix_idx]
                all_times_le[d,:] = this_le
                all_beta = np.array(ds['BTRAN'][:])[0]
                this_beta = all_beta[catch_pix_idx]
                all_times_beta[d,:] = this_beta
                all_infil = np.array(ds['QINFIL'][:])[0]
                this_infil = all_infil[catch_pix_idx]
                all_times_infil[d,:] = this_infil
                all_runoff = np.array(ds['RUNOFF'][:])[0]
                this_runoff = all_runoff[catch_pix_idx]
                all_times_runoff[d,:] = this_runoff
                all_baseflow = np.array(ds['BASEFLOW'][:])[0]
                this_baseflow = all_baseflow[catch_pix_idx]
                all_times_baseflow[d,:] = this_baseflow
                all_laiv = np.array(ds['CNLAI'][:])[0]
                this_laiv = all_laiv[catch_pix_idx]
                all_times_laiv[d,:] = this_laiv
                all_lait = np.array(ds['CNTLAI'][:])[0]
                this_lait = all_lait[catch_pix_idx]
                all_times_lait[d,:] = this_lait
                all_lat = np.array(ds['lat'][:])
                this_lat = all_lat[catch_pix_idx]
                all_times_lat[d,:] = this_lat
                all_lon = np.array(ds['lon'][:])
                this_lon = all_lon[catch_pix_idx]
                all_times_lon[d,:] = this_lon
                if read_met_forcing:
                    all_tair = np.array(ds_force['Tair'][:])[0]
                    this_tair = all_tair[catch_pix_idx]
                    all_times_tair[d,:] = this_tair
                    all_qair = np.array(ds_force['Qair'][:])[0]
                    this_qair = all_qair[catch_pix_idx]
                    all_times_qair[d,:] = this_qair
                    all_lwdown = np.array(ds_force['LWdown'][:])[0]
                    this_lwdown = all_lwdown[catch_pix_idx]
                    all_times_lwdown[d,:] = this_lwdown
                    all_swdown = np.array(ds_force['SWdown'][:])[0]
                    this_swdown = all_swdown[catch_pix_idx]
                    all_times_swdown[d,:] = this_swdown
                    all_psurf = np.array(ds_force['Psurf'][:])[0]
                    this_psurf = all_psurf[catch_pix_idx]
                    all_times_psurf[d,:] = this_psurf
                    all_rainf = np.array(ds_force['Rainf'][:])[0]
                    this_rainf = all_rainf[catch_pix_idx]
                    all_times_rainf[d,:] = this_rainf
                    all_rainfsnowf = np.array(ds_force['RainfSnowf'][:])[0]
                    this_rainfsnowf = all_rainfsnowf[catch_pix_idx]
                    all_times_rainfsnowf[d,:] = this_rainfsnowf
                # update the day and return to the top of the loop
                curr += delta
                d += 1
            # if the model thinks there is negative le, just turn it to
            # zero
            all_times_le = np.where(all_times_le<0,0,all_times_le)
            # add each of the variables to out_dict
            # let's make these np arrays into dataframes
            trns_df = pd.DataFrame(
                data=all_times_trns,
                index=time,
                columns=list(self.pixels)
            )
            trns_df.index.name = 'time'
            ave_sm_df = pd.DataFrame(
                data=all_times_ave_sm,
                index=time,
                columns=list(self.pixels)
            )
            ave_sm_df.index.name = 'time'
            root_sm_df = pd.DataFrame(
                data=all_times_root_sm,
                index=time,
                columns=list(self.pixels)
            )
            root_sm_df.index.name = 'time'
            surf_sm_df = pd.DataFrame(
                data=all_times_surf_sm,
                index=time,
                columns=list(self.pixels)
            )
            surf_sm_df.index.name = 'time'
            le_df = pd.DataFrame(
                data=all_times_le,
                index=time,
                columns=list(self.pixels)
            )
            le_df.index.name = 'time'
            beta_df = pd.DataFrame(
                data=all_times_beta,
                index=time,
                columns=list(self.pixels)
            )
            beta_df.index.name = 'time'
            beta_df = pd.DataFrame(
                data=all_times_beta,
                index=time,
                columns=list(self.pixels)
            )
            beta_df.index.name = 'time'
            infil_df = pd.DataFrame(
                data=all_times_infil,
                index=time,
                columns=list(self.pixels)
            )
            infil_df.index.name = 'time'
            runoff_df = pd.DataFrame(
                data=all_times_runoff,
                index=time,
                columns=list(self.pixels)
            )
            runoff_df.index.name = 'time'
            baseflow_df = pd.DataFrame(
                data=all_times_baseflow,
                index=time,
                columns=list(self.pixels)
            )
            baseflow_df.index.name = 'time'
            laiv_df = pd.DataFrame(
                data=all_times_laiv,
                index=time,
                columns=list(self.pixels)
            )
            laiv_df.index.name = 'time'
            lait_df = pd.DataFrame(
                data=all_times_lait,
                index=time,
                columns=list(self.pixels)
            )
            lait_df.index.name = 'time'
            lat_df = pd.DataFrame(
                data=all_times_lat,
                index=time,
                columns=list(self.pixels)
            )
            lat_df.index.name = 'time'
            lon_df = pd.DataFrame(
                data=all_times_lon,
                index=time,
                columns=list(self.pixels)
            )
            lon_df.index.name = 'time'
            if read_met_forcing:
                tair_df = pd.DataFrame(
                    data=all_times_tair,
                    index=time,
                    columns=list(self.pixels)
                )
                tair_df.index.name = 'time'
                qair_df = pd.DataFrame(
                    data=all_times_qair,
                    index=time,
                    columns=list(self.pixels)
                )
                qair_df.index.name = 'time'
                lwdown_df = pd.DataFrame(
                    data=all_times_lwdown,
                    index=time,
                    columns=list(self.pixels)
                )
                lwdown_df.index.name = 'time'
                swdown_df = pd.DataFrame(
                    data=all_times_swdown,
                    index=time,
                    columns=list(self.pixels)
                )
                swdown_df.index.name = 'time'
                psurf_df = pd.DataFrame(
                    data=all_times_psurf,
                    index=time,
                    columns=list(self.pixels)
                )
                psurf_df.index.name = 'time'
                rainf_df = pd.DataFrame(
                    data=all_times_rainf,
                    index=time,
                    columns=list(self.pixels)
                )
                rainf_df.index.name = 'time'
                rainfsnowf_df = pd.DataFrame(
                    data=all_times_rainfsnowf,
                    index=time,
                    columns=list(self.pixels)
                )
                rainfsnowf_df.index.name = 'time'


            out_dict = {
                'trns':trns_df,
                'ave_sm':ave_sm_df,
                'root_sm':root_sm_df,
                'surf_sm':surf_sm_df,
                'le':le_df,
                'beta':beta_df,
                'infil':infil_df,
                'runoff':runoff_df,
                'baseflow':baseflow_df,
                'laiv':laiv_df,
                'lait':lait_df,
                'lat':lat_df,
                'lon':lon_df
            }
            if read_met_forcing:
                out_dict['tair'] = tair_df
                out_dict['qair'] = qair_df
                out_dict['lwdown'] = lwdown_df
                out_dict['swdown'] = swdown_df
                out_dict['psurf'] = psurf_df
                out_dict['rainf'] = rainf_df
                out_dict['rainfsnowf'] = rainfsnowf_df
            # let's create 
            if save_bool == True:
                save_dir = os.path.join(
                    timeseries_dir,exp_name+'.pkl'
                )
                with open(save_dir,'wb') as f:
                    pkl.dump(out_dict,f)
                print('saved experiment to {}.'.format(save_dir))
        #otherwise load it
        else:
            # get where it was saved
            load_dir = os.path.join(
                timeseries_dir,exp_name+'.pkl'
            )
            print('Loading experiment from {}'.format(load_dir))
            # open the file
            with open(load_dir,'rb') as f:
                out_dict = pkl.load(f)
            # let's only take the times that we need
            times = list(out_dict['le'].index)
            times = pd.to_datetime(times)
            # let's see where is our start and end idx
            for t,ti in enumerate(times):
                if ti.date() == self.start:
                    start_idx = t
                if ti.date() == self.end:
                    end_idx = t
            for ke in out_dict.keys():
                out_dict[ke] = out_dict[ke].iloc[start_idx:end_idx+1]
        return out_dict
    def get_fluxcom_timeseries(self,fluxcom_dir,start_err,end_err):
        # load the fluxcom data from previous step
        fluxcom_data = pd.read_csv(fluxcom_dir)
        # set the start and end dates
        start_format = start_err.strftime('%Y-%m-%d')
        end_format = end_err.strftime('%Y-%m-%d')
        # find which rows this is in fluxcom
        fluxcom_time = np.array(fluxcom_data['time'])
        start_idx = np.where(fluxcom_time == start_format)[0][0]
        end_idx = np.where(fluxcom_time == end_format)[0][0]
        # trim to these start and end dates
        fluxcom_trimmed = fluxcom_data.iloc[start_idx:end_idx+1]
        # set the fluxcom columns to be ints instead of strings, aligning with
        # outputs from catchment-CN
        cols = np.array(fluxcom_trimmed.columns)
        cols_no_time = cols[1:]
        new_cols = copy.deepcopy(cols)
        cols_int = cols_no_time.astype(int)
        new_cols[1:] = cols_int
        fluxcom_trimmed.columns = new_cols
        return fluxcom_trimmed
    def get_strm_timeseries(self,strm_dir):
        # load streamflow from previous step
        strm_data = pd.read_csv(strm_dir)
        strm_data = strm_data.set_index('time')
        # trim using start and end dates
        # format the start and end dates
        start_fmt = self.start.strftime('%Y%m%d')
        end_fmt = self.end.strftime('%Y%m%d')
        strm_data = strm_data.loc[start_fmt:end_fmt]
        return strm_data
