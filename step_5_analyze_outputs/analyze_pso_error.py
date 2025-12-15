import sys
sys.path.append('/shared/pso/other_analyses/pft_prevalence')
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import os
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import stats
import pandas as pd
from statsmodels.formula.api import ols
import copy
from get_pft_distribution import pft_dist

class analyze_pso:
    def __init__(self,optimization_type):
        self.optimization_type = optimization_type
    def get_all_info(self,all_info_fname):
        with open(all_info_fname,'rb') as f:
            all_info = pkl.load(f)
        return all_info
    def get_parameter_names(self):
        if self.optimization_type == 'pft':
            parameter_names = [
                'a0_needleleaf_trees', # 1
                'a0_broadleaf_trees', # 2
                'a0_shrub', # 3
                'a0_c3_grass', # 4
                'a0_c4_grass', # 5
                'a0_crop' # 6
            ]
        elif self.optimization_type == 'ef':
            parameter_names = [
                'b_needleleaf_trees', # 1
                'b_broadleaf_trees', # 2
                'b_shrub', # 3
                'b_c3_grass', # 4
                'b_c4_grass', # 5
                'b_crop', # 6
                'a0', # 7
                'a1', # 8
                'a2' # 9
            ]
        return parameter_names
    def plot_parameter_convergence(self,all_info,parameter_names,plots_dir,
                                   exp_names):
        this_exp_name = exp_names[2]
        # get the iterations that the pso ran
        iteration_str = list(all_info.keys())
        #num_iterations = len(iteration_str)
        iterations_considered = 5
        num_iterations = iterations_considered
        iterations = np.zeros(iterations_considered)
        for t in range(iterations_considered):
            this_str = iteration_str[t]
            for c,char in enumerate(this_str):
                if char == '_':
                    this_iter = int(this_str[c+1:])
                    iterations[t] = this_iter
        # for each of these iterations, get the parameter values
        ex_params = all_info['iteration_1']['positions']
        num_particles,num_params = np.shape(ex_params)
        all_particles = np.arange(num_particles)
        all_positions = np.zeros((num_iterations,num_particles,num_params))
        all_obj = np.zeros((num_iterations,num_particles))
        for i in range(iterations_considered):
            it = iteration_str[i]
            all_positions[i,:,:] = all_info[it]['positions']
            all_obj[i,:] = all_info[it]['obj_out_norm']
        # make a plot of the progression of each parameter over time
        for i,param in enumerate(parameter_names):
            plt.figure()
            for j,part in enumerate(all_particles):
                this_param_this_part = all_positions[:,j,i]
                plt.plot(
                    iterations,this_param_this_part,
                    label='Particle {}'.format(part)
                )
            #plt.legend()
            this_savename = os.path.join(
                plots_dir,
                'pso_{param}_movement_{exp}.png'.format(
                    param=param,exp=this_exp_name)
            )
            plt.savefig(this_savename)
            plt.close()
        # make a plot of the progression of objective values over time
        plt.figure()
        for p,part in enumerate(all_particles):
            this_part_obj = all_obj[:,p]
            plt.plot(
                iterations,this_part_obj,
                label='Particle {}'.format(part)
            )
        #plt.legend()
        this_savename = os.path.join(
            plots_dir,
            'pso_objective_value_movement_{exp}.png'.format(
                exp=this_exp_name
            )
        )
        plt.savefig(this_savename)
        plt.close()
    def plot_parameter_maps(self,all_info,tiles,precip_fname,lai_fname,sand_fname,
                            k0_fname,pft_info,default_df,plots_dir,exp_names,
                            gldas_precip_dir,gldas_temp_dir,gldas_pet_dir,
                            gldas_et_dir,canopy_height_dir,returned_rmse_dfs,
                            canopy_map_fname,intersection_info,
                            plot_pso_scatter,
                            all_info_compare='none'):
        this_exp_name = exp_names[2]
        # get env covariates maps
        precip_ds = nc.Dataset(precip_fname)
        norm_precip_vals = np.array(precip_ds['vals'])
        canopy_ds = nc.Dataset(canopy_map_fname)
        norm_canopy_vals = np.array(canopy_ds['vals'])
        #lai_ds = nc.Dataset(lai_fname)
        #lai_vals = np.array(lai_ds['lai'])
        #sand_ds = nc.Dataset(sand_fname)
        #sand_vals = np.array(sand_ds['sand_perc'])
        k0_ds = nc.Dataset(k0_fname)
        k0_vals = np.array(k0_ds['k_sat'])
        gldas_precip_df = pd.read_csv(gldas_precip_dir)
        gldas_precip_df = gldas_precip_df.set_index('time')
        #gldas_temp_df = pd.read_csv(gldas_temp_dir)
        #gldas_temp_df = gldas_temp_df.set_index('time')
        gldas_pet_df = pd.read_csv(gldas_pet_dir)
        gldas_pet_df = gldas_pet_df.set_index('time')
        #gldas_et_df = pd.read_csv(gldas_et_dir)
        #gldas_et_df = gldas_et_df.set_index('time')
        canopy_height_df = pd.read_csv(canopy_height_dir)
        canopy_height_df = canopy_height_df.set_index('time')
        # testing stuff
        #print(precip_fname)
        #print(norm_precip_vals)
        #print(gldas_precip_dir)
        #print(gldas_precip_df)
        #print(canopy_map_fname)
        #print(norm_canopy_vals)
        #print(canopy_height_dir)
        #print(canopy_height_df)
        #sys.exit()
        # let's make a scatter plot of precip values versus default model RMSE
        # to see if they co-vary
        # let's get the model error varaibles
        model_rmse = np.array(default_df.loc['le_rmse'])
        model_rmse = model_rmse[:-1]
        model_bias = np.array(default_df.loc['ave_le_diff'])
        model_bias = model_bias[:-1]
        model_r2 = np.array(default_df.loc['le_r2'])
        model_r2 = model_r2[:-1]
        model_corr = np.array(default_df.loc['le_corr'])
        model_corr = model_corr[:-1]
        model_ubrmse = np.array(default_df.loc['le_ubrmse'])
        model_ubrmse = model_ubrmse[:-1]
        model_bias = np.array(default_df.loc['ave_le_diff'])
        model_bias = model_bias[:-1]
        not_nan_idx = np.where(
            np.isnan(model_rmse) == False
        )
        nan_idx = np.where(
            np.isnan(model_rmse) == True
        )
        # let's get the covariates. start with precip
        gldas_precip_avg = np.array(gldas_precip_df.mean())
        #gldas_temp_avg = np.array(gldas_temp_df.mean())
        gldas_pet_avg = np.array(gldas_pet_df.mean())
        #gldas_et_avg = np.array(gldas_et_df.mean())
        canopy_height_avg = np.array(canopy_height_df.mean())
        pixels_str = gldas_precip_df.columns
        pixels = [
            int(p) for p in pixels_str
        ]
        pixels = np.array(pixels)
        nan_pix = pixels[nan_idx]
        not_nan_pix = pixels[not_nan_idx]
        le_obs = np.array(default_df.loc['le_obs'])
        le_obs = le_obs[not_nan_idx]
        # let's convert precip to mm/day
        gldas_precip_avg = gldas_precip_avg*86400
        # let's convert ET to mm/day
        gldas_pet_avg = gldas_pet_avg/28.94
        #gldas_et_avg = gldas_et_avg/28.94
        # let's get the streamflow errors
        default_strm = returned_rmse_dfs[0]
        strm_rmse = np.array(default_strm.loc['strm_rmse'])
        strm_rmse = strm_rmse[:-1]
        strm_avg_diff = np.array(default_strm.loc['strm_avg_diff'])
        strm_avg_diff = strm_avg_diff[:-1]
        strm_nse = np.array(default_strm.loc['strm_nse'])
        strm_nse = strm_nse[:-1]
        strm_corr = np.array(default_strm.loc['strm_corr'])
        strm_corr = strm_corr[:-1]
        strm_ubrmse = np.array(default_strm.loc['strm_ubrmse'])
        strm_ubrmse = strm_ubrmse[:-1]
        strm_bias = np.array(default_strm.loc['strm_avg_diff'])
        strm_bias = strm_bias[:-1]
        strm_obs = np.array(default_strm.loc['strm_obs'])
        strm_obs = strm_obs[:-1]
        strm_mae = np.array(default_strm.loc['strm_mae'])
        strm_mae = strm_mae[:-1]
        not_nan_idx_strm = np.where(
            np.isnan(strm_rmse) == False
        )
        # get the watersheds we are going to work with
        watersheds = np.array(default_strm.columns)
        watersheds = watersheds[:-1]
        # let's get the predictors in terms of watersheds
        gldas_precip_strm = np.zeros(len(watersheds))
        #gldas_temp_strm = np.zeros(len(watersheds))
        gldas_pet_strm = np.zeros(len(watersheds))
        #gldas_et_strm = np.zeros(len(watersheds))
        canopy_height_strm = np.zeros(len(watersheds))
        for w,wat in enumerate(watersheds):
            this_perc = intersection_info[wat][1]
            this_tiles = intersection_info[wat][0]
            precip_vals = np.zeros(len(this_tiles))
            #temp_vals = np.zeros(len(this_tiles))
            pet_vals = np.zeros(len(this_tiles))
            #et_vals = np.zeros(len(this_tiles))
            canopy_vals = np.zeros(len(this_tiles))
            percentages = np.zeros(len(this_tiles))
            for t,ti in enumerate(this_tiles):
                this_tile_idx = np.where(
                    pixels == ti
                )
                precip_vals[t] = gldas_precip_avg[this_tile_idx]
                #temp_vals[t] = gldas_temp_avg[this_tile_idx]
                pet_vals[t] = gldas_pet_avg[this_tile_idx]
                #et_vals[t] = gldas_et_avg[this_tile_idx]
                canopy_vals[t] = canopy_height_avg[this_tile_idx]
            gldas_precip_strm[w] = np.average(
                precip_vals,weights=this_perc
            )
            #gldas_temp_strm[w] = np.average(
            #    temp_vals,weights=this_perc
            #)
            gldas_pet_strm[w] = np.average(
                pet_vals,weights=this_perc
            )
            #gldas_et_strm[w] = np.average(
            #    et_vals,weights=this_perc
            #)
            canopy_height_strm[w] = np.average(
                canopy_vals,weights=this_perc
            )
        
        # get the ones where we need to do math
        p_pet = gldas_precip_avg/gldas_pet_avg
        #pet_et = gldas_pet_avg - gldas_et_avg
        #inv_canopy_height = 1/canopy_height_avg
        p_pet_strm = gldas_precip_strm/gldas_pet_strm
        #pet_et_strm = gldas_pet_strm - gldas_et_strm
        #inv_canopy_height_strm = 1/canopy_height_strm
        # get rid of nans everywhere
        model_rmse = model_rmse[not_nan_idx]
        model_r2 = model_r2[not_nan_idx]
        model_corr = model_corr[not_nan_idx]
        model_ubrmse = model_ubrmse[not_nan_idx]
        model_bias = model_bias[not_nan_idx]
        gldas_precip_avg = gldas_precip_avg[not_nan_idx]
        #gldas_temp_avg = gldas_temp_avg[not_nan_idx]
        #inv_canopy_height = inv_canopy_height[not_nan_idx]
        canopy_height = canopy_height_avg[not_nan_idx]
        p_pet = p_pet[not_nan_idx]
        #pet_et = pet_et[not_nan_idx]
        # get rid of nans everywhere for streamflow
        strm_rmse = strm_rmse[not_nan_idx_strm]
        strm_avg_diff = strm_avg_diff[not_nan_idx_strm]
        strm_nse = strm_nse[not_nan_idx_strm]
        strm_corr = strm_corr[not_nan_idx_strm]
        strm_ubrmse = strm_ubrmse[not_nan_idx_strm]
        strm_bias = strm_bias[not_nan_idx_strm]
        strm_obs = strm_obs[not_nan_idx_strm]
        strm_mae = strm_mae[not_nan_idx_strm]
        gldas_precip_strm = gldas_precip_strm[not_nan_idx_strm]
        #gldas_temp_strm = gldas_temp_strm[not_nan_idx_strm]
        #inv_canopy_height_strm = inv_canopy_height_strm[not_nan_idx_strm]
        canopy_height_strm = canopy_height_strm[not_nan_idx_strm]
        p_pet_strm = p_pet_strm[not_nan_idx_strm]
        #pet_et_strm = pet_et_strm[not_nan_idx_strm]
        # now for where we want to look at the multiple linear regression
        can_mean = np.mean(canopy_height)
        can_std = np.std(canopy_height)
        can_norm = (canopy_height - can_mean)/can_std
        #can_max = np.max(canopy_height)
        #can_norm = canopy_height/can_max
        gldas_precip_mean = np.mean(gldas_precip_avg)
        gldas_precip_std = np.std(gldas_precip_avg)
        gldas_precip_norm = (
            (
                gldas_precip_avg - gldas_precip_mean
            )/gldas_precip_std
        )
        #gldas_precip_max = np.max(gldas_precip_avg)
        #gldas_precip_norm = gldas_precip_avg/gldas_precip_max
        p_pet_mean = np.mean(p_pet)
        p_pet_std = np.std(p_pet)
        p_pet_norm = (p_pet - p_pet_mean)/p_pet_std
        #p_pet_max = np.max(p_pet)
        #p_pet_norm = p_pet/p_pet_max
        can_mean_strm = np.mean(canopy_height_strm)
        can_std_strm = np.std(canopy_height_strm)
        can_norm_strm = (
            (
                canopy_height_strm - can_mean_strm
            )/can_std_strm
        )
        #can_max_strm = np.max(canopy_height_strm)
        #can_norm_strm = canopy_height_strm/can_max_strm
        gldas_precip_mean_strm = np.mean(gldas_precip_strm)
        gldas_precip_std_strm = np.std(gldas_precip_strm)
        gldas_precip_norm_strm = (
            (
                gldas_precip_strm - gldas_precip_mean_strm
            )/gldas_precip_std_strm
        )
        #gldas_precip_max_strm = np.max(gldas_precip_strm)
        #gldas_precip_norm_strm = gldas_precip_strm/gldas_precip_max_strm
        p_pet_mean_strm = np.mean(p_pet_strm)
        p_pet_std_strm = np.std(p_pet_strm)
        p_pet_norm_strm = (
            (
                p_pet_strm - p_pet_mean_strm
            )/p_pet_std_strm
        )
        #p_pet_max_strm = np.max(p_pet_strm)
        #p_pet_norm_strm = p_pet_strm/p_pet_max_strm
        data_precip_strm = pd.DataFrame(
            {
                'precip':gldas_precip_norm_strm,
                'canopy_height':can_norm_strm,
                'strm_rmse':strm_rmse
            }
        )
        model_precip_strm = ols(
            "strm_rmse ~ precip + canopy_height",data_precip_strm
        ).fit()
        model_precip_coefs_strm = model_precip_strm._results.params
        precip_a0_strm = model_precip_coefs_strm[0]
        precip_a1_strm = model_precip_coefs_strm[1]
        precip_a2_strm = model_precip_coefs_strm[2]
        precip_r2_strm = model_precip_strm.rsquared
        data_ppet_strm = pd.DataFrame(
            {
                'p_pet':p_pet_norm_strm,
                'canopy_height':can_norm_strm,
                'strm_rmse':strm_rmse
            }
        )
        model_ppet_strm = ols(
            "strm_rmse ~ p_pet + canopy_height",data_ppet_strm
        ).fit()
        model_ppet_coefs_strm = model_ppet_strm._results.params
        ppet_a0_strm = model_ppet_coefs_strm[0]
        ppet_a1_strm = model_ppet_coefs_strm[1]
        ppet_a2_strm = model_ppet_coefs_strm[2]
        ppet_r2_strm = model_ppet_strm.rsquared
        data_precip = pd.DataFrame(
            {
                'precip':gldas_precip_norm,
                'canopy_height':can_norm,
                'le_ubrmse':model_ubrmse
            }
        )
        model_precip = ols("le_ubrmse ~ precip + canopy_height",data_precip).fit()
        model_precip_coefs = model_precip._results.params
        precip_a0 = model_precip_coefs[0]
        precip_a1 = model_precip_coefs[1]
        precip_a2 = model_precip_coefs[2]
        precip_r2 = model_precip.rsquared
        data_ppet = pd.DataFrame(
            {
                'p_pet':p_pet_norm,
                'canopy_height':can_norm,
                'le_ubrmse':model_ubrmse
            }
        )
        model_ppet = ols("le_ubrmse ~ p_pet + canopy_height",data_ppet).fit()
        model_ppet_coefs = model_ppet._results.params
        ppet_a0 = model_ppet_coefs[0]
        ppet_a1 = model_ppet_coefs[1]
        ppet_a2 = model_ppet_coefs[2]
        ppet_r2 = model_ppet.rsquared
        data_precip_bias = pd.DataFrame(
            {
                'precip':gldas_precip_norm,
                'canopy_height':can_norm,
                'le_bias':model_bias
            }
        )
        model_precip_bias = ols(
            "le_bias ~ precip + canopy_height",data_precip_bias
        ).fit()
        model_precip_bias_coefs = model_precip_bias._results.params
        precip_bias_a0 = model_precip_bias_coefs[0]
        precip_bias_a1 = model_precip_bias_coefs[1]
        precip_bias_a2 = model_precip_bias_coefs[2]
        precip_bias_r2 = model_precip_bias.rsquared
        data_precip_bias_strm = pd.DataFrame(
            {
                'precip':gldas_precip_norm_strm,
                'canopy_height':can_norm_strm,
                'strm_bias':strm_bias
            }
        )
        model_precip_bias_strm = ols(
            "strm_bias ~ precip + canopy_height",data_precip_bias_strm
        ).fit()
        model_precip_bias_coefs_strm = model_precip_bias_strm._results.params
        precip_bias_a0_strm = model_precip_bias_coefs_strm[0]
        precip_bias_a1_strm = model_precip_bias_coefs_strm[1]
        precip_bias_a2_strm = model_precip_bias_coefs_strm[2]
        precip_bias_r2_strm = model_precip_bias_strm.rsquared
        # and now for all of the cases where we try to see how much value
        # canopy height adds
        m_ppet,b_ppet,na,nan,nana = (
            stats.linregress(p_pet,model_bias)
        )
        remain_res_ppet = (
            model_bias -
            m_ppet*p_pet -
            b_ppet
        )
        m_precip,b_precip,na,nan,nana = (
            stats.linregress(gldas_precip_avg,model_bias)
        )
        remain_res_precip = (
            model_bias -
            m_precip*gldas_precip_avg -
            b_precip
        )
        m_precip_strm,b_precip_strm,na,nan,nana = (
            stats.linregress(gldas_precip_strm,strm_bias)
        )
        remain_res_precip_strm = (
            strm_bias -
            m_precip_strm*gldas_precip_strm -
            b_precip_strm
        )
        m_ppet_strm,b_ppet_strm,na,nan,nana = (
            stats.linregress(p_pet_strm,strm_bias)
        )
        remain_res_ppet_strm = (
            strm_bias -
            m_ppet_strm*p_pet_strm -
            b_ppet_strm
        )
        # let's make sure default streamflow error is normally distributed
        strm_rmse_normalized = strm_rmse# - np.mean(
        #    strm_rmse
        #)
        #strm_rmse_normalized = strm_rmse_normalized/np.std(
        #    strm_rmse_normalized
        #)
        strm_rmse_norm = strm_rmse/strm_obs
        strm_rmse_norm_normalized = strm_rmse_norm# - np.mean(
        #    strm_rmse_norm
        #)
        #strm_rmse_norm_normalized = strm_rmse_norm_normalized/np.std(
        #    strm_rmse_norm_normalized
        #)
        binwidth = 0.15
        plt.figure()
        plt.hist(
            strm_rmse_normalized,
            bins=np.arange(
                min(strm_rmse_normalized),
                max(strm_rmse_normalized) + binwidth,
                binwidth
            )
        )
        save_name = os.path.join(
            plots_dir,
            'strm_rmse_hist.png'
        )
        plt.savefig(save_name)
        plt.close()
        plt.figure()
        plt.hist(
            strm_rmse_norm_normalized,
            bins=np.arange(
                min(strm_rmse_norm_normalized),
                max(strm_rmse_norm_normalized) + binwidth,
                binwidth
            )
        )
        save_name = os.path.join(
            plots_dir,
            'strm_rmse_norm_hist.png'
        )
        plt.savefig(save_name)
        plt.close()
        too_small_idx = np.where(strm_obs < 0.05)
        big_enough_idx = np.where(strm_obs >= 0.05)
        print('number of catchments removed:')
        print(np.shape(too_small_idx))
        strm_rmse_big = strm_rmse[big_enough_idx]
        strm_obs_big = strm_obs[big_enough_idx]
        strm_rmse_big_normalized = strm_rmse_big# - np.mean(
        #    strm_rmse_big
        #)
        #strm_rmse_big_normalized = strm_rmse_big_normalized/np.std(
        #    strm_rmse_big_normalized
        #)
        strm_rmse_big_norm = strm_rmse_big/strm_obs_big
        strm_rmse_big_norm_normalized = strm_rmse_big_norm# - np.mean(
        #    strm_rmse_norm
        #)
        #strm_rmse_big_norm_normalized = strm_rmse_big_norm_normalized/np.std(
        #    strm_rmse_big_norm_normalized
        ##)
        plt.figure()
        plt.hist(
            strm_rmse_big_normalized,
            bins=np.arange(
                min(strm_rmse_big_normalized),
                max(strm_rmse_big_normalized) + binwidth,
                binwidth
            )
        )
        save_name = os.path.join(
            plots_dir,
            'strm_rmse_big_hist.png'
        )
        plt.savefig(save_name)
        plt.close()
        plt.figure()
        plt.hist(
            strm_rmse_big_norm_normalized,
            bins=np.arange(
                min(strm_rmse_big_norm_normalized),
                max(strm_rmse_big_norm_normalized) + binwidth,
                binwidth
            )
        )
        save_name = os.path.join(
            plots_dir,
            'strm_rmse_big_norm_hist.png'
        )
        plt.savefig(save_name)
        plt.close()
        strm_mae_normalized = strm_mae# - np.mean(
        #    strm_mae
        #)
        #strm_mae_normalized = strm_mae_normalized/np.std(
        #    strm_mae_normalized
        #)
        strm_mae_norm = strm_mae/strm_obs
        strm_mae_norm_normalized = strm_mae_norm# - np.mean(
        #    strm_mae_norm
        #)
        #strm_mae_norm_normalized = strm_mae_norm_normalized/np.std(
        #    strm_mae_norm_normalized
        #)
        plt.figure()
        plt.hist(
            strm_mae_normalized,
            bins=np.arange(
                min(strm_mae_normalized),
                max(strm_mae_normalized) + binwidth,
                binwidth
            )
        )
        save_name = os.path.join(
            plots_dir,
            'strm_mae_hist.png'
        )
        plt.savefig(save_name)
        plt.close()
        plt.figure()
        plt.hist(
            strm_mae_norm_normalized,
            bins=np.arange(
                min(strm_mae_norm_normalized),
                max(strm_mae_norm_normalized) + binwidth,
                binwidth
            )
        )
        save_name = os.path.join(
            plots_dir,
            'strm_mae_norm_hist.png'
        )
        plt.savefig(save_name)
        plt.close()
        too_small_idx = np.where(strm_obs < 0.05)
        big_enough_idx = np.where(strm_obs >= 0.05)
        print('number of catchments removed:')
        print(np.shape(too_small_idx))
        strm_mae_big = strm_mae[big_enough_idx]
        strm_obs_big = strm_obs[big_enough_idx]
        strm_mae_big_normalized = strm_mae_big #- np.mean(
        #    strm_mae_big
        #)
        #strm_mae_big_normalized = strm_mae_big_normalized/np.std(
        #    strm_mae_big_normalized
        #)
        strm_mae_big_norm = strm_mae_big/strm_obs_big
        strm_mae_big_norm_normalized = strm_mae_big_norm# - np.mean(
        #    strm_mae_norm
        #)
        #strm_mae_big_norm_normalized = strm_mae_big_norm_normalized/np.std(
        #    strm_mae_big_norm_normalized
        #)
        plt.figure()
        plt.hist(
            strm_mae_big_normalized,
            bins=np.arange(
                min(strm_mae_big_normalized),
                max(strm_mae_big_normalized) + binwidth,
                binwidth
            )
        )
        save_name = os.path.join(
            plots_dir,
            'strm_mae_big_hist.png'
        )
        plt.savefig(save_name)
        plt.close()
        plt.figure()
        plt.hist(
            strm_mae_big_norm_normalized,
            bins=np.arange(
                min(strm_mae_big_norm_normalized),
                max(strm_mae_big_norm_normalized) + binwidth,
                binwidth
            )
        )
        save_name = os.path.join(
            plots_dir,
            'strm_mae_big_norm_hist.png'
        )
        plt.savefig(save_name)
        plt.close()
        # let's just look at within-PFT variability
        p_d = pft_dist()
        pft_distribution = p_d.get_pso_pfts(pft_info,intersection_info)
        # at the pixel scale, we have enough >80% pixels to make this plot for
        # Needleleaf, Broadleaf, Shrub, Cool c3 grass, crop
        # at the watershed scale we have Needleleaf, Broadleaf, warm c4 grass
        need_pix = pft_distribution['pix_80']['pixels'].loc['Needleleaf']
        num_del = 0
        for p,pix in enumerate(need_pix):
            if pix in nan_pix:
                need_pix = np.delete(need_pix,p - num_del)
                num_del += 1
        need_le_ubrmse = np.zeros(len(need_pix))
        need_le_obs = np.zeros(len(need_pix))
        need_precip = np.zeros(len(need_pix))
        need_canopy = np.zeros(len(need_pix))
        for p,pix in enumerate(need_pix):
            this_idx = np.where(not_nan_pix == pix)
            need_le_ubrmse[p] = model_ubrmse[this_idx]
            need_le_obs[p] = le_obs[this_idx]
            need_precip[p] = gldas_precip_avg[this_idx]
            need_canopy[p] = canopy_height_avg[this_idx]
        need_le_ubrmse_norm = need_le_ubrmse/need_le_obs
        broad_pix = pft_distribution['pix_80']['pixels'].loc['Broadleaf']
        num_del = 0
        for p,pix in enumerate(broad_pix):
            if pix in nan_pix:
                broad_pix = np.delete(broad_pix,p - num_del)
                num_del += 1
        broad_le_ubrmse = np.zeros(len(broad_pix))
        broad_le_obs = np.zeros(len(broad_pix))
        broad_precip = np.zeros(len(broad_pix))
        broad_canopy = np.zeros(len(broad_pix))
        for p,pix in enumerate(broad_pix):
            this_idx = np.where(not_nan_pix == pix)
            broad_le_ubrmse[p] = model_ubrmse[this_idx]
            broad_le_obs[p] = le_obs[this_idx]
            broad_precip[p] = gldas_precip_avg[this_idx]
            broad_canopy[p] = canopy_height_avg[this_idx]
        broad_le_ubrmse_norm = broad_le_ubrmse/broad_le_obs
        shrub_pix = pft_distribution['pix_80']['pixels'].loc['Shrub']
        num_del = 0
        for p,pix in enumerate(shrub_pix):
            if pix in nan_pix:
                shrub_pix = np.delete(shrub_pix,p)
                num_del += 1
        shrub_le_ubrmse = np.zeros(len(shrub_pix))
        shrub_le_obs = np.zeros(len(shrub_pix))
        shrub_precip = np.zeros(len(shrub_pix))
        shrub_canopy = np.zeros(len(shrub_pix))
        for p,pix in enumerate(shrub_pix):
            this_idx = np.where(not_nan_pix == pix)
            shrub_le_ubrmse[p] = model_ubrmse[this_idx]
            shrub_le_obs[p] = le_obs[this_idx]
            shrub_precip[p] = gldas_precip_avg[this_idx]
            shrub_canopy[p] = canopy_height_avg[this_idx]
        shrub_le_ubrmse_norm = shrub_le_ubrmse/shrub_le_obs
        c3grass_pix = pft_distribution['pix_80']['pixels'].loc['Cool c3 grass']
        num_del = 0
        for p,pix in enumerate(c3grass_pix):
            if pix in nan_pix:
                c3grass = np.delete(c3grass,p)
                num_del += 1
        c3grass_le_ubrmse = np.zeros(len(c3grass_pix))
        c3grass_le_obs = np.zeros(len(c3grass_pix))
        c3grass_precip = np.zeros(len(c3grass_pix))
        c3grass_canopy = np.zeros(len(c3grass_pix))
        for p,pix in enumerate(c3grass_pix):
            this_idx = np.where(not_nan_pix == pix)
            c3grass_le_ubrmse[p] = model_ubrmse[this_idx]
            c3grass_le_obs[p] = le_obs[this_idx]
            c3grass_precip[p] = gldas_precip_avg[this_idx]
            c3grass_canopy[p] = canopy_height_avg[this_idx]
        c3grass_le_ubrmse_norm = c3grass_le_ubrmse/c3grass_le_obs
        crop_pix = pft_distribution['pix_80']['pixels'].loc['Crop']
        num_del = 0
        for p,pix in enumerate(crop_pix):
            if pix in nan_pix:
                crop = np.delete(crop,p)
                num_del += 1
        crop_le_ubrmse = np.zeros(len(crop_pix))
        crop_le_obs = np.zeros(len(crop_pix))
        crop_precip = np.zeros(len(crop_pix))
        crop_canopy = np.zeros(len(crop_pix))
        for p,pix in enumerate(crop_pix):
            this_idx = np.where(not_nan_pix == pix)
            crop_le_ubrmse[p] = model_ubrmse[this_idx]
            crop_le_obs[p] = le_obs[this_idx]
            crop_precip[p] = gldas_precip_avg[this_idx]
            crop_canopy[p] = canopy_height_avg[this_idx]
        crop_le_ubrmse_norm = crop_le_ubrmse/crop_le_obs
        # let's do this for stream
        need_wat = pft_distribution['wat_80']['watersheds'].loc[
            'Needleleaf'
        ]
        need_strm_rmse = np.zeros(len(need_wat))
        need_strm_obs = np.zeros(len(need_wat))
        need_strm_precip = np.zeros(len(need_wat))
        need_strm_canopy = np.zeros(len(need_wat))
        for w,wat in enumerate(need_wat):
            this_idx = np.where(watersheds == wat)
            need_strm_rmse[w] = strm_rmse[this_idx]
            need_strm_obs[w] = strm_obs[this_idx]
            need_strm_precip[w] = gldas_precip_strm[this_idx]
            need_strm_canopy[w] = canopy_height_strm[this_idx]
        need_strm_rmse_norm = need_strm_rmse/need_strm_obs
        broad_wat = pft_distribution['wat_80']['watersheds'].loc[
            'Broadleaf'
        ]
        broad_strm_rmse = np.zeros(len(broad_wat))
        broad_strm_obs = np.zeros(len(broad_wat))
        broad_strm_precip = np.zeros(len(broad_wat))
        broad_strm_canopy = np.zeros(len(broad_wat))
        for w,wat in enumerate(broad_wat):
            this_idx = np.where(watersheds == wat)
            broad_strm_rmse[w] = strm_rmse[this_idx]
            broad_strm_obs[w] = strm_obs[this_idx]
            broad_strm_precip[w] = gldas_precip_strm[this_idx]
            broad_strm_canopy[w] = canopy_height_strm[this_idx]
        broad_strm_rmse_norm = broad_strm_rmse/broad_strm_obs

        # let's make the plots where we just look within PFTs
        x_names = [
            'need_mean_annual_precip',
            'need_canopy_height',
            'broad_mean_annual_precip',
            'broad_canopy_height',
            'shrub_mean_annual_precip',
            'shrub_canopy_height',
            'c3grass_mean_annual_precip',
            'c3grass_canopy_height',
            'crop_mean_annual_precip',
            'crop_canopy_height',
            'need_strm_mean_annual_precip',
            'need_strm_canopy_height',
            'broad_strm_mean_annual_precip',
            'broad_strm_canopy_height'
        ]
        x_vals = [
            need_precip,
            need_canopy,
            broad_precip,
            broad_canopy,
            shrub_precip,
            shrub_canopy,
            c3grass_precip,
            c3grass_canopy,
            crop_precip,
            crop_canopy,
            need_strm_precip,
            need_strm_canopy,
            broad_strm_precip,
            broad_strm_canopy
        ]
        x_units = [
            'mm/day',
            'm',
            'mm/day',
            'm',
            'mm/day',
            'm',
            'mm/day',
            'm',
            'mm/day',
            'm',
            'mm/day',
            'm',
            'mm/day',
            'm'
        ]
        y_names = [
            'need_le_ubrmse_norm',
            'need_le_ubrmse_norm',
            'broad_le_ubrmse_norm',
            'broad_le_ubrmse_norm',
            'shrub_le_ubrmse_norm',
            'shrub_le_ubrmse_norm',
            'c3grass_le_ubrmse_norm',
            'c3grass_le_ubrmse_norm',
            'crop_le_ubrmse_norm',
            'crop_le_ubrmse_norm',
            'need_strm_rmse_norm',
            'need_strm_rmse_norm',
            'broad_strm_rmse_norm',
            'broad_strm_rmse_norm'
        ]
        y_vals = [
            need_le_ubrmse_norm,
            need_le_ubrmse_norm,
            broad_le_ubrmse_norm,
            broad_le_ubrmse_norm,
            shrub_le_ubrmse_norm,
            shrub_le_ubrmse_norm,
            c3grass_le_ubrmse_norm,
            c3grass_le_ubrmse_norm,
            crop_le_ubrmse_norm,
            crop_le_ubrmse_norm,
            need_strm_rmse_norm,
            need_strm_rmse_norm,
            broad_strm_rmse_norm,
            broad_strm_rmse_norm
        ]
        y_units = [
            '-',
            '-',
            '-',
            '-',
            '-',
            '-',
            '-',
            '-',
            '-',
            '-',
            '-',
            '-',
            '-',
            '-',
            '-'
        ]
        for x,x_nam in enumerate(x_names):
            slope,intercept,r_value,na,nan = (
                stats.linregress(x_vals[x],y_vals[x])
            )
            save_name = os.path.join(
                plots_dir,
                'scatter_{x}_vs_{y}.png'.format(
                    x=x_nam,y=y_names[x]
                )
            )
            plt.figure()
            plt.scatter(x_vals[x],y_vals[x],s=1)
            plt.plot(
                x_vals[x],
                intercept+slope*x_vals[x],
                label='best fit line'
            )
            plt.annotate(
                'r = {:.2f}'.format(r_value),
                xy=(0.05,0.95),xycoords='axes fraction'
            )
            plt.xlabel(
                '{thing} ({unit})'.format(
                    thing=x_nam,unit=x_units[x]
                )
            )
            plt.ylabel(
                '{thing} ({unit})'.format(
                    thing=y_names[x],unit=y_units[x]
                )
            )
            plt.legend(loc=4)
            #this_y_bounds = y_axes[y_nam]
            #plt.ylim(this_y_bounds[0],this_y_bounds[1])
            plt.savefig(save_name)
            plt.close()
        # make all of the plots for le at teh pixel scale
        # decide what plots to make. code will make all possible combinations
        # of provided x and y combinations
        if plot_pso_scatter:
            x_names = [
                'mean_annual_precip',
                'precip_over_potential_et',
                'canopy_height'
            ]
            x_vals = [
                gldas_precip_avg,
                p_pet,
                canopy_height
            ]
            x_units = [
                'mm/day',
                '-',
                'm'
            ]
            y_names = [
                'catchcn_model_le_bias_gleam',
                'rem_resid_precip',
                'rem_resid_ppet'
            ]
            y_vals = [
                model_bias,remain_res_precip,remain_res_ppet
            ]
            y_units = [
                'W/m2','W/m2','W/m2'
            ]
            y_axes = {
                y_names[0]:[0,60],
                y_names[1]:[-40,40],
                y_names[2]:[-40,40]
            }
            for x,x_nam in enumerate(x_names):
                for y,y_nam in enumerate(y_names):
                    slope,intercept,r_value,na,nan = (
                        stats.linregress(x_vals[x],y_vals[y])
                    )
                    save_name = os.path.join(
                        plots_dir,
                        'scatter_{x}_vs_{y}.png'.format(
                            x=x_nam,y=y_nam
                        )
                    )
                    plt.figure()
                    plt.scatter(x_vals[x],y_vals[y],s=1)
                    plt.plot(
                        x_vals[x],
                        intercept+slope*x_vals[x],
                        label='best fit line'
                    )
                    plt.annotate(
                        'r = {:.2f}'.format(r_value),
                        xy=(0.05,0.95),xycoords='axes fraction'
                    )
                    plt.xlabel(
                        '{thing} ({unit})'.format(
                            thing=x_nam,unit=x_units[x]
                        )
                    )
                    plt.ylabel(
                        '{thing} ({unit})'.format(
                            thing=y_nam,unit=y_units[y]
                        )
                    )
                    plt.legend(loc=4)
                    this_y_bounds = y_axes[y_nam]
                    plt.ylim(this_y_bounds[0],this_y_bounds[1])
                    plt.savefig(save_name)
                    plt.close()
            # make all plots for streamflow at the watershed scale
            x_names = [
                'mean_annual_precip',
                'precip_over_potential_et',
                'canopy_height'
            ]
            x_vals = [
                gldas_precip_strm,
                p_pet_strm,
                canopy_height_strm
            ]
            x_units = [
                'mm/day',
                '-',
                'm'
            ]
            y_names = [
                'catchcn_strm_bias_camels',
                'rem_resid_precip_strm',
                'rem_resid_ppet_strm'
            ]
            y_vals = [
                strm_bias,
                remain_res_precip_strm,
                remain_res_ppet_strm
            ]
            y_units = [
                'mm/day','mm/day','mm/day'
            ]
            y_axes = {
                y_names[0]:[0,3],
                y_names[1]:[-2,2],
                y_names[2]:[-2,2]
            }
            for x,x_nam in enumerate(x_names):
                for y,y_nam in enumerate(y_names):
                    slope,intercept,r_value,na,nan = (
                        stats.linregress(x_vals[x],y_vals[y])
                    )
                    save_name = os.path.join(
                        plots_dir,
                        'scatter_{x}_vs_{y}.png'.format(
                            x=x_nam,y=y_nam
                        )
                    )
                    plt.figure()
                    plt.scatter(x_vals[x],y_vals[y],s=1)
                    plt.plot(
                        x_vals[x],
                        intercept+slope*x_vals[x],
                        label='best fit line'
                    )
                    plt.annotate(
                        'r = {:.2f}'.format(r_value),
                        xy=(0.05,0.95),xycoords='axes fraction'
                    )
                    plt.xlabel(
                        '{thing} ({unit})'.format(
                            thing=x_nam,unit=x_units[x]
                        )
                    )
                    plt.ylabel(
                        '{thing} ({unit})'.format(
                            thing=y_nam,unit=y_units[y]
                        )
                    )
                    plt.legend(loc=4)
                    this_y_bounds = y_axes[y_nam]
                    plt.ylim(this_y_bounds[0],this_y_bounds[1])
                    plt.savefig(save_name)
                    plt.close()
            intercepts = [
                precip_a0,
                ppet_a0,
                precip_a0_strm,
                ppet_a0_strm,
                precip_bias_a0,
                precip_bias_a0_strm
            ]
            a1s = [
                precip_a1,
                ppet_a1,
                precip_a1_strm,
                ppet_a1_strm,
                precip_bias_a1,
                precip_bias_a1_strm
            ]
            a2s = [
                precip_a2,
                ppet_a2,
                precip_a2_strm,
                ppet_a2_strm,
                precip_bias_a2,
                precip_bias_a2_strm
            ]
            a1_vals = [
                gldas_precip_norm,
                p_pet_norm,
                gldas_precip_norm_strm,
                p_pet_norm_strm,
                gldas_precip_norm,
                gldas_precip_norm_strm
            ]
            a2_vals = [
                can_norm,
                can_norm,
                can_norm_strm,
                can_norm_strm,
                can_norm,
                can_norm_strm
            ]
            a1_names = [
                'precip',
                'p_pet',
                'precip_strm',
                'p_pet_strm',
                'precip_bias',
                'precip_bias_strm'
            ]
            a2_names = [
                'canopy_height',
                'canopy_height',
                'canopy_height_strm',
                'canopy_height_strm',
                'canopy_height_bias',
                'canopy_height_bias_strm'
            ]
            r2_vals = [
                precip_r2,
                ppet_r2,
                precip_r2_strm,
                ppet_r2_strm,
                precip_bias_r2,
                precip_bias_r2_strm
            ]
            ys = [
                model_ubrmse,
                model_ubrmse,
                strm_rmse,
                strm_rmse,
                model_bias,
                strm_bias
            ]
            x_labs = [
                'predicted LE ubRMSE (W/m2)',
                'predicted LE ubRMSE (W/m2)',
                'predicted streamflow RMSE (mm/day)',
                'predicted streamflow RMSE (mm/day)',
                'predicted LE bias (W/m2)',
                'predicted streamflow bias (mm/day)'
            ]
            y_labs = [
                'actual LE ubRMSE (W/m2)',
                'actual LE ubRMSE (W/m2)',
                'actual streamflow RMSE (mm/day)',
                'actual streamflow RMSE (mm/day)',
                'actual LE bias (W/m2)',
                'actual streamflow bias (mm/day)'
            ]
            types = [
                'et',
                'et',
                'strm',
                'strm',
                'et',
                'strm'
            ]
            for i,inter in enumerate(intercepts):
                save_name = os.path.join(
                    plots_dir,
                    'scatter_a1_{}_a2_{}.png'.format(
                        a1_names[i],
                        a2_names[i]
                    )
                )
                this_x = (
                    inter +
                    a1s[i]*a1_vals[i] +
                    a2s[i]*a2_vals[i]
                )
                this_y = ys[i]
                x_min = np.min(this_x)
                x_max = np.max(this_x)
                y_min = np.min(this_y)
                y_max = np.max(this_y)
                if types[i] == 'et':
                    y_min_min = y_min - 10
                    y_max_max = y_max + 10
                elif types[i] == 'strm':
                    y_min_min = y_min - 1
                    y_max_max = y_max + 1
                x_lab = x_labs[i]
                y_lab = y_labs[i]
                plt.figure()
                plt.scatter(this_x,this_y,s=1)
                plt.plot(
                    [y_min_min,y_max_max],
                    [y_min_min,y_max_max]
                )
                plt.annotate(
                    'R2 = {:.2f}'.format(r2_vals[i]),
                    xy=(0.05,0.95),xycoords='axes fraction'
                )
                plt.annotate(
                    '{} = {:.2f}'.format(a1_names[i],a1s[i]),
                    xy=(0.05,0.9),xycoords='axes fraction'
                )
                plt.annotate(
                    '{} = {:.2f}'.format(a2_names[i],a2s[i]),
                    xy=(0.05,0.85),xycoords='axes fraction'
                )
                plt.xlim(y_min_min,y_max_max)
                plt.ylim(y_min_min,y_max_max)
                plt.xlabel(x_labs[i])
                plt.ylabel(y_labs[i])
                plt.savefig(save_name)
        # let's run my fun experiment for estimating a0, a1, and a2
        # possibilities
        #print(default_df)
        #le_obs = np.array(default_df.loc['le_obs'])
        #le_obs = le_obs[not_nan_idx]
        #g1_copied = copy.deepcopy(le_obs)
        #max_le = np.max(g1_copied)
        #min_le = np.min(g1_copied)
        #g1_init = (g1_copied - min_le)/(max_le - min_le)
        #g1_init = g1_init*3.5
        #g1_init = g1_init + 0.5
        #too_small_idx = np.where(g1_init < 0.5)
        #g1_init[too_small_idx] = 0.5
        #num_loops = 1000
        #a0s = np.zeros(num_loops)
        #a1s = np.zeros(num_loops)
        #a2s = np.zeros(num_loops)
        #r2s = np.zeros(num_loops)
        #min_g1 = np.zeros(num_loops)
        #max_g1 = np.zeros(num_loops)
        #mean_g1 = np.zeros(num_loops)
        #std_g1 = np.zeros(num_loops)
        #mu = 0
        #sigma = 1
        #num_g1s = len(le_obs)
        #for l in range(num_loops):
        #    if l % 10 == 0:
        #        print('iteration {}'.format(l))
        #    noise = np.random.normal(mu,sigma,num_g1s)
        #    this_g1 = g1_init + noise
        #    #this_g1 = np.random.rand(len(g1_init))
        #    #this_g1 = this_g1*4.5
        #    #this_g1 = this_g1 + 0.5
        #    too_small_idx = np.where(this_g1 < 0.5)
        #    this_g1[too_small_idx] = 0.5
        #    min_g1[l] = np.min(this_g1)
        #    max_g1[l] = np.max(this_g1)
        #    mean_g1[l] = np.mean(this_g1)
        #    std_g1[l] = np.std(this_g1)
        #    this_data = pd.DataFrame(
        #        {
        #            'precip':gldas_precip_norm,
        #            'canopy_height':can_norm,
        #            'g1':this_g1
        #        }
        #    )
        #    model_g1 = ols("g1 ~ precip + canopy_height",this_data).fit()
        #    this_coefs = model_g1._results.params
        #    a0s[l] = this_coefs[0]
        #    a1s[l] = this_coefs[1]
        #    a2s[l] = this_coefs[2]
        #    r2s[l] = model_g1.rsquared
        #min_a0 = np.min(a0s)
        #max_a0 = np.max(a0s)
        #mean_a0 = np.mean(a0s)
        #std_a0 = np.std(a0s)
        #ovr_min_g1 = np.mean(min_g1)
        #ovr_max_g1 = np.mean(max_g1)
        #ovr_mean_g1 = np.mean(mean_g1)
        #ovr_std_g1 = np.mean(std_g1)
        #min_a1 = np.min(a1s)
        #max_a1 = np.max(a1s)
        #mean_a1 = np.mean(a1s)
        #std_a1 = np.std(a1s)
        #min_a2 = np.min(a2s)
        #max_a2 = np.max(a2s)
        #mean_a2 = np.mean(a2s)
        #std_a2 = np.std(a2s)
        #min_r2 = np.min(r2s)
        #max_r2 = np.max(r2s)
        #mean_r2 = np.mean(r2s)
        #std_r2 = np.std(r2s)
        #print('|  var  |  min  |  max  |  mean |  std  |')
        #print('|-------|-------|-------|-------|-------|')
        #print(
        #    '| {} | {} | {} | {} | {} |'.format(
        #        ' g1  ',
        #        str(ovr_min_g1)[:5].ljust(5),
        #        str(ovr_max_g1)[:5],
        #        str(ovr_mean_g1)[:5],
        #        str(ovr_std_g1)[:5]
        #    )
        #)
        #print(
        #    '| {} | {} | {} | {} | {} |'.format(
        #        ' a0  ',
        #        str(min_a0)[:5],
        #        str(max_a0)[:5],
        #        str(mean_a0)[:5],
        #        str(std_a0)[:5]
        #    )
        #)
        #print(
        #    '| {} | {} | {} | {} | {} |'.format(
        #        ' a1  ',
        #        str(min_a1)[:5],
        #        str(max_a1)[:5],
        #        str(mean_a1)[:5],
        #        str(std_a1)[:5]
        #    )
        #)
        #print(
        #    '| {} | {} | {} | {} | {} |'.format(
        #        ' a2  ',
        #        str(min_a2)[:5],
        #        str(max_a2)[:5],
        #        str(mean_a2)[:5],
        #        str(std_a2)[:5]
        #    )
        #)
        #print(
        #    '| {} | {} | {} | {} | {} |'.format(
        #        ' r2  ',
        #        str(min_r2)[:5],
        #        str(max_r2)[:5],
        #        str(mean_r2)[:5],
        #        str(std_r2)[:5]
        #    )
        #)
        #sys.exit()
        # get the global best parameters from the PSO
        iteration_keys = list(all_info.keys())
        iterations_considered = 5
        this_it_key = iteration_keys[iterations_considered-1]
        best_positions = all_info[this_it_key]['global_best_positions']
        # assign these positions to their relevant variables
        if self.optimization_type == 'pft':
            a0_g1_dict = {
                'Needleleaf evergreen temperate tree':best_positions[0],
                'Needleleaf evergreen boreal tree':best_positions[0],
                'Needleleaf deciduous boreal tree':best_positions[0],
                'Broadleaf evergreen tropical tree':best_positions[1],
                'Broadleaf evergreen temperate tree':best_positions[1],
                'Broadleaf deciduous tropical tree':best_positions[1],
                'Broadleaf deciduous temperate tree':best_positions[1],
                'Broadleaf deciduous boreal tree':best_positions[1],
                'Broadleaf evergreen temperate shrub':best_positions[2],
                'Broadleaf deciduous boreal shrub':best_positions[2],
                'Broadleaf deciduous temperate shrub':best_positions[2],
                'Broadleaf deciduous temperate shrub[moisture +deciduous]': (
                    best_positions[2]
                ),
                'Broadleaf deciduous temperate shrub[moisture stress only]': (
                    best_positions[2]
                ),
                'Arctic c3 grass':best_positions[3],
                'Cool c3 grass':best_positions[3],
                'Cool c3 grass [moisture + deciduous]':best_positions[3],
                'Cool c3 grass [moisture stress only]':best_positions[3],
                'Warm c4 grass [moisture + deciduous]':best_positions[4],
                'Warm c4 grass [moisture stress only]':best_positions[4],
                'Warm c4 grass':best_positions[4],
                'Crop':best_positions[5],
                'Crop [moisture + deciduous]':best_positions[5],
                'Crop [moisture stress only]':best_positions[5],
                '(Spring temperate cereal)':best_positions[5],
                '(Irrigated corn)':best_positions[5],
                '(Soybean)':best_positions[5],
                '(Corn)':best_positions[5]
            }
        default_g1_dict = {
            'Needleleaf evergreen temperate tree':2.3,
            'Needleleaf evergreen boreal tree':2.3,
            'Needleleaf deciduous boreal tree':2.3,
            'Broadleaf evergreen tropical tree':4.1,
            'Broadleaf evergreen temperate tree':4.1,
            'Broadleaf deciduous tropical tree':4.4,
            'Broadleaf deciduous temperate tree':4.4,
            'Broadleaf deciduous boreal tree':4.4,
            'Broadleaf evergreen temperate shrub':4.7,
            'Broadleaf deciduous boreal shrub':4.7,
            'Broadleaf deciduous temperate shrub':4.7,
            'Broadleaf deciduous temperate shrub[moisture + deciduous]':4.7,
            'Broadleaf deciduous temperate shrub[moisture stress only]':4.7,
            'Arctic c3 grass':2.2,
            'Cool c3 grass':5.3,
            'Cool c3 grass [moisture + deciduous]':5.3,
            'Cool c3 grass [moisture stress only]':5.3,
            'Warm c4 grass [moisture + deciduous]':1.6,
            'Warm c4 grass [moisture stress only]':1.6,
            'Warm c4 grass':1.6,
            'Crop':5.79,
            'Crop [moisture + deciduous]':5.79,
            'Crop [moisture stress only]':5.79,
            '(Spring temperate cereal)':5.79,
            '(Irrigated corn)':5.79,
            '(Soybean)':5.79,
            '(Corn)':5.79
        }
        if self.optimization_type == 'ef':
            bi_g1_dict = {
                'Needleleaf evergreen temperate tree':best_positions[0],
                'Needleleaf evergreen boreal tree':best_positions[0],
                'Needleleaf deciduous boreal tree':best_positions[0],
                'Broadleaf evergreen tropical tree':best_positions[1],
                'Broadleaf evergreen temperate tree':best_positions[1],
                'Broadleaf deciduous tropical tree':best_positions[1],
                'Broadleaf deciduous temperate tree':best_positions[1],
                'Broadleaf deciduous boreal tree':best_positions[1],
                'Broadleaf evergreen temperate shrub':best_positions[2],
                'Broadleaf deciduous boreal shrub':best_positions[2],
                'Broadleaf deciduous temperate shrub':best_positions[2],
                'Broadleaf deciduous temperate shrub[moisture +deciduous]': (
                    best_positions[2]
                ),
                'Broadleaf deciduous temperate shrub[moisture stress only]': (
                    best_positions[2]
                ),
                'Arctic c3 grass':best_positions[3],
                'Cool c3 grass':best_positions[3],
                'Cool c3 grass [moisture + deciduous]':best_positions[3],
                'Cool c3 grass [moisture stress only]':best_positions[3],
                'Warm c4 grass [moisture + deciduous]':best_positions[4],
                'Warm c4 grass [moisture stress only]':best_positions[4],
                'Warm c4 grass':best_positions[4],
                'Crop':best_positions[5],
                'Crop [moisture + deciduous]':best_positions[5],
                'Crop [moisture stress only]':best_positions[5],
                '(Spring temperate cereal)':best_positions[5],
                '(Irrigated corn)':best_positions[5],
                '(Soybean)':best_positions[5],
                '(Corn)':best_positions[5]
            }
            a0 = best_positions[6]
            a1 = best_positions[7]
            a2 = best_positions[8]
        if all_info_compare != 'none':
            compare_iterations_considered = 5
            compare_keys = list(all_info_compare.keys())
            this_compare_key = compare_keys[compare_iterations_considered-1]
            best_positions_compare = (
                all_info_compare[this_compare_key]['global_best_positions']
            )
            aj_g1_dict = {
                'Needleleaf evergreen temperate tree':best_positions_compare[0],
                'Needleleaf evergreen boreal tree':best_positions_compare[0],
                'Needleleaf deciduous boreal tree':best_positions_compare[0],
                'Broadleaf evergreen tropical tree':best_positions_compare[1],
                'Broadleaf evergreen temperate tree':best_positions_compare[1],
                'Broadleaf deciduous tropical tree':best_positions_compare[1],
                'Broadleaf deciduous temperate tree':best_positions_compare[1],
                'Broadleaf deciduous boreal tree':best_positions_compare[1],
                'Broadleaf evergreen temperate shrub':best_positions_compare[2],
                'Broadleaf deciduous boreal shrub':best_positions_compare[2],
                'Broadleaf deciduous temperate shrub':best_positions_compare[2],
                'Broadleaf deciduous temperate shrub[moisture +deciduous]': (
                    best_positions_compare[2]
                ),
                'Broadleaf deciduous temperate shrub[moisture stress only]': (
                    best_positions_compare[2]
                ),
                'Arctic c3 grass':best_positions_compare[3],
                'Cool c3 grass':best_positions_compare[3],
                'Cool c3 grass [moisture + deciduous]':best_positions_compare[3],
                'Cool c3 grass [moisture stress only]':best_positions_compare[3],
                'Warm c4 grass [moisture + deciduous]':best_positions_compare[4],
                'Warm c4 grass [moisture stress only]':best_positions_compare[4],
                'Warm c4 grass':best_positions_compare[4],
                'Crop':best_positions_compare[5],
                'Crop [moisture + deciduous]':best_positions_compare[5],
                'Crop [moisture stress only]':best_positions_compare[5],
                '(Spring temperate cereal)':best_positions_compare[5],
                '(Irrigated corn)':best_positions_compare[5],
                '(Soybean)':best_positions_compare[5],
                '(Corn)':best_positions_compare[5]
        }
        #ksat_alpha = best_positions[5]
        ksat_alpha = 4.5
        #ksat_beta = best_positions[6]
        ksat_beta = 5
        #ksat_const_1 = best_positions[7]
        ksat_const_1 = 3.5
        #ksat_const_2 = best_positions[8]
        ksat_const_2 = 1.5
        #ksat_sand_exp = best_positions[9]
        ksat_sand_exp = 0.13

        # calculate the parameter at each location
        # first for default ksat
        default_ksat = k0_vals
        # then for pso ksat
        #ks_max = k0_vals*(10**(ksat_const_1 - ksat_const_2*(sand_vals**(ksat_sand_exp))))
        #pso_ksat = ks_max - (
        #    (ks_max - k0_vals)/(1 + ((lai_vals/ksat_alpha)**ksat_beta))
        #)
        # then for default g1
        default_g1_init = np.zeros(len(norm_precip_vals))
        default_g1 = np.zeros(len(default_g1_init))
        for g,g1 in enumerate(default_g1_init):
            this_perc = [
                pft_info['pft_1_perc'].loc[tiles[g]],
                pft_info['pft_2_perc'].loc[tiles[g]],
                pft_info['pft_3_perc'].loc[tiles[g]],
                pft_info['pft_4_perc'].loc[tiles[g]]
            ]
            this_pfts = [
                pft_info['pft_1_name'].loc[tiles[g]],
                pft_info['pft_2_name'].loc[tiles[g]],
                pft_info['pft_3_name'].loc[tiles[g]],
                pft_info['pft_4_name'].loc[tiles[g]]
            ]
            effective_g1 = 0
            for p,perc in enumerate(this_perc):
                effective_g1 += (perc/100)*default_g1_dict[this_pfts[p]]
            default_g1[g] = effective_g1
        # then for pso g1
        #pso_g1_init = -0.163747 + 0.025*precip_vals
        if self.optimization_type == 'ef':
            pso_g1_init = np.zeros(len(norm_precip_vals))
            pso_g1 = np.zeros(len(pso_g1_init))
            a0_term = np.zeros(len(pso_g1_init))
            a1_term = np.zeros(len(pso_g1_init))
            a2_term = np.zeros(len(pso_g1_init))
            a1_a2_term = np.zeros(len(pso_g1_init))
            a0_a1_a2_term = np.zeros(len(pso_g1_init))
            bi_term = np.zeros(len(pso_g1_init))
            for g,g1 in enumerate(pso_g1_init):
                this_perc = [
                    pft_info['pft_1_perc'].loc[tiles[g]],
                    pft_info['pft_2_perc'].loc[tiles[g]],
                    pft_info['pft_3_perc'].loc[tiles[g]],
                    pft_info['pft_4_perc'].loc[tiles[g]]
                ]
                this_pfts = [
                    pft_info['pft_1_name'].loc[tiles[g]],
                    pft_info['pft_2_name'].loc[tiles[g]],
                    pft_info['pft_3_name'].loc[tiles[g]],
                    pft_info['pft_4_name'].loc[tiles[g]]
                ]
                effective_bi = 0
                for p,perc in enumerate(this_perc):
                    effective_bi += (perc/100)*bi_g1_dict[this_pfts[p]]
                pso_g1[g] = (
                    effective_bi*(
                        a0 + a1*norm_precip_vals[g] +
                        a2*norm_canopy_vals[g]
                    )
                )
                bi_term[g] = effective_bi
                a0_term[g] = a0
                a1_term[g] = a1*norm_precip_vals[g]
                a2_term[g] = a2*norm_canopy_vals[g]
                a1_a2_term[g] = a1_term[g] + a2_term[g]
                a0_a1_a2_term[g] = (
                    a0_term[g] + a1_term[g] + a2_term[g]
                )
            pso_g1 = np.where(pso_g1 < 0.5, 0.5, pso_g1)
        if self.optimization_type == 'pft':
            pso_g1_init = np.zeros(len(norm_precip_vals))
            pso_g1 = np.zeros(len(pso_g1_init))
            for g,g1 in enumerate(pso_g1_init):
                this_perc = [
                    pft_info['pft_1_perc'].loc[tiles[g]],
                    pft_info['pft_2_perc'].loc[tiles[g]],
                    pft_info['pft_3_perc'].loc[tiles[g]],
                    pft_info['pft_4_perc'].loc[tiles[g]]
                ]
                this_pfts = [
                    pft_info['pft_1_name'].loc[tiles[g]],
                    pft_info['pft_2_name'].loc[tiles[g]],
                    pft_info['pft_3_name'].loc[tiles[g]],
                    pft_info['pft_4_name'].loc[tiles[g]]
                ]
                effective_a0 = 0
                for p,perc in enumerate(this_perc):
                    effective_a0 += (perc/100)*a0_g1_dict[this_pfts[p]]
                pso_g1[g] = (
                    effective_a0
                )
            pso_g1 = np.where(pso_g1 < 0.5, 0.5, pso_g1)
        if all_info_compare != 'none':
            compare_g1_init = np.zeros(len(norm_precip_vals))
            compare_g1 = np.zeros(len(compare_g1_init))
            for g,g1 in enumerate(compare_g1_init):
                this_perc = [
                    pft_info['pft_1_perc'].loc[tiles[g]],
                    pft_info['pft_2_perc'].loc[tiles[g]],
                    pft_info['pft_3_perc'].loc[tiles[g]],
                    pft_info['pft_4_perc'].loc[tiles[g]]
                ]
                this_pfts = [
                    pft_info['pft_1_name'].loc[tiles[g]],
                    pft_info['pft_2_name'].loc[tiles[g]],
                    pft_info['pft_3_name'].loc[tiles[g]],
                    pft_info['pft_4_name'].loc[tiles[g]]
                ]
                effective_aj = 0
                for p,perc in enumerate(this_perc):
                    effective_aj += (perc/100)*aj_g1_dict[this_pfts[p]]
                compare_g1[g] = effective_aj
            compare_g1 = np.where(compare_g1 < 0.5, 0.5, compare_g1)
            compare_g1_diff = pso_g1 - compare_g1
        


        #sys.exit()
        # make a histogram of the default ksat versus final ksat
        #bins = np.arange(0,0.001+0.00005,0.00005)
        #plt.figure()
        #plt.hist(
        #    k0_vals,
        #    alpha=0.5,
        #    label='default k_sat',
        #    bins=bins
        #)
        #plt.hist(
        #    pso_ksat,
        #    alpha=0.5,
        #    label='pso k_sat',
        #    bins=bins
        #)
        #plt.legend()
        #savename = os.path.join(
        #    plots_dir,
        #    'pso_vs_default_ksat_histogram.png'
        #)
        #plt.savefig(savename)
        #plt.close()
        # lets get diff between default and ksat for all
        diff_pso_g1 = pso_g1 - default_g1
        perc_diff_pso_g1 = diff_pso_g1/default_g1
        #diff_pso_ksat = pso_ksat - default_ksat
        #perc_diff_pso_ksat = diff_pso_ksat/default_ksat
        # now get lats/lons and plot
        # define the lats and lons for the points
        lons = default_df.loc['lon']
        lons = lons.drop(labels=['all'])
        lats = default_df.loc['lat']
        lats = lats.drop(labels=['all'])
        # what are we plotting in each of these maps?
        if self.optimization_type == 'ef' and all_info_compare != 'none':
            names = [
                'default_g1','pso_g1',
                'diff_g1',
                'perc_diff_pso_g1',
                'a0_term_map',
                'a1_term_map','a2_term_map',
                'a1_a2_term_map','exp_diff_g1',
                'a0_a1_a2_term','bi_term'
            ]
            vals = [
                default_g1,pso_g1,diff_pso_g1,
                perc_diff_pso_g1,a0_term,a1_term,
                a2_term,a1_a2_term,compare_g1_diff,
                a0_a1_a2_term,bi_term
            ]
            plot_type = [
                'g1','g1','diff_g1',
                'perc_diff','a0','a1','a2',
                'a1a2','diff_g1_compare','a0a1a2',
                'bi'
            ]
        elif self.optimization_type == 'ef' and all_info_compare == 'none':
            names = [
                'default_g1','pso_g1',
                'diff_g1',
                'perc_diff_pso_g1',
                'a0_term_map',
                'a1_term_map','a2_term_map',
                'a1_a2_term_map','a0_a1_a2_term',
                'bi_term'
            ]
            vals = [
                default_g1,pso_g1,diff_pso_g1,
                perc_diff_pso_g1,a0_term,a1_term,
                a2_term,a1_a2_term,a0_a1_a2_term,
                bi_term
            ]
            plot_type = [
                'g1','g1','diff_g1',
                'perc_diff','a0','a1','a2',
                'a1a2','a0a1a2','bi'
            ]
        elif self.optimization_type == 'pft':
            names = [
                'default_g1','pso_g1',
                'diff_g1',
                'perc_diff_pso_g1'
            ]
            vals = [
                default_g1,pso_g1,diff_pso_g1,
                perc_diff_pso_g1
            ]
            plot_type = [
                'g1','g1','diff_g1','perc_diff'
            ]
        cmaps = {
            'g1':'rainbow',
            'diff_g1':'bwr',
            'perc_diff':'bwr',
            'a0':'bwr',
            'a1':'bwr',
            'a2':'bwr',
            'a1a2':'bwr',
            'a0a1a2':'bwr',
            'diff_g1_compare':'bwr',
            'bi':'bwr'
        }
        vmins = {
            'g1':0,
            'diff_g1':-6,
            'perc_diff':-6,
            'a0':-10,
            'a1':-5,
            'a2':-5,
            'a1a2':-5,
            'a0a1a2':-10,
            'diff_g1_compare':-5,
            'bi':-1
        }
        vmaxs = {
            'g1':3,
            'diff_g1':6,
            'perc_diff':5,
            'a0':10,
            'a1':5,
            'a2':5,
            'a1a2':5,
            'a0a1a2':10,
            'diff_g1_compare':5,
            'bi':1
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
            val = vals[p]
            # make the scatterplot
            scatter = plt.scatter(
                lons,lats,marker='s',s=1,c=val,transform=ccrs.PlateCarree(),
                cmap=cmaps[plot_type[p]],vmin=vmins[plot_type[p]],vmax=vmaxs[plot_type[p]]
            )
            # add colorbar
            plt.colorbar(scatter)
            # save
            savename = '{name}_map_{exp}.png'.format(
                name=names[p],exp=this_exp_name
            )
            savename = os.path.join(
                plots_dir,savename
            )
            plt.savefig(savename,dpi=300,bbox_inches='tight')
            plt.close()
        print('created maps of parameters')
