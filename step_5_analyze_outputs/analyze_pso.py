import pickle as pkl
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import stats

class analyze_pso:
    def __init__(self):
        pass
    def get_all_info(self,all_info_fname):
        with open(all_info_fname,'rb') as f:
            all_info = pkl.load(f)
        return all_info
    def get_parameter_names(self):
        parameter_names = [
            #'g1_aj_forests',
            #'g1_aj_croplands',
            #'g1_aj_grasslands',
            #'g1_aj_savannas',
            #'g1_aj_shrublands',
            #'ksat_alpha',
            #'ksat_beta',
            #'ksat_const_1',
            #'ksat_const_2',
            #'ksat_sand_exp'
            'g1_a1_forests',
            'g1_a1_shrublands',
            'g1_a1_savannas',
            'g1_a1_grasslands',
            'g1_a1_croplands'
            'g1_a0_forests',
            'g1_a0_shrublands',
            'g1_a0_savannas',
            'g1_a0_grasslands',
            'g1_a0_croplands'
        ]
        return parameter_names
    def plot_parameter_convergence(self,all_info,parameter_names,plots_dir,
                                   exp_names):
        this_exp_name = exp_names[2]
        # get the iterations that the pso ran
        iteration_str = list(all_info.keys())
        num_iterations = len(iteration_str)
        iterations = np.zeros(len(iteration_str))
        for t,this_str in enumerate(iteration_str):
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
        for i,it in enumerate(iteration_str):
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
            plt.legend()
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
        plt.legend()
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
                            all_info_compare='none'):
        this_exp_name = exp_names[2]
        # get env covariates maps
        precip_ds = nc.Dataset(precip_fname)
        precip_vals = np.array(precip_ds['precip'])
        lai_ds = nc.Dataset(lai_fname)
        lai_vals = np.array(lai_ds['lai'])
        sand_ds = nc.Dataset(sand_fname)
        sand_vals = np.array(sand_ds['sand_perc'])
        k0_ds = nc.Dataset(k0_fname)
        k0_vals = np.array(k0_ds['k_sat'])
        # let's compare using default bonnetti vals, just to see if we are
        # doing this correctly
        ks_max = (
            k0_vals*(10**(3.5-(1.5*(sand_vals**0.13))))
        )
        k_sat = (
            ks_max - ((ks_max - k0_vals)/(1 + ((lai_vals/4.5)**5)))
        )
        diff_bonetti_ksat = k_sat - k0_vals
        perc_diff_bonetti_ksat = (k_sat - k0_vals)/k0_vals
        # make a histogram to compare these distributions
        bins = np.arange(0,0.001+0.00005,0.00005)
        plt.figure()
        plt.hist(
            k0_vals,
            alpha=0.5,
            label='default k_sat',
            bins=bins
        )
        plt.hist(
            k_sat,
            alpha=0.5,
            label='bonnetti k_sat',
            bins=bins
        )
        plt.legend()
        savename = os.path.join(
            plots_dir,
            'bonnetti_vs_default_ksat_histogram.png'
        )
        plt.savefig(savename)
        plt.close()
        # let's make a scatter plot of precip values versus default model RMSE
        # to see if they co-vary
        print(default_df)
        model_rmse = np.array(default_df.loc['le_rmse'])
        model_rmse = model_rmse[:-1]
        rmse_slope,rmse_intercept,rmse_r_value,na,nan = (
            stats.linregress(precip_vals,model_rmse)
        )
        model_bias = np.array(default_df.loc['ave_le_diff'])
        model_bias = model_bias[:-1]
        bias_slope,bias_intercept,bias_r_value,na,nan = (
            stats.linregress(precip_vals,model_bias)
        )
        precip_save_name = os.path.join(
            plots_dir,
            'scatter_precip_vals_vs_default_catchmentcn45_le_rmse.png'
        )
        precip_bias_save_name = os.path.join(
            plots_dir,
            'scatter_precip_vals_vs_default_catchmentcn45_le_bias.png'
        )
        plt.figure()
        plt.scatter(precip_vals,model_rmse,s=1)
        plt.plot(
            precip_vals,
            rmse_intercept+rmse_slope*precip_vals,
            label='best fit line'
        )
        plt.text(200,35,'r = {:.2f}'.format(rmse_r_value))
        plt.xlabel('mean annual precipitation (mm/month)')
        plt.ylabel('Catchment-CN4.5 LE RMSE 2007 (W/m2)')
        plt.legend()
        plt.savefig(precip_save_name)
        plt.close()
        plt.figure()
        plt.scatter(precip_vals,model_bias,s=1)
        plt.plot(
            precip_vals,
            bias_intercept+bias_slope*precip_vals,
            label='best fit line'
        )
        plt.text(200,35,'r = {:.2f}'.format(bias_r_value))
        plt.xlabel('mean annual precipitation (mm/month)')
        plt.ylabel('Catchment-CN4.5 LE bias 2007 (W/m2)')
        plt.legend()
        plt.savefig(precip_bias_save_name)
        plt.close()
        # get the global best parameters from the PSO
        iteration_keys = list(all_info.keys())
        this_it_key = iteration_keys[-1]
        best_positions = all_info[this_it_key]['global_best_positions']
        # assign these positions to their relevant variables
        a1_g1_dict = {
            'Forest':best_positions[0],
            'Shrublands':best_positions[1],
            'Savannas':best_positions[2],
            'Grasslands':best_positions[3],
            'Croplands':best_positions[4]
        }
        a0_g1_dict = {
            'Forest':best_positions[5],
            'Shrublands':best_positions[6],
            'Savannas':best_positions[7],
            'Grasslands':best_positions[8],
            'Croplands':best_positions[9]
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
        if all_info_compare != 'none':
            best_positions_compare = (
                all_info_compare[this_it_key]['global_best_positions']
            )
            aj_g1_dict = {
                'Forest':best_positions_compare[0],
                'Shrublands':best_positions_compare[1],
                'Savannas':best_positions_compare[2],
                'Grasslands':best_positions_compare[3],
                'Croplands':best_positions_compare[4]
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
        ks_max = k0_vals*(10**(ksat_const_1 - ksat_const_2*(sand_vals**(ksat_sand_exp))))
        pso_ksat = ks_max - (
            (ks_max - k0_vals)/(1 + ((lai_vals/ksat_alpha)**ksat_beta))
        )
        # then for default g1
        default_g1_init = np.zeros(len(precip_vals))
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
        pso_g1_init = precip_vals
        pso_g1 = np.zeros(len(pso_g1_init))
        int_term = np.zeros(len(pso_g1_init))
        slope_term = np.zeros(len(pso_g1_init))
        for g,g1 in enumerate(pso_g1_init):
            this_perc = [
                pft_info['pft_1_perc'].loc[tiles[g]],
                pft_info['pft_2_perc'].loc[tiles[g]],
                pft_info['pft_3_perc'].loc[tiles[g]],
                pft_info['pft_4_perc'].loc[tiles[g]]
            ]
            this_pfts = [
                pft_info['pft_1_simple'].loc[tiles[g]],
                pft_info['pft_2_simple'].loc[tiles[g]],
                pft_info['pft_3_simple'].loc[tiles[g]],
                pft_info['pft_4_simple'].loc[tiles[g]]
            ]
            effective_a1 = 0
            effective_a0 = 0
            for p,perc in enumerate(this_perc):
                effective_a1 += (perc/100)*a1_g1_dict[this_pfts[p]]
                effective_a0 += (perc/100)*a0_g1_dict[this_pfts[p]]
            pso_g1[g] = effective_a0 + effective_a1*g1
            #pso_g1[g] = effective_a1
            int_term[g] = effective_a0
            slope_term[g] = effective_a1*g1
        if all_info_compare != 'none':
            compare_g1_init = precip_vals
            compare_g1 = np.zeros(len(compare_g1_init))
            for g,g1 in enumerate(compare_g1_init):
                this_perc = [
                    pft_info['pft_1_perc'].loc[tiles[g]],
                    pft_info['pft_2_perc'].loc[tiles[g]],
                    pft_info['pft_3_perc'].loc[tiles[g]],
                    pft_info['pft_4_perc'].loc[tiles[g]]
                ]
                this_pfts = [
                    pft_info['pft_1_simple'].loc[tiles[g]],
                    pft_info['pft_2_simple'].loc[tiles[g]],
                    pft_info['pft_3_simple'].loc[tiles[g]],
                    pft_info['pft_4_simple'].loc[tiles[g]]
                ]
                effective_aj = 0
                for p,perc in enumerate(this_perc):
                    effective_aj += (perc/100)*aj_g1_dict[this_pfts[p]]
                compare_g1[g] = effective_aj
            pso_g1 = np.where(pso_g1 < 0.5, 0.5, pso_g1)
            compare_g1 = np.where(compare_g1 < 0.5, 0.5, compare_g1)
            compare_g1_diff = pso_g1 - compare_g1
        # make a histogram of the default ksat versus final ksat
        bins = np.arange(0,0.001+0.00005,0.00005)
        plt.figure()
        plt.hist(
            k0_vals,
            alpha=0.5,
            label='default k_sat',
            bins=bins
        )
        plt.hist(
            pso_ksat,
            alpha=0.5,
            label='pso k_sat',
            bins=bins
        )
        plt.legend()
        savename = os.path.join(
            plots_dir,
            'pso_vs_default_ksat_histogram.png'
        )
        plt.savefig(savename)
        plt.close()
        # lets get diff between default and ksat for all
        diff_pso_g1 = pso_g1 - default_g1
        perc_diff_pso_g1 = diff_pso_g1/default_g1
        diff_pso_ksat = pso_ksat - default_ksat
        perc_diff_pso_ksat = diff_pso_ksat/default_ksat
        # now get lats/lons and plot
        # define the lats and lons for the points
        lons = default_df.loc['lon']
        lons = lons.drop(labels=['all'])
        lats = default_df.loc['lat']
        lats = lats.drop(labels=['all'])
        # what are we plotting in each of these maps?
        names = [
            'default_g1','pso_g1','default_ksat','pso_ksat',
            'diff_g1','diff_pso_ksat','diff_bonetti_ksat',
            'perc_diff_bonetti_ksat','perc_diff_pso_g1',
            'perc_diff_pso_ksat','intercept_term_map',
            'slope_term_map','pso_ef_g1_vs_pso_pft_g1'
        ]
        vals = [
            default_g1,pso_g1,default_ksat,pso_ksat,diff_pso_g1,diff_pso_ksat,
            diff_bonetti_ksat,perc_diff_bonetti_ksat,
            perc_diff_pso_g1,perc_diff_pso_ksat,int_term,slope_term,
            compare_g1_diff
        ]
        plot_type = [
            'g1','g1','ksat','ksat','diff_g1','diff_ksat','diff_ksat',
            'perc_diff','perc_diff','perc_diff','g1','g1','diff_g1'
        ]
        cmaps = {
            'g1':'winter',
            'ksat':'winter',
            'diff_g1':'bwr',
            'diff_ksat':'bwr',
            'perc_diff':'bwr'
        }
        vmins = {
            'g1':0,
            'ksat':0,
            'diff_g1':-1,
            'diff_ksat':-0.00005,
            'perc_diff':-5
        }
        vmaxs = {
            'g1':.5,
            'ksat':0.00005,
            'diff_g1':1,
            'diff_ksat':0.00005,
            'perc_diff':5
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
