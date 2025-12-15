import sys
import geopandas as gpd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

class compare_experiments:
    def __init__(self):
        pass
    def plot_obj_options(self,exp_1_names,exp_2_names,exp_1_pix_df,
                         exp_2_pix_df,exp_3_pix_df,exp_4_pix_df,
                         exp_1_wat_df,exp_2_wat_df,exp_3_wat_df,
                         exp_4_wat_df,plots_dir,intersection_info):
        camels_sheds = np.sort(list(intersection_info.keys()))
        # weights for the two different terms
        weights = [0.5,0.5]
        # get the different possible error metrics
        # we will do this for five different cases
        # first, the case where we calculate the different objective function
        # options at all watersheds
        # second, where we calculate where both ET and streamflow are strongly
        # under-predicted (Washington watersheds 12040500, 12189500)
        # third, where there ET is under-predicted and streamflow is
        # over-predicted (souther GA, 02314500)
        # finally, a watershed typical of the midwest, ET over-predicted and
        # streamflow pretty accurate (south datkota, 06360500)
        # for exp 1
        exp_1_rmse_et = weights[0]*(
            exp_1_pix_df['all'].loc['le_rmse']/
            exp_1_pix_df['all'].loc['le_obs']
        )
        exp_1_ubrmse_et = weights[0]*(
            exp_1_pix_df['all'].loc['le_ubrmse']/
            exp_1_pix_df['all'].loc['le_obs']
        )
        exp_1_rmse_strm = weights[1]*(
            exp_1_wat_df['all'].loc['strm_rmse']/
            exp_1_wat_df['all'].loc['strm_obs']
        )
        exp_1_ubrmse_strm = weights[1]*(
            exp_1_wat_df['all'].loc['strm_ubrmse']/
            exp_1_wat_df['all'].loc['strm_obs']
        )
        exp_1_bias_strm = weights[1]*np.abs(
            exp_1_wat_df['all'].loc['strm_avg_diff']/
            exp_1_wat_df['all'].loc['strm_obs']
        )
        exp_1_corr_strm = weights[1]*(
            exp_1_wat_df['all'].loc['strm_corr']
        )
        exp_1_le_ubrmse = np.array(exp_1_pix_df.loc['le_ubrmse'])
        exp_1_le_obs = np.array(exp_1_pix_df.loc['le_obs'])
        exp_1_ubrmse_norm_et = weights[0]*(
            np.nanmean(
                exp_1_le_ubrmse[:-1]/exp_1_le_ubrmse[:-1]
            )
        )
        exp_1_strm_rmse = np.array(exp_1_wat_df.loc['strm_rmse'])
        exp_1_strm_obs = np.array(exp_1_wat_df.loc['strm_obs'])
        exp_1_rmse_norm_strm = weights[1]*(
            np.nanmean(
                exp_1_strm_rmse[:-1]/exp_1_strm_obs[:-1]
            )
        )
        exp_1_strm_mae = np.array(exp_1_wat_df.loc['strm_mae'])
        exp_1_mae_norm_strm = weights[1]*(
            np.nanmean(
                exp_1_strm_mae[:-1]/exp_1_strm_mae[:-1]
            )
        )
        exp_1_mae_norm_def_strm = weights[1]*(
            np.nanmean(
                exp_1_strm_mae[:-1]/exp_1_strm_mae[:-1]
            )
        )
        ## for exp 2
        exp_2_rmse_et = weights[0]*(
            exp_2_pix_df['all'].loc['le_rmse']/
            exp_2_pix_df['all'].loc['le_obs']
        )
        exp_2_ubrmse_et = weights[0]*(
            exp_2_pix_df['all'].loc['le_ubrmse']/
            exp_2_pix_df['all'].loc['le_obs']
        )
        exp_2_rmse_strm = weights[1]*(
            exp_2_wat_df['all'].loc['strm_rmse']/
            exp_2_wat_df['all'].loc['strm_obs']
        )
        exp_2_ubrmse_strm = weights[1]*(
            exp_2_wat_df['all'].loc['strm_ubrmse']/
            exp_2_wat_df['all'].loc['strm_obs']
        )
        exp_2_bias_strm = weights[1]*np.abs(
            exp_2_wat_df['all'].loc['strm_avg_diff']/
            exp_2_wat_df['all'].loc['strm_obs']
        )
        exp_2_corr_strm = weights[1]*(
            exp_2_wat_df['all'].loc['strm_corr']
        )
        exp_2_le_ubrmse = np.array(exp_2_pix_df.loc['le_ubrmse'])
        exp_2_le_obs = np.array(exp_2_pix_df.loc['le_obs'])
        exp_2_ubrmse_norm_et = weights[0]*(
            np.nanmean(
                exp_2_le_ubrmse[:-1]/exp_1_le_ubrmse[:-1]
            )
        )
        exp_2_strm_rmse = np.array(exp_2_wat_df.loc['strm_rmse'])
        exp_2_strm_obs = np.array(exp_2_wat_df.loc['strm_obs'])
        exp_2_rmse_norm_strm = weights[1]*(
            np.nanmean(
                exp_2_strm_rmse[:-1]/exp_2_strm_obs[:-1]
            )
        )
        exp_2_strm_mae = np.array(exp_2_wat_df.loc['strm_mae'])
        exp_2_mae_norm_strm = weights[1]*(
            np.nanmean(
                exp_2_strm_mae[:-1]/exp_2_strm_mae[:-1]
            )
        )
        exp_2_mae_norm_def_strm = weights[1]*(
            np.nanmean(
                exp_2_strm_mae[:-1]/exp_1_strm_mae[:-1]
            )
        )
        ## for exp 3
        exp_3_rmse_et = weights[0]*(
            exp_3_pix_df['all'].loc['le_rmse']/
            exp_3_pix_df['all'].loc['le_obs']
        )
        exp_3_ubrmse_et = weights[0]*(
            exp_3_pix_df['all'].loc['le_ubrmse']/
            exp_3_pix_df['all'].loc['le_obs']
        )
        exp_3_rmse_strm = weights[1]*(
            exp_3_wat_df['all'].loc['strm_rmse']/
            exp_3_wat_df['all'].loc['strm_obs']
        )
        exp_3_ubrmse_strm = weights[1]*(
            exp_3_wat_df['all'].loc['strm_ubrmse']/
            exp_3_wat_df['all'].loc['strm_obs']
        )
        exp_3_bias_strm = weights[1]*np.abs(
            exp_3_wat_df['all'].loc['strm_avg_diff']/
            exp_3_wat_df['all'].loc['strm_obs']
        )
        exp_3_corr_strm = weights[1]*(
            exp_3_wat_df['all'].loc['strm_corr']
        )
        exp_3_le_ubrmse = np.array(exp_3_pix_df.loc['le_ubrmse'])
        exp_3_le_obs = np.array(exp_3_pix_df.loc['le_obs'])
        exp_3_ubrmse_norm_et = weights[0]*(
            np.nanmean(
                exp_3_le_ubrmse[:-1]/exp_1_le_ubrmse[:-1]
            )
        )
        exp_3_strm_rmse = np.array(exp_3_wat_df.loc['strm_rmse'])
        exp_3_strm_obs = np.array(exp_3_wat_df.loc['strm_obs'])
        exp_3_rmse_norm_strm = weights[1]*(
            np.nanmean(
                exp_3_strm_rmse[:-1]/exp_3_strm_obs[:-1]
            )
        )
        exp_3_strm_mae = np.array(exp_3_wat_df.loc['strm_mae'])
        exp_3_mae_norm_strm = weights[1]*(
            np.nanmean(
                exp_3_strm_mae[:-1]/exp_3_strm_mae[:-1]
            )
        )
        exp_3_mae_norm_def_strm = weights[1]*(
            np.nanmean(
                exp_3_strm_mae[:-1]/exp_1_strm_mae[:-1]
            )
        )
        ## for exp 4
        exp_4_rmse_et = weights[0]*(
            exp_4_pix_df['all'].loc['le_rmse']/
            exp_4_pix_df['all'].loc['le_obs']
        )
        exp_4_ubrmse_et = weights[0]*(
            exp_4_pix_df['all'].loc['le_ubrmse']/
            exp_4_pix_df['all'].loc['le_obs']
        )
        exp_4_rmse_strm = weights[1]*(
            exp_4_wat_df['all'].loc['strm_rmse']/
            exp_4_wat_df['all'].loc['strm_obs']
        )
        exp_4_ubrmse_strm = weights[1]*(
            exp_4_wat_df['all'].loc['strm_ubrmse']/
            exp_4_wat_df['all'].loc['strm_obs']
        )
        exp_4_bias_strm = weights[1]*np.abs(
            exp_4_wat_df['all'].loc['strm_avg_diff']/
            exp_4_wat_df['all'].loc['strm_obs']
        )
        exp_4_corr_strm = weights[1]*(
            exp_4_wat_df['all'].loc['strm_corr']
        )
        exp_4_le_ubrmse = np.array(exp_4_pix_df.loc['le_ubrmse'])
        exp_4_le_obs = np.array(exp_4_pix_df.loc['le_obs'])
        exp_4_ubrmse_norm_et = weights[0]*(
            np.nanmean(
                exp_4_le_ubrmse[:-1]/exp_1_le_ubrmse[:-1]
            )
        )
        exp_4_strm_rmse = np.array(exp_4_wat_df.loc['strm_rmse'])
        exp_4_strm_obs = np.array(exp_4_wat_df.loc['strm_obs'])
        exp_4_rmse_norm_strm = weights[1]*(
            np.nanmean(
                exp_4_strm_rmse[:-1]/exp_4_strm_obs[:-1]
            )
        )
        exp_4_strm_mae = np.array(exp_4_wat_df.loc['strm_mae'])
        exp_4_mae_norm_strm = weights[1]*(
            np.nanmean(
                exp_4_strm_mae[:-1]/exp_4_strm_mae[:-1]
            )
        )
        exp_4_mae_norm_def_strm = weights[1]*(
            np.nanmean(
                exp_4_strm_mae[:-1]/exp_1_strm_mae[:-1]
            )
        )
        # lets make the plot
        et_rmse = [
            exp_1_rmse_et,
            exp_2_rmse_et,
            exp_3_rmse_et,
            exp_4_rmse_et
        ]
        et_ubrmse = [
            exp_1_ubrmse_et,
            exp_2_ubrmse_et,
            exp_3_ubrmse_et,
            exp_4_ubrmse_et
        ]
        strm_rmse = [
            exp_1_rmse_strm,
            exp_2_rmse_strm,
            exp_3_rmse_strm,
            exp_4_rmse_strm
        ]
        strm_ubrmse = [
            exp_1_ubrmse_strm,
            exp_2_ubrmse_strm,
            exp_3_ubrmse_strm,
            exp_4_ubrmse_strm
        ]
        strm_bias = [
            np.abs(exp_1_bias_strm),
            np.abs(exp_2_bias_strm),
            np.abs(exp_3_bias_strm),
            np.abs(exp_4_bias_strm)
        ]
        et_ubrmse_norm = [
            exp_1_ubrmse_norm_et,
            exp_2_ubrmse_norm_et,
            exp_3_ubrmse_norm_et,
            exp_4_ubrmse_norm_et
        ]
        strm_rmse_norm = [
            exp_1_rmse_norm_strm,
            exp_2_rmse_norm_strm,
            exp_3_rmse_norm_strm,
            exp_4_rmse_norm_strm
        ]
        strm_mae_norm = [
            exp_1_mae_norm_strm,
            exp_2_mae_norm_strm,
            exp_3_mae_norm_strm,
            exp_4_mae_norm_strm
        ]
        strm_mae_def_norm = [
            exp_1_mae_norm_def_strm,
            exp_2_mae_norm_def_strm,
            exp_3_mae_norm_def_strm,
            exp_4_mae_norm_def_strm
        ]
        # another way to make the plot
        exp_1_et = [
            exp_1_rmse_et,
            exp_1_ubrmse_et,
            exp_1_ubrmse_et,
            exp_1_ubrmse_et
        ]
        exp_1_strm = [
            exp_1_rmse_strm,
            exp_1_rmse_strm,
            exp_1_bias_strm,
            exp_1_ubrmse_strm
        ]
        exp_2_et = [
            exp_2_rmse_et,
            exp_2_ubrmse_et,
            exp_2_ubrmse_et,
            exp_2_ubrmse_et
        ]
        exp_2_strm = [
            exp_2_rmse_strm,
            exp_2_rmse_strm,
            exp_2_bias_strm,
            exp_2_ubrmse_strm
        ]
        exp_3_et = [
            exp_3_rmse_et,
            exp_3_ubrmse_et,
            exp_3_ubrmse_et,
            exp_3_ubrmse_et
        ]
        exp_3_strm = [
            exp_3_rmse_strm,
            exp_3_rmse_strm,
            exp_3_bias_strm,
            exp_3_ubrmse_strm
        ]
        exp_4_et = [
            exp_4_rmse_et,
            exp_4_ubrmse_et,
            exp_4_ubrmse_et,
            exp_4_ubrmse_et
        ]
        exp_4_strm = [
            exp_4_rmse_strm,
            exp_4_rmse_strm,
            exp_4_bias_strm,
            exp_4_ubrmse_strm
        ]
        error_metrics = [
            'rmse_et_rmse_strm',
            'ubrmse_et_rmse_strm',
            'ubrmse_et_bias_strm',
            'ubrmse_et_ubrmse_strm'
        ]
        # make the plots
        N = len(et_rmse)
        ind = np.arange(N)
        width = 0.35
        # for rmse rmse
        fig = plt.subplots()
        p1 = plt.bar(ind,et_rmse,width)
        p2 = plt.bar(ind,strm_rmse,width,bottom=et_rmse)
        plt.ylabel('J')
        plt.title('J = rmse(ET) + rmse(strm)')
        plt.xticks(
            ind,(
                exp_1_names[0],
                exp_1_names[2],
                exp_2_names[0],
                exp_2_names[2]
            ),rotation=45
        )
        plt.ylim(0,2)
        plt.legend((p1[0],p2[0]),('rmse(ET)','rmse(strm)'))
        save_name = os.path.join(
            plots_dir,'obj_function_rmse_rmse_{}.png'.format(
                'all'
            )
        )
        plt.savefig(save_name,bbox_inches='tight')
        plt.close()
        # for ubrmse rmse
        fig = plt.subplots()
        p1 = plt.bar(ind,et_ubrmse,width)
        p2 = plt.bar(ind,strm_rmse,width,bottom=et_ubrmse)
        plt.ylabel('J')
        plt.title('J = ubrmse(ET) + rmse(strm)')
        plt.xticks(
            ind,(
                exp_1_names[0],
                exp_1_names[2],
                exp_2_names[0],
                exp_2_names[2]
            ),rotation=45
        )
        plt.ylim(0,2)
        plt.legend((p1[0],p2[0]),('ubrmse(ET)','rmse(strm)'))
        save_name = os.path.join(
            plots_dir,'obj_function_ubrmse_rmse_{}.png'.format(
                'all'
            )
        )
        plt.savefig(save_name,bbox_inches='tight')
        plt.close()
        # for ubrmse bias
        fig = plt.subplots()
        p1 = plt.bar(ind,et_ubrmse,width)
        p2 = plt.bar(ind,strm_bias,width,bottom=et_ubrmse)
        plt.ylabel('J')
        plt.title('J = ubrmse(ET) + bias(strm)')
        plt.xticks(
            ind,(
                exp_1_names[0],
                exp_1_names[2],
                exp_2_names[0],
                exp_2_names[2]
            ),rotation=45
        )
        plt.ylim(0,2)
        plt.legend((p1[0],p2[0]),('ubrmse(ET)','bias(strm)'))
        save_name = os.path.join(
            plots_dir,'obj_function_ubrmse_bias_{}.png'.format(
                'all'
            )
        )
        plt.savefig(save_name,bbox_inches='tight')
        plt.close()
        # for ubrmse ubrmse
        fig = plt.subplots()
        p1 = plt.bar(ind,et_ubrmse,width)
        p2 = plt.bar(ind,strm_ubrmse,width,bottom=et_ubrmse)
        plt.ylabel('J')
        plt.title('J = ubrmse(ET) + ubrmse(strm)')
        plt.xticks(
            ind,(
                exp_1_names[0],
                exp_1_names[2],
                exp_2_names[0],
                exp_2_names[2]
            ),rotation=45
        )
        plt.ylim(0,2)
        plt.legend((p1[0],p2[0]),('ubrmse(ET)','ubrmse(strm)'))
        save_name = os.path.join(
            plots_dir,'obj_function_ubrmse_ubrmse_{}.png'.format(
                'all'
            )
        )
        plt.savefig(save_name,bbox_inches='tight')
        plt.close()
        # for ubrmse_norm rmse_norm
        fig = plt.subplots()
        p1 = plt.bar(ind,et_ubrmse_norm,width)
        p2 = plt.bar(ind,strm_rmse_norm,width,bottom=et_ubrmse_norm)
        plt.ylabel('J')
        plt.title('J = ubrmse_norm(ET) + rmse_norm(strm)')
        plt.xticks(
            ind,(
                exp_1_names[0],
                exp_1_names[2],
                exp_2_names[0],
                exp_2_names[2]
            ),rotation=45
        )
        plt.ylim(0,2)
        plt.legend((p1[0],p2[0]),('ubrmse_norm(ET)','rmse_norm(strm)'))
        save_name = os.path.join(
            plots_dir,'obj_function_rmse_norm_ubrmse_norm_{}.png'.format(
                'all'
            )
        )
        plt.savefig(save_name,bbox_inches='tight')
        plt.close()
        # for ubrmse_norm mae_norm
        fig = plt.subplots()
        p1 = plt.bar(ind,et_ubrmse_norm,width)
        p2 = plt.bar(ind,strm_mae_norm,width,bottom=et_ubrmse_norm)
        plt.ylabel('J')
        plt.title('J = ubrmse_norm(ET) + mae_norm(strm)')
        plt.ylim([0,1])
        plt.xticks(
            ind,(
                exp_1_names[0],
                exp_1_names[2],
                exp_2_names[0],
                exp_2_names[2]
            ),rotation=45
        )
        plt.legend((p1[0],p2[0]),('ubrmse_norm(ET)','mae_norm(strm)'))
        save_name = os.path.join(
            plots_dir,'obj_function_mae_norm_ubrmse_norm_{}.png'.format(
                'all'
            )
        )
        plt.savefig(save_name,bbox_inches='tight')
        plt.close()
        # for exp 1
        fig = plt.subplots()
        p1 = plt.bar(ind,exp_1_et,width)
        p2 = plt.bar(ind,exp_1_strm,width,bottom=exp_1_et)
        plt.ylabel('J')
        plt.title(exp_1_names[0])
        plt.xticks(
            ind,(
                error_metrics[0],
                error_metrics[1],
                error_metrics[2],
                error_metrics[3]
            ),rotation=45
        )
        plt.ylim(0,2)
        plt.legend((p1[0],p2[0]),('error(ET)','error(strm)'))
        save_name = os.path.join(
            plots_dir,'obj_function_{}_{}.png'.format(
                exp_1_names[0],'all'
            )
        )
        plt.savefig(save_name,bbox_inches='tight')
        plt.close()
        # for exp 2
        fig = plt.subplots()
        p1 = plt.bar(ind,exp_2_et,width)
        p2 = plt.bar(ind,exp_2_strm,width,bottom=exp_2_et)
        plt.ylabel('J')
        plt.title(exp_1_names[2])
        plt.xticks(
            ind,(
                error_metrics[0],
                error_metrics[1],
                error_metrics[2],
                error_metrics[3]
            ),rotation=45
        )
        plt.ylim(0,2)
        plt.legend((p1[0],p2[0]),('error(ET)','error(strm)'))
        save_name = os.path.join(
            plots_dir,'obj_function_{}_{}.png'.format(
                exp_1_names[2],'all'
            )
        )
        plt.savefig(save_name,bbox_inches='tight')
        plt.close()
        # for exp 3
        fig = plt.subplots()
        p1 = plt.bar(ind,exp_3_et,width)
        p2 = plt.bar(ind,exp_3_strm,width,bottom=exp_3_et)
        plt.ylabel('J')
        plt.title(exp_2_names[0])
        plt.xticks(
            ind,(
                error_metrics[0],
                error_metrics[1],
                error_metrics[2],
                error_metrics[3]
            ),rotation=45
        )
        plt.ylim(0,2)
        plt.legend((p1[0],p2[0]),('error(ET)','error(strm)'))
        save_name = os.path.join(
            plots_dir,'obj_function_{}_{}.png'.format(
                exp_2_names[0],'all'
            )
        )
        plt.savefig(save_name,bbox_inches='tight')
        plt.close()
        # for exp 4
        fig = plt.subplots()
        p1 = plt.bar(ind,exp_4_et,width)
        p2 = plt.bar(ind,exp_4_strm,width,bottom=exp_4_et)
        plt.ylabel('J')
        plt.title(exp_2_names[2])
        plt.xticks(
            ind,(
                error_metrics[0],
                error_metrics[1],
                error_metrics[2],
                error_metrics[3]
            ),rotation=45
        )
        plt.ylim(0,2)
        plt.legend((p1[0],p2[0]),('error(ET)','error(strm)'))
        save_name = os.path.join(
            plots_dir,'obj_function_{}_{}.png'.format(
                exp_2_names[2],'all'
            )
        )
        plt.savefig(save_name,bbox_inches='tight')
        plt.close()
        # now for specific watersheds
        selected_sheds = [
            12040500,12189500,2314500,6360500
        ]
        for s,shed in enumerate(selected_sheds):
            exp_1_rmse_et = weights[0]*(
                exp_1_wat_df[shed].loc['le_rmse']/
                exp_1_wat_df[shed].loc['le']
            )
            exp_1_ubrmse_et = weights[0]*(
                exp_1_wat_df[shed].loc['le_ubrmse']/
                exp_1_wat_df[shed].loc['le']
            )
            exp_1_rmse_strm = weights[1]*(
                exp_1_wat_df[shed].loc['strm_rmse']/
                exp_1_wat_df[shed].loc['strm']
            )
            exp_1_ubrmse_strm = weights[1]*(
                exp_1_wat_df[shed].loc['strm_ubrmse']/
                exp_1_wat_df[shed].loc['strm']
            )
            exp_1_bias_strm = weights[1]*np.abs(
                exp_1_wat_df[shed].loc['strm_avg_diff']/
                exp_1_wat_df[shed].loc['strm']
            )
            exp_1_corr_strm = weights[1]*(
                exp_1_wat_df[shed].loc['strm_corr']
            )
            ## for exp 2
            exp_2_rmse_et = weights[0]*(
                exp_2_wat_df[shed].loc['le_rmse']/
                exp_2_wat_df[shed].loc['le']
            )
            exp_2_ubrmse_et = weights[0]*(
                exp_2_wat_df[shed].loc['le_ubrmse']/
                exp_2_wat_df[shed].loc['le']
            )
            exp_2_rmse_strm = weights[1]*(
                exp_2_wat_df[shed].loc['strm_rmse']/
                exp_2_wat_df[shed].loc['strm']
            )
            exp_2_ubrmse_strm = weights[1]*(
                exp_2_wat_df[shed].loc['strm_ubrmse']/
                exp_2_wat_df[shed].loc['strm']
            )
            exp_2_bias_strm = weights[1]*np.abs(
                exp_2_wat_df[shed].loc['strm_avg_diff']/
                exp_2_wat_df[shed].loc['strm']
            )
            exp_2_corr_strm = weights[1]*(
                exp_2_wat_df[shed].loc['strm_corr']
            )
            ## for exp 3
            exp_3_rmse_et = weights[0]*(
                exp_3_wat_df[shed].loc['le_rmse']/
                exp_3_wat_df[shed].loc['le']
            )
            exp_3_ubrmse_et = weights[0]*(
                exp_3_wat_df[shed].loc['le_ubrmse']/
                exp_3_wat_df[shed].loc['le']
            )
            exp_3_rmse_strm = weights[1]*(
                exp_3_wat_df[shed].loc['strm_rmse']/
                exp_3_wat_df[shed].loc['strm']
            )
            exp_3_ubrmse_strm = weights[1]*(
                exp_3_wat_df[shed].loc['strm_ubrmse']/
                exp_3_wat_df[shed].loc['strm']
            )
            exp_3_bias_strm = weights[1]*np.abs(
                exp_3_wat_df[shed].loc['strm_avg_diff']/
                exp_3_wat_df[shed].loc['strm']
            )
            exp_3_corr_strm = weights[1]*(
                exp_3_wat_df[shed].loc['strm_corr']
            )
            ## for exp 4
            exp_4_rmse_et = weights[0]*(
                exp_4_wat_df[shed].loc['le_rmse']/
                exp_4_wat_df[shed].loc['le']
            )
            exp_4_ubrmse_et = weights[0]*(
                exp_4_wat_df[shed].loc['le_ubrmse']/
                exp_4_wat_df[shed].loc['le']
            )
            exp_4_rmse_strm = weights[1]*(
                exp_4_wat_df[shed].loc['strm_rmse']/
                exp_4_wat_df[shed].loc['strm']
            )
            exp_4_ubrmse_strm = weights[1]*(
                exp_4_wat_df[shed].loc['strm_ubrmse']/
                exp_4_wat_df[shed].loc['strm']
            )
            exp_4_bias_strm = weights[1]*np.abs(
                exp_4_wat_df[shed].loc['strm_avg_diff']/
                exp_4_wat_df[shed].loc['strm']
            )
            exp_4_corr_strm = weights[1]*(
                exp_4_wat_df[shed].loc['strm_corr']
            )
            # lets make the plot
            et_rmse = [
                exp_1_rmse_et,
                exp_2_rmse_et,
                exp_3_rmse_et,
                exp_4_rmse_et
            ]
            et_ubrmse = [
                exp_1_ubrmse_et,
                exp_2_ubrmse_et,
                exp_3_ubrmse_et,
                exp_4_ubrmse_et
            ]
            strm_rmse = [
                exp_1_rmse_strm,
                exp_2_rmse_strm,
                exp_3_rmse_strm,
                exp_4_rmse_strm
            ]
            strm_ubrmse = [
                exp_1_ubrmse_strm,
                exp_2_ubrmse_strm,
                exp_3_ubrmse_strm,
                exp_4_ubrmse_strm
            ]
            strm_bias = [
                np.abs(exp_1_bias_strm),
                np.abs(exp_2_bias_strm),
                np.abs(exp_3_bias_strm),
                np.abs(exp_4_bias_strm)
            ]
            # another way to make the plot
            exp_1_et = [
                exp_1_rmse_et,
                exp_1_ubrmse_et,
                exp_1_ubrmse_et,
                exp_1_ubrmse_et
            ]
            exp_1_strm = [
                exp_1_rmse_strm,
                exp_1_rmse_strm,
                exp_1_bias_strm,
                exp_1_ubrmse_strm
            ]
            exp_2_et = [
                exp_2_rmse_et,
                exp_2_ubrmse_et,
                exp_2_ubrmse_et,
                exp_2_ubrmse_et
            ]
            exp_2_strm = [
                exp_2_rmse_strm,
                exp_2_rmse_strm,
                exp_2_bias_strm,
                exp_2_ubrmse_strm
            ]
            exp_3_et = [
                exp_3_rmse_et,
                exp_3_ubrmse_et,
                exp_3_ubrmse_et,
                exp_3_ubrmse_et
            ]
            exp_3_strm = [
                exp_3_rmse_strm,
                exp_3_rmse_strm,
                exp_3_bias_strm,
                exp_3_ubrmse_strm
            ]
            exp_4_et = [
                exp_4_rmse_et,
                exp_4_ubrmse_et,
                exp_4_ubrmse_et,
                exp_4_ubrmse_et
            ]
            exp_4_strm = [
                exp_4_rmse_strm,
                exp_4_rmse_strm,
                exp_4_bias_strm,
                exp_4_ubrmse_strm
            ]
            error_metrics = [
                'rmse_et_rmse_strm',
                'ubrmse_et_rmse_strm',
                'ubrmse_et_bias_strm',
                'ubrmse_et_ubrmse_strm'
            ]
            # make the plots
            N = len(et_rmse)
            ind = np.arange(N)
            width = 0.35
            # for rmse rmse
            fig = plt.subplots()
            p1 = plt.bar(ind,et_rmse,width)
            p2 = plt.bar(ind,strm_rmse,width,bottom=et_rmse)
            plt.ylabel('J')
            plt.title('J = rmse(ET) + rmse(strm)')
            plt.xticks(
                ind,(
                    exp_1_names[0],
                    exp_1_names[2],
                    exp_2_names[0],
                    exp_2_names[2]
                ),rotation=45
            )
            plt.ylim(0,2)
            plt.legend((p1[0],p2[0]),('rmse(ET)','rmse(strm)'))
            save_name = os.path.join(
                plots_dir,'obj_function_rmse_rmse_{}.png'.format(
                    shed
                )
            )
            plt.savefig(save_name,bbox_inches='tight')
            plt.close()
            # for ubrmse rmse
            fig = plt.subplots()
            p1 = plt.bar(ind,et_ubrmse,width)
            p2 = plt.bar(ind,strm_rmse,width,bottom=et_ubrmse)
            plt.ylabel('J')
            plt.title('J = ubrmse(ET) + rmse(strm)')
            plt.xticks(
                ind,(
                    exp_1_names[0],
                    exp_1_names[2],
                    exp_2_names[0],
                    exp_2_names[2]
                ),rotation=45
            )
            plt.ylim(0,2)
            plt.legend((p1[0],p2[0]),('ubrmse(ET)','rmse(strm)'))
            save_name = os.path.join(
                plots_dir,'obj_function_ubrmse_rmse_{}.png'.format(
                    shed
                )
            )
            plt.savefig(save_name,bbox_inches='tight')
            plt.close()
            # for ubrmse bias
            fig = plt.subplots()
            p1 = plt.bar(ind,et_ubrmse,width)
            p2 = plt.bar(ind,strm_bias,width,bottom=et_ubrmse)
            plt.ylabel('J')
            plt.title('J = ubrmse(ET) + bias(strm)')
            plt.xticks(
                ind,(
                    exp_1_names[0],
                    exp_1_names[2],
                    exp_2_names[0],
                    exp_2_names[2]
                ),rotation=45
            )
            plt.ylim(0,2)
            plt.legend((p1[0],p2[0]),('ubrmse(ET)','bias(strm)'))
            save_name = os.path.join(
                plots_dir,'obj_function_ubrmse_bias_{}.png'.format(
                    shed
                )
            )
            plt.savefig(save_name,bbox_inches='tight')
            plt.close()
            # for ubrmse ubrmse
            fig = plt.subplots()
            p1 = plt.bar(ind,et_ubrmse,width)
            p2 = plt.bar(ind,strm_ubrmse,width,bottom=et_ubrmse)
            plt.ylabel('J')
            plt.title('J = ubrmse(ET) + ubrmse(strm)')
            plt.xticks(
                ind,(
                    exp_1_names[0],
                    exp_1_names[2],
                    exp_2_names[0],
                    exp_2_names[2]
                ),rotation=45
            )
            plt.ylim(0,2)
            plt.legend((p1[0],p2[0]),('ubrmse(ET)','ubrmse(strm)'))
            save_name = os.path.join(
                plots_dir,'obj_function_ubrmse_ubrmse_{}.png'.format(
                    shed
                )
            )
            plt.savefig(save_name,bbox_inches='tight')
            plt.close()
            # for exp 1
            fig = plt.subplots()
            p1 = plt.bar(ind,exp_1_et,width)
            p2 = plt.bar(ind,exp_1_strm,width,bottom=exp_1_et)
            plt.ylabel('J')
            plt.title(exp_1_names[0])
            plt.xticks(
                ind,(
                    error_metrics[0],
                    error_metrics[1],
                    error_metrics[2],
                    error_metrics[3]
                ),rotation=45
            )
            plt.ylim(0,2)
            plt.legend((p1[0],p2[0]),('error(ET)','error(strm)'))
            save_name = os.path.join(
                plots_dir,'obj_function_{}_{}.png'.format(
                    exp_1_names[0],shed
                )
            )
            plt.savefig(save_name,bbox_inches='tight')
            plt.close()
            # for exp 2
            fig = plt.subplots()
            p1 = plt.bar(ind,exp_2_et,width)
            p2 = plt.bar(ind,exp_2_strm,width,bottom=exp_2_et)
            plt.ylabel('J')
            plt.title(exp_1_names[2])
            plt.xticks(
                ind,(
                    error_metrics[0],
                    error_metrics[1],
                    error_metrics[2],
                    error_metrics[3]
                ),rotation=45
            )
            plt.ylim(0,2)
            plt.legend((p1[0],p2[0]),('error(ET)','error(strm)'))
            save_name = os.path.join(
                plots_dir,'obj_function_{}_{}.png'.format(
                    exp_1_names[2],shed
                )
            )
            plt.savefig(save_name,bbox_inches='tight')
            plt.close()
            # for exp 3
            fig = plt.subplots()
            p1 = plt.bar(ind,exp_3_et,width)
            p2 = plt.bar(ind,exp_3_strm,width,bottom=exp_3_et)
            plt.ylabel('J')
            plt.title(exp_2_names[0])
            plt.xticks(
                ind,(
                    error_metrics[0],
                    error_metrics[1],
                    error_metrics[2],
                    error_metrics[3]
                ),rotation=45
            )
            plt.ylim(0,2)
            plt.legend((p1[0],p2[0]),('error(ET)','error(strm)'))
            save_name = os.path.join(
                plots_dir,'obj_function_{}_{}.png'.format(
                    exp_2_names[0],shed
                )
            )
            plt.savefig(save_name,bbox_inches='tight')
            plt.close()
            # for exp 4
            fig = plt.subplots()
            p1 = plt.bar(ind,exp_4_et,width)
            p2 = plt.bar(ind,exp_4_strm,width,bottom=exp_4_et)
            plt.ylabel('J')
            plt.title(exp_2_names[2])
            plt.xticks(
                ind,(
                    error_metrics[0],
                    error_metrics[1],
                    error_metrics[2],
                    error_metrics[3]
                ),rotation=45
            )
            plt.ylim(0,2)
            plt.legend((p1[0],p2[0]),('error(ET)','error(strm)'))
            save_name = os.path.join(
                plots_dir,'obj_function_{}_{}.png'.format(
                    exp_2_names[2],shed
                )
            )
            plt.savefig(save_name,bbox_inches='tight')
            plt.close()

    def plot_diff(self,exp_1_names,exp_2_names,exp_1_pix_pso_df,exp_2_pix_pso_df,
                  exp_1_wat_pso_df,exp_2_wat_pso_df,geojson_fname,states_shp,
                  plots_dir,skinny_plot):
        # lets make some watershed scale difference plots
        # perc change in le rmse between experiments
        exp_perc_change_le_rmse = (
            exp_1_pix_pso_df.loc['perc_change_le_rmse'] -
            exp_2_pix_pso_df.loc['perc_change_le_rmse']
        )
        avg_exp_perc_change_le_rmse = exp_perc_change_le_rmse['all']
        # change in le rmse between experiments
        exp_change_le_rmse = (
            exp_1_pix_pso_df.loc['change_le_rmse'] -
            exp_2_pix_pso_df.loc['change_le_rmse']
        )
        avg_exp_change_le_rmse = exp_change_le_rmse['all']
        # change in le r2 between experiments
        exp_change_le_r2 = (
            exp_1_pix_pso_df.loc['change_le_r2'] -
            exp_2_pix_pso_df.loc['change_le_r2']
        )
        avg_exp_change_le_r2 = exp_change_le_r2['all']
        # change in corr between experiments
        exp_change_le_corr = (
            exp_1_pix_pso_df.loc['change_le_corr'] -
            exp_2_pix_pso_df.loc['change_le_corr']
        )
        avg_exp_change_le_corr = exp_change_le_corr['all']
        # perc change in le ubrmse between experiments
        exp_perc_change_le_ubrmse = (
            exp_1_pix_pso_df.loc['perc_change_le_ubrmse'] -
            exp_2_pix_pso_df.loc['perc_change_le_ubrmse']
        )
        avg_exp_perc_change_le_ubrmse = exp_perc_change_le_ubrmse['all']
        # change in le ubrmse between experiments
        exp_change_le_ubrmse = (
            exp_1_pix_pso_df.loc['change_le_ubrmse'] -
            exp_2_pix_pso_df.loc['change_le_ubrmse']
        )
        avg_exp_change_le_ubrmse = exp_change_le_ubrmse['all']
        # difference of average le
        exp_diff_ave_le = (
            exp_1_pix_pso_df.loc['ave_le'] -
            exp_2_pix_pso_df.loc['ave_le']
        )
        avg_exp_diff_ave_le = exp_diff_ave_le['all']
        # difference in normalized le ubrmse
        exp_change_le_ubrmse_norm = (
            exp_1_pix_pso_df.loc['change_le_ubrmse_norm'] -
            exp_2_pix_pso_df.loc['change_le_ubrmse_norm']
        )
        avg_exp_change_le_ubrmse_norm = (
            exp_change_le_ubrmse_norm['all']
        )
        # put values to be plotted into list for plotting
        vals = [
            exp_diff_ave_le,exp_change_le_ubrmse,
            exp_change_le_ubrmse_norm
        ]
        # put the averages that correspond to these values
        avgs = [
            avg_exp_diff_ave_le,avg_exp_change_le_ubrmse,
            avg_exp_change_le_ubrmse_norm
        ]
        # put the name the corresponds to each of these values
        names = [
            'exp_change_ave_le','exp_change_le_ubrmse',
            'exp_change_le_ubrmse_norm'
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
        types = names
        cmaps = {
            'exp_change_ave_le':'bwr',
            'exp_change_le_ubrmse':'bwr',
            'exp_change_le_ubrmse_norm':'bwr'
        }
        vmins = {
            'exp_change_ave_le':-10,
            'exp_change_le_ubrmse':-10,
            'exp_change_le_ubrmse_norm':-.25
        }
        vmaxs = {
            'exp_change_ave_le':10,
            'exp_change_le_ubrmse':10,
            'exp_change_le_ubrmse_norm':.25
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
            lons = exp_1_pix_pso_df.loc['lon']
            lats = exp_1_pix_pso_df.loc['lat']
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
            savename = '{name}_{exp_1}_vs_{exp_2}_optimization_pixels.png'.format(
                name=names[p],exp_1=exp_1_names[2],exp_2=exp_2_names[2]
            )
            savename = os.path.join(
                plots_dir,savename
            )
            plt.savefig(savename,dpi=300,bbox_inches='tight')
            plt.close()
        # now lets do this for our streamflow friends
        # percent change in strm rmse
        # change in strm perc rmse between experiments
        exp_perc_change_strm_rmse = (
            exp_1_wat_pso_df.loc['perc_change_strm_rmse'] -
            exp_2_wat_pso_df.loc['perc_change_strm_rmse']
        )
        avg_exp_perc_change_strm_rmse = exp_perc_change_strm_rmse['all']
        exp_perc_change_strm_rmse = np.array(exp_perc_change_strm_rmse)
        exp_perc_change_strm_rmse = exp_perc_change_strm_rmse[:-1]
        # change in strm rmse between experiments
        exp_change_strm_rmse = (
            exp_1_wat_pso_df.loc['change_strm_rmse'] -
            exp_2_wat_pso_df.loc['change_strm_rmse']
        )
        avg_exp_change_strm_rmse = exp_change_strm_rmse['all']
        exp_change_strm_rmse = np.array(exp_change_strm_rmse)
        exp_change_strm_rmse = exp_change_strm_rmse[:-1]
        # change in normalized strm rmse between experiments
        exp_change_strm_rmse_norm = (
            exp_1_wat_pso_df.loc['change_strm_rmse_norm'] -
            exp_2_wat_pso_df.loc['change_strm_rmse_norm']
        )
        avg_exp_change_strm_rmse_norm = exp_change_strm_rmse_norm['all']
        exp_change_strm_rmse_norm = np.array(exp_change_strm_rmse_norm)
        exp_change_strm_rmse_norm = exp_change_strm_rmse_norm[:-1]
        # change in individualyk normalized strm rmse between experimetns
        exp_change_strm_rmse_norm_ind = (
            exp_1_wat_pso_df.loc['change_strm_rmse_norm_ind'] -
            exp_2_wat_pso_df.loc['change_strm_rmse_norm_ind']
        )
        avg_exp_change_strm_rmse_norm_ind = (
            exp_change_strm_rmse_norm_ind['all']
        )
        exp_change_strm_rmse_norm_ind = np.array(
            exp_change_strm_rmse_norm_ind
        )
        exp_change_strm_rmse_norm_ind = (
            exp_change_strm_rmse_norm_ind[:-1]
        )
        # change in strm r2 between experiments
        #exp_change_strm_r2 = (
        #    exp_1_wat_pso_df.loc['change_strm_r2'] -
        #    exp_2_wat_pso_df.loc['change_strm_r2']
        #)
        #avg_exp_change_strm_r2 = exp_change_strm_r2['all']
        #exp_change_strm_r2 = np.array(exp_change_strm_r2)
        #exp_change_strm_r2 = exp_change_strm_r2[:-1]
        ## change in strm corr between experiments
        #exp_change_strm_corr = (
        #    exp_1_wat_pso_df.loc['change_strm_corr'] -
        #    exp_2_wat_pso_df.loc['change_strm_corr']
        #)
        #avg_exp_change_strm_corr = exp_change_strm_corr['all']
        #exp_change_strm_corr = np.array(exp_change_strm_corr)
        #exp_change_strm_corr = exp_change_strm_corr[:-1]
        ## change in strm perc ubrmse between experiments
        #exp_perc_change_strm_ubrmse = (
        #    exp_1_wat_pso_df.loc['perc_change_strm_ubrmse'] -
        #    exp_2_wat_pso_df.loc['perc_change_strm_ubrmse']
        #)
        #avg_exp_perc_change_strm_ubrmse = exp_perc_change_strm_ubrmse['all']
        #exp_perc_change_strm_ubrmse = np.array(exp_perc_change_strm_ubrmse)
        #exp_perc_change_strm_ubrmse = exp_perc_change_strm_ubrmse[:-1]
        # change in strm nse between experiments
        exp_change_strm_nse = (
            exp_1_wat_pso_df.loc['change_strm_nse'] -
            exp_2_wat_pso_df.loc['change_strm_nse']
        )
        avg_exp_change_strm_nse = exp_change_strm_nse['all']
        exp_change_strm_nse = np.array(exp_change_strm_nse)
        exp_change_strm_nse = exp_change_strm_nse[:-1]
        # change in strm nse between experiments
        exp_change_strm_nse_ind = (
            exp_1_wat_pso_df.loc['change_strm_nse_ind'] -
            exp_2_wat_pso_df.loc['change_strm_nse_ind']
        )
        avg_exp_change_strm_nse_ind = exp_change_strm_nse_ind['all']
        exp_change_strm_nse_ind = np.array(exp_change_strm_nse_ind)
        exp_change_strm_nse_ind = exp_change_strm_nse_ind[:-1]
        # change in individaully averaged nse between experiments
        # change in strm between experiments
        exp_diff_strm = (
            exp_1_wat_pso_df.loc['diff_strm'] -
            exp_2_wat_pso_df.loc['diff_strm']
        )
        avg_exp_diff_strm = exp_diff_strm['all']
        exp_diff_strm = np.array(exp_diff_strm)
        exp_diff_strm = exp_diff_strm[:-1]
        # change in normalized strm rmse between experiments
        exp_change_strm_rmse_norm = (
            exp_1_wat_pso_df.loc['diff_strm_rmse_norm'] -
            exp_2_wat_pso_df.loc['diff_strm_rmse_norm']
        )
        avg_exp_change_strm_rmse_norm = exp_change_strm_rmse_norm['all']
        exp_change_strm_rmse_norm = np.array(exp_change_strm_rmse_norm)
        exp_change_strm_rmse_norm = exp_change_strm_rmse_norm[:-1]
        # and plot this data
        # now let's get the shapes that we need for plotting
        huc6s = gpd.read_file(geojson_fname)
        # now let's get everythin in arrays for proper plotting
        names = [
            'exp_diff_strm',
            'exp_change_strm_rmse_norm'
        ]
        vals = [
            exp_diff_strm,
            exp_change_strm_rmse_norm_ind
        ]
        avgs = [
            avg_exp_diff_strm,
            avg_exp_change_strm_rmse_norm_ind
        ]
        types = names
        cmaps = {
            'exp_change_strm_rmse':'bwr',
            'exp_change_strm_rmse_norm':'bwr',
            'exp_change_strm_rmse_norm_ind':'bwr',
            'exp_change_strm_nse':'bwr',
            'exp_change_strm_nse_ind':'bwr',
            'exp_diff_strm':'bwr'
        }
        vmins = {
            'exp_change_strm_rmse':-.4,
            'exp_change_strm_rmse_norm':-.4,
            'exp_change_strm_rmse_norm_ind':-1,
            'exp_change_strm_nse':-1,
            'exp_change_strm_nse_ind':-1,
            'exp_diff_strm':-1
        }
        vmaxs = {
            'exp_change_strm_rmse':.4,
            'exp_change_strm_rmse_norm':.4,
            'exp_change_strm_rmse_norm_ind':1,
            'exp_change_strm_nse':1,
            'exp_change_strm_nse_ind':1,
            'exp_diff_strm':1
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
        all_hucs = np.array(exp_1_wat_pso_df.columns)
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
                    ax.fill(xs,ys,fc='none',ec='k',linewidth=.2)
                except:
                    for geom in this_geom.geoms:
                        xs,ys = geom.exterior.xy
                        ax.fill(xs,ys,fc='none',ec='k',linewidth=.2)
            ## get a list of all the hucs
            #all_hucs = list(huc6s['huc6'])
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
                '{name}_{exp_1}_vs_{exp_2}_huc6.png'.format(
                    name=names[n],exp_1=exp_1_names[2],exp_2=exp_2_names[2]
                )
            )
            plt.savefig(this_savename,dpi=350,bbox_inches='tight')
            plt.close()
