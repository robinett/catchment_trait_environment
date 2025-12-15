import sys
import spei
import scipy.stats as scs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

class drought:
    def __init__(self):
        pass
    def find_drought(self,fluxcom_timeseries,timeseries_w_precip,roll_len,
                     plots_dir,plot_spei,save_spei,load_spei,
                     spei_fname,is_drought_fname):
        if not load_spei:
            precip_kg_m2_s = timeseries_w_precip['rainfsnowf']
            precip_mm_day = precip_kg_m2_s*86400
            cols_all = list(fluxcom_timeseries.columns)
            idx_all = list(precip_kg_m2_s.index)
            idx_spei = idx_all[roll_len:]
            pixels_all = cols_all[1:]
            # create the dataframes where we will store things
            spei_df = pd.DataFrame(index=idx_all)
            is_drought_df = pd.DataFrame(index=idx_all)
            for p,pix in enumerate(pixels_all):
                print(pix)
                this_precip_mm_day = precip_mm_day[pix]
                correct_index = list(this_precip_mm_day.index)
                this_et_w_m2 = fluxcom_timeseries[pix]
                this_et_mm_day = this_et_w_m2/28.94
                this_et_mm_day.index = correct_index
                p_minus_et = this_precip_mm_day - this_et_mm_day
                if not np.isnan(p_minus_et[0]):
                    spei_calc = spei.spei(p_minus_et,dist=scs.fisk)
                    spei_calc = spei_calc.rolling(roll_len).mean()
                    if plot_spei:
                        tmin,tmax = pd.to_datetime(['2006-03-01','2006-12-31'])
                        plt.figure(figsize=(8,4))
                        spei.plot.si(spei_calc)
                        plt.xlim(tmin,tmax)
                        plt.grid()
                        plt.title('SPEI for pixel {}'.format(pix))
                        savename = os.path.join(
                            plots_dir,
                            'spei_pix_{}_smooth_{}.png'.format(pix,roll_len)
                        )
                        plt.savefig(savename)
                        plt.close()
                    cons_days = 0
                    is_drought = np.repeat(False,len(spei_calc))
                    for s,sp in enumerate(spei_calc):
                        if sp < -1 and not np.isnan(sp):
                            cons_days += 1
                        else:
                            cons_days = 0
                        if cons_days >= 90:
                            is_drought[s-89:s+1] = True
                else:
                    spei_calc = np.repeat(np.nan,len(p_minus_et))
                    is_drought = np.repeat(np.nan,len(p_minus_et))
                spei_df[pix] = spei_calc
                spei_df = spei_df.copy()
                is_drought_df[pix] = is_drought
                is_drought_df = is_drought_df.copy()
            if save_spei:
                spei_df.to_csv(spei_fname)
                is_drought_df.to_csv(is_drought_fname)
        elif load_spei:
            spei_df = pd.read_csv(spei_fname)
            is_drought_df = pd.read_csv(is_drought_fname)
        return [spei_df,is_drought_df]
    def trim_to_drought(self,data_df,is_drought_df):
        print(is_drought_df)
        print(data_df)
