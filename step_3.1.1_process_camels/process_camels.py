import pandas as pd
import numpy as np
import sys
import os
import copy
import datetime
import geopandas as gpd
from dateutil.relativedelta import relativedelta

class camels:
    def __init__(self,start,end):
        self.start = start
        self.end  = end
    def get_camels(self,camels_fname):
        # get the tiles that we are going to run as defined in step 1
        camels = pd.read_csv(camels_fname)
        camels = camels.set_index('camel')
        return camels
    def get_streamflow(self,camels_fname,chosen_camels_df):
        # let's get the camels that we are looking for
        chosen_camels = np.array(chosen_camels_df.index)
        cols = chosen_camels
        all_dates = []
        curr = copy.deepcopy(self.start)
        while curr <= self.end:
            all_dates.append(curr.strftime('%Y%m%d'))
            curr += relativedelta(years=1)
        camels_truth = pd.DataFrame(columns=cols,index=all_dates)
        camels_truth.index.name = 'time'
        data_cols = [
            'gauge','year','mon','day','strm','qc'
        ]
        for c,cam in enumerate(chosen_camels):
            print(
                'processing camel {} out of {}'.format(
                    c,len(chosen_camels)
                )
            )
            camel_str = str(cam).zfill(8)
            file_found = False
            dir_guess = 0
            while not file_found:
                dir_guess_str = str(dir_guess).zfill(2)
                this_fname = os.path.join(
                    camels_fname,
                    dir_guess_str,
                    '{}_streamflow_qc.txt'.format(
                        camel_str
                    )
                )
                file_found = os.path.isfile(this_fname)
                dir_guess += 1
            this_data = pd.read_csv(
                this_fname,
                delim_whitespace=True,
                header=None,
                names=data_cols
            )
            start_of_year = copy.deepcopy(self.start)
            end_of_year = (
                start_of_year +
                relativedelta(years=1) -
                relativedelta(days=1)
            )
            yearly_avg = np.zeros(0)
            while end_of_year <= self.end:
                beginning_aligns = True
                try:
                    start_idx = np.where(
                        (this_data['year'] == start_of_year.year) &
                        (this_data['mon'] == start_of_year.month) &
                        (this_data['day'] == start_of_year.day)
                    )[0][0]
                except:
                    beginning_aligns = False
                if not beginning_aligns:
                    print(this_data)
                    diff = end_of_year - start_of_year
                    days_diff = diff.days + 1
                    this_strm = np.repeat(np.nan,days_diff)
                else:
                    try:
                        end_idx = np.where(
                            (this_data['year'] == end_of_year.year) &
                            (this_data['mon'] == end_of_year.month) &
                            (this_data['day'] == end_of_year.day)
                        )[0][0]
                        missing_end = False
                    except:
                        print(this_data)
                        data_end_date = datetime.date(
                            this_data['year'].iloc[-1],
                            this_data['mon'].iloc[-1],
                            this_data['day'].iloc[-1]
                        )
                        end_idx = np.where(
                            (this_data['year'] == data_end_date.year) &
                            (this_data['mon'] == data_end_date.month) &
                            (this_data['day'] == data_end_date.day)
                        )[0][0]
                        diff = end_of_year - data_end_date
                        days_diff = diff.days
                        missing_end = True
                    this_data_trim = this_data.iloc[
                        start_idx:end_idx+1
                    ]
                    this_strm = np.array(this_data_trim['strm'])
                    nan_idx = np.where(this_strm == -999)
                    this_strm[nan_idx] = np.nan
                    if missing_end:
                        for d in range(days_diff):
                            this_strm = np.append(this_strm,np.nan)
                num_days_missing = 0
                for s in this_strm:
                    if np.isnan(s) == True:
                        num_days_missing += 1
                if num_days_missing > 30:
                    this_year_avg = np.nan
                else:
                    this_year_avg = np.nanmean(this_strm)
                yearly_avg = np.append(yearly_avg,this_year_avg)
                start_of_year += relativedelta(years=1)
                end_of_year += relativedelta(years=1)
            # lets convert
            this_area = chosen_camels_df['area'].loc[cam]
            # strm data oringally in cf/s
            yearly_avg = yearly_avg/(3280.84**3) # km3/s
            yearly_avg = yearly_avg/this_area # km/s
            yearly_avg = yearly_avg*1e6 # mm/s
            yearly_avg = yearly_avg*86400 # mm/d
            camels_truth[cam] = yearly_avg
            camels_truth = camels_truth.copy()
        return camels_truth
    def save_truth(self,camels_truth,out_dir):
        fname = os.path.join(
            out_dir,
            'camels_truth_yearly_{st}_{end}_mm_day.csv'.format(
                st = self.start,
                end = self.end
            )
        )
        print(camels_truth)
        camels_truth.to_csv(fname)












