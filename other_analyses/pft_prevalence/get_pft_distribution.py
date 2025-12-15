import pandas as pd
import numpy as np
import pickle
import sys

class pft_dist:
    def get_info(self,fname):
        data = pd.read_csv(fname)
        data = data.set_index('tile')
        return data
    def get_intersection_info(self,fname):
        with open(fname,'rb') as f:
            data = pickle.load(f)
        return data
    def get_pso_pfts(self,pft_data,intersection_info):
        tiles = np.array(pft_data.index)
        our_pfts = [
            'Needleleaf',
            'Broadleaf',
            'Shrub',
            'Cool c3 grass',
            'Warm c4 grass',
            'Crop'
        ]
        pft_converter = {
            'Needleleaf evergreen temperate tree':our_pfts[0],
            'Needleleaf evergreen boreal tree':our_pfts[0],
            'Needleleaf deciduous boreal tree':our_pfts[0],
            'Broadleaf evergreen tropical tree':our_pfts[1],
            'Broadleaf evergreen temperate tree':our_pfts[1],
            'Broadleaf deciduous tropical tree':our_pfts[1],
            'Broadleaf deciduous temperate tree':our_pfts[1],
            'Broadleaf deciduous boreal tree':our_pfts[1],
            'Broadleaf evergreen temperate shrub':our_pfts[2],
            'Broadleaf deciduous temperate shrub':our_pfts[2],
            'Broadleaf deciduous temperate shrub[moisture stress only]':our_pfts[2],
            'Broadleaf deciduous boreal shrub':our_pfts[2],
            'Arctic c3 grass':our_pfts[3],
            'Cool c3 grass':our_pfts[3],
            'Cool c3 grass [moisture stress only]':our_pfts[3],
            'Warm c4 grass':our_pfts[4],
            'Warm c4 grass [moisture stress only]':our_pfts[4],
            'Crop':our_pfts[5],
            'Crop [moisture stress only]':our_pfts[5],
            '(Corn)':our_pfts[5],
            '(Irrigated corn)':our_pfts[5],
            '(Spring temperate cereal)':our_pfts[5],
            '(Irrigated spring temperate cereal)':our_pfts[5],
            '(winter temperate cereal)':our_pfts[5],
            '(Irrigated winter temperate cereal)':our_pfts[5],
            '(Soybean)':our_pfts[5],
            '(Irrigated Soybean)':our_pfts[5],
        }
        pft_1_perc = np.array(pft_data['pft_1_perc'])
        pft_2_perc = np.array(pft_data['pft_2_perc'])
        pft_3_perc = np.array(pft_data['pft_3_perc'])
        pft_4_perc = np.array(pft_data['pft_4_perc'])
        pft_percentages = [
            pft_1_perc,
            pft_2_perc,
            pft_3_perc,
            pft_4_perc
        ]
        their_pft_1 = np.array(pft_data['pft_1_name'])
        their_pft_2 = np.array(pft_data['pft_2_name'])
        their_pft_3 = np.array(pft_data['pft_3_name'])
        their_pft_4 = np.array(pft_data['pft_4_name'])
        our_pft_df = pd.DataFrame(columns=['num','avg_perc'])
        pft_80_pix_df = pd.DataFrame(columns=['num','pixels'])
        pft_85_pix_df = pd.DataFrame(columns=['num','pixels'])
        pft_90_pix_df = pd.DataFrame(columns=['num','pixels'])
        pft_95_pix_df = pd.DataFrame(columns=['num','pixels'])
        for p,pft in enumerate(our_pfts):
            print(pft)
            pft_exists = np.zeros(len(tiles))
            pft_perc = np.zeros(len(tiles))
            pft_80 = 0
            pft_80_pixels = np.zeros(0)
            pft_85 = 0
            pft_85_pixels = np.zeros(0)
            pft_90 = 0
            pft_90_pixels = np.zeros(0)
            pft_95 = 0
            pft_95_pixels = np.zeros(0)
            for t,ti in enumerate(tiles):
                these_pfts = [
                    their_pft_1[t],
                    their_pft_2[t],
                    their_pft_3[t],
                    their_pft_4[t]
                ]
                for tp,these in enumerate(these_pfts):
                    this_our_pft = pft_converter[these]
                    if this_our_pft == pft:
                        pft_exists[t] = 1
                        pft_perc[t] += pft_percentages[tp][t]
                if pft_perc[t] >= 80:
                    pft_80 += 1
                    pft_80_pixels = np.append(pft_80_pixels,ti)
                if pft_perc[t] >= 85:
                    pft_85 += 1
                    pft_85_pixels = np.append(pft_85_pixels,ti)
                if pft_perc[t] >= 90:
                    pft_90 += 1
                    pft_90_pixels = np.append(pft_90_pixels,ti)
                if pft_perc[t] >= 95:
                    pft_95 += 1
                    pft_95_pixels = np.append(pft_95_pixels,ti)
            total_num = np.sum(pft_exists)
            average_perc_idx = np.where(pft_perc != 0)
            average_perc = np.average(pft_perc[average_perc_idx])
            our_pft_df.loc[pft] = [total_num,average_perc]
            pft_80_pix_df.loc[pft] = [pft_80,pft_80_pixels]
            pft_85_pix_df.loc[pft] = [pft_85,pft_85_pixels]
            pft_90_pix_df.loc[pft] = [pft_90,pft_90_pixels]
            pft_95_pix_df.loc[pft] = [pft_95,pft_95_pixels]
        print('all')
        print(our_pft_df)
        print('pix_80')
        print(pft_80_pix_df)
        watersheds = list(intersection_info.keys())
        pft_80_wat_df = pd.DataFrame(columns=['num','watersheds'])
        pft_85_wat_df = pd.DataFrame(columns=['num','watersheds'])
        pft_90_wat_df = pd.DataFrame(columns=['num','watersheds'])
        pft_95_wat_df = pd.DataFrame(columns=['num','watersheds'])
        for p,pft in enumerate(our_pfts):
            avg_in_wat = np.zeros(len(watersheds))
            num_above_80 = 0
            wat_above_80 = np.zeros(0)
            for w,wat in enumerate(watersheds):
                this_tiles = intersection_info[wat][0]
                this_perc = intersection_info[wat][1]
                pft_perc = np.zeros(len(this_tiles))
                for t,ti in enumerate(this_tiles):
                    this_idx = np.where(tiles == ti)[0][0]
                    these_pfts = [
                        their_pft_1[this_idx],
                        their_pft_2[this_idx],
                        their_pft_3[this_idx],
                        their_pft_4[this_idx]
                    ]
                    for tp,these in enumerate(these_pfts):
                        this_our_pft = pft_converter[these]
                        if this_our_pft == pft:
                            pft_perc[t] += pft_percentages[tp][this_idx]
                this_avg_in_wat = np.average(
                    pft_perc,weights=this_perc
                )
                avg_in_wat[w] = this_avg_in_wat
                if this_avg_in_wat > 80:
                    num_above_80 += 1
                    wat_above_80 = np.append(wat_above_80,wat)
            pft_80_wat_df.loc[pft] = [num_above_80,wat_above_80]
        print('wat_80')
        print(pft_80_wat_df)
        to_return = {
            'all':our_pft_df,
            'pix_80':pft_80_pix_df,
            'wat_80':pft_80_wat_df
        }
        return to_return







