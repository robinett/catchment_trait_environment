import sys
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/pso/step_5.1_analyze_outputs/funcs'
)
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large'
)
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/pso/step_2x_env_covariates'
)
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/pso/step_5_analyze_outputs'
)
import os
import datetime
import numpy as np
import pickle as pickle
import pandas as pd
import netCDF4 as nc
import statsmodels.formula.api as sm
import copy
import matplotlib as mpl
import xarray as xr
from tqdm import tqdm
import matplotlib.colors as mcolors
from get_env_covariates import get_covariates
from extrapolation import extrapolate
from general_functions import gen_funcs
from choose_tiles_camels import choose
from get_timeseries import get_timeseries
from get_averages_and_error import averages_and_error
from plot_watersheds import plot_wat
from plot_other import plot_other
from plot_pixels import plot_pixels
from analyze_timeseries_pixel_scale import analyze_pix

def main():
    # where are we?
    base_dir = '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate'
    # where do we put outputs?
    out_dir = os.path.join(
        base_dir,
        'outputs'
    )
    # where do we put plots?
    plots_dir = os.path.join(
        base_dir,
        'plots'
    )
    # where is the data
    data_dir = os.path.join(
        base_dir,
        'data'
    )
    # where is the optimization PFT info?
    opt_pft_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_5_analyze_outputs/' +
        'outputs/pft_distribution.csv'
    )
    # where is the extrapolation PFT info?
    ext_pft_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/' +
        'outputs/processed_pft_info_extrap.csv'
    )
    # where is the potimization canopy height information?
    opt_canopy_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_2x_env_covariates/outputs/' +
        'canopy_height_normalized.nc4'
    )
    # where is the extrapolation canopy height info?
    ext_canopy_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/' +
        'outputs/canopy_height_normalized.nc4'
    )
    # where is the optimization MAP information?
    opt_precip_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_2x_env_covariates/outputs/' +
        'gldas_avg_precip_normalized.nc4'
    )
    # where is the extrapolation canopy height info?
    ext_precip_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/' +
        'outputs/gldas_avg_precip_normalized.nc4'
    )
    # where is the optimization pixel information located?
    # where are the analysis pixels located?
    opt_pixels_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_1x_choose_tiles_large' +
        '/outputs/intersecting_catch_tiles.csv'
    )
    # where is the extrapolation pixel information located?
    ext_pixels_fname = (
        '/discover/nobackup/trobinet/from_aws/pso/step_6x_extrapolate/' +
        'outputs/extrapolation_pix.csv'
    )
    positions_pft_fname = (
        '/discover/nobackup/trobinet/GEOSldas_pso_final_pft/' +
        'GEOSldas/positions.csv'
    )
    positions_ef_fname = (
        '/discover/nobackup/trobinet/GEOSldas_pso_final_ef/' +
        'GEOSldas/positions.csv'
    )
    gen = gen_funcs()
    # load the optimization pixels
    opt_pixels = gen.get_pixels(opt_pixels_fname)
    # load the extrapolation
    ext_pixels = gen.get_pixels(ext_pixels_fname)
    # merge to get all CONUS pixels
    conus_pixels = np.append(opt_pixels,ext_pixels)
    conus_pixels = np.sort(conus_pixels)
    #u,c = np.unique(conus_pixels,return_counts=True)
    #dups = u[c>1]
    # merge PFT info
    opt_pft = pd.read_csv(opt_pft_fname)
    opt_pft = opt_pft.set_index('tile')
    # select PFT info for all CONUS pixels
    ext_pft = pd.read_csv(ext_pft_fname)
    ext_pft = ext_pft.set_index('tile')
    # merge these two
    conus_pft = pd.concat([opt_pft,ext_pft])
    conus_pft = conus_pft[~conus_pft.index.duplicated(keep='first')]
    # get opt environmental coavriates
    opt_h = xr.open_dataset(opt_canopy_fname)
    opt_h = opt_h.rename({'vals':'h'})
    opt_map = xr.open_dataset(opt_precip_fname)
    opt_map = opt_map.rename({'vals':'map'})
    opt_env = xr.merge([opt_h,opt_map])
    # get ext environmental covariates
    ext_h = xr.open_dataset(ext_canopy_fname)
    ext_h = ext_h.rename({'vals':'h'})
    ext_map = xr.open_dataset(ext_precip_fname)
    ext_map = ext_map.rename({'vals':'map'})
    ext_env = xr.merge([ext_h,ext_map])
    # get conus environmental covariates by merging the two
    conus_env = xr.concat([opt_env,ext_env],dim='tile_num')
    # get the total pft area for both opt and conus
    # first for opt
    simple_map = {
        'Needleleaf evergreen temperate tree':'needleleaf_trees',
        'Needleleaf evergreen boreal tree':'needleleaf_trees',
        'Needleleaf deciduous boreal tree':'broadleaf_trees',
        'Broadleaf evergreen tropical tree':'broaleaf_trees',
        'Broadleaf evergreen temperate tree':'broadleaf_trees',
        'Broadleaf deciduous tropical tree':'broadleaf_trees',
        'Broadleaf deciduous temperate tree':'broadleaf_trees',
        'Broadleaf deciduous boreal tree':'broadleaf_trees',
        'Broadleaf evergreen temperate shrub':'shrub',
        'Broadleaf deciduous boreal shrub':'shrub',
        'Broadleaf deciduous temperate shrub':'shrub',
        'Broadleaf deciduous temperate shrub[moisture +deciduous]': (
            'shrub'
        ),
        'Broadleaf deciduous temperate shrub[moisture stress only]': (
            'shrub'
        ),
        'Arctic c3 grass':'c3_grass',
        'Cool c3 grass':'c3_grass',
        'Cool c3 grass [moisture + deciduous]':'c3_grass',
        'Cool c3 grass [moisture stress only]':'c3_grass',
        'Warm c4 grass [moisture + deciduous]':'c4_grass',
        'Warm c4 grass [moisture stress only]':'c4_grass',
        'Warm c4 grass':'c4_grass',
        'Crop':'crop',
        'Crop [moisture + deciduous]':'crop',
        'Crop [moisture stress only]':'crop',
        '(Spring temperate cereal)':'crop',
        '(Irrigated corn)':'crop',
        '(Soybean)':'crop',
        '(Corn)':'crop',
        '(Irrigated spring temperate cereal)':'crop'
    }
    # overall areas for each PFT
    opt_need = 0
    opt_broad = 0
    opt_shrub = 0
    opt_c3 = 0
    opt_c4 = 0
    opt_crop = 0
    # climate distributions for each PFT
    opt_map_vals = opt_env['map'].values
    opt_h_vals = opt_env['h'].values
    opt_env_pix = opt_env['tile'].values
    opt_need_h = np.zeros(0)
    opt_need_map = np.zeros(0)
    opt_broad_h = np.zeros(0)
    opt_broad_map = np.zeros(0)
    opt_shrub_h = np.zeros(0)
    opt_shrub_map = np.zeros(0)
    opt_c3_h = np.zeros(0)
    opt_c3_map = np.zeros(0)
    opt_c4_h = np.zeros(0)
    opt_c4_map = np.zeros(0)
    opt_crop_h = np.zeros(0)
    opt_crop_map = np.zeros(0)
    for p,pix in enumerate(opt_pixels):
        for i in range(4):
            if i == 0:
                used_need = False
                used_broad = False
                used_shrub = False
                used_c3 = False
                used_c4 = False
                used_crop = False
            # add up the percentages of each PFT in that pixels
            # also get the MAP and h for each PFT that exists in that pixel
            j = i+1
            this_pft = opt_pft['pft_{}_name'.format(j)].loc[pix]
            this_pft_simple = simple_map[this_pft]
            this_perc = opt_pft['pft_{}_perc'.format(j)].loc[pix]
            this_env_idx = np.where(opt_env_pix == pix)[0][0]
            this_h = opt_h_vals[this_env_idx]
            this_map = opt_map_vals[this_env_idx]
            if this_pft_simple == 'needleleaf_trees' and this_perc > 0:
                opt_need += this_perc
                if used_need == False:
                    opt_need_h = np.append(opt_need_h,this_h)
                    opt_need_map = np.append(opt_need_map,this_map)
                    used_need = True
            elif this_pft_simple == 'broadleaf_trees' and this_perc > 0:
                opt_broad += this_perc
                if used_broad == False:
                    opt_broad_h = np.append(opt_broad_h,this_h)
                    opt_broad_map = np.append(opt_broad_map,this_map)
                    used_broad = True
            elif this_pft_simple == 'shrub' and this_perc > 0:
                opt_shrub += this_perc
                if used_shrub == False:
                    opt_shrub_h = np.append(opt_shrub_h,this_h)
                    opt_shrub_map = np.append(opt_shrub_map,this_map)
                    used_shrub = True
            elif this_pft_simple == 'c3_grass' and this_perc > 0:
                opt_c3 += this_perc
                if used_c3 == False:
                    opt_c3_h = np.append(opt_c3_h,this_h)
                    opt_c3_map = np.append(opt_c3_map,this_map)
                    used_c3 = True
            elif this_pft_simple == 'c4_grass' and this_perc > 0:
                opt_c4 += this_perc
                if used_c4 == False:
                    opt_c4_h = np.append(opt_c4_h,this_h)
                    opt_c4_map = np.append(opt_c4_map,this_map)
                    used_c4 = True
            elif this_pft_simple == 'crop' and this_perc > 0:
                opt_crop += this_perc
                if used_crop == False:
                    opt_crop_h = np.append(opt_crop_h,this_h)
                    opt_crop_map = np.append(opt_crop_map,this_map)
                    used_crop == True
    total_perc = np.sum([
        opt_need,
        opt_broad,
        opt_shrub,
        opt_c3,
        opt_c4,
        opt_crop
    ])
    opt_need_perc = opt_need/total_perc
    opt_broad_perc = opt_broad/total_perc
    opt_shrub_perc = opt_shrub/total_perc
    opt_c3_perc = opt_c3/total_perc
    opt_c4_perc = opt_c4/total_perc
    opt_crop_perc = opt_crop/total_perc
    # make a pie chart of the perc for each 
    p_o = plot_other()
    # let's get the colors we want
    ## Define start and end colors
    #start_color = "#8fe9ff"
    #end_color = "#614ba5"
    ## Generate 6 evenly spaced colors
    #colors = mcolors.LinearSegmentedColormap.from_list(
    #    "custom_gradient", [start_color, end_color]
    #)
    #color_list = [
    #    mcolors.to_hex(colors(i)) for i in np.linspace(0, 1, 6)
    #]
    #def reorder_list(original_list):
    #    from random import shuffle
    #    # Initialize variables
    #    n = len(original_list)
    #    if n <= 1:
    #        return original_list
    #    # Track adjacency
    #    adjacency = set()
    #    for i in range(n - 1):
    #        adjacency.add((original_list[i], original_list[i + 1]))
    #        adjacency.add((original_list[i + 1], original_list[i]))
    #    # Shuffle to start with
    #    shuffled_list = original_list[:]
    #    shuffle(shuffled_list)
    #    # Try to reorder
    #    def is_valid(new_list):
    #        for i in range(len(new_list) - 1):
    #            if (new_list[i], new_list[i + 1]) in adjacency:
    #                return False
    #        return True
    #    # Perform reordering
    #    for _ in range(1000):  # Number of attempts to find a valid permutation
    #        shuffle(shuffled_list)
    #        if is_valid(shuffled_list):
    #            return shuffled_list
    #    # If no valid permutation found, return the shuffled list
    #    return shuffled_list
    #color_list = reorder_list(color_list)
    color_list = [
        '#8e5e96',
        '#86c36c',
        '#cb9ac2',
        '#75c7ed',
        '#050708',
        '#ef3f29'
    ]
    p_o.plot_pie(
        ['need','broad','shrub','c3','c4','crop'],
        [
            opt_need_perc,opt_broad_perc,opt_shrub_perc,
            opt_c3_perc,opt_c4_perc,opt_crop_perc
        ],
        color_list,
        os.path.join(
            plots_dir,
            'opt_pie.svg'
        )
    )
    # now for conus
    conus_need = 0
    conus_broad = 0
    conus_shrub = 0
    conus_c3 = 0
    conus_c4 = 0
    conus_crop = 0
    # climate distributions for each PFT
    conus_map_vals = conus_env['map'].values
    conus_h_vals = conus_env['h'].values
    conus_env_pix = conus_env['tile'].values
    conus_need_h = np.zeros(0)
    conus_need_map = np.zeros(0)
    conus_broad_h = np.zeros(0)
    conus_broad_map = np.zeros(0)
    conus_shrub_h = np.zeros(0)
    conus_shrub_map = np.zeros(0)
    conus_c3_h = np.zeros(0)
    conus_c3_map = np.zeros(0)
    conus_c4_h = np.zeros(0)
    conus_c4_map = np.zeros(0)
    conus_crop_h = np.zeros(0)
    conus_crop_map = np.zeros(0)
    for p,pix in enumerate(conus_pixels):
        for i in range(4):
            if i == 0:
                used_need = False
                used_broad = False
                used_shrub = False
                used_c3 = False
                used_c4 = False
                used_crop = False
            # add up the percentages of each PFT in that pixels
            # also get the MAP and h for each PFT that exists in that pixel
            j = i+1
            this_pft = conus_pft['pft_{}_name'.format(j)].loc[pix]
            this_pft_simple = simple_map[this_pft]
            this_perc = conus_pft['pft_{}_perc'.format(j)].loc[pix]
            this_env_idx = np.where(conus_env_pix == pix)[0][0]
            this_h = conus_h_vals[this_env_idx]
            this_map = conus_map_vals[this_env_idx]
            if this_pft_simple == 'needleleaf_trees' and this_perc > 0:
                conus_need += this_perc
                if used_need == False:
                    conus_need_h = np.append(conus_need_h,this_h)
                    conus_need_map = np.append(conus_need_map,this_map)
                    used_need = True
            elif this_pft_simple == 'broadleaf_trees' and this_perc > 0:
                conus_broad += this_perc
                if used_broad == False:
                    conus_broad_h = np.append(conus_broad_h,this_h)
                    conus_broad_map = np.append(conus_broad_map,this_map)
                    used_broad = True
            elif this_pft_simple == 'shrub' and this_perc > 0:
                conus_shrub += this_perc
                if used_shrub == False:
                    conus_shrub_h = np.append(conus_shrub_h,this_h)
                    conus_shrub_map = np.append(conus_shrub_map,this_map)
                    used_shrub = True
            elif this_pft_simple == 'c3_grass' and this_perc > 0:
                conus_c3 += this_perc
                if used_c3 == False:
                    conus_c3_h = np.append(conus_c3_h,this_h)
                    conus_c3_map = np.append(conus_c3_map,this_map)
                    used_c3 = True
            elif this_pft_simple == 'c4_grass' and this_perc > 0:
                conus_c4 += this_perc
                if used_c4 == False:
                    conus_c4_h = np.append(conus_c4_h,this_h)
                    conus_c4_map = np.append(conus_c4_map,this_map)
                    used_c4 = True
            elif this_pft_simple == 'crop' and this_perc > 0:
                conus_crop += this_perc
                if used_crop == False:
                    conus_crop_h = np.append(conus_crop_h,this_h)
                    conus_crop_map = np.append(conus_crop_map,this_map)
                    used_crop == True
    total_perc = np.sum([
        conus_need,
        conus_broad,
        conus_shrub,
        conus_c3,
        conus_c4,
        conus_crop
    ])
    conus_need_perc = conus_need/total_perc
    conus_broad_perc = conus_broad/total_perc
    conus_shrub_perc = conus_shrub/total_perc
    conus_c3_perc = conus_c3/total_perc
    conus_c4_perc = conus_c4/total_perc
    conus_crop_perc = conus_crop/total_perc
    # make a pie chart of the perc for each 
    p_o = plot_other()
    p_o.plot_pie(
        ['need','broad','shrub','c3','c4','crop'],
        [
            conus_need_perc,conus_broad_perc,conus_shrub_perc,
            conus_c3_perc,conus_c4_perc,conus_crop_perc
        ],
        color_list,
        os.path.join(
            plots_dir,
            'conus_pie.svg'
        )
    )
<<<<<<< HEAD
    # plot this as a double bar too since it might be a better way
    # to show it
    p_o.multibar_plot(
        [
            [
                opt_need_perc,
                opt_broad_perc,
                opt_shrub_perc,
                opt_c3_perc,
                opt_c4_perc,
                opt_crop_perc
            ],
            [
                conus_need_perc,
                conus_broad_perc,
                conus_shrub_perc,
                conus_c3_perc,
                conus_c4_perc,
                conus_crop_perc
            ]
        ],
        ['need','broad','shrub','c3','c4','crop'],
        ['opt','conus'],
        plots_dir,
        'landcover_rep_doublebar.svg',
        'PFT',
        '% of area',
        bar_colors=['#8e5e96','#75c7ed']
    )
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    # and extrapolation only for later
    ext_need = 0
    ext_broad = 0
    ext_shrub = 0
    ext_c3 = 0
    ext_c4 = 0
    ext_crop = 0
    # climate distributions for each PFT
    ext_map_vals = ext_env['map'].values
    ext_h_vals = ext_env['h'].values
    ext_env_pix = ext_env['tile'].values
    ext_need_h = np.zeros(0)
    ext_need_map = np.zeros(0)
    ext_broad_h = np.zeros(0)
    ext_broad_map = np.zeros(0)
    ext_shrub_h = np.zeros(0)
    ext_shrub_map = np.zeros(0)
    ext_c3_h = np.zeros(0)
    ext_c3_map = np.zeros(0)
    ext_c4_h = np.zeros(0)
    ext_c4_map = np.zeros(0)
    ext_crop_h = np.zeros(0)
    ext_crop_map = np.zeros(0)
    for p,pix in enumerate(ext_pixels):
        for i in range(4):
            if i == 0:
                used_need = False
                used_broad = False
                used_shrub = False
                used_c3 = False
                used_c4 = False
                used_crop = False
            # add up the percentages of each PFT in that pixels
            # also get the MAP and h for each PFT that exists in that pixel
            j = i+1
            this_pft = ext_pft['pft_{}_name'.format(j)].loc[pix]
            this_pft_simple = simple_map[this_pft]
            this_perc = ext_pft['pft_{}_perc'.format(j)].loc[pix]
            this_env_idx = np.where(ext_env_pix == pix)[0][0]
            this_h = ext_h_vals[this_env_idx]
            this_map = ext_map_vals[this_env_idx]
            if this_pft_simple == 'needleleaf_trees' and this_perc > 0:
                ext_need += this_perc
                if used_need == False:
                    ext_need_h = np.append(ext_need_h,this_h)
                    ext_need_map = np.append(ext_need_map,this_map)
                    used_need = True
            elif this_pft_simple == 'broadleaf_trees' and this_perc > 0:
                ext_broad += this_perc
                if used_broad == False:
                    ext_broad_h = np.append(ext_broad_h,this_h)
                    ext_broad_map = np.append(ext_broad_map,this_map)
                    used_broad = True
            elif this_pft_simple == 'shrub' and this_perc > 0:
                ext_shrub += this_perc
                if used_shrub == False:
                    ext_shrub_h = np.append(ext_shrub_h,this_h)
                    ext_shrub_map = np.append(ext_shrub_map,this_map)
                    used_shrub = True
            elif this_pft_simple == 'c3_grass' and this_perc > 0:
                ext_c3 += this_perc
                if used_c3 == False:
                    ext_c3_h = np.append(ext_c3_h,this_h)
                    ext_c3_map = np.append(ext_c3_map,this_map)
                    used_c3 = True
            elif this_pft_simple == 'c4_grass' and this_perc > 0:
                ext_c4 += this_perc
                if used_c4 == False:
                    ext_c4_h = np.append(ext_c4_h,this_h)
                    ext_c4_map = np.append(ext_c4_map,this_map)
                    used_c4 = True
            elif this_pft_simple == 'crop' and this_perc > 0:
                ext_crop += this_perc
                if used_crop == False:
                    ext_crop_h = np.append(ext_crop_h,this_h)
                    ext_crop_map = np.append(ext_crop_map,this_map)
                    used_crop == True
    # plot the pdf for MAP and h for both optimization and CONUS
    # for each PFT
    # for need
<<<<<<< HEAD
    line_styles = ['-','--','-','--']
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
    all_need_vals = [
        opt_need_map,conus_need_map,
        opt_need_h,conus_need_h
    ]
    labels = [
        'opt_need_map','conus_need_map',
        'opt_need_h','conus_need_h'
    ]
    colors = [
        '#854e8d','#71c6ed',
        '#ef3f29','#f0a861'
    ]
    kernel_widths = [1,1,1,1]
    line_widths = [1.5,1.5,1.5,1.5]
    fills = [False,False,False,False]
    save_name = os.path.join(
        plots_dir,'need_env_pdf.svg'
    )
    xlim = [-3,4]
    p_o.plot_pdf(
        all_need_vals,
        labels,
        colors,
        kernel_widths,
        line_widths,
        fills,
<<<<<<< HEAD
        line_styles,
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
        save_name,
        xlim=xlim
    )
    # for broad
    all_broad_vals = [
        opt_broad_map,conus_broad_map,
        opt_broad_h,conus_broad_h
    ]
    labels = [
        'opt_broad_map','conus_broad_map',
        'opt_broad_h','conus_broad_h'
    ]
    colors = [
        '#854e8d','#71c6ed',
        '#ef3f29','#f0a861'
    ]
    kernel_widths = [1,1,1,1]
    line_widths = [1.5,1.5,1.5,1.5]
    fills = [False,False,False,False]
    save_name = os.path.join(
        plots_dir,'broad_env_pdf.svg'
    )
    p_o.plot_pdf(
        all_broad_vals,
        labels,
        colors,
        kernel_widths,
        line_widths,
        fills,
<<<<<<< HEAD
        line_styles,
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
        save_name,
        xlim=xlim
    )
    # for shrub
    all_shrub_vals = [
        opt_shrub_map,conus_shrub_map,
        opt_shrub_h,conus_shrub_h
    ]
    labels = [
        'opt_shrub_map','conus_shrub_map',
        'opt_shrub_h','conus_shrub_h'
    ]
    colors = [
        '#854e8d','#71c6ed',
        '#ef3f29','#f0a861'
    ]
    kernel_widths = [1,1,1,1]
    line_widths = [1.5,1.5,1.5,1.5]
    fills = [False,False,False,False]
    save_name = os.path.join(
        plots_dir,'shrub_env_pdf.svg'
    )
    p_o.plot_pdf(
        all_shrub_vals,
        labels,
        colors,
        kernel_widths,
        line_widths,
        fills,
<<<<<<< HEAD
        line_styles,
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
        save_name,
        xlim=xlim
    )
    # for c3
    all_c3_vals = [
        opt_c3_map,conus_c3_map,
        opt_c3_h,conus_c3_h
    ]
    labels = [
        'opt_c3_map','conus_c3_map',
        'opt_c3_h','conus_c3_h'
    ]
    colors = [
        '#854e8d','#71c6ed',
        '#ef3f29','#f0a861'
    ]
    kernel_widths = [1,1,1,1]
    line_widths = [1.5,1.5,1.5,1.5]
    fills = [False,False,False,False]
    save_name = os.path.join(
        plots_dir,'c3_env_pdf.svg'
    )
    p_o.plot_pdf(
        all_c3_vals,
        labels,
        colors,
        kernel_widths,
        line_widths,
        fills,
<<<<<<< HEAD
        line_styles,
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
        save_name,
        xlim=xlim
    )
    # for c4
    all_c4_vals = [
        opt_c4_map,conus_c4_map,
        opt_c4_h,conus_c4_h
    ]
    labels = [
        'opt_c4_map','conus_c4_map',
        'opt_c4_h','conus_c4_h'
    ]
    colors = [
        '#854e8d','#71c6ed',
        '#ef3f29','#f0a861'
    ]
    kernel_widths = [1,1,1,1]
    line_widths = [1.5,1.5,1.5,1.5]
    fills = [False,False,False,False]
    save_name = os.path.join(
        plots_dir,'c4_env_pdf.svg'
    )
    p_o.plot_pdf(
        all_c4_vals,
        labels,
        colors,
        kernel_widths,
        line_widths,
        fills,
<<<<<<< HEAD
        line_styles,
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
        save_name,
        xlim=xlim
    )
    # for crop
    all_crop_vals = [
        opt_crop_map,conus_crop_map,
        opt_crop_h,conus_crop_h
    ]
    labels = [
        'opt_crop_map','conus_crop_map',
        'opt_crop_h','conus_crop_h'
    ]
    colors = [
        '#854e8d','#71c6ed',
        '#ef3f29','#f0a861'
    ]
    kernel_widths = [1,1,1,1]
    line_widths = [1.5,1.5,1.5,1.5]
    fills = [False,False,False,False]
    save_name = os.path.join(
        plots_dir,'crop_env_pdf.svg'
    )
    p_o.plot_pdf(
        all_crop_vals,
        labels,
        colors,
        kernel_widths,
        line_widths,
        fills,
<<<<<<< HEAD
        line_styles,
=======
>>>>>>> 7585e74 (adding rest of updates from transition back to Discover)
        save_name,
        xlim=xlim
    )
    ## plot the pdf for MAP for both optimization and CONUS
    #opt_map_vals = opt_env['map'].values
    #conus_map_vals = conus_env['map'].values
    #all_map_vals = [opt_map_vals,conus_map_vals]
    #labels = ['opt','conus']
    #colors = ['#854e8d','#71c6ed']
    #kernel_widths = [1,1]
    #line_widths = [1.5,1.5]
    #save_name = os.path.join(
    #    plots_dir,'map_pdf_100.svg'
    #)
    #p_o.plot_pdf(
    #    all_map_vals,
    #    labels,
    #    colors,
    #    kernel_widths,
    #    line_widths,
    #    save_name
    #)
    ## plot the pdf for h for both optimizations
    #opt_h_vals = opt_env['h'].values
    #conus_h_vals = conus_env['h'].values
    #all_h_vals = [opt_h_vals,conus_h_vals]
    #save_name = os.path.join(
    #    plots_dir,'h_pdf_100.svg'
    #)
    #p_o.plot_pdf(
    #    all_h_vals,
    #    labels,
    #    colors,
    #    kernel_widths,
    #    line_widths,
    #    save_name
    #)
    ## get g1 by landcover for both opt and extrap
    #a_pix = analyze_pix()
    #g = get_timeseries()
    ##pft_info_opt = a_pix.get_pft_info(raw_pft_fname,processed_pft_fname,opt_pixels)
    ## start with positions
    pos_pft = pd.read_csv(positions_pft_fname,header=None)
    pos_pft = np.array(pos_pft)[0]
    ## now put it in the weird format so that our g1 decoder can accept it
    #position_pft_names = [
    #    'a0_needleaf_trees',
    #    'a0_broadleaf_trees',
    #    'a0_shrub',
    #    'a0_c3_grass',
    #    'a0_c4_grass',
    #    'a0_crop'
    #]
    #pos_pft_final = {}
    #for p,pft in enumerate(position_pft_names):
    #    this_pos = pos_pft[p]
    #    this_pos = [this_pos]
    #    this_pos = [this_pos]
    #    pos_pft_final[pft] = this_pos
    ## now we need a fake objective value to pass in this weird format
    #fake_obj = {}
    #fake_obj['all'] = [[1]]
    #g1_pft_opt_df = g.get_g1_map(
    #    pos_pft_final,
    #    fake_obj,
    #    opt_precip_fname,
    #    opt_canopy_fname,
    #    opt_pft_fname,
    #    'pft',0
    #)
    #g1_pft_opt = g1_pft_opt_df.loc['g1']
    ## now let's do this for ef
    pos_ef = pd.read_csv(positions_ef_fname,header=None)
    pos_ef = np.array(pos_ef)[0]
    ## now put it in the weird format so that our g1 decoder can accept it
    #position_ef_names = [
    #    'b_needleleaf_trees',
    #    'b_broadleaf_trees',
    #    'b_shrub',
    #    'b_c3_grass',
    #    'b_c4_grass',
    #    'b_crop',
    #    'a0_intercept',
    #    'a1_precip_coef',
    #    'a2_canopy_coef'
    #]
    #pos_ef_final = {}
    #for p,pft in enumerate(position_ef_names):
    #    this_pos = pos_ef[p]
    #    this_pos = [this_pos]
    #    this_pos = [this_pos]
    #    pos_ef_final[pft] = this_pos
    #g1_ef_opt_df = g.get_g1_map(
    #    pos_ef_final,
    #    fake_obj,
    #    opt_precip_fname,
    #    opt_canopy_fname,
    #    opt_pft_fname,
    #    'ef',0
    #)
    #g1_ef_opt = g1_ef_opt_df.loc['g1']
    ## make ext g1 pdf that also has line for opt
    ## start with positions
    ## now we need a fake objective value to pass in this weird format
    #g1_pft_ext_df = g.get_g1_map(
    #    pos_pft_final,
    #    fake_obj,
    #    ext_precip_fname,
    #    ext_canopy_fname,
    #    ext_pft_fname,
    #    'pft',0
    #)
    #g1_pft_ext = g1_pft_ext_df.loc['g1']
    ## now let's do this for ef
    #g1_ef_ext_df = g.get_g1_map(
    #    pos_ef_final,
    #    fake_obj,
    #    ext_precip_fname,
    #    ext_canopy_fname,
    #    ext_pft_fname,
    #    'ef',0
    #)
    #g1_ef_ext = g1_ef_ext_df.loc['g1']
    # make the extrapolation figure
    #g1_by_pft_opt = {
    #    'needleleaf_trees':pd.DataFrame(
    #        {'type':['pft','ef','perc']}
    #    ).set_index('type'),
    #    'broadleaf_trees':pd.DataFrame(
    #        {'type':['pft','ef','perc']}
    #    ).set_index('type'),
    #    'shrub':pd.DataFrame(
    #        {'type':['pft','ef','perc']}
    #    ).set_index('type'),
    #    'c3_grass':pd.DataFrame(
    #        {'type':['pft','ef','perc']}
    #    ).set_index('type'),
    #    'c4_grass':pd.DataFrame(
    #        {'type':['pft','ef','perc']}
    #    ).set_index('type'),
    #    'crop':pd.DataFrame(
    #        {'type':['pft','ef','perc']}
    #    ).set_index('type'),
    #}
    #print(g1_by_pft_opt)
    ## for each pixel
    #tile_pft_info_opt = pd.read_csv(opt_pft_fname)
    #tile_pft_info_opt = tile_pft_info_opt.set_index('tile')
    #tile_pft_info_opt_df = pd.DataFrame(columns=opt_pixels)
    #tile_pft_info_opt_df['pfts'] = [
    #    'needleleaf_trees','broadleaf_trees','shrub','c3_grass','c4_grass',
    #    'crop'
    #]
    #tile_pft_info_opt_df = tile_pft_info_opt_df.set_index('pfts')
    #for p,pix in enumerate(opt_pixels):
    #    for pft in range(4):
    #        p = pft + 1
    #        this_pft = tile_pft_info_opt['pft_{}_name'.format(p)].loc[pix]
    #        this_pft_simple = simple_map[this_pft]
    #        this_perc = tile_pft_info_opt['pft_{}_perc'.format(p)].loc[pix]
    #        if pd.isnull(tile_pft_info_opt_df[pix].loc[this_pft_simple]) == True:
    #            tile_pft_info_opt_df[pix].loc[this_pft_simple] = this_perc
    #        else:
    #            tile_pft_info_opt_df[pix].loc[this_pft_simple] += this_perc
    #for p,pix in enumerate(opt_pixels):
    #    # for each pft
    #    for p,pft in enumerate(tile_pft_info_opt_df.index.to_list()):
    #        # if not nan, calculate ef-based g1 and pft-based g1 for
    #        # that pft at that pixel
    #        this_val = tile_pft_info_opt_df[pix].loc[pft]
    #        if pd.isnull(this_val) == False:
    #            # add to the pft-based series or g1-based series
    #            this_pft_g1 = pos_pft[pft]
    #            this_ef_g1 = (
    #                positions_ef[pft]*(
    #                    pos_ef['b0'] +
    #                    positions_ef['b1']*precip_df[pix].loc['precip'] +
    #                    positions_ef['b2']*canopy_df[pix].loc['canopy']
    #                )
    #            )
    #            this_pft_g1 = list(this_pft_g1)[0]
    #            this_ef_g1  = list(this_ef_g1)[0]
    #            if this_pft_g1 < 0.5:
    #                this_pft_g1 = 0.5
    #            if this_ef_g1 < 0.5:
    #                this_ef_g1 = 0.5
    #            comb_g1 = [this_pft_g1,this_ef_g1,this_val]
    #            g1_by_pft[pft][pix] = comb_g1
    pfts = [
        'needleleaf_trees','broadleaf_trees','shrub','c3_grass','c4_grass',
        'crop'
    ]
    default_g1s = {
        'needleleaf_trees':2.3,
        'broadleaf_trees':4.25,
        'shrub':4.7,
        'c3_grass':5.3,
        'c4_grass':1.6,
        'crop':5.79
    }
    lin_g1s = {
        'needleleaf_trees':2.35,
        'broadleaf_trees':3.97,
        'shrub':4.22,
        'c4_grass':1.62,
        'c3_grass':4.5,
        'crop':5.79
    }
    pos_pft = pd.read_csv(positions_pft_fname,header=None)
    pos_pft = np.array(pos_pft)[0]
    pos_ef = pd.read_csv(positions_ef_fname,header=None)
    pos_ef = np.array(pos_ef)[0]
    position_pft_names = [
        'a0_needleaf_trees',
        'a0_broadleaf_trees',
        'a0_shrub',
        'a0_c3_grass',
        'a0_c4_grass',
        'a0_crop'
    ]
    position_ef_names = [
        'b_needleleaf_trees',
        'b_broadleaf_trees',
        'b_shrub',
        'b_c3_grass',
        'b_c4_grass',
        'b_crop',
        'a0_intercept',
        'a1_precip_coef',
        'a2_canopy_coef'
    ]
    g1_pft_opt_need = pos_pft[0]
    g1_pft_opt_broad = pos_pft[1]
    g1_pft_opt_shrub = pos_pft[2]
    g1_pft_opt_c3 = pos_pft[3]
    g1_pft_opt_c4 = pos_pft[4]
    g1_pft_opt_crop = pos_pft[5]
    g1_ef_opt_need = pos_ef[0]*(
        pos_ef[6] +
        pos_ef[7]*opt_need_map +
        pos_ef[8]*opt_need_h
    )
    g1_ef_opt_broad = pos_ef[1]*(
        pos_ef[6] +
        pos_ef[7]*opt_broad_map +
        pos_ef[8]*opt_broad_h
    )
    g1_ef_opt_shrub = pos_ef[2]*(
        pos_ef[6] +
        pos_ef[7]*opt_shrub_map +
        pos_ef[8]*opt_shrub_h
    )
    g1_ef_opt_c3 = pos_ef[3]*(
        pos_ef[6] +
        pos_ef[7]*opt_c3_map +
        pos_ef[8]*opt_c3_h
    )
    g1_ef_opt_c4 = pos_ef[4]*(
        pos_ef[6] +
        pos_ef[7]*opt_c4_map +
        pos_ef[8]*opt_c4_h
    )
    g1_ef_opt_crop = pos_ef[5]*(
        pos_ef[6] +
        pos_ef[7]*opt_crop_map +
        pos_ef[8]*opt_crop_h
    )
    g1_pft_opts = {
        'needleleaf_trees':g1_pft_opt_need,
        'broadleaf_trees':g1_pft_opt_broad,
        'shrub':g1_pft_opt_shrub,
        'c3_grass':g1_pft_opt_c3,
        'c4_grass':g1_pft_opt_c4,
        'crop':g1_pft_opt_crop
    }
    g1_ef_opts = {
        'needleleaf_trees':g1_ef_opt_need,
        'broadleaf_trees':g1_ef_opt_broad,
        'shrub':g1_ef_opt_shrub,
        'c3_grass':g1_ef_opt_c3,
        'c4_grass':g1_ef_opt_c4,
        'crop':g1_ef_opt_crop
    }
    g1_ef_ext_need = pos_ef[0]*(
        pos_ef[6] +
        pos_ef[7]*ext_need_map +
        pos_ef[8]*ext_need_h
    )
    g1_ef_ext_broad = pos_ef[1]*(
        pos_ef[6] +
        pos_ef[7]*ext_broad_map +
        pos_ef[8]*ext_broad_h
    )
    g1_ef_ext_shrub = pos_ef[2]*(
        pos_ef[6] +
        pos_ef[7]*ext_shrub_map +
        pos_ef[8]*ext_shrub_h
    )
    g1_ef_ext_c3 = pos_ef[3]*(
        pos_ef[6] +
        pos_ef[7]*ext_c3_map +
        pos_ef[8]*ext_c3_h
    )
    g1_ef_ext_c4 = pos_ef[4]*(
        pos_ef[6] +
        pos_ef[7]*ext_c4_map +
        pos_ef[8]*ext_c4_h
    )
    g1_ef_ext_crop = pos_ef[5]*(
        pos_ef[6] +
        pos_ef[7]*ext_crop_map +
        pos_ef[8]*ext_crop_h
    )
    g1_ef_exts = {
        'needleleaf_trees':g1_ef_ext_need,
        'broadleaf_trees':g1_ef_ext_broad,
        'shrub':g1_ef_ext_shrub,
        'c3_grass':g1_ef_ext_c3,
        'c4_grass':g1_ef_ext_c4,
        'crop':g1_ef_ext_crop
    }
    tot_ext = len(ext_pixels)
    for p,pft in enumerate(pfts):
        this_default_g1 = default_g1s[pft]
        this_lin_g1 = lin_g1s[pft]
        this_ef_opt_g1 = g1_ef_opts[pft]
        this_ef_ext_g1 = g1_ef_exts[pft]
        this_pft_g1 = g1_pft_opts[pft]
        this_pix = len(this_ef_ext_g1)
        this_perc = len(this_ef_ext_g1)/tot_ext*100
        this_ef_g1_avg = np.mean(this_ef_ext_g1)
        this_ef_g1_median = np.median(this_ef_ext_g1)
        print('total pix for {} is {} at {}%'.format(pft,this_pix,this_perc))
        save_name = os.path.join(
            plots_dir,
            'ext_{}_g1_pdf.svg'
        )
        p_o.plot_pdf(
            [this_ef_opt_g1,this_ef_ext_g1],
            ['ef_opt_g1','ef_ext_g1'],
            ['#0048c8','#87c56d'],
            [0.5,0.5],
            [3,3],
            [False,True],
            save_name.format(pft),
            vert_lines=[
                this_ef_g1_median,this_pft_g1,
                this_default_g1,this_lin_g1
            ],
            vert_line_labels=[
                'ef_median','pft_avg',
                'default','lin'
            ],
            vert_line_colors=[
                '#358414','#8e5e96',
                '#080a0b','#ef3f29'
            ],
            xlim=[0,8],
            min_val=[0.5],
            fill_alphas=[1,0.5]
        )




















if __name__ == '__main__':
    main()
