#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
main_reviewer_responses.py

Goal:
Compare 3 model configs for:
  a) soil moisture skill (vs ESA CCI ACTIVE)
     using CDF-matched soil moisture
  b) carbon cycle realism: how well does
     GPP track SIF?
     (pixelwise anomaly correlation)

This is structured to mirror the style of
main_final_paper_figures.py :contentReference[oaicite:0]{index=0}
so it drops summary printouts, hist-style
plots, and saves outputs for the paper
/ reviewer response bundle.
"""

import sys
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/' +
    'pso/step_5.1_analyze_outputs/funcs'
)
sys.path.append(
    '/discover/nobackup/trobinet/from_aws/' +
    'pso/step_5.1_analyze_outputs/tasks'
)

import os
import datetime
import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt

# we will reuse plotting helpers so figures
# look like the rest of the paper
from plot_other import plot_other

########################################
# helper funcs
########################################

def load_truth_csv(fname,
                   date_col='date',
                   pixel_col='pixel',
                   val_col='val'):
    """
    Generic loader for truth csvs expected to
    be long-form:
        date, pixel, val
    Returns:
        df with columns:
        [date, pixel, val]
        date is datetime.date
    """

    df = pd.read_csv(fname)

    # try to coerce date
    if df[date_col].dtype != 'datetime64[ns]':
        df[date_col] = pd.to_datetime(
            df[date_col]
        ).dt.date

    # standardize column names
    df = df.rename(
        columns={
            date_col: 'date',
            pixel_col: 'pixel',
            val_col: 'val'
        }
    )

    # sort for sanity
    df = df.sort_values(
        by=['pixel', 'date']
    ).reset_index(drop=True)

    return df


def load_model_csv(fname,
                   date_col='date',
                   pixel_col='pixel',
                   val_col='val'):
    """
    Same idea as load_truth_csv, but for a
    model output file. Assumes long-form
    (date,pixel,val).

    You will call this for:
      - model soil moisture (sm)
      - model GPP (gpp)

    Returns standardized df with:
    [date, pixel, val]
    """
    return load_truth_csv(
        fname,
        date_col=date_col,
        pixel_col=pixel_col,
        val_col=val_col
    )


def cdf_match_model_to_truth(
        model_series,
        truth_series
    ):
    """
    Do empirical CDF matching pixelwise.

    Inputs:
        model_series : 1D np.array of model sm
        truth_series : 1D np.array of truth sm
    Both should already be aligned in time
    (same days) and have no NaNs.

    Returns:
        matched_model : 1D np.array
        with model distribution mapped
        onto truth distribution.

    Method:
    1. sort unique model vals
    2. get CDF rank for each unique val
    3. get the quantile of the truth at
       those same ranks
    4. interpolate back to each model val
    """

    # unique sorted model values
    m_sorted_unique = np.sort(
        np.unique(model_series)
    )

    # percentile rank of those values
    # rank in [0,1]
    m_ranks = np.argsort(
        np.argsort(m_sorted_unique)
    ).astype(float)
    if len(m_ranks) > 1:
        m_ranks = (
            m_ranks /
            (len(m_sorted_unique) - 1)
        )
    else:
        # edge case single value
        m_ranks = np.array([0.5])

    # now map ranks -> quantiles of truth
    # using np.quantile
    t_quants = np.quantile(
        truth_series,
        m_ranks
    )

    # build interp from model_value
    # -> matched_value
    matched_vals = np.interp(
        model_series,
        m_sorted_unique,
        t_quants
    )

    return matched_vals


def pixelwise_join(a_df,
                   b_df,
                   pixel,
                   start_date=None,
                   end_date=None):
    """
    Grab one pixel from a_df and b_df,
    inner join on date,
    optionally subset by date window.

    Expects columns: date, pixel, val
    Returns:
        adf, bdf merged on date
        (as plain np arrays)
    """
    a_sub = a_df[a_df['pixel'] == pixel]
    b_sub = b_df[b_df['pixel'] == pixel]

    if start_date is not None:
        a_sub = a_sub[
            a_sub['date'] >= start_date
        ]
        b_sub = b_sub[
            b_sub['date'] >= start_date
        ]
    if end_date is not None:
        a_sub = a_sub[
            a_sub['date'] <= end_date
        ]
        b_sub = b_sub[
            b_sub['date'] <= end_date
        ]

    merged = pd.merge(
        a_sub[['date','val']],
        b_sub[['date','val']],
        on='date',
        how='inner',
        suffixes=('_a','_b')
    ).dropna()

    return merged['val_a'].to_numpy(), \
           merged['val_b'].to_numpy()


def mae(a, b):
    """mean abs error, ignoring NaN"""
    diff = np.abs(a - b)
    return np.nanmean(diff)


def corr(a, b):
    """
    Pearson r, ignoring NaN.
    If <2 finite samples, return nan.
    """
    good = np.isfinite(a) & np.isfinite(b)
    if np.sum(good) < 2:
        return np.nan
    return np.corrcoef(
        a[good],
        b[good]
    )[0,1]


########################################
# analysis funcs
########################################

def evaluate_soil_moisture(
    truth_sm_df,
    model_sm_dfs,
    pixels=None,
    start_date=None,
    end_date=None,
    min_len=30
):
    """
    For each pixel and each model:
    1. intersect time series with truth
    2. CDF-match model->truth
    3. compute MAE between
       CDF-matched model and truth

    Inputs:
        truth_sm_df : df with cols
                      [date,pixel,val]
        model_sm_dfs: dict of
                      {model_name: df}
                      each df same cols
        pixels      : list of pixels
                      to analyze.
                      If None, use
                      union of truth.
        start_date, end_date:
                      datetime.date or None
        min_len     : min # of days to call
                      this pixel "valid"

    Returns:
        sm_skill_df : DataFrame
            rows   = pixel
            cols   = model_name
            val    = MAE after CDF-match
        summary    : dict with model means
    """

    if pixels is None:
        pixels = sorted(
            list(
                truth_sm_df['pixel'].unique()
            )
        )

    # per-pixel MAE per model
    out_df = pd.DataFrame(
        index=pixels,
        columns=list(model_sm_dfs.keys()),
        dtype=float
    )

    for pix in pixels:
        # truth series for this pixel
        # (we only pull once below)
        for mname, mdf in model_sm_dfs.items():
            # join model + truth
            mod_vals, tru_vals = pixelwise_join(
                mdf,
                truth_sm_df,
                pix,
                start_date,
                end_date
            )

            # skip short joins
            if len(mod_vals) < min_len:
                out_df.loc[pix, mname] = np.nan
                continue

            # CDF match model->truth
            matched_mod = cdf_match_model_to_truth(
                mod_vals,
                tru_vals
            )

            # score
            this_mae = mae(
                matched_mod,
                tru_vals
            )
            out_df.loc[pix, mname] = this_mae

    # summarize mean MAE across pixels
    summary = {}
    for mname in model_sm_dfs.keys():
        summary[mname] = np.nanmean(
            out_df[mname].to_numpy()
        )

    return out_df, summary


def evaluate_sif_gpp_corr(
    truth_sif_df,
    model_gpp_dfs,
    pixels=None,
    start_date=None,
    end_date=None,
    min_len=30
):
    """
    For each pixel and each model:
    1. join daily SIF and daily GPP
    2. subtract per-pixel mean from
       BOTH SIF and GPP to get
       anomalies
    3. compute Pearson r between
       SIF_anom and GPP_anom

    Inputs:
        truth_sif_df : df cols
                       [date,pixel,val]
        model_gpp_dfs: dict of
                       {model_name: df}
        pixels       : list of pixels
        start_date, end_date : same
        min_len      : min overlap days

    Returns:
        corr_df : DataFrame
            rows   = pixel
            cols   = model_name
            val    = anomaly corr
        summary : dict mean corr
    """

    if pixels is None:
        pixels = sorted(
            list(
                truth_sif_df['pixel'].unique()
            )
        )

    out_df = pd.DataFrame(
        index=pixels,
        columns=list(model_gpp_dfs.keys()),
        dtype=float
    )

    for pix in pixels:
        for mname, gpp_df in model_gpp_dfs.items():

            gpp_vals, sif_vals = pixelwise_join(
                gpp_df,
                truth_sif_df,
                pix,
                start_date,
                end_date
            )

            if len(gpp_vals) < min_len:
                out_df.loc[pix, mname] = np.nan
                continue

            # anomaly = value - pixel mean
            gpp_anom = (
                gpp_vals - np.nanmean(gpp_vals)
            )
            sif_anom = (
                sif_vals - np.nanmean(sif_vals)
            )

            out_df.loc[pix, mname] = corr(
                gpp_anom,
                sif_anom
            )

    summary = {}
    for mname in model_gpp_dfs.keys():
        summary[mname] = np.nanmean(
            out_df[mname].to_numpy()
        )

    return out_df, summary


def make_pdf_plot(
    arrays_list,
    labels_list,
    colors_list,
    save_name,
    xlim=None
):
    """
    Thin wrapper around plot_other.plot_pdf
    to match fig style used elsewhere.

    arrays_list : list of 1D arrays
    labels_list : list of str
    colors_list : list of color str
    save_name   : str, path to save
    xlim        : 2-elem list or None
    """

    # default styling borrowed from
    # main_final_paper_figures.py
    kernels = [1 for _ in arrays_list]
    lw      = [1 for _ in arrays_list]
    fills   = [False for _ in arrays_list]
    styles  = ['-' for _ in arrays_list]

    p_o = plot_other()
    p_o.plot_pdf(
        arrays_list,
        labels_list,
        colors_list,
        kernels,
        lw,
        fills,
        styles,
        save_name,
        xlim=xlim
    )


########################################
# main
########################################

def main():

    ####################################
    # paths / config
    ####################################

    base_dir = (
        '/discover/nobackup/trobinet/' +
        'from_aws/pso/step_7_reviews'
    )

    plots_dir = os.path.join(
        base_dir,
        'plots'
    )
    out_dir = os.path.join(
        base_dir,
        'outputs'
    )

    os.makedirs(plots_dir,
                exist_ok=True)
    os.makedirs(out_dir,
                exist_ok=True)

    # reviewer windows:
    # soil moisture truth is 1991-10-01
    # thru 2014-12-31
    # SIF truth is 2001-01-01
    # thru 2014-12-31
    sm_start = datetime.date(
        1991,10,1
    )
    sm_end   = datetime.date(
        2014,12,31
    )
    sif_start = datetime.date(
        2001,1,1
    )
    sif_end   = datetime.date(
        2014,12,31
    )

    # truth CSVs you gave
    sm_truth_fname = os.path.join(
        base_dir,
        'outputs',
        'sm_truth_cci_active_pct_' +
        '1991-10-01_2014-12-31.csv'
    )

    sif_truth_fname = os.path.join(
        base_dir,
        'outputs',
        'sif_truth_oco2_all_daily_' +
        '2001-01-01_2014-12-31.csv'
    )

    ####################################
    # model info
    ####################################

    # We assume 3 model configs that we
    # want to compare, analogous to
    # "default", "pft", "ef" in the main
    # paper workflow.
    #
    # TODO:
    # fill in actual file paths.
    #
    # Assumed long-form columns:
    #   date,pixel,val
    #
    # Soil moisture model outputs
    # after running land model.
    # GPP model outputs from same
    # three experiments.

    model_info = {
        'default':{
            'soil_moisture_csv': (
                'TODO/path/default_sm.csv'
            ),
            'gpp_csv':(
                'TODO/path/default_gpp.csv'
            )
        },
        'pft':{
            'soil_moisture_csv': (
                'TODO/path/pft_sm.csv'
            ),
            'gpp_csv':(
                'TODO/path/pft_gpp.csv'
            )
        },
        'ef':{
            'soil_moisture_csv': (
                'TODO/path/ef_sm.csv'
            ),
            'gpp_csv':(
                'TODO/path/ef_gpp.csv'
            )
        }
    }

    ####################################
    # load truth
    ####################################

    sm_truth_df = load_truth_csv(
        sm_truth_fname,
        date_col='date',
        pixel_col='pixel',
        val_col='val'  # % soil moisture
    )

    sif_truth_df = load_truth_csv(
        sif_truth_fname,
        date_col='date',
        pixel_col='pixel',
        val_col='val'  # SIF
    )

    ####################################
    # load model data
    ####################################

    model_sm_dfs = {}
    model_gpp_dfs = {}

    for mname,info in model_info.items():

        this_sm = load_model_csv(
            info['soil_moisture_csv'],
            date_col='date',
            pixel_col='pixel',
            val_col='val'   # model SM
        )

        this_gpp = load_model_csv(
            info['gpp_csv'],
            date_col='date',
            pixel_col='pixel',
            val_col='val'   # model GPP
        )

        # store
        model_sm_dfs[mname] = this_sm
        model_gpp_dfs[mname] = this_gpp

    ####################################
    # figure A:
    # soil moisture skill
    ####################################

    sm_skill_df, sm_summary = (
        evaluate_soil_moisture(
            truth_sm_df=sm_truth_df,
            model_sm_dfs=model_sm_dfs,
            pixels=None,
            start_date=sm_start,
            end_date=sm_end,
            min_len=30
        )
    )

    # print summary MAE per model
    print('Soil moisture MAE (CDF-' +
          'matched) per model:')
    for mname,val in sm_summary.items():
        print('  {}: {}'.format(
            mname,val
        ))

    # relative improvement numbers,
    # like we did with J improvements.
    # We'll treat "ef" as the fancy
    # model and compare vs pft/default
    if 'ef' in sm_summary and \
       'pft' in sm_summary:
        diff_ef_pft = (
            sm_summary['ef'] -
            sm_summary['pft']
        )
        perc_imp_ef_pft = (
            diff_ef_pft /
            sm_summary['pft']
        ) * -100.0
        print(
            'EF improved SM MAE ' +
            'over PFT by {}%'.format(
                perc_imp_ef_pft
            )
        )

    if 'ef' in sm_summary and \
       'default' in sm_summary:
        diff_ef_def = (
            sm_summary['ef'] -
            sm_summary['default']
        )
        perc_imp_ef_def = (
            diff_ef_def /
            sm_summary['default']
        ) * -100.0
        print(
            'EF improved SM MAE ' +
            'over default by {}%'.format(
                perc_imp_ef_def
            )
        )

    # plot PDF of per-pixel MAE for
    # each model to show distribution
    sm_mae_arrays = []
    sm_labels = []
    sm_colors = [
        '#000000',
        '#dc70ab',
        '#86c148'
    ]
    # make sure order matches color
    desired_order = ['default',
                     'pft',
                     'ef']
    for m_i,mname in enumerate(
        desired_order
    ):
        if mname in sm_skill_df.columns:
            sm_mae_arrays.append(
                sm_skill_df[mname].to_numpy()
            )
            sm_labels.append(mname)

    sm_pdf_name = os.path.join(
        plots_dir,
        'soil_moisture_mae_pdf.svg'
    )

    make_pdf_plot(
        arrays_list=sm_mae_arrays,
        labels_list=sm_labels,
        colors_list=sm_colors[
            0:len(sm_labels)
        ],
        save_name=sm_pdf_name,
        xlim=[0,0.3]  # TODO tweak
    )

    ####################################
    # figure B:
    # SIF-vs-GPP correlation
    ####################################

    sif_corr_df, sif_summary = (
        evaluate_sif_gpp_corr(
            truth_sif_df=sif_truth_df,
            model_gpp_dfs=model_gpp_dfs,
            pixels=None,
            start_date=sif_start,
            end_date=sif_end,
            min_len=30
        )
    )

    print('Mean pixelwise anomaly ' +
          'corr(SIF,GPP) per model:')
    for mname,val in sif_summary.items():
        print('  {}: {}'.format(
            mname,val
        ))

    if 'ef' in sif_summary and \
       'pft' in sif_summary:
        diff = (
            sif_summary['ef'] -
            sif_summary['pft']
        )
        print(
            'EF minus PFT mean ' +
            'corr: {}'.format(diff)
        )

    if 'ef' in sif_summary and \
       'default' in sif_summary:
        diff = (
            sif_summary['ef'] -
            sif_summary['default']
        )
        print(
            'EF minus default mean ' +
            'corr: {}'.format(diff)
        )

    # PDF / histogram-style plot
    # of correlation values
    corr_arrays = []
    corr_labels = []
    corr_colors = [
        '#000000',
        '#dc70ab',
        '#86c148'
    ]
    for m_i,mname in enumerate(
        desired_order
    ):
        if mname in sif_corr_df.columns:
            corr_arrays.append(
                sif_corr_df[mname].to_numpy()
            )
            corr_labels.append(mname)

    corr_pdf_name = os.path.join(
        plots_dir,
        'sif_gpp_corr_pdf.svg'
    )

    make_pdf_plot(
        arrays_list=corr_arrays,
        labels_list=corr_labels,
        colors_list=corr_colors[
            0:len(corr_labels)
        ],
        save_name=corr_pdf_name,
        xlim=[-1,1]
    )

    ####################################
    # save numeric outputs we may want
    # to cite directly in the response
    ####################################

    sm_skill_out = os.path.join(
        out_dir,
        'soil_moisture_mae_by_pixel.csv'
    )
    sm_skill_df.to_csv(
        sm_skill_out
    )

    sif_corr_out = os.path.join(
        out_dir,
        'sif_gpp_corr_by_pixel.csv'
    )
    sif_corr_df.to_csv(
        sif_corr_out
    )

    sm_summary_out = os.path.join(
        out_dir,
        'soil_moisture_mae_summary.csv'
    )
    pd.Series(sm_summary).to_csv(
        sm_summary_out
    )

    sif_summary_out = os.path.join(
        out_dir,
        'sif_gpp_corr_summary.csv'
    )
    pd.Series(sif_summary).to_csv(
        sif_summary_out
    )

    print('Wrote:')
    print('  ' + sm_skill_out)
    print('  ' + sif_corr_out)
    print('  ' + sm_summary_out)
    print('  ' + sif_summary_out)


if __name__ == '__main__':
    main()

