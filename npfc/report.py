"""
Module report
================

A module with helper functions for computing pre-defined plots for the analysis
of fragment combinations.
"""

# standard
import warnings
import logging
import argparse
import sys
from shutil import rmtree
from datetime import datetime
import re
from pathlib import Path
from collections import OrderedDict
# data handling
import numpy as np
import json
import pandas as pd
from pandas import DataFrame
import networkx as nx
# data visualization
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker
# from pylab import savefig
from adjustText import adjust_text
# chemoinformatics
import rdkit
from rdkit import Chem
from rdkit.Chem import Mol
# docs
from typing import List
from typing import Tuple
from rdkit import RDLogger
# custom libraries
import npfc
from npfc import utils
from npfc import load
from npfc import save
from npfc import draw
from npfc import fragment_combination
from multiprocessing import Pool


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


DEFAULT_PALETTE = {'gray': '#808080',
                   'green': '#2CA02C',
                   'blue': '#378EBF',
                   'red': '#EB5763',
                   'orange': '#EBC81A',
                   }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def init_report_globals():
    """This function is used to init a reporting session by setting the default
    parameters for displaying information.
    """
    pd.set_option('display.max_columns', 20)
    pd.set_option('display.max_rows', 20)
    pd.set_option('max_colwidth', 70)


def print_title(title: str, level: int = 2, pad: int = 160):
    """Print a title with padding.

    :param title: title to print
    :param level: level of importance of the title (lower level is more important, min=1)
    :param pad: a padding added for highlighting the title ("="*pad)
    """
    if level == 1:
        # super title
        logging.info("=" * pad)
        logging.info(title.upper().center(pad))
        logging.info("=" * pad)
    elif level == 2:
        # normal title
        logging.info("*" * pad)
        logging.info(title.upper().center(pad))
        logging.info("*" * pad)
    elif level == 3:
        # subtitle
        logging.info("-" * pad)
        logging.info(title.upper().center(pad))
        logging.info("-" * pad)
    elif level == 4:
        # sub-subtitle
        logging.info("." * pad)
        logging.info(title.upper().center(pad))
        logging.info("." * pad)
    else:
        raise ValueError(f"ERROR! UNKNOWN LEVEL! ({level})")


def _get_chunks(WD: str, pattern: str) -> List[str]:
    """This function returns the list of all chunks found in a given directory.
    It does not use the _000_ notation because it can also be applied to fragments.

    :param WD: the WD where all subfolders for each step are present (i.e. 'natural/dnp/data')
    :param pattern: the pattern that desired chunks contain (i.e. '1_input/data/*([0-9][0-9][0-9])?.sdf.gz')
    :return: the list of chunks
    """
    WD = Path(WD)
    pattern = re.compile(pattern)
    chunks = [str(x) for x in list(WD.glob("*"))]
    return sorted(list(filter(pattern.match, chunks)))





def save_barplot(df: DataFrame,
                 output_plot: str,
                 x_name: str,
                 y_name: str,
                 title: str,
                 color: str,
                 x_label: str = None,
                 y_label: str = None,
                 rotate_x: int = 0,
                 perc_labels: str = None,
                 perc_label_size: int = 15,
                 fig_size: Tuple[int] = (24, 12),
                 force_order: bool = True,
                 print_rank: bool = False,
                 ):
    """This function helps for computing automated barplots using seaborn.
    It sets up somewhat standardized figure output for a harmonized rendering.

    :param df: the DataFrame with data to plot
    :param output_plot: the output plot full file name
    :param x_name: DF column name to use for x-axis
    :param y_name: DF column name to use for y-axis
    :param x_label: the name to display on the plot for x-axis
    :param y_label: the name to display on the plot for y-axis
    :param rotate_x: rotate the x-axis ticks anticlock-wise
    :param perc_labels: DF column name to use display percentage labels above bars
    :param perc_label_size: annotation text size for perc_labels
    :param color: color to use for bars, theoritically could also be a list of colors
    :param fig_size: tuple of integers defining the plot dimensions (x, y)
    :param force_order: force the plot to display the bars in the provided order. I do not know why sometimes this is needed, sometimes not depending on the plot.
    :param print_rank: display the rank as #i above the bar label
    :return: the figure in searborn format
    """
    # detect format from file extension
    format = Path(output_plot).suffix[1:].lower()
    if format != 'svg' and format != 'png':
        raise ValueError(f"ERROR! UNKNOWN PLOT FORMAT! ('{format}')")
    logging.debug(f"FORMAT FOR PLOT: '{format}'")

    # delete existing file for preventing stacking of plots
    p = Path(output_plot)
    if p.exists():
        p.unlink()

    # general style for plotting
    # sns.set(rc={'figure.figsize': fig_size})
    plt.figure(figsize=fig_size)
    sns.set_style('whitegrid', {'axes.edgecolor': '0.2'})
    sns.set_context("paper", font_scale=2)

    # barplot
    if force_order:
        ax = sns.barplot(x=df[x_name], y=df[y_name], color=color, order=df[x_name])
    else:
        ax = sns.barplot(x=df[x_name], y=df[y_name], color=color)
    ax.set_title(title, fontsize=24, y=1.02)
    ax.tick_params(labelsize=20)
    ax.tick_params(axis='x', rotation=rotate_x)
    # ax.set_xticklabels(df[x_name])
    ax.set_xlabel(x_label, fontsize=25, labelpad=20)
    ax.set_ylabel(y_label, fontsize=25, labelpad=20)
    # if no entry, then y ticks get all confused
    if df[y_name].sum() == 0:
        ax.set_yticks((0, 1, 2, 3, 4, 5))
        ax.set_ylim((0, 5))

    # format y labels
    label_format = '{:,.0f}'
    ticks_loc = ax.get_yticks().tolist()
    ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    ax.set_yticklabels([label_format.format(x) for x in ticks_loc])

    # add percentage annotations
    if perc_labels is not None:
        for a, i in zip(ax.patches, range(len(df.index))):
            row = df.iloc[i]
            ax.text(row.name, a.get_height(), row[perc_labels], color='black', ha='center', fontdict={'fontsize': perc_label_size})

    # save
    figure = ax.get_figure()
    figure.subplots_adjust(bottom=0.2)
    figure.savefig(output_plot, format=format, dpi=600)
    plt.clf()
    plt.close()

    return figure


def save_lineplot(df: DataFrame,
                  output_plot: str,
                  x_name: str,
                  y_name: str,
                  title: str,
                  color: str,
                  x_label: str = None,
                  y_label: str = None,
                  rotate_x: int = 0,
                  perc_labels: str = None,
                  perc_label_size: int = 15,
                  fig_size: Tuple[int] = (24, 12),
                  fill: bool = False,
                  ):
    """This function helps for computing automated lineplot using seaborn.
    It sets up somewhat standardized figure output for a harmonized rendering.

    :param df: the DataFrame with data to plot
    :param output_plot: the output plot full file name
    :param x_name: DF column name to use for x-axis
    :param y_name: DF column name to use for y-axis
    :param x_label: the name to display on the plot for x-axis
    :param y_label: the name to display on the plot for y-axis
    :param rotate_x: rotate the x-axis ticks anticlock-wise
    :param perc_labels: DF column name to use display percentage labels above bars
    :param perc_label_size: annotation text size for perc_labels
    :param color: color to use for bars, theoritically could also be a list of colors
    :param fig_size: tuple of integers defining the plot dimensions (x, y)
    :return: the figure in searborn format
    """
    # detect format from file extension
    format = Path(output_plot).suffix[1:].lower()
    if format != 'svg' and format != 'png':
        raise ValueError(f"ERROR! UNKNOWN PLOT FORMAT! ('{format}')")
    logging.debug(f"FORMAT FOR PLOT: '{format}'")

    # delete existing file for preventing stacking of plots
    p = Path(output_plot)
    if p.exists():
        p.unlink()

    # general style for plotting
    # sns.set(rc={'figure.figsize': fig_size})
    plt.figure(figsize=fig_size)
    sns.set_style('whitegrid', {'axes.edgecolor': '0.2'})
    sns.set_context("paper", font_scale=2)

    # lineplot
    fig, ax = plt.subplots(figsize=fig_size)
    ax = sns.scatterplot(x=x_name, y=y_name, data=df, ax=ax, s=20)
    ax = sns.lineplot(x=df[x_name], y=df[y_name], color=color, ax=ax)
    # fix xticks that would contain floats instead of ints
    plt.xticks(np.arange(min(df[x_name]), max(df[x_name]) + 1, 1))  # https://stackoverflow.com/questions/30327153/seaborn-pylab-changing-xticks-from-float-to-int
    # update all lines width
    for l in ax.lines:
        plt.setp(l, linewidth=3)
    ax.margins(0, 0)
    ax.set_title(title, fontsize=24, y=1.02)
    ax.tick_params(labelsize=20)
    ax.tick_params(axis='x', rotation=rotate_x)
    # ax.set_xticklabels(df[x_name])
    ax.set_xlabel(x_label, fontsize=25, labelpad=20)
    ax.set_ylabel(y_label, fontsize=25, labelpad=20)
    # if no entry, then y ticks get all confused
    if df[y_name].sum() == 0:
        ax.set_yticks((0, 1, 2, 3, 4, 5))
        ax.set_ylim((0, 5))

    # format y labels
    label_format = '{:,.0f}'
    ticks_loc = ax.get_yticks().tolist()
    ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    ax.set_yticklabels([label_format.format(x) for x in ticks_loc])

    # add percentage annotations
    if perc_labels is not None:
        for i in range(len(df.index)):
            row = df.iloc[i]
            ax.text(row[x_name], row[y_name], row[perc_labels], color='black', ha='center', fontdict={'fontsize': perc_label_size})

    if fill:
        plt.fill_between(df[x_name].values, df[y_name].values)

    # adjust overlapping text annotations (percentages with axes, etc.)
    texts = [x for x in ax.texts]
    adjust_text(texts)

    # save
    figure = ax.get_figure()
    figure.subplots_adjust(bottom=0.2)
    figure.savefig(output_plot, format=format, dpi=600)
    plt.clf()
    plt.close()

    return figure


def save_kdeplot(df: DataFrame,
                 output_plot: str,
                 x_name: str,
                 title: str,
                 color: str,
                 x_label: str = None,
                 y_label: str = None,
                 normalize_x: bool = True,
                 fig_size: Tuple[int] = (24, 12),
                 ):
    """This function helps for computing automated kdeplots using seaborn.
    It sets up somewhat standardized figure output for a harmonized rendering.

    :param df: the DataFrame with data to plot
    :param output_plot: the output plot full file name
    :param x_name: DF column name to use for x-axis
    :param x_label: the name to display on the plot for x-axis
    :param y_label: the name to display on the plot for y-axis
    :param color: color to use for bars, theoritically could also be a list of colors
    :param fig_size: tuple of integers defining the plot dimensions (x, y)
    :return: the figure in searborn format
    """
    # detect format from file extension
    format = Path(output_plot).suffix[1:].lower()
    if format != 'svg' and format != 'png':
        raise ValueError(f"ERROR! UNKNOWN PLOT FORMAT! ('{format}')")
    logging.debug(f"FORMAT FOR PLOT: '{format}'")

    # delete existing file for preventing stacking of plots
    p = Path(output_plot)
    if p.exists():
        p.unlink()

    # general style for plotting
    sns.set(rc={'figure.figsize': fig_size})
    sns.set_style('whitegrid', {'axes.edgecolor': '0.2'})
    sns.set_context("paper", font_scale=2)

    ax = sns.kdeplot(df[x_name], shade=True, label='', color=color)
    ax.set_title(title, fontsize=24, y=1.02)
    ax.tick_params(labelsize=20)
    ax.tick_params(axis='x', rotation=0)
    ax.set_xlim(0, 1)
    # ax.set_xticklabels(df[x_name])
    label_format = '{:,.0%}'
    ticks_loc = ax.get_xticks().tolist()
    ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    ax.set_xticklabels([label_format.format(x) for x in ticks_loc])

    #ax.set_xticklabels(['{:,.0%}'.format(x) for x in ax.get_xticks()])
    ax.set_xlabel(x_label, fontsize=25, labelpad=20)
    ax.set_ylabel(y_label, fontsize=25, labelpad=20)

    # save
    figure = ax.get_figure()
    figure.savefig(output_plot, dpi=600)
    plt.clf()
    plt.close()

    return figure


def save_distplot(df: DataFrame,
                  output_plot: str,
                  x_name: str,
                  title: str,
                  color: str,
                  x_label: str = None,
                  y_label: str = None,
                  normalize_x: bool = True,
                  fig_size: Tuple[int] = (24, 12),
                  bins=None
                  ):
    """This function helps for computing automated distplots using seaborn.
    It sets up somewhat standardized figure output for a harmonized rendering.

    :param df: the DataFrame with data to plot
    :param output_plot: the output plot full file name
    :param x_name: DF column name to use for x-axis
    :param x_label: the name to display on the plot for x-axis
    :param y_label: the name to display on the plot for y-axis
    :param color: color to use for bars, theoritically could also be a list of colors
    :param fig_size: tuple of integers defining the plot dimensions (x, y)
    :param bins: number of bins
    :return: the figure in searborn format
    """
    # detect format from file extension
    format = Path(output_plot).suffix[1:].lower()
    if format != 'svg' and format != 'png':
        raise ValueError(f"ERROR! UNKNOWN PLOT FORMAT! ('{format}')")
    logging.debug(f"FORMAT FOR PLOT: '{format}'")

    # delete existing file for preventing stacking of plots
    p = Path(output_plot)
    if p.exists():
        p.unlink()

    # general style for plotting
    sns.set(rc={'figure.figsize': fig_size})
    sns.set_style('whitegrid', {'axes.edgecolor': '0.2'})
    sns.set_context("paper", font_scale=2)
    vals = []
    for i in range(len(df)):
        row = df.iloc[i]
        vals += [row[x_name]] * row['Count']
    ax = sns.distplot(vals, kde=False, label='', color=color, bins=bins, hist_kws=dict(alpha=1))
    # ax.set_xlim(0, df[x_name].max())
    ax.set_title(title, fontsize=24, y=1.02)
    ax.tick_params(labelsize=20)
    ax.tick_params(axis='x', rotation=0)

    label_format = '{:,}'
    ticks_loc = ax.get_yticks().tolist()
    ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    ax.set_yticklabels([label_format.format(x) for x in ticks_loc])

    # ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks()])
    ax.set_xlabel(x_label, fontsize=25, labelpad=20)
    ax.set_ylabel(y_label, fontsize=25, labelpad=20)

    # save
    figure = ax.get_figure()
    figure.savefig(output_plot, dpi=600)
    plt.clf()
    plt.close()

    return figure

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
