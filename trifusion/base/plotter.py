#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  Copyright 2012 Unknown <diogo@arch>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

import matplotlib
# Using the agg backend prevent a segmentation fault when quiting the
# application
matplotlib.use("agg")

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.stats import spearmanr
from scipy.interpolate import spline
import numpy as np
from itertools import chain
from collections import Counter
import seaborn as sns
import os

# Set of 10 easily distinguishable colors that will be used when generating plots
# to ensure consistency in color usage. In case more than 10 colors are required
# (ugh...) they will be randomly generated henceforth.

clr_list = [[0, .53, .66],  # light blue
            [1, .16, .16],     # light red
            [.21, .78, .22],  # light green
            [.83, .55, .37],  # light brown
            [.55, .37, .83],  # purple
            [1, .4, 0],     # light orange
            [.77, .77, .77],  # light grey
            [.78, .22, .22],  # dark red
            [.1, .1, .1],     # dark grey
            [0, .66, 0]]      # green

two_clr_list = [[0, .53, .66],
                [.62, .80, .85]]


class SetProps(object):
    """
    This decorator is use to automatically set several plot properties based
    on the arguments provided to the plotting functions. In this way,
    there is no need to write repetitive code in each function. The supported
    arguments are:

    .: ax_names: list, first element x-axis label; second element y-axis label
    .: title: string, title of the plot
    """

    def __init__(self, func):
        self.func = func

    def __call__(self, *args, **kwargs):

        plt.clf()
        plt.close()

        res = self.func(*args, **kwargs)

        # Set axis names
        override = {
            "fontsize": 14,
            "fontweight": "bold",
            "color": "grey"
        }
        if "ax_names" in kwargs:
            if kwargs["ax_names"][0]:
                plt.xlabel(kwargs["ax_names"][0], labelpad=15, **override)
            if kwargs["ax_names"][1]:
                plt.ylabel(kwargs["ax_names"][1], labelpad=15, **override)

        # Set title
        if "title" in kwargs:
            plt.title(kwargs["title"], y=1.05, **override)

        return res


def _add_labels(ax_names):
    """
    Wrapper to the add labels routine
    :param ax_names: List of two items. [x-axis name, y-axis name]
    """

    if ax_names[0]:
        plt.xlabel(ax_names[0])

    if ax_names[1]:
        plt.ylabel(ax_names[1])


@SetProps
def scatter_plot(data, correlation=False, ax_names=None, table_header=None,
                 title=None):
    """
    Builds a scatter plot from a 2D data list. Also calculates the
    correlation coefficient is requested.
    :param data: 2D list, containing the x and y data points
    :param correlation: Boolean. If True, the spearman's rank correlation
    coefficient is also calculated
    :param ax_names: list, with two string elements for the x and y axis
    :param table_header: list, with table header. Each element represents a
    column
    """

    plt.rcParams["figure.figsize"] = (8, 6)

    # Use ggpot style
    plt.style.use("ggplot")

    fig, ax = plt.subplots()

    # Generate scatter plot
    ax.scatter(data[0], data[1], alpha=0.5, edgecolors=None)

    # Plot best fit line
    ax.plot(data[0], np.poly1d(np.polyfit(data[0], data[1], 1))(data[0]))

    # Set plot limits
    max1 = max(data[0])
    max2 = max(data[1])
    ax.set_xlim(0, max1 + (max1 * .1))
    ax.set_ylim(0, max2 + (max2 * .1))

    # Calculate spearman's rank correlation coefficient and add it to the plot
    if correlation:
        rho, pval = spearmanr(data[0], data[1])

        ax.text((ax.get_xlim()[1] - ax.get_xlim()[0]) * .05,
                ax.get_ylim()[1] - (ax.get_ylim()[1] * .05),
                r"$\rho$ = {}; p-value = {}".format(round(rho, 4),
                                                    round(pval, 4)))

    # Generate table data
    if table_header:
        table = [table_header]
    else:
        table = []

    for i, j in zip(*data):
        table.append([i, j])

    return fig, None, table


@SetProps
def bar_plot(data, labels=None, title=None, ax_names=None,
             lgd_list=None, reverse_x=False, table_header=None):
    """
    Builds simple bar plot from a data_list
    :param data: list with data to be plotted.
    :param labels: list with x axis labels
    :param title: string, plot title
    :param ax_names: list. Names of the axis [yaxis_name, xaxis_name]
    :param reverse_x: Boolean, determines whether the x-axis is reversed or not
    """

    if len(labels) > 10:
        plt.rcParams["figure.figsize"] = (len(labels) / 3, 6)
    else:
        plt.rcParams["figure.figsize"] = (8, 6)

    # Use ggpot style
    plt.style.use("ggplot")

    fig, ax = plt.subplots()

    # Create plot
    for i, d in enumerate(data):
        # Get color from 10 color list. If more than 10 colors are required,
        # randomly generate new ones
        if len(data) > 2:
            try:
                clr = clr_list[i]
            except IndexError:
                clr = np.random.rand(3, 1)
        else:
            clr = two_clr_list[i]

        # Create plot for first entry
        if d == data[0]:
            ax.bar(np.arange(len(d)), d, align="center", color=clr,
                   label=lgd_list[i] if lgd_list else None)
        # Add plots on top
        else:
            ax.bar(np.arange(len(d)), d, align="center", color=clr,
                   bottom=data[i - 1],
                   label=lgd_list[i] if lgd_list else None)

    # Set legend
    if lgd_list:
        lgd = plt.legend(bbox_to_anchor=(1, .5), loc="center left",
                         fancybox=True, shadow=True, framealpha=0.8)
    else:
        lgd = None

    # Set labels at the center of bars
    if labels:
        ax.set_xticks(np.arange(len(labels)))
        # Determine rotation and alignmnet from labels
        if max([len(x) for x in labels]) >= 3:
            ax.set_xticklabels(labels, ha="right", rotation=45)
        else:
            ax.set_xticklabels(labels, ha="center")

        # Setting the x range to avoid large spaces between the axis and bars
        plt.xlim([min(np.arange(len(labels))) - .5,
                  max(np.arange(len(labels))) + .5])

    # Invert x-axis so that higher taxa prevalence is shown in the left
    if reverse_x:
        plt.gca().invert_xaxis()

    # Generating table structure
    if table_header:
        table = [table_header]
    else:
        table = []

    for l, x in zip(labels, data[0]):
        table.append([l, x])

    return fig, lgd, table


@SetProps
def multi_bar_plot(data_list, labels=None, lgd_list=None, title=None,
                   ax_names=None):
    """
    General purpose multiple bar plot with custom layout. Returns a pyplot
    object.
    :param data_list: list with data to be plotted in list/array type in each
    entry. Data with two groups would be like [[1.23, .53], [1.55, .12]]
    :param labels: list, containing the labels
    :param lgd_list: list, The legend string for each data set in data_list
    """

    plt.rcParams["figure.figsize"] = (8, 6)

    # Use ggpot style
    plt.style.use("ggplot")

    fig, ax = plt.subplots()

    # Set margin between bar groups
    margin = 0.1

    # Determine bar group width according to the number of data lists
    w = (1. - 2. * margin) / len(data_list)
    # Determine max plt height to calculate text label space
    max_height = max((x for y in data_list for x in y))

    # Create bar plots
    for i, d in enumerate(data_list):
        # Get color from 10 color list. If more than 10 colors are required,
        # randomly generate new ones
        try:
            clr = clr_list[i]
        except IndexError:
            clr = np.random.rand(3, 1)

            # Determine position of xdata
        xdata = np.arange(len(d)) + margin + (i * w)
        bplt = ax.bar(xdata, d, w, color=clr, label=lgd_list[i])

    # Add labels
    if labels:
        ax.set_xticks(np.arange(len(labels)) + 0.5)
        ax.set_xticklabels(labels, rotation=45, ha="center", fontweight="bold")

    # Set label colors
    ax.tick_params(axis="x", colors=[.4, .4, .4])
    ax.tick_params(axis="y", colors=[.3, .3, .3])

    # Add legend
    lgd = plt.legend(bbox_to_anchor=(1, .5), loc="center left",
                     fancybox=True, shadow=True, framealpha=0.8)

    # Add text to the top of bar with the corresponding value
    for b in bplt:

        # If the bar width is large enough to accommodate the text with
        # horizontal orientation
        if w >= 0.8:
            xpos = b.get_x() + b.get_width() / 2
            ypos = b.get_height() * .90
            ax.text(xpos, ypos, b.get_height(), horizontalalignment="center",
                    weight='bold', color="white", size=14)
        else:
            if b.get_height() / max_height >= 0.2:
                xpos = b.get_x() + b.get_width() / 2 * 1.3
                ypos = b.get_height() * .95
                ax.text(xpos, ypos, b.get_height(),
                        horizontalalignment="center", color="white",
                        verticalalignment="top", weight='bold',
                        size=14, rotation=90)
            else:
                xpos = b.get_x() + b.get_width() / 2
                ypos = b.get_height() * 1.1
                ax.text(xpos, ypos, b.get_height(),
                        horizontalalignment="center", color="black",
                        verticalalignment="bottom", weight='bold',
                        size=14, rotation=90)

    # Automatically adjust figure size to accommodate all labels
    # plt.tight_layout()

    return fig, lgd


@SetProps
def interpolation_plot(data, title=None, ax_names=None):
    """
    Creates a black and white interpolation plot from data, which must consist
    of a 0/1 matrix for absence/presence of taxa in genes
    :param data: numpy array of variable shape.
    """

    plt.rcParams["figure.figsize"] = (8, 6)

    # Use ggpot style
    plt.style.use("ggplot")

    fig, ax = plt.subplots()

    # Setting the aspect ratio proportional to the data
    ar = float(len(data[0])) / float(len(data)) * .2

    ax.imshow(data, interpolation="none", cmap="Greys", aspect=ar)

    ax.grid(False)

    return fig, None, None


@SetProps
def stacked_bar_plot(data, labels, legend=None, table_header=None, title=None,
                     ax_names=None, normalize=False, normalize_factor=None):
    """
    Creates a tight stacked bar plot
    :param data: list, data for 2 groups should be like [[1,2], [2,3]]
    :param labels: list, should match the number of items in data.
    """

    plt.rcParams["figure.figsize"] = (8, 6)

    if normalize:

        data_original = np.copy(data)
        data = np.array([[y / normalize_factor for y in x] for x in data])

    if len(labels) > 10:
        plt.rcParams["figure.figsize"] = (len(labels) / 3, 6)
    else:
        plt.rcParams["figure.figsize"] = (8, 6)

    # Use ggpot style
    plt.style.use("ggplot")

    fig, ax = plt.subplots()

    # Get bottom positions for stacked bar
    bottoms = np.cumsum(np.vstack((np.zeros(data.shape[1]), data)),
                        axis=0)[:-1]

    w = .8

    # Get x positions
    xpos = [x + (w / 2) for x in range(len(labels))]

    for c, d in enumerate(data):

        if len(data) <= 10:
            c1 = c if c < 9 else c - 10
            clr = cm.Vega10(c1, 1)
        else:
            c1 = c if c < 19 else c - 20
            clr = cm.Vega20c(c1, 1)

        if c == 0:
            bplot = ax.bar(xpos, d, w, color=clr, label=legend[c], alpha=.9)
        else:
            bplot = ax.bar(xpos, d, w, color=clr, label=legend[c], alpha=.9,
                           bottom=bottoms[c])

    # Set x labels
    plt.xticks([x + (w / 2) for x in xpos], labels, ha="right", rotation=45)

    # Set legend
    if legend:
        if len(legend) <= 4:
            cols = 1
        else:
            cols = len(legend) if len(legend) < 3 else 3
        borderpad = cols * -6
        lgd = plt.legend(loc=7, fancybox=True,
                         shadow=True, framealpha=.8, ncol=cols,
                         borderaxespad=borderpad)

    if ax_names:
        _add_labels(ax_names)

    # Generate table structure
    if table_header:
        table = [table_header]
    else:
        table = []
    if normalize:
        for i, lbl in enumerate(labels):
            table.append([lbl] + list(chain.from_iterable((int(x[i]), y[i])
                                      for x, y in zip(*[data_original, data]))))
    else:
        for i, lbl in enumerate(labels):
            table.append([lbl] + [int(x[i]) for x in data])

    return fig, lgd, table


@SetProps
def box_plot(data, labels=None, title=None, ax_names=None):
    """
    Creates a boxplot from a data series
    :param data: list, data to be plotted
    :param labels: list, x-axis labels
    :param title: string, plot title
    :param ax_names: list, first element for x-axis, second for y-axis
    """

    plt.rcParams["figure.figsize"] = (8, 6)

    if len(labels) > 10:
        plt.rcParams["figure.figsize"] = (len(labels) / 3, 6)
    else:
        plt.rcParams["figure.figsize"] = (8, 6)

    # Use ggpot style
    plt.style.use("ggplot")

    fig, ax = plt.subplots()

    bplot = ax.boxplot(data, patch_artist=True)

    # Change outline color, fill color and line width of the boxes
    for box in bplot["boxes"]:
        box.set(color="#7570b3", linewidth=2)
        box.set(facecolor="#1b9e77")

    # Change the color and line width of the whiskers
    for whisker in bplot["whiskers"]:
        whisker.set(color="#7570b3", linewidth=2)

    # Change color and line width of the caps
    for cap in bplot["caps"]:
        cap.set(color="#7570b3", linewidth=2)

    # Change color and line width of the medias
    for media in bplot["medians"]:
        media.set(color="#b2df8a", linewidth=2)

    # Change the style of the fliers and their fill
    for flier in bplot["fliers"]:
        flier.set(marker="o", color="#e7298a", alpha=0.5)

    if labels:
        ax.set_xticklabels(labels, rotation=45, ha="right")

    # Generate table structure
    table = [["", "1st Quartile", "Median", "3rd Quartile"]]

    for x, d in zip(labels, data):
        table.append([x, np.percentile(d, 25), np.median(d),
                      np.percentile(d, 75)])

    return fig, None, table


def histogram_smooth(data, ax_names=None, table_header=None, title=None,
                     legend=None):
    """
    Creates a smooth line plot with colored areas with the same distribution
    as an histogram
    :param data: list, with histogram data. To support multiple plots,
    provide each plot data as a list. If you want to produce two smooth
    histograms, provide a data list with two lists.
    :param title: string, title for the plot
    :param ax_names: list, names for the x-axis and y-axis, respectively
    :param table_header: list, header for the table
    """

    from scipy.interpolate import UnivariateSpline

    plt.rcParams["figure.figsize"] = (8, 6)

    plt.style.use("ggplot")

    fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True)

    # Start building table
    if table_header:
        table = [table_header]
    else:
        table = []

    for (p, d), axn, l in zip(enumerate(data), ax, legend):
        # Convert data into histogram
        y, x = np.histogram(d, 50)
        x = x[:-1] + (x[1] - x[0]) / 2
        # Fit a one-dimensional smoothing spline to git a given set of points
        # from the histogram
        f = UnivariateSpline(x, y, s=200)
        # Set title for subplot
        axn.set_title(l, color=clr_list[p], fontweight="bold")
        # Create subplot
        axn.plot(x, f(x), color=clr_list[p])
        # Color area below line lot
        axn.fill_between(x, f(x), alpha=0.5, color=clr_list[p])

    # If axis names are provided, add them to figure
    if ax_names:
        fig.text(0.5, 0.02, ax_names[0], ha="center", fontweight="bold",
                 color="grey", fontsize=14)
        fig.text(0.04, 0.5, ax_names[1], va="center", rotation="vertical",
                 fontweight='bold', color="grey", fontsize=14)

    # Build histogram with fixed bin range for table construction
    hist_data = []
    for d in data:
        y, x = np.histogram(d, 100, (0., 1.))
        hist_data.append(y)

    # Populate table
    for i, d in zip(np.arange(0, 1, 0.01), zip(*hist_data)):
        table.append([i] + list(d))

    return fig, None, table


@SetProps
def histogram_plot(data, title=None, ax_names=None, table_header=None,
                   real_bin_num=False):
    """
    Creates an histogram from data
    :param data: list
    :param title: string
    :param ax_names: list, first element for x-axis, second for y-axis
    :param table_header: list, containing the header of the table as the
    first element
    :param real_bin_num: boolean, if true, the histogram bins will be real
    numbers
    """

    # Use ggpot style
    plt.style.use("ggplot")

    plt.rcParams["figure.figsize"] = (8, 6)

    fig, ax = plt.subplots()

    c_data = Counter(data)
    max_data = max(c_data)
    min_data = min(c_data)
    if len(c_data) > 50:
        bins = int(len(c_data) / 10)
    else:
        bins = len(c_data)

    vals, b, _ = plt.hist(data, bins, histtype="stepfilled",
                                color=clr_list[0])

    plt.axvline(np.mean(data), linewidth=2, color="r", alpha=.8,
                linestyle="--")

    # Add cutom artist for legend
    mean_artist = plt.Line2D((0, 1), (0, 1), color="r", linestyle="--")

    lgd = ax.legend([mean_artist], ["Mean"], loc=0, frameon=True, fancybox=True,
                    shadow=True, framealpha=.8, fontsize="large")
    lgd.get_frame().set_facecolor("white")

    # Generate table structure
    if table_header:
        table = [table_header]
    else:
        table = []

    if real_bin_num:
        # If real_bin_num was set, then the bins in the table should be real
        # numbers
        vals, b = np.histogram(data, bins, (min_data, max_data + 1))

        c = 0
        for b, v in zip(b, vals):
            table.append(["{} - {}".format(int(c), int(b)), v])
            c = b + 1

        table.append(["{} - {}".format(int(c), max_data), vals[-1]])

    else:
        for p, val in zip(b, vals):
            table.append([p, val])

    return fig, lgd, table


@SetProps
def triangular_heat(data, labels, color_label=None, title=None):
    """
    Creates a triangular heatmap plot based on an array of triangular shape,
    or with a masked triangle
    :return:
    """

    plt.style.use("ggplot")

    plt.rcParams["figure.figsize"] = (8, 6)

    fig, ax = plt.subplots()

    cmap = cm.get_cmap("jet", 100)
    cmap.set_bad("w")

    heat = ax.imshow(data, interpolation="nearest", cmap=cmap)

    # Set ticks for label positioning
    ax.set_xticks(np.arange(data.shape[1]), minor=False)
    ax.set_yticks(np.arange(data.shape[0]), minor=False)
    # Set minor ticks that will make the grid
    ax.set_xticks(np.arange(data.shape[1]) + .5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]) - .5, minor=True)
    # Set x axis labels on top and y axis labels on right
    ax.xaxis.tick_top()
    # Set axis labels
    # Remove first entry of xlabel
    xlabel = ["" if x == 0 else labels[x] for x in range(len(labels))]
    ax.set_xticklabels(xlabel, rotation=45, ha="left")
    # Remove last entry of ylabel
    ylabel = ["" if x == len(labels) - 1 else labels[x] for x in
              range(len(labels))]
    ax.set_yticklabels(ylabel)
    # Remove major ticks
    plt.tick_params(axis="both", which="major", top="off", right="off",
                    left="off", labelsize=5)
    # Remove minor ticks
    plt.tick_params(axis="y", which="minor", right="off")

    plt.grid(True, which="minor")
    plt.grid(False, which="major")

    plt.gca().invert_xaxis()

    cbar = plt.colorbar(heat)
    if not color_label:
        cbar.set_label("Similarity proportion")
    else:
        cbar.set_label(color_label)

    # Generate table
    table = [[""] + [x for x in labels[::-1]]]
    for p, sp in enumerate(labels):
        table.append([sp] + list(data[p])[::-1])

    return fig, None, table


@SetProps
def outlier_densisty_dist(data, outliers, outliers_labels=None, ax_names=None,
                          title=None):
    """
    Creates a density distribution for data and highlights outliers
    :param data: 1D array containing data points
    :param outliers: 1D array containing the outliers
    :param outliers_labels: 1D array containing the name of the outliers.
    Must match outliers in size
    :param ax_names: list, with names for x and y axis, respectively
    """

    plt.rcParams["figure.figsize"] = (8, 6)

    fig, ax = plt.subplots()

    # Create density function
    sns.distplot(data, rug=True, hist=False, color="black")

    # Plot outliers
    ax.plot(outliers, np.zeros_like(outliers), "ro", clip_on=False,
            label="Outliers")

    # Create legend
    lgd = ax.legend(frameon=True, loc=2, fancybox=True, shadow=True,
                    framealpha=.8, prop={"weight": "bold"})

    if outliers_labels:
        table = [[os.path.basename(x)] for x in outliers_labels]
    else:
        table = None

    return fig, lgd, table


@SetProps
def sliding_window(data, window_size, ax_names=None, table_header=None,
                   title=None):
    """
    Creates a sliding window plot
    :param data:
    :param window_size:
    :param ax_names:
    :param table_header:
    """

    plt.style.use("ggplot")

    plt.rcParams["figure.figsize"] = (8, 6)

    data = np.array(data)
    x = np.arange(0, len(data) * window_size, window_size)

    table_data = []

    fig, ax = plt.subplots()

    ax.set_axis_bgcolor("#f2f2f2")

    xnew = np.linspace(x.min(), x.max(), 500)
    xsmooth = spline(x, data, xnew)

    # p = ax.plot(xnew, xsmooth)

    for i in range(0, len(xsmooth)):
        ax.plot(xnew[i:i + 2], xsmooth[i:i + 2], color=cm.rainbow(1 - (xsmooth[i] / float(max(
            xsmooth))), 1), linewidth=2.5)
        table_data.append([xnew[i], xsmooth[i]])

    ax.set_ylim([min(data) - min(data) * .1, max(data) + max(data) * .1])

    if ax_names:
        _add_labels(ax_names)

    if table_header:
        table = [table_header] + table_data

    else:
        table = table_data

    return fig, None, table


__author__ = "Diogo N. Silva"
__credits__ = ["Diogo N. Silva", "Tiago F. Jesus"]
