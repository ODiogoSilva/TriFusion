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
import functools

# Set of 10 easily distinguishable colors that will be used when generating
# plots to ensure consistency in color usage. In case more than 10 colors
#  are required (ugh...) they will be randomly generated henceforth.

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


def set_props(func):
    """Decorator used to set general plot attributes

    This decorator is use to automatically set several plot properties based
    on the arguments provided to the plotting functions. In this way,
    there is no need to write repetitive code in each function. The accepted
    arguments are:

      - ax_names : list.
        Axis names. First element is x-axis label, second is y-axis label
      - title : str.
        Title of the plot
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        plt.clf()
        plt.close()

        res = func(*args, **kwargs)

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

    return wrapper


@set_props
def scatter_plot(data, correlation=False, ax_names=None, table_header=None,
                 title=None):
    """Generates a scatter plot from a 2D array.

    Builds a scatter plot from a 2D array. Also calculates the
    correlation coefficient if requested by the `correlation` argument.

    Parameters
    ----------
    data : numpy.array
        2D array containing the x and y data points.
    correlation : bool
        If True, the spearman's rank correlation coefficient is calculated
        and added to the plot as an annotation
    ax_names : list
        List with the labels for the x-axis (first element) and y-axis
        (second element).
    table_header : list
        List with the header of the table object. Each element represents
        a column.
    title : str
        Title of the plot.

    Returns
    -------
    fig : matplotlib.pyplot.Figure
        Figure object of the plot.
    _ : None
        Placeholder for the legend object. Not used here but assures
        consistency across other plotting methods.
    table : list
        Table data in list format. Each item in the list corresponds to a
        table row.
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


@set_props
def bar_plot(data, labels=None, title=None, ax_names=None,
             lgd_list=None, reverse_x=False, table_header=None):
    """Creates a bar plot from a `data` array.

    If a multi-dimensional array is provided, it will create a stacked bar
    plot.

    Parameters
    ----------
    data : numpy.array
        Single or multi-dimensional array with plot data.
    labels : list
        List of xtick labels. Should have the same length as `data`.
    title : str
        Title of the plot.
    ax_names : list
        List with the labels for the x-axis (first element) and y-axis
        (second element).
    lgd_list : list
        For categorical plots, provide the label of each category.
    reverse_x : bool
        If True, reverse the x-axis orientation.
    table_header : list
        List with the header of the table object. Each element represents
        a column.

    Returns
    -------
    fig : matplotlib.Figure
        Figure object of the plot.
    lgd : matplotlib.Legend
        Legend object of the plot.
    table : list
        Table data in list format. Each item in the list corresponds to a
        table row.
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


@set_props
def multi_bar_plot(data, labels=None, lgd_list=None, title=None,
                   ax_names=None):
    """Creates a multiple bar plot.

    Creates a multiple bar plot from a multi-dimensional array. The bar
    plots from each array will be grouped and displayed side by side.

    Parameters
    ----------
    data : numpy.array
        Single or multi-dimensional array.
    labels : list
        List of xtick labels. Should have the same length as `data`.
    title : str
        Title of the plot.
    lgd_list : list
        For categorical plots, provide the label of each category.
    ax_names : list
        List with the labels for the x-axis (first element) and y-axis
        (second element).

    Returns
    -------
    fig : matplotlib.Figure
        Figure object of the plot.
    lgd : matplotlib.Legend
        Legend object of the plot.
    table : list
        Table data in list format. Each item in the list corresponds to a
        table row.
    """

    plt.rcParams["figure.figsize"] = (8, 6)

    # Use ggpot style
    plt.style.use("ggplot")

    fig, ax = plt.subplots()

    # Set margin between bar groups
    margin = 0.1

    # Determine bar group width according to the number of data lists
    w = (1. - 2. * margin) / len(data)
    # Determine max plt height to calculate text label space
    max_height = max((x for y in data for x in y))

    # Create bar plots
    for i, d in enumerate(data):
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


@set_props
def interpolation_plot(data, title=None, ax_names=None):
    """Creates black and white interpolation plot

    Creates a black and white interpolation plot from data, which must consist
    of a 0/1 matrix for absence/presence of taxa in genes.

    Parameters
    ----------
    data : numpy.array
        Single or multi-dimensional array with plot data.
    title : str
        Title of the plot.
    ax_names : list
        List with the labels for the x-axis (first element) and y-axis
        (second element).

    Returns
    -------
    fig : matplotlib.pyplot.Figure
        Figure object of the plot.
    _ : None
        Placeholder for the legend object. Not used here but assures
        consistency across other plotting methods.
    _ : None
        Placeholder for the table header list. Not used here but assures
        consistency across other plotting methods.
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


@set_props
def stacked_bar_plot(data, labels, legend=None, table_header=None, title=None,
                     ax_names=None, normalize=False, normalize_factor=None):
    """Creates a stacked bar plot.

    Parameters
    ----------
    data : numpy.array
        Multi-dimensional array.
    labels : list
        List of xtick labels. Should have the same length as `data`.
    title : str
        Title of the plot.
    table_header : list
        List with the header of the table object. Each element represents
        a column.
    title : str
        Title of the plot.
    ax_names : list
        List with the labels for the x-axis (first element) and y-axis
        (second element).
    normalize : bool
        If True, values of the `data` array will be normalized by the
        `normalize_factor`
    normalize_factor : int or float
        Number used to normalize values of the `data` array.

    Returns
    -------
    fig : matplotlib.Figure
        Figure object of the plot.
    lgd : matplotlib.Legend
        Legend object of the plot.
    table : list
        Table data in list format. Each item in the list corresponds to a
        table row.
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


@set_props
def box_plot(data, labels=None, title=None, ax_names=None):
    """Creates box (whisker) plot.

    Parameters
    ----------
    data : numpy.array
        Single or multi-dimensional array with plot data.
    labels : list
        List of xtick labels. Should have the same length as `data`.
    title : str
        Title of the plot.
    ax_names : list
        List with the labels for the x-axis (first element) and y-axis
        (second element).

    Returns
    -------
    fig : matplotlib.pyplot.Figure
        Figure object of the plot.
    _ : None
        Placeholder for the legend object. Not used here but assures
        consistency across other plotting methods.
    table : list
        Table data in list format. Each item in the list corresponds to a
        table row.
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
    """Creates a smooth histogram-like plot.

    Creates a smooth line plot with colored areas with the same distribution
    as an histogram. It supports a multi-dimensional array, in which case
    each array will be used to create vertically stacked subplots.

    Parameters
    ----------
    data : numpy.array
        Array with plot data.
    ax_names : list
        List with the labels for the x-axis (first element) and y-axis
        (second element).
    table_header : list
        List with the header of the table object. Each element represents
        a column.
    title : str
        Title of the plot.
    legend : list
        If using a multi-dimensional array, provide the name of each
        subplot.

    Returns
    -------
    fig : matplotlib.pyplot.Figure
        Figure object of the plot.
    _ : None
        Placeholder for the legend object. Not used here but assures
        consistency across other plotting methods.
    table : list
        Table data in list format. Each item in the list corresponds to a
        table row.
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


@set_props
def histogram_plot(data, title=None, ax_names=None, table_header=None,
                   real_bin_num=False):
    """Creates an histogram from data.

    Parameters
    ----------
    data : numpy.array
        Array with plot data.
    title : str
        Title of the plot.
    ax_names : list
        List with the labels for the x-axis (first element) and y-axis
        (second element).
    table_header : list
        List with the header of the table object. Each element represents
        a column.
    real_bin_num : bool
        If True, then the table data will be forced to be in real numbers.

    Returns
    -------
    fig : matplotlib.Figure
        Figure object of the plot.
    lgd : matplotlib.Legend
        Legend object of the plot.
    table : list
        Table data in list format. Each item in the list corresponds to a
        table row.
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


@set_props
def triangular_heat(data, labels, color_label=None, title=None):
    """Creates a triangular heatmap plot.

    Parameters
    ----------
    data : numpy.array
        Triangular array with plot data.
    labels : list
        List of xtick labels. Should have the same length as `data`.
    title : str
        Title of the plot.
    color_label : str
        Label for colorbar.

    Returns
    -------
    fig : matplotlib.pyplot.Figure
        Figure object of the plot.
    _ : None
        Placeholder for the legend object. Not used here but assures
        consistency across other plotting methods.
    table : list
        Table data in list format. Each item in the list corresponds to a
        table row.
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


@set_props
def outlier_densisty_dist(data, outliers, outliers_labels=None, ax_names=None,
                          title=None):
    """Creates a density distribution for outlier plots.

    Parameters
    ----------
    data : numpy.array
        1D array containing data points.
    outliers : numpy.array
        1D array containing the outliers.
    outliers_labels : list or numpy.array
        1D array containing the labels for each outlier.
    title : str
        Title of the plot.

    Returns
    -------
    fig : matplotlib.Figure
        Figure object of the plot.
    lgd : matplotlib.Legend
        Legend object of the plot.
    table : list
        Table data in list format. Each item in the list corresponds to a
        table row.
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


@set_props
def sliding_window(data, window_size, ax_names=None, table_header=None,
                   title=None):
    """Creates a sliding window plot.

    Parameters
    ----------
    data : numpy.array
        Array with plot data.
    window_size : int
        Length of window size for sliding window.
    table_header : list
        List with the header of the table object. Each element represents
        a column.
    title : str
        Title of the plot.
    ax_names : list
        List with the labels for the x-axis (first element) and y-axis
        (second element).

    Returns
    -------
    fig : matplotlib.pyplot.Figure
        Figure object of the plot.
    _ : None
        Placeholder for the legend object. Not used here but assures
        consistency across other plotting methods.
    table : list
        Table data in list format. Each item in the list corresponds to a
        table row.
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

    if table_header:
        table = [table_header] + table_data

    else:
        table = table_data

    return fig, None, table


__author__ = "Diogo N. Silva"
__credits__ = ["Diogo N. Silva", "Tiago F. Jesus"]
