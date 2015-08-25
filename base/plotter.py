#!/usr/bin/env python3
# -*- coding: utf-8 -*-#
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
#
#  Author: Diogo N. Silva
#  Version:
#  Last update:

import matplotlib
# Using the agg backend prevent a segmentation fault when quiting the
# application
matplotlib.use("agg")

import matplotlib.pyplot as plt
import numpy as np


"""
Set of 10 easily distinguishable colors that will be used when generating plots
to ensure consistency in color usage. In case more than 10 colors are required
(ugh...) they will be randomly generated henceforth.
"""
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


def bar_plot(data_list, labels=None, title=None, ax_names=None,
             lgd_list=None, reverse_x=False):
    """
    Builds simple bar plot from a data_list
    :param data_list: list with data to be plotted.
    :param labels: list with x axis labels
    :param title: string, plot title
    :param ax_names: list. Names of the axis [yaxis_name, xaxis_name]
    :param reverse_x: Boolean, determines whether the x-axis is reversed or not
    """

    if len(labels) > 10:
        plt.rcParams["figure.figsize"] = (8 * len(labels) / 10, 6)
    else:
        plt.rcParams["figure.figsize"] = (8, 6)

    # Use ggpot style
    plt.style.use("ggplot")

    fig, ax = plt.subplots()

    # Create plot
    for i, d in enumerate(data_list):
        # Get color from 10 color list. If more than 10 colors are required,
        # randomly generate new ones
        if len(data_list) > 2:
            try:
                clr = clr_list[i]
            except IndexError:
                clr = np.random.rand(3, 1)
        else:
            clr = two_clr_list[i]

        # Create plot for first entry
        if d == data_list[0]:
            ax.bar(np.arange(len(d)), d, align="center", color=clr,
                   label=lgd_list[i] if lgd_list else None)
        # Add plots on top
        else:
            ax.bar(np.arange(len(d)), d, align="center", color=clr,
                   bottom=data_list[i - 1],
                   label=lgd_list[i] if lgd_list else None)

    # Set legend
    if lgd_list:
        lgd = plt.legend(bbox_to_anchor=(1, .5), loc="center left",
                         fancybox=True, shadow=True, framealpha=0.8)
    else:
        lgd = None

    # Set axys names
    if ax_names[0]:
        plt.xlabel(ax_names[0])
    if ax_names[1]:
        plt.ylabel(ax_names[1])

    # Set title
    if title:
        plt.title(title)

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

    return plt, lgd


def multi_bar_plot(data_list, labels=None, lgd_list=None):
    """
    General purpose multiple bar plot with custom layout. Returns a pyplot
    object.
    :param data_list: list with data to be plotted in list/array type in each
    entry. Data with two groups would be like [[1.23, .53], [1.55, .12]]
    :param labels: list, containing the labels
    :param lgd_list: list, The legend string for each data set in data_list
    """

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

    return plt, lgd


def interpolation_plot(data):
    """
    Creates a black and white interpolation plot from data, which must consist
    of a 0/1 matrix for absence/presence of taxa in genes
    :param data: numpy array of variable shape.
    """

    # Use ggpot style
    plt.style.use("ggplot")

    fig, ax = plt.subplots()

    # Setting the aspect ratio proportional to the data
    ar = len(data[0]) / 300

    ax.imshow(data, interpolation="none", cmap="Greys", aspect=ar)

    return plt


__author__ = 'diogo'
