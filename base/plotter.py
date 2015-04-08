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

import matplotlib.pyplot as plt
import numpy as np


def bar_plot(data, labels=None, title=None):
    """
    General purpose bar plot with custom layout. Returns a pyplot object
    :param data: list/array with data to be plotted
    :param labels: list, containing the labels
    """

    # Use ggpot style
    plt.style.use("ggplot")

    fig, ax = plt.subplots()

    # Create bar pot
    bplt = ax.bar(np.arange(len(data)), data, align="center")

    # Add labels
    if labels:
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="center")

    # Set label colors
    ax.tick_params(axis="x", colors=[.5, .5, .5])
    ax.tick_params(axis="y", colors=[.5, .5, .5])

    # Add title
    plt.title(title, size=23, color=[.4, .4, .4], loc="center")

    # Add text to the top of bar with the corresponding value
    for b in bplt:

        xpos = b.get_x() + b.get_width() / 2
        ypos = b.get_height() * .92
        ax.text(xpos, ypos, b.get_height(), horizontalalignment="center",
                weight='bold', color="white", size=14)

    # Automatically adjust figure size to accomodate all labels
    plt.tight_layout()

    return plt

__author__ = 'diogo'
