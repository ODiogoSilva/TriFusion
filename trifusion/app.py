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

"""
Nomenclature guide for :class:`.TriFusionApp`.
==============================================

Given the large number of methods needed to give functionality to
the app, this nomenclature guide was created to aid in the naming
of new methods so that the code can be more easily browsed and
understood. Note that this guide only targets methods that
perform similar tasks and, therefore, can be grouped by a common
prefix name. Other methods that perform more unique operations
may have different names.

Method's names will be given based on their main operation and
specific task. For example, a method in charge of toggle the side
panel, should be named "toggle_sidepanel", being "toggle" the
common prefix and "sidepanel" the keyword linked ot the specific
task.

1. Toggles.

"toggle_[specific_task]", e.g. "toggle_sidepanel"

Methods use to toggle certain widgets or values/attributes.

2. Dialogues.

"dialog_[specific_task]", e.g. "dialog_format"

Methods that generate dialogues throughout the app, usually in
the form of popups

3. Populating methods.

"populate_[specific_task]", e.g., "populate_input_files"

Methods that populate certain widgets, usually gridlayouts, with
other widgets

4. Add/Remove

"add_[specific_task]", e.g., "add_bookmark"
"remove_[specific_task]", e.g., "remove_taxa_group"

Methods that add or remove widgets, usually buttons/togglebuttons,
from other widgets

5. Saves.

"save_[specific_task]", e.g., "save_file"

Methods that save specific settings from the several options of the
sapp

6. Updates.

"update_[specific_task]", e.g., "update_tabs"

Wrapper methods used to update several attributes or widgets of the
app

7. Checks.

"check_[specific_task]", e.g., "check_filters"

Methods that perform some kind of sanity checks to user input data

8. Unique operations

[specific_task]_[unique_operation], e.g., "sidepanel_animation"

When the method performs a unique operations, the specific_task
should prefix the name of the method.
"""

# Standard libraries imports
from os.path import dirname, exists, splitext
from collections import OrderedDict
from functools import partial
from copy import deepcopy
from scipy import misc
import multiprocessing
import matplotlib.patches as patches
import subprocess
import platform
import logging
import ctypes
import pickle
import urllib
import string
import signal
import shutil
import Queue
import time
import stat
import sys
import os
import threading

# Move to Application's directory. This is a way of avoiding encoding
# issues when the full path to the application's directory contains
# non-ASCII characters. This way, the cur_dir attribute will be always
# set to "."
called_dir = os.getcwd()
if getattr(sys, "frozen", False):
    os.chdir(dirname(sys.executable))
elif __file__:
    os.chdir(dirname(os.path.abspath(__file__)))

# freeze_support must be called here so that multiprocessing work
# correctly on windows
multiprocessing.freeze_support()

# Bypass default argument handling from kivy
os.environ["KIVY_NO_ARGS"] = "1"

# Comment this section until the EventLoop for sphinx building.
# Kivy imports
from kivy.config import Config

# Sets some kivy configurations before creating main window
Config.set("kivy", "desktop", 1)
Config.set("kivy", "exit_on_escape", 0)
Config.set("graphics", "resizable", 1)
Config.set("graphics", "fullscreen", 0)
Config.set("graphics", "minimum_height", 700)
Config.set("graphics", "minimum_width", 1000)
Config.set("graphics", "height", 700)
Config.set("graphics", "width", 1000)
Config.set("graphics", "multisamples", 4)
Config.set("input", "mouse", "mouse, disable_multitouch")

# Force creation of main window
from kivy.base import EventLoop
EventLoop.ensure_window()

from kivy.app import App
from kivy.animation import Animation
from kivy.uix.checkbox import CheckBox
from kivy.uix.scrollview import ScrollView
from kivy.lang import Builder
from kivy.properties import ListProperty, DictProperty
from kivy.clock import Clock
from kivy.uix.treeview import TreeViewLabel
from kivy.logger import Logger

# Local TriFusion imports
try:
    from __init__ import __version__, __build__
    import data.resources.theme.default as tm
    from ortho import protein2dna
    from process.base import Base
    from process.data import Partitions
    from process.error_handling import EmptyAlignment, KillByUser
    from process.sequence import AlignmentList
    from data.resources.info_data import orthology_plots, \
        informative_storage
    from data.resources.background_tasks import remove_tmp, \
        load_proc, get_stats_summary, background_process, \
        background_export_groups, orto_execution, process_execution, \
        orto_update_filters, get_active_group, get_orto_data, get_stats_data
    from data.resources.custom_widgets import *
    from data.resources.stats import *
    from base.plotter import scatter_plot, bar_plot, \
        interpolation_plot, stacked_bar_plot, box_plot, histogram_smooth, \
        histogram_plot, triangular_heat, outlier_densisty_dist, \
        sliding_window
    from base.html_creator import HtmlTemplate
    from ortho.OrthomclToolbox import MultiGroupsLight
except ImportError:
    from trifusion import __version__, __build__
    import trifusion.data.resources.theme.default as tm
    from trifusion.ortho import protein2dna
    from trifusion.process.base import Base
    from trifusion.process.data import Partitions
    from trifusion.process.error_handling import EmptyAlignment, \
        KillByUser
    from trifusion.process.sequence import AlignmentList
    from trifusion.data.resources.info_data import orthology_plots, \
        informative_storage
    from trifusion.data.resources.background_tasks import remove_tmp, \
        load_proc, get_stats_summary, background_process, \
        background_export_groups, orto_execution, process_execution, \
        orto_update_filters, get_active_group, get_orto_data, get_stats_data
    from trifusion.data.resources.custom_widgets import *
    from trifusion.data.resources.stats import *
    from trifusion.base.plotter import scatter_plot, bar_plot, \
        interpolation_plot, stacked_bar_plot, box_plot, histogram_smooth, \
        histogram_plot, triangular_heat, outlier_densisty_dist, \
        sliding_window
    from trifusion.base.html_creator import HtmlTemplate
    from trifusion.ortho.OrthomclToolbox import MultiGroupsLight

###################################
# Modifications to kivy source code
###################################

# MEMO 1
# scatter.py on_touch_up function was modified to prevent a bug from
# crashing  the app. Before the fix, when right-clicking in a
# ScatterLayout and then  left-clicking in the generate ball would crash
# the app with a KeyError  exception. The fix on lines 595-598 handles
# this exception by not deleting a touch from _last_touch_pos
# dictionary, since the key is not there anymore.
# MEMO 2
# scrollview.py _change_touch_mode function was modified to prevent a bug
# from # crashing the app. Before the fix, sometimes the app would crash
# due to a KeyError exception on line 920. I did not detect any
# specific pattern for this  error and I could never replicate the bug,
# but every now and then the app  would crash. The fix on lines 920-923
# handles the KeyError exception by  returning the function.
# MEMO 3
# To avoid issues with non-ascii characters on the path of the binary
# file, some changes were performed on matplotlib and kivy files.
# .:Matplotlib:.
#   .: matplotlib/__init__.py :.
#      [l324](wrap()) Wrap the self.report() method around a try
#      statement and handle all exceptions attributing False to spoke
#      [l822](matplotlib_fname) hardcode the join() method with "mpl-data"
#   .: matplotlib/font_manager.py:.
#      [l1026](__init__()) Replace the "datapath" value of the rcParams
#      dictionary wit a relative path ->
#      rcParams["datapath"] = os.path.relpath(rcParams["datapath"])
#   .: matplotlib/mathtext.py:.
#      [l1013](StandardPsFonts()) Wrap the get_data_path() function call
#      with a relative path ->
#      basepath = os.path.join( os.path.relpath(get_data_path()),
#                              'fonts', 'afm' )
#   .: matplotlib/style/core.py:.
#      [l30] Same as previous ->
#      BASE_LIBRARY_PATH = os.path.join(
#          os.path.relpath(mpl.get_data_path()), 'stylelib')
# .:Kivy:.
#   .: kivy/atlas.py :.
#      [l206] Wrap the dirname(filename) in a relative path ->
#      d = os.path.relpath(dirname(filename))
#    .: kivy/core/text/text_sdl2.py :.
#      [l21](_get_font_id()) Replace the "font_name_r" value of the
#      self.options dictionary with a relative path ->
#      self.options['font_name_r'] = os.path.relpath(
#          self.options['font_name_r'])
#    .: kivy/core/window/window_sdl2.py :.
#       [l339](set_icon()) Wrap the filename variable in a relative path ->
#       self._win.set_window_icon(str(os.path.relpath(filename)))

# ==============================================================================
#                                  EXCEPTIONS
# ==============================================================================

MCL_FILE = None


class InputTypeError(Exception):
    pass


class TriFusionApp(App):

    #######################
    #
    # GUI RELATED VARIABLES
    #
    #######################

    # Referecence to blue and red colors
    _blue = tm.c_app_blue
    _red = tm.c_app_red

    show_side_panel = BooleanProperty(False)
    """
    bool controlling the toggling of main headers. Is set to True
    when the sidepanel is open and False when not.
    """
    lock_side_panel = BooleanProperty(False)
    """bool that prevents the animation of the sidepanel when True."""

    screen = None
    """
    Currently open Screen object.
    """

    # Variable containing screen names
    screen_names = ListProperty()
    """List of basename screen names."""
    available_screens = ListProperty()
    """List of paths to screen files."""
    loaded_screens = {}
    """
    Dictionary that maps the screen name to a Screen object. Entries are
    only populated when Screens are loaded for the first time. This
    allows the state of each screen to be saved.
    """
    plot_screens = ListProperty()
    """List of screens dedicated for plots."""

    # Attributes to know current and previous screen
    current_screen = StringProperty()
    """String with the name of the current Screen object."""
    previous_screen = StringProperty()
    """String with the name of the last Screen object."""

    index = NumericProperty(-1)
    """Integer with the index of the current screen. The index is
    retrieved from the :attr:`~.TriFusionApp.screen_names` attribute."""

    temp_dir = StringProperty()
    """String with path to temporary directory."""
    log_file = StringProperty()
    """String with path to TriFusion's log file."""
    bm_file = StringProperty()
    """String with path to bookmarks file."""
    projects_file = StringProperty()
    """String with path to projects file."""

    cur_dir = "."
    """
    String with the current working directory. Used to fetch the kivy
    screen files.
    """

    file_list = []
    """
    List with the paths of alignment files loaded into TriFusion.
    This attribute should not be modified unless the files are removed
    from the application.
    """
    active_file_list = ListProperty()
    """
    List with the paths of the active alignment files loaded into TriFusion.
    Can be changed during the session.
    """

    filename_map = DictProperty()
    """
    Dictionary mapping the basename of alignment files to their full path.
    Some parts of TriFusion need only the basename of the alignment files,
    so this attributes allows the conversion.
    """

    active_taxa_list = ListProperty()
    """
    List with the active taxa loaded into TriFusion. Can be changed
    during the session.
    """

    home_path = unicode(expanduser("~"))
    """String with path to home directory."""

    bookmarks = [[], {}]
    """
    List with information of the bookmarks for the filechoosers. The first
    element is a list with the full path of the bookmarks. The second
    element is a dictionary mapping the full path to the bookmark name.
    """
    bookmarks_bt = []
    """
    List of bookmark Button widgets. Used for mouse over purposes.
    """

    drag_files = []
    """
    List that stores the paths of files drag and dropped into TriFusion's
    window.
    """
    drag_c = 0
    """
    Integer that works as a counter of the number of files that were dropped
    into TriFusion's window.
    """

    _popup = ObjectProperty(None)
    """Reference to Popup object, if open."""
    _subpopup = ObjectProperty(None)
    """Reference to sub Popup object, if open."""
    _exit_popup = ObjectProperty(None)
    """Reference to exit Popup object, if open"""

    # Dictionary containing the values for the main process operations
    main_operations = DictProperty({
        "concatenation": False,
        "conversion": False,
        "reverse_concatenation": False})
    """
    Dictionary mapping the main operations to a bool value. The active
    main operation, if any, has a True value. Only one True value can
    exist for the main operations.
    """

    # Dictionary containing all values of the switches and checkboxes in the
    # process screen
    secondary_operations = OrderedDict([("collapse", False),
                                        ("filter", False),
                                        ("gcoder", False),
                                        ("consensus", False)])
    """
    Dictionary mapping the secondary operations to a bool value. Active
    secondary operations have a True value.
    """

    secondary_options = DictProperty([("interleave", False),
                                      ("zorro", False),
                                      ("collapse_file", False),
                                      ("collapse_filter", False),
                                      ("filter_file", False),
                                      ("taxa_filter", False),
                                      ("codon_filter", False),
                                      ("gap_filter", False),
                                      ("minimum_taxa_filter", False),
                                      ("variation_filter", False),
                                      ("gcoder_file", False),
                                      ("consensus_file", False),
                                      ("consensus_single", False)])
    """
    Dictionary mapping secondary options to a bool value. Active secondary
    options have a True value.
    """

    file_overwrite = None
    """
    Attribute controlling whether output files should be overwritten
    or skipped. Has a None value when inactive, "overwrite" when there is
    an instruction to overwrite, and "skip" when there is a instruction to
    skip.
    """
    file_apply_all = False
    """
    bool attribute that determines whether the behaviour of
    :attr:`~.TriFusionApp.file_overwrite` should be applied to all files
    (True) or not (False).
    """

    ortho_search_options = None
    """
    Reference to the
    :class:`trifusion.data.resources.custom_widgets.OrthologySearchGrid
    object, containing the additional option's widgets of the Orthology
    Search screen.
    """
    orto_search_height = None
    """
    Integer with the height property of the
    :attr:`~.TriFusionApp.ortho_search_options` object.
    """

    # Attribute for the gridlayout widget that will contain all main options
    # for the process module
    process_grid_wgt = None
    """
    Reference to the
    :class:`trifusion.data.resources.custom_widgets.ProcessGeneral object
    with the general options widgets of the Process screen.
    """
    process_options = None
    """
    Reference to the 
    :class:`trifusion.data.resources.custom_widgets.AdditionalProcessContents
    object with the additional options widgets of the Process screen.
    """
    process_height = None
    """
    Integer with the height property of the
    :attr:`~.TriFusionApp.process_options` object.
    """

    stats_summary = None
    """
    Reference to the
    :class:`trifusion.data.resources.custom_widgets.StatsSummary object
    containing the widgets of the overall summary statistics in the Statistics
    screen.
    """
    stats_table = None
    """
    Reference to the 
    :class:`trifusion.data.resources.custom_widgets.GeneTable object
    containing the widgets of the gene table summary statistics in the
    Statistics screen.
    """
    MAX_TABLE_N = 50
    """
    Integer with the maximum number of rows allowed for the gene table
    in the :class:`trifusion.data.resources.custom_widgets.GeneTable object.
    """
    gene_table_selection = []
    """
    Pandas DataFrame with the selection of gene summary statistics for the
    :class:`trifusion.data.resources.custom_widgets.GeneTable object.
    """
    gene_master_table = []
    """
    List with the table data from the
    :class:`trifusion.data.resources.custom_widgets.GeneTable object.
    """

    # Attribute for the widget containing the treeview showing the
    # operations queue
    operation_tv = ObjectProperty(None)
    """TreeView object containing the queue of operations."""
    main_nodes = DictProperty()
    """
    Dictionary that maps string identifiers of main nodes to
    those Node objects.
    """

    # Attributes containing plot related elements
    plot_backups = {}
    """
    Dictionary that maps Statistics plot analyses identifiers to their
    respective table data in a list object. This object is only populated
    when those analyses are executed.
    """
    #TODO: Check memory usage of plot_backups
    current_plt_idx = ListProperty(None)
    """
    List that contains two elements. The first is the string identifier
    of the current plot analyses. The second is the list of arguments
    that should be passed to the respective plotting method.
    """
    current_plot = ObjectProperty(None, allownone=True)
    """
    matplotlib.Figure object of the current plot.
    """
    current_lgd = None
    """matplotlib.Legend object associated to the current plot. Can be None."""
    current_table = ObjectProperty(None, allownone=True)
    """List with table data associated with the current plot."""

    # Patch attribute
    plt_patch = None
    """Stores matplotlib.patches objects that may be added to plots."""
    orto_plt_method = {
        "Taxa distribution":
        [bar_plot, "Species_distribution.png"],
        "Gene copy distribution":
        [bar_plot, "Gene_copy_distribution.png"],
        "Taxa coverage":
        [bar_plot, "Species_coverage.png"],
        "Taxa gene copies":
        [bar_plot, "Species_copy_number.png"]
    }
    """
    Dictionary that maps the string identifier of orthology plot analyses
    to a list. This list has two elements. The first is the plotting
    method defined in :mod:`trifusion.base.plotter`, and the second is
    the name of the temporary figure file that will be generated.
    """
    stats_plt_method = {
            "Gene occupancy":
            [interpolation_plot, "gene_occupancy.png"],
            "Distribution of missing data sp":
            [stacked_bar_plot, "missing_data_distribution_sp.png"],
            "Distribution of missing data":
            [histogram_smooth, "missing_data_distribution.png"],
            "Distribution of missing orthologs":
            [bar_plot, "missing_gene_distribution.png"],
            "Distribution of missing orthologs avg":
            [histogram_plot, "missing_gene_distribution_avg.png"],
            "Cumulative distribution of missing genes":
            [bar_plot, "cumulative_distribution_missing_genes.png"],
            "Distribution of sequence size":
            [box_plot, "avg_seqsize_species.png"],
            "Distribution of sequence size all":
            [histogram_plot, "avg_seqsize.png"],
            "Proportion of nucleotides or residues":
            [bar_plot, "char_proportions.png"],
            "Proportion of nucleotides or residues sp":
            [stacked_bar_plot, "char_proportions_sp.png"],
            "Pairwise sequence similarity":
            [histogram_plot, "similarity_distribution.png"],
            "Pairwise sequence similarity sp":
            [triangular_heat, "similarity_distribution_sp.png"],
            "Pairwise sequence similarity gn":
            [sliding_window, "similarity_distribution_gn.png"],
            "Segregating sites":
            [histogram_plot, "segregating_sites.png"],
            "Segregating sites prop":
            [histogram_plot, "segregating_sites_prop.png"],
            "Segregating sites sp":
            [triangular_heat, "segregating_sites_sp.png"],
            "Segregating sites gn":
            [sliding_window, "segregating_sites_gn.png"],
            "Alignment length/Polymorphism correlation":
            [scatter_plot, "length_polymorphism_correlation.png"],
            "Distribution of taxa frequency":
            [histogram_plot, "distribution_taxa_frequency.png"],
            "Allele Frequency Spectrum":
            [histogram_plot, "allele_frequency_spectrum.png"],
            "Allele Frequency Spectrum prop":
            [histogram_plot, "allele_frequency_spectrum_prop.png"],
            "Allele Frequency Spectrum gn":
            [histogram_plot, "allele_frequency_spectrum_gn.png"],
            "Missing data outliers":
            [outlier_densisty_dist, "Missing_data_outliers.png"],
            "Missing data outliers sp":
            [outlier_densisty_dist, "Missing_data_outliers_sp.png"],
            "Segregating sites outliers":
            [outlier_densisty_dist, "Segregating_sites_outliers.png"],
            "Segregating sites outliers sp":
            [outlier_densisty_dist, "Segregating_sites_outliers_sp.png"],
            "Sequence size outliers":
            [outlier_densisty_dist, "Sequence_size_outliers.png"],
            "Sequence size outliers sp":
            [outlier_densisty_dist, "Sequence_size_outliers_sp.png"]
        }
    """
    Dictionary that maps the string identifier of statistics plot analyses
    to a list. This list has two elements. The first is the plotting
    method defined in :mod:`trifusion.base.plotter`, and the second is
    the name of the temporary figure file that will be generated.
    """

    previous_sets = {"Files": [], "Taxa": [], "Stats": []}
    """
    Dictionary storing previous data sets defined in TriFusion. This is
    used to evaluate whether the plot methods should be executed or ignored.
    The Stats key refers to the previous active data sets when the stats
    summary statistics overview was performed. It should contain two lists,
    the first with the file set and the second with the taxa set.

        previous_sets = {"Files": [], "Taxa": [], "Stats": []}
    """

    previous_stats_toggle = None
    """
    Reference to
    :class:`~trifusion.data.resources.custom_widgets.StatsToggleWgt object.
    It's used to restore the widget to the Statistics screen when exiting
    and entering it when a plot is loaded.
    """

    # Attributes for storing taxa and file buttons for side panel. These
    #  will be used when search for files/taxa and for loading only
    # button subsets for very large data sets. Each list element
    # pertains a single file/taxon and it will be a tupple containing the
    #  main button, information button and remove button.
    sp_file_bts = ListProperty()
    """
    List of widgets from individual files in the "Files" tab of the
    sidepanel. This attribute is used when searching for files and for
    loading only Button subsets of very large data sets. Each element
    will be a tuple containing the main Button, the information Button,
    and the remove Button widgets.
    """
    sp_taxa_bts = ListProperty()
    """
    List of widgets from individual taxa in the "Taxa" tab of the
    sidepanel. This attribute is used when searching for taxa and for
    loading only Button subsets of very large data sets. Each element
    will be a tuple containing the main Button, the information Button,
    and the remove Button widgets.
    """
    sp_partition_bts = ListProperty()
    """
    List of widgets from individual partitions in the "Partitions" tab of the
    sidepanel. This attribute is used when searching for partitions and for
    loading only Button subsets of very large data sets. Each element
    will be a tuple containing the main Button, the counter Button,
    and the edit Button widgets.
    """

    # Attributes that control the amount of taxa/file buttons showing at the
    # side panel. To avoid staggering the app with tons of buttons, a
    # maximum number of buttons showing initially is set. More buttons
    # can be later added.
    MAX_FILE_BUTTON = NumericProperty(20)
    """
    Integer with the maximum number of file Buttons allowed in the "Files"
    tab of the sidepanel. Prevents the loading of all file Buttons in
    very large datasets.
    """
    count_files = NumericProperty(0)
    """
    Integer that servers as a counter of the number of file Buttons currently
    loaded into the "Files" tab of the sidepanel.
    """
    MAX_PARTITION_BUTTON = NumericProperty(20)
    """
    Integer with the maximum number of partition Buttons allowed in the
    "Partitions" tab of the sidepanel. Prevents the loading of all partition
    Buttons in very large datasets.
    """
    count_partitions = NumericProperty(0)
    """
    Integer that servers as a counter of the number of partition Buttons
    currently loaded into the "Partitions" tab of the sidepanel.
    """

    mouse_over_bts = DictProperty({
        "Files": [],
        "Taxa": [],
        "Partitions": []})
    """
    Attributes storing the ToggleButton objects from Taxa/File/Partitions
    tabs in the sidepanel. Used mostly for mouse over events.
    """

    previous_mouse_over = StringProperty("")
    """
    The Button text property of the previous mouse over event. This will
    allow the assessment of whether the current mouse collision is for the
    same button (in which case the mouse over will not be triggered) or for
    a different button (in which case the mouse over is triggered).
    """

    mouse_over_ready = BooleanProperty(True)
    """
    bool attribute that is a locking mechanism of the mouse over event.
    When there is a scheduled event for a mouse over this attribute is set
    to False, which prevents further events from being scheduled in the
    meantime. When the scheduled event is dispatched, the lock is released
    and it returns to True.
    """

    old_mouse_over = None
    """
    Label widget of the previous mouse over. Used so that this label can be
    removed.
    """
    fancy_bt = None
    """
    Reference to the 
    :class:`~trifusion.data.resources.custom_widgets.FancyButton` object
    from the current mouse over.
    """
    touch = None
    """
    Touch object. Can be used to query the type of touch on certain events.
    """

    removable_media = []
    """
    List attribute that stores paths of currently active removable media.
    """

    arrow_block = BooleanProperty(False)
    """
    bool attribute that ocks arrow keys keybindings when a text input is
    focused.
    """

    # Whether SidePanel's More options dialog is active or not
    sp_moreopts = BooleanProperty(False)
    """
    bool attribute of whether the side panel's MoreOptions dialog is active
     (True) or not (False).
    """

    # Whether a FancyDropDown is... well... dropped
    fancy_dropped = BooleanProperty(False)
    """
    bool attribute of whether the
    :class:'trifusion.data.resources.custom_widgets.FancyDropDown object
    if active (True) or not (False).
    """

    last_sp_bt = {"Files": None, "Taxa": None, "Partitions": None}
    """
    Dictionary that stores the last selected Button object in the sidepanel
    for Files, Taxa and Partitions.
    
        last_sp_bt = {"Files": None, "Taxa": None, "Partitions": None}
    """

    is_control_pressed = BooleanProperty(False)
    """
    bool attribute that stores information on whether the control key
    is being pressed (True) or not (False).
    """

    is_shift_pressed = BooleanProperty(False)
    """
    bool attribute that stores information on whether the shift key
    is being pressed (True) or not (False).
    """

    prev_tb = None
    """
    Attribute that stores the previous Button for general purpose multi
    selection of toggle buttons using shift clicking
    """

    mouse_position = ListProperty()
    """
    List of the mouse position coordinates in [xpos, ypos].
    """

    lock_stats = BooleanProperty(False)
    """
    bool attribute that determines whether statistics background processes
    should be interrupted (True) or not (False).
    """

    # Attributes used to determine is a specific background process
    # should be terminated by the user. All values are set to True by
    # default and their values only change on the start up of the
    # appropriate function.

    terminate_stats = BooleanProperty(True)
    """
    bool attribute used to determine if a summary statists background
    process should be terminated by the user. Value is True by default and
    they only change on the start up of the appropriate background function.
    """

    terminate_process_exec = BooleanProperty(True)
    """
    bool attribute used to determine if a process execution background
    process should be terminated by the user. Value is True by default and
    they only change on the start up of the appropriate background function.
    """

    terminate_background = BooleanProperty(True)
    """
    bool attribute used to determine if a generic background
    process should be terminated by the user. Value is True by default and 
    they only change on the start up of the appropriate background function.
    """

    terminate_load_files = BooleanProperty(True)
    """
    bool attribute used to determine if a file loading background
    process should be terminated by the user. Value is True by default and 
    they only change on the start up of the appropriate background function.
    """

    terminate_orto_search = BooleanProperty(True)
    """
    bool attribute used to determine if a orthology search background
    process should be terminated by the user. Value is True by default and 
    they only change on the start up of the appropriate background function.
    """

    terminate_group_export = BooleanProperty(True)
    """
    bool attribute used to determine if a ortholog export background
    process should be terminated by the user. Value is True by default and 
    they only change on the start up of the appropriate background function.
    """

    ################################
    #
    # CORE PROGRAM RELATED VARIABLES
    #
    ################################

    sqldb = StringProperty("")
    """String with path to sqlite database."""

    ortho_dir = StringProperty("")
    """String with path to OrthoMCL output directory."""

    orto_export_dir = StringProperty("")
    """String with path to directory for ortholog exports."""

    usearch_file = StringProperty("")
    """String with path to usearch executable."""
    usearch_db = StringProperty("goodProteins_db")
    """String with name of usearch's database file."""
    usearch_output = StringProperty("AllVsAll.out")
    """String with name of the usearch output file."""
    usearch_evalue = StringProperty("0.00001")
    """String with evalue for usearch execution."""

    # MCL/Groups attributes
    mcl_file = StringProperty("")
    """String with path to mcl executable."""
    ortholog_prefix = StringProperty("My_group")
    """String with prefix for ortholog groups."""
    group_prefix = StringProperty("group")
    """String with prefix for orthology search group files."""
    mcl_inflation = ListProperty(["3"])
    """List with inflation values, as strings, for the mcl execution."""

    # Protein quality filters
    protein_min_len = NumericProperty(10)
    """Integer with the minimum protein length for OrthoMCL."""
    protein_max_stop = NumericProperty(20)
    """Integer with the maximum percentage of stops for OrthoMCL."""

    # Orthology cluster filters
    orto_max_gene = NumericProperty(1)
    """
    Integer with maximum number of gene copies allowed for each
    ortholog group.
    """
    orto_min_sp = NumericProperty(3)
    """
    Integer with the minimum taxa representation allowed for each
    ortholog group..
    """

    # Attributes for exporting groups as protein/nucleotides
    protein_db = StringProperty("")
    """String with path to protein database file."""
    cds_db = ListProperty("")
    """List with paths to the CDS files corresponding to proteomes."""

    # Attribute containing the path to the proteome files
    proteome_files = []
    """List of paths to proteome files loaded into TriFusion."""
    active_proteome_files = ListProperty()
    """List of active proteome files."""

    # Attribute containing the orthology group files
    ortho_group_files = ListProperty()
    """List with paths to ortholog group files loaded into TriFusion."""
    ortho_groups = None
    """
    Reference to a :class:`~trifusion.ortho.OrthomclToolbox.MultiGroupLight`
    object.
    """
    active_group = None
    """
    Reference to a :class:`~trifusion.ortho.OrthomclToolbox.GroupLight`
    object, currently active.
    """
    active_group_name = None
    """
    String with name of the currently active
    :class:`~trifusion.ortho.OrthomclToolbox.GroupLight`.
    """

    alignment_list = None
    """
    Reference to a :class:`~trifusion.process.sequence.AlignmentList` object.
    """
    active_partitions = ListProperty()
    """
    List of currently active partitions.
    """

    trigger_taxa_update = False
    """
    bool attribute that determined whether taxa information needs to be
    updated (True) or not (False).
    """

    original_tx_inf = DictProperty()
    """
    Dictionary that maps a taxon string to the several informative
    properties for the entire loaded data set, shown when clicking the
    information button in the sidepanel.
    """
    active_tx_inf = DictProperty()
    """
    Dictionary that maps a taxon string to the several informative
    properties for the active data set, shown when clicking the
    information button in the sidepanel.
    """
    original_file_inf = DictProperty()
    """
    Dictionary that maps a filename string to the several informative
    properties for the entire loaded data set, shown when clicking the
    information button in the sidepanel. 
    """
    active_file_inf = DictProperty()
    """
    Dictionary that maps a filename string to the several informative 
    properties for the active data set, shown when clicking the
    information button in the sidepanel. 
    """
    #TODO: must check RAM usage of these dictionaries

    export_mode = None
    """
    Tuple attribute that determines the export mode in the
    sidepanel. The first element can be either 'taxa' or 'file'. The
    second element can be 'all' or 'select'.
    """

    sequence_types = StringProperty()
    """
    String attribute storing the sequence type of alignment files
    currently loaded.
    """

    # Attribute for taxa and file groups
    taxa_groups = DictProperty()
    """
    Dictionary storing custom taxa groups defined in TriFusion. Keys are
    group names and values are a list of taxon names.
    """
    file_groups = DictProperty()
    """
    Dictionary storing custom file groups defined in TriFusion. Keys are
    group names and values are a list of filenames.
    """
    dataset_file = None
    """
    String attribute with path to a data set file.
    """

    output_file = StringProperty("")
    """
    String with path to output file for Process execution.
    """
    output_dir = StringProperty("")
    """String with path to output directory for Process execution."""
    conversion_suffix = StringProperty("")
    """
    String with suffix to add for conversion operation of Process
    execution.
    """

    output_formats = ListProperty(["fasta"])
    """
    List of output currently selected output formats for Process execution.
    """

    ld_hat = BooleanProperty(False)
    """
    bool attribute that determines whether the Fasta output format shoud
    be compliant with LDhat (True) or a regular fasta (False).
    """
    # Determines whether the part.File associated with phylip format is
    # created
    create_partfile = BooleanProperty(True)
    """
    bool attribute that determines whether the partition file associated with
    the phylip format should be created (True) or not (False).
    """

    use_nexus_partitions = BooleanProperty(True)
    """
    bool attribute that determines whether the charset partitions block
    in nexus format should be written (True) or not (False).
    """

    use_nexus_models = BooleanProperty(True)
    """
    bool attribute that determines whether the substitution models
    associated to alignment partitions should be written to the nexus
    file (True) or not (False).
    """

    phylip_truncate_name = BooleanProperty(False)
    """
    bool attribute that determines whether taxa names should be truncated
    to 10 characters (True) or not (False) in a phylip output file.
    """

    ima2_options = ListProperty([None, None, None, None])
    """
    Four element list of additional options for IMa2 format:
    
        1. str - population file path
        2. str - Population tree string
        3. list - mutation model for each partition (one element if
    applies to all)
        4. list - inheritance scalaer for each partition (one element if
     applies to all)
    """

    missing_filter_settings = ListProperty([(True, 25, 50), (True, 0)])
    """
    Two element list attribute storing the missing data filter settings.
    The first element is a three element tuple:
    
        1. bool. Whether the missing data filter is active (True) or not
           (False).
        2. int. Percentage of gaps allowed.
        3. int, Percentage of missing data allowed.
    
    The second element is a two element tuple:
        
        1. bool. Whether the minimum taxa representation filter is active
           (True) or not (False).
        2. int. Percentage of the min taxa representation
    """

    taxa_filter_settings = ListProperty([None, None])
    """
    List attribute storing the taxa filter settings. The first element of
    the list should be the filter mode (either "Contain" or "Exclude")
    and the second element should be a string with the name of the taxa
    group (from the :attr:`~TriFusionApp.taxa_group` attribute).
    """

    codon_filter_settings = ListProperty([True, True, True])
    """
    List attribute storing the alignment filter settings. This will determine
    which codon positions will be written to the output (only for DNA
    sequences), so this will consist of a list containing 3 elements that
    correspond to each position. Positions will be saved or filtered
    depending on the boolean value of the list position. Ex. [True,
    True, True] will save all positions, whereas [True, True, False]
    will only save the first two positions.
    """

    variation_filter = ListProperty([None, None, None, None])
    """
    List attribute storing the alignment variation filter settings. The
    elements correspond to:
    
        1. Minimum variable sites.
        2. Maximum variable sites.
        3. Minimum informative sites.
        4. Maximum informative sites.
    """

    use_app_partitions = BooleanProperty(False)
    """
    bool atribute determining whether reverse concatenation will use a
    partition file (False) or the partitions defined in the app (True).
    """

    partitions_file = StringProperty("")
    """
    String with path to partitions file.
    """
    
    rev_infile = StringProperty("")
    """
    String with nameo the alignment file to be reverse concatenated.
    """

    zorro_suffix = StringProperty("")
    """String with the suffix of zorro files."""
    zorro_dir = StringProperty("")
    """String with the directory containing the zorro files."""

    hap_prefix = StringProperty("Hap")
    """String with the haplotype prefix for collapse operation."""

    def build(self):

        # Setting main window title
        self.title = "TriFusion - Streamline phylogenomics"

        # Set icon
        self.icon = join(self.cur_dir, "data", "backgrounds",
                         "trifusion-icon-512.png")

        # Set temporary directory path
        self.temp_dir = join(self.user_data_dir, "tmp")
        # Set directory for error records
        self.log_file = join(self.user_data_dir, "log", "error.out")
        # Set path for file hosting filechooser bookmarks
        self.bm_file = join(self.user_data_dir, "bookmarks")
        # Set path for file hosting the user projects
        self.projects_file = join(self.user_data_dir, "projects")
        # Set path for sequence data sqlite database
        self.sqldb = join(self.temp_dir, "sequencedb")

        logging.basicConfig(filename=self.log_file, level=logging.DEBUG,)

        # Setting available screens
        self.available_screens = ["main", "Orthology", "Process",
                                  "Statistics", "fc", "group_compare",
                                  "plot", "orto_plot"]
        self.screen_names = self.available_screens

        # Transforming screen names into complete paths to be loaded by kivy
        self.available_screens = [join(self.cur_dir, "data", "screens",
                                       "{}.kv".format(screen)) for screen in
                                  self.available_screens]

        # Store screen names specifically designed for plot display
        self.plot_screens = ["group_compare", "plot", "orto_plot",
                             "Statistics"]

        self.loaded_screens = dict((sc, None) for sc in
                                   self.available_screens)

        # First thing is go to main screen
        self.go_screen(0)

        # Set method for closing side panel when touching outside
        Window.bind(on_touch_up=lambda x, y: self.sidepanel_on_touch(y))

        # Listen to keybindings
        Window.bind(on_key_down=self._on_keyboard_events)

        Window.bind(on_key_up=self._release_events)

        # Execute cleaning function when exiting app
        Window.bind(on_request_close=lambda x:
                    self.check_action("Are you sure you want to quit?",
                                      self._exit_clean, popup_level=3))

        Window.bind(on_stop=lambda x: self._exit_clean())

        Window.bind(on_motion=self.mouse_zoom)

        Window.bind(on_dropfile=self.load_files_dragndrop)

        Window.clearcolor = tm.c_main_background

        # Orthology widgets
        self.ortho_search_options = OrthologySearchGrid()
        # Check for orthology programs
        self._check_orto_programs()

        # Process widgets
        # Creating GridLayout instance for general options of Process
        self.process_grid_wgt = ProcessGeneral()
        # Create GridLayout instance for additional options of Process.
        self.process_options = AdditionalProcessContents()
        # Create Summary stats box instance for Statistics screen
        self.stats_summary = StatsSummary()
        self.stats_table = GeneTable()

        # Initialize operation queue treeview in side panel
        self.operation_queue_init()

        # Initialize projects
        self.projects_init()

        self.fancy_bt = ObjectProperty(FancyButton())

        # Set schedule for mouse over events on side panel
        Clock.schedule_interval(lambda x: self._on_mouseover_tabs(), .1)

        # Sets schedule for removable media automatic detection
        Clock.schedule_interval(lambda x: self._check_removable_media(), 2)

        # Set path to sqlite database
        self.sqldb = join(self.temp_dir, "trifusion.sql3")

        self._start_clean()

        # Now that the sqlite database path is defined, initialize the
        # AlignmentList object
        self.alignment_list = AlignmentList([], sql_db=self.sqldb)

        # If arguments were provided at the command line, load data into
        # the app assuming they are alignments
        if sys.argv[1:]:
            Clock.schedule_once(
                lambda dt: self.load_files_startup(sys.argv[1:]), .1)

    def mouse_zoom(self, *vals):
        """
        :param vals: touch event list
        """

        # Only perform any actions in plot screens and pop is not active
        if self.screen.name in self.plot_screens and \
                self._popup not in self.root_window.children and \
                self._subpopup not in self.root_window.children:

            # Ignore if in statistics screen and sidepanel expanded
            if self.screen.name == "Statistics" and \
                        self.screen.ids.stats_panel.width == 410:
                return

            # Add exceptions here
            try:
                if isinstance(self.screen.ids.plot_content.children[
                              0].children[0], StatsSummary) or \
                        isinstance(self.screen.ids.plot_content.children[
                                   0].children[0], GeneTable):
                    return
            except (IndexError, ValueError):
                pass

            motion = vals[2]

            # Check if motion is mouse scroll
            if motion.is_mouse_scrolling:

                # Only perform an action at the begining of the motion
                if vals[1] == "begin" and self.is_control_pressed:

                    if motion.button == "scrollup":

                        self.screen.ids.plot_content.scale -= .1

                    elif motion.button == "scrolldown":

                        self.screen.ids.plot_content.scale += .1

    def _start_clean(self):
        """
        In the event of unexpected exits, clean the tmp directory on app
        start
        :return:
        """

        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

        os.makedirs(self.temp_dir)

    def _check_exec(self, path, program):
        """
        Checks the execution of the MCL software from a given path. If the MCL
        is correct, returns True
        """

        msg = {"mcl": b"mcl",
               "usearch": "usearch"}

        # If the OS is windows and the program is usearch, first check if the
        # necessary DLL vcomp100.dll is present in the app's directory
        if sys.platform in ["win32", "cygwin"] and program == "usearch":
            flist = [x.lower() for x in os.listdir(self.user_data_dir)]
            if "vcomp100.dll" not in flist:
                self.dismiss_subpopup()
                return self.dialog_warning(
                    "Missing a DLL file",
                    "You are providing the USEARCH executable but the"
                    " vcomp100.dll file is missing. Please download it from"
                    " the link below and place the DLL in"
                    " C:\Users\<USER>\AppData\Roaming\\trifusion. Then"
                    " try again.\n\n"
                    "[ref=http://drive5.com/usearch/vcomp100.dll][b]"
                    "[color=#37abc8ff]http://drive5.com/usearch/vcomp100.dll"
                    "[/b][/color][/ref]")

        # First, try to execute the file directly from the path
        try:
            res, _ = subprocess.Popen([path, "--version"],
                                      stdout=subprocess.PIPE).communicate()

            return res.startswith(msg[program])

        except OSError as e:
            # This exception may occur because the file is not executable
            if e.errno == 13:

                try:
                    # Make the file executable
                    st = os.stat(path)
                    os.chmod(path, st.st_mode | stat.S_IEXEC)

                    res, _ = subprocess.Popen(
                        [path, "--version"],
                        stdout=subprocess.PIPE).communicate()

                    return res.startswith(msg[program])

                except OSError:
                    # The provided file is not executable and TriFusion does
                    # not have permissions to change that.
                    return False
            else:
                # The first OSError exception could not be handled. Fail the
                # executable check
                return False

    def _check_orto_programs(self):
        """
        Checks if the USEARCH and MCL programs are installed and reachable.
        1. Check if executables are in the app_dir
        2. Check if executables are callable from subprocess
        3. If both checks fail, modify ortho_search_grid to highlight the
        problem.
        """

        dll_path = None

        # The mcl executable requires the .exe extension in windows to
        # work properly
        if sys.platform in ["win32", "cygwin"]:
            mcl_string = "mcl.exe"
        else:
            mcl_string = "mcl"

        # Check if MCL_FILE has been set. This happens when running TriFusion
        # from the executable binary.
        if MCL_FILE:
            self.mcl_file = MCL_FILE
        # Check MCL in app_dir
        elif os.path.exists(join(self.user_data_dir, mcl_string)):
            self.mcl_file = join(self.user_data_dir, mcl_string)
        # Check MCL in package directory, if running from source
        elif os.path.exists(join(os.getcwd(), "data", "resources", "mcl")):
            mcl_dir = join(os.getcwd(), "data", "resources", "mcl")
            # For Windows
            if sys.platform in ["win32", "cygwin"]:
                # For Windows 64bit
                if platform.architecture()[0] == "64bit":
                    mcl_path = join(mcl_dir, "windows", "64bit", "mcl64.exe")
                    dll_path = join(mcl_dir, "windows", "64bit",
                                    "cygwin1.dll")
                # For Windows 32bit
                else:
                    mcl_path = join(mcl_dir, "windows", "32bit", "mcl32.exe")
                    dll_path = join(mcl_dir, "windows", "32bit",
                                    "cygwin1.dll")

            # For MacOS
            elif sys.platform == "darwin":
                mcl_path = join(mcl_dir, "MacOS", mcl_string)

            # For Linux
            else:
                mcl_path = join(mcl_dir, "linux", mcl_string)

            # Copy mcl executable to trifusion dir
            shutil.copyfile(mcl_path, join(self.user_data_dir,
                                           mcl_string))
            mcl_file = join(self.user_data_dir, mcl_string)
            # Copy dll necessary in windows version
            if dll_path:
                shutil.copyfile(dll_path, join(self.user_data_dir,
                                               "cygwin1.dll"))
            # Make mcl executable
            st = os.stat(mcl_file)
            os.chmod(mcl_file, st.st_mode | stat.S_IEXEC)
            # Test mcl executable
            if self._check_exec(mcl_file, "mcl"):
                self.mcl_file = mcl_file
            else:
                # Remove bad executable from app dir
                os.remove(join(self.user_data_dir, mcl_string))
                if dll_path:
                    os.remove(join(self.user_data_dir, "cygwin1.dll"))

        else:
            # If not in app_dir check is its reachable system-wide by
            # subprocess
            try:
                subprocess.call(["mcl"])
                self.mcl_file = mcl_string
            except OSError:
                pass

        # If the mcl file could not be found through any of the previous
        #  methods, modify the attributes of the orto_search_grid.
        if not self.mcl_file:
            self.ortho_search_options.ids.mcl_check.background_src = \
                "data/backgrounds/red_noise.png"
            self.ortho_search_options.ids.mcl_check_ico.icon_src = \
                "data/backgrounds/warning_icon.png"
            self.ortho_search_options.ids.mcl_check_ico.icon_size = \
                (20, 20)
            self.ortho_search_options.ids.mcl_check_lbl.text = \
                "MCL is not installed or reachable"
            fix_bt = FixButton()
            fix_bt.bind(on_release=lambda x: self.dialog_mcl_fix())
            self.ortho_search_options.ids.mcl_fix_box.add_widget(fix_bt)

        # Check USEARCH DLL in app_dir for windows only
        if sys.platform in ["win32", "cygwin"] and \
                not os.path.exists(join(self.user_data_dir, "vcomp100.dll")):
            # Copy necessary DLL to app dir if it doesn't exist
            try:
                if platform.architecture()[0] == "64bit":
                    shutil.copyfile(join(os.getcwd(), "data", "resources",
                                         "usearch", "64bit", "vcomp100.dll"),
                                    join(self.user_data_dir, "vcomp100.dll"))
                else:
                    shutil.copyfile(join(os.getcwd(), "data", "resources",
                                         "usearch", "32bit", "vcomp100.dll"),
                                    join(self.user_data_dir, "vcomp100.dll"))
            except IOError:
                pass

        # CHeck USEARCH in app_dir
        if os.path.exists(join(self.user_data_dir, "usearch")):
            self.usearch_file = join(self.user_data_dir, "usearch")
        else:
            # If not in app_dir check is its reachable system-wide by
            # subprocess
            try:
                subprocess.call(["usearch"])
                self.usearch_file = "usearch"
            # MCL software not reachable. Modify attributes and
            # orto_search_grid
            except OSError:
                self.ortho_search_options.ids.usearch_check.background_src = \
                    "data/backgrounds/red_noise.png"
                self.ortho_search_options.ids.usearch_check_ico.icon_src = \
                    "data/backgrounds/warning_icon.png"
                self.ortho_search_options.ids.usearch_check_ico.icon_size = \
                    (20, 20)
                self.ortho_search_options.ids.usearch_check_lbl.text = \
                    "USEARCH is not installed or reachable"
                fix_bt = FixButton()
                fix_bt.bind(on_release=lambda x: self.dialog_usearch_fix())
                self.ortho_search_options.ids.usearch_fix.add_widget(fix_bt)

    def dialog_mcl_fix(self):
        """
        Opens a dialog with information for fixing the missing MCL executable
        """

        content = DialogMCLFix(cancel=self.dismiss_popup)

        self.show_popup(title="MCL troubleshooting...", content=content,
                        size=(400, 330), close_bt=True)

    def dialog_usearch_fix(self):
        """
        Opens a dialog with information for fixing the missing USEARCH
        executable
        """

        content = DialogUSEARCHFix(cancel=self.dismiss_popup)

        self.show_popup(title="USEARCH troubleshooting...", content=content,
                        size=(400, 330), close_bt=True)

    def load_files_startup(self, file_list):
        """
        Gives support for loading data into the app when launching the
        app via command line and providing the files as arguments. Issues
        the same dialog as the load_files_dragndrop to let the user
        choose whether the input files are alignments or proteomes
        :param file_list: list, with file paths as strings or unicode
        """

        content = InputType(cancel=self.dismiss_popup)
        content.files = file_list

        self.show_popup(title="", content=content,
                        size=(350, 160),
                        separator_color=(0, 0, 0, 0),
                        close_bt=True)

    def load_files_dragndrop(self, *args):
        """
        This function gives functionality to the drag and drop feature
        that automatically loads files when droped into the app window.
        Note that this method will issue different data loding methods
        depending on the active screen. Proteome alignments will be opened
        in the Orthology screen, while sequence alignment will be opened
        in the Process and Statistics screens.
        :param args: list, first element is the SDL2 object, second
        element is the file path
        """

        def drag_check(dt):
            """
            In order to support drag and drop of multiple files, a periodic
            check is made to assess if the self.drag_files object is
            still being populated, or if it reached a static state. Only
            when it reaches a static state, will the app issue the
            corresponding methods to load all files into the app.
            :param dt:
            """

            # If the drag_files attribute continues to be populated,
            # its lenght will be higher than the current counter. In such
            #  case, update the counter
            if len(self.drag_files) > self.drag_c:
                self.drag_c = len(self.drag_files)
                # In case only a single file was dragged and dropped,
                # this additional scheduled event will allow the single
                # file to the loaded. Otherwise, it will do nothing.
                Clock.schedule_once(drag_check, .1)
            # In this case, the drag_files attribute has reached a static
            #  size and the loading methods can be issued
            elif len(self.drag_files) == self.drag_c:
                # Issue methods only when drag_files is populated
                if self.drag_files:

                    content = InputType(cancel=self.dismiss_popup)
                    content.files = self.drag_files

                    self.show_popup(title="", content=content,
                                    size=(350, 160),
                                    separator_color=(0, 0, 0, 0),
                                    border_color=tm.c_popup_border,
                                    close_bt=True)

                    # Load proteomes
                    # if self.screen.name == "Orthology":
                    #     self.load_proteomes(self.drag_files)
                    # if self.screen.name in ["Process", "Statistics",
                    #                         "main"]:
                    #     self.load_files_subproc(self.drag_files)

                    # Reset drag counter and drag_files attributes
                    self.drag_c = 0
                    self.drag_files = []

        # This functionality should only be triggered when in one of the
        # main screens

        self.dismiss_all_popups()

        self.drag_files.append(unicode(args[1], "utf8"))

        Clock.schedule_once(drag_check, .1)

    def _exit_clean(self):
        """
        This method is issued when the application is closed and performs
        any necessary clean up operations
        """

        self.terminate_stats = True
        self.terminate_load_files = True
        self.terminate_background = True

        self.run_in_background(remove_tmp, self.stop,
                               [self.temp_dir, self.alignment_list.con],
                               no_arg2=True, cancel=False,
                               msg="Cleaning temporary files and "
                                   "exiting...")

    def _update_path(self, path):
        """
        This method updates the filechooser path when clicking on the path
        label

        :param path: string, with destination path
        """

        # If in filechooser, do this
        try:
            self.screen.ids.path_bx.children[0].text = path
            self.screen.ids.icon_view_tab.path = path
        # When in popup, do thissd_filechooser
        except AttributeError:
            try:
                self._popup.content.ids.path_bx.children[0].text = path
                self._popup.content.ids.sd_filechooser.path = path
            except AttributeError:
                self._subpopup.content.ids.path_bx.children[0].text = path
                self._subpopup.content.ids.sd_filechooser.path = path

    def _release_events(self, *vals):
        """
        Method that releases keyboard events that are triggered in
        _on_keyboard_events when the key is released
        :param vals: input list from on_key_up
        """

        key_code = vals[1]

        if key_code == 305:
            self.is_control_pressed = False

        if key_code == 304:
            self.is_shift_pressed = False

    def _on_keyboard_events(self, *vals):
        """
        Methods that listens to keyboard input and triggers events or
        changes properties when acceptable keyboard shortcuts are entered
        :param vals: input list from on_key_down function
        """

        # Get modifier (usually ctrl or shift)
        # TODO: The modifier in MacOS is different. Must check on this.
        modifier = "".join(vals[-1])
        key_code = vals[1]

        # Switch between Ctrl for Windows/Linux and Cmd for MacOS as the
        # preferred modifier
        if sys.platform == "darwin":
            mod_key = "meta"
        else:
            mod_key = "ctrl"

        if key_code == 305:
            self.is_control_pressed = True

        if key_code == 304:
            self.is_shift_pressed = True

        # ==================================================================
        # Popup keybindings
        # ==================================================================

        def popup_keys(backn, backd, bt1, bt2, bt3=None):
            """
            Wrapper function that provides functionality to arrow keys for
            navigating through selection buttons in popups
            :param backn: string for background normal
            :param backd: string for background down
            :param bt1: Button widget one (usually for Ok buttons)
            :param bt2: Button widget two (usually for Cancel buttons)
            :param bt3: Optional button widget three.
            """

            if not self.arrow_block:
                # This will deal with cases with only two buttons to cycle
                if not bt3:
                    # if left arrow key
                    if key_code == 276:
                        bt1.background_normal = backn
                        bt2.background_normal = backd
                    # if right arrow key
                    if key_code == 275:
                        bt1.background_normal = backd
                        bt2.background_normal = backn
                    # if enter key. Dispatch the events of the focused
                    # button
                    if key_code == 13:
                        if bt1.background_normal == backn:
                            bt1.dispatch("on_release")
                        else:
                            bt2.dispatch("on_release")

                # This will cycle through three buttons
                else:
                    bt_list = [bt1, bt2, bt3]
                    idx = [i.background_normal for i
                           in bt_list].index(backn)
                    if key_code == 276 and idx > 0:
                        idx -= 1
                    if key_code == 275 and idx < 2:
                        idx += 1

                    for bt in bt_list:
                        if bt_list.index(bt) == idx:
                            bt_list[idx].background_normal = backn
                        else:
                            bt.background_normal = backd

                    if key_code == 13:
                        bt_on = [i for i in bt_list if
                                 i.background_normal == backn][0]
                        bt_on.dispatch("on_release")

        # ==================================================================
        # Popup keybindings
        # ==================================================================

        if self._exit_popup in self.root_window.children:
            if "check_ok" in self._exit_popup.content.ids:
                bn = "data/backgrounds/check_ok.png"
                bd = "data/backgrounds/check_cancel.png"
                ok_bt = self._exit_popup.content.ids["check_ok"]
                cancel_bt = self._exit_popup.content.ids["check_cancel"]
                popup_keys(bn, bd, ok_bt, cancel_bt)

        if self._exit_popup in self.root_window.children:
            return

        if self._popup in self.root_window.children:
            if "find_bx" in self._popup.content.ids:
                if modifier == mod_key and key_code == 102:
                    if self._popup.content.ids.find_bx.height == 30:
                        Animation(height=0, d=.3, t="out_quart").start(
                            self._popup.content.ids.find_bx)
                        self._popup.content.ids.sd_filechooser.f = "*"
                        self._popup.content.ids.text_filter.text = ""
                        self._popup.content.ids.text_filter.focus = False
                    else:
                        Animation(height=30, d=.3, t="out_quart").start(
                            self._popup.content.ids.find_bx)
                        self._popup.content.ids.text_filter.focus = True

        if self._subpopup in self.root_window.children:
            if "find_bx" in self._subpopup.content.ids:
                if modifier == mod_key and key_code == 102:
                    if self._subpopup.content.ids.find_bx.height == 30:
                        Animation(height=0, d=.3,
                                  t="out_quart").start(
                            self._subpopup.content.ids.find_bx)
                        self._subpopup.content.ids.sd_filechooser.f = "*"
                        self._subpopup.content.ids.text_filter.text = ""
                        self._subpopup.content.ids.text_filter.focus = False
                    else:
                        Animation(height=30, d=.3,
                                  t="out_quart").start(
                            self._subpopup.content.ids.find_bx)
                        self._subpopup.content.ids.text_filter.focus = True

        # Toggle path editor
        if self._popup in self.root_window.children:
            if "path_toggle" in self._popup.content.ids:
                if modifier == mod_key and key_code == 108:
                    st = self._popup.content.ids.path_toggle.state
                    self._popup.content.ids.path_toggle.state = \
                        "down" if st != "down" else "normal"
                    self._popup.content.ids.path_toggle.dispatch("on_release")

        if self._popup in self.root_window.children:
            if "text_filter" in self._popup.content.ids:
                # Ctrl + f toggles focus on find text input field
                if modifier == mod_key and key_code == 102:
                    self._popup.content.ids.text_filter.focus = True

        if self._subpopup in self.root_window.children:
            if "ok_bt" in self._subpopup.content.ids:
                bn = "data/backgrounds/bt_process.png"
                bd = "data/backgrounds/bt_process_off.png"
                ok_bt = self._subpopup.content.ids["ok_bt"]
                cancel_bt = self._subpopup.content.ids["cancel_bt"]
                popup_keys(bn, bd, ok_bt, cancel_bt)

        elif self._popup in self.root_window.children:
            if "check_ok" in self._popup.content.ids:
                bn = "data/backgrounds/check_ok.png"
                bd = "data/backgrounds/check_cancel.png"
                ok_bt = self._popup.content.ids["check_ok"]
                cancel_bt = self._popup.content.ids["check_cancel"]
                popup_keys(bn, bd, ok_bt, cancel_bt)
            # In this case there are three buttons to cicle
            elif "apply_bt" in self._popup.content.ids:
                bn = "data/backgrounds/bt_process.png"
                bd = "data/backgrounds/bt_process_off.png"
                ok_bt = self._popup.content.ids["ok_bt"]
                cancel_bt = self._popup.content.ids["cancel_bt"]
                apply_bt = self._popup.content.ids["apply_bt"]
                popup_keys(bn, bd, ok_bt, apply_bt, cancel_bt)
            # Only two buttons to cicle
            elif "ok_bt" in self._popup.content.ids:
                bn = "data/backgrounds/bt_process.png"
                bd = "data/backgrounds/bt_process_off.png"
                ok_bt = self._popup.content.ids["ok_bt"]
                cancel_bt = self._popup.content.ids["cancel_bt"]
                popup_keys(bn, bd, ok_bt, cancel_bt)

            elif "close_bt" in self._popup.content.ids:
                if key_code == 13:
                    self._popup.content.ids.close_bt.dispatch("on_release")

        # ======================================================================
        # Filechooser keybindings
        # ======================================================================

        if self.screen.name == "fc":
            # Keybinding ctrl+f that brings focus to the "Find" field in the
            # Filechooser screen
            if modifier == mod_key and key_code == 102:
                self.screen.ids.text_filter.focus = True

            # Keybinding ctrl+backspace to clear selection
            if modifier == mod_key and key_code == 8:
                self.screen.ids.clear_s.dispatch("on_release")

            # Add bookmarks with ctrl+d
            if modifier == mod_key and key_code == 100:
                self.screen.ids.add_bk_bt.dispatch("on_release")

            # Toggle manual path writing with ctrl+l
            if modifier == mod_key and key_code == 108:
                if self.screen.ids.path_toggle.state == "down":
                    self.screen.ids.path_toggle.state = "normal"
                else:
                    self.screen.ids.path_toggle.state = "down"
                self.screen.ids.path_toggle.dispatch("on_release")

            # Select all files with ctrl+a
            if modifier == mod_key and key_code == 97:
                self.screen.ids.icon_view_tab.selection = \
                    [x for x in self.screen.ids.icon_view_tab.files if not
                    os.path.isdir(x)]

            if self._popup not in self.root_window.children:
                # Use arrow keys and enter to navigate through open/cancel
                # buttons and selecting them
                bn = "data/backgrounds/bt_process.png"
                bd = "data/backgrounds/bt_process_off.png"
                open_close_bt = self.screen.ids.open_close_bt
                open_bt = self.screen.ids.open_bt
                cancel_bt = self.screen.ids.cancel_bt
                popup_keys(bn, bd, open_bt, open_close_bt, cancel_bt)

        # ==================================================================
        # General keybindings
        # ==================================================================

        if not [x for x in self.root_window.children
                if isinstance(x, CustomPopup)]:

            # Keybinding ctrl+o that opens the Filechooser screen
            if modifier == mod_key and key_code == 111:
                self.go_screen(self.screen_names.index("fc"))

            # Changing main screens between Orthology, Process and
            # Statistics
            if modifier == mod_key and key_code == 49:
                self.root.ids.h_ortho.dispatch("on_release")
                self.root.ids.h_ortho.state = "down"
            if modifier == mod_key and key_code == 50:
                self.root.ids.h_process.dispatch("on_release")
                self.root.ids.h_process.state = "down"
            if modifier == mod_key and key_code == 51:
                self.root.ids.h_stat.dispatch("on_release")
                self.root.ids.h_stat.state = "down"

            # Toggle side panel (slash)
            if key_code == 92:
                self.root.ids.ap.dispatch("on_release")

        # ==================================================================
        # Text input autocompletion
        # ==================================================================

        # Use tab for auto completion when textinput is focused
        path_wgt = None
        if key_code == 9:
            if "path_bx" in self.screen.ids:
                path_wgt = self.screen.ids.path_bx.children[0]
            elif self._popup:
                if "path_bx" in self._popup.content.ids:
                    path_wgt = self._popup.content.ids.path_bx.children[0]

            if isinstance(path_wgt, TextInput):
                path = path_wgt.text
                s = self._auto_completion(path)
                path_wgt.text = s

    def _common_path(self, file_list=None):
        """
        Returns a string with the longest common path contained in the
        self.file_list dir. This is used by FileChoosers to open the nearest
        directory to their input files
        """

        if file_list:
            fl = file_list
        else:
            fl = self.file_list

        # Return None if there are no input files
        if not fl:
            return

        # Get common path
        common_path = ""
        for char in zip(*fl):
            if len(set(char)) == 1:
                common_path += "".join(set(char))
            else:
                break

        # Get nearest directory from common path
        while not os.path.isdir(common_path):
            common_path = sep.join(common_path.split(sep)[:-1])

        return common_path

    def _auto_completion(self, path):
        """
        Method used for providing auto completion for text input widgets
        navigating the os file system
        """

        # Check if current path exists. If it does, there is no path to auto
        # complete and no action is required.
        if os.path.exists(path) or path == "":
            return path
        # If the path does not exist, get the nearest parent directory and
        # use the final string of the path to text for auto completion
        else:
            original_path = path
            s = path.split(sep)[-1]
            path = sep.join(path.split(sep)[:-1])

        # Get list of contents that may
        dirlist = [x for x in os.listdir(path) if x.startswith(s) and
                   os.path.isdir(join(path, x))]

        # If there is only one match in dirlist, return that match
        if len(dirlist) == 1:
            return join(path, dirlist[0])

        # If there are multiple matches in dirlist, return the longest
        # common substring
        elif len(dirlist) > 1:
            return join(path, os.path.commonprefix(dirlist))

        else:
            self.dialog_floatcheck("Path does not exist",
                                   t="warning")
            return original_path

    @staticmethod
    def _determine_collision(wgt, mp):
        """
        Provided a widget, this function determines whether the mouse is
        colliding with its window coordinates
        :param wgt: ToggleButton widget inside a relative layout
        :return: Boolean. True is mouse position collides with wgt
        """

        # Retrieving window widget coordinates
        window_pos = wgt.to_window(wgt.pos[0], wgt.pos[1])
        # Creating dummy widget to determine collision
        dummy_wgt = Widget(pos=window_pos, size_hint=(None, None),
                           size=wgt.size)

        return dummy_wgt.collide_point(mp[0], mp[1])

    def _check_removable_media(self):
        """
        Checks, when a Filechooser dialog is active, for the insertion or
        removal of external storage devices, and updates the bookmark
        buttons accordingly
        """

        def detect_changes():
            """
            Detects changes between the app attribute self.removable_media
            and the removable media detected in the system. If
            differences are found, returns the operation (remove ou add)
            and a list of paths
            :return: Tuple, 1? element is the operation (add/remove); 2?
            element is a list of paths with which to perform the operation
            """

            # For linux
            if sys.platform in ["linux", "linux2"]:

                # Get current removable media
                i = subprocess.Popen(["mount | grep -e /media -e /mnt | "
                    "awk '{print $3}'"], shell=True, stdout=subprocess.PIPE)
                removable_media = i.communicate()[0].split("\n")[:-1]

            # For MacOS
            elif sys.platform == "darwin":

                removable_media = [join("/Volumes", x) for x in
                                   os.listdir("/Volumes/")]

            # For Windows
            elif sys.platform in ["win32", "cygwin"]:

                removable_media = re.findall(
                    r"[A-Z]+:.*$", os.popen("mountvol /").read(), re.MULTILINE)

                # Prevent no disk on drive popup error in windows
                SEM_FAILCRITICALERRORS = 1
                SEM_NOOPENFILEERRORBOX = 0x8000
                SEM_FAIL = SEM_NOOPENFILEERRORBOX | SEM_FAILCRITICALERRORS

                kernel32 = ctypes.WinDLL('kernel32')
                oldmode = ctypes.c_uint()
                kernel32.SetThreadErrorMode(SEM_FAIL, ctypes.byref(oldmode))

            else:
                return

            cur, act = set(removable_media), set(self.removable_media)

            if cur == act:
                return False

            elif cur.difference(act):
                return "add", [p for p in cur.difference(act)]

            elif act.difference(cur):
                return "remove", [p for p in act.difference(cur)]

        def update_media(operation, paths, wgt, fc_wgt):
            """
            Updates media, according to the operation and list of paths
            :param operation: string, Either "add" or "remove"
            :param paths: list, containing the paths to the added or removed
            :param wgt: Widget, preferentially a gridlayout where the
            bookmark buttons will be added/removed
            :param fc_wgt: FileChooser widget in which bookmark operations
            will be performed
            """

            if operation == "add":

                for p in paths:

                    # First, check if the provided path is a directory and
                    # is reachable. If not, we do not add it as a bookmark
                    if not os.path.isdir(p):
                        return

                    # For windows, get the name of the drive manually
                    if sys.platform in ["win32", "cygwin"]:
                        name = os.path.splitdrive(p)[0]
                    else:
                        name = None

                    self.add_bookmark_bt(p, wgt, fc_wgt, rm_bt=False,
                                         name=name)
                    self.removable_media.append(p)

            elif operation == "remove":

                for bt in wgt.children:
                    if bt.id in paths:
                        self.removable_media.remove(bt.id)
                        wgt.remove_widget(bt)

        if self.screen.name == "fc":

            res = detect_changes()

            if res:
                update_media(res[0], res[1], self.screen.ids.sv_mycomp,
                             self.screen.ids.icon_view_tab)

        elif [x for x in self.root_window.children if
                isinstance(x, CustomPopup)]:
            # Get popup content
            wgt = [x for x in self.root_window.children if
                   isinstance(x, CustomPopup)][0].content

            if isinstance(wgt, SaveDialog) or \
                    isinstance(wgt, LoadMultipleDialog):

                res = detect_changes()

                if res:
                    update_media(res[0], res[1], wgt.ids.sv_mycomp,
                                 wgt.ids.sd_filechooser)

    def _on_mouseover_tabs(self):
        """
        Provides mouse over events throughout the app. In order to reduce
        computations in each cycle, the mouse over events should only be
        tested when the appropriate screen/widget is visible
        """

        try:
            self.mouse_position = self.root_window.mouse_pos
        except ReferenceError:
            pass

        # Get mouse position coordinates
        mp = self.root_window.mouse_pos
        # Set collision attribute
        collision = False
        # Set side button list
        sidebt_list = [x for x in self.root.ids.side_bt.children if
                       isinstance(x, ToggleButton)]

        def show_label(mouse, wgt):
            """
            Use this function with a Clock schedule to delay the
            introduction of the label widget. Otherwise, it could become
            cumbersome to have an label appearing right after the mouse
            colliding with the button
            :param mouse: mouse position tuple
            :param wgt: Layout widget, containing a descriptive label
            """

            # Checking if the current mouse position is the same as the
            #  mouse position when the mouse over was triggered. This
            # ensures that if the user changes the mouse position while
            # this event is scheduled, the label will not be added to the
            #  root_window but the lock in self.mouse_over_ready is removed
            if self.root_window.mouse_pos == mouse:
                # Add widget to root layout
                self.root_window.add_widget(wgt)
                # Update old label widget
                self.old_mouse_over = wgt
                # Update old label text
                self.previous_mouse_over = wgt.text

            # Unlocking mouse over
            self.mouse_over_ready = True

        def create_sidebt_wgt(text, p, s):
            """
            Creates the label for the sidebt mouseover
            :param text: string. Text for the label
            :param p: tuple. position of the wgt
            :param s: tuple. size of the wgt
            """

            clear_mouse_overs()
            side_bt = SideLabel(text=text, pos=p, size_hint=(None, None),
                                size=s, bold=True, border=(0, 0, 0, 0))

            return side_bt

        def create_fancy_label(text, wgt, lbl_height=30, adjust_pos=False,
                               c=(0.216, 0.67, 0.784, 1), wgt_pos=None,
                               wgt_size=None, orientation="horizontal",
                               line_c=(1, 1, 1, 1)):
            """
            Creates a fancy label akin to the mouse over labels in github.
            This method is quite versatile as it is able to calculate the
            available space for the label and determine its orientation
            accordingly. It also starts by removing any other potential
            fancy labels.

            :param text: string, text that will be displayed in the label
            :param wgt: widget that triggered the mouse over
            :param lbl_height: int, height of the label
            :param adjust_pos: Boolean, If True the position of wgt will
            be adjusted to window coordenates; If False, use the origial
            position of the wgt. This is only relevant when wgt_pos=None
            :param c: tuple. Label color
            :param wgt_pos: tuple/list, If not None, this will be the
            position used to create the label. Superseeds wgt.pos
            :param wgt_size: tuple/list, If not None this will be the size
            of the wgt. Superseed wgt.size
            :param orientation: string, whether the label will be display in
            the same horizontal or vertical plane of the widget
            :param line_c: tuple/list, color of the border line and arrow
            """

            # If the current mouse position is no longer the same when the
            # label was issued, then do nothing
            if self.root_window.mouse_pos != mp:
                # Unlocking mouse over
                self.mouse_over_ready = True
                # Cleans any possible mouse overs
                clear_mouse_overs()
                return

            # Cleans any possible mouse overs
            clear_mouse_overs()

            # Determines which coordinates to use
            if wgt_pos:
                wgt_pos = wgt_pos
            else:
                if adjust_pos:
                    wgt_pos = wgt.to_window(wgt.pos[0], wgt.pos[1])
                else:
                    wgt_pos = wgt.pos

            # Determines size attribute
            if wgt_size:
                wgt_size = wgt_size
            else:
                wgt_size = wgt.size

            # Create label
            self.fancy_bt = FancyButton(text=text, height=lbl_height,
                                        id=text, background_color=c)
            # Update label texture size, so that we can evaluate the
            # available space for the label
            self.fancy_bt.texture_update()
            # Set border line color
            self.fancy_bt.line_clr = line_c

            # Determine if there is enough space for  the label to be
            # properly shown. If not, truncante the label width to 70% of
            #  window width
            if (wgt_pos[0] + wgt_size[0] + 5 +
                    self.fancy_bt.texture_size[0] + 50 >
                    self.root.width and wgt_pos[0] - wgt_size[0] - 5 -
                    self.fancy_bt.texture_size[0] - 50 < 0):
                self.fancy_bt.width = self.root.width * .7
                self.fancy_bt.text_size = self.fancy_bt.size
                self.fancy_bt.halign = "center"
                self.fancy_bt.valign = "middle"
                self.fancy_bt.texture_update()
                self.fancy_bt.height = self.fancy_bt.texture_size[1] + 16
            else:
                # Determine label size. Add horizontal margin space (10)
                self.fancy_bt.size = (self.fancy_bt.texture_size[0] +
                                      10, lbl_height)

            # Create horizontal label
            if orientation == "horizontal":
                # Determine if the label has space to the right. If not,
                # flip the orientation
                if (wgt_pos[0] + wgt_size[0] + 5 + self.fancy_bt.width <
                        self.root.width):

                    # Determine position of arrow widget
                    point_pos = wgt_pos[0] + wgt_size[0] + 5, wgt_pos[1] + \
                        wgt_size[1] / 2 - 6

                    # Determine label position
                    self.fancy_bt.pos = point_pos[0] + 7, wgt_pos[1] + \
                        wgt_size[1] / 2 - self.fancy_bt.height / 2

                    # Create arrow widget with left arrow
                    point_wgt = FancyMarker(background_normal=join("data",
                                            "backgrounds",
                                            "box_arrow_right.png"),
                                            pos=point_pos, size=(7, 12),
                                            background_color=line_c)
                # In case this else code is executed, it means there is no
                #  space for the label to be shown in the right,
                # so it will show to the left
                else:
                    # Determine position of arrow widget
                    point_pos = wgt_pos[0] - 10, wgt_pos[1] + \
                        wgt_size[1] / 2 - 6

                    # Determine label position
                    self.fancy_bt.pos = point_pos[0] - self.fancy_bt.width,\
                        wgt_pos[1] + wgt_size[1] / 2 - \
                        self.fancy_bt.height / 2

                    # Create arrow widget with left arrow
                    point_wgt = FancyMarker(background_normal=join("data",
                                            "backgrounds",
                                            "box_arrow_left.png"),
                                            pos=point_pos, size=(7, 12),
                                            id=text,
                                            background_color=line_c)
            # For vertical orientation
            else:
                # For now, show always on top
                # Determine position of arrow
                point_pos = [wgt_pos[0] + (wgt_size[0] / 2),
                             wgt_pos[1] + wgt_size[1] + 5]

                # Determine position of label
                self.fancy_bt.pos = [point_pos[0] -
                    (self.fancy_bt.width / 2) + 6, point_pos[1] + 7]

                # Create arrow widget with down arrow
                point_wgt = FancyMarker(background_normal=join("data",
                                        "backgrounds",
                                        "box_arrow_down.png"),
                                        pos=point_pos, size=(12, 7),
                                        id=text,
                                        background_color=line_c)

            for w in [self.fancy_bt, point_wgt]:
                self.root_window.add_widget(w)

            # Unlocking mouse over
            self.mouse_over_ready = True

        def clear_mouse_overs():
            """
            Clears fancy mouseovers, if any
            """

            for i in [i for i in self.root_window.children
                      if isinstance(i, FancyButton) or
                      isinstance(i, FancyMarker)]:
                self.root_window.remove_widget(i)

        # Only do this routine when the filechooser screen is on
        if self.screen.name == "fc" and self.mouse_over_ready and \
                self.show_side_panel is False:
            case_bt = self.screen.ids.case_bt
            for bt in self.bookmarks_bt + [case_bt]:
                if self._determine_collision(bt, mp):
                    collision = True
                    if bt == case_bt:
                        if "Case sensitive" not in self.previous_mouse_over:
                            if case_bt.state == "down":
                                create_fancy_label("Case sensitive is ON",
                                                   case_bt,
                                                   orientation="vertical")
                            else:
                                create_fancy_label("Case sensitive is OFF",
                                                   case_bt,
                                                   orientation="vertical")

                    else:
                        # Saving relevant attributes, otherwise they would
                        # be lost
                        txt = bt.id
                        pos = bt.to_window(bt.pos[0], bt.pos[1])
                        size = bt.size

                        # If the current collision is different from the
                        # existing fancy lavel, remove it
                        if [txt] != [x.text for x in
                                     self.root_window.children if
                                     isinstance(x, FancyButton)]:
                            clear_mouse_overs()

                        Clock.schedule_once(lambda i: create_fancy_label(
                            txt, bt, adjust_pos=True, wgt_pos=pos,
                            wgt_size=size), .8)
                        self.mouse_over_ready = False

            # If no collision is detected, remove any remaining label
            # widget
            if collision is False and \
                    self.old_mouse_over in self.root_window.children:
                self.root_window.remove_widget(self.old_mouse_over)

        # Only do this routine if the side panel is open
        if self.show_side_panel and self.mouse_over_ready \
                and not self.sp_moreopts:
            # Get active tab in side panel
            active_tab = self.root.ids.main_tp.current_tab.text
            # Get remove all button
            rm_bt = [self.root.ids.rm_all_File, self.root.ids.rm_all_Taxa]
            part_bts = [self.root.ids.merge_part, self.root.ids.split_part,
                       self.root.ids.add_part]

            # Iterate over buttons of active tab
            for bt in self.mouse_over_bts[active_tab] + sidebt_list + \
                    rm_bt + part_bts:
                # Determine if there is a collision with mouse position
                if self._determine_collision(bt, mp) and self._popup not in\
                        self.root_window.children:

                    if bt in self.mouse_over_bts[active_tab]:
                        if self._determine_collision(self.root.ids.sv_file,
                                                     mp)\
                                or self._determine_collision(
                                    self.root.ids.sv_sp, mp):
                            collision = True
                        else:
                            continue
                    else:
                        # Set collision marker to true
                        collision = True
                    # This will determine if a new label button will be
                    # added to the layout, based on the text of the
                    # button. If the text is already in the previous
                    # mouse over, then do nothing. If the text is some
                    # new button, then do something
                    if bt in self.mouse_over_bts[active_tab] + sidebt_list:
                        if bt.text != self.previous_mouse_over:
                            # Check if there is an old label button and
                            # remove it
                            if self.old_mouse_over:
                                self.root_window.remove_widget(
                                    self.old_mouse_over)

                        if bt in sidebt_list:
                            pos = (bt.center_x + bt.width * .5,
                                   bt.center_y - bt.height * .5)
                            size = (200, bt.height)
                            label = create_sidebt_wgt(bt.att, pos, size)
                            show_label(mp, label)

                        elif bt in self.mouse_over_bts[active_tab]:

                            # Saving relevant attributes, otherwise they
                            # would be lost
                            txt = bt.text
                            pos = bt.to_window(bt.pos[0], bt.pos[1])
                            size = bt.size

                            # If the current collision is different from the
                            # existing fancy lavel, remove it
                            if [txt] != [x.text for x in
                                         self.root_window.children if
                                         isinstance(x, FancyButton)]:
                                clear_mouse_overs()

                            # Schedule the introduction of the label widget
                            Clock.schedule_once(lambda i:
                                create_fancy_label(txt, bt, wgt_pos=pos,
                                                   wgt_size=size,
                                                   c=(.3, .3, .3, .95)), 1)

                            # Locking mouse over so that no additional label
                            # widgets are added during the waiting time
                            self.mouse_over_ready = False

                    elif bt in rm_bt:
                        if active_tab == "Taxa" and bt == rm_bt[1]:
                            if "Removes all files and taxa" \
                                    not in self.previous_mouse_over:

                                if not bt.disabled:
                                    create_fancy_label(
                                        "Removes all files and taxa",
                                        rm_bt[1], adjust_pos=True,
                                        c=(1, .33, .33, 1))

                        elif active_tab == "Files" and bt == rm_bt[0]:
                            if "Removes all files and taxa" \
                                    not in self.previous_mouse_over:

                                if not bt.disabled:
                                    create_fancy_label(
                                        "Removes all files and  taxa",
                                        rm_bt[0], adjust_pos=True,
                                        c=(1, .33, .33, 1))

                    elif bt in part_bts:
                        bt_text = {"data/backgrounds/group_bt.png":
                            "Merge selected partitions",
                                   "data/backgrounds/split_bt.png":
                            "Split selected partitions",
                                   "data/backgrounds/add_bt35.png":
                            "Import partition scheme"}

                        if active_tab == "Partitions":
                            create_fancy_label(
                                bt_text[bt.background_normal],
                                bt, adjust_pos=True, orientation="vertical")

            # If no collision is detected, remove any remaining
            # label widget
            if collision is False and \
                    self.old_mouse_over in self.root_window.children:
                self.root_window.remove_widget(self.old_mouse_over)

        # Only do this when Statistics screen is on
        if self.screen.name == "Statistics" and not self.show_side_panel \
                and not self._popup in self.root_window.children:
            if self._determine_collision(self.screen.ids.sum_stats_bt, mp):
                collision = True
                create_fancy_label("Show summary statistics",
                                   self.screen.ids.sum_stats_bt,
                                   adjust_pos=True)

        # Only do this when plot screen is on
        if self.screen.name in self.plot_screens and self._popup not in \
                self.root_window.children:

            # Get PlotToolbar object
            toolbar_wgt = None
            try:
                toolbar_wgt = [x for x in self.root_window.children
                               if isinstance(x, OrtoPlotToolbar)][0]
            except IndexError:
                try:
                    toolbar_wgt = [x for x in self.root_window.children
                                   if isinstance(x, StatsPlotToolbar)][0]
                except IndexError:
                    toolbar_wgt = None

            if toolbar_wgt:
                # For headless plot screens a back button is added to
                # root_window
                if self.screen.name == "plot":
                    # Get back bt
                    back_bt = [x for x in self.root_window.children
                               if isinstance(x, BackButton)][0]

                    if self._determine_collision(back_bt, mp):
                        if back_bt.opacity != 1:
                            Animation(
                                opacity=1, d=.3, t="out_quart").start(back_bt)
                    else:
                        if back_bt.opacity == 1:
                            Animation(
                                opacity=.2, d=.3, t="out_quart").start(back_bt)

                # Change toolbar opacity to become visible when collision is
                # true
                if self._determine_collision(toolbar_wgt, mp):
                    if toolbar_wgt.opacity != 1:
                        Animation(
                            opacity=1, d=.3, t="out_quart").start(toolbar_wgt)

                else:
                    if toolbar_wgt.opacity == 1:
                        Animation(
                            opacity=.2, d=.3, t="out_quart").start(toolbar_wgt)

                if self.screen.name == "Statistics":

                    try:
                        toggle_wgt = [x for x in self.root_window.children
                                      if isinstance(x, StatsToggleWgt)][0]
                        if self._determine_collision(toggle_wgt, mp):
                            if toggle_wgt.opacity != 1:
                                Animation(opacity=1, d=.3, t="out_quart"). \
                                    start(toggle_wgt)
                        else:
                            if toggle_wgt.opacity == 1:
                                Animation(opacity=.2, d=.3, t="out_quart"). \
                                    start(toggle_wgt)
                    except IndexError:
                        pass

                # Check for collision with export figure or export table buttons
                if self._determine_collision(toolbar_wgt.ids.export_fig, mp):
                    collision = True
                    if "Export as graphics" not in [x.id for x in
                                                    self.root_window.children]:
                        # Create fancy label
                        create_fancy_label("Export as graphics",
                                           toolbar_wgt.ids.export_fig,
                                           line_c=(0.216, 0.67, 0.784, 1))

                elif self._determine_collision(toolbar_wgt.ids.export_table,
                                               mp):
                    collision = True
                    if "Export as table" not in [x.id for x in
                                                 self.root_window.children]:
                        # Create fancy label
                        create_fancy_label("Export as table",
                                           toolbar_wgt.ids.export_table,
                                           line_c=(0.216, 0.67, 0.784, 1))

                if self.screen.name != "Statistics":
                    if self._determine_collision(toolbar_wgt.ids.export_group,
                                                 mp):
                        collision = True
                        if "Export group" not in [x.id for x in
                                                  self.root_window.children]:
                            # Create fancy label
                            create_fancy_label("Export group",
                                               toolbar_wgt.ids.export_group,
                                               line_c=(0.216, 0.67, 0.784, 1))

        # Only do this in Orthology screen
        if (self.screen.name == "Orthology" and self.show_side_panel is
                False and self._popup not in self.root_window.children):

            id_to_txt = {"sp_vis": "Species focused exploration",
                         "gn_vis": "Ortholog focused exploration"}

            # Determine collision with add groups button
            if self._determine_collision(self.screen.ids.add_group, mp):
                collision = True
                if "Add group files" not in [x.id for x in
                                             self.root_window.children]:
                    # Create fancy label
                    create_fancy_label("Add group files",
                                       self.screen.ids.add_group,
                                       adjust_pos=True)

            # Determine collision with orthology graph visualization bts
            for bt, bt_id in zip([self.screen.ids.sp_vis,
                                  self.screen.ids.gn_vis],
                                 ["sp_vis", "gn_vis"]):
                if self._determine_collision(bt, mp):
                    collision = True
                    if id_to_txt[bt_id] not in [x.id for x in
                                                self.root_window.children]:
                        # Create fancy label
                        create_fancy_label(id_to_txt[bt_id], bt,
                                           orientation="vertical",
                                           adjust_pos=True)

        if collision is False and self.mouse_over_ready:
            self.previous_mouse_over = ""
            clear_mouse_overs()

    def switch_path_wgt(self, wgt_id, path_bx, fc_wgt):

        def path_updater():
            if os.path.exists(txt.text):
                fc_wgt.path = unicode(txt.text)
                fc_wgt.previous_dir.append(fc_wgt.path)
            else:
                return self.dialog_floatcheck(
                    "Directory does not exist", t="error")

        label = PathLabel()
        txt = PathText(id="path_editor")

        fc_wgt.bind(path=txt.setter("text"))
        fc_wgt.bind(path=label.setter("text"))

        path_bx.clear_widgets()

        if wgt_id == "label":
            label.text = fc_wgt.path
            path_bx.add_widget(label)
        else:
            txt.text = fc_wgt.path
            txt.bind(on_text_validate=lambda x: path_updater())
            path_bx.add_widget(txt)
            path_bx.children[0].focus = True

    def create_folder(self, text):

        path = self._popup.content.ids.sd_filechooser.path
        dir_name = join(path, text)

        if os.path.exists(dir_name):
            return self.dialog_floatcheck(
                "The specified folder already exists", t="error")
        else:
            os.makedirs(dir_name)
            self._popup.content.ids.sd_filechooser.path = unicode(dir_name)
            self.dismiss_subpopup()

    # ######################### SCREEN NAVIGATION ##########################

    def go_screen(self, idx, direct="left"):
        """
        Method used to go to a specific screen by specifying and index and
        transition direction
        :param idx: integer. Index value of the screen from
        self.screen_names
        :param direct: string. The direction of the transition
        """

        if self.screen:
            if self.screen.name not in self.plot_screens or \
                    self.screen.name == "Statistics":
                screen_path = join(self.cur_dir, "data", "screens",
                                   "{}.kv".format(self.screen.name))
                self.loaded_screens[screen_path] = self.screen

            # Automatic removal of plot toolbar when not in a plot screen
            if self.screen_names[idx] not in self.plot_screens:
                self.dismiss_plot_wgt()

            # Automatic removal of fancy dropdown
            self.remove_fancy_dropdown()

            # Removes old toolbar when switching directly from orto plot
            #  widget to Statistics
            if self.screen_names[idx] == "Statistics" and \
                    self.current_screen != "Statistics":
                self.dismiss_plot_wgt()

        self.index = idx

        # Precludes a transition if the current screen is the same as the
        # target screen
        if self.current_screen != self.screen_names[idx]:
            # Update previous screen
            self.previous_screen = self.current_screen
            # Update current screen
            self.current_screen = self.screen_names[idx]
            # Make the switch
            self.root.ids.sm.switch_to(self.load_screen(idx),
                                       direction=direct)

        if self.current_screen == "Statistics":
            # Generate or update statistics screen
            self.statistics_show_summary()

    def go_previous_screen(self):
        """
        Method that returns to the previous screen, set by
        self.previous_screen
        """

        if self.previous_screen != "":
            previous_idx = self.screen_names.index(self.previous_screen)
            self.go_screen(previous_idx, "right")

    def load_screen(self, idx):
        """
        Loads the current screen according to the corresponding kv file
        :param idx: The index of the screen to be loaded
        """

        if self.loaded_screens[self.available_screens[idx]]:
            self.screen = self.loaded_screens[self.available_screens[idx]]
            if self.screen.name == "fc":
                self.screen.ids.icon_view_tab.selection = []
        else:
            self.screen = Builder.load_file(self.available_screens[idx])

            # If the screen to be loaded is the filechooser, set the
            # home path as the default
            if basename(self.available_screens[idx]) == "fc.kv":
                self.screen.ids.icon_view_tab.path = self.home_path
                # Initialize bookmarks
                self.bookmark_init(self.screen.ids.sv_book,
                                   self.screen.ids.sv_mycomp,
                                   self.screen.ids.icon_view_tab)

        # Disengage acion bar toggle buttons when entering main filechooser
        if basename(self.available_screens[idx]) == "fc.kv":
            self.disengage_groups()

        if basename(self.available_screens[idx]) == "Statistics.kv":
            # Add StatsToggleWidget, if present
            try:
                if isinstance(self.screen.ids.plot_content.children[0].
                        children[0], Image):
                    self.show_plot_toolbar(toolbar_type="stats")
                    if self.previous_stats_toggle:
                        self.root_window.add_widget(
                            self.previous_stats_toggle)
            except IndexError:
                pass

        return self.screen

    def go_carousel(self, slide, bt_id):
        """
        Method used by other buttons outside the side buttons of the side
        panel to go to specific slides of the side panel
        :param slide: int, the index of the target slide
        :param bt_id: string, the id of the corresponding button
        :return:
        """

        self.toggle_groups(self.root.ids[bt_id])
        self.root.ids[bt_id].state = "down"

        panel_car = self.root.ids.carousel
        panel_car.load_slide(panel_car.slides[slide])

    @staticmethod
    def toggle_groups(wgt):
        """
        This method generates a desired behaviour for groups of toggle
        buttons By default, when a toggle button is pressed, the state
        will be down and a new screen/slide is presented. However,
        if the same toggle button is pressed again, it's state will
        return to the normal state while the same screen/slide is showed.
        To prevent this behaviour, this method will disable the active
        toggle button in the group and enable any other previously
        disabled button.

        To allow a seamless transition, ensure that background_disabled_dow
        is the same as background_down, and that disabled_color is the
        same as color.

        :param wgt: The toggle button widget. Must belong to a group.
        """

        # Iterating over the children of the parent may not be optimal, but
        # using the get_widgets(groupname) method could result in some
        # issues with garbage collector of kivy. So, for now, this will
        # iterate over all children of the toggle button's parent
        for i in [x for x in wgt.parent.children
                  if isinstance(x, ToggleButton)]:
            if i.disabled:
                i.disabled = False
                i.state = "normal"

        wgt.disabled = True

    def disengage_groups(self):
        """
        Methods used to unselect all module header buttons, in case the user
        enters the main filechooser
        """

        for wgt in [self.root.ids.h_ortho, self.root.ids.h_process,
                    self.root.ids.h_stat]:
            wgt.disabled = False
            wgt.state = "normal"

    def check_action(self, text, func, bt_wgt=None, args=None,
                     popup_level=1, check_wgt=CheckDialog):
        """
        General purpose method that pops a dialog checking if the user
        wants to perform a certain action. This method should be passed
        as a function on the 'on_*' with the final function and original
        widget triggering the event as arguments.

        By default, the check action dialog is a CheckDialog widget,
        but alternative widgets can be used by passing them to the
        check_wgt option. If so, properties such as size, separator
        color, etc must be specified for such cases.

        :param text: string, text to appear in the dialog
        :param func: final function if the users chooses to proceed
        :param bt_wgt: widget where the initial 'on_' event occurred
        :param args: list, of arguments to be passed on to func
        :param popup_level: int, level of popup. 1 for _popup, 2 for
        _subpopup and 3 for exit popup
        :param check_wgt: Widget. Specified the check dialog widget

        Usage example:
        This can be applied to the removal button of the bookmarks. In this
        case, the event of the removal button must be like this:

        remove_button.bind(partial(self.check_action,
                           self.remove_bookmark_bt))

        where, self.check_action is this method, and self.remove_bookmark_bt
        is the function that will actually remove the bookmark button. This
        function is then bound to the "OK" button of the check dialog.
        By default, the last argument is the bt_wgt.
        """

        if popup_level == 1:
            check_content = check_wgt(cancel=self.dismiss_popup)
        elif popup_level == 2:
            check_content = check_wgt(cancel=self.dismiss_subpopup)
        else:
            # Do nothing if exit popup already exists. The return True
            # prevents the app from closing, since in that stage it will
            # be listening to return False to close.
            if self._exit_popup in self.root_window.children:
                return True

            check_content = check_wgt(cancel=self.dismiss_exit)

        if isinstance(check_content, CheckDialog):
            # Set size
            size = (300, 200)
            # Sep separator color
            sep_color = self._red
            # Set popup title
            title = "Warning!"
            # Set text for dialog
            check_content.ids.check_text.text = text
        elif isinstance(check_content, CheckProject):
            title = "Project overview"
            size = (350, 230)
            sep_color = self._blue
            check_content.ids.proj_name.text = text[0]
            check_content.ids.file_num.text = text[1]

        # Set the function to be called when dispatching the on_release
        # event of the check_ok button.
        if bt_wgt and not args:
            check_content.ids.check_ok.bind(
                on_release=lambda val: func(bt_wgt))
        else:
            if args:
                check_content.ids.check_ok.bind(on_release=lambda val:
                func(*args))
            else:
                check_content.ids.check_ok.bind(
                    on_release=lambda val: func())

        if popup_level in [1, 2]:
            self.show_popup(title=title, content=check_content,
                            size=size, popup_level=popup_level,
                            separator_color=sep_color)
        elif popup_level == 3:
            self._exit_popup = CustomPopup(
                title="[b]%s[/b]" % title,
                content=check_content,
                size=size,
                size_hint=(None, None),
                separator_color=sep_color,
                border_color=sep_color,
                separator_height=0.6,
                title_color=sep_color,
                background="data/backgrounds/transparent.png")
            self._exit_popup.open()

        return True

    def check_file(self, path, file_name, idx, bw=None):
        """
        Method used by some filechooser dialogs. Checks whether the provided
        file name already exists. If so, issues a check_action popup. If
        not, proceeds as normal
        :param path: string, complete path
        :param file_name: string, file name
        :param idx: string, operation identifier
        :param bw: bool, argument specific for export_graphic idx.
        Determines whether a figure will be exported in original color or
        grayscale
        """

        # Stores methods. key: idx; first value element, method to apply;
        # second value element, list of arguments; third value element the
        # file_name extension
        methods = {
            "main_output":
            [self.save_file, [path, file_name, idx], ""],
            "export":
            [self.export_names, [path, file_name], ".txt"],
            "export_table":
            [self.export_table, [path, file_name], ".csv"],
            "export_outliers":
            [self.export_table, [path, file_name, ".txt"], ".txt"],
            "export_graphic":
            [self.export_graphic, [path, file_name, bw], ""],
            "group":
            [self.orto_export_groups, ["group", path, file_name], ""]
        }

        # Check if files exists.
        try:
            p = join(unicode(path, "utf8"), file_name + methods[idx][2])
        except TypeError:
            p = join(path, file_name + methods[idx][2])

        if os.path.exists(p) and idx != "main_output":

            self.check_action(
                "The file {} already exists. Overwrite?".format(file_name),
                methods[idx][0],
                **{"args": methods[idx][1], "popup_level": 2})

        else:
            methods[idx][0](*methods[idx][1])

        # Close popup automatically except for the group idx
        if idx != "group":
            self.dismiss_popup()

    # ########################## GENERAL USE ###############################

    def run_in_background(self, func, second_func, args1, args2=None,
                          no_arg2=False, msg="Crunching data...",
                          cancel=True):
        """
        This method is intended to run time/resource consuming operations in
        the background, without freezing the app, and return the final
        result to the main thread. This means that complex methods that
        perform multiple changes to the App's attributes are not suitable
        for this method  (changes in a secondary thread will not change
        the App in the main thread). Therefore, the simplest solution to
        the problem is to isolate the time consuming parts of those
        methods, run them in the background, and get their result. Then,
        their result is piped to the follow-up fuction provided as argument
        :param func: intensive callable bound method to run in the
        background
        :param second_func: Follow-up bound method that will use the value
        returned by func
        :param args1: list, with the arguments for func t method. No
        keywords allowed
        :param args2: list, with arguments for second_func. These will be
        added to the argument list returned by func
        :param no_arg2: Boolean. Whether func will return something to
        second_func
        :param msg: string, message to appear in waiting dialog.
        :param cancel: Boolean, Shows a functional cancel button that
        interrupts the background process, if True.
        """

        def check_process_status(p, second_function, args, man, dt):
            """
            This scheduled function will check the status of the second
            process. When finished, it will dismiss the waiting popup,
            get the value returned by the second process and issue the
            following method

            :param p: Process object.
            :param second_function: function object. If not None, this
            function will be executed in the main process, once the first
            process is finished.
            :param args: list. Contains the argument to be passed at the
            second_function
            :param man: Manager object.
            """

            # Kill switch. Terminating thread by user
            if self.terminate_background:

                # Issuing the kill order to the child thread. This signal
                # will propagate through the worker's methods and graciously
                # terminate
                shared_ns.stop = True

                # This small delay seems to fix manager shutdown issues
                time.sleep(.1)

                # Shutting down manager that provides communication between
                # parent and child thread
                man.shutdown()

                # Close loading popup dialog
                self.dismiss_popup()

                # Unschedule the current function
                Clock.unschedule(check_func)

                # Join child process and exit
                p.join()
                return

            if "img" in self._popup.content.ids:
                self._popup.content.ids.img.rotation -= 10

                if shared_ns.counter:
                    self._popup.content.ids.msg2.text = "{}/{}".format(
                        shared_ns.counter, shared_ns.total)

            if not p.is_alive():

                try:
                    if shared_ns.exception:
                        return self.dialog_floatcheck(
                            "An unexpected error occurred. Check the app"
                            " logs", t="error")
                except Exception as e:
                    Logger.exception(e.message)

                try:
                    val = shared_ns.val
                except AttributeError:
                    val = None

                Clock.unschedule(check_func)
                self.dismiss_popup()

                # Manager must be shutdown before closing the process,
                # otherwise the data pipe will prevent the process from
                # closing
                man.shutdown()
                p.join()

                # Checks if there is a second function to run and whether
                # there are additional arguments for secondary function
                if not no_arg2:
                    if second_func:
                        if args2 and val:
                            val.extend(args)
                        elif args2:
                            val = args
                        second_function(*val)
                elif second_function:
                    second_function()

        # Set up manager and shared name space to share information
        # between background and main processes
        manager = multiprocessing.Manager()
        shared_ns = manager.Namespace()
        shared_ns.counter = 0

        # Set kill switch flag
        shared_ns.stop = False

        # Remove lock from background process
        self.terminate_background = False

        # Create process
        p = threading.Thread(target=background_process,
                             args=(func, shared_ns, args1))
        p.daemon = True
        p.start()

        # Remove any possible previous popups
        self.dismiss_popup()
        # Create waiting dialog
        content = CrunchData()
        # Set label
        content.ids.msg.text = msg

        # Determine whether widget should have cancel option or not
        if cancel:
            # Set functionality of cancel button
            size = (230, 230)
            content.ids.cancel_bt.bind(
                on_release=lambda x: setattr(self, "terminate_background",
                                             True))
        else:
            # Remove cancel button
            size = (230, 200)
            content.ids.cancel_anc.remove_widget(content.ids.cancel_bt)
            content.remove_widget(content.ids.cancel_anc)

        # Create popup with waiting dialog
        self.show_popup(title="", content=content, size=size,
                        separator_color=(0, 0, 0, 0),
                        border_color=tm.c_popup_border,
                        auto_dissmiss=False)

        # Schedule function that checks the process' pulse
        check_func = partial(check_process_status, p,
                             second_func, args2, manager)
        Clock.schedule_interval(check_func, .1)

    # ###################### BOOKMARKS OPERATIONS ##########################

    def bookmark_init(self, wgt, dev_wgt, fc_wgt, popup_level=1):
        """
        This will create a pickle file containing a list with the bookmarks
        for the file chooser menu. If no file exists, it will create an
        empty one. If a file already exists, it will load the available
        bookmarks

        :param wgt: Widget object where the bookmark button will be
        :param dev_wgt: Widget object where the system bookmarks will be
        :param fc_wgt: Filechooser widget associated with the bookmarks
        :param popup_level: int, specifies the level of the check_action popup
         that appears when attempting to remove a bookmark
        """

        if exists(self.bm_file):
            self.bookmarks = pickle.load(open(self.bm_file, "rb"))
            # Retrieving the bookmark path list from the self.bookmarks
            bk_list = self.bookmarks[0]
            for bk in bk_list:
                self.add_bookmark_bt(bk, wgt, fc_wgt, popup_level=popup_level)

        else:
            pickle.dump(self.bookmarks, open(self.bm_file, "wb"))

            # Trying to import FM bookmarks in ~/.config/gtk-3.0/bookmarks.
            # This will happen only the first time the app is executed, when
            # the bookmarks will be saved. From there on, bookmarks will be
            # managed in the app
            if sys.platform in ["linux", "linux2"]:
                if exists(join(self.home_path, ".config", "gtk-3.0",
                               "bookmarks")):
                    with open(join(self.home_path, ".config", "gtk-3.0",
                                   "bookmarks")) as bk_file:
                        for bk_line in bk_file:
                            bk = bk_line.split()[0].replace("file://", "")
                            # urllib will ensure special characters with
                            # punctuation are correctly showed
                            bk = urllib.unquote(bk)
                            self.save_bookmark(bk, wgt, fc_wgt,
                                               popup_level=popup_level)

        # Get main paths for linux
        if sys.platform in ["linux", "linux2"]:

            # System
            self.add_bookmark_bt("/", dev_wgt, fc_wgt, name="System",
                                 rm_bt=False)
            # Home
            self.add_bookmark_bt(self.home_path, dev_wgt, fc_wgt,
                                 name="Home", rm_bt=False)

            # Get removable media
            x = subprocess.Popen(["mount | grep -e /media -e /mnt | awk '{"
                "print $3}'"], shell=True, stdout=subprocess.PIPE)
            removable_media = x.communicate()[0]
            if removable_media:
                for path in removable_media.split("\n")[:-1]:
                    # Add path to removable_media attribute. This attribute
                    # is used when updating removable media real time in
                    # _check_removable_media
                    if path not in self.removable_media:
                        self.removable_media.append(path)

                    name = basename(path)
                    self.add_bookmark_bt(path, dev_wgt, fc_wgt, name,
                                         rm_bt=False)

        # Get main devices for windows
        elif sys.platform in ["win32", "cygwin"]:

            devices = re.findall(
                r"[A-Z]+:.*$", os.popen("mountvol /").read(), re.MULTILINE)

            # Prevent no disk on drive popup error in windows
            SEM_FAILCRITICALERRORS = 1
            SEM_NOOPENFILEERRORBOX = 0x8000
            SEM_FAIL = SEM_NOOPENFILEERRORBOX | SEM_FAILCRITICALERRORS

            kernel32 = ctypes.WinDLL('kernel32')
            oldmode = ctypes.c_uint()
            kernel32.SetThreadErrorMode(SEM_FAIL, ctypes.byref(oldmode))

            for d in devices:
                if exists(d):
                    self.removable_media.append(d)
                    self.add_bookmark_bt(d, dev_wgt, fc_wgt, rm_bt=False,
                                         name=os.path.splitdrive(d)[0])

            # Home directory
            self.add_bookmark_bt(self.home_path, dev_wgt, fc_wgt, rm_bt=False)

        # Get some favorites for MacOS
        elif sys.platform == "darwin":

            # Get mounted volumes
            for d in os.listdir("/Volumes"):
                self.removable_media.append(join("/Volumes", d))
                self.add_bookmark_bt(join("/Volumes", d), dev_wgt, fc_wgt,
                                     name=d, rm_bt=False)

            # Home
            self.add_bookmark_bt(self.home_path, dev_wgt, fc_wgt,
                                 name=basename(self.home_path),
                                 rm_bt=False)

            # Desktop
            self.add_bookmark_bt(join(self.home_path, "Desktop"), wgt,
                                 fc_wgt, "Desktop")

            # Documents
            self.add_bookmark_bt(join(self.home_path, "Documents"), wgt,
                                 fc_wgt, "Documents")

            # Downloads
            self.add_bookmark_bt(join(self.home_path, "Downloads"), wgt,
                                 fc_wgt, "Downloads")


    def save_bookmark(self, path, wgt, fc_wgt, popup_level=1):
        """
        This adds functionality to the FileChooser "Add bookmark" button. It
        will grab the selected path and add it to a storage list that
        will be saved as a pickle object and stored in a file defined in
        self.bm_file.
        :param path: String containing the path of the bookmark
        :param wgt: Widget where the bookmark will be added
        :param fc_wgt: FileChooser widget that the bookmark will affect
        :param popup_level: int, specifies the level of the check_action popup
         that appears when attempting to remove a bookmark
        """

        # Load bookmarks object
        self.bookmarks = pickle.load(open(self.bm_file, "rb"))
        # Check if bookmark already exists. Only add bookmark if it does not
        # exist
        if path not in self.bookmarks[0]:

            # Add bookmarks to the full path list
            self.bookmarks[0].append(path)
            # Add mapping of the full path to the bookmark name
            new_map = {basename(path): path}
            self.bookmarks[1] = dict(list(self.bookmarks[1].items()) +
                                     list(new_map.items()))
            self.add_bookmark_bt(path, wgt, fc_wgt, popup_level=popup_level)
            pickle.dump(self.bookmarks, open(self.bm_file, "wb"))

    def add_bookmark_bt(self, bk, wgt, fc_wgt, name=None, rm_bt=True,
                        popup_level=1):
        """
        This will add a bookmark button, along with its removal button. Only
        a bookmark path will be necessary.

        The path of the bookmark will be associated to the respective button
        by it's id.

        :param bk: string. bookmark file path
        :param wgt: Widget, preferentially a gridlayout where the bookmark
        buttons will be added
        :param fc_wgt: FileChooser widget in which bookmark operations will
        be performed
        :param name: string, optional name for bookmark button instead of
        the basename of the path
        :param rm_bt: Boolean, If True, a removal button will be added with
        the bookmark, else the removal button will not be added. The latter
        case is used for System devices bookmarks.
        :param popup_level: int, specifies the level of the check_action popup
         that appears when attempting to remove a bookmark
        """

        # Ignore if path does not exist / is not a directory
        if not os.path.isdir(bk):
            return

        bookmark_name = basename(bk)
        # Define bookmark button
        bt = TFButton(text=name if name else bookmark_name, id=bk,
                      bold=True, height=30, size_hint=(.8, None),
                      background_normal=join("data", "backgrounds",
                                             "bt_process.png"),
                      background_down=join("data", "backgrounds",
                                           "bt_process_off.png"))
        # Bind to function that loads bookmark path into filechooser
        bt.bind(on_release=lambda x: self.bookmark_load(x, fc_wgt))
        # Add to list for mouse over purposes
        self.bookmarks_bt.append(bt)
        wgt.add_widget(bt)

        if rm_bt:
            # Define bookmark removal button
            xbt = Button(size_hint=(None, None), width=30,
                         height=30, id="%sX" % bk, border=(0, 0, 0, 0),
                         background_normal=join("data", "backgrounds",
                                                "remove_bt.png"),
                         background_down=join("data", "backgrounds",
                                              "remove_bt_down.png"))
            # Bind to function that removes bookmark button as well as the
            # path from self.bm_file
            xbt.bind(on_release=partial(self.check_action,
                                        "Are you sure you want to remove"
                                        " this bookmark?",
                                        self.remove_bookmark_bt,
                                        popup_level=popup_level))
            wgt.add_widget(xbt)

    def bookmark_load(self, value, wgt):
        """
        Provided a bookmark button object, it loads the bookmark file path
        that is stored in the button id.
        :param value: bookmark button object
        :param wgt: Filechooser widget that will show the bookmark
        """

        path = value.id
        if os.path.exists(path):
            try:
                wgt.previous_dir.append(wgt.path)
            except KeyError:
                pass
            wgt.path = unicode(path)
            wgt.selection = []
        else:
            self.dialog_floatcheck(
                "The path to the selected bookmark no longer exists.",
                t="error")

    def remove_bookmark_bt(self, value):
        """
        Adds functionality to the removal button associated with each
        bookmark button. This will not only remove the
        :param value: The removal button widget
        """

        # Get the widget to remove the button
        parent_obj = value.parent

        # Get the bookmark button, using the removal button id
        bk_idx = value.id[:-1]
        bk_bt = [x for x in parent_obj.children if bk_idx == x.id][0]

        # Remove both bookmark and removal buttons
        parent_obj.remove_widget(value)
        parent_obj.remove_widget(bk_bt)

        # Update gridlayout height for scrolling purposes
        parent_obj.height -= self.root.height * 0.07

        # Core changes
        bk_name = basename(bk_idx)
        # Remove bookmark path from list and mapping dictionary
        self.bookmarks[0].remove(bk_idx)
        del self.bookmarks[1][bk_name]
        # Update self.bm_file
        pickle.dump(self.bookmarks, open(self.bm_file, "wb"))

    # ####################### PLOT SCREEN OPERATIONS #######################

    def dialog_export_graphic(self):
        """
        Creates a filechooser dialog for graphics exportation. It differs
        from other filechooser dialogs in the presence of a spinner to
        select the graphical extension
        :return:
        """

        content = ExportGraphics(cancel=self.dismiss_popup,
                                 bookmark_init=self.bookmark_init)

        self.show_popup(title="Export as graphic...", content=content,
                        size_hint=(.9, .9))

    def export_graphic(self, path, file_name, bw):
        """
        Saves the current plot object into a file based on file name and
        extension
        :param path: string, path to final directory
        :param file_name: string, name of graphic file
        :param bw: bool, Determined whether the graphic will be in the
        original colors (False) or grayscale (True)
        """

        temp_f = None
        f_path = join(path, file_name)

        if bw:
            temp_f = join(self.temp_dir, file_name)
            file_name = os.path.splitext(file_name)[0] + "_bw" + \
                os.path.splitext(file_name)[1]
            f_path = join(path, file_name)

        if self.current_lgd:
            self.current_plot.savefig(temp_f if temp_f else f_path,
                bbox_extra_artists=(self.current_lgd,),
                bbox_inches="tight")
        else:
            self.current_plot.savefig(temp_f if temp_f else f_path,
                                      bbox_inches="tight")

        if bw:
            image = misc.imread(temp_f, flatten=True)
            misc.imsave(f_path, image)
            os.remove(temp_f)

        self.dialog_floatcheck("Graphic successfully exported!", t="info")

    def export_table(self, path, file_name, ext=".csv", warning_dlg=True):
        """
        Saves the current_table list attribute to a .csv file.
        :param path: string, path to final directory
        :param file_name: string, name of table file
        :param warning_dlg: boolean. If True, issue a float dialog when
        the operation is done. If false do not issue the float dialog.
        """

        # Create table file object handle
        table_handle = open(join(path, file_name + ext), "w")

        # Writing table. Each entry in self.current_table should represent a
        # line in the table
        for l in self.current_table:
            table_handle.write(";".join([str(x) for x in l]) + "\n")

        # Close file object
        table_handle.close()

        if warning_dlg:
            self.dialog_floatcheck("Table successfully exported!", t="info")

    # ####################### SIDE PANEL OPERATIONS ########################

    def dialog_about(self):
        """
        Dialog with the about information on TriFusion
        """

        content = AboutDialog()

        content.version = __version__
        content.build = __build__

        self.show_popup(title="About TriFusion", content=content,
                        close_bt=True, size=(480, 410))

    def sidepanel_on_touch(self, touch):
        """
        This function is binded to the app Window so that it can handle any
        touch_up events. Once the side panel is open, this allows any mouse
        click outside the panel to close it. It gathers information on the
        mouse and side panel position and evaluates a collision. It will
        trigger the side panel closing only when four conditions are met:

        - When there is a mouse input outside the side panel
        - When the variable controlling the side panel (show_side_panel) is
        True, meaning that the panel is extended
        - When the mouse input is outside the previous button in the action
        bar, which is also used to toggle the side panel. This prevents
        issues of toggling the side panel twice with one mouse input
        - When a popup is not open. There are several buttons in the side
        bar that open whose position is outside the side bar. The user
        should be
        able to click anywhere in the popup without the side panel closing.

        In addition, it will handle the status of the partition box dialog,
        associated with the sidepanel. While this box is active, the
        sidepanel must remain open. However, clicks outside the partition
        box will close it.

        :param touch: Touch event
        """

        # Set touch as a app attribute
        self.touch = touch

        def animate_sidebar():

            # ANIMATIONS with hierarchy
            # Animation of main BoxLayout containing child ScrollViews
            self.sidepanel_animation(width=0, wgt=self.root.ids.main_box)
            # Animation of both scrollviews
            self.sidepanel_animation(width=0, wgt=self.root.ids.sp)
            self.sidepanel_animation(width=0, wgt=self.root.ids.sp_bts)

            self.show_side_panel = not self.show_side_panel

        # Get mouse position
        mp = self.root_window.mouse_pos
        # Get side panel and previous button widgets
        side_panel_wgt = self.root.ids.main_box
        ap = self.root.ids.ap

        # If sidepanel's more options widget is active, check if touch is
        # outside. If so, remove widget
        if self.sp_moreopts:
            wgt = [x for x in self.root_window.children if
                   isinstance(x, SP_MoreOpts_Dialog)][0]

            if not wgt.collide_point(mp[0], mp[1]) and not \
                    self._determine_collision(self.root.ids.file_opt, mp):
                self.sidepanel_remove_moreopts()

        # If fancy dropdown is dropped, check if touch is oustide. If so,
        # remove the widget
        if self.fancy_dropped:
            wgt = [x for x in self.root_window.children if
                   isinstance(x, FancyDropDown)][0]

            if not self._determine_collision(wgt, mp):
                self.remove_fancy_dropdown()

        # Check for existence of a partition dialog box
        partition_box = [x for x in self.root_window.children if
                         isinstance(x, PartitionsDialog)]

        # If the partition box exists and the collision is outside it
        if partition_box and not partition_box[0].collide_point(
                mp[0], mp[1]):
            # Check if spinner is open
            spin1 = partition_box[0].ids.codon_spin.is_open
            spin2 = [x.is_open for x in
                     partition_box[0].ids.model_bx.children]

            # If the spinners are not open, remove
            if True not in spin2 and not spin1:
                rm_bt = [x for x in self.root_window.children if
                         isinstance(x, RemoveFloat)][0]
                self.root_window.remove_widget(partition_box[0])
                self.root_window.remove_widget(rm_bt)

        # Check for conditions to close the side panel.
        # If touch is out of panel; if panel is open; is touch is out of
        # menu button; a popup is not open
        if side_panel_wgt.collide_point(mp[0], mp[1]) is False \
                and self.show_side_panel \
                and ap.collide_point(mp[0], mp[1]) is False \
                and self._popup not in self.root_window.children \
                and not partition_box \
                and not self.sp_moreopts:

            if self.screen.name == "Process" or \
                            self.screen.name == "Orthology":
                queue_bt = self.screen.ids.queue_bt
                if queue_bt.collide_point(mp[0], mp[1]) is False:
                    animate_sidebar()
            else:
                animate_sidebar()

    def toggle_sidepanel(self):
        """
        Method controlling the animation toggling of the side panel
        """

        if self.lock_side_panel:
            return

        # Closes partition box, if open
        self.remove_partition_box()

        # Do not toggle when the following conditions are met
        if self.screen.name == "fc":
            if self.screen.ids.text_filter.focus:
                return
            if isinstance(self.screen.ids.path_bx.children[0], PathText):
                if self.screen.ids.path_bx.children[0].focus:
                    return

        # Toggling the state of the panel. This attribute is the main
        # controller of the side panel state. When its True, the side panel
        #  is extended, otherwise the side panel is hidden
        self.show_side_panel = not self.show_side_panel

        # Lock the sidepanel animation. It will only be released after the
        # animations are over
        self.lock_side_panel = True
        Clock.schedule_once(lambda x: setattr(self, "lock_side_panel", False),
                            .4)

        if self.show_side_panel:

            # Redraw the side panel layout. This will ensure that the widget
            # is always on top of all widgets.
            self.root.ids.bx1.remove_widget(self.root.ids.panel_float)
            self.root.ids.bx1.add_widget(self.root.ids.panel_float)

            # Fixing the width of the side panel
            # Main panel width
            sv_panel_width = 330
            # Side buttons width
            sv_bts_width = 60
        else:
            sv_panel_width, sv_bts_width = 0, 0

        # ANIMATIONS with hierarchy
        # Animation of main BoxLayout containing child ScrollViews
        self.sidepanel_animation(width=sv_panel_width * 1.2,
                                 wgt=self.root.ids.main_box)
        # Animation of both scrollviews
        self.sidepanel_animation(width=sv_panel_width,
                                 wgt=self.root.ids.sp)
        self.sidepanel_animation(width=sv_bts_width,
                                 wgt=self.root.ids.sp_bts)

    @staticmethod
    def sidepanel_animation(width, wgt):

        Animation(width=width, d=.3, t="out_quart").start(wgt)

    def sidepanel_search_bts(self, txt, panel):
        """
        Performs a search for file or taxa buttons in the side panel based
        on the text string provided
        :param txt: string, the expression used for the search
        :param panel: string, the panel to perform the search. Can be either
        'files' or 'taxa'.
        """

        if panel == "files":
            # Setting which original button list
            bt_list = self.file_list
            # Setting which sink grid layout
            gl_wgt = self.root.ids.file_sl
            create_bts_mtd = self.sidepanel_create_bts
        elif panel == "taxa":
            bt_list = sorted(self.alignment_list.taxa_names)
            gl_wgt = self.root.ids.taxa_sl
            create_bts_mtd = self.sidepanel_create_bts
        else:
            bt_list = list(self.alignment_list.partitions.partitions.keys())
            gl_wgt = self.root.ids.partition_sl
            create_bts_mtd = self.sidepanel_create_part_bts

        # Find buttons that match the txt string
        if panel == "files":
            found_bts = [basename(el) for el in bt_list if
                         txt.lower() in basename(el).lower()]
        else:
            found_bts = [el for el in bt_list if txt.lower() in el.lower()]

        # Clear the grid and populate with the found bts
        gl_wgt.clear_widgets()
        mouse_bts = []
        for txt in found_bts:
            if panel == "partitions":
                # Get number of alignment for partition
                fls = self.alignment_list.partitions.\
                    partitions_alignments[txt]
                bt, inf_bt, x_bt = create_bts_mtd([txt, fls])
            else:
                bt, inf_bt, x_bt = create_bts_mtd(txt)
            gl_wgt.add_widget(bt)
            gl_wgt.add_widget(inf_bt)
            gl_wgt.add_widget(x_bt)
            mouse_bts.append(bt)

        if panel == "files":
            self.mouse_over_bts["Files"] = mouse_bts
        elif panel == "taxa":
            self.mouse_over_bts["Taxa"] = mouse_bts
        else:
            self.mouse_over_bts["Partitions"] = mouse_bts

    def sidepanel_clear_search(self, panel):
        """
        Clears previous search string and populates with the original
        buttons
        :param panel: string, the panel to clear the search. Can be either
        'files' or 'taxa'
        """

        if panel == "files":
            # Setting which original button list
            bt_list = self.sp_file_bts
            # Setting which sink grid layout
            gl_wgt = self.root.ids.file_sl
        elif panel == "taxa":
            bt_list = self.sp_taxa_bts
            gl_wgt = self.root.ids.taxa_sl
        else:
            bt_list = self.sp_partition_bts
            gl_wgt = self.root.ids.partition_sl

        # Clear the grid and populate with the found bts
        gl_wgt.clear_widgets()
        mouse_bts = []
        for bt, inf_bt, rm_bt in bt_list:
            # Update button states during search operation
            # Default states is down for Taxa/Files and normal for
            # Partitions
            state = "down" if gl_wgt != self.root.ids.partition_sl else \
                "normal"
            # For files
            if bt.id in self.filename_map:
                if self.filename_map[bt.id] not in self.active_file_list:
                    state = "normal"
            # For taxa
            elif bt.text not in self.active_taxa_list:
                state = "normal"
            bt.state = state
            gl_wgt.add_widget(bt)
            gl_wgt.add_widget(inf_bt)
            gl_wgt.add_widget(rm_bt)
            mouse_bts.append(bt)

        if panel == "files":
            self.mouse_over_bts["Files"] = mouse_bts
        elif panel == "taxa":
            self.mouse_over_bts["Taxa"] = mouse_bts
        else:
            self.mouse_over_bts["Partitions"] = mouse_bts

        if panel != "taxa":
            try:
                if self.file_list[self.count_files + 1]:
                    gl_wgt.add_widget(LoadMoreBt())
            except IndexError:
                return

    def sidepanel_remove_moreopts(self):
        """
        Removes widgets from the more options dialog
        """

        for wgt in [x for x in self.root_window.children if
                    isinstance(x, SP_MoreOpts_Dialog) or
                    isinstance(x, FancyMarkerPersist)]:
            self.root_window.remove_widget(wgt)

        self.sp_moreopts = False

    def remove_fancy_dropdown(self):
        """
        Removes FancyDropDown widgets from the root window
        :return:
        """

        for wgt in [x for x in self.root_window.children if
                    isinstance(x, FancyDropDown) or
                    isinstance(x, FancyMarkerPersist)]:
            self.root_window.remove_widget(wgt)

        self.fancy_dropped = False

    def toggle_fancy_dropdown(self, bt, values, orientation="down",
                              add_ds_type=None):
        """
        Adds functionality to the File/Taxa dropdown menus using a
        FancyDropDown widget
        :param Button bt: Button instance that was pressed
        :param list values: List of values that will populate the dropdown
        menu
        :param str orientation: Whether the dropdown will go 'down' or 'up'
        :param str add_ds_type: If not None, it will provide the additional
        data set buttons of the appropriate ds_type ('files' or 'taxa')
        """

        def set_ds(cur_bt, ds_name, wgt):
            cur_bt.text = ds_name

            self.remove_fancy_dropdown()

        if self.fancy_dropped:
            return self.remove_fancy_dropdown()

        # Get vertical size from the number of entries in values
        ds_len = len(values)
        # If the add_ds_type option is provided, add two more items that will
        # correspond to the "active" and "total" datasets
        if add_ds_type:
            ds_len += 2

        # Since bt may be inside a scrollview, this retrieves the real position
        # of the widget
        real_pos = bt.to_window(bt.pos[0], bt.pos[1])

        width = 200
        height = ds_len * 35 + 19

        # Get position
        wgt_x = real_pos[0] - (width - bt.width) / 2
        if orientation == "down":
            wgt_y = real_pos[1] - 10 - height
        else:
            wgt_y = real_pos[1] + bt.height + 17

        # Correct x position, if outside screen
        if wgt_x + width > self.root.width:
            # Determine the width that is over the window bound
            offset = (wgt_x + width) - self.root.width
            wgt_x -= (offset + 20)

        # Correct x position to the left
        wgt_x = wgt_x if wgt_x >= 20 else 20

        # Correct height and y position, if widget falls outside screen
        if wgt_y <= 0:
            height += (wgt_y - 20)
            wgt_y = 20

        elif wgt_y + height >= self.root.height:
            off_set = (wgt_y + height) - self.root.height
            height -= (off_set + 10)

        # Get final size
        size = [width, height]

        # Generate fancy marker
        bg = join("data", "backgrounds", "box_arrow_up.png") if \
            orientation == "down" else \
            join("data", "backgrounds", "box_arrow_down.png")

        if orientation == "down":
            arrow_pos = (real_pos[0] + (bt.width / 2), real_pos[1] - 10)
        else:
            arrow_pos = (real_pos[0] + (bt.width / 2),
                         real_pos[1] + bt.height + 10)

        point_wgt = FancyMarkerPersist(
            background_normal=bg, pos=arrow_pos, size=(12, 7),
            background_color=(0.216, 0.67, 0.784, 1))

        # Generate DropDown
        dropdown = FancyDropDown(pos=(wgt_x, wgt_y), size=size)

        if add_ds_type:
            for name_tpl in ["All {}", "Active {}"]:
                name = name_tpl.format(add_ds_type)
                bg = "data/backgrounds/transparent.png" if name != bt.text \
                    else "data/backgrounds/bt_process.png"
                ds_bt = TFButton(size_hint_y=None, height=30,
                                 text=name,
                                 background_normal=bg)
                ds_bt.bind(on_release=partial(set_ds, bt, name))
                dropdown.ids.grid_wgt.add_widget(ds_bt)

        # Create Buttons for each custom data set group
        for name in values:
            bg = "data/backgrounds/transparent.png" if name != bt.text else \
                "data/backgrounds/bt_process.png"
            ds_bt = TFButton(size_hint_y=None, height=30, text=name,
                             background_normal=bg)
            ds_bt.bind(on_release=partial(set_ds, bt, name))
            dropdown.ids.grid_wgt.add_widget(ds_bt)

        self.root_window.add_widget(point_wgt)
        self.root_window.add_widget(dropdown)

        self.fancy_dropped = True

    def sidepanel_moreopts_dialog(self, bt):
        """
        """

        if self.sp_moreopts:
            return self.sidepanel_remove_moreopts()

        # Get active tab
        active_tab = self.root.ids.main_tp.current_tab.text

        # Get position
        wgt_x = bt.x + (bt.width * 2) - 135
        wgt_y = bt.y + bt.height + 12

        # Generate fancy marker
        point_wgt = FancyMarkerPersist(
            background_normal=join("data", "backgrounds",
                                   "box_arrow_down.png"),
            pos=(bt.x + (bt.width * 2), (bt.y + bt.height + 5)),
            size=(12, 7), background_color=(0.216, 0.67, 0.784, 1))

        dlg_wgt = SP_MoreOpts_Dialog(ds_type=active_tab, pos=(wgt_x, wgt_y))

        self.root_window.add_widget(dlg_wgt)
        self.root_window.add_widget(point_wgt)

        self.sp_moreopts = True

    def load_proteomes(self, selection):
        """
        Similar to load method, but specific for loading proteome files.
        Given the potential size of these files, they are not stored in
        memory, but instead are processed on the fly

        :param selection: list. Contains complete paths to the proteome
        files
        """

        # Check if alignment files have already been loaded. If so, issue
        # warning, since the two data types cannot be loaded in simultaneous
        if self.file_list:
            return self.dialog_warning("Cannot mix alignment and proteome "
                                       "data",
                                       "TriFusion currently has alignment data"
                                       " loaded. If you want to load proteome"
                                       " data, please remove the current"
                                       " data set.")

        # Collect input files. Search directories, if they are provided
        original_list = []
        for f in selection:
            if os.path.isdir(f):
                original_list.extend([join(f, x) for x in os.listdir(f)
                                  if os.path.isfile(join(f, x))])
            else:
                original_list.append(f)

        # Check if files exist
        file_list = [f for f in original_list if exists(f)]

        # If input is empty, raise warning and return
        if not file_list:
            return self.dialog_floatcheck("The provided proteome files do"
                                          " not exist.", t="error")

        # Stores invalid proteome files
        bad_proteomes = {"invalid": [], "no_fasta": [], "no_protein": []}
        good_proteomes = []

        # Filter non

        for f in file_list:
            b = Base()
            er = b.autofinder(f)
            f_short = basename(f)
            if isinstance(er, Exception):
                bad_proteomes["invalid"].append(f_short)
            elif er[0] != "fasta":
                bad_proteomes["no_fasta"].append(f_short)
            elif er[1][0] != "Protein":
                bad_proteomes["no_protein"].append(f_short)
            else:
                good_proteomes.append(f)

        # If there are proteome files already loaded, extend
        if self.proteome_files:
            for f in good_proteomes:
                if f not in self.proteome_files:
                    self.proteome_files.append(f)
                    self.active_proteome_files.append(f)

        # If there are no proteome files loaded
        else:
            self.proteome_files = good_proteomes
            self.active_proteome_files = deepcopy(self.proteome_files)

        # Update gene filter value to number of proteomes. This will
        # automatically set the minimum number of species to the number of
        # proteome files, each of which should represent a species.
        self.orto_min_sp = len(self.proteome_files)

        # Issue float dialog informing that files have been loaded
        if good_proteomes:
            self.dialog_floatcheck("%s file(s) successfully loaded" %
                                   len(good_proteomes), t="info")

            if set(original_list) - set(file_list):
                missing = len(set(original_list) - set(file_list))
                self.dialog_floatcheck("%s file(s) could not be found" %
                                       missing, t="warning")

            # Update the filename - path mapping attribute
            self.filename_map = dict(list(self.filename_map.items()) +
                                     list((x, y) for x, y in
                                          zip([basename(x) for x in
                                               good_proteomes],
                                              good_proteomes)))

            # Populate file buttons in side panel
            self.original_file_inf = self.get_file_information(
                mode="proteome")
            self.populate_input_files(mode="proteome")

        # Automatically set the default orthology directory as the same
        # directory of the input files
        # Remove potential trailing dir separators that prevent
        # basename from finding the last dir name
        self.ortho_dir = self._common_path(file_list).rstrip(sep)

        if list(bad_proteomes.values()) != [[], [], []]:
            msg = ""
            if bad_proteomes["invalid"]:
                msg += "The following files are in invalid format:\n%s" \
                       "\n\n" % ", ".join(bad_proteomes["invalid"])
            if bad_proteomes["no_fasta"]:
                msg += "The following files are not in FASTA format:\n%s" \
                       "\n\n" % ", ".join(bad_proteomes["no_fasta"])
            if bad_proteomes["no_protein"]:
                msg += "The following files do not contain protein " \
                       "sequences:\n%s\n\n" \
                       % ", ".join(bad_proteomes["no_protein"])

            return self.dialog_warning("Invalid proteome files detected",
                                       msg)

    def update_tabs(self):
        """
        Wrapper that updates the contents of the files and taxa tabs
        """

        self.populate_input_files()
        self.populate_species()
        self.populate_partitions()

    def update_taxa(self):
        """
        This checks whether some taxa that were specific to some file(s)
        were removed when that file is removed.
        """

        # If taxa were removed during the update, remove those buttons too
        current_tx = [x[0].text for x in self.sp_taxa_bts]
        removed_taxa = list(set(current_tx) - set(
            self.alignment_list.taxa_names))

        if removed_taxa:
            for i in removed_taxa:
                # Get the corresponding buttons:
                x_but_txt = "%sX" % i
                bt_obj = [x for x in self.root.ids.taxa_sl.children
                          if x_but_txt == x.id][0]
                self.remove_bt(bt_obj)
                try:
                    self.active_taxa_list.remove(i)
                except ValueError:
                    pass

        #self.active_taxa_list = self.alignment_list.taxa_names

    def update_partitions(self):
        """
        Updates partition buttons following any change to input data
        """

        # Check for missing partitions based on id and remove them
        for bt, _, x_bt in self.sp_partition_bts:
            if bt.id not in self.alignment_list.partitions.partitions:
                self.remove_bt(x_bt, parent_wgt=self.root.ids.partition_sl)

    def update_partition_label(self):
        """
        Setsand updates a label on the Partitions tab of the side panel,
        informing how partitions files are selected out of the total
        partitions
        """

        if self.alignment_list:
            # Get total number of partitions
            total_parts = len(self.alignment_list.partitions.partitions)
            # Get selected partitions
            active_parts = len([x for x in
                self.root.ids.partition_sl.children if
                isinstance(x, TGToggleButton) and x.state == "down"])
            self.root.ids.partition_lab.text = "%s of %s partitions " \
                "selected" % (active_parts, total_parts)

    def update_file_label(self):
        """
        Sets and updates a label on the Files tab of the side panel,
        informing how many files are selected out of the total files
        """

        # Determine which list is used to populate
        lst = self.file_list if self.file_list else self.proteome_files
        active_lst = self.active_file_list if self.active_file_list else \
            self.active_proteome_files

        self.root.ids.file_lab.text = "%s of %s files selected" % (
            len(active_lst), len(lst))

    def update_sp_label(self):
        """
        Sets and updates a label on the Taxa tab of the side panel,
        informing how many taxa are selected out of the total taxa. If
        the taxa list is empty, it disables the select/deselect buttons
        """

        self.root.ids.sp_lab.text = "%s of %s taxa selected" % (
            len(self.active_taxa_list), len(self.alignment_list.taxa_names))

        # Add disabled no taxa button when empty
        if len(self.active_taxa_list) == 0 and \
                len(self.alignment_list.taxa_names) == 0:
            if "species_temp" not in [x.id for x in
                                      self.root.ids.taxa_sl.children]:
                no_bt = Button(id="species_temp", text="No species loaded",
                               size_hint_y=None, height=40, disabled=True)
                self.root.ids["species_temp"] = no_bt
                self.root.ids.taxa_sl.add_widget(no_bt)

    def sidepanel_create_bts(self, idx):

        # Determine state based on active_file_list
        if self.filename_map:
            # For files
            if idx in self.filename_map:
                if self.file_list:
                    state = "down" if self.filename_map[idx] in \
                        self.active_file_list else "normal"
                elif self.proteome_files:
                    state = "down" if self.filename_map[idx] in \
                                      self.active_proteome_files else "normal"
            # For taxa
            else:
                state = "down" if idx in self.active_taxa_list else "normal"
        else:
            state = "down"

        bt = TGToggleButton(text=idx, id=idx, state=state, height=30,
                            size_hint_x=.8, shorten=True,
                            shorten_from="right")

        # Setting horizontal text size for shortening
        bt.text_size[0] = bt.size[0] * 2
        # Binding functionality to toggle button
        bt.bind(on_release=self.toggle_selection)

        # Set Information button and add the widget
        inf_bt = Button(size_hint=(None, None), width=30,
                        height=30, id="%s?" % idx,
                        border=(0, 0, 0, 0),
                        background_normal=join("data", "backgrounds",
                                               "info_bt_down.png"),
                        background_down=join("data", "backgrounds",
                                             "info_bt.png"))
        inf_bt.bind(on_release=self.get_popup_info)

        # Set remove button with event binded and add the widget
        x_bt = Button(size_hint=(None, None), width=30,
                      height=30, id="%sX" % idx,
                      border=(0, 0, 0, 0),
                      background_normal=join("data", "backgrounds",
                                             "remove_bt.png"),
                      background_down=join("data", "backgrounds",
                                           "remove_bt_down.png"))
        x_bt.bind(on_release=partial(self.check_action,
                                     "Are you sure you want to remove"
                                     " this item?",
                                     self.remove_bt))

        return bt, inf_bt, x_bt

    def populate_input_files(self, mode="alignment"):
        """
        This method grabs the input files that were selected in the
        FileChooser widget and populates the File tab in the main side panel
        with toggle and remove buttons for each file

        :param mode: string. Determines which list will be used to populate.
        If alignment, it will use file_list; If proteome, it will use
        proteome_list
        """

        # Determine which list is used to populate
        if mode == "alignment":
            lst = self.file_list
        else:
            lst = self.proteome_files

        # Remove the initial disabled button, if it's still there
        if "file_temp" in self.root.ids.keys():
            self.root.ids.file_sl.remove_widget(self.root.ids.file_temp)
            del self.root.ids["file_temp"]
            self.root.ids.file_sl.height = 5

        # Add a label at the end of the file list informing how many files
        # are currently selected out of the total files
        self.update_file_label()
        self.update_partition_label()

        for infile in lst:

            if self.count_files <= self.MAX_FILE_BUTTON:

                self.count_files += 1
                file_name = basename(infile)
                self.sidepanel_add_bts(file_name, "Files")

            else:
                # Check if morebt is already present
                if not [x for x in self.root.ids.file_sl.children if
                        isinstance(x, LoadMoreBt)]:
                    self.root.ids.file_sl.add_widget(LoadMoreBt())
                return

    def sidepanel_add_bts(self, idx, tab_name, partition_fls=None):
        """
        Wrapper that adds buttons to the corresponding File, Taxa or
        Partitions panel. Also prevents duplicates from being entered and
        adds buttons for mouse over events.
        :param idx: string. Identifier of the current item. It will be the
        name of the button
        :param tab_name: string. Identifies in which side panel gridlayout
        the button will be added
        :param partition_fls: string. Only for Partitions tab_name. Includes
        full path of the alignment for the partition.
        """

        # Set attributes to be added
        if tab_name == "Files":
            grid_wgt = self.root.ids.file_sl
            bt_list = self.sp_file_bts

        elif tab_name == "Taxa":
            grid_wgt = self.root.ids.taxa_sl
            bt_list = self.sp_taxa_bts

        else:
            grid_wgt = self.root.ids.partition_sl
            bt_list = self.sp_partition_bts

        # This prevents duplicate entries from being added
        if idx not in [x.id for x in grid_wgt.children]:

            # Create buttons
            if tab_name == "Partitions":
                bt, inf_bt, x_bt = self.sidepanel_create_part_bts(
                    [idx, partition_fls])
            else:
                bt, inf_bt, x_bt = self.sidepanel_create_bts(idx)

            # Add button to storage for mouse over events
            self.mouse_over_bts[tab_name].append(bt)

            # Adds buttons to gridlayout
            grid_wgt.add_widget(bt)
            grid_wgt.add_widget(inf_bt)
            grid_wgt.add_widget(x_bt)

            # Add all three buttons to the storage attribute
            bt_list.append((bt, inf_bt, x_bt))

    def sidepanel_load_more_filebts(self):

        max_buttons = self.MAX_FILE_BUTTON + self.count_files

        self.root.ids.file_sl.remove_widget(
            self.root.ids.file_sl.children[0])

        for _ in range(self.count_files, max_buttons):

            self.count_files += 1

            try:
                infile = basename(self.file_list[self.count_files])
                self.sidepanel_add_bts(infile, "Files")
            except IndexError:
                return

        self.root.ids.file_sl.add_widget(LoadMoreBt())

    def populate_species(self):
        """
        This method grabs the taxa names from the input files that were
        selected in the FileChooser widget and populates the Taxa tab in the
        main side panel with toggle and remove buttons for each taxon
        """

        # Remove the initial disabled button if it's still there
        if "species_temp" in self.root.ids.keys():
            self.root.ids.taxa_sl.remove_widget(self.root.ids.species_temp)
            del self.root.ids["species_temp"]
            self.root.ids.taxa_sl.height = 5

        # Add a label at the end of the taxa list informing how many taxa
        # are currently selected out of the total taxa
        self.update_sp_label()

        for tx in sorted(self.active_taxa_list):

            # Prevents duplicate taxa from being entered
            if tx not in [x.id for x in self.root.ids.taxa_sl.children]:

                self.sidepanel_add_bts(tx, "Taxa")

    def repopulate_partitions(self):
        """
        Wrapper method that re-populates partitions after changes
        """

        # Clear partitions
        self.root.ids.partition_sl.clear_widgets()
        self.count_partitions = 0
        # Re-populate partitions
        self.populate_partitions()

    def populate_partitions(self):
        """
        Populates the partitions tab in the side bar from the partitions
        object associated with alignment objects.

        This method is used when input files are loaded into the program,
        which means there will be no issue with multiple files being
        associated with the same partitions. This kind of change is done
        a posteriori when importing partition files or setting the
        partitions manually.
        """

        # Remove initial disabled button, if it's still there
        if "partition_temp" in self.root.ids.keys():
            try:
                self.root.ids.partition_sl.remove_widget(
                    self.root.ids.partition_temp)
                del self.root.ids["partition_temp"]
            except ReferenceError:
                pass

        for partition, fls in self.alignment_list.partitions.iter_files():

            if self.count_partitions <= self.MAX_PARTITION_BUTTON:

                self.count_partitions += 1

                # Create partition buttons
                self.sidepanel_add_bts(partition, "Partitions",
                                       partition_fls=fls)

            else:
                self.root.ids.partition_sl.add_widget(LoadMoreBt())
                return

    def partitions_merge_dialog(self):
        """
        Dialog that appears when clicking merge partitions asking for the
        name of the new partition
        """

        content = InputTextDialog(cancel=self.dismiss_popup,
                                  action=self.partitions_merge)

        self.show_popup(title="Choose name for new partition",
                        content=content,
                        size=(300, 153.5))

    def partitions_merge(self, name):
        """
        Merge active partitions

        :param name: string. Name of the new partition
        """

        if name in self.alignment_list.partitions.partitions:
            return self.dialog_floatcheck(
                "A partition named %s already exists." % name,
                t="error")

        self.alignment_list.partitions.merge_partitions(
            self.active_partitions, name)

        # Resets active partitions
        self.active_partitions = []

        # Close popup
        self.dismiss_popup()
        self.repopulate_partitions()

    def partitions_split(self, new_range=None, new_names=None):
        """
        Split an active partition

        :param new_range: Optional. Tuple containing the new ranges for the
        new partitions
        :param new_names: Optional. Tuple containing the new names for each
        new partition
        """

        # Check whether any of the new partition names already exist. Issue
        # an error if so.
        if new_names:
            if [x for x in new_names if x in
                    self.alignment_list.partitions.partitions]:
                return self.dialog_floatcheck("Duplicate partition names "
                                              "deteced: {}".format(
                    " ".join([x for x in new_names if x if
                    self.alignment_list.partitions.partitions])), t="error")

        # Get active partition
        active_partition = [x.text for x in
                            self.root.ids.partition_sl.children
                            if isinstance(x, ToggleButton) and
                            x.state == "down"][0]

        if new_range:
            self.alignment_list.partitions.split_partition(
                active_partition, new_range, new_names)

        else:
            self.alignment_list.partitions.split_partition(active_partition)

        # Resets ative partitions
        self.active_partitions = []

        # Close popup
        self.dismiss_popup()
        self.repopulate_partitions()

    def dialog_partitions_split(self):

        content = SplitPartitions(cancel=self.dismiss_popup)

        # Get active partition
        active_partition = [x.text for x in
                            self.root.ids.partition_sl.children
                            if isinstance(x, ToggleButton) and
                            x.state == "down"][0]

        # Get partition range
        part_range = self.alignment_list.partitions.\
            partitions[active_partition]

        # Disable manual split if the current partition has a fragmented
        # range
        if not isinstance(part_range[0], tuple):
            for wgt in [content.ids.ok_bt, content.ids.part1,
                        content.ids.part2]:
                wgt.disabled = True
            content.ids.manual_slider.max = 0
        else:
            content.ids.manual_slider.min = int(part_range[0][0]) + 1
            content.ids.manual_slider.max = int(part_range[0][1]) + 1

        # If partition contains only one file, disable automatic split by
        #  files
        if len(self.alignment_list.partitions.partitions_alignments[
                active_partition]) == 1:
            content.ids.auto_split_bt.disabled = True

        self.show_popup(title="Split partition", content=content,
                        size=(420, 320))

        if not isinstance(part_range[0], tuple):
            self.dialog_floatcheck(
                "Manual split in unavailable for partitions with "
                "fragmented ranges.",
                t="warning")

    def partitions_change_name(self, partition_name, new_name):
        """
        Changes name of a partition
        :param partition_name: string, Original partition name
        :param new_name: string, new partition name
        :return:
        """

        # Change partition name
        self.alignment_list.partitions.change_name(partition_name, new_name)

        # Update button in side panel
        bt = [x for x in self.root.ids.partition_sl.children
              if isinstance(x, ToggleButton) and
              x.text == partition_name][0]
        bt.text = new_name

    def partitions_import_scheme(self, partition_file):
        """
        Imports partitions in partition_file and applies to the
        alignment_list partition object. It applies to the entire
        alignment_list for now.
        :param partition_file: string, path to partition file
        """

        self.alignment_list.partitions.read_from_file(partition_file)

        # Clear partitions
        self.root.ids.partition_sl.clear_widgets()
        #C Clear partition button storage
        self.sp_partition_bts = []
        # Re-populate partitions
        self.populate_partitions()
        # Update partitions label
        self.update_partition_label()

        self.dismiss_popup()

    def sidepanel_create_part_bts(self, idx):
        """
        Creates buttons for each partition
        :param idx: string. unique identifier of partition
        """

        part_name = idx[0]
        fl_num = str(len(idx[1]))

        # Create main button
        bt = TGToggleButton(text=part_name, state="normal", id=part_name,
                            size_hint_y=.8, shorten=True, height=30,
                            shorten_from="right")

        bt.bind(on_release=self.toggle_selection)

        # Setting horizontal text size for shortening
        bt.text_size[0] = bt.size[0] * 2

        # Create file counter button. This button will display the number
        # of alignments included in this partition as its text. The
        # on_release event will show a popup with a list of the alignment
        # files contained
        c_bt = Button(size_hint=(None, None), width=30, text=fl_num,
                      height=30, id="%sC" % part_name, bold=True,
                      border=(0, 0, 0, 0))
        c_bt.background_normal = join("data", "backgrounds",
                                      "bt_process.png")
        c_bt.background_down = join("data", "backgrounds", "count_bt.png")
        c_bt.bind(on_release=lambda x: self.dialog_partition_files(bt.text))

        # Create edition button
        ed_bt = Button(size_hint=(None, None), width=30,
                       height=30, id="%s?" % part_name, border=(0, 0, 0, 0))
        ed_bt.background_normal = join("data", "backgrounds",
                                       "edit_bt_down.png")
        ed_bt.background_down = join("data", "backgrounds", "edit_bt.png")
        ed_bt.bind(on_release=self.dialog_partitions)

        return bt, c_bt, ed_bt

    def set_codon_model(self, codon_partition, wgt=None):
        """
        Changes the model spinners when changing the codon partitioning

        :param codon_partition: string. The codon partition string that
        should correspond to a key in partition_model
        :param wgt: Widget object containing the partitions dialog
        """

        first_background = join("data", "backgrounds", "model_bt1.png")
        second_background = join("data", "backgrounds", "model_bt2.png")
        third_background = join("data", "backgrounds", "model_bt3.png")

        partition_model = {
            "[color=ff5555ff]1[/color] + [color=37abc8ff]2[/color] + "
            "[color=71c837ff]3[/color]":
            [ModelSpinner(background_normal=first_background, id="1"),
            ModelSpinner(background_normal=second_background, id="2"),
            ModelSpinner(background_normal=third_background, id="3")],
            "[color=ff5555ff](1 + 2)[/color] + [color=37abc8ff]3[/color]":
            [ModelSpinner(background_normal=first_background, id="12"),
            ModelSpinner(background_normal=second_background, id="3")],
            "[color=ff5555ff]1[/color] + [color=37abc8ff](2 + 3)[/color]":
            [ModelSpinner(background_normal=first_background, id="1"),
            ModelSpinner(background_normal=second_background, id="23")],
            "[color=ff5555ff](1 + 3)[/color] + [color=37abc8ff]2[/color]":
            [ModelSpinner(background_normal=first_background, id="13"),
            ModelSpinner(background_normal=second_background, id="2")],
            "No partitions": [ModelSpinner(id="0")]}

        if wgt:
            partitions_wgt = wgt
        else:
            partitions_wgt = [x for x in self.root_window.children if
                              isinstance(x, PartitionsDialog)][0]

        partitions_wgt.ids.model_bx.clear_widgets()

        for model in partition_model[codon_partition]:
            partitions_wgt.ids.model_bx.add_widget(model)

    def save_model(self, part_name, partition_wgt, apply_all=False):
        """
        Saves the model currently set in the partitions dialog.
        :param part_name: string, name of partition
        :param partition_wgt: Widget of the Partitions dialog
        :param apply_all: boolean, whether the current model will be
        applied to all partitions or not
        """

        model_spiners = [x for x in partition_wgt.ids.model_bx.children]

        if len(model_spiners) == 1:
            self.alignment_list.partitions.set_model(
                part_name, [model_spiners[0].text], apply_all=apply_all)
        else:
            models = []
            links = []
            for wgt in model_spiners:
                models.extend([wgt.text] * len(wgt.id))
                links.append(wgt.id)

            links, models = [list(x) for x in zip(
                *sorted(zip(links, models), key=lambda pair: pair[0]))]
            self.alignment_list.partitions.set_model(
                part_name, models, links, apply_all=apply_all)

    def remove_partition_box(self):
        """
        Removes a currently active partition box when clicking the X button
        """

        # Gathers active partition box and remove button
        partition_box = [x for x in self.root_window.children if
                         isinstance(x, PartitionsDialog) or
                         isinstance(x, RemoveFloat)]

        # Removes widgets:
        for wgt in partition_box:
            self.root_window.remove_widget(wgt)

    def dialog_partition_files(self, partition_name):
        """
        Shows a popup listing the files associated with a given partition
        :param partition_name: string, name of partition
        :return:
        """

        content = BtList(cancel=self.dismiss_popup)

        plist = self.alignment_list.partitions.partitions_alignments[
            partition_name]

        for f in plist:
            bt = TFButton(text=basename(f), height=40)
            content.ids.rev_inlist.add_widget(bt)

        self.show_popup(title="Alignment files in %s" % partition_name,
                        content=content, size_hint=(.6, .8))

    def dialog_select_taxa_group(self):
        """
        Shows a subpopup listing the taxa groups that have already been
        created. Each taxa group will be a Button widget
        """

        content = BtList(cancel=self.dismiss_subpopup)

        if self.taxa_groups:
            for nm in self.taxa_groups:
                bt = TFButton(text=nm, size_hint_y=None, height=30)
                bt.bind(on_release=self.select_taxa_group)
                content.ids.rev_inlist.add_widget(bt)

        else:
            bt = TFButton(text="No groups defined", disabled=True)
            content.ids.rev_inlist.add_widget(bt)

        self.show_popup(title="Select taxa group", content=content,
                        size_hint=(.4, .8), popup_level=2)

    def select_taxa_group(self, bt):
        """
        Gives functionality to the buttons in the dialog_select_taxa_group.
        Saves the taxa group name and closes subpopup

        :param bt: Button object
        """

        if self.taxa_filter_settings:
            self.taxa_filter_settings[1] = bt.text
        else:
            self.taxa_filter_settings = ["Contain", bt.text]

        self.dismiss_subpopup()

    def dialog_create_group_from_file(self, ds_type):
        """
        Creates a filechooser dialog  to select a file containing a
        taxa/file list that will be used to generate a data set group
        :param ds_type: string. Identifies the data set type. Either taxa or
        files
        """

        content = SaveDialog(cancel=self.dismiss_popup,
                             bookmark_init=self.bookmark_init)

        content.ids.txt_box.clear_widgets()
        content.ids.txt_box.height = 0
        title = "Choose text file to import"
        content.ids.sd_filechooser.text = ds_type

        self.show_popup(title=title, content=content)

    def dialog_select_from_file(self):
        """
        Creates a filechooser dialog to select a text file containing a list
        of files or taxa names to be selected in the side panel

        Calls: select_bt_from_file
        """

        content = SaveDialog(cancel=self.dismiss_popup,
                             bookmark_init=self.bookmark_init)

        content.ids.txt_box.clear_widgets()
        content.ids.txt_box.height = 0
        title = "Choose text file to import"
        content.ids.sd_filechooser.text = "select_from_file"

        self.show_popup(title=title, content=content)

    def dialog_remove_from_file(self):
        """
        Creates a filechooser dialog to select a text file containing a list
        of files or taxa names to be removed in the side panel.

        Then calls: self.remove_bt_from_file
        """

        content = SaveDialog(cancel=self.dismiss_popup,
                             bookmark_init=self.bookmark_init)

        content.ids.txt_box.clear_widgets()
        content.ids.txt_box.height = 0
        title = "Choose text file to import"
        content.ids.sd_filechooser.text = "remove_from_file"

        self.show_popup(title=title, content=content)

    def dialog_import_partitions(self):
        """
        Creates a filechooser dialog to select a partition file and import
        its scheme to the current partition. If one or more partitions
        are active, ask the user if we wants to import the partition
        scheme to the selected partitions or to the whole dataset
        :return:
        """

        content = SaveDialog(cancel=self.dismiss_popup,
                             bookmark_init=self.bookmark_init)

        content.ids.txt_box.clear_widgets()
        content.ids.txt_box.height = 0
        title = "Choose partition scheme file"
        content.ids.sd_filechooser.text = "import_partitions"

        self.show_popup(title=title, content=content)

    def dialog_partitions(self, btx):
        """
        Shows a small widget with partition information

        :param btx: Button widget that will be used to determine which
        partition is being viewed and the position of the partitions dialog
        """

        def flatter(s):
            """
            Creates a flat iterator of tuples. If s is [[(1,2), (2,3)],
            (4,5)] this will yield ((1,2), (2,3), (4,5))

            :param s: list.
            """
            for y in s:
                if isinstance(y, tuple):
                    yield y
                else:
                    for j in y:
                        yield j

        partition_model = {
            "1,2,3":
                "[color=ff5555ff]1[/color] + [color=37abc8ff]2[/color]"
                " + [color=71c837ff]3[/color]",
            "12,3":
                "[color=ff5555ff](1 + 2)[/color] + [color=37abc8ff]3"
                "[/color]",
            "1,23":
                "[color=ff5555ff]1[/color] + [color=37abc8ff](2 + 3)"
                "[/color]",
            "13,2":
                "[color=ff5555ff](1 + 3)[/color] + [color=37abc8ff]2"
                "[/color]"}

        # Get position of partition edit button:
        ed_pos = btx.to_window(btx.pos[0], btx.pos[1])

        # Set position for partitions dialog
        size = (240, 260)
        pos = [ed_pos[0] + btx.width,
               ed_pos[1] + (btx.height / 2) - (size[1] / 2)]

        content = PartitionsDialog(pos=pos, size=size,
                                   size_hint=(None, None))
        rm_wgt = RemoveFloat(pos=[pos[0] + size[0] - 20,
                                  pos[1] + size[1] - 20])

        # Set partition object and partition name

        # Since partition names can be changed and I can only get the
        # partition name from he edition button id (which does not
        # change), this iteration over all three partition buttons for
        # each partition will retrieve the correct partition name
        displayed_partitions = (
            x for x in self.root.ids.partition_sl.children if not
            isinstance(x, LoadMoreBt))
        part_name = [bt.text for ebt, _, bt in
                     zip(*[iter(displayed_partitions)] * 3)
                     if ebt.id == btx.id][0]
        part_obj = self.alignment_list.partitions
        content.ids.partition_name.text = part_name
        content.original_name = part_name

        # Get partition length
        part_range = (y[0] for x, y in self.alignment_list.partitions
                      if x == part_name)
        part_len = sum([x[1] - x[0] + 1 for x in flatter(part_range)])
        content.ids.partition_lenght.text = "{}bp".format(part_len)

        # If there are codon partitions
        if part_obj.partitions[part_name][1]:
            if not part_obj.models[part_name][2]:
                content.ids.codon_spin.text = content.ids.codon_spin.\
                    values[1]
                self.set_codon_model(
                    content.ids.codon_spin.values[1], content)
            else:
                m_key = ",".join(part_obj.models[part_name][2])
                content.ids.codon_spin.text = partition_model[m_key]
                self.set_codon_model(partition_model[m_key], content)
            if part_obj.models[part_name][0][0]:
                for i in range(len(part_obj.models[part_name][0])):
                    params = part_obj.models[part_name][0][i]
                    model = part_obj.get_model_name(params)
                    content.ids.model_bx.children[i].text = model
            else:
                for p, m in enumerate(part_obj.models[part_name][1]):
                    if m:
                        content.ids.model_bx.children[::-1][p].text = m
        elif part_obj.models[part_name][0][0]:
            params = part_obj.models[part_name][0][0]
            model = part_obj.get_model_name(params)
            content.ids.model_bx.children[0].text = model
        elif part_obj.models[part_name][1][0]:
            model = part_obj.models[part_name][1][0]
            content.ids.model_bx.children[0].text = model

        # Give functionality to remove button
        rm_wgt.bind(on_release=lambda y: self.remove_partition_box())

        self.root_window.add_widget(content)
        self.root_window.add_widget(rm_wgt)

    def get_popup_info(self, value):
        """
        This function determines whether the taxa information has changed
        since it has been update (when files are removed). If so,
        it updates the taxa information and issues popup_info. Else,
        just issues popup_info
        :param value: Button object that issued the function
        """

        if self.trigger_taxa_update:
            self.trigger_taxa_update = False
            self.run_in_background(self.get_taxa_information,
                                   self.popup_info,
                                   None,
                                   [value],
                                   cancel=False,
                                   msg="Updating taxa information")
        else:
            self.popup_info(value)

    def popup_info(self, value):
        """
        Generates the pop up information content for the pressed taxa or
        file button
        :param value: the button object is provided when binding
        """

        # Determining if the request comes from file or taxa tab
        if value.parent == self.root.ids.taxa_sl:

            # Get taxa name
            tx = value.id[:-1]

            if tx in self.active_taxa_list:

                # Get the information from the content list. This is done
                # when calling the popup to avoid repeating this
                # operation every time taxa  or files are added/removed.
                self.active_tx_inf = self.get_taxon_information(
                    tx, self.alignment_list)

                content = BoxLayout(orientation="vertical", padding=10,
                                    spacing=10)
                sv = ScrollView(scroll_type=["bars"], bar_width=10)
                all_ds = BoxLayout(orientation="vertical",
                                   height=2 * (30 * 7) + 10,
                                   size_hint_y=None)
                total_ds = TaxaPopup(height=30 * 7)
                total_ds.ids.dataset_label.text = "Complete data set"
                active_ds = TaxaPopup(height=30 * 7)
                active_ds.ids.dataset_label.text = "Active data set"

                # Populate complete data set contents
                total_ds.ids.seq_len.text = "%s" % \
                    self.original_tx_inf[tx]["length"]
                total_ds.ids.indels.text = "%s" % \
                    self.original_tx_inf[tx]["indel"]
                total_ds.ids.missing.text = "%s" % \
                    self.original_tx_inf[tx]["missing"]
                total_ds.ids.ef_seq_len.text = ("%s (%s%%)" % (
                    self.original_tx_inf[tx]["effective_len"],
                    self.original_tx_inf[tx]["effective_len_per"]))
                total_ds.ids.file_cov.text = ("%s (%s%%)" % (
                    self.original_tx_inf[tx]["fl_coverage"],
                    self.original_tx_inf[tx]["fl_coverage_per"]))

                # Populate active data set contents
                active_ds.ids.seq_len.text = "%s" % \
                    self.active_tx_inf["length"]
                active_ds.ids.indels.text = "%s" % \
                    self.active_tx_inf["indel"]
                active_ds.ids.missing.text = "%s" % \
                    self.active_tx_inf["missing"]
                active_ds.ids.ef_seq_len.text = ("%s (%s%%)" % (
                    self.active_tx_inf["effective_len"],
                    self.active_tx_inf["effective_len_per"]))
                active_ds.ids.file_cov.text = ("%s (%s%%)" % (
                    self.active_tx_inf["fl_coverage"],
                    self.active_tx_inf["fl_coverage_per"]))

                close_bl = CloseBox(cancel=self.dismiss_popup)

                all_ds.add_widget(total_ds)
                all_ds.add_widget(active_ds)
                sv.add_widget(all_ds)
                content.add_widget(sv)
                content.add_widget(close_bl)

                self.show_popup(title="Taxon: %s" % value.id[:-1],
                                content=content, size=(450, 400))

        elif value.parent == self.root.ids.file_sl:

            # Get file name
            file_name = value.id[:-1]

            if self.filename_map[file_name] in self.active_file_list:

                content = FilePopup(cancel=self.dismiss_popup)

                # Get the information from the content list. This is done
                # when calling the popup to avoid repeating this
                # operation every time taxa  or files are added/removed.
                self.active_file_inf = self.get_file_information(file_name)

                content.ids.in_format.text = "%s" % \
                    self.active_file_inf["aln_format"]
                content.ids.seq_type.text = "%s" % \
                    self.active_file_inf["seq_type"]
                content.ids.seq_size.text = "%s" % \
                    self.active_file_inf["aln_len"]
                content.ids.n_taxa.text = "%s" % \
                    self.active_file_inf["n_taxa"]

                self.show_popup(title="File: %s" % value.id[:-1],
                                content=content, size=(400, 320))

            if self.filename_map[file_name] in self.proteome_files:

                content = ProteomePopup(cancel=self.dismiss_popup)

                content.ids.n_seq.text = "%s" % \
                    self.original_file_inf[file_name]["n_seq"]
                content.ids.n_res.text = "%s" % \
                    self.original_file_inf[file_name]["n_res"]

                self.show_popup(title="File %s" % value.id[:-1],
                                content=content, size=(400, 200))

    def change_taxa_name(self, old_name, new_name):
        """
        Changes the taxa name on a double tap
        :param old_name: string, original taxon name
        :param new_name: string, new taxon name
        """

        self.alignment_list.change_taxon_name(old_name, new_name)

        # Change active taxa list
        self.active_taxa_list = [new_name if x == old_name else x for x
                                 in self.active_taxa_list]

        # Change sp_taxa_bts attribute
        for taxa, inf, rm in self.sp_taxa_bts:
            if taxa.text == old_name:
                taxa.text = new_name
                taxa.id = new_name
                inf.id = "{}?".format(new_name)
                rm.id = "{}X".format(new_name)

        # Change tx_info attributes
        self.original_tx_inf = dict((new_name, y) if x == old_name else
            (x, y) for x, y in self.original_tx_inf.items())

        # Change mouser_over_bts attribute
        for taxa in [x for x in self.mouse_over_bts["Taxa"]
                     if isinstance(x, ToggleButton)]:
            if taxa.text == old_name:
                taxa.text = new_name

        self.sidepanel_clear_search("taxa")

    def export_names(self, path, file_name):
        """
        Export the names of buttons in the corresponding tab in the side
        panel It listens to the self.export_mode attribute, which is a
        tuple object with the first element being either "file" or "taxa"
        and the second element as "all" or "selected".

        :param path: string. Path to the output file.
        :param file_name. Name of the output file.
        """

        # Create file object
        export_file = open(join(path, file_name) + ".txt", "w")

        if self.export_mode[0] == "files":
            # Export all files
            if self.export_mode[1] == "all":
                fl = self.file_list if self.file_list else self.filename_map
                for x in fl:
                    short_name = basename(x)
                    export_file.write(short_name + "\n")
            # Export selected files
            elif self.export_mode[1] == "selected":
                fl = self.active_file_list if self.active_file_list else \
                    self.active_proteome_files
                for x in fl:
                    short_name = basename(x)
                    export_file.write(short_name + "\n")

        elif self.export_mode[0] == "taxa":
            # Export all taxa
            if self.export_mode[1] == "all":
                for x in self.alignment_list.taxa_names:
                    export_file.write(x + "\n")
            # Export selected taxa
            elif self.export_mode[1] == "selected":
                for x in self.active_taxa_list:
                    export_file.write(x + "\n")

        # Close file handle
        export_file.close()

        # Issue success message
        self.dialog_floatcheck("Successfully exported {} {}.".format(
            self.export_mode[1], self.export_mode[0]),
                               t="info")

    def toggle_multi_selection(self, value):
        """
        Adds multiple selection using shift+click for general purpose lists
        of toggle buttons in data set creator
        """

        def get_selection(bt, bt_list):

            if self.is_shift_pressed:

                if not self.prev_tb:
                    return [bt]

                # For possible exceptions where prev_tb is not in button
                # list, return only the selected button and reset the prev_tb
                # attribute
                elif self.prev_tb not in bt_list:
                    self.prev_tb = None
                    return [bt]

                start = bt_list.index(self.prev_tb)
                stop = bt_list.index(bt)

                if start < stop:
                    return bt_list[start:stop + 1]
                else:
                    return bt_list[stop:start + 1]

            else:
                return [bt]

        # Get the parent layout object
        parent_obj = value.parent

        bt_list = parent_obj.children[::-1]

        sel = get_selection(value, bt_list)

        for bt in sel:

            if value.state == "normal":
                bt.state = "normal"

            elif value.state == "down":
                bt.state = "down"

        # Update labels for data set creator dialog
        # Get dialog
        if isinstance(self._popup.content, TaxaGroupDialog):
            dlg = self._popup.content
        else:
            dlg = self._subpopup.content

        if parent_obj.text == "all":
            lbl = dlg.ids.all_lbl
        else:
            lbl = dlg.ids.select_lbl

        lbl.text = "{} of {} selected".format(len([x for x in bt_list if
                                                   x.state == "down"]),
                                              len(bt_list))

        self.prev_tb = value

    def sidepanel_invert_selection(self, panel):
        """
        
        Parameters
        ----------
        panel

        Returns
        -------

        """

        # Set gridlayout and taxa/file list of the appropraite panel
        if panel == "Files":
            grid_wgt = self.root.ids.file_sl
            lst = self.file_list if self.file_list else self.proteome_files
            active_lst = self.active_file_list if self.file_list \
                else self.active_proteome_files

        else:
            grid_wgt = self.root.ids.taxa_sl
            lst = self.alignment_list.taxa_names
            active_lst = self.active_taxa_list

        # Storage of active itens
        active_items = []

        # Invert selection and get the active elements
        for item in lst:

            if item not in active_lst:
                active_items.append(item)

        # Update panel with new active elements
        for wgt in grid_wgt.children:

            # Get only toggle buttons
            if isinstance(wgt, ToggleButton):

                # Get text. File togglebuttons have the basename of the
                # active list
                txt = self.filename_map[wgt.text] if panel == "Files" \
                    else wgt.text

                if txt in active_items:
                    wgt.state = "down"
                
                else:
                    wgt.state = "normal"

        # Core changes to taxa
        if panel == "Taxa":

            self.active_taxa_list = active_items
            self.update_sp_label()

        elif panel == "Files":

            # Update when alignment files are loaded
            if self.file_list:

                # update active file list
                self.active_file_list = active_items

                # Update AlignmentList object
                self.alignment_list.update_active_alignments(
                    active_items
                )
                
            # Update when proteome files are loaded
            if self.proteome_files:

                self.active_proteome_files = active_items

            self.update_file_label()

    def toggle_selection(self, value):
        """
        Adds functionality for the file and taxa toggle buttons in the side
        panel. It adds or removes the selected taxa from the active lists

        :param value: Button widget.
        """

        def get_selection(bt, tab):
            """
            To support multiple selection using shift+click,
            this function will return the list of buttons to be modified,
            even if it is a single one (when shift is not being pressed)
            :param bt: ToggleButton object from sidepanel
            :param tab: string, identifier of current tab. Either File,
            Taxa or Partitions
            """

            if self.is_shift_pressed:

                # If there is no previous pressed button in the current
                # the behaviour is as if shift is not being pressed
                if not self.last_sp_bt[tab]:
                    return [bt]

                start = self.mouse_over_bts[tab].index(self.last_sp_bt[tab])
                stop = self.mouse_over_bts[tab].index(bt)

                if start < stop:
                    return self.mouse_over_bts[tab][start + 1:stop + 1]
                else:
                    return self.mouse_over_bts[tab][stop:start]

            else:
                return [bt]

        # Get the parent layout object
        parent_obj = value.parent

        if self.touch.is_double_tap and parent_obj == self.root.ids.taxa_sl:
            self.dialog_text(
                "Change taxon name", "change_taxon", value.text)

        # determine active file list
        act_lst = self.active_file_list if self.file_list else \
            self.active_proteome_files

        # Changes concerning the files tab
        if parent_obj == self.root.ids.file_sl:

            sel = get_selection(value, "Files")

            for b in sel:

                # When button is normal (unselected) remove from active list
                if value.state == "normal":

                    try:

                        if self.active_file_list:
                            aln_path = self.filename_map[b.id]
                            self.alignment_list.update_active_alignment(
                                aln_path, "shelve")

                        act_lst.remove(self.filename_map[b.id])
                        b.state = "normal"

                    except ValueError:
                        pass

                # When button is down (selected) add to active list
                elif value.state == "down":
                    if self.filename_map[b.id] not in act_lst:

                        act_lst.append(self.filename_map[b.id])
                        b.state = "down"

                        if self.active_file_list:
                            aln_path = self.filename_map[b.id]
                            self.alignment_list.update_active_alignment(
                                aln_path, "active")

            # Update label
            self.update_file_label()

            # Update last pressed button
            self.last_sp_bt["Files"] = value

        # Changes concerning the taxa tab
        elif parent_obj == self.root.ids.taxa_sl:

            sel = get_selection(value, "Taxa")

            for b in sel:
                # When button is normal (unselected) remove from active list
                if value.state == "normal":
                    try:
                        self.active_taxa_list.remove(b.text)
                        b.state = "normal"
                    except ValueError:
                        pass
                # When button is down (selected) add to active
                elif value.state == "down":
                    if b.text not in self.active_taxa_list:
                        self.active_taxa_list.append(b.text)
                        b.state = "down"

            # Update label
            self.update_sp_label()

            # Update last pressed button
            self.last_sp_bt["Taxa"] = value

        elif parent_obj == self.root.ids.partition_sl:

            sel = get_selection(value, "Partitions")

            for b in sel:
                if value.state == "normal":
                    try:
                        self.active_partitions.remove(b.text)
                        b.state = "normal"
                    except ValueError:
                        pass
                else:
                    if b.text not in self.active_partitions:
                        self.active_partitions.append(b.text)
                        b.state = "down"

            # Update last pressed button
            self.last_sp_bt["Partitions"] = value

            self.update_partition_label()

    def remove_all(self):
        """
        Functionality for the remove all button for taxa and file buttons
        in the side panel. This method will remove all files and taxa
        from the program
        """

        # App changes
        # Clear widgets in side panel
        for panel in [self.root.ids.file_sl, self.root.ids.taxa_sl,
                      self.root.ids.partition_sl]:
            panel.clear_widgets()

        self.root.ids.sp_lab.text = ""
        self.root.ids.file_lab.text = ""
        self.root.ids.partition_lab.text = ""
        self.count_files = 0
        self.count_partitions = 0

        self.clear_process_input()
        self.clear_orto_input()

        # Add disabled no changes button
        if "file_temp" not in [x.id for x in
                               self.root.ids.file_sl.children]:
            no_bt = Button(id="file_temp", text="No files loaded",
                           size_hint_y=None, height=40, disabled=True)
            self.root.ids["file_temp"] = no_bt
            self.root.ids.file_sl.add_widget(no_bt)

        # Add disabled no changes button
        if "species_temp" not in [x.id for x in
                                  self.root.ids.taxa_sl.children]:
            no_bt = Button(id="species_temp", text="No files loaded",
                           size_hint_y=None, height=40, disabled=True)
            self.root.ids["species_temp"] = no_bt
            self.root.ids.taxa_sl.add_widget(no_bt)

    def clear_orto_input(self):
        """
        Clears any input for the orthology screen and related variables and
        attributes
        """

        self.proteome_files = []
        self.active_proteome_files = []
        self.orto_min_sp = 3

    def clear_process_input(self):
        """
        Clears any input for the process/statistics screen and related
        variables and attributes
        """

        # Stop other background processes
        self.terminate_stats = True
        self.terminate_background = True
        self.terminate_load_files = True

        # Give time for processes to end
        time.sleep(.1)

        self.alignment_list.clear_alignments()
        self.original_tx_inf.clear()
        self.active_tx_inf.clear()
        self.original_file_inf.clear()
        self.active_file_inf.clear()
        self.active_taxa_list = []
        self.filename_map.clear()
        self.file_list = []
        self.active_file_list = []
        self.sequence_types = ""
        self.MAX_FILE_BUTTON = 20
        self.MAX_PARTITION_BUTTON = 20
        self.mouse_over_bts.clear()
        self.mouse_over_bts = {"Files": [], "Taxa": [], "Partitions": []}
        self.previous_sets = {"Files": [], "Taxa": [], "Stats": []}
        self.sp_taxa_bts = []
        self.sp_file_bts = []
        self.sp_partition_bts = []
        self.previous_mouse_over = ""
        self.current_plt_idx = []
        self.current_plot = None
        self.current_table = None
        self.output_dir = ""
        self.output_file = ""
        self.conversion_suffix = ""
        self.gene_table_selection = []
        self.gene_master_table = []

        # Clear Statistics screen scatter, if screen is active
        self.dismiss_stats_toggle()
        self.previous_stats_toggle = None
        if self.screen.name == "Statistics":
            try:
                if not isinstance(self.screen.ids.plot_content.children[
                                  0].children[0], NoDataLabel):
                    self.screen.ids.plot_content.clear_widgets()
                    self.screen.ids.plot_content.add_widget(NoDataLabel())
            except IndexError:
                pass
            self.screen.ids.gene_num.text = \
                "Genes: [color=37abc8ff]N/A[/color]"
            self.screen.ids.taxa_num.text = \
                "Taxa: [color=37abc8ff]N/A[/color]"
        else:
            self.loaded_screens[self.available_screens[3]] = None

    def remove_all_groups(self):
        """
        Removes all loaded orthology groups
        """

        # Clear gridlayout contents
        for gl in self.screen.ids.orto_group_glbx.children:
            gl.clear_widgets()

        # Clear orthology cards
        self.screen.ids.card_gl.clear_widgets()

        # Resets ortho_groups object
        # Formal clearing of ortho_groups and temporary files
        self.ortho_groups.clear_groups()
        # Resets the ortho_groups variable so that new loadings will start anew
        self.ortho_groups = None
        self.ortho_group_files = []
        #self.proteome_files = []

    def remove_groups(self, value):
        """
        Removes orthology group buttons
        :param value: Instance of remove button
        """

        # Remove group from MultiGroup object
        self.ortho_groups.remove_group(value.id)

        # Remove from ortho group file container
        self.ortho_group_files.remove(value.id)

        # If all groups have been removed, reset ortho_groups attribute
        if not self.ortho_group_files:
            self.ortho_groups = None

        # Get box container of all gridlayouts
        gl_bx = value.parent.parent

        for gl in gl_bx.children:
            # Remove appropriate item, according to id, from its gridlayout
            gl.remove_widget([x for x in gl.children if
                              x.id == value.id][0])

        # If no group button is active, dispatch the first
        if not [x for x in self.screen.ids.group_gl.children
                if x.state == "down"] and self.screen.ids.group_gl.children:
            self.screen.ids.group_gl.children[-1].dispatch("on_release")
            self.screen.ids.group_gl.children[-1].state = "down"

        if not self.screen.ids.group_gl.children:
            self.screen.ids.card_gl.clear_widgets()

    def remove_bt(self, value, parent_wgt=None):
        """
        Functionality for the "X" remove buttons in the side panel. It
        removes button pairs with similar id's and can be used in both files
        and taxa tabs
        :param value: Button widget to be removed
        :param parent_wgt: Button widget contained. If provided, it
        overrides value.parent
        """

        # APP CHANGES
        # Get the parent layout object from where the widget will be removed
        if parent_wgt:
            parent_obj = parent_wgt
        else:
            parent_obj = value.parent

        # If the button is the last file or taxa element, issue a remove_all
        if parent_obj == self.root.ids.file_sl:
            if len(parent_obj.children) == 3 and len(self.file_list) == 1:
                return self.remove_all()
        elif parent_obj == self.root.ids.taxa_sl:
            if len(parent_obj.children) == 3 and \
                    len(self.alignment_list.taxa_names) == 1:
                return self.remove_all()

        # Get button widgets to be removed
        bt_idx = value.id[:-1]
        inf_idx = value.id[:-1] + "?"
        c_idx = value.id[:-1] + "C"

        # Remove button widgets (name button, info button and remove button)
        try:
            bt = [x for x in parent_obj.children if bt_idx == x.id][0]
            parent_obj.remove_widget(bt)

            # Removes reference to this file/taxa in button and mouse over
            #  vars
            if parent_obj == self.root.ids.file_sl:
                self.sp_file_bts = [x for x in self.sp_file_bts
                                    if x[0].text != bt.text]
                self.mouse_over_bts["Files"].remove(bt)
            elif parent_obj == self.root.ids.taxa_sl:
                self.sp_taxa_bts = [x for x in self.sp_taxa_bts
                                    if x[0].text != bt.text]
                self.mouse_over_bts["Taxa"].remove(bt)

        except (IndexError, ValueError):
            pass
        try:
            inf_bt = [x for x in parent_obj.children if inf_idx == x.id][0]
            parent_obj.remove_widget(inf_bt)
        except IndexError:
            pass
        try:
            cbt = [x for x in parent_obj.children if c_idx == x.id][0]
            parent_obj.remove_widget(cbt)
        except IndexError:
            pass

        # Remove widgets
        parent_obj.remove_widget(value)

        # CORE CHANGES
        # Get the parent tab
        if parent_obj == self.root.ids.file_sl:
            # Set the update_taxa attribute to True so that it updates
            # the next time the taxa information is queried
            self.trigger_taxa_update = True
            # Update file list
            file_path = self.filename_map[bt_idx]
            del self.filename_map[bt_idx]
            # If alignment remove from file list
            if self.file_list:
                self.file_list.remove(file_path)
                # Update alignment object list
                self.alignment_list.remove_file([file_path])
            # If proteome remove, from proteome file
            if self.proteome_files:
                self.proteome_files.remove(file_path)
                # In proteome files, if the total number of proteomes is lower
                # than the set min taxa representation filter, update the
                # filter
                if len(self.proteome_files) < self.orto_min_sp:
                    self.orto_min_sp = len(self.proteome_files)
            # Update active file list. If the file has been removed from the
            # active list, this will handle the exception
            try:
                self.active_file_list.remove(file_path)
            except ValueError:
                pass
            try:
                self.active_proteome_files.remove(file_path)
            except ValueError:
                pass

            # Update active taxa list. This must be executed before calling
            # self.get_taxa_information since this method relies on an
            # updated active taxa list
            self.update_taxa()

            # Updates the partition list
            self.update_partitions()

            # Updates labels
            self.update_file_label()
            self.update_sp_label()
            self.update_partition_label()

        if parent_obj == self.root.ids.taxa_sl:
            self.alignment_list.remove_taxa([bt_idx])
            try:
                self.active_taxa_list.remove(bt_idx)
            except ValueError:
                pass
            # Updates label
            self.update_sp_label()

        if not self.file_list and not self.proteome_files:
            self.clear_process_input()

    def remove_bt_from_selection(self, idx, txt_file=None):
        """
        Adds functionality to the dropdown button options for removing file
        or taxa buttons contained in a text file
        :param idx: string, either 'Files' or 'Taxa'.
        :param txt_file: string, path to txt file containing the files/taxa
        names to be selected. If this is not provided, then the selected
        files/taxa will be removed instead
        """

        self.dismiss_popup()

        selection = []
        selection_idx = []

        # Get sidepanel widget where items will be removed
        wgt = self.root.ids.taxa_sl if idx == "Taxa" else self.root.ids.file_sl

        # Grab selection from text file
        if txt_file:

            with open(txt_file) as fh:

                for line in fh:
                    if line.strip() != "":
                        selection.append(line.strip())
                        selection_idx.append(line.strip() + "X")

        # Grab selection from selected items in the sidepanel
        else:
            selection = self.active_taxa_list if idx == "Taxa" else \
                [basename(x) for x in self.active_file_list]
            selection_idx = [x + "X" for x in selection]

        if idx == "Taxa":
            missing_elements = [x for x in selection if x not in
                                self.alignment_list.taxa_names]
            # All taxa were selected, so remove all
            if not [x for x in self.alignment_list.taxa_names
                    if x not in selection]:
                return self.remove_all()

        else:
            # If a file does not exist in the app, add it here
            missing_elements = [x for x in selection
                                if x not in self.filename_map]
            # All files were selected, so remove all
            if not [x for x in self.filename_map if x not in selection]:
                return self.remove_all()

        # Removes buttons from side panel. In the case of Files, some may
        # not be loaded as a button in the side panel (due to the limit
        # of 20 files by default). These files have to be removed from
        # other structures in the code below
        it = iter([(x, x.id) for x in wgt.children])
        for bt, idx in it:
            if idx in selection_idx:
                self.remove_bt(bt, parent_wgt=wgt)
                selection_idx.remove(idx)

        for idx in selection_idx:
            if idx[:-1] in self.filename_map:
                # Set the update_taxa attribute to True so that it updates
                # the next time the taxa information is queried
                self.trigger_taxa_update = True
                file_path = self.filename_map[idx[:-1]]
                del self.filename_map[idx[:-1]]
                self.file_list.remove(file_path)
                try:
                    self.active_file_list.remove(file_path)
                except ValueError:
                    pass
                self.alignment_list.remove_file([file_path])

        self.update_taxa()
        self.update_partitions()
        self.update_file_label()
        self.update_sp_label()
        self.update_partition_label()

        # If some of the files did not exist in the app, issue a warning
        if missing_elements:
            self.dialog_floatcheck("Warning: {} items could not be removed "
                "because they were not loaded into the "
                "app".format(len(missing_elements)), t="error")

    def select_bt_from_file(self, idx, txt_file):
        """
        Adds functionality to the dropdown button options for selecting
        file or taxa buttons contained in a text file
        :param idx: string, either 'Files' or 'Taxa'.
        :param txt_file: string, path to txt file containing the files/taxa
        names to be selected
        """

        self.dismiss_popup()

        selection = []

        with open(txt_file) as fh:

            for line in fh:
                if line.strip() != "":
                    selection.append(line.strip())

        if idx == "Taxa":

            for bt in [x for x in self.root.ids.taxa_sl.children if
                       isinstance(x, ToggleButton)]:
                if bt.text in selection:
                    bt.state = "down"
                else:
                    bt.state = "normal"

            self.active_taxa_list = [x for x in selection if x in
                                     self.alignment_list.taxa_names]

            self.update_sp_label()

        elif idx == "Files":

            for bt in [x for x in self.root.ids.file_sl.children if
                       isinstance(x, ToggleButton)]:
                if bt.text in selection:
                    bt.state = "down"
                else:
                    bt.state = "normal"

                if self.file_list:

                    self.active_file_list = [
                        self.filename_map[x] for x in
                        selection if x in self.filename_map]
                    self.alignment_list.update_active_alignments(
                        [self.filename_map[x] for x in selection if
                         x in self.filename_map])

                elif self.proteome_files:

                    self.active_proteome_files = [
                        self.filename_map[x] for x in
                        selection if x in self.filename_map]

            self.update_file_label()

    def select_bt(self, value):
        """
        Functionality to the Select All/Deselect All buttons of the side
        panel. The method was made in such a way that it could be of general
        use for buttons in the files and taxa tabs

        :param value: Button widget.
        """

        sv_parent = [x for x in value.parent.parent.parent.children if
                     isinstance(x, ScrollView)][0]

        # This will iterate over the first child of the parent scrollview.
        # Since scroll view only supports one child, this should be fine
        for i in sv_parent.children[0].children:
            # Skips the X buttons
            if isinstance(i, ToggleButton):

                if value.text == "Select All":
                    # App related action
                    i.state = "down"

                elif value.text == "Deselect All":
                    # App related action
                    i.state = "normal"

        # Core changes to files
        if (sv_parent == self.root.ids.sv_file and
                value.text == "Select All"):

            if self.file_list:
                self.active_file_list = self.file_list[:]
                self.alignment_list.update_active_alignments(
                    [x for x in self.file_list])
            elif self.proteome_files:
                self.active_proteome_files = self.proteome_files[:]

        # Core changes to taxa
        if sv_parent == self.root.ids.sv_sp and value.text == "Select All":
            self.active_taxa_list = deepcopy(
                self.alignment_list.taxa_names)

        # Core changes to files
        if (sv_parent == self.root.ids.sv_file and
                value.text == "Deselect All"):
            if self.file_list:
                self.active_file_list = []
                self.alignment_list.update_active_alignments([])
            elif self.proteome_files:
                self.active_proteome_files = []
        # Core changes to taxa
        if (sv_parent == self.root.ids.sv_sp and
                value.text == "Deselect All"):
            self.active_taxa_list = []

        # Core changes to partitions
        if (sv_parent == self.root.ids.sv_partition and
                value.text == "Select All"):
            self.active_partitions = list(
                self.alignment_list.partitions.partitions)
        else:
            self.active_partitions = []

        if (sv_parent == self.root.ids.sv_sp or sv_parent ==
                self.root.ids.sv_file):
            # Updates labels
            self.update_sp_label()
            self.update_file_label()
            self.update_partition_label()

    def dialog_dataset_creator(self, ds_type, popup_level=1):
        """
        Creates a dialog to choose between creating a data set from a file
        or manually

        :param ds_type: string. Identifier of the data type. Can be either
        taxa or files
        :param popup_level: integer, The level of the dialog popup. Can be
        either 1 (creates _popup instance) or 2 (creates _subpopup instance)
        """

        content = DataSetTriageDialog(cancel=self.dismiss_popup)
        content.ds_type = ds_type
        content.popup_level = popup_level

        self.show_popup(title="Choose data set group creation method",
                        content=content, size=(360, 200))

    def dialog_taxagroup(self, ds_type, popup_level=1):
        """
        Creates the layout for the taxa group creation popup.
        :param ds_type: string. Data set type. It may be either "taxa" or
        "files"
        :param popup_level: integer, The level of the dialog popup. Can be
        either 1 (creates _popup instance) or 2 (creates _subpopup instance)
        """

        # Initializing instance for taxa group dialog
        if popup_level == 1:
            content = TaxaGroupDialog(cancel=self.dismiss_popup)
        else:
            content = TaxaGroupDialog(cancel=self.dismiss_subpopup)

        if ds_type == "taxa":
            bt_list = sorted(self.alignment_list.taxa_names, reverse=True)
            title = "Create taxa groups"
            group_list = sorted(self.taxa_groups.keys())
        else:
            bt_list = sorted([basename(x) for x in self.file_list],
                             reverse=True)
            title = "Create file groups"
            group_list = sorted(self.file_groups.keys())

        # Reset previous button for multi selection support
        self.prev_tb = None

        # Populate the gridlayout for all entries
        for i in bt_list:
            # Create togglebutton for each entry
            bt = ToggleButton(text=i, size_hint_y=None, height=30)
            bt.bind(on_release=self.toggle_multi_selection)
            self.add_dataset_bt(bt, content.ids.all_grid, ds_type)

        # Populate created groups, if any
        if group_list:
            content.ids.group_list.clear_widgets()

        for i in group_list:

            self.taxagroups_add_group(i, content.ids.group_list, ds_type)

        content.ds_type = ds_type

        # Show dialog

        self.show_popup(title=title, content=content, size=(900, 600),
                        popup_level=popup_level)

    def add_dataset_bt(self, bt, wgt, ds_type):
        """
        Method for addition of a button to a widget. This method was created
        for the automatic upated of the widgets height when moving buttons
        in the taxa group creation dialog
        :param bt: The button widget
        :param wgt: The sink widget
        :param ds_type: string. Data set type. It may be either "taxa" or
        "files"
        """

        if ds_type == "taxa":
            bt_list = sorted(self.alignment_list.taxa_names, reverse=True)
        else:
            bt_list = sorted([basename(x) for x in self.file_list],
                             reverse=True)

        wgt.add_widget(bt, bt_list.index(bt.text))
        wgt.height += 30

    @staticmethod
    def remove_taxa_bt(bt, wgt):
        """
        Method for addition of a button to a widget. This method was created
        for the automatic upated of the widgets height when moving buttons
        in the taxa group creation dialog
        :param bt: The button widget
        :param wgt: The source widget
        """
        wgt.remove_widget(bt)
        wgt.height -= 30

    def taxagroup_move_taxa(self, source_wgt, sink_wgt, all_taxa, ds_type):
        """
        Method that adds functionality to the addition/removal buttons
        (<<, <, >>, >) in the taxa group dialog.
        :param source_wgt: widget, the gridlayout from where the buttons
        will be moved
        :param sink_wgt: widget, the gridlayout to where buttons will be
        moved
        :param all_taxa: Boolean, if True its as if alsa taxa were selected
        to be moved
        :param ds_type: string. Data set type. It may be either "taxa" or
        "files"
        """

        if ds_type == "taxa":
            bt_list = sorted(self.alignment_list.taxa_names, reverse=True)
        else:
            bt_list = sorted([basename(x) for x in self.file_list],
                             reverse=True)

        # In case all taxa are to be moved
        if all_taxa:
            # Ensures that only toggle buttons are moved
            for tx in bt_list:
                try:
                    bt = [x for x in source_wgt.children if tx == x.text][0]
                    self.remove_taxa_bt(bt, source_wgt)
                    bt.state = "normal"
                    self.add_dataset_bt(bt, sink_wgt, ds_type)
                except IndexError:
                    pass
        else:
            # This workaround is used to add some buttons from the source
            #  to the sink widgets while maintaining their original
            # order. The z-index of widgets is not working quite as I
            # expected, so for now this will produce the desired behaviour
            sink_bts = []
            # Remove buttons from the sink widget and store their names into
            # a storage list
            for bt in sink_wgt.children[::-1]:
                self.remove_taxa_bt(bt, sink_wgt)
                sink_bts.append(bt.text)

            # Gather the buttons from the source widget to be transferred
            # while removing them from the sink
            for bt in source_wgt.children[::-1]:
                if bt.state == "down":
                    self.remove_taxa_bt(bt, source_wgt)
                    sink_bts.append(bt.text)

            # Add buttons to the sink widget in the desired order
            for tx in sorted(sink_bts, reverse=True):
                bt = ToggleButton(text=tx, size_hint_y=None, height=30)
                bt.bind(on_release=self.toggle_multi_selection)
                self.add_dataset_bt(bt, sink_wgt, ds_type)

    def taxagroups_show_taxa(self, name_wgt):
        """
        Creates a popup listing the taxa included in a taxa group given by
        name
        :param name_wgt: widget, widget containing the name of the group as
        text
        """

        # Create root boxlayout
        content = BoxLayout(orientation="vertical", padding=10, spacing=10)
        # Create scroll view in which the gridlayout will be inserted
        sv = ScrollView(bar_width=10, scroll_type=["bars"])
        # Create close button for the popup
        close_bt = Button(text="Close", size_hint_y=None, height=30,
            bold=True,
            background_normal="data/backgrounds/bt_process.png",
            bakcground_down="data/backgrounds/bt_process_off.png")
        # Add functionality to the close button
        close_bt.bind(on_release=self.dismiss_popup)
        # Create gridlayout that will store the buttons with taxa names
        gl = GridLayout(cols=1, size_hint_y=None, height=0, spacing=5)

        # Get ds_type
        parent_wgt = name_wgt.parent

        # If its taxa data set
        if parent_wgt == self.root.ids.taxa_group_grid:
            bt_list = self.taxa_groups[name_wgt.text]
        else:
            bt_list = self.file_groups[name_wgt.text]

        # Create buttons
        for tx in bt_list:
            bt = Button(text=tx, size_hint_y=None, height=30)
            gl.add_widget(bt)
            gl.height += 35

        # Create widget tree
        sv.add_widget(gl)
        content.add_widget(sv)
        content.add_widget(close_bt)

        # Show dialog
        self.show_popup(title="Taxa group: %s" % name_wgt.text,
                        content=content, size_hint=(.3, .7))

    def remove_taxa_group(self, rm_wgt):
        """
        Removes the data set group button from the app list and
        corresponding data set group attribute
        :param rm_wgt: widget, widget of the removal button
        """

        # Remove from app
        parent_wgt = rm_wgt.parent

        bt_idx = rm_wgt.id[:-1]

        try:
            bt = [x for x in parent_wgt.children if bt_idx == x.id][0]
            parent_wgt.remove_widget(bt)
            parent_wgt.remove_widget(rm_wgt)
        except IndexError:
            bt = [x for x in parent_wgt.children if bt_idx == x.text][0]
            # Remove from the dataset dialog
            parent_wgt.remove_widget(bt)
            parent_wgt.remove_widget(rm_wgt)
            # Remove from the sidepanel
            if parent_wgt.ds == "taxa":
                for i in [x for x in self.root.ids.taxa_group_grid.children
                          if x.id == rm_wgt.id or x.id == bt_idx]:
                    self.root.ids.taxa_group_grid.remove_widget(i)
            else:
                for i in [x for x in self.root.ids.file_group_grid.children
                          if x.id == rm_wgt.id or x.id == bt_idx]:
                    self.root.ids.file_group_grid.remove_widget(i)

        # Remove from program attribute
        if parent_wgt == self.root.ids.taxa_group_grid or \
                parent_wgt.ds == "taxa":
            # Remove group from core attribute
            del self.taxa_groups[bt_idx]
            # Remove button from dropdown menu
            # Since the children of a dropdown widget are a gridlayout and
            # not the actual buttons contained in the dropdown menu,
            # this will search for the children of the gridlayout
            for i in [x for x in
                      self.process_grid_wgt.ids.taxa_dropdown.
                      children[0].children if x.text == bt_idx]:

                self.process_grid_wgt.ids.taxa_dropdown.remove_widget(i)
            # Remove button from sidepanel
            for i in [x for x in self.root.ids.taxa_group_grid.children
                      if x.text == bt_idx]:
                self.root.ids.taxa_group_grid.remove_widget(i)
        if parent_wgt == self.root.ids.file_group_grid or \
                parent_wgt.ds == "files":
            # Remove group from core attribute
            del self.file_groups[bt_idx]
            # Remove button from dropdown menu
            for i in [x for x in
                      self.process_grid_wgt.ids.file_dropdown.
                      children[0].children if x.text == bt_idx]:

                self.process_grid_wgt.ids.file_dropdown.remove_widget(i)
            # Remove button from sidepanel
            for i in [x for x in self.root.ids.file_group_grid.children
                      if x.text == bt_idx]:
                self.root.ids.file_group_grid.remove_widget(i)

    def taxagroups_add_group(self, name, wgt, ds_type):
        """
        Adds a dataset button, and corresponding removal button, to the
        group list gridlayut of the Dataset dialog.
        :param name: string, name of the group
        :param wgt: GridLayout widget where the buttons will be added
        :param ds_type: string, dataset type, whether 'taxa' or 'files'
        """

        bt = ShortenToggleButton(text=name, size_hint_y=None, height=30,
            state="normal", disabled_color=(1, 1, 1, 1))

        bt.bind(on_release=self.toggle_groups)
        bt.bind(on_release=lambda x: self.taxagroups_display_group(name,
                                                                   ds_type))
        rm_bt = Button(size_hint=(None, None), width=30,
            height=30, id="{}X".format(name), border=(0, 0, 0, 0),
            background_normal=join("data", "backgrounds", "remove_bt.png"),
            background_down=join("data", "backgrounds",
                                 "remove_bt_down.png"))

        rm_bt.bind(on_release=partial(
            self.check_action,
            "Are you sure you want to remove this group?",
            self.remove_taxa_group,
            popup_level=2))

        wgt.add_widget(bt)
        wgt.add_widget(rm_bt)

    def taxagroups_display_group(self, name, ds_type):
        """
        Gives functionality to the group buttons in the dataset dialog. It
        displays the group constitution in the dialog
        :param name:
        :param ds_type:
        """

        # Get dataset dialog
        root_wgts = self.root_window.children
        if [x for x in root_wgts if isinstance(x, Popup)]:
            dataset_wgt = [x for x in root_wgts if
                           isinstance(x, Popup)][0].content
        elif [x for x in root_wgts if isinstance(x, CustomPopup)]:
            dataset_wgt = [x for x in root_wgts if
                           isinstance(x, CustomPopup)][0].content
        else:
            dataset_wgt = None

        # Reset dataset creator source and sink gridlayouts
        self.taxagroup_move_taxa(dataset_wgt.ids.select_grid,
                                 dataset_wgt.ids.all_grid,
                                 True,
                                 dataset_wgt.ds_type)

        # Get group names list
        if ds_type == "taxa":
            group_lst = self.taxa_groups[name]
        else:
            group_lst = self.file_groups[name]

        # Display group in the selected gridlayout
        for i in group_lst:
            bt = [x for x in dataset_wgt.ids.all_grid.children
                  if i == x.text][0]
            self.remove_taxa_bt(bt, dataset_wgt.ids.all_grid)
            bt.state = "normal"
            self.add_dataset_bt(bt, dataset_wgt.ids.select_grid, ds_type)

        # Change dataset group name in the dialog
        dataset_wgt.ids.group_name.text = name

    def save_dataset_group(self, source_wgt, name, ds_type,
                           group_file=False):
        """
        Adds a taxa group declared using the taxa group creator popup to the
        list of taxa groups in the side panel
        :param source_wgt, gridlayout of the selected items
        :param name: string, name of the group
        :param ds_type: string. Data set type. It may be either "taxa" or
        "files"
        :param group_file: boolean, If True get the group items from the
        self.dataset_file file.
        """

        # Determine if the Dataset dialog is still active
        root_wgts = self.root_window.children
        if [x for x in root_wgts if isinstance(x, Popup)]:
            dataset_wgt = [x for x in root_wgts if
                           isinstance(x, Popup)][-1].content
        elif [x for x in root_wgts if isinstance(x, CustomPopup)]:
            dataset_wgt = [x for x in root_wgts if
                           isinstance(x, CustomPopup)][-1].content
        else:
            dataset_wgt = None

        # This resets the dataset_wgt if its not the correct instance. This
        # happens, for example, when setting data set from file
        if not isinstance(dataset_wgt, TaxaGroupDialog):
            dataset_wgt = None

        if ds_type == "taxa":
            # Make core changes by populating self.taxa_groups dictionary
            self.taxa_groups[name] = []
            group_list = self.taxa_groups[name]
            # This list will be used to ensure that when groups are provided
            # via text file, only compliant entries are saved
            check_list = self.alignment_list.taxa_names
            # Set the grid layout where the group button is to be added
            grid_layout = self.root.ids.taxa_group_grid
            # Set the current dataset group as default for taxa filter if
            # it has not been defined
            self.taxa_filter_settings[1] = name

        else:
            # Make core changes by populating self.file_groups dictionary
            self.file_groups[name] = []
            group_list = self.file_groups[name]
            # This list will be used to ensure that when groups are provided
            # via text file, only compliant entries are saved
            check_list = list(self.filename_map.keys())
            # Set the grid layout where the group button is to be added
            grid_layout = self.root.ids.file_group_grid

        if group_file:
            with open(self.dataset_file) as fh:
                for line in fh:
                    if line.strip() in check_list:
                        group_list.append(line.strip())
        else:
            for bt in source_wgt.children:
                group_list.append(bt.text)

        # If dataset dialog is still active, add the new group
        if dataset_wgt:

            # Remove original button when not groups have been previously
            # added
            if len(dataset_wgt.ids.group_list.children) == 1:
                dataset_wgt.ids.group_list.clear_widgets()

            # Add group button to layout
            if name not in [x.text for x in
                            dataset_wgt.ids.group_list.children]:
                self.taxagroups_add_group(name, dataset_wgt.ids.group_list,
                                          ds_type)

            # Reset dataset creator source and sink gridlayouts
            self.taxagroup_move_taxa(dataset_wgt.ids.select_grid,
                                     dataset_wgt.ids.all_grid,
                                     True,
                                     dataset_wgt.ds_type)

            for bt in [x for x in dataset_wgt.ids.group_list.children
                       if isinstance(x, TGToggleButton)]:
                bt.state = "normal"
                bt.disabled = False

        # App changes by adding two buttons for the taxa group
        # Taxa button itself
        if name not in [x .text for x in grid_layout.children]:
            bt = ShortenButton(text=name, size_hint_y=None,
                                     height=30, id=name)
            bt.bind(on_release=self.taxagroups_show_taxa)
            # Removal button
            x_bt = Button(size_hint=(None, None), width=30,
                          border=(0, 0, 0, 0), height=30, id="%sX" % name,
                          background_normal=join("data", "backgrounds",
                                                 "remove_bt.png"))
            x_bt.bind(on_release=partial(self.check_action,
                                         "Are you sure you want to remove this"
                                         " group?",
                                         self.remove_taxa_group))

            # Add buttons to gridlayout
            for i in [bt, x_bt]:
                grid_layout.add_widget(i)

    def projects_init(self):
        """
        Initializes the projects attribute
        """

        if os.path.exists(self.projects_file):

            with open(self.projects_file, "rb") as projects_fh:
                # Attempts to read the projects file. If the file is
                # somewhat corrupted, remove it and issue a warning
                try:
                    projects_dic = pickle.load(projects_fh)
                except EOFError:
                    os.remove(self.projects_file)
                    self.projects_init()
                # Get ordered dict
                projects_dic = OrderedDict(sorted(projects_dic.items(),
                                                  key=lambda x: (x[1][1],
                                                                 x[0])))

                # Populate sidepanel
                for name, p in projects_dic.items():
                    project_grid = self.root.ids.project_grid
                    self.add_project_bt(name, len(p[0]), p[1], project_grid)

                    # Populates main screen if active
                    if "home_projects_grid" in self.screen.ids:
                        project_grid = self.screen.ids.home_projects_grid
                        self.add_project_bt(name, len(p[0]), p[1],
                                            project_grid)

        # This handles cases where the projects file has
        # not been populated yet
        else:
            projects_dic = OrderedDict()
            with open(self.projects_file, "wb") as projects_fh:
                pickle.dump(projects_dic, projects_fh)

    def save_project(self, name):
        """
        Saves the current alignment_list or proteome_list for quick access
        to the data sets. It automatically detects whether the current
        data set is from the orthology or process modules.

        :param name: string, name of the project
        """

        # Get projects var
        with open(self.projects_file, "rb") as projects_fh:
            projects_dic = pickle.load(projects_fh)

        # Check for project name duplicates
        if name in projects_dic:
            return self.dialog_floatcheck(
                "Project with the same name is already present",
                t="warning")

        with open(self.projects_file, "wb") as projects_fh:
            if self.file_list:
                projects_dic[name] = [self.file_list, "process"]
                pickle.dump(projects_dic, projects_fh)
                ds_type = "process"

            else:
                projects_dic[name] = [self.proteome_files, "orthology"]
                pickle.dump(projects_dic, projects_fh)
                ds_type = "orthology"

        project_grid = self.root.ids.project_grid
        self.add_project_bt(name, len(self.file_list), ds_type,
                            project_grid)

        if "home_projects_grid" in self.screen.ids:
            project_grid = self.screen.ids.home_projects_grid
            self.add_project_bt(name, len(self.file_list), ds_type,
                                project_grid)

    def add_project_bt(self, name, file_num, ds_type, grid_wgt):
        """
        Wrapper that adds a project button to the sidepanel
        :param name: string, name of the project
        :param file_num: int, number of files associated with project
        :param ds_type: string, type of dataset. Can be either 'orthology'
        or 'process'
        :param grid_wgt: GridLayut widget where the project is to be added
        """

        if len(grid_wgt.children) == 1:
            grid_wgt.clear_widgets()

        if ds_type == "orthology":
            ds_bt = ProjectOrtoBt(id=name)
        else:
            ds_bt = ProjectProcBt(id=name)

        bt = TFButton(text=name, size_hint_y=None, height=35, bold=True,
                      id=name)

        # Set markup message for project name
        msg1 = "[b][color=ccccccff]Project name:[/color][/b] {}".format(
            name.encode("utf8"))

        # Set markup message for number of files
        msg2 = "[b][color=ccccccff]Number of files:[/color][/b] {}".format(
            str(file_num))

        if ds_bt == "orthology":
            pass
        else:
            bt.bind(on_release=partial(
                self.check_action,
                [msg1, msg2],
                self.open_project,
                **{"args": [name],
                   "check_wgt": CheckProject}))

        rm_bt = Button(size_hint=(None, None), size=(35, 35),
            id=name, border=(0, 0, 0, 0),
            background_normal=join("data", "backgrounds",
                                   "remove_bt35.png"),
            background_down=join("data", "backgrounds",
                                 "remove_bt35_down.png"))
        rm_bt.bind(on_release=partial(
            self.check_action,
            "Are you sure you want to remove this project?",
            self.remove_project))

        for i in [ds_bt, bt, rm_bt]:
            grid_wgt.add_widget(i)

    def remove_project(self, wgt):
        """
        Gives functionality to the remove project button
        :param wgt: remove button
        """

        # Remove from projects file
        with open(self.projects_file, "rb") as project_fh:
            project_dic = pickle.load(project_fh)

        del project_dic[wgt.id]

        with open(self.projects_file, "wb") as project_fh:
            pickle.dump(project_dic, project_fh)

        # Remove from sidepanel
        for i in [bt for bt in self.root.ids.project_grid.children
                  if bt.id == wgt.id]:
            self.root.ids.project_grid.remove_widget(i)

        if len(self.root.ids.project_grid.children) == 0:
            self.root.ids.project_grid.add_widget(Button(
                text="No Saved Projects", disabled=True, size_hint_y=None,
                height=35))

        # Remove from home screen, if present
        if "home_projects_grid" in self.screen.ids:
            for i in [bt for bt in
                      self.screen.ids.home_projects_grid.children
                      if bt.id == wgt.id]:
                self.screen.ids.home_projects_grid.remove_widget(i)

            if len(self.screen.ids.home_projects_grid.children) == 0:
                self.screen.ids.home_projects_grid.add_widget(Button(
                    text="No Saved Projects", disabled=True,
                    size_hint_y=None, height=35))

    def open_project(self, name):
        """
        Closes the current data set and opens the project identified by name
        :param name: string, name of the project
        """

        with open(self.projects_file, "rb") as projects_fh:
            project_dic = pickle.load(projects_fh)

        if name in project_dic:
            if project_dic[name][1] == "process":
                # Closes current data sets
                self.remove_all()
                # Check if files are present
                if [x for x in project_dic[name][0] if not os.path.isfile(
                        x)]:
                    self.dialog_floatcheck("Some project "
                        "files no longer exist", t="warning")
                # Opens new dataset
                self.load_files_subproc(project_dic[name][0])

            if project_dic[name][1] == "orthology":
                # Closes current data sets
                self.remove_all()
                # Opens new data set
                self.load_proteomes(project_dic[name][0])

    def more_stats_table_bts(self):
        """
        Increases the MAX_TABLE_N attribute by 50 and adds more buttons
        to the summary statistics gene table view
        """

        # Get results
        with open(join(self.temp_dir, "table.pc"), "rb") as fh:
            table = pickle.load(fh)

        # Remove previous bt
        self.stats_table.ids.table_grid.remove_widget(
            self.stats_table.ids.table_grid.children[0])

        # If the current table has a search operation performed, show more
        # lines from that searched table. If not, show more from the original
        # table
        if any(self.gene_table_selection):
            table = self.gene_table_selection
        else:
            table = self.gene_table

        start = self.MAX_TABLE_N
        self.MAX_TABLE_N += 25

        self.search_add_gene_table_line(table, start)

    def statistics_show_summary(self, force=False, tp="stats"):
        """
        Creates a subprocess for calculating summary statistics for the
        curently acitve data set. It will show either stats or table
        mode, depending on the option provided in the 't' argument.
        :param force: Bool, Clear all widgets from Stats scatter and show
        summary
        :param tp: string, Show overall summary statistics ('stats') or
        summary statistics for each gene in table view ('table')
        """

        def prepare_display(plt_wgt, scatter):

            if self.screen.name == "Statistics":
                # Clear main widget
                plt_wgt.clear_widgets()

                # Remove any potential widgets associated with stats plots
                self.screen.ids.footer_box.clear_widgets()
                self.dismiss_stats_toggle()
                self.dismiss_plot_wgt()

                # Remove any potential widgets associated with stats plots
                self.screen.ids.footer_box.clear_widgets()
                self.dismiss_stats_toggle()
                self.dismiss_plot_wgt()

                # Disable scale, translation and rotation of scatter
                scatter.do_scale = False
                scatter.do_translation = False
                scatter.do_rotation = False
                scatter.scale = 1
                scatter.pos = 0, 0

        def display_stats(plt_wgt):

            # Get results
            with open(join(self.temp_dir, "stats.pc"), "rb") as fh:
                stats, table = pickle.load(fh)

            # Check for exceptions
            if "exception" in stats:
                return self.dialog_floatcheck(stats["exception"], t="error")

            # Populate cards content
            for k, v in stats.items():
                self.stats_summary.ids[k].text = "%s" % v

            # Add stats widget if in Screen
            if self.screen.name == "Statistics":
                plt_wgt.add_widget(self.stats_summary)

            # Set table attribute
            self.current_table = table

        def display_table(plt_wgt):

            # Force statistics side panel close
            self.toggle_stats_panel(force_close=True)

            # Clear Table grid
            self.stats_table = GeneTable()

            # Get results
            with open(join(self.temp_dir, "table.pc"), "rb") as fh:
                table, td = pickle.load(fh)

            # Check for exceptions
            if "exception" in table:
                return self.dialog_floatcheck(table["exception"], t="error")

            # Issue warning for large data sets
            if len(table) > 50:
                self.dialog_floatcheck("Large number of genes to "
                    "display. Showing first 50 for now. For larger data "
                    "sets it is recommended to use the 'export table' "
                    "option", t="warning", dl=10)

            self.search_add_gene_table_line(table)

            # Add Table widget if in screen
            if self.screen.name == "Statistics":
                plt_wgt.add_widget(self.stats_table)

            # Set table attribute
            self.current_table = td
            self.gene_table = table

        def check_process(p, ldg_wgt, plt_wgt, ns, dt):

            try:
                ldg.ids.msg.text = "{}/{}".format(ns.counter, ns.total)
            except AttributeError:
                pass

            # KILL SWITCH. Terminating thread by user
            if self.terminate_stats:

                # Unschedule the current function
                Clock.unschedule(check_func)

                # Issuing the kill order to the child thread. This signal
                # will propagate through the worker's methods and graciously
                # terminate
                ns.stop = True

                # Clear the loading widgets of the Statistics screen
                plt_wgt.clear_widgets()
                plt_wgt.add_widget(NoDataLabel())

                # Remove lock from get_summary_stats method, so that this
                # method can be repeated
                self.lock_stats = False

                # Shutting down manager that provides communication between
                # parent and child thread
                manager.shutdown()

                # Join child process and exit
                p.join()
                return

            # Provides animation for loading widget
            if self.current_screen == "Statistics":

                if self.screen.ids.stats_panel.width == 0 or \
                        self.screen.ids.stats_panel.width == 410:

                    # Sets animation speed
                    ldg_wgt.ids.img.rotation -= 10

            # When worker has finished
            if not p.is_alive():

                # Set the previous sets as the ones that just finished without
                # errors
                self.previous_sets["Stats"] = [file_set, taxa_set]

                # Display stats depending on the type
                prepare_display(plt_wgt, scatter_wgt)
                if tp == "stats":
                    display_stats(plt_wgt)
                elif tp == "table":
                    display_table(plt_wgt)

                # Unschedule the current function
                Clock.unschedule(check_func)

                # Shutting down manager that provides communication between
                # parent and child thread
                manager.shutdown()

                # Joining child thread
                p.join()

                # Remove lock from get_summary_stats()
                self.lock_stats = False

                # Exit
                return

        # Get Scatter widgets
        scatter_wgt = self.screen.ids.plot_content
        plot_wgt = self.screen.ids.plot_content.children[0]

        # If no files have been loaded, ignore
        if not self.file_list:
            return

        # If a process is already running, ignore
        if self.lock_stats:
            return

        # Check if plot has been previously loaded. If yes, ignore.
        # If there are no children, proceed as usual
        try:
            if isinstance(plot_wgt.children[0], Image) and not force:
                return
        except IndexError:
            pass

        # Check if any input alignments are active. If not, depends on
        # whether there is a summary statistics widget already in
        # display. If there is, ignore. If not, use the latest sum stats
        # pickle object to quickly display the stats
        if not self.active_file_list and self.file_list:
            return self.dialog_floatcheck("No active files have been "
                                          "selected", t="error")

        # Get currently active data sets
        file_set_name = self.screen.ids.active_file_set.text
        taxa_set_name = self.screen.ids.active_taxa_set.text

        file_set, taxa_set = self.get_active_sets(file_set_name,
                                                  taxa_set_name)

        # Check if currently active data sets have changed since the last
        #  calculation of summary statistics
        if self.previous_sets["Stats"] == [file_set, taxa_set]:
            # Active data sets are the same. If summary stats are already
            # displayed, ignore
            if isinstance(plot_wgt.children[0], StatsSummary) and \
                    tp == "stats":
                return
            elif isinstance(plot_wgt.children[0], GeneTable) and \
                    tp == "table":
                return
            # If summary stats are not displayed, then retrieve
            # information from the pickle objects and draw the widgets
            else:
                if tp == "stats":
                    prepare_display(plot_wgt, scatter_wgt)
                    display_stats(plot_wgt)
                elif tp == "table":
                    prepare_display(plot_wgt, scatter_wgt)
                    display_table(plot_wgt)
                return

        # Prepare screen by removing any potential plot widgets and
        # changing some properties of the scatter layout
        prepare_display(plot_wgt, scatter_wgt)

        # Disable stats subprocess interruption. Being True by default
        # and explicitly setting is to False when calculating summary
        # stats makes it much more easier to interrupt by other
        # actions/widgets
        self.terminate_stats = False
        self.lock_stats = True

        manager = multiprocessing.Manager()
        shared_ns = manager.Namespace()
        shared_ns.stop = False
        shared_ns.counter = 0

        p = threading.Thread(target=get_stats_summary,
                             kwargs={"aln_list": self.alignment_list,
                                     "dest": self.temp_dir,
                                     "active_file_set": file_set,
                                     "active_taxa_set": taxa_set,
                                     "ns": shared_ns})

        p.daemon = True
        p.start()

        # Add loading widget
        ldg = StatsLoading()
        scatter_wgt.add_widget(ldg)

        check_func = partial(check_process, p, ldg,
                             self.screen.ids.plot_content.children[0],
                             shared_ns)
        Clock.schedule_interval(check_func, .1)

    def search_add_gene_table_line(self, gene_table, start=0):
        """
        Wraps the creation of gene table lines into the table_grid widget.
        :param gene_table
        :return:
        """

        cols = ["nsites", "taxa", "var", "inf", "gap", "missing"]

        for p, row in enumerate(gene_table.itertuples()):

            # Ignore entries before the provided starting point
            if p <= start:
                pass

            if p > self.MAX_TABLE_N:
                bt = MoreTableBt()
                self.stats_table.ids.table_grid.add_widget(bt)
                break

            # Create and populate TableLine
            x = TableLine()
            x.ids.gn_name.text = "{}. {}".format(row.Index + 1, row.genes)
            x.ids.gn_name.size_hint_x = 0.5

            for idx, val in zip(cols, row[2:]):
                x.ids[idx].text = "{}".format(int(val))

            # Add TableLine to main grid
            self.stats_table.ids.table_grid.add_widget(x)

    def search_statistics_gene_table(self, s):
        """
        Searches the summary statistic gene table view, according to the
        provided string.
        :param s: string. Used to search alignment names.
        """

        # Store the master gene table on another attribute, so that it can be
        # later retrieved
        self.gene_master_table = self.current_table[:]

        # Info for new table data
        table = [["Gene name", "Number of sites", "Number of taxa",
                  "Variable sites", "Informative sites", "Gaps",
                  "Missing data"]]

        # clear gene table widget
        self.stats_table.ids.table_grid.clear_widgets()

        # Search the gene table object for gene names matching the provided
        # string
        self.gene_table_selection = self.gene_table[
            self.gene_table["genes"].str.contains(s)]

        self.search_add_gene_table_line(self.gene_table_selection)

        # Create new table data
        for k in self.gene_table_selection.itertuples():
            # Add table line
            table.append(list(k[1:]))

        self.current_table = table

    def search_statistics_clear(self):
        """
        Clears the search in the gene table summary statistics
        """

        # Only do this if there is an active selection being displayed
        if any(self.gene_table_selection):
            # Clear gene_table widget
            self.stats_table.ids.table_grid.clear_widgets()

            # Re-populate with original genes
            self.search_add_gene_table_line(self.gene_table)

            # Reset gene_table_selection
            self.gene_table_selection = []

            # Reset search field hint text
            self.stats_table.ids.gn_txt.text = ""

            # Set the current_table attribute to the total gene list
            self.current_table = self.gene_master_table[:]

    def search_sort_gene_table(self, sortby, ascending=True):
        """
        Sorts the gene table by any column header.
        :param str sortby: Column header name to perform the sort
        :param bool ascending: Whether the sorting will be ascending or
        descending
        """

        col_map = {"Gene name": "genes",
                   "N": "nsites",
                   "Taxa": "taxa",
                   "V": "var",
                   "P": "inf",
                   "Gap": "gap",
                   "M": "missing"}

        # Info for new table data
        table = [["Gene name", "Number of sites", "Number of taxa",
                  "Variable sites", "Informative sites", "Gaps",
                  "Missing data"]]

        # clear gene table widget
        self.stats_table.ids.table_grid.clear_widgets()

        self.gene_table.sort_values([col_map[sortby]],
                                   ascending=ascending,
                                   inplace=True)

        self.search_add_gene_table_line(self.gene_table)

        # Create new table data
        for k in self.gene_table.itertuples():
            # Add table line
            table.append(list(k[1:]))

        self.current_table = table

    def statistics_populate_groups(self, ds_type):
        """
        This method is called when the dataset selection buttons in the
        Statistics sidepanel are pressed. They populate the respective
        dropdown menu with the currently set groups
        :param ds_type: string, data set type. Can be either 'taxa' or
        'files'
        """

        def set_dd_txt(txt, wgt):
            """
            Using this function instead of lamdba to set the callback on the
            on_release event of the dropdown buttons ensures that each button
            has each unique callback and it's not overwritte (as in the lamdba)
            :param txt: string, text provided to select() method
            :param wgt: dummy argument, required because of how the bind()
            method works
            """
            dd_wgt.select(txt)

        if ds_type == "taxa":
            dd_wgt = self.screen.ids.taxa_dropdown
            grid_children = \
                self.screen.ids.taxa_dropdown.children[0].children
            group_atr = self.taxa_groups
        else:
            dd_wgt = self.screen.ids.file_dropdown
            grid_children = \
                self.screen.ids.file_dropdown.children[0].children
            group_atr = self.file_groups

        # Remove discarded groups
        for bt in grid_children:
            if bt.text in ["All taxa", "Active taxa", "All files",
                           "Active files"]:
                pass
            elif bt.text not in group_atr:
                dd_wgt.remove_widget(bt)

        # Add new groups
        current_groups = [x.text for x in grid_children]
        for g in group_atr:
            if g not in current_groups:
                dd_bt = Button(
                    text=g, size_hint_y=None, height=40, bold=True,
                    background_normal=join("data", "backgrounds",
                                           "spinner_opt.png"))
                dd_bt.bind(on_release=partial(set_dd_txt, g))
                dd_wgt.add_widget(dd_bt)

    def dialog_general_info(self, idx):
        """
        Generates the popup with information for several components of the
        application
        :param idx: string. Identifier of the informative content to be
        shown. It must be present in the dictionary keys of the
        informative_storage variable in data/resources/info_data.py
        """

        content = InfoPopup(cancel=self.dismiss_popup)

        # Retrieve title and body text
        title_str, body_str = informative_storage[idx]

        # Add text body
        content.ids.content.text = body_str

        self.show_popup(title=title_str, content=content,
                        size=(500, 400), close_bt=True)

    def operation_queue_init(self):
        """
        This will create the skeleton for the tree view of the operations
        queued for execution in an appropriate screen of the side panel.
        """

        # Create tree view instance
        self.operation_tv = TreeView(hide_root=True, size_hint_y=None,
                                     height=269)

        # Create main nodes for each module
        ortho_node = self.operation_tv.add_node(
            TreeViewLabel(text="Orthology Operations", bold=True,
                          font_size=20, color=(1, 0.3, 0.3, .2)))
        proc_node = self.operation_tv.add_node(
            TreeViewLabel(text="Process Operations", bold=True,
                          font_size=20, color=(.3, .3, 1, 1)))

        # Main subnodes for Process
        main_op_node = self.operation_tv.add_node(TreeViewLabel(
            text="Main Operation", bold=True, font_size=15, opacity=.2),
            proc_node)
        secondary_op_node = self.operation_tv.add_node(TreeViewLabel(
            text="Secondary Operations", bold=True, font_size=15,
            opacity=.2),
            proc_node)
        format_node = self.operation_tv.add_node(TreeViewLabel(
            text="Output Formats", bold=True, font_size=15, opacity=.2),
            proc_node)
        main_file = self.operation_tv.add_node(TreeViewLabel(
            text="Output File", bold=True, font_size=15, opacity=.2),
            proc_node)

        # Save main nodes
        self.main_nodes = {"ortho": ortho_node, "proc": proc_node,
                           "proc_main": main_op_node,
                           "proc_form": format_node,
                           "proc_sec": secondary_op_node,
                           "main_file": main_file}

        self.root.ids.op_sv.add_widget(self.operation_tv)

    def save_operation_queue(self):
        """
        This populate the operations queue tree view with the operations and
        options selected by the user

        As of now, it listens to information stored in:

            - self.main_operations, to gather the main operation
            - self.process_switches, to gather secondary operations
            - self.output_formats, to gather information on the output
            formats
            - self.
        """

        def clear_nodes(parent):
            old_nodes = []
            for node in self.operation_tv.iterate_all_nodes(parent):
                if node.text != parent.text:
                    old_nodes.append(node)
            for node in old_nodes:
                self.operation_tv.remove_node(node)
                self.operation_tv.height -= 24

        def add_node(text, parent):
            self.operation_tv.add_node(
                TreeViewLabel(text=text, font_size=16), parent)
            self.operation_tv.height += 24
            if parent.opacity != 1:
                parent.opacity = 1

        # PROCESS NODES
        # Main operation
        # Clear old nodes
        clear_nodes(self.main_nodes["proc_main"])
        try:
            main_op = [nm for nm, bl in self.main_operations.items()
                       if bl is True][0]
            add_node("%s" % main_op, self.main_nodes["proc_main"])
            # Open Process node
            if self.main_nodes["proc"].is_open is False:
                self.operation_tv.toggle_node(self.main_nodes["proc"])
            # Open main operation node if closed
            if self.main_nodes["proc_main"].is_open is False:
                self.operation_tv.toggle_node(self.main_nodes["proc_main"])
        except IndexError:
            self.main_nodes["proc_main"].opacity = .2

        # Output format
        # Clear old nodes
        clear_nodes(self.main_nodes["proc_form"])
        if self.output_formats:
            for ft in self.output_formats:
                add_node("%s" % ft, self.main_nodes["proc_form"])
            # Open output format node is closed
            if self.main_nodes["proc_form"].is_open is False:
                self.operation_tv.toggle_node(self.main_nodes["proc_form"])
        else:
            self.main_nodes["proc_form"].opacity = .2

        # Secondary operations
        # Clear old nodes
        clear_nodes(self.main_nodes["proc_sec"])
        secondary_op = [nm for nm, bl in self.secondary_operations.items()
                        if bl is True]
        if secondary_op:
            for op in secondary_op:
                add_node("%s" % op, self.main_nodes["proc_sec"])

            if self.main_nodes["proc_sec"].is_open is False:
                self.operation_tv.toggle_node(self.main_nodes["proc_sec"])
        else:
            self.main_nodes["proc_sec"].opacity = .2

        # Output file
        clear_nodes(self.main_nodes["main_file"])
        # for conversion
        if self.main_operations["conversion"]:
            add_node("[Based on input] (main)",
                     self.main_nodes["main_file"])
            if self.main_nodes["main_file"].is_open is False:
                self.operation_tv.toggle_node(self.main_nodes["main_file"])
            # Output files from secondary operations
            if secondary_op:
                for op in secondary_op:
                    if self.secondary_options["%s_file" % op]:
                        add_node("*_%s (%s)" % (op, op),
                                 self.main_nodes["main_file"])
        # for concatenation
        elif self.main_operations["concatenation"]:
            if self.output_file == "":
                add_node("[empty] (main)", self.main_nodes["main_file"])
            else:
                add_node("%s (main)" % basename(self.output_file),
                         self.main_nodes["main_file"])

            # Output files from secondary operations
            if secondary_op:
                for op in secondary_op:
                    if self.secondary_options["%s_file" % op]:
                        add_node("%s_%s (%s)" % (basename(self.output_file),
                                                 op, op),
                                 self.main_nodes["main_file"])
            if self.main_nodes["main_file"].is_open is False:
                self.operation_tv.toggle_node(self.main_nodes["main_file"])
        else:
            self.main_nodes["main_file"].opacity = .2

    def process_clear_options(self):

        # CORE CHANGES
        # Clear main operations
        self.main_operations = dict((op, False) for op in
                                    self.main_operations)

        # Clear active data set
        self.process_grid_wgt.ids.active_taxa_set.text = "Active taxa"
        self.process_grid_wgt.ids.active_file_set.text = "Active files"

        # Clear output formats
        self.output_formats = ["fasta"]

        # Clear filters, consensus, haplotype name and zorro suffix
        self.missing_filter_settings = [(True, 25, 50), (True, 0)]
        self.taxa_filter_settings = [None, None]
        self.codon_filter_settings = [True, True, True]
        self.hap_prefix = "Hap"
        self.zorro_suffix = ""
        self.zorro_dir = ""
        self.process_options.ids.consensus_mode.text = \
            "IUPAC" if self.sequence_types != "Protein" else "Soft mask"

        # Clear output file
        self.output_file = ""

        # Clear secondary operations
        self.secondary_operations = dict((op, False) for op in
                                         self.secondary_operations)

        # Clear secondary options
        self.secondary_options = dict((op, False) for op in
                                      self.secondary_options)

        # APP CHANGES
        # Deselect main operation
        for bt in ["conv", "conc"]:
            self.screen.ids[bt].state = "normal"
            self.screen.ids[bt].disabled = False

        # Changes in buttons with dynamic text
        # Output format
        self.process_grid_wgt.ids.conv_formatbt.text = "Fasta"
        # Output file
        self.process_grid_wgt.ids.conv.text = "Select..."
        # Zorro settings
        self.process_options.ids.zorro.background_normal = \
            "data/backgrounds/bt_process_off.png"
        self.process_options.ids.zorro.text = "Off"

        # Turn switches off
        for switch in self.secondary_operations:
            try:
                self.process_options.ids[switch].active = False
            except KeyError:
                pass

        for switch in self.secondary_options:
            try:
                self.process_options.ids[switch].active = False
            except KeyError:
                pass

    def orthology_clear_options(self):
        """
        Resets orthology search options to default values
        """

        self.ortho_dir = ""
        self.orto_export_dir = ""

        self.usearch_db = "goodProteins_db"
        self.usearch_output = "AllVsAll.out"
        self.usearch_evalue = "0.00001"

        self.ortholog_prefix = "MyGroup"
        self.group_prefix = "group"
        self.mcl_inflation = ["3"]

        self.orto_max_gene = 1
        self.orto_min_sp = 3

        self.screen.ids.usearch_threads.text = "1"

        self.ortho_search_options.ids.inflation_bt.text = "['3']"

    # ########################### PLOT SCREENS #############################

    def show_stats_toggle(self, args1, args2, active_bt, single_gene=None):
        """
        Adds a toggle widget to some Statistics plots that allow the user to
        toggle plots between the whole data set and species perspectives
        :param args1: dictionary, key must be "plt_idx", and value the
        plot index string for the species plot type
        :param args2: dictionary, key must be "plt_idx", and value the
        plot index string for the average plot type
        :param active_bt: string, identifier of the active button. Can be
        either "sp" (Species), "avg" (Average) or "gene" (Single gene)
        :param single_gene: dictionary, key must be "plt_idx", and value the
        plot index string for the single gene plot type
        """

        content = StatsToggleWgt()

        if not single_gene and active_bt != "gene":
            content.remove_widget(content.ids.gene)
            content.remove_widget(content.ids.gene_sep)
        else:
            content.gene_args = single_gene
            if active_bt == "gene":
                content.ids.gene.background_normal = \
                    "data/backgrounds/bt_focus.png"
                content.ids.gene.text = "Change gene"
            else:
                content.ids.gene.text = "Single gene"

        for wgt_name in ["avg", "sp"]:
            if wgt_name == active_bt:
                content.ids[wgt_name].state = "down"
                content.ids[wgt_name].disabled = True
            else:
                content.ids[wgt_name].state = "normal"
                content.ids[wgt_name].disabled = False

        if args2:
            content.args2 = args2
        else:
            content.ids.avg.disabled = True

        if args1:
            content.args1 = args1
        else:
            content.ids.sp.disabled = True

        self.previous_stats_toggle = content

        self.root_window.add_widget(content)

    def show_plot_toolbar(self, toolbar_type="orto"):
        """
        Adds a PlotToolbar BoxLayout to self.root_window. This is meant to
        be an auxiliary toolbar for specific operations related to plots.
        :param toolbar_type: string, determines whether an orto plot toolbar
        is displayed ('orto') or a stats plot toolbar ('stats')
        """

        if toolbar_type == "orto":
            # Determine position
            content = OrtoPlotToolbar()
        else:
            content = StatsPlotToolbar()

        self.root_window.add_widget(content)

    def show_back_bt(self):
        """
        Adds a back button to self.root_window, which will navigate to the
        previous screen. This is meant for headless plot screens
        """

        bt = BackButton()

        self.root_window.add_widget(bt)

    def dismiss_stats_toggle(self):
        """
        Removes the stats toggle widget
        """

        try:
            wgt = [x for x in self.root_window.children if
                   isinstance(x, StatsToggleWgt)][0]
            self.root_window.remove_widget(wgt)
        except IndexError:
            pass

    def dismiss_plot_wgt(self):
        """
        Removes plot widgets from the root window
        """

        try:
            for wgt in [x for x in self.root_window.children if
                        isinstance(x, OrtoPlotToolbar) or
                        isinstance(x, BackButton) or
                        isinstance(x, StatsPlotToolbar) or
                        isinstance(x, StatsToggleWgt)]:
                self.root_window.remove_widget(wgt)
        except IndexError:
            pass

    # ######################### ORTHOLOGY SCREEN ###########################

    def toggle_orto_soptions(self):
        """
        Controls the toggling of the GridLayout with the advanced options
        for the Orthology screen, Ortholog search slide
        """

        if not self.orto_search_height:
            self.orto_search_height = self.screen.ids.gl_orto_search.height

        if (self.screen.ids.adv_search_options.text ==
                "Show additional options"):
            # Add widget to main grid
            self.screen.ids.gl_orto_search.add_widget(
                self.ortho_search_options)
            # Animate widget entrance
            Animation(opacity=1, d=.5, t="in_quad").start(
                self.ortho_search_options)
            # Update button text
            self.screen.ids.adv_search_options.text = \
                "Hide additional options"

            self.screen.ids.gl_orto_search.height = \
                self.orto_search_height + sum(x.height + 5 for x in
                    self.ortho_search_options.ids.mcl_grid.children) + 30

        elif self.screen.ids.adv_search_options.text == \
                "Hide additional options":
            self.screen.ids.gl_orto_search.height = self.orto_search_height
            # Remove widget from main grid
            self.screen.ids.gl_orto_search.remove_widget(
                self.ortho_search_options)
            # Update button text
            self.screen.ids.adv_search_options.text = \
                "Show additional options"

    def dialog_search_report(self, stat_storage, groups):
        """
        Creates the dialog that reports the results of the Orthology search
        :param stat_storage: dictionary. Each entry corresponds to an
        inflation value, which will have a list as a value. The list will
        contain:

           [total_orts, species_compliant_orts, gene_compliant_orts,
           final_orts]
        """

        content = OrthoReportDialog(cancel=self.dismiss_popup)

        for inf in sorted(self.mcl_inflation):
            stats = stat_storage[inf]

            # Creating graphical report
            report_wgt = OrthoGraphicReport()

            # Setting inflation value attribute
            report_wgt.inf = inf

            # Setting total orthologs
            report_wgt.ids.total_ort.text = str(stats[0])

            # Setting gene compliant
            report_wgt.ids.gf_txt.text = str(stats[2])
            report_wgt.ids.gf_box.size_hint_x = float(stats[2]) / \
                float(stats[0])

            # Setting species compliant
            report_wgt.ids.sf_txt.text = str(stats[3])
            report_wgt.ids.sf_box.size_hint_x = float(stats[3]) / \
                float(stats[0])

            # Setting final orthologs
            report_wgt.ids.final_txt.text = str(stats[4])
            report_wgt.ids.final_box.size_hint_x = float(stats[4]) / \
                float(stats[0])

            # Adding widget to carousel
            content.ids.report_car.add_widget(report_wgt)

        content.ids.ok_bt.bind(on_release=lambda x: self.load_groups(
            groups, groups.filters))

        self.show_popup(title="Orthology search report", content=content,
                        size=(400, 470))

    def dialog_import_groups(self):

        """
        Creates filechooser dialog to select group files to be imported
        """

        content = LoadMultipleDialog(cancel=self.dismiss_popup,
                                     bookmark_init=self.bookmark_init)

        self.show_popup(title="Choose group file(s) to import",
                        content=content, size_hint=(.9, .9))

    def dialog_protein_filter(self):

        content = ProteinFilterDialog(cancel=self.dismiss_popup)

        self.show_popup(title="Protein filter settings", content=content,
                        size=(350, 200))

    def dialog_inflation(self):
        """
        Creates dialog for inflation values selection
        """

        content = InflationDialog(cancel=self.dismiss_popup)
        # Updated dialog with the selected inflation values
        for i in content.ids.inflation_bx.children:
            if i.text in self.mcl_inflation:
                i.state = "down"
            else:
                i.state = "normal"

        self.show_popup(title="MCL inflation settings", content=content,
                        size=(300, 220))

    def dialog_orto_setfilter(self, group_name):
        """
        A similar dialog to dialog_ortho_filter but for the explore screen.
        Contains an additional option of applying the specified filters
        to all group files
        :param group_name: string. name for the group object
        """

        content = OrtoSetFiltersDialog(cancel=self.dismiss_popup)
        content.group_name = group_name

        sp_list = self.ortho_groups.taxa_list[self.active_group_name]
        excluded_tx = self.ortho_groups.excluded_taxa[self.active_group_name]
        content.ids.exclude_bt.text = "{} included".format(
            len(sp_list))
        content.excluded_tx = excluded_tx

        self.show_popup(title="Set/change ortholog filters for %s file" %
                              group_name, content=content, size=(400, 390))

    def dialog_ortho_filter(self):
        """
        Creates dialog for orthology cluster filters
        """

        content = OrtoFilterDialog(cancel=self.dismiss_popup)

        self.show_popup(title="Ortholog filters", content=content,
                        size=(400, 200))

    def save_ortho_filters(self, gene_filt, sp_filt):
        """
        Save orthology clusters filters
        :param gene_filt: int. Integer for gene filter threshold
        :param sp_filt: int. Integer for species filter threshold
        """

        try:

            if int(gene_filt) < 0:
                return self.dialog_floatcheck(
                    "Maximum number of gene copies must be higher "
                    "than 0", t="error")

            self.orto_max_gene = int(gene_filt)

        except ValueError:
            return self.dialog_floatcheck(
                "Maximum number of gene copies '{}' must be a real "
                "numbers".format(gene_filt), t="error")

        try:

            if float(sp_filt) < 0:
                return self.dialog_floatcheck(
                    "Minimum number of taxa must be a positive "
                    "value", t="error")

            if 0 < float(sp_filt) < 1:
                orto_min_sp = float(sp_filt)
            else:
                orto_min_sp = int(sp_filt)

        except ValueError:
            return self.dialog_floatcheck(
                "Minimum number of taxa '{}' must be a number".format(
                    sp_filt), t="error")

        # Add check for min species. If this filter is set to a number
        # greater that the number of proteome input files (which should
        # represent a single species each) this will issue a warning.
        if self.proteome_files and \
                        int(sp_filt) > len(self.active_proteome_files):
            return self.dialog_floatcheck("Minimum number of "
                                          "species larger than the provided"
                                          " proteomes", t="warning")
        else:
            self.orto_min_sp = orto_min_sp

        self.dismiss_popup()

        return self.orto_max_gene, self.orto_min_sp

    def save_inflation(self, inflation_wgt):
        """
        Save inflation values
        :param inflation_wgt: Widget. Widget containing the inflation values
        """

        for wgt in inflation_wgt.children:
            if wgt.state == "down" and wgt.text not in self.mcl_inflation:
                self.mcl_inflation.append(wgt.text)
            elif wgt.state == "normal" and wgt.text in self.mcl_inflation:
                self.mcl_inflation.remove(wgt.text)

        self.ortho_search_options.ids.inflation_bt.text = \
            str(sorted(self.mcl_inflation))

    def save_protein_filters(self, min_len, max_stop):
        """
        Saves protein length and stop percentage filters
        :param min_len: int. Minimum sequence length
        :param max_stop: int. Maximum percentage of stop codons
        """

        # Check arguments for compliance. Must be int

        try:
            # Check if max_stop value is between 0 and 100
            if 0 < float(max_stop) > 100:
                return self.dialog_floatcheck(
                    "Maximum codon stops value must be between 0 "
                    "and 100", t="error")

            self.protein_min_len = abs(int(min_len))
            self.protein_max_stop = abs(int(float(max_stop)))
            self.dismiss_popup()
        except ValueError:
            return self.dialog_floatcheck(
                "{} and {} must be "
                "numbers".format(min_len, max_stop), t="error")

    def dialog_export_groups(self):
        """
        Dialog for group exportation.
        """

        content = ExportGroupDialog(cancel=self.dismiss_all_popups)
        self.show_popup(title="Export group as...", content=content,
                        size=(600, 360))

    def dialog_export_groups_filechooser(self, idx):
        """
        When clicking Export groups in the main dialog for group
        exportation, the user is redirected to a filechooser to choose
        the output directory and output file name, in the case of group
        exportation.
        :param idx: string. Determines the export mode: if "protein" or
        "nucleotide", it will export sequence files, if "group", it will
        export to another group file
        """

        content = SaveDialog(cancel=self.dismiss_popup,
                             bookmark_init=self.bookmark_init)

        # Set the path from previously imported groups, if any
        content.ids.sd_filechooser.path = unicode(self.orto_export_dir) if \
            self.orto_export_dir else unicode(self.home_path)

        if idx == "nucleotide" or idx == "protein":
            content.ids.txt_box.clear_widgets()
            content.ids.txt_box.height = 0
            title = "Export sequences to directory ..."

        else:
            title = "Export group file to directory ..."

        content.ids.sd_filechooser.text = idx

        self.show_popup(title=title, content=content)

    def orto_export_groups(self, export_idx, output_dir=None,
                           output_name=None):
        """
        This will handle the group exportation of the orthology screen. The
        intensive export methods will run in the background while updating
        the main process of their process.
        :param export_idx: string, with the identifier of exportation. Can
        be either 'group', 'protein' or 'dna' the database from which
        sequences will be retrieved
        :param output_dir: string. Path to output directory
        :param output_name: string, for group exportation, provide the name
        of the output filtered file
        """

        def check_process(p, dt):

            # Update dialog text
            try:
                content.ids.msg.text = shared_ns.act
            except AttributeError:
                pass

            # Update progress bar
            try:
                content.ids.pb.value = shared_ns.progress
            except AttributeError:
                pass

            # Update maximum progress bar value
            try:
                content.ids.pb.max = shared_ns.max_pb
            except AttributeError:
                pass

            if self.terminate_group_export:

                content.ids.msg.text = "Terminating on you command. " \
                                       "Hang tight!"

                # Issuing the kill order to the child thread. This signal
                # will propagate through the worker's methods and graciously
                # terminate
                shared_ns.stop = True

                # This small delay seems to fix manager shutdown issues
                time.sleep(.1)

                # Shutting down manager that provides communication between
                # parent and child thread
                manager.shutdown()

                # Close loading popup dialog
                self.dismiss_popup()

                # Unschedule the current function
                Clock.unschedule(func)

                # Join child process and exit
                p.join()
                return

            # Check process status
            if not p.is_alive():

                try:
                    if shared_ns.exception:
                        return self.dialog_floatcheck(
                            "An unexpected error occurred when exporting "
                            "orthologs. Check the app logs.", t="error")
                except Exception as e:
                    Logger.exception(e.message)

                Clock.unschedule(func)
                self.dismiss_popup()
                if not shared_ns.missed:
                    self.dialog_floatcheck(
                        "%s orthologs successfully exported" % shared_ns.good,
                        t="info")
                else:
                    self.dialog_floatcheck(
                        "%s orthologs exported. However, %s sequences "
                        "could not be retrieved! Missed sequence headers "
                        "were stored in missed_sequences.log" %
                        (shared_ns.good, shared_ns.missed), t="info")

                p.join()

        # Update orthology export directory, if necessary
        if output_dir != self.orto_export_dir:
            self.orto_export_dir = output_dir

        # Get active_group if not already defined by fetching the active
        # group button in the Explore screen
        if not self.active_group:
            group_id = [x.id for x in self.screen.ids.group_gl.children
                        if x.state == "down"][0]
            self.active_group = self.ortho_groups.get_group(group_id)

        # Let the current screen define which filters are going to be used
        # for ortholog export.
        # If the screen is ORTO_PLOT, then use the fitlers specific to the
        # current active group.
        # If the screen is ORTHOLOGY, then use the ortho groups general filters
        if self.screen.name == "Orthology":
            # Update filter values
            self.active_group.update_filters(*self.ortho_groups.filters[
                self.active_group.name])
            self.active_group.exclude_taxa(self.ortho_groups.excluded_taxa[
                self.active_group.name])

        # Set up manager and shared name space to share information
        # between background and main processes
        manager = multiprocessing.Manager()
        shared_ns = manager.Namespace()

        # Initialize stop flag to allow thread to be killed by the user
        shared_ns.stop = False

        method_store = {
            "group":
                [self.active_group.export_filtered_group,
                 [output_name, output_dir]],
            "protein":
                [self.active_group.retrieve_sequences,
                 [self.sqldb, self.protein_db, output_dir]],
            "nucleotide":
                [protein2dna.convert_group,
                 [self.sqldb, self.cds_db, self.protein_db,
                  self.active_group, self.usearch_file, output_dir]]}

        # Get method and args
        m = method_store[export_idx]

        # Remove lock from background process
        self.terminate_group_export = False

        # Create process
        p = threading.Thread(target=background_export_groups,
                             args=(m[0], shared_ns, m[1]))
        p.daemon = True
        p.start()

        # Remove any previous popups
        self.dismiss_all_popups()
        # Create loading process dialog
        content = LoadProgressDialog()
        content.ids.pb.max = self.active_group.all_compliant
        # Give functionality to cancel button
        content.ids.close_bt.bind(
            on_release=lambda i: setattr(self, "terminate_group_export",
                                         True))

        self.show_popup(title="Exporting...", content=content,
                        size=(400, 250), auto_dissmiss=False)

        # Schedule function that checks the process' pulse
        func = partial(check_process, p)
        Clock.schedule_interval(func, 1)

    def orto_report_dialog(self):
        """
        Generates a filechooser upon clicking on the "Generate full
        report" button in Orthology Explore. The filechooser will provide
        the directory where the report will be generated
        """

        content = SaveDialog(cancel=self.dismiss_popup,
                             bookmark_init=self.bookmark_init)

        # Set the path from previously imported groups, if any
        content.ids.sd_filechooser.path = unicode(self.orto_export_dir) if \
            self.orto_export_dir else self.home_path

        title = "Generate full report in directory..."

        # Remove file name text input
        content.ids.txt_box.height = 0
        content.ids.txt_box.clear_widgets()

        # Set identifier for filechooser
        content.ids.sd_filechooser.text = "orto_report"

        self.show_popup(title=title, content=content)

    def get_active_group_light(self):
        if not self.active_group:
            group_id = [x.id for x in self.screen.ids.group_gl.children if x.state == "down"][0]
            self.active_group = self.ortho_groups.get_group(group_id)
        return self.active_group

    def orto_generate_report(self, out_dir, ns=None):
        """
        Generates full orthology report on the specified directory.
        :param out_dir: string, path to directory where the report will be
        generated
        """

        # Create directory that will store figures
        fig_dir = join(out_dir, "Figures")
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)

        if ns:
            ns.files = len(MultiGroupsLight.calls)

        active_group_light = self.get_active_group_light()
        for command in MultiGroupsLight.calls:

            if ns:
                if ns.stop:
                    raise KillByUser("")
                ns.counter += 1

            data = getattr(active_group_light, command)(fig_dir)
            fig, lgd, _ = bar_plot(**data)

            plot_file = join(fig_dir, command + ".png")
            if lgd:
                fig.savefig(plot_file, bbox_extra_artists=(lgd,),
                            bbox_inches="tight", dpi=200)
            else:
                fig.savefig(plot_file, bbox_inches="tight", dpi=200)

        html = HtmlTemplate(out_dir, "Orthology report", orthology_plots)
        html.write_file()

        self.dialog_floatcheck("Orthology automatic report successfully "
                               "generated.", t="info")

    def orto_compare_groups(self, groups_objs=None, selected_groups=None):
        """
        Switches to the orthology group comparison screen and presents the
        initial plot comparing total orthologs across group files
        :param groups_objs: MultiGroupLight object. Provide only when
        updating filters in the plot screen
        :param selected_groups: list. If provided, should contain the name
        of the groups that should be plotted.
        """

        # Displays correspondence
        displays = {"total_ort": "1", "sp_ort": "2", "gn_ort": "3",
                    "final_ort": "4"}

        # Determine MultiGroupLight object
        if groups_objs:
            groups_objs = groups_objs
        else:
            groups_objs = self.ortho_groups

        # Update slider max values
        self.screen.ids.gn_spin.max = \
            max(groups_objs.max_extra_copy.values())
        self.screen.ids.sp_spin.max = groups_objs.species_number[0]

        # Update initial slider values
        if not self.screen.ids.header_content.original_filt:
            self.screen.ids.gn_spin.value = self.orto_max_gene
            self.screen.ids.sp_spin.value = self.orto_min_sp if \
                self.orto_min_sp != 3 else groups_objs.species_number[0]

        self.screen.ids.header_content.original_filt = \
            [self.screen.ids.gn_spin.value, self.screen.ids.sp_spin.value]
        self.screen.ids.header_content.excluded_taxa = \
            groups_objs.excluded_taxa[self.active_group_name]

        # Set group object for screen. This property will be used when
        # changing  which filters should be displayed in the compare
        # plot screen
        self.screen.group_obj = groups_objs

        # Get active displays
        stats = "".join([y for x, y in displays.items()
                         if self.screen.ids[x].active])

        if stats:
            # Create first comparison plot of total orthologs
            self.current_plot, self.current_lgd, self.current_table = \
                groups_objs.bar_orthologs(group_names=selected_groups,
                                          dest=self.temp_dir, stats=stats)

            # Load plot
            self.load_plot(join(self.temp_dir, "Final_orthologs.png"),
                           self.screen.ids.plot_content)

        else:
            self.screen.ids.plot_content.children[0].clear_widgets()

    def dialog_plot_check_filters(self, previous_filt, prev_tx,
                                  current_filt, current_tx):
        """
        Dialog that is triggered when exiting a plot screen and the current
        plot filters are different from the filters before entering the plot
        screen. Will ask the user if the filters should be updated for the
        remaining app or not.
        """

        # If filters are the same, do nothing and return
        if (previous_filt == current_filt) and (prev_tx == current_tx):
            self.go_previous_screen()
            return

        content = PlotChangeFilters()

        content.ids.prev_gn.text = str(previous_filt[0])
        content.ids.cur_gn.text = str(current_filt[0])

        content.prev_filt = previous_filt
        content.prev_tx = prev_tx
        content.cur_filt = current_filt
        content.cur_tx = current_tx

        if previous_filt[0] != current_filt[0]:
            content.ids.prev_gn.color = self._red
            content.ids.cur_gn.color = self._red

        content.ids.prev_min_tx.text = str(previous_filt[1])
        content.ids.cur_min_tx.text = str(current_filt[1])

        if previous_filt[1] != current_filt[1]:
            content.ids.prev_min_tx.color = self._red
            content.ids.cur_min_tx.color = self._red

        content.ids.prev_ex_tx.text = str(len(prev_tx))
        content.ids.cur_ex_tx.text = str(len(current_tx))

        if sorted(prev_tx) != sorted(current_tx):
            content.ids.prev_ex_tx.color = self._red
            content.ids.cur_ex_tx.color = self._red

        self.show_popup(title="", content=content, size=(520, 400),
                        separator_color=(0, 0, 0, 0))

    def orto_show_plot(self, active_group, plt_idx, filt=None,
                       exclude_taxa=None):
        """
        Loads a orto_plot screen for orthology graphical exploration based
        on the plot index. This method can be called in three ways:

        ..: Orthology Explore screen button, which generates the plot for
            the first time and sets initial attributes for plot screen
            header. This only uses the plt_idx argument.
        ..: Update button in plot screen, used after gene/species filters
            have been changed. This uses the plt_idx and filt arguments.
        ..: Taxa filter button in plot screen, used after changing the taxa
            that should be included in the plot analysis. This uses the
            plt_idx and exclude_taxa arguments

        :param active_group: Group Object.
        :param plt_idx: string, id of the plot in plt_method to issue the
        appropriate method
        :param filt: list, contains the gn and sp filters for group object,
        respectively. If none, this will set a new plot and all plot screen
        attributes. Else, it will only update the plot image and ortholog
        statistics.
        :param exclude_taxa: list, each element should be a taxon name to be
        excluded from the plot
        """

        # Set active group
        if active_group:
            self.active_group = active_group
        else:
            active_group = self.active_group

        # Exclude taxa, if any
        if exclude_taxa:
            # If all taxa were excluded issue a warning and do nothing more
            if set(exclude_taxa) == set(self.active_group.species_list):
                return self.dialog_floatcheck(
                    "At least one taxon must be included.",
                    t="warning")

        # Setting filters for the first time
        if not filt and not exclude_taxa:
            self.screen.ids.gn_spin.value = self.active_group.max_extra_copy
            self.screen.ids.sp_spin.value = 1
            self.screen.ids.header_content.original_filt = \
                [self.active_group.max_extra_copy, 1]
            self.screen.ids.header_content.root_filt = \
                [active_group.gene_threshold, active_group.species_threshold]
            self.screen.ids.header_content.root_excluded = \
                active_group.excluded_taxa

            self.active_group.exclude_taxa([])
            self.active_group.update_filters(
                self.active_group.max_extra_copy, 1, True)

        self.run_in_background(
            get_orto_data,
            self.orto_write_plot,
            [active_group, plt_idx, filt, exclude_taxa],
            [plt_idx]
        )

    def orto_write_plot(self, plot_data, plt_idx):
        """
        Provide with the data structure and a plt_idx string identifier, this
        function will create the plot file and update the plot screen.
        :param plot_data: list/np array, data structure to be used in plot
        construction
        :param plt_idx: string, identification string of the plot.
        """

        if not plot_data["data"]:
            self.screen.ids.plot_content.children[0].clear_widgets()
            return self.dialog_floatcheck("There is not data to plot with"
                                          "the current filters", t="error")

        # Update excluded taxa attribute
        self.screen.ids.header_content.excluded_taxa = \
            self.active_group.excluded_taxa

        # Update orthology filter attributes
        self.screen.ids.header_content.original_filt = \
            [self.active_group.gene_threshold,
             self.active_group.species_threshold]

        # Update orthology filter max values
        self.screen.ids.gn_spin.max = self.active_group.max_extra_copy
        self.screen.ids.sp_spin.max = len(self.active_group.species_list)
        # Update slider values if they are outside bounds
        if self.screen.ids.gn_spin.value > self.screen.ids.gn_spin.max:
            self.screen.ids.gn_spin.value = self.screen.ids.gn_spin.max

        if (self.screen.ids.sp_spin.value >
                len(self.active_group.species_list)):

            self.screen.ids.sp_spin.value = \
                len(self.active_group.species_list)

        # Set the current plt_idx for update reference
        self.screen.ids.header_content.plt_idx = plt_idx

        # Update summary attributes
        self.screen.ids.orto_sum.text = "[size=26][color=71c837ff]%s" \
            "[/color][/size][size=13]/[color=ff5555ff]%s[/color]" \
            "[/size]" % \
            (str(self.active_group.all_compliant),
             str(self.active_group.all_clusters))

        self.screen.ids.taxa_sum.text = "[size=26][color=71c837ff]%s" \
            "[/color][/size][size=13]/[color=ff5555ff]%s[/color]" \
            "[/size]" % \
            (len(self.active_group.species_list),
             len(self.active_group.excluded_taxa) +
             len(self.active_group.species_list))

        # Call corresponding method and catch plot object
        self.current_plot, self.current_lgd, self.current_table = \
            self.orto_plt_method[plt_idx][0](**plot_data)

        # Get plot path
        plot_path = join(self.temp_dir, self.orto_plt_method[plt_idx][1])

        # Save plot figure
        if self.current_lgd:
            self.current_plot.savefig(
                plot_path,
                bbox_extra_artists=(self.current_lgd,),
                bbox_inches="tight", dpi=200)
        else:
            self.current_plot.savefig(
                plot_path,
                bbox_inches="tight", dpi=200)

        # Load plot
        self.load_plot(plot_path,
                       self.screen.ids.plot_content)

    @staticmethod
    def load_plot(file_path, scatter_wgt):
        """
        Loads a new plot into a ScatterLayout. This will clear all previous
        content and load a new image based on the file_path argument.
        This assumes that the current screen is a plot related screen.
        :param file_path: string. Path to the image to be loaded
        :param scatter_wgt: ScatterLayout object, where the plot is to be
        loaded
        :return:
        """

        # Clear previous content
        scatter_wgt.children[0].clear_widgets()

        # Add content
        img_wgt = Image(source=file_path, nocache=True)
        scatter_wgt.children[0].add_widget(img_wgt)

        # Reset position and scale of Scatter
        scatter_wgt.scale = 1
        scatter_wgt.pos = (0, 0)

    def dialog_set_exclude_orto_taxa(self):
        """
        Creates a dialog for the orthology change filters option.
        """

        def set_excluded():
            excluded_taxa = [x.text for x in
                             self._subpopup.content.ids.rev_inlist.children if
                             x.state == "normal"]

            if len(excluded_taxa) == \
                    len(self._subpopup.content.ids.rev_inlist.children):
                return self.dialog_floatcheck("At least one taxon must be "
                                              "included.", t="error")

            self._popup.content.excluded_tx = excluded_taxa
            self._popup.content.ids.exclude_bt.text = "{} included".format(
                len(sp_list) + len(ex_list) - len(excluded_taxa))
            self.dismiss_subpopup()

        sp_list = self.ortho_groups.taxa_list[self.active_group_name]
        ex_list = self.ortho_groups.excluded_taxa[self.active_group_name]

        content = InputList(cancel=self.dismiss_subpopup)
        # Show the select/deselect all buttons
        content.ids.selection_bx.height = 50

        # Add button for each taxon
        for taxon in sorted(sp_list + ex_list):
            bt = TGToggleButton(text=taxon, height=30, state="down")
            # deselect button if taxa is excluded
            if self._popup.content.excluded_tx is not None:
                if taxon in self._popup.content.excluded_tx:
                    bt.state = "normal"
            elif taxon in ex_list:
                bt.state = "normal"
            # Add button to list
            content.ids.rev_inlist.add_widget(bt)

        content.ids.ok_bt.bind(on_release=lambda x: set_excluded())

        self.show_popup(title="Included taxa", content=content,
                        size_hint=(.3, .8), popup_level=2)

    def dialog_exclude_orto_taxa(self, plt_idx):

        content = InputList(cancel=self.dismiss_popup)
        # Show the select/deselect all buttons
        content.ids.selection_bx.height = 50

        # Add button for each taxon
        for taxon in sorted(self.active_group.species_list +
                            self.active_group.excluded_taxa):
            bt = TGToggleButton(text=taxon, height=30, state="down")
            # deselect button if taxa is excluded
            if taxon in self.active_group.excluded_taxa:
                bt.state = "normal"
            # Add button to list
            content.ids.rev_inlist.add_widget(bt)

        # Add bindings to Ok button
        content.ids.ok_bt.bind(on_release=lambda j:
            self.run_in_background(
                get_active_group,
                self.orto_show_plot,
                [self.ortho_groups, self.active_group,
                 str(self.active_group_name)],
                [str(plt_idx), [int(self.screen.ids.gn_spin.value),
                 int(self.screen.ids.sp_spin.value)], [x.text for x in
                     content.ids.rev_inlist.children if x.state == "normal"]],
                False,
                cancel=True))

        self.show_popup(title="Included taxa", content=content,
                        size_hint=(.3, .8))

    def orto_change_state(self):
        """
        Toggle selection or deselection of group checkboxes
        """

        # In case all are already selected, deselect all
        if len(self.screen.ids.group_check.children) == \
                len([x for x in self.screen.ids.group_check.children if
                     x.active]):
            for chk in self.screen.ids.group_check.children:
                chk.active = False
        else:
            for chk in self.screen.ids.group_check.children:
                chk.active = True

    def orto_check_state(self):
        """
        sets the "Compare" button disabled attribute according to the
        number of active check boxes
        """

        # Change Compare button state
        if (len([x for x in self.screen.ids.group_check.children if
                 x.active]) >= 2):
            self.screen.ids.compare_group_bt.disabled = False
        else:
            self.screen.ids.compare_group_bt.disabled = True

        # Check master check box state
        if all([x.active for x in self.screen.ids.group_check.children]):
            self.screen.ids.master_chk.active = True
        else:
            self.screen.ids.master_chk.active = False

    def load_groups(self, groups_obj, default_filters):
        """
        Loads the group files generated by the Orthology search or manually
        imported into the app. This method only accepts MultiGroup objects.
        Loading groups from files is more computationally intensive and
        should be done in the background using the run_in_background and
        load_group_files methods.

        :param groups_obj: MultiGroup object
        :param default_filters: Tuple. With filters for MultiGroups object
        """

        if groups_obj:

            if self.ortho_groups:
                self.ortho_groups.add_multigroups(groups_obj)
                new_groups = [x for x in groups_obj.groups if
                              x not in self.ortho_groups.duplicate_groups]
            else:
                self.ortho_groups = groups_obj
                self.ortho_groups.filters = default_filters
                new_groups = self.ortho_groups.groups.keys()

            self.ortho_group_files.extend(list(groups_obj.groups.keys()))

            # Check if any group file is duplicate. If so, issue a warning
            if groups_obj.duplicate_groups or groups_obj.bad_groups:

                msg = ""

                if groups_obj.duplicate_groups:
                    msg += "The following group files were found to be" \
                           " duplicate and were not loaded:\n\n[b]%s[/b]" %\
                           "\n".join(basename(x) for x in
                                     groups_obj.duplicate_groups)

                    # Reset duplicate_groups attribute
                    self.ortho_groups.duplicate_groups = []

                if groups_obj.bad_groups:
                    msg += "The following group files could not be parsed "\
                           "as group files:\n\n[b]%s[/b]" % \
                           "\n".join(basename(x) for x in
                                     groups_obj.bad_groups)

                    # Reset bad_groups attribute
                    self.ortho_groups.bad_groups = []

                self.dialog_warning("Invalid group files detected", msg)

            if groups_obj.groups:

                # Removes "No groups loaded" button if it still exists
                try:
                    self.screen.ids.group_check.remove_widget(
                        self.screen.ids.no_bt)
                except ReferenceError:
                    pass

                # Populate the app gridlayout with group buttons
                for gname in sorted(groups_obj.groups):

                    # If group name contains full path, get only file name
                    gname_short = basename(gname)

                    if gname not in self.ortho_groups.duplicate_groups:
                        # Create check box for multiple group selection
                        chk = CheckBox(id=gname, size_hint=(.1, None),
                                       height=30)
                        chk.bind(active=lambda i, y:
                            self.orto_check_state())

                        self.screen.ids.group_check.add_widget(chk)

                        # Create group button
                        bt = ToggleButton(text=gname_short, id=gname,
                            group="group_bts", size_hint_y=None,
                            height=30, shorten=True, shorten_from="right",
                            halign="center", bold=True,
                            background_down=join("data", "backgrounds",
                                                 "bt_process.png"),
                            background_normal=join("data", "backgrounds",
                                                   "bt_process_off.png"),
                            background_disabled_down=join("data",
                                                          "backgrounds",
                                                          "bt_process.png"),
                            disabled_color=(1, 1, 1, 1))

                        # Apparently I need to use partial instead of lambda
                        # in order to provide a diferent group object as
                        # argument. Using lambda will overwrite the group
                        # objects of all buttons with the last group of the
                        # iteration.
                        bt.bind(on_release=partial(self.orthology_card,
                                                   gname))

                        # Add box to gridlayout
                        self.screen.ids.group_gl.add_widget(bt)

                        # Create removal button
                        x_bt = Button(size_hint=(None, None), width=30,
                            height=30, id=gname, border=(0, 0, 0, 0),
                            background_normal=join("data", "backgrounds",
                                                   "remove_bt.png"),
                            background_down=join("data", "backgrounds",
                                                 "remove_bt_down.png"))
                        x_bt.bind(on_release=partial(
                            self.check_action,
                            "Are you sure you want to remove this group?",
                            self.remove_groups))

                        self.screen.ids.group_rm.add_widget(x_bt)

                # If last group name contains a directory, set it as the
                # default export dir
                try:
                    path = dirname(gname)
                    if os.path.exists(path):
                        self.orto_export_dir = path
                except AttributeError:
                    pass

                self.run_in_background(
                    orto_update_filters,
                    self.orthology_card,
                    [self.ortho_groups, None, None, [],
                     new_groups, True],
                    None,
                    False,
                    msg="Setting up filters...",
                    cancel=False)

    def orthology_card(self, group_name=None, bt=None):
        """
        Generates the descriptive cards with general information for a group
        file.
        :param group_name
        :param bt: ToggleButton instance
        """

        # If no group button is active, dispatch the first
        # The comparison is done with __class__.__name__ to avoid
        # issues of the isinstance() method comparing the same class
        # instance from different imports (trifusion.ortho vs. ortho)
        if group_name and group_name.__class__.__name__ == "MultiGroupsLight":
            try:
                self.ortho_groups = group_name
            except Exception as e:
                Logger.exception(e.message)

        if group_name and not \
                group_name.__class__.__name__ == "MultiGroupsLight":
            pass

        elif (not [x for x in self.screen.ids.group_gl.children
                   if x.state == "down"] and
                self.screen.ids.group_gl.children):
            self.screen.ids.group_gl.children[-1].state = "down"
            self.screen.ids.group_gl.children[-1].disabled = True
            group_name = self.screen.ids.group_gl.children[-1].id

        else:
            group_name = [x.id for x in self.screen.ids.group_gl.children
                          if x.state == "down"][0]

        self.active_group_name = group_name

        # Create desired behaviour for group toggle buttons
        if bt:
            self.toggle_groups(bt)

        # Get statistics from group object
        stats = self.ortho_groups.groups_stats[group_name]["stats"]

        # Create cards
        cards = DescriptionBox(opacity=0)

        # Populate card with group information
        cards.prot_txt = str(stats[1])
        cards.ortholog_txt = str(stats[0])
        cards.taxa_txt = str(len(self.ortho_groups.groups_stats[group_name]
                                 ["species"]))
        cards.group_name = basename(group_name)

        # Create gauge plots, if there are any filtered groups
        if self.ortho_groups.filters[group_name][0] or \
                self.ortho_groups.filters[group_name][1] or \
                self.ortho_groups.filters[group_name] == (0, 0):

            # Create species filter plot and add to box
            sp_filter_plot = GaugePlot()
            sp_filter_plot.txt = "After species filter"
            sp_filter_plot.proportion = float(stats[3]) / float(stats[0])
            sp_filter_plot.ortholog_num = "{0:,}".format(stats[3])
            cards.ids.gauge_bx.add_widget(sp_filter_plot)

            # Create gene filter plot and add to box
            gn_filter_plot = GaugePlot()
            gn_filter_plot.txt = "After gene filter"
            gn_filter_plot.proportion = float(stats[2]) / float(stats[0])
            gn_filter_plot.ortholog_num = "{0:,}".format(stats[2])
            cards.ids.gauge_bx.add_widget(gn_filter_plot)

            # Create final ortholog plot
            final_ortholog_plot = GaugePlot()
            final_ortholog_plot.txt = "Final orthologs"
            final_ortholog_plot.proportion = float(stats[4]) / \
                float(stats[0])
            final_ortholog_plot.ortholog_num = "{0:,}".format(stats[4])
            cards.ids.gauge_bx.add_widget(final_ortholog_plot)

        # Add button to generate full report
        full_rep_bt = Button(text="Generate full report", bold=True,
                             background_normal=join("data", "backgrounds",
                                                    "check_ok.png"),
                             background_down=join("data", "backgrounds",
                                                  "check_cancel.png"),
                             size=(180, 40), size_hint=(None, None),
                             font_size=17)
        full_rep_bt.bind(on_release=lambda i: self.orto_report_dialog())
        # Anchor layout that will hold full report button
        anc = AnchorLayout(anchor_y="center", achor_x="center",
                           size_hint_y=None, height=100)
        anc.add_widget(full_rep_bt)
        cards.add_widget(anc)

        self.screen.ids.card_gl.clear_widgets()

        # Add card
        Clock.schedule_once(
            lambda y: self.screen.ids.card_gl.add_widget(cards), .3)
        Clock.schedule_once(
            lambda y: Animation(opacity=1, d=.3, t="out_quart").start(
                cards), .3)

    # ########################## POPUP OPS #################################

    def show_popup(self, title, content, size_hint=(.9, .9), size=None,
                   separator_color=None, close_bt=None, auto_dissmiss=True,
                   popup_level=1, border_color=None):
        """
        General purpose method to create a popup widget

        :param str title: Title of the popup
        :param widget content: The contents of the popup widget
        :param tuple size_hint: Size hint for the widget
        :param tuple size: The absolute size for the popup. If this
        argument is used, the size_hint will be ignored
        :param list separator_color: List with rgb color of popup separator
        :param bool close_bt: If True, Adds a (X) button to the top
        right of the popup that dismisses the popup
        :param bool auto_dismiss: Sets whether the popup should be
        automatically dismissed on pressing Esc
        :param int popup_level: Determines the level of the popup. Two
        values are currently supported. 1 assigns the popup to _popup; 2
        assigns the popup to _subpopup
        """

        # This prevents defining a mutable argument.
        if not separator_color:
            separator_color = [47 / 255., 167 / 255., 212 / 255., 1.]

        if not border_color:
            border_color = separator_color

        # Ignore size_hint is absolute size is provided
        if size:
            popup_obj = CustomPopup(title="[b]%s[/b]" % title,
                content=content, size=size,
                size_hint=(None, None),
                auto_dismiss=auto_dissmiss,
                separator_color=separator_color,
                border_color=border_color,
                separator_height=0.6,
                title_color=separator_color,
                background="data/backgrounds/transparent.png")
        else:
            popup_obj = CustomPopup(title="[b]%s[/b]" % title,
                content=content, size_hint=size_hint,
                separator_color=separator_color,
                border_color=border_color,
                auto_dismiss=auto_dissmiss,
                title_color=separator_color,
                separator_height=0.6,
                background="data/backgrounds/transparent.png")

        popup_obj.bind(on_dismiss=self.dismiss_close_float)

        popup_obj.open()

        if popup_level == 1:
            self._popup = popup_obj
        else:
            self._subpopup = popup_obj

        if close_bt:
            pos = ((self.root.width / 2) + (self._popup.size[0] / 2) - 25,
                   (self.root.height / 2) + (self._popup.size[1] / 2) - 25)
            rm_wgt = CloseFloat(pos=pos, opacity=1)
            rm_wgt.bind(on_release=self.dismiss_popup if popup_level == 1 else
                        self.dismiss_subpopup)
            self.root_window.add_widget(rm_wgt)

    def dismiss_all_popups(self, *args):
        """
        Method that force closes all popups in thre screen
        """

        for wgt in (x for x in self.root_window.children
                    if isinstance(x, CustomPopup)):
            wgt.dismiss()

        try:
            rm_wgt = [x for x in self.root_window.children if
                      isinstance(x, CloseFloat)][0]
            self.root_window.remove_widget(rm_wgt)
        except IndexError:
            pass

    def dismiss_close_float(self, *args):
        """
        Use to automatically remove any possibly CloseFloat widgets from the
        root window, when a Popup is automatically dismissed
        """

        try:
            rm_wgt = [x for x in self.root_window.children if
                      isinstance(x, CloseFloat)][0]
            self.root_window.remove_widget(rm_wgt)
        except IndexError:
            pass

    def dismiss_popup(self, *args):
        """
        General purpose method to close popups from the screen
        """
        if self._popup:
            self._popup.dismiss()
            try:
                rm_wgt = [x for x in self.root_window.children if
                          isinstance(x, CloseFloat)][0]
                self.root_window.remove_widget(rm_wgt)
            except IndexError:
                pass

    def dismiss_subpopup(self, *args):
        """
        General purpose method to close sub-popups from the screen
        """
        self._subpopup.dismiss()

        try:
            rm_wgt = [x for x in self.root_window.children if
                      isinstance(x, CloseFloat)][0]
            self.root_window.remove_widget(rm_wgt)
        except IndexError:
            pass

    def dismiss_exit(self, *args):
        """
        Dismiss function for exit popup
        """

        if self._exit_popup:
            self._exit_popup.dismiss()

    # ########################## PROCESS SCREEN ############################

    def save_file(self, path, file_name=None, idx=None, auto_close=True,
                  suffix=""):
        """
        Adds functionality to the save button in the output file chooser. It
        gathers information on the specified path through filechooser, file
        name through textinput and the widget text when called.

        For now, only one main output file can be provided, so the its path
        is stored in a string attribute.

        :param path: string. complete path
        :param file_name: string. file name only
        :param idx: string. An id of where the filechooser is calling. This
        allows the addition of custom behaviours for different dialogs
        :param auto_close: Boolean. When True, the dismiss_popup method is called in
        the end of the function
        :param suffix: string. Allow the addition of a suffix to the file name
        or to be save in an attribute
        """

        if idx == "main_output":
            if self.main_operations["concatenation"]:
                # Adds output file to storage
                try:
                    self.output_file = join(unicode(path, "utf8"),
                                            file_name)
                except TypeError:
                    self.output_file = join(path, file_name)
                self.output_dir = path
                # Renames the output file button text
                self.process_grid_wgt.ids.conv.text = file_name

            else:
                try:
                    self.output_dir = unicode(path, "utf8")
                except TypeError:
                    self.output_dir = path
                self.process_grid_wgt.ids.conv.text = basename(path)
                self.conversion_suffix = suffix

        elif idx == "ortho_dir":
            self.ortho_dir = path
            self.screen.ids.orto_dir.text = basename(path)

        elif idx == "protein_db":
            self.protein_db = path[0]

        elif idx == "cds_db":
            self.cds_db = path

        elif idx == "orto_export_dir":
            self.orto_export_dir = path

        elif idx == "ima2_popfile":
            self.ima2_options[0] = path

        elif idx == "zorro_dir":
            self.zorro_dir = path

        elif idx == "mcl_fix":

            # The mcl executable requires the .exe extension in windows to
            # work properly
            if sys.platform in ["win32", "cygwin"]:
                mcl_string = "mcl.exe"
            else:
                mcl_string = "mcl"

            # Copy new mcl_file to app_dir
            if self._check_exec(path, "mcl"):
                # Copy mcl executable to app dir
                shutil.copyfile(path, join(self.user_data_dir, mcl_string))
                self.mcl_file = join(self.user_data_dir, mcl_string)
                # Make it executable
                st = os.stat(self.mcl_file)
                os.chmod(self.mcl_file, st.st_mode | stat.S_IEXEC)
                self.dismiss_all_popups()
                # Set the orthology fields to green
                self.ortho_search_options.ids.mcl_check.background_src = \
                    "data/backgrounds/green_noise.png"
                self.ortho_search_options.ids.mcl_check_ico.icon_src = \
                    "data/backgrounds/check_icon.png"
                self.ortho_search_options.ids.mcl_check_ico.icon_size = \
                    (20, 16)
                self.ortho_search_options.ids.mcl_check_lbl.text = \
                    "MCL is installed and reachable"
                self.ortho_search_options.ids.mcl_fix_box.clear_widgets()
            else:
                self.dialog_floatcheck("The provided MCL executable does "
                                       "not seem to be correct.", t="error")

        elif idx == "usearch_fix":
            # Copy new usearch_file to app_dir
            shutil.copyfile(path, join(self.user_data_dir, "usearch"))
            self.usearch_file = join(self.user_data_dir, "usearch")
            # Make it executable
            st = os.stat(self.usearch_file)
            os.chmod(self.usearch_file, st.st_mode | stat.S_IEXEC)
            if self._check_exec(self.usearch_file, "usearch"):
                # Copy mcl executable to app dir
                self.dismiss_all_popups()
                # Set the orthology fields to green
                self.ortho_search_options.ids.usearch_check.background_src =\
                    "data/backgrounds/green_noise.png"
                self.ortho_search_options.ids.usearch_check_ico.icon_src = \
                    "data/backgrounds/check_icon.png"
                self.ortho_search_options.ids.usearch_check_ico.icon_size = \
                    (20, 16)
                self.ortho_search_options.ids.usearch_check_lbl.text = \
                    "USEARCH is installed and reachable"
                self.ortho_search_options.ids.usearch_fix.clear_widgets()
            else:
                os.remove(self.usearch_file)
                self.usearch_file = ""
                self.dialog_floatcheck("The provided USEARCH executable"
                                       " does "
                                       "not seem to be correct.", t="error")

        if auto_close:
            self.dismiss_popup()

    def save_format(self, value):
        """
        Method that stores the output formats specified through the formats
        dialog in the Process screen.

        The active formats are stored in a self.output_formats list

        :param value: widget object.
        """

        # Add the selected output formats to the storage list and remove
        # deselected formats if they have been selected before
        for idx, wgt in value.ids.items():
            if isinstance(wgt, ToggleButton):
                if wgt.state == "down" and idx not in self.output_formats:
                    self.output_formats.append(idx)
                elif wgt.state == "normal" and idx in self.output_formats:
                    self.output_formats.remove(idx)

        if not self.output_formats:
            return self.dialog_floatcheck(
                "Please choose at least one output format",
                t="warning")

        self.dismiss_popup()

        # If IMa2 is among the chosen ouptut formats, issue its
        # additional options dialog.
        if "ima2" in self.output_formats:
            # Only show dialog when the additional options have not been set
            if None in self.ima2_options:
                self.dialog_ima2_extra()

    def save_gapfilter(self, filter_act, within_chk, among_chk, gap_val,
                       mis_val, min_tx_val):
        """
        Stores the information of the FilterDialog
        :param filter_act: Boolean, whether the filter is active (True) or
        not (False)
        :param within_chk: Boolean, whether the within alignment filter
        is active (True) or not (False)
        :param among_chk: Boolean, whether the multiple alignments filter
        is active (True) or not (False)
        :param gap_val: integer, proportion of gap threshold
        :param mis_val: integer, proportion of missing data threshold
        :param min_tx_val: integer, proportion of minimum taxa
        representation
        """

        if filter_act and not within_chk and not among_chk:
            return self.dialog_floatcheck("Warning: No filters were "
                "checked. Specify at least one to active the filter.",
                                          t="error")

        self.secondary_options["gap_filter"] = filter_act

        # Save only when the filter is set to active
        if filter_act:
            self.missing_filter_settings = [(within_chk, gap_val, mis_val),
                                            (among_chk, min_tx_val)]

        self.dismiss_popup()

    def save_taxafilter(self, filter_act, filter_mode, taxa_group):
        """
        Stores the information of the taxa filter dialog of the process
        screen
        :param filter_mode: string, determines the taxa filtering mode. Can
        be either 'Contain' or 'Exclude'
        :param taxa_group: string, with the name of the taxa group that
        should be present in the taxa_group attribute
        """

        # Check if a taxa group has been selected. If not issue a warning
        # and do not close popup
        if taxa_group == "No group selected":
            return self.dialog_floatcheck("Warning: A taxa group must be "
                "specified to activate the filter", t="error")

        self.secondary_options["taxa_filter"] = filter_act

        # Save only when the filter is set to active
        if filter_act:
            self.taxa_filter_settings = [filter_mode, taxa_group]

        self.dismiss_all_popups()

    def save_codonfilter(self, filter_act, position_list):
        """
        Stores the information of the alignment filter dialog of the process
        screen
        :param filter_act: Boolean, whether the filter is active (True) or
        not (False)
        :param position_list: A list of three elements, containing which
        positions should be saved (True) or filtered (False)
        """

        self.secondary_options["codon_filter"] = filter_act

        if filter_act:
            self.codon_filter_settings = position_list

        self.dismiss_popup()

    def save_variationfilter(self, filter_act, *pargs):
        """
        Stores the information of the variation filter dialog.
        :param filter_act: Boolean, whether the filter is active (True) or
        not (False)
        :param pargs, list of tuples, with first element being the
        checkbox active property and the second element the value of said
        parameter
        """

        # Check if any parameter was checked. If not, issue warning and
        # do not close popup
        if not any([p[0] for p in pargs]) and filter_act:
            return self.dialog_floatcheck("Warning: No parameters were "
                "checked. Specify at least one to activate the filter.",
                                          t="error")

        # Check if all checked parameters have values. If not,
        # issue warning and do not close popup
        for p in pargs:
            if p[0] and not p[1]:
                return self.dialog_floatcheck("Warning: Checked "
                    "parameters must have values specified.", t="error")

        # Check for consistency between min and max values
        for i in range(0, len(pargs), 2):
            try:
                vals = (int(pargs[i][1]), int(pargs[i + 1][1]))
                if vals[0] and vals[1] and vals[0] > vals[1]:
                    return self.dialog_floatcheck("Warning: Minimum values "
                        "cannot be greater than maximum values", t="error")
            except ValueError:
                pass

        self.secondary_options["variation_filter"] = filter_act

        if filter_act:
            for i, p in enumerate(pargs):
                p = int(p[1]) if p[0] else None
                self.variation_filter[i] = p

        self.dismiss_popup()

    def save_ima2_opts(self, pop_string, mutation, inheritance):
        """
        Check each text input parameter given by the user and, if they
        checkout, saves them.
        :param pop_string: string, with population tree
        :param mutation:  string, with mutation model
        :param inheritance: string, with inheritance scalar
        """

        for i, msg in zip([self.ima2_options[0], pop_string, mutation,
                           inheritance],
                          ["population file",
                           "population tree string",
                           "mutation model",
                           "inheritance scalar"]):
            if not i:
                return self.dialog_floatcheck("Please specify a "
                                              "{}".format(msg), t="warning")

        self.ima2_options[1] = pop_string
        self.ima2_options[2] = mutation
        self.ima2_options[3] = inheritance

        self.dismiss_all_popups()

    def dialog_ima2_extra(self, popup_level=1):
        """
        Dialog with extra options for IMa2 output format
        """

        if popup_level == 1:
            content = IMa2Extra(cancel=self.dismiss_popup)
        else:
            content = IMa2Extra(cancel=self.dismiss_subpopup)

        if self.ima2_options[0]:
            content.ids.popfile.text = basename(self.ima2_options[0])

        for i, idx in zip(self.ima2_options[1:], ["pop_string",
                                                  "mutation",
                                                  "inheritance"]):
            if i:
                content.ids[idx].text = i

        title = "IMa2 additional options"

        self.show_popup(title=title, content=content, size=(380, 310),
                        popup_level=popup_level)

    def dialog_nexus_extra(self):
        """
        Dialog with extra options for nexus output format
        """

        content = NexusExtra(cancel=self.dismiss_subpopup)

        content.ids.nexus_check.active = self.use_nexus_partitions
        content.ids.nexus_model_check.active = self.use_nexus_models

        self.show_popup(title="Nexus extra options", content=content,
                        size=(500, 210), popup_level=2)

    def dialog_phylip_extra(self):
        """
        Dialog with extra options for phylip output format
        """

        content = PhylipExtra(cancel=self.dismiss_subpopup)

        content.ids.part_check.active = self.create_partfile
        content.ids.trunc_names_check.active = self.phylip_truncate_name

        self.show_popup(title="Phylip extra options",
                        content=content, size=(400, 230),
                        popup_level=2)

    def dialog_fasta_extra(self):
        """
        Dialog with extra options for fasta output format
        """

        content = FastaExtra(cancel=self.dismiss_subpopup)

        content.ids.ldhat_check.active = self.ld_hat

        self.show_popup(title="Fasta extra options", content=content,
                        size=(400, 160), popup_level=2)

    def save_zorro_settings(self, suffix):
        """
        Handles the information provided by the user in the ZorroDialog
        :param suffix: string, suffix of the ZORRO files
        """

        # Only perform file checking is switch is active
        if self._popup.content.ids.zorro_switch.active:

            # Check if the zorro files exist
            for f in self.active_file_list:
                # Check if a specific directory for the ZORRO files has been
                # specified. If so, get only the file name of the
                # active_file_list and use the zorro_dir to provide the path
                if self.zorro_dir:
                    f = splitext(f.split(sep)[-1])[0]
                    f = "{}{}.txt".format(join(self.zorro_dir, f), suffix)
                # If zorro_dir is not specified, assume that the ZORRO files
                # are in he same working directory of the input files
                else:
                    f = sep.join(f.split(".")[0:-1])
                    f = "%s%s.txt" % (f, suffix)
                if not os.path.isfile(f):
                    return self.dialog_floatcheck(
                        "File %s does not exist" % f, t="error")

        self.update_process_switch(
            "zorro", self._popup.content.ids.zorro_switch.active)

        if self.secondary_options["zorro"]:
            # Save auxiliary files suffix
            self.zorro_suffix = suffix
            # Change background and text of the Zorro button in process
            # options
            self.process_options.ids.zorro.background_normal = \
                "data/backgrounds/bt_process.png"
            self.process_options.ids.zorro.text = "Active"

        else:
            self.process_options.ids.zorro.background_normal = \
                "data/backgrounds/bt_process_off.png"
            self.process_options.ids.zorro.text = "OFF"

        self.dismiss_popup()

    def save_reverseconc_partfile(self, partfile):
        """
        Save partition file for reverse concatenation an update button in
        reverse concatenation dialog
        :param partfile: string. Path to partition file
        """

        self.dismiss_subpopup()

        self.partitions_file = partfile

        # Update partition file button
        self._popup.content.ids.part_file.text = basename(partfile)
        self._popup.content.ids.part_file.background_normal = \
            "data/backgrounds/bt_process.png"

        # Deactivate user defined partitions
        self._popup.content.ids.use_parts.state = "normal"
        self.use_app_partitions = False

        # Change background of used defined partition button
        self._popup.content.ids.use_parts.state = "normal"

    def dialog_floatcheck(self, text, t, dl=5):
        """
        Creates a floating label with informative text on the right upper
        corner of the app. This is used for showing errors, warnings and
        general informative messages that fade in and fade out after a time
        :param text: string, text to appear in the label
        :param t: string, with type of check. Can be either error or info
        :param dl: integer, duration of the float
        """

        def fade_in():
            Animation(opacity=1, d=.5, t="out_quart").start(check_wgt)
            Animation(opacity=1, d=.5, t="out_quart").start(rm_wgt)
            Animation(opacity=1, d=.5, t="out_quart").start(ico_wgt)

        def fade_out():
            # Get height of warning
            h = check_wgt.height
            pos_h = check_wgt.pos[1]

            for wgt in [check_wgt, rm_wgt, ico_wgt]:
                Animation(opacity=0, d=.2, t="out_quart").start(wgt)

            Clock.schedule_once(
                lambda dt: self.root_window.remove_widget(check_wgt), .2)
            Clock.schedule_once(
                lambda dt: self.root_window.remove_widget(rm_wgt), .2)
            Clock.schedule_once(
                lambda dt: self.root_window.remove_widget(ico_wgt), .2)

            # Check if any warnings are still present. If so, update vertical
            # position
            current_floats = [i for i in self.root_window.children
                              if isinstance(i, WarningFloat) or
                              isinstance(i, RemoveFloat) or
                              isinstance(i, IconFloat)]

            if len(current_floats) > 1:

                for wgt in current_floats:
                    if pos_h > wgt.pos[1]:
                        pos = (wgt.pos[0], wgt.pos[1] + 15 + h)
                        Animation(pos=pos, t="out_quart").start(wgt)

        # The spacing variable will be updated depending on the number of
        # warning floats present in the root window.
        spacing = 0

        # Check if float is already present.
        try:
            wf = [x for x in self.root_window.children
                  if isinstance(x, WarningFloat)]

            # If a warning with the same text is already present, ignore
            if [x for x in wf if x.text == text]:
                return

            # Else, update the spacing according to the number of floats
            # present
            spacing += sum([x.height + 15 for x in wf])

        except IndexError:
            pass

        # Get position from root window
        x, y = self.root_window.size

        # Create widget
        check_wgt = WarningFloat(text=text, opacity=0, markup=True)
        check_wgt.texture_update()
        check_wgt.pos = x - check_wgt.size[0] - 20, \
            y - check_wgt.size[1] - 60 - spacing
        # Create remove button
        rm_wgt = RemoveFloat(pos=(x - 38, y - 75 - spacing), opacity=0)
        rm_wgt.bind(on_release=lambda arg: fade_out())
        # Create icon float
        ico_wgt = IconFloat(pos=(check_wgt.pos[0] - 20,
                                 check_wgt.pos[1] + (check_wgt.height / 2)
                                 - 10), opacity=0)

        # Determine background color
        if t in ["error", "warning"]:
            check_wgt.background_src = "data/backgrounds/red_noise.png"
            ico_wgt.icon_src = "data/backgrounds/warning_icon.png"
        else:
            check_wgt.background_src = "data/backgrounds/green_noise.png"
            ico_wgt.icon_src = "data/backgrounds/check_icon.png"
            ico_wgt.size = (20, 16)
            ico_wgt.pos[1] += 2

        # Add widget
        self.root_window.add_widget(check_wgt)
        self.root_window.add_widget(rm_wgt)
        self.root_window.add_widget(ico_wgt)

        # Set animations
        fade_in()
        Clock.schedule_once(lambda arg: fade_out(), dl)

    def check_partitions_file(self):
        """
        This will make some checks on the partitions file provided by the
        user. It will check for errors in the format of the file itself,
        and whether the partitions are correctly defined
        """

        aln_obj = self.alignment_list.retrieve_alignment(self.rev_infile)

        part_obj = Partitions()
        part_obj.set_length(aln_obj.locus_length)

        er = part_obj.read_from_file(self.partitions_file)

        if not er:
            return

        elif isinstance(er, Exception):
            self.dialog_warning("Invalid partitions file", er.value)
            return er

        # Check for the validity of the partitions file
        # if isinstance(er, InvalidPartitionFile):
        #     return self.dialog_floatcheck(
        #         "The provided partitions file is invalid. Please check "
        #         "the file or replace with an appropriate one.", t="error")
        #
        # # Check for the validity of the partitions file
        # if isinstance(aln_er, InvalidPartitionFile):
        #     return self.dialog_floatcheck(str(aln_er), t="error")
        # else:
        #     return True

    def save_reverseconc_settings(self, use_parts=False):
        """
        Handles the information provided by the LoadDialog with settings
        for the reverse concatenation
        :param use_parts: boolean. If True, use a partition file. Else
        Use user defined partitions
        """

        if not self.partitions_file and not self.rev_infile and \
                not self.use_app_partitions and \
                self.main_operations["reverse_concatenation"]:
            return self.dialog_floatcheck("Please provide a partition "
                    "file and file to reverse concatenate OR use defined "
                    "partitions defined in side panel.", t="error")

        if not self.use_app_partitions:

            if not self.partitions_file:
                # Check if partition file was selected. If not, return
                # warning
                return self.dialog_floatcheck(
                    "Please provide a partitions file.", t="error")

            if not self.rev_infile:
                # Check if file to reverse concatenate was specified. If
                # not, return warning
                return self.dialog_floatcheck(
                    "Please provide a file to reverse concatenate.",
                    t="error")

            # Check for the validity of the partitions file
            er = self.check_partitions_file()

            # Prevent further action is the partitions file is invalid
            if er:
                return

        if self.main_operations["reverse_concatenation"]:
            self.screen.ids.rev_conc.background_normal = \
                "data/backgrounds/bt_process.png"
            self.screen.ids.rev_conc.text = "Active"

        else:
            self.screen.ids.rev_conc.background_normal = \
                "data/backgrounds/bt_process_off.png"
            self.screen.ids.rev_conc.text = "OFF"

        self.dismiss_popup()

    def dialog_format(self):
        """
        Creates the dialog containing the buttons to select output formats.
        """

        # Inherits the layout defined in the .kv file under <FormatDialog>
        content = FormatDialog(cancel=self.dismiss_popup)

        # This will mark the respective buttons for each format as active
        # or not depending on whether they have been previously selected
        # or not. It allows the selected button states to be persistent
        # when visiting the popup multiple times
        for idx, wgt in content.ids.items():
            if idx in self.output_formats:
                wgt.state = "down"
            else:
                wgt.state = "normal"

        # Show popup
        self.show_popup(title="Choose output format", content=content,
                        size=(300, 500))

    def dialog_reverse_inlist(self):

        content = InputList(cancel=self.dismiss_subpopup)

        def set_infile(txt, wgt):

            # Update file to reverse concatenate button
            self.rev_infile = txt
            self._popup.content.ids.rev_infile.background_normal = \
                "data/backgrounds/bt_process.png"
            self._popup.content.ids.rev_infile.text = basename(txt)
            # Deactivate use defined partitions
            self._popup.content.ids.use_parts.state = "down"
            self.use_app_partitions = False

            self.dismiss_subpopup()

        if self.file_list:
            for infile in self.file_list:
                bt = Button(text=basename(infile), size_hint_y=None,
                            height=30, bold=True, shorten=True,
                            shorten_from="right", halign="center",
                            valign="middle")
                bt.text_size[0] = bt.size[0] * 3.5
                bt.bind(on_release=partial(set_infile, infile))
                content.ids.rev_inlist.add_widget(bt)

        self.show_popup(title="Choose input file",
                        content=content,
                        size_hint=(.5, .8),
                        popup_level=2)

    def dialog_reverse_concatenation(self, title="Choose input file"):
        """
        Generates a general purpose file chooser to request additional data
        :param title: string, A custom title for the load data dialog popup
        """

        content = RevConcDialog(cancel=self.dismiss_popup)
        content.ids.rev_conc.active = \
            self.main_operations["reverse_concatenation"]

        # Check if use app partitions button has been selected
        if self.use_app_partitions:
            content.ids.use_parts.state = "down"

        # Check if partitions file was already selected. If so, update the
        # corresponding button
        if self.partitions_file:
            content.ids.part_file.background_normal = \
                "data/backgrounds/bt_process.png"
            content.ids.part_file.text = basename(self.partitions_file)
            content.ids.defined.collapse = True
            content.ids.manual.collapse = False

        # Check if input file to reverse concatenate was already selected.
        #  If so, update the corresponding button
        if self.rev_infile:
            content.ids.rev_infile.background_normal = \
                "data/backgrounds/bt_process.png"
            content.ids.rev_infile.text = basename(self.rev_infile)

        self.show_popup(title=title, content=content, size=(450, 360))

        if self.use_app_partitions:
            content.ids.defined.collapse = False
            content.ids.manual.collapse = True
        else:
            content.ids.manual.collapse = False
            content.ids.defined.collapse = True

    def dialog_load_partfile(self):

        content = SaveDialog(cancel=self.dismiss_subpopup,
                             bookmark_init=self.bookmark_init)

        # If input files have already been provided, use their directory
        # as a starting point for the partition file chooser. Otherwise,
        # use the home path
        if self.file_list:
            content.ids.sd_filechooser.path = unicode(
                dirname(self.file_list[0]))
        else:
            content.ids.sd_filechooser.path = self.home_path

        content.ids.sd_filechooser.text = "partition_file"
        content.ids.txt_box.height = 0
        content.ids.txt_box.clear_widgets()

        self.show_popup(title="Choose partition file",
                        content=content, size_hint=(.9, .9),
                        popup_level=2)

    def dialog_filechooser(self, idx=None, popup_level=1):
        """
        Generates a file chooser popup for the user to select an output file

        :param idx: string. An id of where the filechooser is calling. This
        allows the addition of custom behaviours for different dialogs
        :param popup_level: int, level of popup. 1 for _popup, 2 for
        _subpopup and 3 for exit popup
        """

        def update_bw_box(ext_wgt, bw_wgt):
            """
            Function used to update the disable status of the black and white
            check box depending on the plot extensions that is selected
            """
            if ext_wgt.text != ".png":
                bw_wgt.ids.bw.disabled = True
                bw_wgt.ids.bw.active = False
                bw_wgt.ids.lbl.color = (.7, .7, .7, .7)
            else:
                bw_wgt.ids.bw.disabled = False
                bw_wgt.ids.lbl.color = (1, 1, 1, 1)

        # Determines cancel action depending on popup_level
        cancel = self.dismiss_popup if popup_level == 1 else \
            self.dismiss_subpopup

        # Lists the idx that require the selection of file extensions
        idx_with_ext = ["export_graphic"]

        # Lists the idx that do not required file name
        idx_no_file = ["ortho_dir", "zorro_dir", "protein_db", "mcl_fix",
                       "usearch_fix", "cds_db", "ima2_popfile"]

        # Maps idx for which an extension label is provided in the filechooser
        # The key is the idx, the value is the extension to appear in the label
        idx_ext_label = {"export_outliers": ".txt",
                         "export_table": ".csv",
                         "export": ".txt"}

        # Inherits the layout defined in the .kv file under <SaveDialog>
        content = SaveDialog(cancel=cancel,
                             bookmark_init=self.bookmark_init)

        # Maps the filechooser idx to the corresponding title
        title_map = {
            "main_output":
                "Choose destination directory of output file(s)",
            "export":
                "Choose file for exportation",
            "export_table":
                "Export as table...",
            "export_graphic":
                "Export as graphic...",
            "zorro_dir":
                "Choose directory with ZORRO weight files",
            "protein_db":
                "Choose protein sequence database file",
            "cds_db":
                "Choose cds/transcript sequence database file(s)",
            "export_outliers":
                "Export outliers...",
            "ortho_dir":
                "Choose destination directory for OrthoMCL output files...",
            "mcl_fix":
                "Select the MCL executable...",
            "usearch_fix":
                "Select the USEARCH executable...",
            "ima2_popfile":
                " Select the population file for IMa2 format"
        }

        # Add extension selection spinner, if idx in idx_with_ext
        if idx in idx_with_ext:
            ext_spinner = ExtSpinner()
            content.ids.txt_box.ids["ext"] = ext_spinner
            bw_box = BWCheck()
            ext_spinner.bind(text=lambda x, y: update_bw_box(ext_spinner,
                                                             bw_box))
            content.ids.txt_box.add_widget(ext_spinner)
            content.ids.txt_box.add_widget(bw_box)
            content.ids.txt_box.ids["bw"] = bw_box

        # Remove file name text input if idx in idx_no_file
        if idx in idx_no_file:
            content.ids.txt_box.clear_widgets()
            content.ids.txt_box.height = 0

        # Add extension label if idx in idx_ext_label
        if idx in idx_ext_label:
            ext_label = ExtLabel()
            ext_label.text = idx_ext_label[idx]
            content.ids.txt_box.add_widget(ext_label)

        # Custom behaviour for main output file chooser dialog
        if idx == "main_output":
            # If the main operation is conversion or reverse concatenation,
            # remove the TextInput asking for the file name
            if self.main_operations["conversion"] or \
                    self.main_operations["reverse_concatenation"]:
                content.ids.txt_label.text = "Suffix:"
                if self.conversion_suffix:
                    content.ids.text_input.text = self.conversion_suffix
                else:
                    content.ids.text_input.hint_text = \
                        "Suffix that will be appended to each output file " \
                        "name. You may leave empty."
            else:
                title_map["main_output"] = "Choose output file"

        # Custom behaviour for orthology output directory
        elif idx == "ortho_dir":
            if self.ortho_dir:
                content.ids.sd_filechooser.path = unicode(self.ortho_dir)

        if idx == "cds_db":
            content.ids.sd_filechooser.multiselect = True

        # Save output file for conversion/concatenation purposes
        # Providing this operation will allow the filechooser widget to
        # know which output file is this
        content.ids.sd_filechooser.text = idx

        self.show_popup(title=title_map[idx], content=content,
                        popup_level=popup_level)

    def dialog_taxafilter(self):
        """
        Generates dialog for taxa filter in the additional options of the
        process screen
        """

        content = TaxaFilterDialog(cancel=self.dismiss_all_popups)

        self.show_popup(title="Advanced taxa filter", content=content,
                        size=(350, 400))

    def dialog_codonfilter(self):
        """
        Generates dialog for alignment filter in the additional options of
        the process screen
        """

        content = CodonFilterDialog(cancel=self.dismiss_popup)

        self.show_popup(title="Advanced alignment filter", content=content,
                        size=(330, 280))

    def dialog_filter(self):
        """
        Generates the settings popup for filtering options
        """

        content = FilterDialog(cancel=self.dismiss_popup)
        # Update filter values if they were changed
        if self.missing_filter_settings:
            # For within alignment
            content.ids.within_chk.active = \
                self.missing_filter_settings[0][0]
            content.ids.gap_sli.value = self.missing_filter_settings[0][1]
            content.ids.mis_sli.value = self.missing_filter_settings[0][2]

            # For minimum taxa
            content.ids.among_chk.active = \
                self.missing_filter_settings[1][0]
            content.ids.min_taxa.value = self.missing_filter_settings[1][1]

        content.ids.gap_filter.active = self.secondary_options["gap_filter"]

        self.show_popup(title="Set filter thresholds", content=content,
                        size=(400, 490))

    def dialog_variationfilter(self):
        """
        Generates the settings popup for variation filter options
        """

        p_list = [("var_min", "var_min_val"), ("var_max", "var_max_val"),
                  ("inf_min", "inf_min_val"), ("inf_max", "inf_max_val")]

        content = VariationFilterDialog(cancel=self.dismiss_popup)

        content.ids.variation_filter.active = self.secondary_options[
            "variation_filter"]

        # Populate content with previous filters
        for i, p in enumerate(self.variation_filter):
            content.ids[p_list[i][0]].active = True if p else False
            content.ids[p_list[i][1]].text = str(p) if p else ""

        self.show_popup(title="Set variation filter options",
                        content=content, size=(250, 440))

    def dialog_filter_report(self, filtered_stats, nalns):
        """
        Generates graphic report at the end of a processing operation
        with active filters
        :param filtered_stats: Dictionary with the values of the filtered
         alignments for each filter type
        :param nalns: integer, total number of processed alignments
        """

        content = FilterGraphicReport(cancel=self.dismiss_popup)

        # Set total number of alignments
        content.ids.total_alns.text = str(nalns)

        # Set final alignments
        final_alns = sum([x for x in filtered_stats.values() if x])
        content.ids.final_box.size_hint_x = float(final_alns) / float(nalns)
        content.bar_color = (.44, .78, .22, 1)

        if final_alns:
            content.ids.final_txt.text = str(final_alns)
        else:
            content.ids.final_box.size_hint_x = 0.00001
            content.ids.final_txt2.text = "0"

        # Set graphics for each filter
        for f, val in filtered_stats.items():

            filt_wgt = IndividualFilterWgt()
            filt_wgt.ids.main_lbl.text = f
            if val:
                filt_wgt.ids.gf_box.size_hint_x = float(val) / float(nalns)
                filt_wgt.ids.gf_txt.text = str(val)
            elif val == 0:
                filt_wgt.ids.gf_box.size_hint_x = 0.00001
                filt_wgt.ids.gf_txt2.text = "0"
            else:
                filt_wgt.ids.gf_box.size_hint_x = 0.00001
                filt_wgt.ids.gf_txt2.text = "Not applied"

            content.ids.filter_bx.add_widget(filt_wgt)

        self.show_popup("Filter report", content=content, size=(400, 520))

    @staticmethod
    def check_variation_filters(value):
        """
        Same function as check_filters method, but checks only for
        integer compliance
        :param value: text_input.text
        """

        x = value.replace(",", ".")

        try:
            x = int(x)
        except ValueError:
            # Check if value can be converted to float by removing
            # all non digits. To avoid problems like converting 2.23
            # to 223, the first step is to get the first value before
            #  a ".", if any
            try:
                all_d = string.maketrans("", "")
                nodigs = all_d.translate(all_d, string.digits)
                x = x.encode("ascii", "ignore")
                x = int(x.translate(all_d, nodigs))
            except ValueError:
                return False

        return True, abs(x)

    @staticmethod
    def check_filters(value):
        """
        Method that validates the input of the text input in filter
        settings. It handles common mistakes, such as using "," instead
        of "." for decimal places and truncates values between the range
        of 0 and 100. If the text input cannot be converted to float,
        it will return false and the slider value will not change
        :param value: text_input.text
        """

        try:
            # Check if value can be converted to float just by replacing
            # "," to "."
            x = float(value.replace(",", "."))
        except ValueError:
            try:
                # Check if value can be converted to float by removing
                # all non digits. To avoid problems like converting 2.23
                # to 223, the first step is to get the first value before
                #  a ".", if any
                x = value.split(".")[0]
                all_d = string.maketrans("", "")
                nodigs = all_d.translate(all_d, string.digits)
                x = x.encode("ascii", "ignore")
                x = float(x.translate(all_d, nodigs))
            except ValueError:
                return False

        if x > 100:
            corrected_val = 100
        elif x < 0:
            corrected_val = 0
        else:
            corrected_val = int(x)
        return True, corrected_val

    @staticmethod
    def check_partition_split(value, prange):
        """
        Method that validates the text input for manual split of _partitions.
        Ensures that the input can be converted to int, and whether it is
        inside the available range.
        """

        try:
            # Check if value can be converted to int
            x = int(value)
        except ValueError:
            try:
                # Check if value can be converted to int by removing
                # all non digits.
                x = value
                all_d = string.maketrans("", "")
                nodigs = all_d.translate(all_d, string.digits)
                x = x.encode("ascii", "ignore")
                x = int(x.translate(all_d, nodigs))
            except ValueError:
                return False

        if x > prange[1]:
            corrected_val = prange[1]
        elif x < prange[0]:
            corrected_val = prange[0]
        else:
            corrected_val = x

        return True, corrected_val

    def dialog_execution(self):
        """
        Generates the dialog for Process execution. It also preforms several
        sanity checks before issuing the dialog
        """

        if not self.active_file_list:
            return self.dialog_floatcheck(
                "No input files were loaded or are active",
                t="error")

        content = ExecutionDialog(cancel=self.dismiss_popup)

        # Get number of files to be processed
        fs_name = self.process_grid_wgt.ids.active_file_set.text

        # Get total number of files loaded
        if fs_name == "All files":
            file_num = len(self.file_list)

        # Get active files
        elif fs_name == "Active files":
            file_num = len(self.active_file_list)

        # Get files from custom file set
        else:
            file_num = len(self.file_groups[fs_name])

        # Perform pre-execution checks

        # Get main operation
        try:
            main_op = [nm for nm, bl in self.main_operations.items()
                       if bl is True][0]
            content.ids.main_op.text = "[b][size=18][color=37abc8ff]Main " \
                "operation:[/color][/size][/b] %s" % main_op
        except IndexError:
            return self.dialog_floatcheck(
                "Please select a main operation", t="error")

        # Get secondary operations
        secondary_op = [nm for nm, bl in self.secondary_operations.items()
                        if bl is True]
        if secondary_op:
            content.ids.sec_op.text = "[b][size=18][color=37abc8ff]" \
                "Secondary operation(s):[/color][/size][/b] %s" %\
                ", ".join(secondary_op)
        else:
            content.ids.sec_op.text = "[b][size=18][color=37abc8ff]" \
                "Secondary operation(s):[/color][/size][/b] None"

        # Get output formats
        content.ids.out_form.text = "[b][size=18][color=37abc8ff]Output " \
            "format(s):[/color][/size][/b] %s" % ", ".join(
                self.output_formats)

        # Get output files
        # In case concatenation
        if main_op == "concatenation":
            # Check if an output directory has been selected
            if self.output_file == "":
                return self.dialog_floatcheck(
                    "No output file has been selected", t="error")

            out_file = basename(self.output_file)
            add_files = [out_file + "_" + nm for nm, bl in
                         self.secondary_operations.items() if bl and
                         self.secondary_options["%s_file" % nm]]
            content.ids.out_files.text = "[b][size=18][color=37abc8ff]" \
                "Output file(s):[/color][/size][/b] (%s) %s, %s" % \
                (len(add_files) + 1, out_file, ", ".join(add_files))

        # In case conversion
        if main_op == "conversion":
            # Check if an output file has been selected
            if self.output_dir == "":
                return self.dialog_floatcheck(
                    "No output directory has been selected",
                    t="error")
            try:
                # Check for additional files
                add_files = [nm for nm, bl in self.secondary_operations.items()
                             if bl and self.secondary_options["%s_file" % nm]]
                content.ids.out_files.text = "[b][size=18][color=37abc8ff]"\
                    "Output file(s):[/color][/size][/b] %s converted " \
                    "file(s)" % (file_num + file_num * len(add_files))
            # In case aln_obj has not being defined, probably because there
            # are no input files
            except AttributeError:
                return self.dialog_floatcheck(
                    "No input files havebeen selected", t="error")

        if main_op == "reverse_concatenation":
            if self.output_dir == "":
                return self.dialog_floatcheck(
                    "No output directory has been selected",
                    t="error")

        if main_op == "reverse_concatenation" and not self.rev_infile and \
                len(self.file_list) > 1:
            return self.dialog_floatcheck(
                "Reverse concatenation using _partitions defined in "
                "the app requires only one input alignment. Please select"
                " a single file to reverse concatenate in the Reverse "
                "concatenate settings", t="error")

        try:
            self.show_popup(
                title="Process execution summary - Processing %s file(s)" %
                      file_num,
                content=content,
                size=(600, 350))

        except AttributeError:
            return self.dialog_floatcheck(
                "No input files have been selected", t="error")

        # Check if main operation is reverse concatenation and if the active
        # taxa set smaller than the complete set. If so, issue a warning
        # that individual _partitions that do not contain any of the selected
        # taxa will not be written
        if main_op == "reverse_concatenation" and self.active_taxa_list != \
                self.alignment_list.taxa_names:
            self.dialog_warning("Reverse concatenation warning",
                                "A data set is being reverse "
                                "concatenated with a taxa subset "
                                "specified. Individual alignments that"
                                " do not contain the selected/active"
                                " taxa will not be written")

    def dialog_text(self, title, idx, msg=None):
        """
        Generates a simple text dialog to capture text input
        :param title: string. Title of the popup
        :param idx: string. Identifier of the dialog
        :param msg: string. Option message to appear in the dialog input
        text
        """

        if idx == "new_folder":
            content = TextDialog(cancel=self.dismiss_subpopup)
        else:
            content = TextDialog(cancel=self.dismiss_popup)

        content.text = idx

        if idx == "haplotype":
            content.ids.txt_dlg.text = self.hap_prefix

        elif idx == "usearch_db":
            content.ids.txt_dlg.text = self.usearch_db

        elif idx == "usearch_out":
            content.ids.txt_dlg.text = self.usearch_output

        elif idx == "evalue":
            content.ids.txt_dlg.text = self.usearch_evalue

        elif idx == "orto_group":
            content.ids.txt_dlg.text = self.ortholog_prefix

        elif idx == "groups":
            content.ids.txt_dlg.text = self.group_prefix

        elif idx == "change_taxon":
            content.ids.txt_dlg.text = msg
            content.old_name = msg

        elif idx == "project":

            # This will generate a default name for the project, based on
            # untitled and a number. To avoid duplications, this increases
            # the counter c while there are untitled_[0-9] projects
            c = 1
            with open(self.projects_file, "rb") as project_fh:
                project_dic = pickle.load(project_fh)
                while "untitled_{}".format(c) in project_dic:
                    c += 1

            content.ids.txt_dlg.text = "untitled_{}".format(c)
            content.ids.txt_dlg.select_all()

        if idx == "new_folder":
            popup_level = 20
        else:
            popup_level = 1

        self.show_popup(title=title, content=content, popup_level=popup_level,
                        size=(400, 150))

    def dialog_warning(self, msg1, msg2):

        content = WarningDialog(cancel=self.dismiss_subpopup)
        content.ids.warning_l.text = \
            "[b][color=#ff5555ff][size=18]%s[/size][/color][/b]\n\n%s" %\
            (msg1, msg2)

        self.show_popup(title="[b]Error![/b]", content=content,
                        size=(550, 300),
                        separator_color=[255 / 255., 85 / 255., 85 / 255., 1.],
                        popup_level=2)

    def dialog_zorro(self):

        content = ZorroDialog(cancel=self.dismiss_popup)

        content.ids.zorro_switch.active = self.secondary_options["zorro"]
        content.ids.zorro_txt.text = self.zorro_suffix

        self.show_popup(title="ZORRO support", content=content,
                        size=(400, 260))

    def save_text(self, text, idx):
        """
        Saves a text input, whose attribute depends on the idx argument
        :param text: string. The text to be assigned to the attribute
        :param idx: string. This will determine which attribute will be used
        property
        """

        if idx == "haplotype":
            self.hap_prefix = text
            self.process_options.ids.hapbt.text = text

        elif idx == "usearch_db":
            self.usearch_db = text

        elif idx == "usearch_out":
            self.usearch_output = text

        elif idx == "evalue":
            try:
                float(text)
                self.usearch_evalue = text
            except ValueError:
                return self.dialog_floatcheck(
                    "E-value must be a number", t="error")

        elif idx == "orto_group":
            self.ortholog_prefix = text

        elif idx == "groups":
            self.group_prefix = text

    def update_main_operations(self, op):
        """
        Updates the app attribute containing the main operations of the
        Process screen, self.main_operations. Only one main operation can
        be active.
        :param op: The name of the operation to turn on (all others will be
        disabled)
        """

        # The text of the Output file/directory field changes depending on
        # whether the main operation is a concatenation or a conversion
        file_text = "[size=18][b]Output file[/b][/size]\n[size=13]Save " \
                    "output file to the selected file.[/size]"
        dir_text = "[size=18][b]Output directory[/b][/size]\n[size=13]" \
                   "Save output file(s) to the selected directory.[/size]"

        self.main_operations = {k: True if k == op else False for k in
                                self.main_operations}

        # Disables output file button and other conversion/concatenation
        # specific buttons
        if op == "conversion":
            if self.output_dir == "":
                self.process_grid_wgt.ids.conv.text = "Select..."
            else:
                self.process_grid_wgt.ids.conv.text = \
                    basename(self.output_dir)
            self.process_grid_wgt.ids.output_label.text = dir_text
            Animation(height=0, d=.32, t="in_quad").start(
                self.screen.ids.sub_conc)
        elif op == "concatenation":
            if self.output_file == "":
                self.process_grid_wgt.ids.conv.text = "Select..."
            else:
                self.process_grid_wgt.ids.conv.text = \
                    basename(self.output_file)
            self.process_grid_wgt.ids.output_label.text = file_text
            Animation(height=50, d=.32, t="in_quad").start(
                self.screen.ids.sub_conc)
        elif op == "reverse_concatenation":
            if self.output_dir == "":
                self.process_grid_wgt.ids.conv.text = "Select..."
            else:
                self.process_grid_wgt.ids.conv.text = \
                    basename(self.output_dir)
            self.process_grid_wgt.ids.output_label.text = dir_text

    def save_main_operation(self, op):
        """
        This controls the appearance of the general options after the user
        chooses the main operation (Conversion vs. Concatenation). When the
        main operation is chosen for the first time, the general options
        are introduced, but further selection of the main operation only
        changes the state of the operations button.
        :param op: string, the main operations. values are "concatenation"
         and "conversion"
        """

        # If the general options widget is not shown yet, show them
        if self.process_grid_wgt.opacity == 0:
            # Add to Process screen, via the appropriate scrollview
            self.screen.ids.process_sv.add_widget(self.process_grid_wgt)

            # Store the original height of the gridlayout containing the
            # general options so that it can be updated
            self.process_height = self.process_grid_wgt.height

            # Animate the appearance of the general options via changes in
            # opacity to give a fade in effect
            Animation(opacity=1, d=.32, t="in_quad").start(
                self.process_grid_wgt)

        # Store the active main operation
        self.update_main_operations(op)

    def toggle_process_options(self):
        """
        Controls the toggling of the GridLayout with the additional options
        for the process screen.
        """

        # Shows additional options
        if self.process_grid_wgt.ids.opt_bt.state == "down":

            # Adds widget to the Process widget tree an animates its
            # appearance via changes in opacity to give impression of fade
            # in
            self.process_grid_wgt.add_widget(self.process_options)

            Animation(
                opacity=1, d=.5, t="in_quad").start(self.process_options)

            # Update the height of the GridLayout to allow scrolling
            self.process_grid_wgt.height = self.process_height + \
                sum([x.height for x in
                     self.process_options.ids.filter_grid.children]) + 100

            # Change text in the toggle button
            self.process_grid_wgt.ids.opt_bt.text = \
                "Hide additional options"

        # Hides additional options
        elif self.process_grid_wgt.ids.opt_bt.state == "normal":

            # Removes widget from widget tree. However, all settings that
            #  were change while the widget was visible will be stored
            self.process_grid_wgt.remove_widget(self.process_options)

            # Change text in the toggle button
            self.process_grid_wgt.ids.opt_bt.text = \
                "Show additional options"

    # ######################### STATISTICS SCREEN ##########################

    def toggle_stats_panel(self, force_close=None, force_open=None):
        """
        Controls the animation of the statistics panel
        :param force_close: Boolean. If True, close the stats panel
        regardless of current position
        :param force_open: Boolean. If True, open the stats panel
        regardless of current position
        """

        expanded_width = 410

        if self.screen.ids.stats_panel.width == expanded_width and \
                not force_open:
            self.sidepanel_animation(width=0,
                wgt=self.screen.ids.stats_panel)
            self.sidepanel_animation(width=40,
                wgt=self.screen.ids.sidepanel_container)

        elif self.screen.ids.stats_panel.width == 0 and not force_close:
            self.sidepanel_animation(width=expanded_width,
                wgt=self.screen.ids.stats_panel)
            self.sidepanel_animation(width=expanded_width + 40,
                wgt=self.screen.ids.sidepanel_container)

    def toggle_stats_options(self, bt, idx):
        """
        Toggles the main data exploration analyses options in the Statistics
        screen
        """

        def transfer_wgts(source_wgts, sink_gl):

            for wgt in source_wgts:
                sink_gl.add_widget(wgt)

            sink_gl.active_grid = True
            bt.disabled = False

        bt.disabled = True

        # Storage of Options buttons separated by major analysis types
        wgts = {
            "General information":
            [self.screen.ids.general_information, [SizeDistribution(),
                                                   NucAAProportions(),
                                                   TaxaDistribution()]],
            "Missing Data":
            [self.screen.ids.missing_data_opts, [GeneOccupancy(),
                                                 MissingOrto(),
                                                 MissingData(),
                                                 CumulativeMissingOrto()]],
            "Polymorphism and Variation":
            [self.screen.ids.polymorphism_data_opts, [SequenceSimilarity(),
                                                      SegregatingSites(),
                                                      LVCorrelation(),
                                                      AFS()]],
            "Outlier detection":
            [self.screen.ids.outlier_opts, [OutlierMissing(),
                                            OutlierSegregating(),
                                            OutlierSize()]]}

        # Get active type
        main_gl = self.screen.ids.main_stats_opts
        for gl in [x for x in main_gl.children if isinstance(x, OptsGrid)]:
            if gl.active_grid:
                active_gl = gl

                if active_gl.grid_name == idx:
                    bt.disabled = False
                    return
                else:
                    active_gl.clear_widgets()
                    Animation(height=5, d=.32, t="in_quad").start(active_gl)
                    active_gl.active_grid = False

        # Update sink grid height
        gl_height = sum(x.height + 10 for x in wgts[idx][1]) + 5
        Animation(height=gl_height, d=.32, t="in_quad").start(wgts[idx][0])

        Clock.schedule_once(
            lambda x: transfer_wgts(wgts[idx][1], wgts[idx][0]), .32)

    def stats_write_plot(self, plot_data, footer, plt_idx):
        """
        Provided with the data structure and a plt_idx string identifier,
        this function will create the plot file and app variable,
        and load it into the Statistics screen.

        To configure the Gene/Species/Average plot trigger, this method
        impors the stats_compliant object from data/resources/stats.py

        :param plot_data: list/np array, data structure to be used in plot
        construction
        :param footer: List, containing information on the number of
        active genes and taxa
        :param plt_idx: string, identification string of the plot. Usually
        is the text property of the issuing button.
        """

        if isinstance(plot_data, EmptyAlignment):
            return self.dialog_floatcheck(
                "Active alignment is empty", t="error")

        if plot_data:

            if "exception" in plot_data:
                if plot_data["exception"] == "empty_data":
                    self.screen.ids.plot_content.clear_widgets()
                    self.dialog_floatcheck(
                        "No data available for plotting", t="error")
                elif plot_data["exception"] == "single_alignment":
                    self.screen.ids.plot_content.clear_widgets()
                    self.dialog_floatcheck(
                        "Selected plot cannot be executed on single "
                        "alignments", t="error")

                self.screen.ids.footer_box.clear_widgets()
                self.populate_stats_footer(["NA", "NA"])
                return self.dismiss_plot_wgt()

            # Set new plot attributes
            self.current_plot, self.current_lgd, self.current_table = \
                self.stats_plt_method[plt_idx][0](**plot_data)

            # Save plot elements in a backup. This backup can then be
            # accessed when fast switching plots.
            self.plot_backups[plt_idx] = self.current_table

            pickle.dump(self.current_plot,
                        open(join(self.temp_dir, plt_idx.replace("/", "")),
                             "wb"))

            if self.current_lgd:
                self.current_plot.savefig(
                    join(self.temp_dir, self.stats_plt_method[plt_idx][1]),
                    bbox_extra_artists=(self.current_lgd,),
                    bbox_inches="tight", dpi=200)
            else:
                self.current_plot.savefig(
                    join(self.temp_dir, self.stats_plt_method[plt_idx][1]),
                    bbox_inches="tight", dpi=200)

        # Dismiss stats toggle widget, if present
        self.dismiss_stats_toggle()
        if plt_idx in stats_compliant:
            # Show toggle widget
            self.show_stats_toggle(**stats_compliant[plt_idx])
        else:
            self.previous_stats_toggle = None

        # List of plots for which an horizontal separator is available:
        hseparator_plots = ["Pairwise sequence similarity gn"]

        # List of plots with a switcher between absolute values and
        # proportions
        prop_plots = [["Segregating sites",
                       "Allele Frequency Spectrum"],
                      ["Segregating sites prop",
                       "Allele Frequency Spectrum prop"]]

        # List of plots with outlier footer
        outlier_plots = ["Missing data outliers",
                         "Missing data outliers sp",
                         "Segregating sites outliers",
                         "Segregating sites outliers sp",
                         "Sequence size outliers",
                         "Sequence size outliers sp"]

        # Adds or removes the horizontal threshold option slider from the
        # Screen footer
        if plt_idx in hseparator_plots:
            self.screen.ids.footer_box.clear_widgets()
            hwgt = HseparatorFooter()
            ylims = self.current_plot.gca().get_ylim()
            hwgt.ids.slider.min, hwgt.ids.slider.max = \
                [int(x) for x in ylims]
            hwgt.plt_file = self.stats_plt_method[plt_idx][1]
            self.screen.ids.footer_box.add_widget(hwgt)

        # Adds or removes the Absolute/Proportion switcher from the footer
        elif plt_idx in prop_plots[0] + prop_plots[1]:

            swt = PropSwitcher()

            try:
                prop_idx = prop_plots[1][prop_plots[0].index(plt_idx)]
                abs_idx = plt_idx
                swt.ids.abs.state = "down"
            except ValueError:
                abs_idx = prop_plots[0][prop_plots[1].index(plt_idx)]
                prop_idx = plt_idx
                swt.ids.proportion.state = "down"

            # Clear footer box
            self.screen.ids.footer_box.clear_widgets()

            # Create PropSwitcher instance and provide plt_idx to checkboxes
            swt.ids.abs.plt_idx = abs_idx
            swt.ids.proportion.plt_idx = prop_idx

            self.screen.ids.footer_box.add_widget(swt)

        # Adds or removes outlier footer
        elif plt_idx in outlier_plots:

            self.screen.ids.footer_box.clear_widgets()
            outlier_wgt = OutlierFooter()

            # Outliers are present when self.current_table is not None
            if self.current_table:
                outlier_num = len(self.current_table)
                # Set label with number of outliers
                outlier_wgt.ids.outlier_number.text = \
                    "{} outliers found".format(outlier_num)
                # Set the data set type based on the ending of the
                # plt_idx. Species plots always end with " sp"
                if plt_idx.endswith(" sp"):
                    outlier_wgt.ds_type = "Taxa"
                else:
                    outlier_wgt.ds_type = "Files"
            else:
                # If no outliers were found, disabled all buttons
                outlier_wgt.ids.o_rm.disabled = True
                outlier_wgt.ids.o_exp.disabled = True
                outlier_wgt.ids.o_view.disabled = True
                outlier_wgt.ids.outlier_number.text = "No outliers found"

            self.screen.ids.footer_box.add_widget(outlier_wgt)

        else:
            self.screen.ids.footer_box.clear_widgets()

        # Add plot toolbar
        self.show_plot_toolbar(toolbar_type="stats")

        self.load_plot(
                join(self.temp_dir, self.stats_plt_method[plt_idx][1]),
                self.screen.ids.plot_content)

        if footer:
            self.populate_stats_footer(footer)

    def stats_sethline(self, val, plt_file, inverted=False):
        """
        Sets an horizontal threshold bar to the current plot object
        :param val: integer, The y-axis value for the horizontal line
        :param plt_file: string, The path to the plot file
        :param inverted: Boolean, Determine whether the shade will be
        above or below the plot
        """

        # Get plot limits for patch position and size
        xlim = self.current_plot.gca().get_xlim()
        ylim = self.current_plot.gca().get_ylim()

        # Removes previous patch, if present
        try:
            self.plt_patch.remove()
        except AttributeError:
            pass

        # Get position and size
        if inverted:
            xy = (xlim[0], val)
            height = ylim[1] - val
        else:
            xy = (xlim[0], ylim[0])
            height = val - ylim[0]

        self.plt_patch = patches.Rectangle(
            xy, xlim[1], height, edgecolor="grey", fc="#dcdcdc",
            alpha=.7, zorder=10)

        self.current_plot.gca().add_patch(self.plt_patch)

        self.current_plot.savefig(join(self.temp_dir, plt_file),
                                  bbox_inches="tight", dpi=200)

        self.load_plot(join(self.temp_dir, plt_file),
                       self.screen.ids.plot_content)

    def stats_select_plot(self, gn_idx, sp_idx, avg_idx):
        """
        Dialog used for plot type selection when choosing an option from the
        stats panel
        """

        # Check whether there is data loaded
        if not self.file_list:
            return self.dialog_floatcheck("Warning: No data has been "
                                          "loaded into the app", t="error")

        # Check if active file list is empty
        if not self.active_file_list:
            return self.dialog_floatcheck("Warning: There are no active "
                                          "files selected", t="error")

        wgts = ["gene", "sp", "avg"]

        content = PlotTriageDialog(cancel=self.dismiss_popup)

        for w, f in zip(wgts, [gn_idx, sp_idx, avg_idx]):
            if not f:
                content.ids[w].disabled = True
            else:
                content.ids[w].plt_idx = f

        self.show_popup(title="Select plot type", content=content,
                        size=(400, 180))

    def stats_handle_outliers(self, operation, ds_type):
        """
        Function that handles the processing of outliers in the stats plots
        :param operation: string, three operations are available:
            .:remove: Removes outliers from data set
            .:export: Export outliers to text file
            .:view: Opens a new dialog with the list of outliers
        :param ds_type: string, two ds_types are available:
            .:Taxa: Removes taxa
            .:Files: Removes files
        """

        # Remove elements from app
        if operation == "remove":

            # The removal of taxa/files uses the export_table function to
            # dump the outliers into a file and then the
            # remove_bt_from_file function to remove them form the app
            self.export_table(self.temp_dir, "outlier_temp", warning_dlg=False)

            self.remove_bt_from_selection(ds_type,
                                     txt_file=join(self.temp_dir,
                                                   "outlier_temp.csv"))

            # Remove temporary outlier text file
            os.remove(join(self.temp_dir, "outlier_temp.csv"))

            # Refresh plot after outlier removal
            self.stats_show_plot(self.current_plt_idx[0],
                                 self.current_plt_idx[1])

            self.dialog_floatcheck(
                "{} outliers where removed".format(len(self.current_table)),
                t="info")

        if operation == "export":
            self.dialog_filechooser("export_outliers")

        if operation == "view":

            content = InputList(cancel=self.dismiss_popup)
            content.ids.ok_bx.remove_widget(content.ids.ok_bt)
            content.ids.ok_bx.width = 100
            for outlier in self.current_table:
                bt = TFButton(text=outlier[0], height=30)
                content.ids.rev_inlist.add_widget(bt)

            self.show_popup(
                title="{} {} outlier(s)".format(len(self.current_table),
                                                ds_type),
                content=content, size=(450, 500))

    def populate_stats_footer(self, footer):
        """
        Populates the footer of the Statistics screen with information on
        the active number of genes and taxa
        :param footer: list, first element contains the number of genes, the
        second element contains the number of taxa
        """

        self.screen.ids.gene_num.text = \
            "Genes: [color=37abc8ff]{}[/color]".format(footer[0])

        self.screen.ids.taxa_num.text = \
            "Taxa: [color=37abc8ff]{}[/color]". format(footer[1])

    def get_active_sets(self, file_set_name=None, taxa_set_name=None):
        """
        Returns a tuple with the file set list as first element and taxa
        set list as second element. List sets are only return for non
        None arguments. If, for example, only file_set_name is provided,
        then the return tuple will be ([file list], None)
        :param file_set_name: string, name of file set
        :param taxa_set_name: string, name of file set
        """

        if file_set_name:
            if file_set_name == "All files":
                file_set = [x for x in self.file_list]
            elif file_set_name == "Active files":
                file_set = [x for x in self.active_file_list]
            else:
                file_set = [self.filename_map[x] for x in
                            self.file_groups[file_set_name]]
        else:
            file_set = None

        if taxa_set_name:
            if taxa_set_name == "Active taxa":
                taxa_set = [x for x in self.active_taxa_list]
            elif taxa_set_name == "All taxa":
                taxa_set = [x for x in self.alignment_list.taxa_names]
            else:
                taxa_set = self.taxa_groups[taxa_set_name]
        else:
            taxa_set = None

        return file_set, taxa_set

    def stats_show_plot(self, plt_idx, additional_args=None):
        """
        Wrapper that executes plot data gathering and execution. The method
        that gathers the data for plot production runs in the background.
        Once it's finished, that data is fed to the stats_write_plot
        method that will create the plot file and load it into the program.

        :param plt_idx: string, identification string of the plot. Usually
        is the text property of the issuing button.
        :param additional_args:
        """

        # Interrupt summary statistics, if running
        self.terminate_stats = True

        # Reset scatter properties for plots
        self.screen.ids.plot_content.do_scale = True
        self.screen.ids.plot_content.do_translation = True
        self.screen.ids.plot_content.do_rotation = True

        # Check whether there is data loaded
        if not self.file_list:
            return self.dialog_floatcheck("Warning: No data has been "
                                          "loaded into the app", t="error")

        # Check if active file list is empty
        if not self.active_file_list:
            return self.dialog_floatcheck("Warning: There are no active "
                                          "files selected", t="error")

        # Set active file and taxa sets
        file_set_name = self.screen.ids.active_file_set.text
        taxa_set_name = self.screen.ids.active_taxa_set.text

        # Note: When setting file_set, the objects file_list,
        # active_taxa_list, etc should not be directly assigned since
        # changes to these objects will also change the taxa_set,
        # and consequently, the self.previous_sets
        file_set, taxa_set = self.get_active_sets(file_set_name,
                                                  taxa_set_name)

        # List of gene specific plots. These are always removed
        gene_specific = {"Pairwise sequence similarity gn":
            "similarity_distribution_gn.png",
            "Segregating sites gn": "segregating_sites_gn.png",
            "Allele Frequency Spectrum gn":
            "allele_frequency_spectrum_gn.png"}

        # Remove gene specific plots if they exist
        if plt_idx in gene_specific:
            try:
                os.remove(join(self.temp_dir, gene_specific[plt_idx]))
            except EnvironmentError:
                pass

        # This will check if the current data sets are different from the
        # previous. If not so, it will then check if there is a temporary
        # plot file for the current plt_idx. If so, do not run
        # get_stats_data and show the previous plot instead.
        if file_set == self.previous_sets["Files"] and \
                taxa_set == self.previous_sets["Taxa"]:

            # Checks if the figure file for the selected plot already exists
            if os.path.exists(join(self.temp_dir,
                                   self.stats_plt_method[plt_idx][1])):

                self.toggle_stats_panel(force_close=True)

                # Get the current plot from the backup
                self.current_table = self.plot_backups[plt_idx]
                self.current_plot = pickle.load(
                    open(join(self.temp_dir, plt_idx.replace("/", "")), "rb"))

                return self.stats_write_plot(None, None, plt_idx)
        else:
            # Update previous sets
            self.previous_sets["Files"] = file_set
            self.previous_sets["Taxa"] = taxa_set

            # Remove previous plots when data sets are modified
            for p in self.stats_plt_method.values():
                try:
                    os.remove(join(self.temp_dir, p[1]))
                except OSError:
                    pass

        self.run_in_background(
            func=get_stats_data,
            second_func=self.stats_write_plot,
            args1=[self.alignment_list, plt_idx, file_set,
                   list(taxa_set), additional_args,
                   "use_ns"],
            args2=[plt_idx])

        self.toggle_stats_panel(force_close=True)

        # Sets the current_plt_idx so that plots can be updated and
        # reload the current plot
        self.current_plt_idx = [plt_idx, additional_args]

    def dialog_select_gene(self, plt_idx):
        """
        Generates dialog for selecting gene for single gene plot creation
        """

        no_window = ["Allele Frequency Spectrum gn"]

        content = SelectGeneDialog(cancel=self.dismiss_popup)

        # Remove sliding window parameter if plt_idx does not use it
        if plt_idx in no_window:
            content.remove_widget(content.ids.sliding_window)

        content.plt_idx = plt_idx

        # By default show up to 20 files at first
        for i in range(20):
            try:
                f = basename(self.file_list[i])

                bt = TGToggleButton(text=f, id=f, state="normal", height=30,
                    shorten=True, shorten_from="right",
                    disabled_color=(1, 1, 1, 1),
                    background_disabled_down=join("data", "backgrounds",
                                                  "bt_process.png"))

                bt.text_size[0] = bt.size[0] * 3
                bt.bind(on_release=self.toggle_groups)
                content.ids.gn_grid.add_widget(bt)
                content.gene_counter += 1

            except IndexError:
                break

        try:
            if self.file_list[content.gene_counter + 1]:
                content.ids.gn_grid.add_widget(StatsMoreBt())
        except IndexError:
            pass

        self.show_popup(title="Select gene for sliding window analysis...",
                        content=content, size_hint=(.4, .9))

    def stats_search_genes(self, txt):
        """
        Searches loaded genes for the single gene display of the Statistics
        screen
        """

        # When empty search, clear grid layout
        if txt == "":
            return self._popup.content.ids.gn_grid.clear_widgets()

        self._popup.content.ids.gn_grid.clear_widgets()

        found_bts = [basename(x) for x in self.file_list if
                     txt.lower() in basename(x).lower()]

        for f in found_bts:

            bt = TGToggleButton(text=f, id=f, state="normal", height=30,
                shorten=True, shorten_from="right",
                disabled_color=(1, 1, 1, 1),
                background_disabled_down=join("data", "backgrounds",
                                              "bt_process.png"))

            bt.text_size[0] = bt.size[0] * 3
            bt.bind(on_release=self.toggle_groups)
            self._popup.content.ids.gn_grid.add_widget(bt)

    def stats_load_more_genes(self):
        """
        Functionality to the load next 20 genes in the stats gene selection
        dialog
        """

        # Remove previous StatsMorebt
        self._popup.content.ids.gn_grid.remove_widget(
            self._popup.content.ids.gn_grid.children[0])

        for i in range(self._popup.content.gene_counter,
                       self._popup.content.gene_counter + 20):

            try:
                f = basename(self.file_list[i])

                bt = TGToggleButton(text=f, id=f, state="normal", height=30,
                    shorten=True, shorten_from="right",
                    disabled_color=(1, 1, 1, 1),
                    background_disabled_down=join("data", "backgrounds",
                                                  "bt_process.png"))

                bt.text_size[0] = bt.size[0] * 3
                bt.bind(on_release=self.toggle_groups)
                self._popup.content.ids.gn_grid.add_widget(bt)

            except IndexError:
                break

        try:
            if self.file_list[self._popup.content.gene_counter + 1]:
                self._popup.content.ids.gn_grid.add_widget(StatsMoreBt())
        except IndexError:
            pass

    def stats_clear_search(self):
        """
        Clears stats gene search to default buttons
        """

        self._popup.content.ids.gn_grid.clear_widgets()

        for i in range(20):
            try:
                f = basename(self.file_list[i])

                bt = TGToggleButton(text=f, id=f, state="normal", height=30,
                    shorten=True, shorten_from="right",
                    disabled_color=(1, 1, 1, 1),
                    background_disabled_down=join("data", "backgrounds",
                                                  "bt_process.png"))
                bt.text_size[0] = bt.size[0] * 3
                bt.bind(on_release=self.toggle_groups)
                self._popup.content.ids.gn_grid.add_widget(bt)

            except IndexError:
                break

    # ##################################
    #
    # CORE RELATED METHODS AND FUNCTIONS
    #
    # ##################################

    def load_files_subproc(self, files):

        multiprocessing.freeze_support()

        def check_proc(p, dt):

            try:
                content.ids.pb.value = shared_ns.progress
                content.ids.msg.text = shared_ns.m
            except AttributeError:
                pass

            if not p.is_alive():

                Clock.unschedule(func)
                self.dismiss_popup()

                try:
                    if shared_ns.exception == "multiple_type":
                        return self.dialog_warning(
                            "Multiple sequence types detected",
                            "The selected input alignments contain more "
                            "than one sequence type (DNA, RNA, Protein). "
                            "Please select input files of the  same "
                            "sequence type")
                    if shared_ns.exception:
                        return self.dialog_floatcheck(
                            "Unexpected error when loading input data. "
                            "Check app logs", t="error")
                except Exception as e:
                    pass

                # Get the alignment object from the child thread and load
                # it into the app
                aln_obj = queue.get()
                self.load_files(file_list, aln_obj)

                # Shutting down manager that provides communication between
                # parent and child thread
                manager.shutdown()

                # Join child process and exit
                p.join()
                return

            # Kill switch
            if self.terminate_load_files:
                content.ids.msg.text = "Canceling..."

                # Issuing the kill order to the child thread. This signal
                # will propagate through the worker's methods and graciously
                # terminate
                shared_ns.stop = True

                # This small delay seems to fix manager shutdown issues
                time.sleep(.1)

                # Shutting down manager that provides communication between
                # parent and child thread
                manager.shutdown()

                # Close loading popup dialog
                self.dismiss_popup()

                # Unschedule the current function
                Clock.unschedule(func)

                # Clean up alignments that were already loaded into the
                # AlignmentList object
                self.alignment_list.remove_file(file_list)

                # Join child process and exit
                p.join()
                return

        # Check if proteome files have already been loaded. If so, issue
        # warning, since the two data types cannot be loaded in simultaneous
        if self.proteome_files:
            return self.dialog_warning("Cannot mix alignment and proteome "
                                       "data",
                                       "TriFusion currently has proteome data"
                                       " loaded. If you want to load alignment"
                                       " data, please remove the current"
                                       " data set.")

        # To support for opening all files in one or more directories, all
        # entries in files will be checked if they are directories. If so,
        # all files in that directory will be appended to file_list instead
        file_list = []
        for i in files:
            if os.path.isdir(i):
                file_list.extend([join(i, x) for x in
                                 os.listdir(i)
                                 if os.path.isfile(join(i, x))])
            elif os.path.isfile(i):
                file_list.append(i)
            else:
                # If the file has not been found, try joining the called_dir
                # path with the file name
                mod_path = join(called_dir, os.path.split(i)[1])
                if os.path.isfile(mod_path):
                    file_list.append(mod_path)

        file_list = sorted(file_list)

        if not file_list:
            return self.dialog_floatcheck("No valid input files were"
                                          " provided", t="error")

        # Set up manager and shared name space to share information
        # between background and main processes
        manager = multiprocessing.Manager()
        shared_ns = manager.Namespace()
        queue = Queue.Queue()

        # Initialize stop flag to allow thread to be killed by the user
        shared_ns.stop = False

        # Remove lock from background process
        self.terminate_load_files = False

        # Create process
        p = threading.Thread(
                target=load_proc,
                args=(self.alignment_list, file_list, shared_ns,
                      queue))

        p.daemon = True
        p.start()

        # Remove any possible previous popups
        self.dismiss_popup()

        # Create loading process dialog
        content = LoadProgressDialog()
        content.ids.msg.text = "Initializing"
        content.ids.pb.max = len(file_list)

        # Give functionality to cancel button
        content.ids.close_bt.bind(
            on_release=lambda j: setattr(self, "terminate_load_files",
                                         True))

        self.show_popup(title="Loading files", content=content,
                        size=(400, 300), auto_dissmiss=False)

        # Schedule function that checks the process' pulse
        func = partial(check_proc, p)
        Clock.schedule_interval(func, .01)

    def load_files(self, selection=None, aln_list=None):
        """
        Loads the selected input files into the program using the
        AlignmentList object provided by aln_list. The loading process is
        divided in phases:

        .: Load the AlignmentList object from a pickle object

        .: Check for invalid alignments (duplicates, invalid input
        formats, sequences with unequal length) and issue warning message

        .: If any valid input files passed all checks, then the app
        structures are updated

        To allow different files to be loaded in different occasions,
        all checks are performed on the aln_list object, and not in the
        app attribute self.alignment_list

        :param aln_pickle: string, path to pickle file with the
        AlignmentObject
        :param selection: list, with the path of all files provided to
        the app
        """

        # Check for consistency in sequence type across alignments
        if self.sequence_types:
            current_seq_type = set([self.sequence_types] +
                aln_list.format_list())
        else:
            current_seq_type = aln_list.format_list()

        if not self.sequence_types and current_seq_type:
            self.sequence_types = list(current_seq_type)[0]

        # Set the new alignment_list attribute
        self.alignment_list = aln_list
        # Updating active alignment list
        self.active_taxa_list = self.alignment_list.taxa_names

        # removes bad alignment files from selection list
        selection = [path for path in selection if path not in
                     aln_list.bad_alignments + aln_list.non_alignments]

        # If duplicate alignments were loaded, issue a warning
        if self.alignment_list.duplicate_alignments:
            self.dialog_floatcheck(
                "Duplicate input alignments detected and ignored",
                t="error")
            # Reset the duplicate alignment storage, so that it doesn't
            # issue the warning every time data is loaded
            self.alignment_list.duplicate_alignments = []

        # Checking if there are invalid input alignments
        if aln_list.bad_alignments or aln_list.non_alignments:
            msg = ""
            if aln_list.bad_alignments:
                msg += "The following input file(s) could not be open:"\
                    "\n\n[b]%s[/b]\n\n" % "\n".join(basename(x) for x in
                                                    aln_list.bad_alignments)

                # Reset the bad_alignments attribute to avoid triggering
                # this in subsequent loading operations
                self.alignment_list.bad_alignments = []

            if aln_list.non_alignments:
                msg += "The following input file(s) contain(s) " \
                    "sequences of unequal length:\n\n[b]%s[/b]" % \
                    "\n".join(basename(x) for x in aln_list.non_alignments)

                # Reset the non_alignments attribute to avoid triggering
                # this in subsequent loading operations
                self.alignment_list.non_alignments = []

            # After gathering all information on bad input files and
            # non-aligned sequences, issue the error message
            self.dialog_warning("Invalid input file(s) detected", msg)

        # If new alignments have been provided update app structures. If
        # all provided alignments were bad, ignore this
        if selection:
            # If data has been previously loaded, updated these attributes
            if self.file_list:
                # Updating complete and active file lists
                self.file_list.extend(selection)
                self.active_file_list.extend(selection)
                # Update the filename - path mapping attribute
                self.filename_map = dict(
                    list(self.filename_map.items()) + list((x, y) for x, y in
                    zip([basename(x) for x in selection], selection)))

            # If no data has been previously loaded, set the attributed
            else:
                # Set an attribute with the input file list
                self.file_list = selection
                # Setting active file list and path list
                self.active_file_list = deepcopy(self.file_list)
                # Sett the filename - path mapping attribute
                self.filename_map = dict((x, y) for x, y in zip(
                    [basename(x) for x in selection], selection))

            # Update active taxa list
            self.update_taxa()
            # Populates files and taxa contents
            self.update_tabs()
            # Gathers taxa  and file information
            self.get_taxa_information()

            # Issue final dialog with the number of files successfully
            # loaded
            self.dialog_floatcheck(
                "%s file(s) successfully loaded" % len(selection),
                t="info")

            # If on statistics screen, issue the summary statistics overview
            if self.current_screen == "Statistics":
                self.statistics_show_summary()

    @staticmethod
    def get_taxon_information(tx, aln_list):
        """
        Akin to the get_taxa_information method, but it only looks for
        the information of a single taxon.

        :param tx: string, taxon name
        :param aln_list: AlignmentList object from where the information
        will be retrieved

        """

        tx_inf = {}
        # Counter for alignment missing the taxa
        tx_missing = 0

        # Get full sequence
        sequence = []
        # This assures that the information will only be gathered if the
        # active data set is not empty
        if aln_list.alignments:
            for aln in aln_list:
                if tx in aln.taxa_idx:
                    sequence.append(aln.get_sequence(tx))
                else:
                    tx_missing += 1

            # Retrieve missing data symbol
            missing_symbol = aln.sequence_code[1]

            sequence = "".join(sequence)

            # Get sequence length
            seq_len = len(sequence)
            tx_inf["length"] = seq_len
            # Get indel number
            tx_inf["indel"] = sequence.count("-")
            # Get missing data
            tx_inf["missing"] = sequence.count(missing_symbol)
            # Get effective sequence length in absolute and percentage
            tx_inf["effective_len"] = seq_len - \
                (tx_inf["indel"] + tx_inf["missing"])

            tx_inf["effective_len_per"] = \
                round((tx_inf["effective_len"] * 100) / seq_len, 2)

            # Get number of files containing the taxa in absolute and
            # percentage
            tx_inf["fl_coverage"] = len(aln_list.alignments) -\
                tx_missing
            tx_inf["fl_coverage_per"] = \
                round(((tx_inf["fl_coverage"] * 100) /
                       len(aln_list.alignments)), 2)

        else:
            # This handles the case where the active data set is empty
            tx_inf["length"] = "NA"
            tx_inf["indel"] = "NA"
            tx_inf["missing"] = "NA"
            tx_inf["effective_len"] = "NA"
            tx_inf["effective_len_per"] = "NA"
            tx_inf["fl_coverage"] = "NA"
            tx_inf["fl_coverage_per"] = "NA"

        return tx_inf

    def get_taxa_information(self):
        """
        This method will gather all available information for all taxa and
        sets a number of related attributes. This methods should only be
        used to generate data for the complete data set. The active data
        set is usually only gathered for a single taxon and upon demand.

        Sets the self.original_tx_inf attribute. Contains the relevant
        information for the taxa popup. All corresponding values are the
        result of additions across all active input alignments. Contains
        the following keys:

            - length: Total alignment length including missing data and gaps
            - indel: The number columns containing indel/gaps
            - missing: The number of columns containing missing data
            - effective_len: Alignment length excluding gaps and missing
            data
            - effective_len_per: Same as above but in percentage
            - fl_coverage: The number of files containing the focal taxon
            - fl_coverage_per: Same as above but in percentage
        """

        # Reset trigger for taxa update
        self.trigger_taxa_update = False

        for tx in self.active_taxa_list:
            # Add entry to storage dictionary
            self.original_tx_inf[tx] = \
                self.get_taxon_information(tx, self.alignment_list)

    def get_file_information(self, file_name=None, mode="alignment"):
        """
        Similar to get_taxa_information, but generating information for the
        files in the file tab.
        :param mode: string. The type of file information to retrieve. May
        be 'alignment' or 'proteome'

        :return: file_inf (dictionary). Contains all relevant content for
        the file popup. It contains the following keys:

            ..:mode: alignment:
                - aln_format: The format of the input file
                - seq_type: The sequence type. If DNA, RNA, Protein.
                - n_taxa: Number of taxa
                - aln_len: Length of the alignment
                - model: The model of sequence evolution, if applicable.
                 This is
                usually only present on Nexus input format
            ..mode: proteome:
                - n_seq: Number of sequences
                - n_res: Number of residues
        """

        # main storage
        file_inf = {}

        if mode == "alignment":
            # Iterating over file path since the name of the Alignment
            # objects contain the full path
            if self.alignment_list:
                aln_path = self.filename_map[file_name]
                aln = self.alignment_list.retrieve_alignment(aln_path)

                # Get input format
                file_inf["aln_format"] = aln.input_format

                # Get sequence type
                file_inf["seq_type"] = aln.sequence_code[0]

                # Get number of species
                file_inf["n_taxa"] = len([x for x in
                    aln.taxa_idx if x in self.active_taxa_list])

                # Get length of largest sequence if not aligned, or
                # alignment length
                file_inf["aln_len"] = aln.locus_length

        elif mode == "proteome":
            for p in self.proteome_files:
                p_name = basename(p)
                file_inf[p_name] = {}
                n_seq, n_res = self.get_proteome_information(p)

                file_inf[p_name]["n_seq"] = n_seq
                file_inf[p_name]["n_res"] = n_res

        return file_inf

    @staticmethod
    def get_proteome_information(proteome_file):
        """
        Returns informative stats for a given proteome_file
        :param proteome_file: string. Path to proteome file
        """

        file_handle = open(proteome_file)

        n_seq = 0
        n_res = 0

        for line in file_handle:
            if line.startswith(">"):
                n_seq += 1
            else:
                n_res += len(line.strip())

        file_handle.close()

        return n_seq, n_res

    def update_process_switch(self, switch_id, state):
        """
        Listens and updates the attribute process_switches when their state
        changes.
        :param switch_id: string, name of the switch according to the keys
        in process_switches
        :param state: Boolean, current state of the corresponding switch
        """

        if switch_id in self.secondary_operations:
            self.secondary_operations[switch_id] = state
        else:
            self.secondary_options[switch_id] = state

    def dialog_orto_execution(self):
        """
        Creates and populates the pre-execution dialog for orthology search
        """

        # Check for input proteomes
        if not self.proteome_files:
            return self.dialog_floatcheck(
                "Please provide proteome files as input data", t="error")

        # Check for output directory
        if self.ortho_dir == "":
            return self.dialog_floatcheck(
                "Please specify an output directory for orthology results",
                t="error")

        # Check for the USEARCH executable
        if not self.usearch_file:
            return self.dialog_usearch_fix()

        # check for the MCL executable
        if not self.mcl_file:
            return self.dialog_mcl_fix()

        content = OrtoExecutionDialog(cancel=self.dismiss_popup)

        self.show_popup(
            title="Orthology search execution summary - Processing %s"
                  " file(s)" % len(self.active_proteome_files),
            content=content,
            size=(550, 350))

        # Correct minimum species filter
        if self.orto_min_sp > len(self.active_proteome_files):
            self.orto_min_sp = len(self.active_proteome_files)
            self.dialog_floatcheck(
                "Corrected the minimum number of taxa per cluster value "
                "to the maximum number allowed by the active proteomes.",
                t="warning")

        content.ids.gene_filter.text = \
            "[b][size=18][color=37abc8ff]Maximum number of gene copies " \
            "per cluster:[/color][/size][/b] %s" %\
            self.orto_max_gene

        content.ids.sp_filter.text = \
            "[b][size=18][color=37abc8ff]Minimum number of taxa per " \
            "cluster:[/color][/size][/b] %s" %\
            self.orto_min_sp

        content.ids.eval.text =\
            "[b][size=18][color=37abc8ff]USEARCH e-value threshold:" \
            "[/color][/size][/b] %s" %\
            self.usearch_evalue

        content.ids.inflation.text = \
            "[b][size=18][color=37abc8ff]MCL inflation value(s):" \
            "[/color][/size][/b] %s" %\
            ", ".join(self.mcl_inflation)

        content.ids.threads.text = \
            "[b][size=18][color=37abc8ff]Threads :[/color][/size][/b] %s" %\
            self.screen.ids.usearch_threads.text

    def orthology_search_exec(self):
        """
        Main function that executes all queued procedures of the orthology
        module
        """

        tasks_map = {"schema": "Database setup",
                     "adjust": "Adjusting Fasta files",
                     "filter": "Filtering Fasta files",
                     "usearch": "Running USearch",
                     "parse": "Parsing USearch output",
                     "pairs": "Obtaining pairs",
                     "mcl": "Running MCL",
                     "dump": "Dumping groups",
                     "filter_groups": "Filtering group files"}

        wgt_ref = {}

        def finish_op(wgt, op_name):
            wgt.ids.load_box.clear_widgets()
            check_wgt = ProgressFinished()
            wgt.ids.load_box.add_widget(check_wgt)
            wgt.ids.main_lbl.font_size = 18
            wgt.ids.main_lbl.color = (.7, .7, .7, 1)
            wgt.ids.main_lbl.text = tasks_map[op_name]
            wgt.ids.secondary_lbl.text = ""
            wgt.ids.secondary_lbl.color = (.7, .7, .7, 1)
            closed_ops.append(op_name)

        def start_op(wgt):
            wgt.ids.load_box.clear_widgets()
            spinner = LoadSpinner()
            counter = LoadCounter()
            wgt_ref["spinner"] = spinner
            wgt_ref["counter"] = counter
            wgt.ids.load_box.add_widget(counter)
            wgt.ids.load_box.add_widget(spinner)
            wgt.ids.main_lbl.font_size = 18
            wgt.ids.main_lbl.color = self._blue
            wgt.ids.secondary_lbl.color = self._blue

        def check_process(p, dt):
            """
            Checks the status of the background process "p" and updates
            the progress dialog label
            """

            # Listens for cancel signal
            if self.terminate_orto_search:
                # Issuing the kill order to the child thread. This signal
                # will propagate through the worker's methods and graciously
                # terminate
                shared_ns.stop = True

                # If a subprocess has been issued, kill it
                try:
                    if shared_ns.subp:
                        os.kill(shared_ns.subp, signal.SIGINT)
                        shared_ns.subp = None
                except AttributeError:
                    pass

                # This small delay seems to fix manager shutdown issues
                time.sleep(.1)

                # Shutting down manager that provides communication between
                # parent and child thread
                manager.shutdown()

                # Close loading popup dialog
                self.dismiss_popup()

                # Unschedule the current function
                Clock.unschedule(func)

                # Clear sqlitedb
                os.remove(join(self.temp_dir, "orthoDB.db"))

                # Join child process and exit
                p.join()
                return

            # Update progress dialog
            try:

                if shared_ns:
                    # If current_op is not an empty list, it means that this
                    # is not the first task.
                    if current_op:
                        # Check if the task defined in the previous cycle is
                        # the same as the current task. If not, do this
                        if shared_ns.task != current_op[0]:
                            # Check if all finished tasks where closed
                            for op in [x for x in shared_ns.finished_tasks
                                       if x not in closed_ops]:
                                # This will get previously finished tasks
                                # that were completed sooner
                                if op != current_op[0]:
                                    pwgt = content.ids[op]
                                    finish_op(pwgt, op)
                            # Get the previous task progress widget
                            pwgt = content.ids[current_op[0]]
                            # Finish animation for previous task
                            finish_op(pwgt, current_op[0])
                            # Set current task progress as the current_op
                            current_op[0] = shared_ns.task
                            # Get the current task progress widget
                            cwgt = content.ids[shared_ns.task]
                            # Start progress animations for current task
                            start_op(cwgt)

                        else:
                            pass
                    # This means that this is the first task
                    else:
                        # Add the first task to current_op
                        current_op.append(shared_ns.task)
                        # Get current progress task widget
                        cwgt = content.ids[shared_ns.task]
                        # Start progress animations for current task
                        start_op(cwgt)

                    # Update the current spinner... spin
                    wgt_ref["spinner"].ids.img.rotation -= 10

                    # Update the counter percentage
                    if shared_ns.total and shared_ns.counter and \
                            current_op[0] != "usearch" and \
                                    current_op[0] != "mcl":

                        # Get percentage
                        perc = int((float(shared_ns.counter) /
                                    float(shared_ns.total)) * 100.)

                        # Set percentage on progress dialog
                        wgt_ref["counter"].text = "{}%".format(perc)

                        # Set secondary progress label, if any
                        if shared_ns.msg:
                            content.ids[shared_ns.task].\
                                ids.secondary_lbl.text = shared_ns.msg

                    if current_op[0] == "usearch" or current_op[0] == "mcl":

                        # Provide uncertain counter for usearch and mcl
                        # operations.
                        wgt_ref["counter"].text = "??"

            except AttributeError:
                pass

            if not p.is_alive():

                Clock.unschedule(func)
                self.dismiss_popup()

                try:
                    if shared_ns.exception:
                        self.dialog_floatcheck(shared_ns.exception,
                                               t="error")
                except Exception as e:
                    Logger.exception(e.message)

                try:
                    # Set the protein database file
                    self.protein_db = join(self.ortho_dir,
                                           "backstage_files",
                                           self.usearch_db)
                    self.dialog_search_report(shared_ns.stats,
                                              shared_ns.groups)
                except Exception as e:
                    Logger.exception(e.message)
                    if not shared_ns.exception:
                        self.dialog_floatcheck("Unexpected error when "
                                               "search orthologs. Check "
                                               "app logs.", t="error")

                manager.shutdown()
                p.join()

        # Create directory that will store intermediate files during
        # orthology search
        int_dir = join(self.ortho_dir, "backstage_files")
        if not os.path.exists(int_dir):
            os.makedirs(int_dir)

        # Set up manager and shared name space to share information
        # between background and main processes
        manager = multiprocessing.Manager()
        shared_ns = manager.Namespace()

        shared_ns.stop = False
        shared_ns.exception = None
        shared_ns.total = 0
        shared_ns.counter = 0
        shared_ns.msg = None
        current_op = []
        closed_ops = []

        # Remove lock from background process
        self.terminate_orto_search = False

        # Create Process instance
        p = threading.Thread(
            target=orto_execution,
            args=(
                shared_ns,
                self.temp_dir,
                self.active_proteome_files,
                self.protein_min_len,
                self.protein_max_stop,
                self.usearch_file,
                self.usearch_evalue,
                self.screen.ids.usearch_threads.text,
                self.usearch_output,
                self.mcl_file,
                self.mcl_inflation,
                self.ortholog_prefix,
                self.group_prefix,
                self.orto_max_gene,
                self.orto_min_sp,
                self.sqldb,
                self.ortho_dir,
                self.usearch_db))

        p.daemon = True
        p.start()

        # Remove any possible previous popups
        self.dismiss_popup()
        # Create Progression dialog
        content = OrtoProgressDialog()

        # Populate progress dialog content
        for op, txt in tasks_map.items():
            content.ids[op].ids.main_lbl.text = txt

        self.show_popup(title="Running Orthology Search", content=content,
                        size=(340, 600), auto_dissmiss=False)

        # Schedule function that checks the process' pulse
        func = partial(check_process, p)
        Clock.schedule_interval(func, .2)

    def dialog_fileoverwrite(self, msg):
        """
        General dialog for querying whether an existing file should be
        overwritten or skipped
        """

        content = FileOverwriteDialog(cancel=self.dismiss_subpopup)

        content.ids.dlg_text.text = msg

        self.show_popup(title="File already exists...", content=content,
                        size=(450, 220), popup_level=2)

    def process_exec(self):
        """
        Main function that executes all queued procedures of the process
        module
        """

        ops_map = {"reverse_concatenation": "Reverse concatenation",
                   "collapse_filter": "Filter (Collapse)",
                   "filter": "Filter",
                   "concatenation": "Concatenation",
                   "collapse": "Collapse",
                   "gcoder": "Gap coding",
                   "consensus": "Consensus",
                   "write": "Writing output"}

        def build_progress_dialog(wgt):

            # Get the main and additional file operations
            main_ops = []
            additional_ops = []

            # Iterate over a fixed order list of possible ops
            for op in ["reverse_concatenation",
                       "collapse_filter",
                       "filter",
                       "concatenation",
                       "collapse",
                       "gcoder",
                       "consensus"]:

                try:
                    # Add to main ops if op is a main operation and set to True
                    if self.main_operations[op]:
                        main_ops.append(op)
                except KeyError:
                    pass

                try:
                    if self.secondary_operations[op]:
                        try:
                            # Add to additional ops if op is True for secondary
                            # operations and true for *_file
                            if self.secondary_options["{}_file".format(op)]:
                                additional_ops.append(op)
                            else:
                                main_ops.append(op)
                        except KeyError:
                            # Add to main ops if op is True for secondary
                            # operations but not True for *_file
                            main_ops.append(op)
                except KeyError:
                    pass

            # Set width of dialog based on the presence or not of additional
            # files
            if additional_ops:
                width = 650
            else:
                width = 350
                wgt.ids.additional_outputs.clear_widgets()
                wgt.ids.additional_outputs.size_hint_x = None
                wgt.ids.additional_outputs.width = 0
                wgt.remove_widget(wgt.ids.additional_outputs)

            # Set height of dialog based on the number of main/additional ops
            ml = len(main_ops)
            al = len(additional_ops)
            max_ops = ml if ml > al else al

            height = 200 + (max_ops * 50)

            wgt_ref = {}

            # Add main ops
            for op in main_ops:
                p = ProgressBox()
                p.ids.main_lbl.text = ops_map[op]
                wgt.ids.main_box.add_widget(p)
                wgt_ref[op] = p

            # Create write files operation
            p = ProgressBox()
            p.ids.main_lbl.text = ops_map["write"]
            wgt.ids.main_box.add_widget(p)
            wgt_ref["write"] = p

            # Add additional ops
            for op in additional_ops:
                p = ProgressBox()
                p.ids.main_lbl.text = ops_map[op]
                wgt.ids.ad_box.add_widget(p)
                wgt_ref[op] = p

            # Give functionality to cancel button
            content.ids.cancel_bt.bind(
                on_release=lambda x: setattr(self,
                                             "terminate_process_exec",
                                             True))
            self.show_popup(title="Process execution...", content=content,
                            size=(width, height), auto_dissmiss=False)

            return wgt_ref

        def finish_op(wgt, op_name):
            wgt.ids.load_box.clear_widgets()
            check_wgt = ProgressFinished()
            wgt.ids.load_box.add_widget(check_wgt)
            wgt.ids.main_lbl.font_size = 18
            wgt.ids.main_lbl.color = (.7, .7, .7, 1)
            wgt.ids.main_lbl.text = ops_map[op_name]
            wgt.ids.secondary_lbl.text = ""
            wgt.ids.secondary_lbl.color = (.7, .7, .7, 1)

        def start_op(wgt):
            wgt.ids.load_box.clear_widgets()
            spinner = LoadSpinner()
            counter = LoadCounter()
            wgt_ref["spinner"] = spinner
            wgt_ref["counter"] = counter
            wgt.ids.load_box.add_widget(counter)
            wgt.ids.load_box.add_widget(spinner)
            wgt.ids.main_lbl.font_size = 18
            wgt.ids.main_lbl.color = self._blue
            wgt.ids.secondary_lbl.color = self._blue

        def check_process(p, man, dt):

            # Update progress dialog
            try:
                # Start only when the first task has been initiated
                if shared_ns.task:
                    # If current_op is not empty, it means some task has
                    # already been previously provided.
                    if current_op:
                        # Check if the previously provided task does not match
                        # with the current task. If True, modify the progress
                        # dialog so that the previous task is finished and
                        # start the current task
                        if shared_ns.task != current_op[0]:
                            # Get widget of previous task
                            pwgt = wgt_ref[current_op[0]]
                            # Modify progress dialog of previous task so that
                            # it appears as finished
                            finish_op(pwgt, current_op[0])
                            # Set the previous task as the current task for
                            # future comparisons
                            current_op[0] = shared_ns.task
                            # Get widget of current task
                            cwgt = wgt_ref[shared_ns.task]
                            # Start loading progress of current task
                            start_op(cwgt)
                        else:
                            cwgt = wgt_ref[shared_ns.task]

                    # In this case, this is the first task being performed.
                    # We only need to start the loading process of the current
                    # task
                    else:
                        # Set the previous task as the current task
                        current_op.append(shared_ns.task)
                        # Get widget of current task
                        cwgt = wgt_ref[shared_ns.task]
                        # Start loading progress of current task
                        start_op(cwgt)

                    wgt_ref["spinner"].ids.img.rotation -= 10

                    if shared_ns.total and shared_ns.counter:
                        # Get percentage
                        try:
                            perc = int((float(shared_ns.counter) /
                                        float(shared_ns.total)) * 100.)
                        except TypeError:
                            perc = 0

                        # Set percentage on progress dialog
                        wgt_ref["counter"].text = "{}%".format(perc)

                        # Set main label, if any
                        if shared_ns.main_msg:
                            cwgt.ids.main_lbl.text = shared_ns.main_msg

                        # Set secondary label, if any
                        if shared_ns.msg:
                            try:
                                cwgt.ids.secondary_lbl.text = shared_ns.msg
                            except ValueError:
                                pass

            except AttributeError:
                pass

            # Interrupt subporcess on user demand
            if self.terminate_process_exec:
                shared_ns.stop = True
                man.shutdown()
                Clock.unschedule(check_func)
                self.dismiss_all_popups()
                # Removes all temporary database tables
                self.alignment_list.remove_tables()
                return

            try:
                content.ids.msg.text = shared_ns.msg
            except AttributeError:
                pass

            # Checks for signal that the current file being processed
            # already exists
            try:
                if shared_ns.file_dialog:
                    self.dialog_fileoverwrite("File %s already exists" %
                                              basename(shared_ns.file_dialog))
                    # This will prevent the dialog from being
                    # re-generated while the user chooses the option
                    shared_ns.file_dialog = None
            except AttributeError:
                pass

            # Transmit to the subprocess the option of overwrite/skip
            try:
                if self.file_overwrite:
                    shared_ns.status = self.file_overwrite
                    # If the user does not check the apply to all,
                    # reset the file_overwrite attribute for the next file
                    if not self.file_apply_all:
                        self.file_overwrite = None
                        shared_ns.status = None
                    else:
                        shared_ns.apply_all = True
            except AttributeError:
                pass

            if not p.is_alive():
                # Resets file_apply_all and file_overwrite attributes
                self.file_apply_all = False
                self.file_overwrite = None
                Clock.unschedule(check_func)
                self.dismiss_all_popups()

                # Removes all temporary database tables
                self.alignment_list.remove_tables()

                p.join()

                # If process execution ended with an error, issue warning.
                try:

                    self.dialog_warning(shared_ns.exception["exception"][0],
                                        shared_ns.exception["exception"][1])
                    man.shutdown()
                except:
                    if shared_ns.proc_files == 1:
                        self.dialog_floatcheck(
                            "All Done! %s file was successfully processed"
                            % shared_ns.proc_files, t="info")
                    else:
                        self.dialog_floatcheck(
                            "All Done! %s files were successfully processed"
                            % shared_ns.proc_files, t="info")

                    if shared_ns.filtered_alns:
                        self.dialog_filter_report(shared_ns.filtered_alns,
                                                  shared_ns.proc_files)
                    man.shutdown()

        # Set up manager and shared name space to share information
        # between background and main processes
        manager = multiprocessing.Manager()
        shared_ns = manager.Namespace()
        shared_ns.status = None
        shared_ns.apply_all = False

        shared_ns.stop = False
        shared_ns.task = shared_ns.total = shared_ns.counter = \
            shared_ns.msg = shared_ns.main_msg = None
        shared_ns.finished_tasks = []
        current_op = []

        # Packing arguments to background process
        process_kwargs = {"aln_list": self.alignment_list,
            "file_set_name": self.process_grid_wgt.ids.active_file_set.text,
            "file_list": list(self.file_list),
            "filename_map": dict(self.filename_map),
            "file_groups": dict(self.file_groups),
            "taxa_set_name": self.process_grid_wgt.ids.active_taxa_set.text,
            "active_taxa_list": list(self.active_taxa_list),
            "ns": shared_ns,
            "taxa_groups": dict(self.taxa_groups),
            "hap_prefix": str(self.hap_prefix),
            "secondary_operations": self.secondary_operations,
            "secondary_options": dict(self.secondary_options),
            "missing_filter_settings": list(self.missing_filter_settings),
            "taxa_filter_settings": list(self.taxa_filter_settings),
            "codon_filter_settings": list(self.codon_filter_settings),
            "variation_filter_settings": list(self.variation_filter),
            "output_file": self.output_file,
            "conversion_suffix": self.conversion_suffix,
            "rev_infile": self.rev_infile,
            "main_operations": dict(self.main_operations),
            "zorro_suffix": str(self.zorro_suffix),
            "partitions_file": self.partitions_file,
            "output_formats": list(self.output_formats),
            "create_partfile": bool(self.create_partfile),
            "use_nexus_partitions": bool(self.use_nexus_partitions),
            "use_nexus_models": bool(self.use_nexus_models),
            "phylip_truncate_name": bool(self.phylip_truncate_name),
            "output_dir": self.output_dir,
            "use_app_partitions": bool(self.use_app_partitions),
            "consensus_type": self.process_options.ids.consensus_mode.text,
            "ld_hat": bool(self.ld_hat),
            "ima2_params": list(self.ima2_options)}

        # Remove lock from background process
        self.terminate_process_exec = False

        # Create process
        p = threading.Thread(target=process_execution,
                     kwargs=process_kwargs)
        p.start()

        # Remove any possible previous popups
        self.dismiss_popup()
        # Create waiting dialog
        content = ProcessExecutionProgress()
        wgt_ref = build_progress_dialog(content)

        # Schedule function that checks the process' pulse
        check_func = partial(check_process, p, manager)
        Clock.schedule_interval(check_func, .1)


def main():

    multiprocessing.freeze_support()

    global MCL_FILE

    # This will handle the pointer to the kv file whether TriFusion is a
    # one-file executable, a one-dir executable, or running from source.
    if getattr(sys, "frozen", False):

        # Get relative path of mcl binary
        mcl_path = None
        if sys.platform == "linux2":
            mcl_path = os.path.join("data", "resources", "mcl", "linux", "mcl")
        elif sys.platform in ["win32", "cygwin"]:
            if platform.architecture()[0] == "64bit":
                mcl_path = os.path.join("data", "resources", "mcl", "windows",
                                        "64bit", "mcl64.exe")
            else:
                mcl_path = os.path.join("data", "resources", "mcl", "windows",
                                        "32bit", "mcl32.exe")
        elif sys.platform == "darwin":
            mcl_path = os.path.join("data", "resources", "mcl", "MacOS", "mcl")

        # One-file
        try:
            # Get main kv file
            kv_file = os.path.join(sys._MEIPASS, "trifusion.kv")
            os.chdir(sys._MEIPASS)
            # Get mcl executable path
            MCL_FILE = os.path.join(sys._MEIPASS, mcl_path)
        # One-dir
        except AttributeError:
            kv_file = "trifusion.kv"
            MCL_FILE = mcl_path

    # Source
    elif __file__:
        kv_file = "trifusion.kv"

    TriFusionApp(kv_file="trifusion.kv").run()

if __name__ == "__main__":
    main()
