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

# Kivy imports
from kivy.app import App
from kivy.uix.togglebutton import ToggleButton
from kivy.uix.button import Button
from kivy.animation import Animation
from kivy.uix.image import Image
from kivy.uix.popup import Popup
from kivy.uix.widget import Widget
from kivy.uix.label import Label
from kivy.core.window import Window
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.gridlayout import GridLayout
from kivy.uix.tabbedpanel import TabbedPanel
from kivy.uix.floatlayout import FloatLayout
from kivy.uix.textinput import TextInput
from kivy.uix.anchorlayout import AnchorLayout
from kivy.uix.spinner import Spinner
from kivy.uix.filechooser import FileChooserListView, FileChooserIconView
from kivy.uix.checkbox import CheckBox
from kivy.lang import Builder
from kivy.properties import NumericProperty, StringProperty, BooleanProperty,\
    ListProperty, ObjectProperty, DictProperty
from kivy.uix.screenmanager import Screen
from kivy.config import Config
from kivy.clock import Clock
from kivy.uix.treeview import TreeView, TreeViewLabel
from kivy.graphics import Color, Rectangle
from kivy.uix.scrollview import ScrollView

# Main program imports
from ortho import orthomcl_pipeline as opipe
from ortho import OrthomclToolbox as ot
from process.sequence import AlignmentList
from process.base import Base
from process import data
from data.resources.info_data import informative_storage
from data.resources.db_tools import *

# Other imports
import os
from os.path import dirname, join, exists, abspath, pardir, expanduser
from os import sep
from collections import OrderedDict
from copy import deepcopy
from functools import partial
import pickle
import multiprocessing
import time
import tempfile
import shutil
import re

Config.set("kivy", "log_level", "warning")
Config.set("kivy", "desktop", 1)
Config.set("kivy", "exit_on_escape", 0)
Config.set("graphics", "resizable", 0)
Config.set("graphics", "fullscreen", 0)

#### MEMO
# scatter.py on_touch_up function was modified to prevent a bug from crashing
# the app. Before the fix, when right-clicking in a ScatterLayout and then
# left-clicking in the generate ball would crash the app with a KeyError
# exception. The fix on lines 595-598 handles this exception by not deleting
# a touch from _last_touch_pos dictionary, since the key is not there anymore.

class ShowcaseScreen(Screen):
    fullscreen = BooleanProperty(False)

    def add_widget(self, *args):
        if 'content' in self.ids:
            return self.ids.content.add_widget(*args)
        return super(ShowcaseScreen, self).add_widget(*args)


class StringInput(TextInput):
    """
    Modification of TextInput that only accepts integers
    """

    def insert_text(self, substring, from_undo=False):

        s = re.sub(r"\D", "", substring)

        return super(StringInput, self).insert_text(s, from_undo=from_undo)


class FileChooserL(FileChooserListView):
    """
    Modification of the FileChooserListView widget that fixes an issue of path
    update when clicking in the parent directory
    """

    def __init__(self, **kwargs):
        super(FileChooserL, self).__init__(**kwargs)

    def open_entry(self, entry):
        """
        Modification of the open entry method so that when the entry.path is
        "../", the path is updated to the parent directory, instead of
        appending the entry.path
        """
        try:
            # Just check if we can list the directory. This is also what
            # _add_file does, so if it fails here, it would also fail later
            # on. Do the check here to prevent setting path to an invalid
            # directory that we cannot list.
            self.file_system.listdir(entry.path)
        except OSError:
            entry.locked = True
        else:
            # If entry.path is to jump to previous directory, update path with
            # parent directory
            if entry.path == "../":
                self.path = abspath(join(self.path, pardir))
                self.selection = []
            else:
                self.path = join(self.path, entry.path)
                self.selection = []


class FileChooserM(FileChooserIconView):
    """
    Modification of the FileChooserIconView widget that fixes n issue of path
    update when clicking in the parent directory and provides support for
    multiple file selection using Shift+Click. To achieve this, a few methods
    were added that listen to keyboard input in order to capture when the
    shift key is being pressed. These methods change the new shift attribute of
    the class, which is used when an entry is touched.

    The current Shift+Click implementation supports forward and backward
    selection from the last entry touched.
    """

    shift = False

    def __init__(self, **kwargs):
        super(FileChooserM, self).__init__(**kwargs)
        Window.bind(on_key_down=self.keyboard_listen)
        Window.bind(on_key_up=self.release_shift)

    def keyboard_listen(self, *vals):
        """
        Listens to keyboard when a key is pressed. It is used to set the shift
        attribute to True when Shift is being pressed
        :param vals: keyboard input
        """

        key_code = vals[1:3]
        if key_code == (304, 50):
            self.shift = True

    def release_shift(self, *vals):
        """
        Listens to keyboard when a key is released. It is used to set the
        shift attribute to False when the Shift key is released.
        :param vals:
        :return:
        """

        key_code = vals[1:3]
        if key_code == (304, 50):
            self.shift = False

    def open_entry(self, entry):
        """
        Modification of the open entry method so that when the entry.path is
        "../", the path is updated to the parent directory, instead of
        appending the entry.path
        """
        try:
            # Just check if we can list the directory. This is also what
            # _add_file does, so if it fails here, it would also fail later
            # on. Do the check here to prevent setting path to an invalid
            # directory that we cannot list.
            self.file_system.listdir(entry.path)
        except OSError:
            entry.locked = True
        else:
            # If entry.path is to jump to previous directory, update path with
            # parent directory
            if entry.path == "../":
                self.path = abspath(join(self.path, pardir))
                self.selection = []
            else:
                self.path = join(self.path, entry.path)
                self.selection = []

    def entry_touched(self, entry, touch):
        """
        (internal) This method must be called by the template when an entry
        is touched by the user. Supports Shift+Clicking for multiple selection
        """
        if (
            'button' in touch.profile and touch.button in (
                'scrollup', 'scrolldown', 'scrollleft', 'scrollright')):
            return False

        _dir = self.file_system.is_dir(entry.path)
        dirselect = self.dirselect

        if _dir and dirselect and touch.is_double_tap:
            self.open_entry(entry)
            return

        if self.shift and self.selection:
            # Get index of last selection entry and current entry
            idx_selection = self.files.index(self.selection[-1])
            idx_current = self.files.index(entry.path)

            # If current entry is ahead of last selection, select files
            # going forward
            if idx_selection < idx_current:
                idx_s = idx_selection
                idx_f = idx_current
            else:
                idx_s = idx_current
                idx_f = idx_selection

        if self.multiselect:
            if entry.path in self.selection:
                # This will deselect multiple files when the shift key is down
                # while clicking
                if self.shift and self.selection:
                    for i in range(idx_s + 1, idx_f + 1):
                        f = self.files[i]
                        if f in self.selection:
                            self.selection.remove(f)
                else:
                    self.selection.remove(entry.path)
            else:
                if _dir and not self.dirselect:
                    self.open_entry(entry)
                    return
                # This will select multiple files when the shift key is down
                # while clicking
                if self.shift and self.selection:
                    for i in range(idx_s, idx_f + 1):
                        f = self.files[i]
                        if f not in self.selection:
                            self.selection.append(f)
                else:
                    self.selection.append(entry.path)
        else:
            if _dir and not self.dirselect:
                self.open_entry
                return
            self.selection = [entry.path, ]


class CustomPopup(Popup):
    """
    Modification of Popup class with a few additional feature.

    .: The title does not wrap, but instead is shortened
    .: A custom background may be provided using the custom_background attribute
    """

    def __init__(self, **kwargs):
        super(CustomPopup, self).__init__(**kwargs)
        label = self.children[0].children[-1]
        label.shorten = True
        label.shorten_from = "right"
        label.markup = True

        # New attributes
        try:
            self.custom_background = kwargs["custom_background"]
        except KeyError:
            self.custom_background = None

        # Set custom background
        if self.custom_background:
            gl = self.children[0]
            with gl.canvas.before:
                Color(.7, .7, .7, 1)
                self.rect = Rectangle(
                    source=self.custom_background,
                    pos=self.pos,
                    size=self.size)

                # This will update the background position when the popup is
                # set to the final position
                self.bind(size=self._update_rect, pos=self._update_rect)

    def _update_rect(self, instance, value):
        self.rect.pos = instance.pos
        self.rect.size = instance.size


class AutoCompTextInput(TextInput):
    """
    Modified widget of text input in which the tab key does not introduce a
    tabular space. This is meant to use with _auto_completion, in which the
    tab key serves as a keybinding
    """

    def insert_text(self, substring, from_undo=False):
        if substring == "\t":
            s = ""
        else:
            s = substring
        return super(AutoCompTextInput, self).insert_text(s,
                                                          from_undo=from_undo)


class PartitionsDialog(BoxLayout):
    """
    Custom layout for partition box when editing partitions
    """
    pass


class ModelSpinner(Spinner):
    """
    Custom Spinner that takes a background_normal argument to set the
    background
    """

    def __init__(self, **kwargs):
        super(ModelSpinner, self).__init__(**kwargs)

        try:
            self.background_normal = kwargs["background_normal"]
        except KeyError:
            pass


class ExportGraphics(BoxLayout):
    cancel = ObjectProperty(None)


class PlotToolbar(BoxLayout):
    pass


class OrthologySearchGrid(TabbedPanel):
    pass


class DescriptionBox(BoxLayout):
    # Attribute for number of proteins
    prot_txt = StringProperty()
    # Attribute for number of taxa
    taxa_txt = StringProperty()
    # Attribute for total number of orthologs
    ortholog_txt = StringProperty()


class GaugePlot(BoxLayout):
    # Attribute for Gauge plot top label
    txt = StringProperty()
    # Attribute for proportion for gauge plot. This proportion should range
    # between 0 and 1. It will be automatically adapted in the gauge plot
    proportion = NumericProperty()
    # Attribute for number of orthologs
    ortholog_num = StringProperty()


class OrtoSetFiltersDialog(BoxLayout):
    cancel = ObjectProperty(None)


class OrthoReportDialog(BoxLayout):
    cancel = ObjectProperty(None)


class OrthoGraphicReport(BoxLayout):
    pass


class OrtoFilterDialog(BoxLayout):
    cancel = ObjectProperty(None)


class OrtoExecutionDialog(BoxLayout):
    cancel = ObjectProperty(None)


class OrtoProgressDialog(BoxLayout):
    pass


class ProteinFilterDialog(BoxLayout):
    cancel = ObjectProperty(None)


class MySQLDialog(BoxLayout):
    cancel = ObjectProperty(None)


class InflationDialog(BoxLayout):
    cancel = ObjectProperty(None)


class LoadMultipleDialog(BoxLayout):
    """
    A Filechooser widget in Icon view that allows multpiple file choosing
    """
    cancel = ObjectProperty(None)


class CloseBox(BoxLayout):
    """
    This is part of the taxa information popup. It contains the closing button
    """
    cancel = ObjectProperty(None)


class RemoveFloat(Button):
    """
    Simple (X) float button that can be associated with several root_window
    widgets for closing buttons
    """
    pass


class WarningFloat(Label):
    """
    The general purpose unintruside warning float for errors and informations.
    This dialog is added to the root_window with fade in and fade out animations
    """
    pass


class InfoPopup(BoxLayout):
    """
    Dialog for help texts that can be accessed for several options using the
    "?" button
    """
    cancel = ObjectProperty(None)


class MouseOverLabel(Button):
    """
    General use mouse over label for diverse buttons
    """
    pass


class RevConcDialog(BoxLayout):
    """
    Reverse concatenation dialog
    """
    cancel = ObjectProperty(None)


class InputList(BoxLayout):
    """
    Dialog with list of input files to select file for reverse concatenation
    """
    cancel = ObjectProperty(None)


class SideLabel(Label):
    """
    Mouseover label for side option buttons in the sidepanel
    """
    pass


class LoadMoreBt(AnchorLayout):
    """
    Custom button widget for for the "load more button" in the side panel
    """
    pass


class PathLabel(Label):
    """
    Dialog for the Label with the path for the main file chooser
    """
    pass


class PathText(AutoCompTextInput):
    """
    Dialog for the TextInput for the main file chooser that controls the path
     and inherits from the custom AutoCompTextInput for auto completion
    """
    pass


class ExecutionDialog(BoxLayout):
    """
    The Execution dialog for Process operations
    """
    cancel = ObjectProperty(None)


class ProteomePopup(BoxLayout):
    """
    Informative popup for proteome files when clicking the "i" button in the
    side panel
    """
    cancel = ObjectProperty(None)


class LoadProgressDialog(BoxLayout):
    """
    Dialog for the progression dialog when loading files into the program
    """
    pass


class FilePopup(BoxLayout):
    """
    Class with a custom BoxLayout controlling the informative popup for the
     file buttons in the File tab of the side panel
    """
    cancel = ObjectProperty(None)


class TaxaPopup(BoxLayout):
    """
    Class with a custom BoxLayout controlling the informative popup for the
     taxa buttons in the Taxa tab of the side panel
    """
    pass


class CheckDialog(BoxLayout):
    """
    Class controlling the layout of a general purpose dialog to check if the
    user wants of perform a certain action
    """
    cancel = ObjectProperty(None)


class WarningDialog(BoxLayout):
    """
    Class controlling the layout of a general purpose dialog to warn the user
    of certain events
    """
    cancel = ObjectProperty(None)


class LoadDialog(BoxLayout):
    """
    Class controlling a general purpose layout for loading additional files
    """
    cancel = ObjectProperty(None)


class SaveDialog(FloatLayout):
    """
    Class controlling the layout of the save file dialog in the Process screen
    """
    save = ObjectProperty(None)
    text_input = ObjectProperty(None)
    cancel = ObjectProperty(None)


class FormatDialog(BoxLayout):
    """
    Class controlling the layout of the output format dialog in the Process
    screen
    """
    cancel = ObjectProperty(None)


class FilterDialog(BoxLayout):
    """
    Class controlling the layout of the gap/missing filtering options in the
    Process screen
    """
    cancel = ObjectProperty(None)


class TextDialog(BoxLayout):
    """
    Class controlling a simple text input popup
    """
    cancel = ObjectProperty(None)


class NexusExtra(BoxLayout):
    cancel = ObjectProperty(None)


class PhylipExtra(BoxLayout):
    """
    Class controlling the dialog with extra options for phylip output format
    """
    cancel = ObjectProperty(None)


class ProcessGeneral(GridLayout):
    """
    Class controlling the layout of the general options of the Process screen
    """
    pass


class AdditionalProcessContents(TabbedPanel):
    """
    Class controlling the layout of the additional options of the Process screen
    """
    pass


class TaxaGroupDialog(BoxLayout):
    """
    Class controlling the layout of the taxa group creation dialogue in the
    side panel
    """
    cancel = ObjectProperty(None)


class ZorroDialog(BoxLayout):
    """
    Class controlling the layout of the ZORRO operation dialog
    """
    cancel = ObjectProperty(None)

#===============================================================================
#                                  EXCEPTIONS
#===============================================================================


class InputTypeError(Exception):
    pass


class TriFusionApp(App):

    #######################
    #
    # GUI RELATED VARIABLES
    #
    #######################

    # Setting Boolean controlling the toggling of main headers
    show_side_panel = BooleanProperty(False)

    # Attribute for current screen object
    screen = None

    # Variable containing screen names
    screen_names = ListProperty()
    available_screens = ListProperty()
    loaded_screens = {}
    plot_screens = ListProperty()

    # Attributes to know current and previous screen
    current_screen = StringProperty()
    previous_screen = StringProperty()

    # Attribute to load screens
    index = NumericProperty(-1)

    # Create temporary directory for transient files
    temp_dir = tempfile.TemporaryDirectory()

    # Getting current directory to fetch the screen kv files
    cur_dir = dirname(__file__)

    # Only the original input files. SHOULD NOT BE MODIFIED
    file_list = ListProperty()
    # Dynamic list containing only the activated files
    active_file_list = ListProperty()
    # Dictionary mapping file names to their corresponding full paths. This
    # attribute must exist, because some parts of the code need only the file
    # name instead of the full path, but a connection to the full path must
    # be maintained for future reference
    filename_map = DictProperty()

    # Setting the list of taxa names
    active_taxa_list = ListProperty()

    # Attribute to the home path
    home_path = expanduser("~")

    # Attribute with bookmarks for file chooser. The bookmark attribute is a
    # list containing a list with the full path of the bookmarks as the
    # first element and a dictionary mapping the full path to the bookmark
    # name as the second element
    bookmarks = [[], {}]
    bm_file = join(cur_dir, "data", "resources", "bookmarks")
    # For mouse over purposes
    bookmarks_bt = []

    _popup = ObjectProperty(None)
    _subpopup = ObjectProperty(None)

    # Dictionary containing the values for the main process operations
    main_operations = DictProperty({"concatenation": False, "conversion": False,
                       "reverse_concatenation": False})

    # Dictionary containing all values of the switches and checkboxes in the
    # process screen
    secondary_operations = OrderedDict([("collapse", False),
                                    ("filter", False),
                                    ("gcoder", False)])

    secondary_options = OrderedDict([("interleave", False),
                                    ("zorro", False),
                                    ("collapse_file", False),
                                    ("filter_file", False),
                                    ("gcoder_file", False)])

    # Attributes for the Orthology screen widgets
    ortho_search_options = None
    orto_search_height = None

    # Attribute for the gridlayout widget that will contain all main options
    # for the process module
    process_grid_wgt = None
    process_options = None
    process_height = None

    # Attribute for the widget containing the treeview showing the operations
    # queue
    operation_tv = ObjectProperty(None)
    main_nodes = DictProperty()

    # Attributes containing plot related elements
    current_plot = ObjectProperty(None)
    current_lgd = ObjectProperty(None)

    # Attributes for storing taxa and file buttons for side panel. These will
    # be used when search for files/taxa and for loading only button subsets
    # for very large data sets. Each list element pertains a single file/taxon
    # and it will be a tupple containing the main button, information button
    # and remove button.
    sp_file_bts = ListProperty()
    sp_taxa_bts = ListProperty()
    sp_partition_bts = ListProperty()

    # Attributes that control the amount of taxa/file buttons showing at the
    # side panel. To avoid staggering the app with tons of buttons, a
    # maximum number of buttons showing initially is set. More buttons
    # can be later added.
    MAX_FILE_BUTTON = NumericProperty(20)
    count_files = NumericProperty(0)
    MAX_PARTITION_BUTTON = NumericProperty(20)
    count_partitions = NumericProperty(0)

    # Attributes storing the toggle buttons from Taxa/File panels. Mostly for
    # mouse_over events
    # Contains the button widgets from the Files and Taxa tabs
    mouse_over_bts = DictProperty({"Files": [], "Taxa": [], "Partitions": []})
    # The button text of the previous mouse over event. This will allow the
    # assessment of whether the current mouse collision is for the same button
    # (in which case the mouse over will not be triggered) or for a different
    # button (in which case the mouse over is triggered)
    previous_mouse_over = StringProperty("")
    # This is a locking mechanism of the mouse over event. When there is a
    # scheduled event for a mouse over this attribute is set to False, which
    # prevents further events from being scheduled in the meantime. When the
    # scheduled event is dispatched, the lock is released and it returns to
    # True
    mouse_over_ready = BooleanProperty(True)
    # Stores the previous mouse over label button so that it can be removed
    old_mouse_over = None

    ################################
    #
    # CORE PROGRAM RELATED VARIABLES
    #
    ################################

    # MySQL access
    mysql_pass = StringProperty("")

    # OrthoMCL output directory
    ortho_dir = StringProperty("")

    # USEARCH database
    usearch_db = StringProperty("goodProteins_db")
    usearch_output = StringProperty("AllVsAll.out")
    usearch_evalue = StringProperty("0.00001")

    # MCL/Groups attributes
    ortholog_prefix = StringProperty("My_group")
    group_prefix = StringProperty("group")
    mcl_inflation = ListProperty(["3"])

    # Protein quality filters
    protein_min_len = NumericProperty(10)  # Absolute
    protein_max_stop = NumericProperty(20)  # Percentage

    # Orthology cluster filters
    orto_max_gene = NumericProperty(1)
    orto_min_sp = NumericProperty(3)

    # Attribute containing the path to the proteome files
    proteome_files = ListProperty()
    active_proteome_files = ListProperty()

    # Attribute containing the orthology group files
    ortho_groups = ObjectProperty(None)

    # List storing the original alignment object variables. SHOULD NOT BE
    # MODIFIED
    alignment_list = ObjectProperty(None)
    # List of active alignment object variables.
    active_alignment_list = ObjectProperty(None)

    # Attributes containing the original and active taxa information
    # dictionaries
    original_tx_inf = DictProperty()
    active_tx_inf = DictProperty()
    # Same as before but for file information
    original_file_inf = DictProperty()
    active_file_inf = DictProperty()

    # Export mode. Tuple with first element "taxa" or "file" and second element
    # as "all" or "selected"
    export_mode = None

    # Attribute storing the sequence types currently loaded
    sequence_types = ListProperty()

    # Attribute for taxa and file groups
    taxa_groups = DictProperty()
    file_groups = DictProperty

    # Attribute containing the objects for the several possible output files.
    output_file = StringProperty("")
    output_dir = StringProperty("")

    # Attribute storing active output formats. Fasta is True by default
    output_formats = ListProperty(["fasta"])

    # Attributes for extra options of output formats
    # Determines whether the part.File associated with phylip format is created
    create_partfile = BooleanProperty(True)
    # Determines whether the charset partitions in Nexus input files are to
    # be used in the output file
    use_nexus_partitions = BooleanProperty(True)

    # Attribute storing the filter settings. The list should contain gap
    # threshold as first element and missing data threshold as second element
    filter_settings = ListProperty([25, 50])

    # Partitions file
    partitions_file = StringProperty("")
    # Input file for reverse concatenation
    rev_infile = StringProperty("")

    # Attribute for ZORRO settings
    zorro_suffix = StringProperty("")

    # Attribute storing the haplotype prefix
    hap_prefix = StringProperty("Hap")

    ##################################
    #
    # GUI RELATED METHODS AND FUNCTIONS
    #
    ##################################

    def build(self):

        # Setting main window title
        self.title = "TriFusion - Streamline phylogenomics"

        # Setting available screens
        self.available_screens = ["main", "Orthology", "Process",
                                  "Statistics", "fc", "group_compare"]
        self.screen_names = self.available_screens

        # Transforming screen names into complete paths to be loaded by kivy
        self.available_screens = [join(self.cur_dir, "data", "screens",
                                 "{}.kv".format(screen)) for screen in
                                  self.available_screens]

        # Store screen names specifically designed for plot display
        self.plot_screens = ["group_compare"]

        self.loaded_screens = dict((sc, None) for sc in self.available_screens)

        # First thing is go to main screen
        self.go_screen(0)

        # Set method for closing side panel when touching outside
        Window.bind(on_touch_up=self.sidepanel_on_touch)

        # Listen to keybindings
        Window.bind(on_key_down=self._on_keyboard_events)

        # Execute cleaning function when exiting app
        Window.bind(on_request_close=lambda x: self._exit_clean())

        #### Orthology widgets
        self.ortho_search_options = OrthologySearchGrid()

        #### Process widgets
        # Creating GridLayout instance for general options of Process
        self.process_grid_wgt = ProcessGeneral()
        # Create GridLayout instance for additional options of Process.
        self.process_options = AdditionalProcessContents()

        # Initialize operation queue treeview in side panel
        self.operation_queue_init()

        # Set schedule for mouse over events on side panel
        Clock.schedule_interval(lambda x: self._on_mouseover_tabs(), .1)

        # This corrects a weird bug where the buttons in the dropdown of the
        # taxa tab in the side panel are open by default.
        self.root.ids.taxa_dd.dismiss()
        self.root.ids.file_dd.dismiss()

        """
        ------------------------ METHOD NOMENCLATURE GUIDE ---------------------

        Given the large number of methods needed to give functionality to the
        app, this nomenclature guide was created to aid in the naming of new
        methods so that the code can be more easily browsed and understood. Note
        that this guide only targets methods that perform similar tasks and,
        therefore, can be grouped by a common prefix name. Other methods that
        perform more unique operations may have different names.

        Method's names will be given based on their main operation and specific
        task. For example, a method in charge of toggle the side panel, should
        be named "toggle_sidepanel", being "toggle" the common prefix and
        "sidepanel" the keyword linked ot the specific task.

        1. Toggles.

        "toggle_[specific_task]", e.g. "toggle_sidepanel"

        Methods use to toggle certain widgets or values/attributes.

        2. Dialogues.

        "dialog_[specific_task]", e.g. "dialog_format"

        Methods that generate dialogues throughout the app, usually in the form
        of popups

        3. Populating methods.

        "populate_[specific_task]", e.g., "populate_input_files"

        Methods that populate certain widgets, usually gridlayouts, with other
        widgets

        4. Add/Remove

        "add_[specific_task]", e.g., "add_bookmark"
        "remove_[specific_task]", e.g., "remove_taxa_group"

        Methods that add or remove widgets, usually buttons/togglebuttons, from
        other widgets

        5. Saves.

        "save_[specific_task]", e.g., "save_file"

        Methods that save specific settings from the several options of the app

        6. Updates.

        "update_[specific_task]", e.g., "update_tabs"

        Wrapper methods used to update several attributes or widgets of the app

        7. Checks.

        "check_[specific_task]", e.g., "check_filters"

        Methods that perform some kind of sanity checks to user input data

        8. Unique operations

        [specific_task]_[unique_operation], e.g., "sidepanel_animation"

        When the method performs a unique operations, the specific_task should
        prefix the name of the method.
        """

    def _exit_clean(self):
        """
        This method is issued when the application is closed and performs any
        necessary clean up operations
        """

        shutil.rmtree(self.temp_dir.name)

    def _on_keyboard_events(self, *vals):
        """
        Methods that listens to keyboard input and triggers events or changes
        properties when acceptable keyboard shortcuts are entered
        :param vals: input list from on_key_down function
        """

        # Get modifier (usually ctrl or shift)
        # TODO: The modifier in MacOS is different. Must check on this.
        modifier = "".join(vals[-1])
        # Get key
        key = vals[-2].encode("utf-8")
        key_code = vals[1:3]

        #=======================================================================
        # Popup keybindings
        #=======================================================================

        def popup_keys(backn, backd, bt1, bt2, bt3=None):
            """
            Wrapper function that provides functionality to arrow keys for
            navigating through selection buttons in popups
            :param backn: string for background normal
            :param backd: string for background down
            :param bt1: Button widget one (usually for Ok buttons)
            :param bt2: Button widget two (usually for Cancel buttons)
            """

            # This will deal with cases with only two buttons to cycle
            if not bt3:
                # if left arrow key
                if key_code == (276, 113):
                    bt1.background_normal = backn
                    bt2.background_normal = backd
                # if right arrow key
                if key_code == (275, 114):
                    bt1.background_normal = backd
                    bt2.background_normal = backn
                # if enter key. Dispatch the events of the focused button
                if key_code == (13, 36):
                    if bt1.background_normal == backn:
                        bt1.dispatch("on_release")
                    else:
                        bt2.dispatch("on_release")

            # This will cycle through three buttons
            else:
                bt_list = [bt1, bt2, bt3]
                idx = [x.background_normal for x in bt_list].index(backn)
                if key_code == (276, 113) and idx > 0:
                    idx -= 1
                if key_code == (275, 114) and idx < 2:
                    idx += 1

                for bt in bt_list:
                    if bt_list.index(bt) == idx:
                        bt_list[idx].background_normal = backn
                    else:
                        bt.background_normal = backd

                if key_code == (13, 36):
                    bt_on = [x for x in bt_list if
                             x.background_normal == backn][0]
                    bt_on.dispatch("on_release")

        # Check only when a popup is active
        if self._subpopup in self.root_window.children:
            if "ok_bt" in self._subpopup.content.ids:
                bn = join("data", "backgrounds", "bt_process.png")
                bd = join("data", "backgrounds", "bt_process_off.png")
                ok_bt = self._subpopup.content.ids["ok_bt"]
                cancel_bt = self._subpopup.content.ids["cancel_bt"]
                popup_keys(bn, bd, ok_bt, cancel_bt)

        elif self._popup in self.root_window.children:
            if "check_ok" in self._popup.content.ids:
                bn = join("data", "backgrounds", "check_ok.png")
                bd = join("data", "backgrounds", "check_cancel.png")
                ok_bt = self._popup.content.ids["check_ok"]
                cancel_bt = self._popup.content.ids["check_cancel"]
                popup_keys(bn, bd, ok_bt, cancel_bt)
            # In this case there are three buttons to cicle
            elif "apply_bt" in self._popup.content.ids:
                bn = join("data", "backgrounds", "bt_process.png")
                bd = join("data", "backgrounds", "bt_process_off.png")
                ok_bt = self._popup.content.ids["ok_bt"]
                cancel_bt = self._popup.content.ids["cancel_bt"]
                apply_bt = self._popup.content.ids["apply_bt"]
                popup_keys(bn, bd, ok_bt, apply_bt, cancel_bt)
            # Only two buttons to cicle
            elif "ok_bt" in self._popup.content.ids:
                bn = join("data", "backgrounds", "bt_process.png")
                bd = join("data", "backgrounds", "bt_process_off.png")
                ok_bt = self._popup.content.ids["ok_bt"]
                cancel_bt = self._popup.content.ids["cancel_bt"]
                popup_keys(bn, bd, ok_bt, cancel_bt)

            if "close_bt" in self._popup.content.ids:
                if key_code == (13, 36):
                    self._popup.content.ids.close_bt.dispatch("on_release")

        #=======================================================================
        # Filechooser keybindings
        #=======================================================================

        if self.screen.name == "fc":
            # Keybinding ctrl+f that brings focus to the "Find" field in the
            # Filechooser screen
            if modifier == "ctrl" and key == b'\x06':
                self.screen.ids.text_filter.focus = True

            # Keybinding ctrl+backspace to clear selection
            if modifier == "ctrl" and key_code == (8, 22):
                self.screen.ids.clear_s.dispatch("on_release")

            # Add bookmarks with ctrl+d
            if modifier == "ctrl" and key_code == (100, 40):
                self.screen.ids.add_bk_bt.dispatch("on_release")

            # Toggle manual path writing with ctrl+l
            if modifier == "ctrl" and key_code == (108, 46):
                if self.screen.ids.path_toggle.state == "down":
                    self.screen.ids.path_toggle.state = "normal"
                else:
                    self.screen.ids.path_toggle.state = "down"
                self.screen.ids.path_toggle.dispatch("on_release")

            # Select all files with ctrl+a
            if modifier == "ctrl" and key_code == (97, 38):
                self.screen.ids.icon_view_tab.selection = \
                    [x for x in self.screen.ids.icon_view_tab.files if not
                     os.path.isdir(x)]

            # Use arrow keys and enter to navigate through open/cancel buttons
            # and selecting them
            bn = join("data", "backgrounds", "bt_process.png")
            bd = join("data", "backgrounds", "bt_process_off.png")
            open_close_bt = self.screen.ids.open_close_bt
            open_bt = self.screen.ids.open_bt
            cancel_bt = self.screen.ids.cancel_bt
            popup_keys(bn, bd, open_close_bt, open_bt, cancel_bt)

        #=======================================================================
        # General keybindings
        #=======================================================================

        # Keybinding ctrl+o that opens the Filechooser screen
        if modifier == "ctrl" and key == b'\x0f':
            self.go_screen(self.screen_names.index("fc"))

        # Changing main screens between Orthology, Process and Statistics
        if modifier == "ctrl" and key_code == (49, 10):
            self.root.ids.h_ortho.dispatch("on_release")
            self.root.ids.h_ortho.state = "down"
        if modifier == "ctrl" and key_code == (50, 11):
            self.root.ids.h_process.dispatch("on_release")
            self.root.ids.h_process.state = "down"
        if modifier == "ctrl" and key_code == (51, 12):
            self.root.ids.h_stat.dispatch("on_release")
            self.root.ids.h_stat.state = "down"

        # Toggle side panel (slash)
        if key_code == (92, 49):
            self.root.ids.ap.dispatch("on_release")

        #=======================================================================
        # Text input autocompletion
        #=======================================================================

        # Use tab for auto completion when textinput is focused
        if key_code == (9, 23):
            if "path_bx" in self.screen.ids:
                if isinstance(self.screen.ids.path_bx.children[0], TextInput):
                    path = self.screen.ids.path_bx.children[0].text
                    s = self._auto_completion(path)
                    self.screen.ids.path_bx.children[0].text = s

    def _common_path(self):
        """
        Returns a string with the longest common path contained in the
        self.file_list dir. This is used by FileChoosers to open the nearest
        directory to their input files
        """

        # Return None if there are no input files
        if not self.file_list:
            return

        # Get common path
        common_path = ""
        for char in zip(*self.file_list):
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

        # If there are multiple matches in dirlist, return the longest common
        # substring
        elif len(dirlist) > 1:
            return join(path, os.path.commonprefix(dirlist))

        else:
            self.dialog_floatcheck("WARNING: Path does not exist",
                                   t="error")
            return original_path

    def _on_mouseover_tabs(self):
        """
        Provides mouse over events throughout the app. In order to reduce
        computations in each cycle, the mouse over events should only be tested
        when the appropriate screen/widget is visible
        """

        # Get mouse position coordinates
        mp = self.root_window.mouse_pos
        # Set collision attribute
        collision = False
        # Set side button list
        sidebt_list = [x for x in self.root.ids.side_bt.children if
                       isinstance(x, ToggleButton)]

        def show_label(mouse, wgt):
            """
            Use this function with a Clock schedule to delay the introduction
            of the label widget. Otherwise, it could become cumbersome to have
            an label appearing right after the mouse colliding with the button
            :param wgt: Layout widget, containing a descriptive label
            """

            # Checking if the current mouse position is the same as the mouse
            # position when the mouse over was triggered. This ensures that if
            # the user changes the mouse position while this event is
            # scheduled, the label will not be added to the root_window but the
            # lock in self.mouse_over_ready is removed
            if self.root_window.mouse_pos == mouse:
                # Add widget to root layout
                self.root_window.add_widget(wgt)
                # Update old label widget
                self.old_mouse_over = wgt
                # Update old label text
                self.previous_mouse_over = wgt.text

            # Unlocking mouse over
            self.mouse_over_ready = True

        def determine_collision(wgt):
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

        def create_label(text):
            """
            Creates a general use label
            :param text: Text to appear on the button
            :return: Button object
            """

            label_width = (len(text) + 7) * 7

            if label_width > self.root.width * .7:
                info_bt = MouseOverLabel(text=text, pos=mp,
                                size=(self.root.width * .6, 40),
                                size_hint=(None, None))
            else:
                info_bt = MouseOverLabel(text=text, pos=mp,
                                        size=(label_width, 40),
                            size_hint=(None, None))

            return info_bt

        def create_sidebt_wgt(text, p, s):
            """
            Creates the label for the sidebt mouseover
            :param text: string. Text for the label
            """

            side_bt = SideLabel(text=text, pos=p, size_hint=(None, None),
                             size=s, bold=True, border=(0, 0, 0, 0))

            return side_bt

        def create_fancy_label(text, wgt, lbl_height=30):
            """
            Creates a fancy label. This is akin to the mouse over
            labels in github, for example
            """

            # Create label
            lbl_wgt = Button(text=text, size_hint=(None, None), height=30,
                             background_disabled_normal=join("data",
                                "backgrounds", "fancy_mo_background.png"),
                             disabled=True, disabled_color=(1, 1, 1, 1),
                             bold=True, id=text)
            # Update label texture size, so that we can evaluate the available
            # space for the label
            lbl_wgt.texture_update()

            # Determine label size. Add horizontal margin space (10)
            lbl_wgt.size = (lbl_wgt.texture_size[0] + 10, lbl_height)

            # Determine if the label has space to the right. If not, flip the
            # orientation
            if wgt.pos[0] + lbl_wgt.width < self.root.width:
                # Determine position of arrow widget
                point_pos = wgt.pos[0] + wgt.width, wgt.pos[1] + wgt.height / \
                            2 - 6

                # Determine label position
                lbl_wgt.pos = point_pos[0] + 7,\
                          wgt.pos[1] + wgt.height / 2 - lbl_wgt.height / 2

                # Create arrow widget with left arrow
                point_wgt = Button(background_normal=join("data",
                                                        "backgrounds",
                                                        "box_arrow_right.png"),
                                   pos=point_pos, size=(7, 12), id=text,
                                   size_hint=(None, None),
                                   border=(0, 0, 0, 0))

            else:

                # Determine position of arrow widget
                point_pos = wgt.pos[0] - 10, wgt.pos[1] + wgt.height / 2 - 6

                # Determine label position
                lbl_wgt.pos = point_pos[0] - lbl_wgt.width,\
                          wgt.pos[1] + wgt.height / 2 - lbl_wgt.height / 2

                # Create arrow widget with left arrow
                point_wgt = Button(background_normal=join("data",
                                                          "backgrounds",
                                                          "box_arrow_left.png"),
                                   pos=point_pos, size=(7, 12), id=text,
                                   size_hint=(None, None),
                                   border=(0, 0, 0, 0))

            return point_wgt, lbl_wgt

        # Only do this routine when the filechooser screen is on
        if self.screen.name == "fc" and self.mouse_over_ready and \
                self.show_side_panel is False:
            case_bt = self.screen.ids.case_bt
            for bt in self.bookmarks_bt + [case_bt]:
                if determine_collision(bt):
                    collision = True
                    if bt == case_bt:
                        if "Case sensitive" not in self.previous_mouse_over:
                            if case_bt.state == "down":
                                label = create_label("Case sensitive is ON")
                            else:
                                label = create_label("Case sensitive is OFF")

                            Clock.schedule_once(lambda x: show_label(mp, label),
                                                .1)
                            self.mouse_over_ready = False

                    else:
                        if bt.id != self.previous_mouse_over:
                            if self.old_mouse_over:
                                self.root_window.remove_widget(
                                    self.old_mouse_over)

                            label = create_label(text=bt.id)

                            Clock.schedule_once(lambda x: show_label(mp, label),
                                                .8)
                            self.mouse_over_ready = False
            else:
                # If no collision is detected, remove any remaining label widget
                if collision is False and \
                   self.old_mouse_over in self.root_window.children:
                    self.root_window.remove_widget(self.old_mouse_over)

        # Only do this routine if the side panel is open
        if self.show_side_panel and self.mouse_over_ready:
            # Get active tab in side panel
            active_tab = self.root.ids.main_tp.current_tab.text
            # Get remove all button
            rm_bt = [self.root.ids.rm_all_file, self.root.ids.rm_all_taxa]

            # Iterate over buttons of active tab
            for bt in self.mouse_over_bts[active_tab] + sidebt_list + rm_bt:
                # Determine if there is a collision with mouse position
                if determine_collision(bt):
                    if bt in self.mouse_over_bts[active_tab]:
                        if determine_collision(self.root.ids.sv_file) or \
                                determine_collision(self.root.ids.sv_sp):
                            collision = True
                        else:
                            continue
                    else:
                        # Set collision marker to true
                        collision = True
                    # This will determine if a new label button will be added
                    # to the layout, based on the text of the button. If the
                    # text is already in the previous mouse over, then do
                    # nothing. If the text is some new button, then do something
                    if bt in self.mouse_over_bts[active_tab] + sidebt_list:
                        if bt.text != self.previous_mouse_over:
                            # Check if there is an old label button and remove
                            # it
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
                                # Create label widget
                                label = create_label(text=bt.text)

                                # Schedule the introduction of the label widget
                                Clock.schedule_once(lambda x: show_label(mp,
                                                        label), .8)
                                # Locking mouse over so that no additional label
                                # widgets are added during the waiting time
                                self.mouse_over_ready = False

                    elif bt in rm_bt:
                        if "Removes all files and taxa" \
                                not in self.previous_mouse_over:

                            label = create_label("Removes all files and taxa")
                            Clock.schedule_once(lambda x: show_label(mp, label),
                                                .3)
                            # Locking mouse over so that no additional label
                            # widgets are added during the waiting time
                            self.mouse_over_ready = False

            else:
                # If no collision is detected, remove any remaining label widget
                if collision is False and \
                   self.old_mouse_over in self.root_window.children:
                    self.root_window.remove_widget(self.old_mouse_over)

        # Only do this when plot screen is on
        if self.screen.name in self.plot_screens and self._popup not in \
                self.root_window.children:
            # Get PlotToolbar object
            toolbar_wgt = [x for x in self.root_window.children
                           if isinstance(x, PlotToolbar)][0]

            # Change toolbar opacity to become visible when collision is true
            if determine_collision(toolbar_wgt):
                if toolbar_wgt.opacity != 1:
                    Animation(opacity=1, d=.3, t="out_quart").start(toolbar_wgt)

            else:
                if toolbar_wgt.opacity == 1:
                    Animation(opacity=.2, d=.3, t="out_quart").start(
                        toolbar_wgt)

            # Check for collision with export figure or export table buttons
            if determine_collision(toolbar_wgt.ids.export_fig):
                if "mo1" not in [x.id for x in self.root_window.children]:
                    wgt, wgt2 = create_fancy_label("Export as graphics",
                                    toolbar_wgt.ids.export_fig, (100, 30))
                    self.root_window.add_widget(wgt)
                    self.root_window.add_widget(wgt2)

            else:
                for i in [x for x in self.root_window.children if
                          x.id == "Export as graphics"]:
                    self.root_window.remove_widget(i)

            if determine_collision(toolbar_wgt.ids.export_table):
                if "mo1" not in [x.id for x in self.root_window.children]:
                    wgt, wgt2 = create_fancy_label("Export as table",
                                    toolbar_wgt.ids.export_table, (100, 30))
                    self.root_window.add_widget(wgt)
                    self.root_window.add_widget(wgt2)
            else:
                for i in [x for x in self.root_window.children if
                          x.id == "Export as table"]:
                    self.root_window.remove_widget(i)

        if collision is False:
            self.previous_mouse_over = ""

    def switch_path_wgt(self, wgt_id):

        def path_updater():
            if os.path.exists(txt.text):
                self.screen.ids.icon_view_tab.path = txt.text
            else:
                return self.dialog_floatcheck("ERROR: Directory does not exist",
                                              t="error")

        label = PathLabel()
        txt = PathText()

        fc_wgt = self.screen.ids.icon_view_tab
        fc_wgt.bind(path=txt.setter("text"))
        fc_wgt.bind(path=label.setter("text"))

        self.screen.ids.path_bx.clear_widgets()

        if wgt_id == "label":
            label.text = fc_wgt.path
            self.screen.ids.path_bx.add_widget(label)
        else:
            txt.text = fc_wgt.path
            txt.bind(on_text_validate=lambda x: path_updater())
            self.screen.ids.path_bx.add_widget(txt)
            self.screen.ids.path_bx.children[0].focus = True

    def create_folder(self, text):

        path = self._popup.content.ids.sd_filechooser.path
        dir_name = join(path, text)

        if os.path.exists(dir_name):
            return self.dialog_floatcheck("The specified folder already exists",
                                          t="error")
        else:
            os.makedirs(dir_name)
            self._popup.content.ids.sd_filechooser.path = dir_name
            self.dismiss_subpopup()

    ########################## SCREEN NAVIGATION ###############################

    def go_screen(self, idx, direct="left"):
        """
        Method used to go to a specific screen by specifying and index and
        transition direction
        :param idx: integer. Index value of the screen from self.screen_names
        :param direct: string. The direction of the transition
        """

        if self.screen:
            screen_path = join(self.cur_dir, "data", "screens", "{}.kv".format(
                self.screen.name))
            self.loaded_screens[screen_path] = self.screen

            # Automatic removal of plot toolbar when not in a plot screen
            if self.screen_names[idx] not in self.plot_screens:
                self.dismiss_plot_toolbar()

        self.index = idx

        # Precludes a transition if the current screen is the same as the
        # target screen
        if self.current_screen != self.screen_names[idx]:
            # Update previous screen
            self.previous_screen = self.current_screen
            # Update current screen
            self.current_screen = self.screen_names[idx]
            # Make the switch
            self.root.ids.sm.switch_to(self.load_screen(idx), direction=direct)

    def go_previous_screen(self):
        """
        Method that returns to the previous screen, set by self.previous_screen
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

            # If the screen to be loaded is the filechooser, set the home path
            #  as the default
            if self.available_screens[idx].split("/")[-1] == "fc.kv":
                self.screen.ids.icon_view_tab.path = self.home_path
                # Initialize bookmarks
                self.bookmark_init()
                self.switch_path_wgt("label")

        return self.screen

    def go_carousel(self, slide, bt_id):
        """
        Method used by other buttons outside the side buttons of the side panel
        to go to specific slides of the side panel
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
        This method generates a desired behaviour for groups of toggle buttons
        By default, when a toggle button is pressed, the state will be down and
        a new screen/slide is presented. However, if the same toggle button is
        pressed again, it's state will return to the normal state while the
        same screen/slide is showed. To prevent this behaviour, this method
        will disable the active toggle button in the group and enable any
        other previously disabled button.

        To allow a seamless transition, ensure that background_disabled_dow is
        the same as background_down, and that disabled_color is the same as
        color.

        :param wgt: The toggle button widget. Must belong to a group.
        """

        # Iterating over the children of the parent may not be optimal, but
        # using the get_widgets(groupname) method could result in some issues
        # with garbage collector of kivy. So, for now, this will iterate over
        # all children of the toggle button's parent
        for i in [x for x in wgt.parent.children if isinstance(x,
                                                               ToggleButton)]:
            if i.disabled:
                i.disabled = False
                i.state = "normal"

        wgt.disabled = True

    def check_action(self, text, func, bt_wgt=None):
        """
        General purpose method that pops a dialog checking if the user wants to
        perform a certain action. This method should be passed as a function on
        the 'on_*' with the final function and original widget triggering the
        event as arguments
        :param func: final function if the users chooses to process
        :param bt_wgt: widget where the initial 'on_' event occurred

        Usage example:
        This can be applied to the removal button of the bookmarks. In this
        case, the event of the removal button must be like this:

        remove_button.bind(partial(self.check_action, self.remove_bookmark_bt))

        where, self.check_action is this method, and self.remove_bookmark_bt is
        the function that will actually remove the bookmark button. This
        function is then binder to the "OK" button of the check dialog.
        By default, the last argument is the bt_wgt.
        """

        check_content = CheckDialog(cancel=self.dismiss_popup)
        check_content.ids.check_text.text = text
        if bt_wgt:
            check_content.ids.check_ok.bind(on_release=lambda val: func(bt_wgt))
        else:
            check_content.ids.check_ok.bind(on_release=lambda val: func())

        self.show_popup(title="Warning!", content=check_content,
                        size=(250, 200),
                        separator_color=[255 / 255., 85 / 255., 85 / 255., 1.])

    ####################### BOOKMARKS OPERATIONS ###############################

    def bookmark_init(self):
        """
        This will create a pickle file containing a list with the bookmarks
        for the file chooser menu. If no file exists, it will create an empty
        one. If a file already exists, it will load the available bookmarks
        """

        if exists(self.bm_file):
            self.bookmarks = pickle.load(open(self.bm_file, "rb"))
            # Retrieving the bookmark path list from the self.bookmarks
            bk_list = self.bookmarks[0]
            for bk in bk_list:
                self.add_bookmark_bt(bk)

        else:
            pickle.dump(self.bookmarks, open(self.bm_file, "wb"))

    def save_bookmark(self, path):
        """
        This adds functionality to the FileChooser "Add bookmark" button. It
        will grab the selected path and add it to a storage list that
        will be saved as a pickle object and stored in a file defined in
        self.bm_file.
        :param path: String containing the path of the bookmark
        """

        # Load bookmarks object
        self.bookmarks = pickle.load(open(self.bm_file, "rb"))
        # Check if bookmark already exists. Only add bookmark if it does not
        # exist
        if path not in self.bookmarks[0]:
            # Add bookmarks to the full path list
            self.bookmarks[0].append(path)
            # Add mapping of the full path to the bookmark name
            new_map = {path.split(sep)[-1]: path}
            self.bookmarks[1] = dict(list(self.bookmarks[1].items()) +
                                     list(new_map.items()))
            self.add_bookmark_bt(path)
            pickle.dump(self.bookmarks, open(self.bm_file, "wb"))

    def add_bookmark_bt(self, bk):
        """
        This will add a bookmark button, along with its removal button. Only
        a bookmark path will be necessary.

        The path of the bookmark will be associated to the respective button
        by it's id.

        :param bk: string. bookmark file path
        """
        bookmark_name = bk.split("/")[-1]
        # Define bookmark button
        bt = Button(text=bookmark_name, id=bk,
                    height=30, size_hint=(.8, None),
                    background_normal=join("data", "backgrounds",
                                           "bt_process.png"),
                    background_down=join("data", "backgrounds",
                                         "bt_process_off.png"))
        # Bind to function that loads bookmark path into filechooser
        bt.bind(on_release=self.bookmark_load)
        # Add to list for mouse over purposes
        self.bookmarks_bt.append(bt)
        # Define bookmark removal button
        xbt = Button(size_hint=(None, None), width=30,
                     height=30, id="%sX" % bk, border=(0, 0, 0, 0),
                     background_normal=join("data", "backgrounds",
                                            "remove_bt.png"),
                     background_down=join("data", "backgrounds",
                                            "remove_bt_down.png"))
        # Bind to function that removes bookmark button as well as the path
        # from self.bm_file
        xbt.bind(on_release=partial(self.check_action,
                                    "Are you sure you want to remove"
                                    " this bookmark?",
                                    self.remove_bookmark_bt))
        # Add widgets
        self.screen.ids.sv_book.add_widget(bt)
        self.screen.ids.sv_book.add_widget(xbt)

    def bookmark_load(self, value):
        """
        Provided a bookmark button object, it loads the bookmark file path
        that is stored in the button id.
        :param value: bookmark button object
        """

        path = value.id
        self.screen.ids.icon_view_tab.path = path

    def remove_bookmark_bt(self, value):
        """
        Adds functionality to the removal button associated with each bookmark
        button. This will not only remove the
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
        bk_name = bk_idx.split("/")[-1]
        # Remove bookmark path from list and mapping dictionary
        self.bookmarks[0].remove(bk_idx)
        del self.bookmarks[1][bk_name]
        # Update self.bm_file
        pickle.dump(self.bookmarks, open(self.bm_file, "wb"))

    ######################## PLOT SCREEN OPERATIONS ############################

    def dialog_export_graphic(self):
        """
        Creates a filechooser dialog for graphics exportation. It differs from
        other filechooser dialogs in the presence of a spinner to select
        the graphical extension
        :return:
        """

        content = ExportGraphics(cancel=self.dismiss_popup)

        self.show_popup(title="Export graphic as...", content=content,
                        size_hint=(.9, .9))

    def export_graphic(self, path, file_name, ext):
        """
        Saves the current plot object into a file based on file name and
        extension
        :param path: string, path to final directory
        :param file_name: string, name of graphic file
        :param ext: string, extension of graphic file (e.g. png, svg, pdf, etc)
        """

        self.current_plot.savefig(join(path, file_name + ext),
                                  bbox_extra_artists=(self.current_lgd,),
                                  bbox_inches="tight")

    ######################## SIDE PANEL OPERATIONS #############################

    def sidepanel_on_touch(self, *args):
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
        bar, which is also used to toggle the side panel. This prevents issues
        of toggling the side panel twice with one mouse input
        - When a popup is not open. There are several buttons in the side bar
        that open whose position is outside the side bar. The user should be
        able to click anywhere in the popup without the side panel closing.

        In addition, it will handle the status of the partition box dialog,
        associated with the sidepanel. While this box is active, the sidepanel
        must remain open. However, clicks outside the partition box will close
        it.
        """

        def animate_sidebar():

            ## ANIMATIONS with hierarchy
            # Animation of main BoxLayout containing child ScrollViews
            self.sidepanel_animation(width=0,
                                      wgt=self.root.ids.main_box)
            # Animation of both scrollviews
            self.sidepanel_animation(width=0,
                                      wgt=self.root.ids.sp)
            self.sidepanel_animation(width=0,
                                      wgt=self.root.ids.sp_bts)

            self.show_side_panel = not self.show_side_panel

        # Get mouse position
        mp = self.root_window.mouse_pos
        # Get side panel and previous button widgets
        side_panel_wgt = self.root.ids.main_box
        ap = self.root.ids.ap

        # Check for existence of a partition dialog box
        partition_box = [x for x in self.root_window.children if
                         isinstance(x, PartitionsDialog)]

        # If the partition box exists and the collision is outside it
        if partition_box and not partition_box[0].collide_point(
                mp[0], mp[1]):
            # Check if spinner is open
            spin1 = partition_box[0].ids.codon_spin.is_open
            spin2 = [x.is_open for x in partition_box[0].ids.model_bx.children]

            # If the spinners are not open, remove
            if True not in spin2 and not spin1:
                rm_bt = [x for x in self.root_window.children if
                        isinstance(x, RemoveFloat)][0]
                self.root_window.remove_widget(partition_box[0])
                self.root_window.remove_widget(rm_bt)

        # Check for conditions to close the side panel.
        # If touch is out of panel; if panel is open; is touch is out of menu
        # button; a popup is not open
        if side_panel_wgt.collide_point(mp[0], mp[1]) is False\
                and self.show_side_panel \
                and ap.collide_point(mp[0], mp[1]) is False \
                and self._popup not in self.root_window.children \
                and not partition_box:

            if self.screen.name == "Process":
                queue_bt = self.screen.ids.queue_bt
                if queue_bt.collide_point(mp[0], mp[1]) is False:
                    animate_sidebar()
            else:
                animate_sidebar()

    def toggle_sidepanel(self):
        """
        Method controlling the animation toggling of the side panel
        """

        # Toggling the state of the panel. This attribute is the main
        # controller of the side panel state. When its True, the side panel is
        # extended, otherwise the side panel is hidden
        self.show_side_panel = not self.show_side_panel

        if self.show_side_panel:

            # Redraw the side panel layout. This will ensure that the widget
            # is always on top of all widgets.
            self.root.ids.bx1.remove_widget(self.root.ids.panel_float)
            self.root.ids.bx1.add_widget(self.root.ids.panel_float)

            # Fixing the width of the side panel
            # Main panel width
            sv_panel_width = 300
            # Side buttons width
            sv_bts_width = 60
        else:
            sv_panel_width, sv_bts_width = 0, 0

        ## ANIMATIONS with hierarchy
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
        elif panel == "taxa":
            bt_list = sorted(self.active_taxa_list)
            gl_wgt = self.root.ids.taxa_sl
        else:
            bt_list = list(self.alignment_list.partitions.partitions.keys())
            gl_wgt = self.root.ids.partition_sl

        # Find buttons that match the txt string
        if panel == "files":
            found_bts = [el.split(sep)[-1] for el in bt_list if
                         txt.lower() in el.split(sep)[-1].lower()]
        else:
            found_bts = [el for el in bt_list if txt.lower() in el.lower()]

        # Clear the grid and populate with the found bts
        gl_wgt.clear_widgets()
        mouse_bts = []
        for txt in found_bts:
            bt, inf_bt, x_bt = self.sidepanel_create_bts(txt)
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
        Clears previous search string and populates with the original buttons
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

        try:
            d = self.file_list[self.count_files + 1]
            gl_wgt.add_widget(LoadMoreBt())
        except IndexError:
            return

    def load(self, selection, bad_aln):
        """
        Loads selected input files into the program. This should be switched
        after selecting files in a FileChooser widget. The selected files
        will be parsed and the side panel will be populated with information
        on file names and taxa names.

        To allow different files to be loaded on independent occasions,
        before performing the loading operations, the self.file_list
        attribute that contains the full paths to the input files will be
        checked to see if its empty or populated. When empty, attributes will
        be loaded anew; otherwise, the existing populated attributes will be
        extended.

        :param selection: list. Contains the paths to the selected input files
        """

        # Parses the files into the program
        #bad_aln = self.load_files(selection)

        # Checking for input sequence type inconsistencies. If there are
        # alignments with different sequence types, then issue and error popup
        # and do not load the files
        if isinstance(bad_aln, Exception):
            self.dialog_warning("Multiple sequence types detected", "The "
                                "selected input alignments contain more than "
                                "one sequence type (DNA, RNA, Protein). Please"
                                " select input files of the same sequence type")
            return 0
        else:
            # Checking if there are invalid input alignments
            if bad_aln:
                if len(bad_aln) == 1:
                    self.dialog_warning("Invalid input file detected",
                                        "The following input file is invalid:"
                                        "\n\n[b]%s[/b]" %
                                        bad_aln[0].name.split(sep)[-1])
                else:
                    self.dialog_warning("Invalid input files detected",
                                        "The following input files are invalid:"
                                        "\n\n[b]%s[/b]" %
                                        "\n".join(x.name.split(sep)[-1] for x
                                                  in bad_aln))

            # removes bad alignment files from selection list
            selection = [path for path in selection if
                         path.split(sep)[-1].split(".")[0] not in
                         [x.name for x in bad_aln]]

            # If data has been previously loaded, updated these attributes
            if self.file_list:
                # Updating complete and active file lists
                self.file_list.extend(selection)
                self.active_file_list.extend(selection)
                # Update the filename - path mapping attribute
                self.filename_map = dict(list(self.filename_map.items()) +
                                     list((x, y) for x, y in
                                     zip([x.split("/")[-1] for x in selection],
                                         selection)))

            # If no data has been previously loaded, set the attributed
            else:
                # Set an attribute with the input file list
                self.file_list = selection
                # Setting active file list and path list
                self.active_file_list = deepcopy(self.file_list)
                # Sett the filename - path mapping attribute
                self.filename_map = dict((x, y) for x, y in zip(
                    [x.split("/")[-1] for x in selection], selection))

            # If more than one alignment has been provided
            if self.active_alignment_list:
                # Update active taxa list
                self.update_taxa()
                # Populates files and taxa contents
                self.update_tabs()
                # Gathers taxa  and file information
                self.original_tx_inf = self.get_taxa_information()
                #self.original_file_inf = self.get_file_information()

                # Issue float check
                self.dialog_floatcheck("%s file(s) successfully loaded" %
                                       len(selection), t="info")

    def load_proteomes(self, selection):
        """
        Similar to load method, but specific for loading proteome files. Given
        the potential size of these files, they are not stored in memory, but
        instead are processed on the fly
        """

        # Stores invalid proteome files
        bad_proteomes = {"invalid": [], "no_fasta": [], "no_protein": []}
        good_proteomes = []

        # Check input proteomes
        for f in selection:
            b = Base()
            er = b.autofinder(f)
            f_short = f.split(sep)[-1]
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

            # Update the filename - path mapping attribute
            self.filename_map = dict(list(self.filename_map.items()) +
                list((x, y) for x, y in zip([x.split(sep)[-1] for x in
                good_proteomes], good_proteomes)))

            # Populate file buttons in side panel
            self.original_file_inf = self.get_file_information(mode="proteome")
            self.populate_input_files(mode="proteome")

        if list(bad_proteomes.values()) != [[], [], []]:
            msg = ""
            if bad_proteomes["invalid"]:
                msg += "The following files are in invalid format:\n%s\n\n" % \
                        ", ".join(bad_proteomes["invalid"])
            if bad_proteomes["no_fasta"]:
                msg += "The following files are not in FASTA format:\n%s\n\n" \
                        % ", ".join(bad_proteomes["no_fasta"])
            if bad_proteomes["no_protein"]:
                msg += "The following files do not contain protein sequences:" \
                       "\n%s\n\n" % ", ".join(bad_proteomes["no_protein"])

            return self.dialog_warning("Invalid proteome files detected", msg)

    def update_tabs(self):
        """
        Wrapper that updates the contents of the files and taxa tabs
        """

        self.populate_input_files()
        self.populate_species()
        self.populate_partitions()

    def update_taxa(self):
        """
        This checks whether some taxa that were specific to some file(s) were
        removed when that file is removed.
        """

        # If taxa were removed during the update, remove those buttons too
        removed_taxa = list(set(self.active_taxa_list) - set(
            self.active_alignment_list.taxa_names))
        if removed_taxa:
            for i in removed_taxa:
                # Get the corresponding buttons:
                x_but_txt = "%sX" % i
                bt_obj = [x for x in self.root.ids.taxa_sl.children if x_but_txt
                          == x.id][0]
                self.remove_bt(bt_obj)

        self.active_taxa_list = self.active_alignment_list.taxa_names

    def update_partitions(self):
        """
        Updates partition buttons following any change to input data
        """

        # Check for missing partitions based on id and remove them
        for bt, inf_bt, x_bt in self.sp_partition_bts:
            if bt.id not in self.alignment_list.partitions.partitions:
                self.remove_bt(x_bt)

    def update_file_label(self, mode="alignment"):
        """
        Sets and updates a label on the Files tab of the side panel, informing
        how many files are selected out of the total files

        :param mode: string. Determines which list will be used to populate.
        If alignment, it will use file_list; If proteome, it will use
        proteome_list
        """

        # Determine which list is used to populate
        if mode == "alignment":
            lst = self.file_list
            active_lst = self.active_file_list
        else:
            lst = self.proteome_files
            active_lst = self.active_proteome_files

        self.root.ids.file_lab.text = "%s of %s files selected" % (
                                       len(active_lst),
                                       len(lst))

        # Reset the sequence_types list when all files are removed
        if len(active_lst) == 0 and len(lst) == 0:
            # Core changes
            self.sequence_types = []

            # Add disabled no changes button
            if "species_temp" not in [x.id for x in
                                      self.root.ids.file_sl.children]:
                no_bt = Button(id="file_temp", text="No files loaded",
                               size_hint_y=None, height=40, disabled=True)
                self.root.ids["file_temp"] = no_bt
                self.root.ids.file_sl.add_widget(no_bt)

    def update_sp_label(self):
        """
        Sets and updates a label on the Taxa tab of the side panel, informing
        how many taxa are selected out of the total taxa. If the taxa list
        is empty, it disables the select/deselect buttons
        """

        self.root.ids.sp_lab.text = "%s of %s taxa selected" % (
                                       len(self.active_taxa_list),
                                       len(self.alignment_list.taxa_names))

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

        bt = ToggleButton(text=idx, state="down", id=idx,
                          height=30, size_hint=(.8, None), shorten=True,
                          shorten_from="right", halign="center",
                          bold=True,
                          background_down=join("data", "backgrounds",
                                                 "bt_process.png"),
                          background_normal=join("data", "backgrounds",
                                                 "bt_process_off.png"))

        # Setting horizontal text size for shortening
        bt.text_size[0] = bt.size[0] * 1.3
        # Binding functionality to toggle button
        bt.bind(on_release=self.toggle_selection)

        # Set Information button and add the widget
        inf_bt = Button(size_hint=(None, None), width=30,
                        height=30, id="%s?" % idx,
                        background_normal=join("data", "backgrounds",
                                               "info_bt.png"),
                        background_down=join("data", "backgrounds",
                                               "info_bt_down.png"))
        inf_bt.bind(on_release=self.popup_info)

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

        # Add a label at the end of the file list informing how many files are
        # currently selected out of the total files
        self.update_file_label(mode=mode)

        for infile in lst:

            if self.count_files <= self.MAX_FILE_BUTTON:

                self.count_files += 1
                file_name = infile.split("/")[-1]
                self.sidepanel_add_bts(file_name, "Files")

            else:
                self.root.ids.file_sl.add_widget(LoadMoreBt())
                return

    def sidepanel_add_bts(self, idx, tab_name):

        # Set attributes to be added
        if tab_name == "Files":
            grid_wgt = self.root.ids.file_sl
            bt_list = self.sp_file_bts

        elif tab_name == "Taxa":
            grid_wgt = self.root.ids.taxa_sl
            bt_list = self.sp_taxa_bts

        elif tab_name == "Partitions":
            grid_wgt = self.root.ids.partition_sl
            bt_list = self.sp_partition_bts

        # This prevents duplicate entries from being added
        if idx not in [x.id for x in grid_wgt.children]:

            # Create buttons
            if tab_name == "Partitions":
                bt, inf_bt, x_bt = self.sidepanel_create_part_bts(idx)
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

        MAX_BUTTONS = self.MAX_FILE_BUTTON + self.count_files

        self.root.ids.file_sl.remove_widget(self.root.ids.file_sl.children[0])

        for i in range(self.count_files, MAX_BUTTONS):

            self.count_files += 1

            try:
                infile = self.file_list[self.count_files]
                self.sidepanel_add_filebts(infile)
            except IndexError:
                return

        else:
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

        # Add a label at the end of the taxa list informing how many taxa are
        # currently selected out of the total taxa
        self.update_sp_label()

        for tx in sorted(self.active_taxa_list):

            # Prevents duplicate taxa from being entered
            if tx not in [x.id for x in self.root.ids.taxa_sl.children]:

                self.sidepanel_add_bts(tx, "Taxa")

    def populate_partitions(self):
        """
        Populates the partitions tab in the side bar from the partitions object
        associated with alignment objects.

        This method is used when input files are loaded into the program,
        which means there will be no issue with multiple files being associated
        with the same partitions. This kind of change is done a posteriori
        when importing partition files or setting the partitions manually.
        """

        # Remove initial disabled button, if it's still there
        if "partition_temp" in self.root.ids.keys():
            self.root.ids.partition_sl.remove_widget(
                self.root.ids.partition_temp)
            del self.root.ids["partition_temp"]

        for aln in self.alignment_list:
            for partition, vals in aln.partitions:

                if self.count_partitions <= self.MAX_PARTITION_BUTTON:

                    self.count_partitions += 1
                    # Create partition buttons
                    self.sidepanel_add_bts(partition, "Partitions")

                else:
                    self.root.ids.partition_sl.add_widget(LoadMoreBt())
                    return

    def sidepanel_create_part_bts(self, idx):
        """
        Creates buttons for each partition
        :param idx: string. unique identifier of partition
        """

        # Create main button
        bt = ToggleButton(text=idx, state="normal", id=idx,
                          height=30, size_hint=(.8, None), shorten=True,
                          shorten_from="right", halign="center",
                          bold=True,
                          background_down=join("data", "backgrounds",
                                                 "bt_process.png"),
                          background_normal=join("data", "backgrounds",
                                                 "bt_process_off.png"))

        # Setting horizontal text size for shortening
        bt.text_size[0] = bt.size[0] * 1.3

        # Create edition button
        ed_bt = Button(size_hint=(None, None), width=30,
                        height=30, id="%s?" % idx,
                        background_normal=join("data", "backgrounds",
                                               "edit_bt.png"),
                        background_down=join("data", "backgrounds",
                                               "edit_bt_down.png"))
        ed_bt.bind(on_release=self.dialog_partitions)

        # Create removal button
        x_bt = Button(size_hint=(None, None), width=30,
                      height=30, id="%sX" % idx,
                      border=(0, 0, 0, 0),
                      background_normal=join("data", "backgrounds",
                                             "remove_bt.png"),
                      background_down=join("data", "backgrounds",
                                             "remove_bt_down.png"))

        return bt, ed_bt, x_bt

    def set_codon_model(self, codon_partition, wgt=None):
        """
        Changes the model spinners when changing the codon partitioning
        """

        first_background = join("data", "backgrounds", "model_bt1.png")
        second_background = join("data", "backgrounds", "model_bt2.png")
        third_background = join("data", "backgrounds", "model_bt3.png")

        partition_model = {"[color=ff5555ff]1[/color] + [color=37abc8ff]2"
                           "[/color] + [color=71c837ff]3[/color]":
                [ModelSpinner(background_normal=first_background, id="1"),
                ModelSpinner(background_normal=second_background, id="2"),
                ModelSpinner(background_normal=third_background, id="3")],
                           "[color=ff5555ff](1 + 2)[/color] + [color=37abc8ff]"
                           "3[/color]":
                [ModelSpinner(background_normal=first_background, id="12"),
                 ModelSpinner(background_normal=second_background, id="3")],
                           "[color=ff5555ff]1[/color] + [color=37abc8ff](2 + "
                           "3)[/color]":
                [ModelSpinner(background_normal=first_background, id="1"),
                 ModelSpinner(background_normal=second_background, id="23")],
                           "[color=ff5555ff](1 + 3)[/color] + [color=37abc8ff]"
                           "2[/color]":
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

    def dialog_partitions(self, btx):
        """
        Shows a small widget with partition information
        """

        # Get position of partition edit button:
        ed_pos = btx.to_window(btx.pos[0], btx.pos[1])

        # Set position for partitions dialog
        size = (180, 190)
        pos = [ed_pos[0] + btx.width,
               ed_pos[1] + (btx.height / 2) - (size[1] / 2)]

        content = PartitionsDialog(pos=pos, size=size, size_hint=(None, None))
        rm_wgt = RemoveFloat(pos=[pos[0] + size[0] - 20, pos[1] + size[1] - 20])

        # Set partition object and partition name
        part_obj = self.alignment_list.partitions
        part_name = btx.id[:-1]

        #TODO: For now this assumes all codon partitions are unlinked
        # If there are codon partitions
        if part_obj.partitions[part_name][1]:
            content.ids.codon_spin.text = content.ids.codon_spin.values[1]
            self.set_codon_model(content.ids.codon_spin.values[1], content)
            for i in range(len(part_obj.models[part_name][0])):
                params = part_obj.models[part_name][0][i]
                model = part_obj.get_model_name(params)
                content.ids.model_bx.children[i].text = model

        # Give functionality to remove button
        rm_wgt.bind(on_release=lambda x: self.remove_partition_box())

        self.root_window.add_widget(content)
        self.root_window.add_widget(rm_wgt)

    def popup_info(self, value):
        """
        Generates the pop up information content for the pressed taxa or file
        button
        :param value: the button object is provided when binding
        """

        # Determining if the request comes from file or taxa tab
        if value.parent == self.root.ids.taxa_sl:

            # Get taxa name
            tx = value.id[:-1]

            if tx in self.active_taxa_list:

                # Get the information from the content list. This is done when
                # calling the popup to avoid repeating this operation every time
                # taxa  or files are added/removed.
                self.active_tx_inf = self.get_taxa_information()

                content = BoxLayout(orientation="vertical", padding=10,
                                    spacing=10)
                sv = ScrollView(scroll_type=["bars"], bar_width=10)
                all_ds = BoxLayout(orientation="vertical",
                                   height=2 * (30 * 7) + 10, size_hint_y=None)
                total_ds = TaxaPopup(height=30 * 7)
                total_ds.ids.dataset_label.text = "Complete data set"
                active_ds = TaxaPopup(height=30 * 7)
                active_ds.ids.dataset_label.text = "Active data set"

                # Populate complete data set contents
                total_ds.ids.seq_len.text = "%s" % \
                                            self.original_tx_inf[tx]["length"]
                total_ds.ids.indels.text = "%s" % \
                                           self.original_tx_inf[tx]["indel"]
                total_ds.ids.missing.text = "%s" %\
                                            self.original_tx_inf[tx]["missing"]
                total_ds.ids.ef_seq_len.text = ("%s (%s%%)" % (
                    self.original_tx_inf[tx]["effective_len"],
                    self.original_tx_inf[tx]["effective_len_per"]))
                total_ds.ids.file_cov.text = ("%s (%s%%)" % (
                    self.original_tx_inf[tx]["fl_coverage"],
                    self.original_tx_inf[tx]["fl_coverage_per"]))

                # Populate active data set contents
                active_ds.ids.seq_len.text = "%s" % \
                                             self.active_tx_inf[tx]["length"]
                active_ds.ids.indels.text = "%s" %\
                                            self.active_tx_inf[tx]["indel"]
                active_ds.ids.missing.text = "%s" %\
                                             self.active_tx_inf[tx]["missing"]
                active_ds.ids.ef_seq_len.text = ("%s (%s%%)" % (
                    self.active_tx_inf[tx]["effective_len"],
                    self.active_tx_inf[tx]["effective_len_per"]))
                active_ds.ids.file_cov.text = ("%s (%s%%)" % (
                    self.active_tx_inf[tx]["fl_coverage"],
                    self.active_tx_inf[tx]["fl_coverage_per"]))

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

                # Get the information from the content list. This is done when
                # calling the popup to avoid repeating this operation every time
                # taxa  or files are added/removed.
                self.active_file_inf = self.get_file_information()

                content.ids.in_format.text = "%s" % \
                                self.active_file_inf[file_name]["aln_format"]
                content.ids.seq_type.text = "%s" % \
                                self.active_file_inf[file_name]["seq_type"]
                content.ids.is_aln.text = "%s" % \
                                self.active_file_inf[file_name]["is_aln"]
                content.ids.seq_size.text = "%s" % \
                                self.active_file_inf[file_name]["aln_len"]
                content.ids.n_taxa.text = "%s" % \
                                self.active_file_inf[file_name]["n_taxa"]

                self.show_popup(title="File: %s" % value.id[:-1],
                                content=content, size=(400, 320))

            if self.filename_map[file_name] in self.active_proteome_files:

                content = ProteomePopup(cancel=self.dismiss_popup)

                content.ids.n_seq.text = "%s" % \
                                self.original_file_inf[file_name]["n_seq"]
                content.ids.n_res.text = "%s" % \
                                self.original_file_inf[file_name]["n_res"]

                self.show_popup(title="File %s" % value.id[:-1],
                                content=content, size=(400, 200))

    def export_names(self, path, file_name):
        """
        Export the names of buttons in the corresponding tab in the side panel
        It listens to the self.export_mode attribute, which is a tuple object
        with the first element being either "file" or "taxa" and the second
        element as "all" or "selected".

        :param path: string. Path to the output file.
        :param file_name. Name of the output file.
        """

        # Create file object
        export_file = open(join(path, file_name) + ".txt", "w")

        if self.export_mode[0] == "file":
            # Export all files
            if self.export_mode[1] == "all":
                for x in self.file_list:
                    short_name = x.split(sep)[-1]
                    export_file.write(short_name + "\n")
            # Export selected files
            elif self.export_mode[1] == "selected":
                for x in self.active_file_list:
                    short_name = x.split(sep)[-1]
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

    def toggle_selection(self, value):
        """
        Adds functionality for the file and taxa toggle buttons in the side
        panel. It adds or removes the selected taxa from the active lists
        """

        # Get the parent layout object
        parent_obj = value.parent

        # Changes concerning the files tab
        if parent_obj == self.root.ids.file_sl:

            # When button is normal (unselected) remove from active list
            if value.state == "normal":
                self.active_file_list.remove(self.filename_map[value.id])
                self.active_alignment_list.remove_file(
                    [self.filename_map[value.id]])
            # When button is down (selected) add to active list
            elif value.state == "down":
                self.active_file_list.append(self.filename_map[value.id])
                self.active_alignment_list.add_alignment(
                    self.alignment_list.retrieve_alignment(
                        self.filename_map[value.id]))

            # Update label
            self.update_file_label()

        # Changes concerning the taxa tab
        if parent_obj == self.root.ids.taxa_sl:

            # When button is normal (unselected) remove from active list
            if value.state == "normal":
                self.active_taxa_list.remove(value.id)
            # When button is down (selected) add to active
            elif value.state == "down":
                self.active_taxa_list.append(value.id)

            # Update label
            self.update_sp_label()

    def remove_all(self):
        """
        Functionality for the remove all button for taxa and file buttons in the
        side panel. This method will remove all files and taxa from the program
        """

        # Get file remove button list
        file_bts = [x for x in self.root.ids.file_sl.children if
                    x.background_normal == join("data", "backgrounds",
                                                "remove_bt.png")]

        for i in file_bts:
            self.remove_bt(i)

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
        self.ortho_groups = ObjectProperty(None)

    def remove_groups(self, value):
        """
        Removes orthology group buttons
        :param value: Instance of remove button
        """

        # Remove group from MultiGroup object
        self.ortho_groups.remove_group(value.id)

        # Get box container of all gridlayouts
        gl_bx = value.parent.parent

        for gl in gl_bx.children:
            # Remove appropriate item, according to id, from its gridlayout
            gl.remove_widget([x for x in gl.children if x.id == value.id][0])

        # If no group button is active, dispatch the first
        if not [x for x in self.screen.ids.group_gl.children
                if x.state == "down"] and self.screen.ids.group_gl.children:
            self.screen.ids.group_gl.children[-1].dispatch("on_release")
            self.screen.ids.group_gl.children[-1].state = "down"

        if not self.screen.ids.group_gl.children:
            self.screen.ids.card_gl.clear_widgets()

    def remove_bt(self, value):
        """
        Functionality for the "X" remove buttons in the side panel. It
        removes button pairs with similar id's and can be used in both files
        and taxa tabs
        """

        ####### APP CHANGES
        # Get the parent layout object from where the widget will be removed
        parent_obj = value.parent

        # Get button widgets to be removed
        bt_idx = value.id[:-1]
        inf_idx = value.id[:-1] + "?"
        bt = [x for x in parent_obj.children if bt_idx == x.id][0]
        inf_bt = [x for x in parent_obj.children if inf_idx == x.id][0]

        # Remove widgets
        parent_obj.remove_widget(value)
        parent_obj.remove_widget(bt)
        parent_obj.remove_widget(inf_bt)

        ####### CORE CHANGES
        # Get the parent tab
        if parent_obj == self.root.ids.file_sl:
            # Update file list
            file_path = self.filename_map[bt_idx]
            self.file_list.remove(file_path)
            # Update active file list. If the file has been removed from the
            # active list, this will handle the exception
            try:
                self.active_file_list.remove(file_path)
            except ValueError:
                pass
            # Update alignment object list
            self.alignment_list.remove_file([self.filename_map[bt_idx]])
            self.active_alignment_list.remove_file([self.filename_map[bt_idx]])

            # Update active taxa list. This must be executed before calling
            # self.get_taxa_information since this method relies on an
            # updated active taxa list
            self.update_taxa()

            # Updates the partition list
            self.update_partitions()

            # Update pop up content. Since the file has been removed,
            # it should also be excluded from the complete data set
            self.original_tx_inf = self.get_taxa_information(
                alt_list=self.alignment_list)

            # Updates labels
            self.update_file_label()
            self.update_sp_label()

        if parent_obj == self.root.ids.taxa_sl:
            self.active_alignment_list.remove_taxa([bt_idx])
            self.alignment_list.remove_taxa([bt_idx])
            self.active_taxa_list = self.active_alignment_list.taxa_names
            # Updates label
            self.update_sp_label()

    def select_bt(self, value):
        """
        Functionality to the Select All/Deselect All buttons of the side
        panel. The method was made in such a way that it could be of general
        use for buttons in the files and taxa tabs
        """

        sv_parent = [x for x in value.parent.parent.parent.children if
                     "scrollview" in str(x.__class__)][0]

        # This will iterate over the first child of the parent scrollview.
        # Since scroll view only supports one child, this should be fine
        for i in sv_parent.children[0].children:
            # Skips the X buttons
            if "togglebutton" in str(i.__class__):

                if value.text == "Select All":
                    # App related action
                    i.state = "down"

                elif value.text == "Deselect All":
                    # App related action
                    i.state = "normal"

        # Core changes to files
        if sv_parent == self.root.ids.sv_file and value.text == "Select All":
            self.active_file_list = self.file_list[:]
            self.active_alignment_list = deepcopy(
                self.alignment_list)
            # Update label

        #Core changes to taxa
        if sv_parent == self.root.ids.sv_sp and value.text == "Select All":
            self.active_taxa_list = deepcopy(
                self.alignment_list.taxa_names)

        # Core changes to files
        if sv_parent == self.root.ids.sv_file and value.text == "Deselect All":
            self.active_file_list = []
            self.active_alignment_list.clear_files()
        # Core changes to taxa
        if sv_parent == self.root.ids.sv_sp and value.text == "Deselect All":
            self.active_taxa_list = []

        # Updates labels
        self.update_sp_label()
        self.update_file_label()

    def dialog_taxagroup(self, ds_type):
        """
        Creates the layout for the taxa group creation popup.
        :param ds_type: string. Data set type. It may be either "taxa" or
        "files"
        """

        # Initializing instance for taxa group dialog
        content = TaxaGroupDialog(cancel=self.dismiss_popup)

        if ds_type == "taxa":
            bt_list = sorted(self.alignment_list.taxa_names, reverse=True)
            title = "Create taxa groups"
        elif ds_type == "files":
            bt_list = sorted([x.split(sep)[-1] for x in self.file_list],
                             reverse=True)
            title = "Create file groups"

        # Populate the gridlayout for all entries
        for i in bt_list:
            # Create togglebutton for each entry
            bt = ToggleButton(text=i, size_hint_y=None, height=30)
            self.add_dataset_bt(bt, content.ids.all_grid, ds_type)

        content.ds_type = ds_type

        # Show dialog
        self.show_popup(title=title, content=content,
                        size=(700, 500))

    def add_dataset_bt(self, bt, wgt, ds_type):
        """
        Method for addition of a button to a widget. This method was created
        for the automatic upated of the widgets height when moving buttons in
        the taxa group creation dialog
        :param bt: The button widget
        :param wgt: The sink widget
        :param ds_type: string. Data set type. It may be either "taxa" or
        "files"
        """

        if ds_type == "taxa":
            bt_list = sorted(self.alignment_list.taxa_names, reverse=True)
        elif ds_type == "files":
            bt_list = sorted([x.split(sep)[-1] for x in self.file_list],
                             reverse=True)

        wgt.add_widget(bt, bt_list.index(bt.text))
        wgt.height += 30

    @staticmethod
    def remove_taxa_bt(bt, wgt):
        """
        Method for addition of a button to a widget. This method was created
        for the automatic upated of the widgets height when moving buttons in
        the taxa group creation dialog
        :param bt: The button widget
        :param wgt: The source widget
        """
        wgt.remove_widget(bt)
        wgt.height -= 30

    def taxagroup_move_taxa(self, source_wgt, sink_wgt, all_taxa, ds_type):
        """
        Method that adds functionality to the addition/removal buttons (<<, <,
        >>, >) in the taxa group dialog.
        :param source_wgt: widget, the gridlayout from where the buttons will
        be moved
        :param sink_wgt: widget, the gridlayout to where buttons will be moved
        :param all_taxa: Boolean, if True its as if alsa taxa were selected to
        be moved
        :param ds_type: string. Data set type. It may be either "taxa" or
        "files"
        """

        if ds_type == "taxa":
            bt_list = sorted(self.alignment_list.taxa_names, reverse=True)
        elif ds_type == "files":
            bt_list = sorted([x.split(sep)[-1] for x in self.file_list],
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
            # This workaround is used to add some buttons from the source to the
            # sink widgets while maintaining their original order. The z-index
            # of widgets is not working quite as I expected, so for now this
            # will produce the desired behaviour
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
                self.add_dataset_bt(bt, sink_wgt, ds_type)

    def taxagroups_show_taxa(self, name_wgt):
        """
        Creates a popup listing the taxa included in a taxa group given by name
        :param name_wgt: widget, widget containing the name of the group as text
        """

        # Create root boxlayout
        content = BoxLayout(orientation="vertical", padding=10, spacing=10)
        # Create scroll view in which the gridlayout will be inserted
        sv = ScrollView()
        # Create close button for the popup
        close_bt = Button(text="Close", size_hint_y=None, height=30,
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
        if parent_wgt == self.root.ids.file_group_grid:
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
        self.show_popup(title="Taxa group: %s" % name_wgt.text, content=content,
                        size_hint=(.3, .7))

    def remove_taxa_group(self, rm_wgt):
        """
        Removes the data set group button from the app list and corresponding
        data set group attribute
        :param rm_wgt: widget, widget of the removal button
        :return:
        """

        # Remove from app
        parent_wgt = rm_wgt.parent

        bt_idx = rm_wgt.id[:-1]
        bt = [x for x in parent_wgt.children if bt_idx == x.id][0]

        parent_wgt.remove_widget(bt)
        parent_wgt.remove_widget(rm_wgt)

        # Remove from program attribute
        if parent_wgt == self.root.ids.taxa_group_grid:
            del self.taxa_groups[bt_idx]
        if parent_wgt == self.root.ids.file_group_grid:
            del self.file_groups[bt_idx]

    def save_dataset_group(self, source_wgt, name, ds_type):
        """
        Adds a taxa group declared using the taxa group creator popup to the
        list of taxa groups in the side panel
        :param source_wgt, gridlayout of the selected items
        :param name: string, name of the group
        :param ds_type: string. Data set type. It may be either "taxa" or
        "files"
        """

        if ds_type == "taxa":
            # Make core changes by populating self.taxa_groups dictionary
            self.taxa_groups[name] = []
            group_list = self.taxa_groups[name]
            # Set the grid layout where the group button is to be added
            grid_layout = self.root.ids.taxa_group_grid
            # Set dropdown widget
            dd_wgt = self.process_grid_wgt.ids.taxa_dropdown
        elif ds_type == "files":
            # Make core changes by populating self.file_groups dictionary
            self.file_groups[name] = []
            group_list = self.file_groups[name]
            # Set the grid layout where the group button is to be added
            grid_layout = self.root.ids.file_group_grid
            # Set dropdown widget
            dd_wgt = self.process_grid_wgt.ids.file_dropdown

        for bt in source_wgt.children:
            group_list.append(bt.text)

        # App changes by adding two buttons for the taxa group
        # Taxa button itself
        bt = Button(text=name, size_hint=(.8, None), height=30, id=name,
                    background_normal="data/backgrounds/bt_process.png",
                    background_down="data/backgrounds/bt_process_off.png",
                    bold=True)
        bt.bind(on_release=self.taxagroups_show_taxa)
        # Removal button
        x_bt = Button(size_hint=(None, None), width=30, border=(0, 0, 0, 0),
                        height=30, id="%sX" % name,
                        background_normal=join("data", "backgrounds",
                                               "remove_bt.png"))
        x_bt.bind(on_release=partial(self.check_action,
                                     "Are you sure you want to remove this"
                                     " group?",
                                     self.remove_taxa_group))

        # Add buttons to gridlayout
        for i in [bt, x_bt]:
            grid_layout.add_widget(i)

        # Create separator between dropdown items
        separator = Widget(size_hint_y=None, height=3)
        dd_bt = Button(text=name, size_hint_y=None, height=40,
                       background_normal="data/backgrounds/bt_process.png",
                       background_color=(1, 1, 1, .3), bold=True)
        dd_bt.bind(on_release=lambda x:
                   dd_wgt.select(name))
        dd_wgt.add_widget(separator)
        dd_wgt.add_widget(dd_bt)

        # Update gridlayout height
        grid_layout.height += 40

    def dialog_general_info(self, idx):
        """
        Generates the popup with information for several components of the
        application
        :param idx: string. Identifier of the informative content to be shown.
        It must be present in the dictionary keys of the informative_storage
        variable in data/resources/info_data.py
        """

        content = InfoPopup(cancel=self.dismiss_popup)

        # Retrieve title and body text
        title_str, body_str = informative_storage[idx]

        # Add text body
        content.ids.content.text = body_str

        self.show_popup(title=title_str, content=content, size=(400, 300))

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
            TreeViewLabel(text="Orthology Operations", bold=True, font_size=20,
                          color=(1, 0.3, 0.3, .2)))
        proc_node = self.operation_tv.add_node(
            TreeViewLabel(text="Process Operations", bold=True, font_size=20,
                          color=(.3, .3, 1, 1)))
        stat_node = self.operation_tv.add_node(
            TreeViewLabel(text="Statistics Operations", bold=True, font_size=20,
                          color=(.3, 1, .3, .2)))

        # Main subnodes for Process
        main_op_node = self.operation_tv.add_node(TreeViewLabel(
            text="Main Operation", bold=True, font_size=15, opacity=.2),
            proc_node)
        secondary_op_node = self.operation_tv.add_node(TreeViewLabel(
            text="Secondary Operations", bold=True, font_size=15, opacity=.2),
            proc_node)
        format_node = self.operation_tv.add_node(TreeViewLabel(
            text="Output Formats", bold=True, font_size=15, opacity=.2),
            proc_node)
        main_file = self.operation_tv.add_node(TreeViewLabel(
            text="Output File", bold=True, font_size=15, opacity=.2),
            proc_node)

        # Save main nodes
        self.main_nodes = {"ortho": ortho_node, "proc": proc_node,
                           "stat": stat_node, "proc_main": main_op_node,
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
            - self.output_formats, to gather information on the output formats
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
            self.operation_tv.add_node(TreeViewLabel(text=text, font_size=16),
                                       parent)
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
        ## for conversion
        if self.main_operations["conversion"]:
            add_node("[Based on input] (main)", self.main_nodes["main_file"])
            if self.main_nodes["main_file"].is_open is False:
                self.operation_tv.toggle_node(self.main_nodes["main_file"])
            ## Output files from secondary operations
            if secondary_op:
                for op in secondary_op:
                    if self.secondary_options["%s_file" % op]:
                        add_node("*_%s (%s)" % (op, op),
                                 self.main_nodes["main_file"])
        ## for concatenation
        elif self.main_operations["concatenation"]:
            if self.output_file == "":
                add_node("[empty] (main)", self.main_nodes["main_file"])
            else:
                add_node("%s (main)" % self.output_file.split(sep)[-1],
                         self.main_nodes["main_file"])
            ## Output files from secondary operations
            if secondary_op:
                for op in secondary_op:
                    if self.secondary_options["%s_file" % op]:
                        add_node("%s_%s (%s)" % (
                            self.output_file.split(sep)[-1], op, op),
                                 self.main_nodes["main_file"])
            if self.main_nodes["main_file"].is_open is False:
                self.operation_tv.toggle_node(self.main_nodes["main_file"])
        else:
            self.main_nodes["main_file"].opacity = .2

    def process_clear_options(self):

        #### CORE CHANGES
        # Clear main operations
        self.main_operations = dict((op, False) for op in self.main_operations)

        # Clear active data set
        self.process_grid_wgt.ids.active_taxa_set.text = "All taxa"

        # Clear output formats
        self.output_formats = ["fasta"]

        # Clear filters, haplotype name and zorro suffix
        self.filter_settings = [25, 50]
        self.hap_prefix = "Hap"
        self.zorro_suffix = ""

        # Clear output file
        self.output_file = ""

        # Clear secondary operations
        self.secondary_operations = dict((op, False) for op in
                                         self.secondary_operations)

        # Clear secondary options
        self.secondary_options = dict((op, False) for op in
                                      self.secondary_options)

        #### APP CHANGES
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
            self.process_options.ids[switch].active = False

        for switch in self.secondary_options:
            self.process_options.ids[switch].active = False

    ############################ PLOT SCREENS ##################################

    def show_plot_toolbar(self):
        """
        Adds a PlotToolbar BoxLayout to self.root_window. This is meant to be an
        auxiliary toolbar for specific operations related to plots.
        """

        # Determine position
        pos = self.root.width - 50, self.root.height - 480
        content = PlotToolbar(pos=pos)

        self.root_window.add_widget(content)

    def dismiss_plot_toolbar(self):
        """
        Removes the PlotToolbar from self.root_window
        """

        try:
            wgt = [x for x in self.root_window.children if
                   isinstance(x, PlotToolbar)][0]
            self.root_window.remove_widget(wgt)
        except IndexError:
            pass

    ########################## ORTHOLOGY SCREEN ################################

    def toggle_orto_soptions(self):
        """
        Controls the toggling of the GridLayout with the advanced options for
        the Orthology screen, Ortholog search slide
        """

        if not self.orto_search_height:
            self.orto_search_height = self.screen.ids.gl_orto_search.height

        if self.screen.ids.adv_search_options.text == "Show additional options":
            # Add widget to main grid
            self.screen.ids.gl_orto_search.add_widget(self.ortho_search_options)
            # Animate widget entrance
            Animation(opacity=1, d=.5, t="in_quad").start(
                self.ortho_search_options)
            # Update button text
            self.screen.ids.adv_search_options.text = "Hide additional options"

            self.screen.ids.gl_orto_search.height = self.orto_search_height + \
                sum(x.height + 5 for x in
                    self.ortho_search_options.ids.mcl_grid.children) + 30

        elif self.screen.ids.adv_search_options.text == \
                "Hide additional options":
            self.screen.ids.gl_orto_search.height = self.orto_search_height
            # Remove widget from main grid
            self.screen.ids.gl_orto_search.remove_widget(
                self.ortho_search_options)
            # Update button text
            self.screen.ids.adv_search_options.text = "Show additional options"

    def dialog_search_report(self, stat_storage):
        """
        Creates the dialog that reports the results of the Orthology search
        :param stat_storage: dictionary. Each entry corresponds to an inflation
        value, which will have a list as a value. The list will contain:

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
            report_wgt.ids.sf_txt.text = str(stats[1])
            report_wgt.ids.sf_box.size_hint_x = float(stats[1]) / \
                                                float(stats[0])

            # Setting final orthologs
            report_wgt.ids.final_txt.text = str(stats[3])
            report_wgt.ids.final_box.size_hint_x = float(stats[3]) / \
                                                float(stats[0])

            # Adding widget to carousel
            content.ids.report_car.add_widget(report_wgt)

        self.show_popup(title="Orthology search report", content=content,
                        size=(400, 470))

    def dialog_import_groups(self):
        """
        Creates filechooser dialog to select group files to be imported
        """

        content = LoadMultipleDialog(cancel=self.dismiss_popup)

        self.show_popup(title="Choose group file(s) to import", content=content,
                        size_hint=(.9, .9))

    def dialog_mysql(self):
        """
        Creates dialog for MySQL settings
        """

        content = MySQLDialog(cancel=self.dismiss_popup)
        content.ids.txt_dlg.text = self.mysql_pass

        self.show_popup(title="MySQL root password", content=content,
                        size=(200, 150))

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
        """

        content = OrtoSetFiltersDialog(cancel=self.dismiss_popup)
        content.group_name = group_name

        self.show_popup(title="Set/change ortholog filters for %s file" %
                              group_name, content=content, size=(400, 250))

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
        """

        try:
            self.orto_max_gene = int(gene_filt)

        except ValueError:
            return self.dialog_floatcheck("Invalid filter value: '%s'. Must be "
                                          "integer" % gene_filt, t="error")
        try:
            self.orto_min_sp = int(sp_filt)

        except ValueError:
            return self.dialog_floatcheck("Invalid filter value: '%s'. Must be "
                                          "integer" % sp_filt, t="error")

        self.dismiss_popup()

        # Add check for min species. If this filter is set to a number
        # greater that the number of proteome input files (which should
        # represent a single species each) this will issue a warning.
        if self.proteome_files and int(sp_filt) > len(self.proteome_files):
            return self.dialog_floatcheck("WARNING: Minimum number of "
                                          "species larger than the provided"
                                          " proteomes", t="error")

    def save_inflation(self, inflation_wgt):
        """
        Save inflation values
        """

        for wgt in inflation_wgt.children:
            if wgt.state == "down" and wgt.text not in self.mcl_inflation:
                self.mcl_inflation.append(wgt.text)
            elif wgt.state == "normal" and wgt.text in self.mcl_inflation:
                self.mcl_inflation.remove(wgt.text)

        self.ortho_search_options.ids.inflation_bt.text = \
            str(sorted(self.mcl_inflation))

    def save_mysql_pass(self, txt):
        """
        Saves mysql access for database creation and manipulation
        """

        # Setup mysql configuration
        er = sql_setup(txt)
        if er:
            return self.dialog_warning("MySQL configuration error", "MySQL "
                                       "setup exited with the following error:"
                                       "\n\n%s" % er)

        self.mysql_pass = txt

        if txt != "":
            self.screen.ids.mysql_bt.text = "Password set"
        else:
            self.screen.ids.mysql_bt.text = "Select..."

    def save_protein_filters(self, min_len, max_stop):
        """
        Saves protein length and stop percentage filters
        """

        self.protein_min_len = min_len
        self.protein_max_stop = max_stop

    def orto_compare_groups(self):
        """
        Switches to the orthology group comparison screen and presents the
        initial plot comparing total orthologs across group files
        """

        # Displays correspondence
        displays = {"total_ort": "1", "sp_ort": "2", "gn_ort": "3",
                    "final_ort": "4"}

        # Get MultiGroup object with the selected groups
        active_groups = ot.MultiGroups()
        orto_screen = join(self.cur_dir, "data", "screens", "Orthology.kv")

        for gchk in self.loaded_screens[orto_screen].ids.group_check.children:
            if gchk.active:
                active_groups.add_group(self.ortho_groups.get_group(gchk.id))

        # Get active displays
        stats = "".join([y for x, y in displays.items()
                         if self.screen.ids[x].active])

        if stats:
            # Create first comparison plot of total orthologs
            self.current_plot, self.current_lgd = active_groups.bar_orthologs(
                dest=self.temp_dir.name, stats=stats)

            # Load plot
            self.orto_compare_loadplot(join(self.temp_dir.name,
                                            "Final_orthologs.png"))

        else:
            self.screen.ids.plot_content.children[0].clear_widgets()

    def orto_compare_loadplot(self, file_path):
        """
        Loads a new plot into the ScatterLayout. This will clear all previous
        content and load a new image based on the file_path argument
        :param file_path: string. Path to the image to be loaded
        :return:
        """

        # Clear previous content
        self.screen.ids.plot_content.children[0].clear_widgets()

        # Add content
        img_wgt = Image(source=file_path, nocache=True)
        self.screen.ids.plot_content.children[0].add_widget(img_wgt)

    def orto_change_filters(self, gn_filter, sp_filter, group_name=None,
                            apply_all=False):
        """

        :param group_name: string. Name of a group object
        :param gn_filter: int. Gene filter
        :param sp_filter: int. Species filter
        :param apply_all: Boolean. If True, apply filters to all group objects
        else, apply to the current group object
        """

        # Set change list according to apply_all argument value
        if apply_all:
            chg_list = [x for x in self.ortho_groups]
        else:
            chg_list = [x for x in self.ortho_groups if
                        x.name.split(sep)[-1] == group_name]

        for g_obj in chg_list:
            g_obj.update_filters(gn_filter, sp_filter)

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

        if len([x for x in self.screen.ids.group_check.children if x.active])\
                >= 2:
            self.screen.ids.compare_group_bt.disabled = False
        else:
            self.screen.ids.compare_group_bt.disabled = True

    def load_groups(self, group_obj=None, group_files=None):
        """
        Loads the group files generated by the Orthology search or manually
        imported into the app. Use group_obj argument for MultiGroup object
        and group_files for file names
        :param group_obj: MultiGroup object
        :param group_files: string. File names of one or more group files
        """

        # Only perform actions if any of these arguments have been used
        if group_obj or group_files:

            # Removes "No groups loaded" button if it still exists
            try:
                self.screen.ids.group_gl.remove_widget(self.screen.ids.no_bt)
            except ReferenceError:
                pass

            # If file names were provided, create Group objects and then the
            # MultiGroup. Else, use the groups object provided
            if group_files:
                group_list = [ot.Group(f) for f in group_files]
                groups = ot.MultiGroups(group_list)
            elif group_obj:
                groups = group_obj

            # Create or update self.ortho_groups
            if not self.ortho_groups or isinstance(self.ortho_groups,
                                                   ObjectProperty):
                self.ortho_groups = groups
            else:
                self.ortho_groups.add_multigroups(groups)

            # Check if any group file is duplicate. If so, issue a warning
            if self.ortho_groups.duplicate_groups:
                self.dialog_warning("Duplicate group files detected",
                                    "The following group file(s) were found "
                                    "to be duplicate and were not loaded:\n\n"
                                    "[b]%s[/b]" %
                    "\n".join(x.name.split(sep)[-1] for x in
                              self.ortho_groups.duplicate_groups))

            # Populate the app gridlayout with group buttons
            for g in groups:

                # Automatically update group filters to the default values of
                # 1 copy per cluster and all taxa represented
                g.update_filters(self.orto_max_gene, len(g.species_list))
                self.orto_min_sp = len(g.species_list)

                # If group name contains full path, get only file name
                gname = g.group_name.split(sep)[-1]

                if g not in self.ortho_groups.duplicate_groups:
                    # Create check box for multiple group selection
                    chk = CheckBox(id=gname, size_hint_x=.1)
                    chk.bind(active=lambda x, y: self.orto_check_state())
                    self.screen.ids.group_check.add_widget(chk)

                    # Create group button
                    bt = ToggleButton(text=gname, id=gname, group="group_bts",
                              size_hint_y=None, height=30, shorten=True,
                              shorten_from="right", halign="center", bold=True,
                              background_down=join("data", "backgrounds",
                                                     "bt_process.png"),
                              background_normal=join("data", "backgrounds",
                                                     "bt_process_off.png"),
                              background_disabled_down=join("data",
                                                            "backgrounds",
                                                            "bt_process.png"),
                              disabled_color=(1, 1, 1, 1))
                    # Apparently I need to use partial instead of lambda
                    # in order to provide a diferent group object as argument
                    # Using lambda will overwrite the group objects of all
                    # buttons with the last group of the iteration. Go figure..
                    bt.bind(on_release=partial(self.orthology_card, g))

                    # Add box to gridlayout
                    self.screen.ids.group_gl.add_widget(bt)

                    # Create removal button
                    x_bt = Button(size_hint=(None, None), width=30, height=30,
                                  id=gname, border=(0, 0, 0, 0),
                                  background_normal=join("data", "backgrounds",
                                                         "remove_bt.png"),
                                  background_down=join("data", "backgrounds",
                                                       "remove_bt_down.png"))
                    x_bt.bind(on_release=partial(self.check_action,
                                                 "Are you sure you want to"
                                                 "remove this group?",
                                                 self.remove_groups))
                    self.screen.ids.group_rm.add_widget(x_bt)

        # If no group button is active, dispatch the first
        if not [x for x in self.screen.ids.group_gl.children
                if x.state == "down"] and self.screen.ids.group_gl.children:
            self.screen.ids.group_gl.children[-1].dispatch("on_release")
            self.screen.ids.group_gl.children[-1].state = "down"

    def orthology_card(self, group_obj, bt=None):
        """
        Generates the descriptive cards with general information for a group
        file.
        :param group_obj: Group object.
        :param bt: ToggleButton instance
        """

        # Create desired behaviour for group toggle buttons
        if bt:
            self.toggle_groups(bt)

        # Get statistics from group object
        stats = group_obj.basic_group_statistics()

        # Create cards
        cards = DescriptionBox(opacity=0)

        cards.prot_txt = str(stats[1])
        cards.ortholog_txt = str(stats[0])
        cards.taxa_txt = str(len(group_obj.species_list))
        cards.group_name = group_obj.name.split(sep)[-1]

        # Create gauge plots, if there are any filtered groups
        if group_obj.species_threshold or group_obj.gene_threshold or \
              (group_obj.species_threshold, group_obj.gene_threshold) == (0, 0):
            # Create species filter plot and add to box
            sp_filter_plot = GaugePlot()
            sp_filter_plot.txt = "After species filter"
            sp_filter_plot.proportion = float(stats[3]) / float(stats[0])
            sp_filter_plot.ortholog_num = str(stats[3])
            cards.ids.gauge_bx.add_widget(sp_filter_plot)

            # Create gene filter plot and add to box
            gn_filter_plot = GaugePlot()
            gn_filter_plot.txt = "After gene filter"
            gn_filter_plot.proportion = float(stats[2]) / float(stats[0])
            gn_filter_plot.ortholog_num = str(stats[2])
            cards.ids.gauge_bx.add_widget(gn_filter_plot)

            # Create final ortholog plot
            final_ortholog_plot = GaugePlot()
            final_ortholog_plot.txt = "Final orthologs"
            final_ortholog_plot.proportion = float(stats[4]) / float(stats[0])
            final_ortholog_plot.ortholog_num = str(stats[4])
            cards.ids.gauge_bx.add_widget(final_ortholog_plot)

        elif not group_obj.filtered_groups and not group_obj.species_threshold\
                and not group_obj.gene_threshold:
            lb = Label(text="Please specify gene and species filters",
                       bold=True, color=(0.216, 0.67, 0.784, 1))
            cards.ids.gauge_bx.add_widget(lb)

        # Clear any previous content from card gridlayout holder if any
        if len(self.screen.ids.card_gl.children) == 1 and \
                isinstance(self.screen.ids.card_gl.children[0], Label):
            self.screen.ids.card_gl.clear_widgets()

        else:
            try:
                Animation(opacity=0, d=.2, t="out_quart").start(
                    self.screen.ids.card_gl.children[0])
            except IndexError:
                pass

            Clock.schedule_once(lambda x: self.screen.ids.card_gl.clear_widgets(),
                                .2)

        # Add card
        Clock.schedule_once(lambda x: self.screen.ids.card_gl.add_widget(cards),
                            .3)
        Clock.schedule_once(lambda x: Animation(opacity=1, d=.3,
                                                t="out_quart").start(cards),
                            .3)

    ########################### POPUP OPS ######################################

    def show_popup(self, title, content, size_hint=(.9, .9), size=None,
                   separator_color=[47 / 255., 167 / 255., 212 / 255., 1.],
                   custom_background=None):
        """
        General purpose method to create a popup widget
        :param title: string. Title of the popup
        :param content: widget object. The contents of the popup widget
        :param size_hint: tuple. Size hint for the widget
        :param size: tuple. The absolute size for the popup. If this argument is
        used, the size_hint will be ignored
        :param custom_background: string. Provide the path to a custom
        background image for the popup.
        """

        # Ignore size_hint is absolute size is provided
        if size:
            self._popup = CustomPopup(title="[b]%s[/b]" % title,
                                content=content, size=size,
                                size_hint=(None, None), auto_dismiss=False,
                                separator_color=separator_color,
                                title_color=separator_color,
                                custom_background=custom_background)
        else:
            self._popup = CustomPopup(title="[b]%s[/b]" % title,
                                content=content, size_hint=size_hint,
                                auto_dismiss=False,
                                separator_color=separator_color,
                                title_color=separator_color,
                                custom_background=custom_background)
        self._popup.open()

    def dismiss_popup(self, *args):
        """
        General purpose method to close popups from the process screen
        """
        self._popup.dismiss()

    def dismiss_subpopup(self, *args):
        """
    General purpose method to close sub-popups from the process screen
        """
        self._subpopup.dismiss()

    ########################### PROCESS SCREEN #################################

    def save_file(self, path, file_name=None, idx=None):
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
        """

        if idx == "main_output":
            if self.main_operations["concatenation"]:
                # Adds output file to storage
                self.output_file = join(path, file_name)
                # Renames the output file button text
                self.process_grid_wgt.ids.conv.text = file_name

            else:
                self.output_dir = path
                self.process_grid_wgt.ids.conv.text = path.split(sep)[-1]

        if idx == "ortho_dir":
            self.ortho_dir = path
            self.screen.ids.orto_dir.text = path.split(sep)[-1]

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
            return self.dialog_floatcheck("WARNING: Please choose at least one "
                                   "output format", t="error")

        self.dismiss_popup()

    def save_filter(self, gap_val, mis_val):
        """
        Stores the information of the FilterDialog
        """

        self.filter_settings = [gap_val,
                                mis_val]

        self.dismiss_popup()

    def dialog_nexus_extra(self):

        content = NexusExtra(cancel=self.dismiss_subpopup)

        content.ids.nexus_check.active = self.use_nexus_partitions

        self._subpopup = Popup(title="Nexus extra options", content=content,
                           size=(500, 200), size_hint=(None, None))

        self._subpopup.open()

    def dialog_phylip_extra(self):

        content = PhylipExtra(cancel=self.dismiss_subpopup)

        content.ids.part_check.active = self.create_partfile

        self._subpopup = Popup(title="Phylip extra options", content=content,
                           size=(400, 200), size_hint=(None, None))

        self._subpopup.open()

    def save_zorro_settings(self, suffix):
        """
        Handles the information provided by the user in the ZorroDialog
        :param suffix: string, suffix of the ZORRO files
        """

        # Check if the zorro files exist
        for f in self.active_file_list:
            f = sep.join(f.split(".")[0:-1])
            f = "%s%s.txt" % (f, suffix)
            if not os.path.isfile(f):
                return self.dialog_floatcheck("ERROR: File %s does not"
                                              " exist" % f, t="error")

        self.update_process_switch("zorro",
                                   self._popup.content.ids.zorro_switch.active)

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
        """

        self.partitions_file = partfile

        self._popup.content.ids.part_file.text = partfile.split(sep)[-1]
        self._popup.content.ids.part_file.background_normal = \
            "data/backgrounds/bt_process.png"

    def dialog_floatcheck(self, text, t):
        """
        Creates a floating label with informative text on the right upper
        corner of the app. This is used for showing errors, warnings and
        general informative messages that fade in and fade out after a time
        :param t: string, with type of check. Can be either error or info
        :param text: string, text to appear in the label
        :return:
        """

        def fade_in():
            Animation(opacity=1, d=.5, t="out_quart").start(check_wgt)
            Animation(opacity=1, d=.5, t="out_quart").start(rm_wgt)

        def fade_out():
            Animation(opacity=0, d=.5, t="out_quart").start(check_wgt)
            Animation(opacity=0, d=.5, t="out_quart").start(rm_wgt)
            Clock.schedule_once(
                lambda dt: self.root_window.remove_widget(check_wgt), 1)
            Clock.schedule_once(
                lambda dt: self.root_window.remove_widget(rm_wgt), 1)

        # Get position from root window
        x, y = self.root_window.size

        # Create widget
        check_wgt = WarningFloat(text=text, opacity=0, markup=True)
        check_wgt.root_pos = [x, y]
        # Create remove button
        rm_wgt = RemoveFloat(pos=(x - 38, y - 75), opacity=0)
        rm_wgt.bind(on_release=lambda arg: fade_out())

        # Determine background color
        if t == "error":
            check_wgt.cl = (1, .33, .33, 1)
            check_wgt.line_cl = (1, 0, 0, 1)
        else:
            check_wgt.cl = (.35, .63, .17, 1)
            check_wgt.line_cl = (0.2, 1, 0.2, 1)

        # Add widget
        self.root_window.add_widget(check_wgt)
        self.root_window.add_widget(rm_wgt)

        # Set animations
        fade_in()
        Clock.schedule_once(lambda arg: fade_out(), 5)

    def check_partitions_file(self):
        """
        This will make some checks on the partitions file provided by the user.
        It will check for errors in the format of the file itself, and whether
        the partitions are correctly defined
        """

        part_obj = data.Partitions()
        er = part_obj.read_from_file(self.partitions_file)

        aln_obj = self.alignment_list.retrieve_alignment(self.rev_infile)
        aln_er = aln_obj.set_partitions(part_obj)

        # Check for the validity of the partitions file
        if isinstance(er, data.InvalidPartitionFile):
            return self.dialog_floatcheck("The provided partitions file is "
                       "invalid. Please check the file or replace with an "
                       "appropriate one.", t="error")

        # Check for the validity of the partitions file
        if isinstance(aln_er, data.InvalidPartitionFile):
            return self.dialog_floatcheck("The provided partitions in the "
                   "partition file do not match the selected alignment",
                    t="error")
        else:
            return True

    def save_reverseconc_settings(self):
        """
        Handles the information provided by the LoadDialog with settings for the
        reverse concatenation
        """

        # Check if a partition file has been selected
        if self.partitions_file == "":
            return self.dialog_floatcheck("Please provide a partitions file",
                                          t="error")

        # Check for the validity of the partitions file
        er = self.check_partitions_file()

        if er is True:

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

        # This will mark the respective buttons for each format as active or not
        # depending on whether they have been previously selected or not. It
        # allows the selected button states to be persistent when visiting the
        # popup multiple times
        for idx, wgt in content.ids.items():
            if idx in self.output_formats:
                wgt.state = "down"
            else:
                wgt.state = "normal"

        # Show popup
        self.show_popup(title="Choose output format", content=content,
                        size=(300, 400))

    def dialog_reverse_inlist(self):

        content = InputList(cancel=self.dismiss_subpopup)

        def set_infile(txt, wgt):

            self.rev_infile = txt
            self._popup.content.ids.rev_infile.background_normal = \
                "data/backgrounds/bt_process.png"
            self._popup.content.ids.rev_infile.text = txt.split(sep)[-1]
            self.dismiss_subpopup()

        if self.file_list:
            for infile in self.file_list:
                bt = Button(text=infile.split(sep)[-1], size_hint_y=None,
                            height=30, bold=True, shorten=True,
                            shorten_from="right", halign="center",
                            valign="middle")
                bt.text_size[0] = bt.size[0] * 3.5
                bt.bind(on_release=partial(set_infile, infile))
                content.ids.rev_inlist.add_widget(bt)

        self._subpopup = Popup(title="Choose input file", content=content,
                               size_hint=(.5, .8))

        self._subpopup.open()

    def dialog_reverse_concatenation(self, title="Choose input file"):
        """
        Generates a general purpose file chooser to request additional data
        :param title: string, A custom title for the load data dialog popup
        """

        content = RevConcDialog(cancel=self.dismiss_popup)
        content.ids.rev_conc.active = \
            self.main_operations["reverse_concatenation"]

        # Check if partitions file was already selected. If so, update the
        # corresponding button
        if self.partitions_file:
            content.ids.part_file.background_normal = \
                "data/backgrounds/bt_process.png"
            content.ids.part_file.text = \
                self.partitions_file.split(sep)[-1]

        # Check if input file to reverse concatenate was already selected. If
        # so, update the corresponding button
        if self.rev_infile:
            content.ids.rev_infile.background_normal = \
                "data/backgrounds/bt_process.png"
            content.ids.rev_infile.text = \
                self.rev_infile.split(sep)[-1]

        self.show_popup(title=title, content=content, size=(300, 350))

    def dialog_load_partfile(self):

        content = LoadDialog(cancel=self.dismiss_subpopup)

        # If input files have already been provided, use their directory as a
        # starting point for the partition file chooser. Otherwise, use the
        # home path
        if self.file_list:
            content.ids.ld_filechooser.path = sep.join(
                self.file_list[0].split(sep)[:-1])
        else:
            content.ids.ld_filechooser.path = self.home_path

        self._subpopup = Popup(title="Choose partition file", content=content,
                               size_hint=(.9, .9))

        self._subpopup.open()

    def dialog_filechooser(self, idx=None):
        """
        Generates a file chooser popup for the user to select an output file

        :param idx: string. An id of where the filechooser is calling. This
        allows the addition of custom behaviours for different dialogs
        """

        # Inherits the layout defined in the .kv file under <SaveDialog>
        content = SaveDialog(cancel=self.dismiss_popup)

        # Custom behaviour for main output file chooser dialog
        if idx == "main_output":
            # If the main operation is conversion or reverse concatenation,
            # remove the TextInput asking for the file name
            if self.main_operations["conversion"] or \
                    self.main_operations["reverse_concatenation"]:
                content.ids.txt_box.clear_widgets()
                content.ids.txt_box.height = 0
                title = "Choose destination directory of output file(s)"
            else:
                title = "Choose output file"

        if idx == "export":
            title = "Choose file for exportation"

        # Custom behaviour for orthology output directory
        if idx == "ortho_dir":
            content.ids.txt_box.clear_widgets()
            content.ids.txt_box.height = 0
            title = "Choose destination directory for OrthoMCL output files"

        # Save output file for conversion/concatenation purposes
        # Providing this operation will allow the filechooser widget to
        # know which output file is this
        content.ids.sd_filechooser.text = idx

        self.show_popup(title=title, content=content)

    def dialog_filter(self):
        """
        Generates the settings popup for filtering options
        """

        content = FilterDialog(cancel=self.dismiss_popup)
        # Update filter values if they were changed
        if self.filter_settings:
            content.ids.gap_sli.value = self.filter_settings[0]
            content.ids.mis_sli.value = self.filter_settings[1]

        self.show_popup(title="Set filter thresholds", content=content,
                        size=(400, 300))

    @staticmethod
    def check_filters(value):
        """
        Method that validates the input of the text input in filter settings.
        It handles common mistakes, such as using "," instead of "." for
        decimal places and truncates values between the range of 0 and 100.
        If the text input cannot be converted to float, it will return false
        and the slider value will not change
        :param value: text_input.text
        :return:
        """

        try:
            x = float(value.replace(",", "."))
            if x > 100:
                corrected_val = 100
            elif x < 0:
                corrected_val = 0
            else:
                corrected_val = x
            return True, corrected_val
        except ValueError:
            return False

    def dialog_execution(self):
        """
        Generates the dialog for Process execution
        """

        content = ExecutionDialog(cancel=self.dismiss_popup)
        aln_obj = self.update_active_fileset(deepcopy(self.alignment_list))

        # Get main operation
        try:
            main_op = [nm for nm, bl in self.main_operations.items()
                           if bl is True][0]
            content.ids.main_op.text = "[b][size=18][color=37abc8ff]Main " \
                                       "operation:[/color][/size][/b] %s" % \
                                       main_op
        except IndexError:
            return self.dialog_floatcheck("ERROR: Please select a main "
                                          "operation", t="error")

        # Get secondary operations
        secondary_op = [nm for nm, bl in self.secondary_operations.items()
                       if bl is True]
        if secondary_op:
            content.ids.sec_op.text = "[b][size=18][color=37abc8ff]Secondary " \
                                   "operation(s):[/color][/size][/b] %s" %\
                                  ", ".join(secondary_op)
        else:
            content.ids.sec_op.text = "[b][size=18][color=37abc8ff]Secondary " \
                                      "operation(s):[/color][/size][/b] None"

        # Get output formats
        content.ids.out_form.text = "[b][size=18][color=37abc8ff]Output " \
                                   "format(s):[/color][/size][/b] %s" %\
                                  ", ".join(self.output_formats)

        # Get output files
        # In case concatenation
        if main_op == "concatenation":
            # Check if an output directory has been selected
            if self.output_file == "":
                return self.dialog_floatcheck("ERROR: No output file has been "
                                       "selected", t="error")
            out_file = self.output_file.split(sep)[-1]
            add_files = [out_file + "_" + nm for nm, bl in
                         self.secondary_operations.items() if bl]
            content.ids.out_files.text = "[b][size=18][color=37abc8ff]Output " \
                                       "file(s):[/color][/size][/b] (%s) %s"\
                                       ", %s" % (len(add_files) + 1, out_file,
                                       ", ".join(add_files))
        # In case conversion
        if main_op == "conversion":
            # Check if an output file has been selected
            if self.output_dir == "":
                return self.dialog_floatcheck("ERROR: No output directory has"
                                       " been selected", t="error")
            try:
                # Check for additional files
                add_files = [nm for nm, bl in
                             self.secondary_operations.items() if bl]
                content.ids.out_files.text = "[b][size=18][color=37abc8ff]" \
                            "Output file(s):[/color][/size][/b] %s converted " \
                            "file(s)" % \
                    (len(aln_obj.alignment_object_list) +
                    len(aln_obj.alignment_object_list) * len(add_files))
            # In case aln_obj has not being defined, probably because there
            # are no input files
            except AttributeError:
                return self.dialog_floatcheck("ERROR: No input files have"
                                              "been selected", t="error")

        try:
            self.show_popup(title="Process execution summary - Processing %s "
                "file(s)" % len(aln_obj.alignment_object_list),
                            content=content, size=(550, 350))
        except AttributeError:
            return self.dialog_floatcheck("ERROR: No input files have"
                                              "been selected", t="error")

    def dialog_text(self, title, idx):
        """
        Generates a simple text dialog to capture text input
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

        if idx == "new_folder":
            self._subpopup = Popup(title=title, content=content,
                                   size=(200, 150), size_hint=(None, None))
            self._subpopup.open()
        else:
            self.show_popup(title=title, content=content,
                            size=(200, 150))

    def dialog_warning(self, msg1, msg2):

        content = WarningDialog(cancel=self.dismiss_subpopup)
        content.ids.warning_l.text = "[b][color=#ff5555ff][size=18]%s[/size]" \
                                     "[/color][/b]\n\n%s" % (msg1, msg2)

        self._subpopup = CustomPopup(title="Error!", content=content,
                        separator_color=[255 / 255., 85 / 255., 85 / 255., 1.],
                        size_hint=(None, None), size=(550, 300))

        self._subpopup.open()

    def dialog_zorro(self):

        content = ZorroDialog(cancel=self.dismiss_popup)

        content.ids.zorro_switch.active = self.secondary_options["zorro"]
        content.ids.zorro_txt.text = self.zorro_suffix

        self.show_popup(title="ZORRO support", content=content,
                        size=(350, 200))

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
            self.usearch_evalue = text

        elif idx == "orto_group":
            self.ortholog_prefix = text

        elif idx == "groups":
            self.group_prefix = text

    def update_main_operations(self, op):
        """
        Updates the app attribute containing the main operations of the Process
        screen, self.main_operations. Only one main operation can be active.
        :param op: The name of the operation to turn on (all others will be
        disabled)
        """

        """
        The text of the Output file/directory field changes depending on whether
        the main operation is a concatenation or a conversion
        """
        file_text = "[size=18][b]Output file[/b][/size]\n[size=13]Save "\
                          "output file to the selected file.[/size]"
        dir_text = "[size=18][b]Output directory[/b][/size]\n[size=13]" \
                          "Save output file(s) to the selected directory." \
                          "[/size]"

        self.main_operations = {k: True if k == op else False for k in
                                self.main_operations}

       # Disables output file button and other conversion/concatenation
        # specific buttons
        if op == "conversion":
            if self.output_dir == "":
                self.process_grid_wgt.ids.conv.text = "Select..."
            else:
                self.process_grid_wgt.ids.conv.text = \
                    self.output_dir.split(sep)[-1]
            self.process_grid_wgt.ids.output_label.text = dir_text
            Animation(height=0, d=.32, t="in_quad").start(
                self.screen.ids.sub_conc)
        elif op == "concatenation":
            if self.output_file == "":
                self.process_grid_wgt.ids.conv.text = "Select..."
            else:
                self.process_grid_wgt.ids.conv.text = \
                    self.output_file.split(sep)[-1]
            self.process_grid_wgt.ids.output_label.text = file_text
            Animation(height=50, d=.32, t="in_quad").start(
                self.screen.ids.sub_conc)
        elif op == "reverse_concatenation":
            if self.output_dir == "":
                self.process_grid_wgt.ids.conv.text = "Select..."
            else:
                self.process_grid_wgt.ids.conv.text = \
                    self.output_dir.split(sep)[-1]
            self.process_grid_wgt.ids.output_label.text = dir_text

    def save_main_operation(self, op):
        """
        This controls the appearance of the general options after the user
        chooses the main operation (Conversion vs. Concatenation). When the
        main operation is chosen for the first time, the general options
        are introduced, but further selection of the main operation only
        changes the state of the operations button.
        :param op: string, the main operations. values are "concatenation" and
        "conversion"
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
        Controls the toggling of the GridLayout with the additional options for
        the process screen.
        """

        # Shows additional options
        if self.process_grid_wgt.ids.opt_bt.text == "Show additional options":

            # Adds widget to the Process widget tree an animates its
            # appearance via changes in opacity to give impression of fade in
            self.process_grid_wgt.add_widget(self.process_options)
            Animation(opacity=1, d=.5, t="in_quad").start(self.process_options)

            # Update the height of the GridLayout to allow scrolling
            self.process_grid_wgt.height = self.process_height + (50 * len(
                self.process_options.ids.filter_grid.children))

            # Change text in the toggle button
            self.process_grid_wgt.ids.opt_bt.text = "Hide additional options"

        # Hides additional options
        elif self.process_grid_wgt.ids.opt_bt.text == "Hide additional options":

            # Removes widget from widget tree. However, all settings that were
            # change while the widget was visible will be stored
            self.process_grid_wgt.remove_widget(self.process_options)

            # Change text in the toggle button
            self.process_grid_wgt.ids.opt_bt.text = "Show additional options"

    ###################################
    #
    # CORE RELATED METHODS AND FUNCTIONS
    #
    ###################################

    def load_files_subproc(self, files):

        def load_proc(nm):

            aln_obj = AlignmentList([])
            nm.c = 0
            for f in files:
                nm.c += 1
                nm.m = f.split(sep)[-1]
                aln_obj.add_alignment_file(f)
            else:
                nm.alns = aln_obj
            return

        def check_proc(p, dt):

            try:
                content.ids.pb.value = ns.c
                content.ids.msg.text = "Processing file %s" % ns.m
            except AttributeError:
                pass

            if not p.is_alive():
                alns = ns.alns
                Clock.unschedule(func)
                self.dismiss_popup()
                self.load_files(files, alns)

            if content.proc_kill:
                d.terminate()
                self.dismiss_popup()
                Clock.unschedule(func)

        manager = multiprocessing.Manager()
        ns = manager.Namespace()

        #parent_con, child_con = multiprocessing.Pipe()
        d = multiprocessing.Process(target=load_proc, args=(ns, ))

        d.start()

        time.sleep(.1)

        if not d.is_alive():
            alns = ns.alns
            return self.load_files(files, alns)

        content = LoadProgressDialog()
        content.ids.pb.max = len(files)
        self.show_popup(title="Loading files", content=content,
                        size=(400, 300))

        func = partial(check_proc, d)
        Clock.schedule_interval(func, .1)

    def load_files(self, selection, aln_list):
        """
        Loads the selected input files into the program using the
        AlignmentList object. If there are invalid alignment objects in the
        AlignmentList object, they will be returned by this method for
        handling in the load method

        I also sets the alignment and taxa active lists, which should be used
        to perform subsequent operations. For now, the original taxa and
        alignment lists should remain untouched for future implementation of
        undo/redo functionality and as backups
        :returns: List of invalid/badly formatted alignment objects
        """

        #aln_list = AlignmentList(files)

        # Check for consistency in sequence type across alignments
        current_seq_type = set(self.sequence_types + aln_list.format_list())
        if len(current_seq_type) > 1:
            return InputTypeError("Multiple sequence types detected")
        else:
            self.sequence_types.extend(list(current_seq_type))

        # When creating an AlignmentList object, some input alignment may be
        # invalid, in which case they are removed from the
        # alignment_object_list. This will handle the case where all input files
        # are invalid
        if aln_list.alignment_object_list:
            # In case the alignment list object is already populated with
            # previously loaded files, then add to the object
            if self.alignment_list:
                for aln_obj in aln_list:
                    self.alignment_list.add_alignment(aln_obj)
                    self.active_alignment_list.add_alignment(aln_obj)
                # Update active taxa list
                self.active_taxa_list = self.active_alignment_list.taxa_names

            # In case the alignment list object is empty, load it as is
            else:
                # Populate original alignment list
                self.alignment_list = aln_list
                # Updating active alignment list
                self.active_alignment_list = deepcopy(self.alignment_list)
                self.active_taxa_list = self.active_alignment_list.taxa_names

        self.load(selection, aln_list.bad_alignments)

    def get_taxa_information(self, alt_list=None):
        """
        This method will gather all available information for all taxa and set
        a number of related attributes. The way the method is implemented,
        allow the generation of information for both complete (if the method
        is applied in the original data set) and active (if the method is
        applied to the currently data set) data sets.

        :param: alt_list: This argument provides a way to override the
        self.active_alignment_list that is used by default. This may be helpful
        when the complete/original list has been modified and the contents of
        the taxa popup have to reflect those modifications. If no alternative
        is provided, the method will use the self.active_taxa_list.

        :return: tx_inf (dictionary): Contains the relevant information for
        the taxa popup. All corresponding values are the result of additions
        across all active input alignments. Contains the following keys:

            - length: Total alignment length including missing data and gaps
            - indel: The number columns containing indel/gaps
            - missing: The number of columns containing missing data
            - effective_len: Alignment length excluding gaps and missing data
            - effective_len_per: Same as above but in percentage
            - fl_coverage: The number of files containing the focal taxon
            - fl_coverage_per: Same as above but in percentage
        """

        # main storage defined
        tx_inf = {}

        if alt_list:
            aln_list = alt_list
        else:
            aln_list = self.active_alignment_list

        for tx in self.active_taxa_list:

            # Add entry to storage dictionary
            tx_inf[tx] = {}
            # Counter for alignment missing the taxa
            tx_missing = 0

            # Get full sequence
            sequence = ""
            # This assures that the information will only be gathered if the
            # active data set is not empty
            if aln_list.alignment_object_list:
                for aln in aln_list:
                    if tx in aln.alignment:
                        sequence += aln.alignment[tx]
                    else:
                        tx_missing += 1
                else:
                    # Retrieve missing data symbol
                    missing_symbol = aln.sequence_code[1]

                # Get sequence length
                seq_len = len(sequence)
                tx_inf[tx]["length"] = seq_len

                # Get indel number.
                tx_inf[tx]["indel"] = sequence.count("-")

                # Get missing data
                tx_inf[tx]["missing"] = sequence.count(missing_symbol)

                # Get effective sequence length in absolute and percentage
                tx_inf[tx]["effective_len"] = seq_len - (tx_inf[tx]["indel"] +
                                                          tx_inf[tx]["missing"])
                tx_inf[tx]["effective_len_per"] = round(
                    (tx_inf[tx]["effective_len"] * 100) / seq_len, 2)

                # Get number of files containing the taxa in absolute and
                # percentage
                tx_inf[tx]["fl_coverage"] = len(
                    aln_list.alignment_object_list) - \
                    tx_missing
                tx_inf[tx]["fl_coverage_per"] = round(((
                    tx_inf[tx]["fl_coverage"] * 100) / len(
                    aln_list.alignment_object_list)), 2)

            else:
                # This handles the case where the active data set is empty
                tx_inf[tx]["length"] = "NA"
                tx_inf[tx]["indel"] = "NA"
                tx_inf[tx]["missing"] = "NA"
                tx_inf[tx]["effective_len"] = "NA"
                tx_inf[tx]["effective_len_per"] = "NA"
                tx_inf[tx]["fl_coverage"] = "NA"
                tx_inf[tx]["fl_coverage_per"] = "NA"

        return tx_inf

    def get_file_information(self, mode="alignment"):
        """
        Similar to get_taxa_information, but generating information for the
        files in the file tab.
        :param mode: string. The type of file information to retrieve. May be
        'alignment' or 'proteome'

        :return: file_inf (dictionary). Contains all relevant content for the
        file popup. It contains the following keys:

            ..:mode: alignment:
                - aln_format: The format of the input file
                - seq_type: The sequence type. If DNA, RNA, Protein.
                - n_taxa: Number of taxa
                - aln_len: Length of the alignment
                - is_aln: If the input file is in alignment format of
                non-aligned sequence set format
                - model: The model of sequence evolution, if applicable. This is
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
                for infile in self.file_list:
                    file_name = infile.split(sep)[-1]
                    file_inf[file_name] = {}

                    # Get alignment object
                    aln = self.alignment_list.retrieve_alignment(infile)

                    # Get input format
                    file_inf[file_name]["aln_format"] = aln.input_format

                    # Get sequence type
                    file_inf[file_name]["seq_type"] = aln.sequence_code[0]

                    # Get sequence model
                    # if aln.model:
                    #     file_inf[file_name]["model"] = " ".join(aln.model)
                    # else:
                    #     file_inf[file_name]["model"] = "NA"

                    # Get number of species
                    file_inf[file_name]["n_taxa"] = len([x for x in
                        aln.iter_taxa() if x in self.active_taxa_list])

                    # Get if is alignment
                    file_inf[file_name]["is_aln"] = str(aln.is_alignment)
                    # if aln.model:
                    #     file_inf[file_name]["is_aln"] += " (Concatenated)"

                    # Get length of largest sequence if not aligned, or alignment
                    # length
                    file_inf[file_name]["aln_len"] = aln.locus_length

        elif mode == "proteome":
            for p in self.proteome_files:
                p_name = p.split(sep)[-1]
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
        :param switch_id: string, name of the switch according to the keys in
        process_switches
        :param state: Boolean, current state of the corresponding switch
        """

        if switch_id in self.secondary_operations:
            self.secondary_operations[switch_id] = state
        else:
            self.secondary_options[switch_id] = state

    def update_active_fileset(self, aln_obj):
        """
        This method is similar in purpose and functionality to the
        update_active_taxaset, but it updates the set of files. It should be
        used before the said method.
        :param aln_obj: The alignment object being used during execution of
        Process operatios
        """

        # Determine the selected active taxa set from the dropdown menu
        file_set_name = self.process_grid_wgt.ids.active_file_set.text

        if file_set_name == "All files":
            return aln_obj

        if file_set_name == "Active files":
            return self.active_alignment_list

        else:
            file_set = [self.filename_map[x] for x in
                        self.file_groups[file_set_name]]
            file_set = set(self.file_list) - set(file_set)
            aln_obj.remove_file(file_set)
            return aln_obj

    def update_active_taxaset(self, aln_obj):
        """
        Do not use this method on the original self.alignment_list or
        self.active_alignment list, as it may cause unwanted permanent changes
        to the taxa set.

        This will take the complete taxa set from self.alignment_list.taxa_names
        and the currently active taxa set from self.active_taxa_list and remove
        the all taxa that are not present in the active taxa set from the
        alignment object passed as argument. If the lists are the same, no taxa
        will be removed
        """

        # Determine the selected active taxa set from the dropdown menu
        tx_set_name = self.process_grid_wgt.ids.active_taxa_set.text

        if tx_set_name == "All taxa":
            return aln_obj

        if tx_set_name == "Active taxa":
            tx_set = self.active_taxa_list

        else:
            tx_set = self.taxa_groups[tx_set_name]

        # Remove taxa
        aln_obj.remove_taxa(list(set(self.alignment_list.taxa_names) -
                                 set(tx_set)))

        return aln_obj

    def dialog_orto_execution(self):
        """
        Creates and populates the pre-execution dialog for orthology search
        """

        # Check for input proteomes
        if not self.proteome_files:
            return self.dialog_floatcheck("Please provide proteome files as "
                                          "input data", t="error")

        # Check for mysql pass
        if self.mysql_pass == "":
            return self.dialog_floatcheck("Please provide the access password "
                                          "to MySQL", t="error")

        # Check for output directory
        if self.ortho_dir == "":
            return self.dialog_floatcheck("Please specify an output directory"
                                          "for orthology results", t="error")

        content = OrtoExecutionDialog(cancel=self.dismiss_popup)

        content.ids.gene_filter.text = "[b][size=18][color=37abc8ff]Maximum" \
                                       " number of gene copies per cluster:" \
                                       "[/color][/size][/b] %s" % \
                                       self.orto_max_gene

        content.ids.sp_filter.text = "[b][size=18][color=37abc8ff]Minimum" \
                                       " number of taxa per cluster:[/color]" \
                                       "[/size][/b] %s" % self.orto_min_sp

        content.ids.eval.text = "[b][size=18][color=37abc8ff]USEARCH e-value" \
                                " threshold:[/color][/size][/b] %s" % \
                                self.usearch_evalue

        content.ids.inflation.text = "[b][size=18][color=37abc8ff]MCL " \
                                     "inflation value(s):[/color][/size][/b] "\
                                     "%s" % ", ".join(self.mcl_inflation)

        content.ids.threads.text = "[b][size=18][color=37abc8ff]Threads " \
                                   ":[/color][/size][/b] %s" % \
                                   self.screen.ids.usearch_threads.text

        self.show_popup(title="Orthology search execution summary - Processing"
                              " %s file(s)" % len(self.proteome_files),
                        content=content, size=(550, 350))

    def orthology_search_exec(self):
        """
        Main function that executes all queued procedures of the orthology
        module
        """

        def process_dispatch(nm):
            """
            Executes all pipeline subprocesses sequentially and updates the
            Progess dialog label
            """

            if nm.k:
                nm.t = "Installing schema"
                nm.c = 1
                opipe.install_schema("orthomcl.config")

            if nm.k:
                nm.t = "Adjusting Fasta Files"
                nm.c = 2
                db = opipe.adjust_fasta(self.proteome_files)

            if nm.k:
                nm.t = "Filtering Fasta Files"
                nm.c = 3
                opipe.filter_fasta(self.protein_min_len, self.protein_max_stop)

            if nm.k:
                nm.t = "Running USearch"
                nm.c = 4
                opipe.allvsall_usearch("goodProteins.fasta",
                                  self.usearch_evalue,
                                  self.screen.ids.usearch_threads.text,
                                  self.usearch_output)

            if nm.k:
                nm.t = "Parsing USEARCH output"
                nm.c = 5
                opipe.blast_parser(self.usearch_output)

            if nm.k:
                opipe.remove_duplicate_entries()

            if nm.k:
                nm.t = "Loading USEARCH output to database"
                nm.c = 6
                opipe.load_blast("orthomcl.config")

            if nm.k:
                nm.t = "Obtaining Pairs"
                nm.c = 7
                opipe.pairs("orthomcl.config")

            if nm.k:
                opipe.dump_pairs("orthomcl.config")

            if nm.k:
                nm.t = "Running MCL"
                nm.c = 8
                opipe.mcl(self.mcl_inflation)

            if nm.k:
                nm.t = "Dumping groups"
                nm.c = 9
                opipe.mcl_groups(self.mcl_inflation, self.ortholog_prefix,
                                 "1000", self.group_prefix)

            if nm.k:
                nm.t = "Filtering group files"
                nm.c = 10
                stats, groups_obj = opipe.export_filtered_groups(
                                             self.mcl_inflation,
                                             self.group_prefix,
                                             self.orto_max_gene,
                                             self.orto_min_sp, db)

                # stats is a dictionary containing the inflation value as key
                # and a list with the orthologs as value
                nm.stats = stats
                nm.groups = groups_obj

        def check_process(p, dt):
            """
            Checks the status of the background process "p" and updates
            the progress dialog label
            """

            # Updates progress dialog label
            content.ids.msg.text = ns.t
            # Updates progress bar
            content.ids.pb.value = ns.c

            # When the process finishes, close progress dialog and unschedule
            # this callback
            if not p.is_alive():
                Clock.unschedule(func)
                self.dismiss_popup()
                self.dialog_search_report(ns.stats)
                self.load_groups(ns.groups)

            # Listens for cancel signal
            if content.proc_kill:
                ns.k = False
                self.dismiss_popup()
                Clock.unschedule(func)

        # Create Progression dialog
        content = OrtoProgressDialog()
        self.show_popup(title="Running Orthology Search", content=content,
                        size=(400, 200))

        # Setup multiprocess
        # The manager will allow shared variables between independent processes
        # so that the progress dialog label can be updated with the current
        # pipeline status
        manager = multiprocessing.Manager()
        ns = manager.Namespace()
        ns.k = True

        # Create Process instance
        d = multiprocessing.Process(name="daemon", target=process_dispatch,
                                    args=(ns, ))
        # This will make the process run in the background so that the app
        # doesn't freeze
        d.daemon = True

        # Change working directory
        os.chdir(self.ortho_dir)

        # Create directory that will store intermediate files during orthology
        # search
        int_dir = "backstage_files"
        if not os.path.exists(int_dir):
            os.makedirs(int_dir)

        os.chdir(int_dir)

        # Create orthomcl_config
        create_orthomcl_cfg(".")

        # Start pipeline in the background
        d.start()

        # Check status of pipeline and update progress dialog
        func = partial(check_process, d)
        Clock.schedule_interval(func, .1)

    def process_exec(self):
        """
        Main function that executes all queued procedures of the process module
        """

        write_aln = {}

        # Perform checks
        if self.alignment_list is None or not\
                self.alignment_list.alignment_object_list:
            return self.dialog_warning("No input data", "Use 'Menu > Open "
                                       "file(s)' to load input data")

        #####
        # Perform operations
        #####

        # Setting the alignment to use. A deepcopy of the alignment list
        # is used because it may be possible to do changes in the taxa data set
        # of the AlignmentList object, which should not change the original
        # self.alignment_list. This is because when taxa are removed from
        # the alignment list, there is no way to return those taxa to the
        # object
        # Update active file set of the alignment object
        aln_object = self.update_active_fileset(deepcopy(self.alignment_list))
        # Update active taxa set of the alignment object
        aln_object = self.update_active_taxaset(aln_object)
        proc_files = len(aln_object.alignment_object_list)

        # Concatenation
        if self.main_operations["concatenation"]:
            aln_object = aln_object.concatenate()

            # Concatenation of ZORRO files
            if self.secondary_options["zorro"]:
                zorro_data = data.Zorro(aln_object, self.zorro_suffix)
                zorro_data.write_to_file(self.output_file)

        # Reverse concatenation
        if self.main_operations["reverse_concatenation"]:
            partition_obj = data.Partitions()
            # In case the partitions file is badly formatted or invalid, the
            # exception will be returned by the read_from_file method.
            er = partition_obj.read_from_file(self.partitions_file)
            aln_object = aln_object.retrieve_alignment(self.rev_infile)
            aln_object.set_partitions(partition_obj)
            aln_object = aln_object.reverse_concatenate()

        # Collapsing
        if self.secondary_operations["collapse"]:
            if self.secondary_options["collapse_file"]:
                collapsed_aln_obj = deepcopy(aln_object)
                collapsed_aln_obj.collapse(haplotype_name=self.hap_prefix,
                                           haplotypes_file="_collapsed")
                write_aln[self.output_file + "_collapsed"] = collapsed_aln_obj
            else:
                aln_object.collapse(haplotype_name=self.hap_prefix)

        # Filtering
        if self.secondary_operations["filter"]:
            if self.secondary_options["filter_file"]:
                filtered_aln_obj = deepcopy(aln_object)
                filtered_aln_obj.filter_missing_data(self.filter_settings[0],
                                                 self.filter_settings[1])
                write_aln[self.output_file + "_filtered"] = filtered_aln_obj
            else:
                aln_object.filter_missing_data(self.filter_settings[0],
                                               self.filter_settings[1])

        # Gcoder
        if self.secondary_operations["gcoder"]:
            if self.secondary_options["gcoder_file"]:
                gcoded_aln_obj = deepcopy(aln_object)
                gcoded_aln_obj.code_gaps()
                write_aln[self.output_file + "_coded"] = gcoded_aln_obj
            else:
                aln_object.code_gaps()

        # The output file(s) will only be written after all the required
        # operations have been concluded. The reason why there are two "if"
        # statement for "concatenation" is that the input alignments must be
        # concatenated before any other additional operations. If the first
        # if statement did not exist, then all additional options would have
        # to be manually written for both "conversion" and "concatenation". As
        # it is, when "concatenation", the aln_obj is firstly converted into
        # the concatenated alignment, and then all additional operations are
        # conducted in the same aln_obj
        write_aln[self.output_file] = aln_object
        if self.main_operations["concatenation"]:
            if self.output_file == "":
                return self.dialog_warning("Output file not selected",
                                           "Use the 'Select...' button of "
                                           "'Output file' general option to "
                                           "select an output file name")
            for name, obj in write_aln.items():
                obj.write_to_file(self.output_formats, name,
                                interleave=self.secondary_options["interleave"],
                                partition_file=self.create_partfile,
                                use_charset=self.use_nexus_partitions)
        else:
            if self.output_dir == "":
                return self.dialog_warning("Output directory not specified",
                                           "use the 'Select...' button of "
                                           "'Output directory' general option"
                                           " to specify a destination directory"
                                           " for the output file(s)")
            else:
                for name, obj in write_aln.items():
                    name = name.replace(self.output_file, "")
                    obj.write_to_file(self.output_formats, output_suffix=name,
                                interleave=self.secondary_options["interleave"],
                                partition_file=self.create_partfile,
                                output_dir=self.output_dir,
                                use_charset=self.use_nexus_partitions)

        if proc_files == 1:
            self.dialog_floatcheck("All Done! %s file was successfully "
                                   "processed" % proc_files, t="info")
        else:
            self.dialog_floatcheck("All Done! %s files were successfully "
                                   "processed" % proc_files, t="info")

if __name__ == '__main__':
    TriFusionApp().run()