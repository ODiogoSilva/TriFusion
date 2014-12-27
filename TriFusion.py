#!/usr/bin/env python3
# -*- coding: utf-8 -*-#
#
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
#
# NEEDS FIXING:
# Issue 1. The size of the scrollview for the files and taxa tabs does not
# update perfectly with the addition of taxa/files buttons
#

# Kivy imports
from kivy.app import App
from kivy.uix.togglebutton import ToggleButton
from kivy.uix.button import Button
from kivy.animation import Animation
from kivy.uix.popup import Popup
from kivy.uix.label import Label
from kivy.uix.codeinput import CodeInput
from kivy.uix.textinput import TextInput
from kivy.uix.rst import RstDocument
from kivy.uix.boxlayout import BoxLayout
from kivy.clock import Clock
#from kivy.uix.gridlayout import GridLayout
#from kivy.uix.scrollview import ScrollView
from kivy.factory import Factory
from kivy.uix.floatlayout import FloatLayout
from kivy.lang import Builder
from kivy.properties import NumericProperty, StringProperty, BooleanProperty,\
    ListProperty, ObjectProperty
from kivy.uix.screenmanager import Screen

# Main program imports
from process.sequence import AlignmentList

# Other imports
from os.path import dirname, join, exists
from os.path import expanduser
from copy import deepcopy
import pickle


class ShowcaseScreen(Screen):
    fullscreen = BooleanProperty(False)

    def add_widget(self, *args):
        if 'content' in self.ids:
            return self.ids.content.add_widget(*args)
        return super(ShowcaseScreen, self).add_widget(*args)


class SaveDialog(FloatLayout):
    save = ObjectProperty(None)
    text_input = ObjectProperty(None)
    cancel = ObjectProperty(None)


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

    # Getting current directory to fetch the screen kv files
    cur_dir = dirname(__file__)

    # Setting the list of input files variable
    # The path ONLY attribute
    path = StringProperty()
    # Only the original input files. SHOULD NOT BE MODIFIED
    file_list = ListProperty()
    # Dynamic list containing only the activated files
    active_file_list = ListProperty()
    # Dictionary mapping file names to their corresponding full paths. This
    # attribute must exist, because some parts of the code need only the file
    #  name instead of the full path, but a connection to the full path must
    # be maintained for future reference
    filename_map = {}

    # Setting the list of taxa names
    active_taxa_list = ListProperty()

    # Attributes to know current and previous screen
    current_screen = StringProperty()
    previous_screen = StringProperty()

    # Attribute to load screens
    index = NumericProperty(-1)

    # Attribute to the home path
    home_path = expanduser("~")

    # Attribute with bookmarks for file chooser. The bookmark attribute is a
    # list containing a list with the full path of the bookmarks as the
    # first element and a dictionary mapping the full path to the bookmark
    # name as the second element
    bookmarks = [[], {}]
    bm_file = cur_dir + "/data/resources/bookmarks"

    _popup = ObjectProperty(None)

    ################################
    #
    # CORE PROGRAM RELATED VARIABLES
    #
    ################################

    # List storing the original alignment object variables. SHOULD NOT BE
    # MODIFIED
    alignment_list = None
    # List of active alignment object variables.
    active_alignment_list = None

    # Attributes containing the original and active taxa information
    # dictionaries
    original_tx_inf = None
    active_tx_inf = None
    # Same as before but for file information
    original_file_inf = None
    active_file_inf = None

    # Attribute containing the objects for the several possible output files.
    output_files = {"conversion": None, "collapse": None, "gcoder": None}

    ##################################
    #
    # GUI RELATED METHODS AND FUCTIONS
    #
    ##################################

    def build(self):

        # Setting main window title
        self.title = "TriFusion - Streamline phylogenomics"

        # Setting available screens
        self.available_screens = ["main", "Orthology", "Process",
                                  "Statistics", "fc"]
        self.screen_names = self.available_screens

        # Transforming screen names into complete paths to be loaded by kivy
        self.available_screens = [join(self.cur_dir, "data", "screens",
                                 "{}.kv".format(screen)) for screen in
                                  self.available_screens]

        # First thing is go to main screen
        self.go_screen(0)

    def go_screen(self, idx, direct="left"):
        """
        Method used to go to a specific screen by specifying and index and
        transition direction
        :param idx: integer. Index value of the screen from self.screen_names
        :param direct: string. The direction of the transition
        """

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
        self.screen = Builder.load_file(self.available_screens[idx])

        # If the screen to be loaded is the filechooser, set the home path as
        # the default
        if self.available_screens[idx].split("/")[-1] == "fc.kv":
            self.screen.ids.icon_view_tab.path = self.home_path
            # Initialize bookmarks
            self.init_bookmark()

        return self.screen

    def init_bookmark(self):
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
        will be saved as a pickle object
        :param path: String containing the path of the bookmark
        """

        # Load bookmarks object
        self.bookmarks = pickle.load(open(self.bm_file, "rb"))
        # Add bookmarks to the full path list
        self.bookmarks[0].append(path)
        # Add mapping of the full path to the bookmark name
        new_map = {path.split("/")[-1]: path}
        self.bookmarks[1] = dict(list(self.bookmarks[1].items()) +
                                 list(new_map.items()))
        self.add_bookmark_bt(path)
        pickle.dump(self.bookmarks, open(self.bm_file, "wb"))

    def add_bookmark_bt(self, bk):

        bookmark_name = bk.split("/")[-1]
        bt = Button(text=bookmark_name, id=bk,
                    height=self.root.height * 0.05, size_hint=(.8, None))
        bt.bind(on_release=self.load_bookmark)
        xbt = Button(text="X", size_hint=(.14, None),
                     height=self.root.height * 0.05, id="%sX" % bk,
                     background_color=(255, .9, .9, 1), bold=True)
        xbt.bind(on_release=self.remove_bookmark_bt)
        self.screen.ids.sv_book.add_widget(bt)
        self.screen.ids.sv_book.add_widget(xbt)

    def load_bookmark(self, value):

        path = value.id
        self.screen.ids.icon_view_tab.path = path

    def remove_bookmark_bt(self, value):

        parent_obj = value.parent

        bk_idx = value.id[:-1]
        bk_bt = [x for x in parent_obj.children if bk_idx == x.id][0]

        parent_obj.remove_widget(value)
        parent_obj.remove_widget(bk_bt)

        # Core changes
        bk_name = bk_idx.split("/")[-1]
        self.bookmarks[0].remove(bk_idx)
        del self.bookmarks[1][bk_name]
        pickle.dump(self.bookmarks, open(self.bm_file, "wb"))

    def side_panel_toggle(self):
        """
        Method controlling the animation toggling of the side panel
        """

        # Toggling the state of the panel
        self.show_side_panel = not self.show_side_panel

        if self.show_side_panel:

            # Redraw the side panel layout. This will ensure that the widget
            # is always on top of all widgets.
            self.root.ids.bx1.remove_widget(self.root.ids.panel_float)
            self.root.ids.bx1.add_widget(self.root.ids.panel_float)

            # The width of the side panel contents will be relative to the
            # widget width
            sv_panel_width = self.root.width * .32
            # The width of the side panel buttons will be the same as the
            # actionprevious button of the action bar
            sv_bts_width = self.root.ids.ap.children[0].children[-1].width
        else:
            sv_panel_width, sv_bts_width = 0, 0

        ## ANIMATIONS with hierarchy
        # Animation of main BoxLayout containing child ScrollViews
        Animation(width=sv_panel_width * 1.2, d=.3, t="out_quart").start(
            self.root.ids.main_box)
        # Animation of both scrollviews
        Animation(width=sv_panel_width, d=.3, t="out_quart").start(
            self.root.ids.sp)
        Animation(width=sv_bts_width, d=.3, t="out_quart").start(
            self.root.ids.sp_bts)

    def load(self, selection):
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

        if self.file_list:
            # Updating complete and active file lists
            self.file_list.extend(selection)
            self.active_file_list.extend(selection)
            # Update the filename - path mapping attribute
            self.filename_map = dict(list(self.filename_map.items()) +
                list((x, y) for x, y in zip([x.split("/")[-1] for x in
                                             selection], selection)))

        else:
            # Set an attribute with the input file list
            self.file_list = selection
            # Setting active file list and path list
            self.active_file_list = deepcopy(self.file_list)
            # Sett the filename - path mapping attribute
            self.filename_map = dict((x, y) for x, y in zip(
                [x.split("/")[-1] for x in selection], selection))

        # Parses the files into the program
        self.load_files(selection)
        # Update active taxa list
        self.update_taxa()
        # Populates files and taxa contents
        self.update_tabs()
        # Gathers taxa  and file information
        self.original_tx_inf = self.get_taxa_information()
        self.original_file_inf = self.get_file_information()

    def update_tabs(self):
        """
        Wrapper that updates the contents of the files and taxa tabs
        """

        self.populate_input_files()
        self.populate_species()

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

    def populate_input_files(self):
        """
        This method grabs the input files that were selected in the
        FileChooser widget and populates the File tab in the main side panel
        with toggle and remove buttons for each file
        """

        # Remove the initial disabled button, if it's still there
        if "file_temp" in self.root.ids.keys():
            self.root.ids.file_sl.remove_widget(self.root.ids.file_temp)
            del self.root.ids["file_temp"]

        # Enable selection buttons if file list is not empty
        if self.file_list:
            for i in self.root.ids.sb_file.children:
                i.disabled = False
                i.bind(on_release=self.select_bt)

        for infile in self.file_list:

            file_name = infile.split("/")[-1]
            # This prevents duplicate files from being added
            if file_name not in [x.id for x in self.root.ids.file_sl.children]:

                bt = ToggleButton(text=file_name, state="down", id=file_name,
                                  height=self.root.height * 0.05,
                                  size_hint=(.8, None), shorten=True,
                                  shorten_from="right", halign="center")
                # Setting horizontal text size for shortening
                bt.text_size[0] = bt.size[0] * 1.3
                # Binding functionality to toggle button
                bt.bind(on_release=self.toggle_selection)

                # Adds toggle button with file name
                self.root.ids.file_sl.add_widget(bt)

                # Set Information button and add the widget
                inf_bt = Button(text="?", size_hint=(.14, None),
                                height=self.root.height * 0.05,
                                id="%s?" % file_name, bold=True)
                self.root.ids.file_sl.add_widget(inf_bt)
                inf_bt.bind(on_release=self.popup_info)

                # Set remove button with event binded and add the widget
                x_bt = Button(text="X", size_hint=(.14, None),
                              height=self.root.height * 0.05, id="%sX" %
                              file_name, background_color=(255, .9, .9, 1),
                              bold=True)
                x_bt.bind(on_release=self.remove_bt)
                self.root.ids.file_sl.add_widget(x_bt)

                # Updates the size of the grid layout according to the added
                # buttons
                self.root.ids.file_sl.height += self.root.height * .068

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

        # Enable selection buttons if taxa list is not empty
        if self.active_taxa_list:
            for i in self.root.ids.sb_taxa.children:
                i.disabled = False
                i.bind(on_release=self.select_bt)

        for tx in sorted(self.active_taxa_list):

            # Prevents duplicate taxa from being entered
            if tx not in [x.id for x in self.root.ids.taxa_sl.children]:

                bt = ToggleButton(text=tx, state="down", id=tx,
                                  height=self.root.height * 0.05,
                                  size_hint=(.8, None), shorten=True,
                                  shorten_from="right", halign="center")
                # Setting horizontal text size for shortening
                bt.text_size[0] = bt.size[0] * 1.3
                # Binding functionality to toggle button
                bt.bind(on_release=self.toggle_selection)

                # Add toggle button with taxa name
                self.root.ids.taxa_sl.add_widget(bt)

                # Set Information button and add the widget
                inf_bt = Button(text="?", size_hint=(.14, None),
                                height=self.root.height * 0.05,
                                id="%s?" % tx, bold=True)
                self.root.ids.taxa_sl.add_widget(inf_bt)
                inf_bt.bind(on_release=self.popup_info)

                # Set remove button with event binded and add the widget
                x_bt = Button(text="X", size_hint=(.14, None),
                              height=self.root.height * 0.05, id="%sX" % tx,
                              background_color=(255, .9, .9, 1), bold=True)
                self.root.ids.taxa_sl.add_widget(x_bt)
                x_bt.bind(on_press=self.remove_bt)

                # Updates the size of the grid layout according to the added
                # button
                self.root.ids.taxa_sl.height += self.root.height * 0.068

    def dismiss_popup(self):
        self._popup.dismiss()

    def save_file(self, path, filename, idx):

        file_handle = open(join(path, filename), "w")

        self.output_files[idx] = file_handle

        self.screen.ids.conversion.text = filename

        self.dismiss_popup()

    def popup_filechooser(self, value):
        """
        Generates a file chooser popup for the user to select an output file
        """

        content = SaveDialog(cancel=self.dismiss_popup)
        content.ids.sd_filechooser.path = self.home_path

        if value == "conversion":
            content.ids.sd_filechooser.text = "conversion"

        self._popup = Popup(title="Choose output file",
                 content=content, size_hint=(.9, .9))

        self._popup.open()

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

                # For now, the pop up content will be in a CodeInput widget
                # because it is the only widget (along with TextInput) that
                # allow text selection and it may be even possible to add text
                # formatting using python lexer.
                content = CodeInput(text=" -- Complete data set -- \n"
                                    "Sequence length: %s\n"
                                    "Number of indels: %s\n"
                                    "Number missing data: %s\n"
                                    "Effective sequence length: %s (%s%%)\n"
                                    "File coverage: %s (%s%%)\n\n"
                                    " -- Active data set -- \n"
                                    "Sequence length: %s\n"
                                    "Number of indels: %s\n"
                                    "Number missing data: %s\n"
                                    "Effective sequence length: %s (%s%%)\n"
                                    "File coverage: %s (%s%%)\n" % (
                                 self.original_tx_inf[tx]["length"],
                                 self.original_tx_inf[tx]["indel"],
                                 self.original_tx_inf[tx]["missing"],
                                 self.original_tx_inf[tx]["effective_len"],
                                 self.original_tx_inf[tx]["effective_len_per"],
                                 self.original_tx_inf[tx]["fl_coverage"],
                                 self.original_tx_inf[tx]["fl_coverage_per"],
                                 self.active_tx_inf[tx]["length"],
                                 self.active_tx_inf[tx]["indel"],
                                 self.active_tx_inf[tx]["missing"],
                                 self.active_tx_inf[tx]["effective_len"],
                                 self.active_tx_inf[tx]["effective_len_per"],
                                 self.active_tx_inf[tx]["fl_coverage"],
                                 self.active_tx_inf[tx]["fl_coverage_per"]),
                                    readonly=True)

                popup_wgt = Popup(title="Taxon: %s" % value.id[:-1],
                                  content=content, size_hint=(None, None),
                                  size=(400, 400))

                popup_wgt.open()

        elif value.parent == self.root.ids.file_sl:

            # Get file name
            file_name = value.id[:-1]

            if self.filename_map[file_name] in self.active_file_list:

                # Get the information from the content list. This is done when
                # calling the popup to avoid repeating this operation every time
                # taxa  or files are added/removed.
                self.active_file_inf = self.get_file_information()
                content = CodeInput(text="Input format: %s\n"
                                         "Sequence type: %s\n"
                                         "Alignment: %s\n"
                                         "Sequence size: %s\n"
                                         "Number of taxa: %s\n"
                                         "Model: %s\n" % (
                             self.original_file_inf[file_name]["aln_format"],
                             self.original_file_inf[file_name]["seq_type"],
                             self.original_file_inf[file_name]["is_aln"],
                             self.active_file_inf[file_name]["aln_len"],
                             self.active_file_inf[file_name]["n_taxa"],
                             self.active_file_inf[file_name]["model"]),
                                    read_only=True)

                popup_wgt = Popup(title="File: %s" % file_name,
                                  content=content, size_hint=(None, None),
                                  size=(400, 400))

                popup_wgt.open()

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

        # Changes concerning the taxa tab
        if parent_obj == self.root.ids.taxa_sl:

            # When button is normal (unselected) remove from active list
            if value.state == "normal":
                self.active_taxa_list.remove(value.id)
            # When button is down (selected) add to active
            elif value.state == "down":
                self.active_taxa_list.append(value.id)

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

        # Updates the size of the grid layout according to the removed button
        parent_obj.height -= self.root.height * .0715

        ####### CORE CHANGES
        # Get the parent tab

        if parent_obj == self.root.ids.file_sl:
            # Update file list
            file_path = self.filename_map[bt_idx]
            self.file_list.remove(file_path)
            # Update alignment object list
            self.alignment_list.remove_file([self.filename_map[bt_idx]])
            self.active_alignment_list.remove_file([self.filename_map[bt_idx]])

            # Update active taxa list. This must be executed before calling
            # self.get_taxa_information since this method relies on an
            # updated active taxa list
            self.update_taxa()

            # Update pop up content. Since the file has been removed,
            # it should also be excluded from the complete data set
            self.original_tx_inf = self.get_taxa_information(
                alt_list=self.alignment_list)


        if parent_obj == self.root.ids.taxa_sl:
            self.active_alignment_list.remove_taxa([bt_idx])
            self.active_taxa_list = self.active_alignment_list.taxa_names

    def select_bt(self, value):
        """
        Functionality to the Select All/Deselect All buttons of the side
        panel. The method was made in such a way that it could be of general
        use for buttons in the files and taxa tabs
        """

        sv_parent = [x for x in value.parent.parent.children if "scrollview" in
                     str(x.__class__)][0]

        # This will iterate over the first child of the parent scrollview.
        # Since scroll view only supports one child, this should be fine
        for i in sv_parent.children[0].children:
            # Skips the X buttons
            if "togglebutton" in str(i.__class__):

                if value.text == "Select All":
                    # App related action
                    i.state = "down"

                    # Core changes to files
                    if sv_parent == self.root.ids.sv_file:
                        self.active_file_list = deepcopy(self.file_list)
                        self.active_alignment_list = deepcopy(
                            self.alignment_list)
                    #Core changes to taxa
                    if sv_parent == self.root.ids.sv_sp:
                        self.active_taxa_list = deepcopy(
                            self.alignment_list.taxa_names)

                elif value.text == "Deselect All":
                    # App related action
                    i.state = "normal"

                    # Core changes to files
                    if sv_parent == self.root.ids.sv_file:
                        self.active_file_list = []
                        self.active_alignment_list.clear_files()
                    # Core changes to taxa
                    if sv_parent == self.root.ids.sv_sp:
                        self.active_taxa_list = []

    ###################################
    #
    # CORE RELATED METHODS AND FUNCTIONS
    #
    ###################################

    def load_files(self, files):
        """
        Loads the selected input files into the program using the
        AlignmentList object.

        I also sets the alignment and taxa active lists, which should be used
        to perform subsequent operations. For now, the original taxa and
        alignment lists should remain untouched for future implementation of
        undo/redo functionality and as backups
        """

        aln_list = AlignmentList(files)

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

    def get_file_information(self):
        """
        Similar to get_taxa_information, but generating information for the
        files in the file tab.

        :return: file_inf (dictionary). Contains all relevant content for the
        file popup. It contains the following keys:

            - aln_format: The format of the input file
            - seq_type: The sequence type. If DNA, RNA, Protein.
            - n_taxa: Number of taxa
            - aln_len: Length of the alignment
            - is_aln: If the input file is in alignment format of non-aligned
            sequence set format
            - model: The model of sequence evolution, if applicable. This is
            usually only present on Nexus input format
        """

        # main storage
        file_inf = {}

        # Iterating over file path since the name of the Alignment
        # objects contain the full path

        if self.file_list:
            for infile in self.file_list:
                file_name = infile.split("/")[-1]
                file_inf[file_name] = {}

                # Get alignment object
                aln = self.active_alignment_list.retrieve_alignment(infile)

                # Get input format
                file_inf[file_name]["aln_format"] = aln.input_format

                # Get sequence type
                file_inf[file_name]["seq_type"] = aln.sequence_code[0]

                # Get sequence model
                if aln.model:
                    file_inf[file_name]["model"] = " ".join(aln.model)
                else:
                    file_inf[file_name]["model"] = "NA"

                # Get number of species
                file_inf[file_name]["n_taxa"] = len([x for x in aln.iter_taxa()
                                                if x in self.active_taxa_list])

                # Get if is alignment
                file_inf[file_name]["is_aln"] = str(aln.is_alignment)

                # Get length of largest sequence if not aligned, or alignment
                # length
                file_inf[file_name]["aln_len"] = aln.locus_length

        return file_inf

if __name__ == '__main__':
    TriFusionApp().run()