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

# Kivy imports
from kivy.app import App
from kivy.uix.togglebutton import ToggleButton
from kivy.uix.button import Button
from kivy.animation import Animation
#from kivy.uix.gridlayout import GridLayout
#from kivy.uix.scrollview import ScrollView
from kivy.lang import Builder
from kivy.properties import NumericProperty, StringProperty, BooleanProperty,\
    ListProperty
from kivy.uix.screenmanager import Screen

# Main program imports
from process.sequence import Alignment, AlignmentList

# Other imports
from os.path import dirname, join


class ShowcaseScreen(Screen):
    fullscreen = BooleanProperty(False)

    def add_widget(self, *args):
        if 'content' in self.ids:
            return self.ids.content.add_widget(*args)
        return super(ShowcaseScreen, self).add_widget(*args)


class TriFusionApp(App):

    #######################
    #
    # GUI RELATED VARIABLES
    #
    #######################

    # Setting Boolean controlling the toggling of main headers
    show_side_panel = BooleanProperty(False)

    # Variable containing screen names
    screen_names = ListProperty([])
    available_screens = ListProperty([])

    # Getting current directory to fetch the screen kv files
    cur_dir = dirname(__file__)

    # Setting the list of input files variable
    # Only the file names
    file_list = ListProperty([])
    # The complete path
    file_path_list = ListProperty([])

    # Attributes to know current and previous screen
    current_screen = StringProperty()
    previous_screen = StringProperty()

    # Attribute to load screens
    index = NumericProperty(-1)

    ################################
    #
    # CORE PROGRAM RELATED VARIABLES
    #
    ################################

    # List storing alignment object variables
    alignment_list = []

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
        screen = Builder.load_file(self.available_screens[idx])
        return screen

    def side_panel_toggle(self):
        """
        Method controlling the animation toggling of the side panel
        """

        # Toggling the state of the panel
        self.show_side_panel = not self.show_side_panel

        if self.show_side_panel:
            width = self.root.width * .32
            self.root.ids.sv_but.text = "Open File(s)"
        else:
            width = 0
            self.root.ids.sv_but.text = ""

        Animation(width=width, d=.3, t="out_quart").start(self.root.ids.sp)
        # Animate the button so that the folding of the panel is smoother
        Animation(width=width * .8, d=.3, t="out_quart").start(
            self.root.ids.sv_but)

    def load(self, selection):
        """
        Loads selected input files into the program. This should be switched
        after selecting files in a FileChooser widget. The selected files
        will be parsed and the side panel will be populated with information
        on file names and taxa names
        :param selection: list. Contains the paths to the selected input files
        """

        self.file_path_list = selection

        # Setting a list containing only the file names
        self.file_list = [x.split("/")[-1] for x in selection]

        # Populates the files tab in the side panel
        self.populate_input_files()
        # Parses the files into the program
        self.load_files()
        # Populates the taxa tab in the side panel
        self.populate_species()

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

        for infile in self.file_list:

            # This prevents duplicate files from being added
            if infile not in [x.id for x in self.root.ids.file_sl.children]:

                bt = ToggleButton(text=infile, state="down", id=infile,
                                  height=self.root.height * .05,
                                  size_hint=(.8, None))

                # Updates the size of the grid layout according to the added
                # buttons
                self.root.ids.file_sl.height += self.root.height * .05

                # Adds both a toggle button and a button to remove the file
                self.root.ids.file_sl.add_widget(bt)
                self.root.ids.file_sl.add_widget(ToggleButton(text="X",
                     size_hint_x=.2))

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

        for tx in self.alignment_list.taxa_names:

            # Prevents duplicate taxa from being entered
            if tx not in [x.id for x in self.root.ids.taxa_sl.children]:

                bt = ToggleButton(text=tx, state="down", id=tx,
                                  height=self.root.height * 0.05,
                                  size_hint=(.8, None))

                # Updates the size of the grid layout according to the added
                # button
                self.root.ids.taxa_sl.height += self.root.height * 0.05

                self.root.ids.taxa_sl.add_widget(bt)
                x_bt = Button(text="X", size_hint=(.2, None),
                              height=self.root.height * 0.05, id="%sX" % tx)
                x_bt.bind(on_press=self.remove_bt)
                self.root.ids.taxa_sl.add_widget(x_bt)

    def remove_bt(self, value):

        bt_idx = value.id[:-1]
        bt = [x for x in self.root.ids.taxa_sl.children if bt_idx == x.id][0]
        self.root.ids.taxa_sl.remove_widget(value)
        self.root.ids.taxa_sl.remove_widget(bt)

    def load_files(self):

        self.alignment_list = AlignmentList(self.file_path_list)


if __name__ == '__main__':
    TriFusionApp().run()