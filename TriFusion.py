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

# Kivy imports
from kivy.app import App
from kivy.uix.togglebutton import ToggleButton
from kivy.uix.button import Button
from kivy.animation import Animation
from kivy.uix.popup import Popup
from kivy.uix.scrollview import ScrollView
from kivy.uix.widget import Widget
from kivy.uix.label import Label
from kivy.core.window import Window
from kivy.uix.codeinput import CodeInput
from kivy.uix.textinput import TextInput
from kivy.uix.rst import RstDocument
from kivy.uix.boxlayout import BoxLayout
from kivy.clock import Clock
from kivy.uix.gridlayout import GridLayout
#from kivy.uix.scrollview import ScrollView
from kivy.factory import Factory
from kivy.uix.tabbedpanel import TabbedPanel
from kivy.uix.floatlayout import FloatLayout
from kivy.lang import Builder
from kivy.properties import NumericProperty, StringProperty, BooleanProperty,\
    ListProperty, ObjectProperty
from kivy.uix.screenmanager import Screen
from kivy.config import Config
from kivy.clock import Clock
from kivy.core.text.markup import MarkupLabel
from kivy.uix.treeview import TreeView, TreeViewLabel

# Main program imports
from process.sequence import AlignmentList
from process import data

# Other imports
from os.path import dirname, join, exists
from os import sep
from collections import OrderedDict
from os.path import expanduser
from copy import deepcopy
from functools import partial
import pickle

Config.set("kivy", "log_level", "warning")


class ShowcaseScreen(Screen):
    fullscreen = BooleanProperty(False)

    def add_widget(self, *args):
        if 'content' in self.ids:
            return self.ids.content.add_widget(*args)
        return super(ShowcaseScreen, self).add_widget(*args)


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

    # Getting current directory to fetch the screen kv files
    cur_dir = dirname(__file__)

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
    bm_file = join(cur_dir, "data", "resources", "bookmarks")

    _popup = ObjectProperty(None)
    _subpopup = ObjectProperty(None)

    # Dictionary containing the values for the main process operations
    main_operations = {"concatenation": False, "conversion": False}

    # Dictionary containing all values of the switches and checkboxes in the
    # process screen
    secondary_operations = OrderedDict([("collapse", False),
                                    ("filter", False),
                                    ("gcoder", False),
                                    ("reverse_concatenation", False)])

    secondary_options = OrderedDict([("interleave", False),
                                    ("zorro", False),
                                    ("collapse_file", False),
                                    ("filter_file", False),
                                    ("gcoder_file", False)])

    # Attribute for the gridlayout widget that will contain all main options
    # for the process module
    process_grid_wgt = None
    process_options = None
    process_height = None

    # Attribute for the widget containint the treeview showing the operations
    # queue
    operation_tv = ObjectProperty(None)
    main_nodes = {}

    # Attributes storing the toggle buttons from Taxa/File panels. Mostly for
    # mouse_over events
    # Contains the button widgets from the Files and Taxa tabs
    mouse_over_bts = {"Files": [], "Taxa": []}
    # The button text of the previous mouse over event. This will allow the
    # assessment of whether the current mouse collision is for the same button
    # (in which case the mouse over will not be triggered) or for a different
    # button (in which case the mouse over is triggered)
    previous_mouse_over = ""
    # This is a locking mechanism of the mouse over event. When there is a
    # scheduled event for a mouse over this attribute is set to False, which
    # prevents further events from being scheduled in the meantime. When the
    # scheduled event is dispatched, the lock is released and it returns to
    # True
    mouse_over_ready = True
    # Stores the previous mouse over label button so that it can be removed
    old_mouse_over = None

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

    # Attribute for taxa groups
    taxa_groups = {}

    # Attribute containing the objects for the several possible output files.
    output_file = ""

    # Attribute storing active output formats. Fasta is True by default
    output_formats = ["fasta"]

    # Attributes for extra options of output formats
    # Determines whether the part.File associated with phylip format is created
    create_partfile = True

    # Attribute storing the filter settings. The list should contain gap
    # threshold as first element and missing data threshold as second element
    filter_settings = [25, 50]

    # Attribute for ZORRO settings
    zorro_suffix = ""

    # Attribute storing the haplotype prefix
    hap_prefix = "Hap"

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
                                  "Statistics", "fc"]
        self.screen_names = self.available_screens

        # Transforming screen names into complete paths to be loaded by kivy
        self.available_screens = [join(self.cur_dir, "data", "screens",
                                 "{}.kv".format(screen)) for screen in
                                  self.available_screens]
        self.loaded_screens = dict((sc, None) for sc in self.available_screens)

        # First thing is go to main screen
        self.go_screen(0)

        Window.bind(on_touch_up=self.sidepanel_on_touch)
        Window.bind(on_key_down=self._on_keyboard_events)
        #self.sine_panel_routine()

        # Creating GridLayout instance for general options of Process
        self.process_grid_wgt = ProcessGeneral()

        # Create GridLayout instance for additional options of Process.
        self.process_options = AdditionalProcessContents()

        # Initialize operation queue treeview in side panel
        self.operation_queue_init()

        Clock.schedule_interval(self._on_mouseover_tabs, .1)

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

        # Keybinding ctrl+f that brings focus to the "Find" field in the
        # Filechooser screen
        if modifier == "ctrl" and key == b'\x06' and self.screen.name == "fc":
            self.screen.ids.text_filter.focus = True

        # Keybinding ctrl+o that opens the Filechooser screen
        if modifier == "ctrl" and key == b'\x0f':
            self.go_screen(self.screen_names.index("fc"))

    def _on_mouseover_tabs(self, dt):

        # Get mouse position coordinates
        mp = self.root_window.mouse_pos
        # Set collision attribute
        collision = False

        def show_label(mouse, wgt, *args):
            """
            Use this function with a Clock schedule to delay the introduction
            of the label widget. Otherwise, it could become cumbersome to have
            an label appearing right after the mouse colliding with the button
            :param wgt: FloatLayout widget, containing a descriptive label
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
            window_pos = wgt.to_window(bt.pos[0], bt.pos[1])
            # Creating dummy widget to determine collision
            dummy_wgt = Widget(pos=window_pos, size_hint=(None, None),
                                 size=wgt.size)

            return dummy_wgt.collide_point(mp[0], mp[1])

        def create_label_wgt(text):
            """
            Creates the label to be introduced in the floatlayout.
            :param text, string, text of a button
            """

            info_bt = Button(text=text, pos=mp, size=(len(text) * 8, 40),
                        size_hint=(None, None))

            return info_bt

        # Only do this routine if the side panel is open
        if self.show_side_panel and self.mouse_over_ready:
            # Get active tab in side panel
            active_tab = self.root.ids.main_tp.current_tab.text
            # Iterate over buttons of active tab
            for bt in self.mouse_over_bts[active_tab]:
                # Determine if there is a collision with mouse position
                if determine_collision(bt):
                    # Set collision marker to true
                    collision = True

                    # This will determine if a new label button will be added
                    # to the layout, based on the text of the button. If the
                    # text is already in the previous mouse over, then do
                    # nothing. If the text is some new button, then do something
                    if bt.text != self.previous_mouse_over:
                        # Check if there is an old label button and remove it
                        if self.old_mouse_over:
                            self.root_window.remove_widget(self.old_mouse_over)

                        # Create label widget
                        label = create_label_wgt(text=bt.text)

                        # Schedule the introduction of the label widget
                        Clock.schedule_once(partial(show_label, mp, label), .8)
                        # Locking mouse over so that no additional label widgets
                        # are added during the waiting time
                        self.mouse_over_ready = False

            else:
                # If no collision is detected, remove any remaining label widget
                if collision is False and \
                   self.old_mouse_over in self.root_window.children:
                    self.root_window.remove_widget(self.old_mouse_over)

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
        else:
            self.screen = Builder.load_file(self.available_screens[idx])

            # If the screen to be loaded is the filechooser, set the home path
            #  as the default
            if self.available_screens[idx].split("/")[-1] == "fc.kv":
                self.screen.ids.icon_view_tab.path = self.home_path
                # Initialize bookmarks
                self.bookmark_init()

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
        for i in wgt.parent.children:
            if i.disabled:
                i.disabled = False
                i.state = "normal"

        wgt.disabled = True

    def check_action(self, func, bt_wgt):
        """
        General purpose method that pops a dialog checking if the user wants to
        perform a certain action. This method should be passed as a function on
        the 'on_*' with the final function and original widget triggering the
        event as arguments
        :param func: final function if the users choses to process
        :param bt_wgt: widget where the initial 'on_' event occurred

        Usage example:
        This can be applied to the removal button of the bookmarks. In this
        case, the event of the removal button must be like this:

        remove_button.bind(partial(self.check_action, self.remove_bookmark_bt))

        where, self.check_action is this method, and self.remove_bookmark_bt is
        the function that will actually remove the bookmak button. This function
        is then binder to the "OK" button of the check dialog. By default, the
        last argument is the bt_wgt.

        """

        check_content = CheckDialog(cancel=self.dismiss_popup)
        check_content.ids.check_ok.bind(on_release=lambda val: func(bt_wgt))

        self.show_popup(title="Warning!", content=check_content,
                        size=(200, 150))

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
            new_map = {path.split("/")[-1]: path}
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
                    height=30, size_hint=(.8, None))
        # Bind to function that loads bookmark path into filechooser
        bt.bind(on_release=self.bookmark_load)
        # Define bookmark removal button
        xbt = Button(text="X", size_hint=(.14, None),
                     height=30, id="%sX" % bk,
                     background_color=(255, .9, .9, 1), bold=True)
        # Bind to function that removes bookmark button as well as the path
        # from self.bm_file
        xbt.bind(on_release=partial(self.check_action, self.remove_bookmark_bt))
        # Update gridlayout height
        self.screen.ids.sv_book.height += self.root.height * 0.07
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

    ########################### SIDE PANEL EXP #################################

    def toggle_midpanel(self, width):

        Animation(width=width, d=.32, t="out_quart").start(
            self.root.ids.sv_over)

    def sine_panel_routine(self):

        def toggle_mouse_over(dt):

            mp = self.root_window.mouse_pos

            first_bt = self.root.ids.first_sidebt
            second_bt = self.root.ids.second_sidebt

            if first_bt.collide_point(mp[0], mp[1]):
                over_width = self.root.width * .325

            elif second_bt.collide_point(mp[0], mp[1]):
                over_width = self.root.width * .325

            else:
                over_width = 0

            self.toggle_midpanel(over_width)

        Clock.schedule_interval(toggle_mouse_over, 1)

    ######################## SIDE PANEL OPERATIONS #############################

    def sidepanel_on_touch(self, *args):
        """
        This function is binded to the app Window so that it can handle any
        touch_up events. Once the side panel is open, this allows any mouse
        click outside the panel to close it. It gathers information on the
        mouse and side panel position and evaluates a collision. It will
        trigger the side panel closing only when four conditions are met:

         - When there is a mouse input outside the side panel
         - When the variable controling the side panel (show_side_panel) is
         True, meaning that the panel is extended
         - When the mouse input is outside the previous button in the action
         bar, which is also used to toggle the side panel. This prevents issues
         of toggling the side panel twice with one mouse input
         - When a popup is not open. There are several buttons in the side bar
         that open whose position is outside the side bar. The user should be
         able to click anywhere in the popup without the side panel closing.
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
        mous_pos = self.root_window.mouse_pos
        # Get side panel and previous button widgets
        side_panel_wgt = self.root.ids.main_box
        ap = self.root.ids.ap

        # Check for conditions to close the side panel.
        # If touch is out of panel; if panel is open; is touch is out of menu
        # button; a popup is not open
        if side_panel_wgt.collide_point(mous_pos[0], mous_pos[1]) is False\
                and self.show_side_panel \
                and ap.collide_point(mous_pos[0], mous_pos[1]) is False \
                and self._popup not in self.root_window.children:

            if self.screen.name == "Process":
                queue_bt = self.screen.ids.queue_bt
                if queue_bt.collide_point(mous_pos[0], mous_pos[1]) is False:
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
            sv_bts_width = 50
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

        # Parses the files into the program
        bad_aln = self.load_files(selection)

        # Checking if there are invalid input alignments
        if bad_aln:
            if len(bad_aln) == 1:
                self.dialog_warning("Invalid input file detected",
                                    "The following input file is invalid:\n\n"
                                    "[b]%s[/b]" %
                                    bad_aln[0].name.split(sep)[-1])
            else:
                self.dialog_warning("Invalid input files detected",
                                    "The following input files are invalid:\n\n"
                                    "[b]%s[/b]" %
                                    "\n".join(x.name.split(sep)[-1] for x
                                              in bad_aln))

        # removes bad alignment files from selection list
        selection = [path for path in selection if path not in
                    [x.name for x in bad_aln]]

        if self.file_list:
            # Updating complete and active file lists
            self.file_list.extend(selection)
            self.active_file_list.extend(selection)
            # Update the filename - path mapping attribute
            self.filename_map = dict(list(self.filename_map.items()) +
                                     list((x, y) for x, y in
                                     zip([x.split("/")[-1] for x in selection],
                                         selection)))

        else:
            # Set an attribute with the input file list
            self.file_list = selection
            # Setting active file list and path list
            self.active_file_list = deepcopy(self.file_list)
            # Sett the filename - path mapping attribute
            self.filename_map = dict((x, y) for x, y in zip(
                [x.split("/")[-1] for x in selection], selection))

        if self.active_alignment_list:
            # Update active taxa list
            self.update_taxa()
            # Populates files and taxa contents
            self.update_tabs()
            # Gathers taxa  and file information
            self.original_tx_inf = self.get_taxa_information()
            self.original_file_inf = self.get_file_information()

            # Unlocks options dependent on the availability of input data
            self.root.ids.tx_group_bt.disabled = False

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

    def update_file_label(self):
        """
        Sets and updates a label on the Files tab of the side panel, informing
        how many files are selected out of the total files
        :return:
        """

        self.root.ids.file_lab.text = "%s of %s files selected" % (
                                       len(self.active_file_list),
                                       len(self.file_list))

        # Check if there are 0 out of 0 files. In this case, disabled the
        # select/unselect all buttons
        if len(self.active_file_list) == 0 and len(self.file_list) == 0:
            self.root.ids.sb_file.disabled = True
        else:
            self.root.ids.sb_file.disabled = False

    def update_sp_label(self):
        """
        Sets and updates a label on the Taxa tab of the side panel, informing
        how many taxa are selected out of the total taxa
        :return:
        """
        self.root.ids.sp_lab.text = "%s of %s taxa selected" % (
                                       len(self.active_taxa_list),
                                       len(self.alignment_list.taxa_names))

        # Check if there are 0 out of 0 taxa. In this case, disabled the
        # select/unselect all buttons
        if len(self.active_taxa_list) == 0 and \
                len(self.alignment_list.taxa_names) == 0:
            self.root.ids.sb_taxa.disabled = True
        else:
            self.root.ids.sb_taxa.disabled = False

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
            self.root.ids.file_sl.height = 5

        # Enable selection buttons if file list is not empty
        if self.file_list:
            for i in self.root.ids.sb_file.children:
                i.disabled = False
                i.bind(on_release=self.select_bt)

        # Add a label at the end of the file list informing how many files are
        # currently selected out of the total files
        self.update_file_label()

        for infile in self.file_list:

            file_name = infile.split("/")[-1]
            # This prevents duplicate files from being added
            if file_name not in [x.id for x in self.root.ids.file_sl.children]:

                bt = ToggleButton(text=file_name, state="down", id=file_name,
                                  height=30, size_hint=(.8, None), shorten=True,
                                  shorten_from="right", halign="center")
                # Add button to storage for mouse over events
                self.mouse_over_bts["Files"].append(bt)
                # Setting horizontal text size for shortening
                bt.text_size[0] = bt.size[0] * 1.3
                # Binding functionality to toggle button
                bt.bind(on_release=self.toggle_selection)

                # Adds toggle button with file name
                self.root.ids.file_sl.add_widget(bt)

                # Set Information button and add the widget
                inf_bt = Button(text="?", size_hint=(.14, None),
                                height=30, id="%s?" % file_name, bold=True)
                self.root.ids.file_sl.add_widget(inf_bt)
                inf_bt.bind(on_release=self.popup_info)

                # Set remove button with event binded and add the widget
                x_bt = Button(text="X", size_hint=(.14, None),
                              height=30, id="%sX" % file_name,
                              background_color=(255, .9, .9, 1), bold=True)
                x_bt.bind(on_release=partial(self.check_action, self.remove_bt))
                self.root.ids.file_sl.add_widget(x_bt)

                # Updates the size of the grid layout according to the added
                # buttons
                self.root.ids.file_sl.height += 35

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

        # Enable selection buttons if taxa list is not empty
        if self.active_taxa_list:
            for i in self.root.ids.sb_taxa.children:
                i.disabled = False
                i.bind(on_release=self.select_bt)

        # Add a label at the end of the taxa list informing how many taxa are
        # currently selected out of the total taxa
        self.update_sp_label()

        for tx in sorted(self.active_taxa_list):

            # Prevents duplicate taxa from being entered
            if tx not in [x.id for x in self.root.ids.taxa_sl.children]:

                bt = ToggleButton(text=tx, state="down", id=tx,
                                  height=30, size_hint=(.8, None), shorten=True,
                                  shorten_from="right", halign="center")
                # Add button to storage for mouse over events
                self.mouse_over_bts["Taxa"].append(bt)
                # Setting horizontal text size for shortening
                bt.text_size[0] = bt.size[0] * 1.3
                # Binding functionality to toggle button
                bt.bind(on_release=self.toggle_selection)

                # Add toggle button with taxa name
                self.root.ids.taxa_sl.add_widget(bt)

                # Set Information button and add the widget
                inf_bt = Button(text="?", size_hint=(.14, None),
                                height=30, id="%s?" % tx, bold=True)
                self.root.ids.taxa_sl.add_widget(inf_bt)
                inf_bt.bind(on_release=self.popup_info)

                # Set remove button with event binded and add the widget
                x_bt = Button(text="X", size_hint=(.14, None),
                              height=30, id="%sX" % tx,
                              background_color=(255, .9, .9, 1), bold=True)
                self.root.ids.taxa_sl.add_widget(x_bt)
                x_bt.bind(on_release=partial(self.check_action, self.remove_bt))

                # Updates the size of the grid layout according to the added
                # button
                self.root.ids.taxa_sl.height += 35

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

                close_bt = Button(text="Close", size_hint_y=None, height=40)
                close_bt.bind(on_release=self.dismiss_popup)

                all_ds.add_widget(total_ds)
                all_ds.add_widget(active_ds)
                sv.add_widget(all_ds)
                content.add_widget(sv)
                content.add_widget(close_bt)

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
                                self.original_file_inf[file_name]["aln_format"]
                content.ids.seq_type.text = "%s" % \
                                self.original_file_inf[file_name]["seq_type"]
                content.ids.is_aln.text = "%s" % \
                                self.original_file_inf[file_name]["is_aln"]
                content.ids.seq_size.text = "%s" % \
                                self.active_file_inf[file_name]["aln_len"]
                content.ids.n_taxa.text = "%s" % \
                                self.active_file_inf[file_name]["n_taxa"]

                self.show_popup(title="File: %s" % value.id[:-1],
                                content=content, size=(350, 310))

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
                        # Update label
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

        # Updates labels
        self.update_sp_label()
        self.update_file_label()

    def dialog_taxagroup(self):
        """
        Creates the layout for the taxa group creation popup.
        """

        # Initializing instance for taxa group dialog
        content = TaxaGroupDialog(cancel=self.dismiss_popup)

        # Populate the gridlayout for all taxa
        for i in self.alignment_list.taxa_names:
            # Create togglebutton for each taxa
            bt = ToggleButton(text=i, size_hint_y=None, height=30)
            self.add_taxa_bt(bt, content.ids.all_grid)

        # Show dialog
        self.show_popup(title="Create taxa groups", content=content,
                        size=(700, 500))

    @staticmethod
    def add_taxa_bt(bt, wgt):
        """
        Method for addition of a button to a widget. This method was created
        for the automatic upated of the widgets height when moving buttons in
        the taxa group creation dialog
        :param bt: The button widget
        :param wgt: The sink widget
        """

        wgt.add_widget(bt)
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

    def taxagroup_move_taxa(self, source_wgt, sink_wgt, all_taxa):
        """
        Method that adds functionality to the addition/removal buttons (<<, <,
        >>, >) in the taxa group dialog.
        :param source_wgt: widget, the gridlayout from where the buttons will
        be moved
        :param sink_wgt: widget, the gridlayout to where buttons will be moved
        :param all_taxa: Boolean, if True its as if alsa taxa were selected to
        be moved
        """

        # In case all taxa are to be moved
        if all_taxa:
            # Ensures that only toggle buttons are moved
            for bt in source_wgt.children[::-1]:
                self.remove_taxa_bt(bt, source_wgt)
                bt.state = "normal"
                self.add_taxa_bt(bt, sink_wgt)
        else:
            for bt in source_wgt.children[::-1]:
                if bt.state == "down":
                    self.remove_taxa_bt(bt, source_wgt)
                    bt.state = "normal"
                    self.add_taxa_bt(bt, sink_wgt)

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
        close_bt = Button(text="Close", size_hint_y=.2)
        # Add functionality to the close button
        close_bt.bind(on_release=self.dismiss_popup)
        # Create gridlayout that will store the buttons with taxa names
        gl = GridLayout(cols=1, size_hint_y=None, height=0, spacing=5)

        # Add buttons with taxa names to the gridlayout and update the height
        # of the gridlayout
        for tx in self.taxa_groups[name_wgt.text]:
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
        Removes the taxa group button from the app list and taxa_groups
        attribute
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
        del self.taxa_groups[bt_idx]

    def save_taxa_group(self, source_wgt, name):
        """
        Adds a taxa group declared using the taxa group creator popup to the
        list of taxa groups in the side panel
        :param source_wgt, gridlayout of the selected taxa
        :param name: string, name of the group
        """

        # Make core changes by populating self.taxa_groups dictionary
        self.taxa_groups[name] = []

        for bt in source_wgt.children:
            self.taxa_groups[name].append(bt.text)

        # App changes by adding three buttons for the taxa group
        # Taxa button itself
        bt = Button(text=name, size_hint=(.8, None), height=40, id=name)
        bt.bind(on_release=self.taxagroups_show_taxa)
        # Removal button
        x_bt = Button(text="X", bold=True, size_hint=(.14, None),
                        height=40, id="%sX" % name,
                        background_color=(255, .9, .9, 1))
        x_bt.bind(on_release=partial(self.check_action, self.remove_taxa_group))

        # Add buttons to gridlayout
        for i in [bt, x_bt]:
            self.root.ids.group_grid.add_widget(i)

        # Add this options to the dropdown menu in the general options of the
        # Process screen
        if "empty_taxaset" in self.process_grid_wgt.ids.keys():
            self.process_grid_wgt.ids.dropdown_gl.remove_widget(
                self.process_grid_wgt.ids.empty_taxaset)
            del self.process_grid_wgt.ids["empty_taxaset"]

        # Create separator between dropdown items
        separator = Widget(size_hint_y=None, height=3)
        dd_bt = Button(text=name, size_hint_y=None, height=40,
                       background_normal="data/backgrounds/bt_process.png",
                       background_color=(1, 1, 1, .3))
        dd_bt.bind(on_release=lambda dd_bt:
                   self.process_grid_wgt.ids.dataset_dropdown.select(name))
        self.process_grid_wgt.ids.dataset_dropdown.add_widget(separator)
        self.process_grid_wgt.ids.dataset_dropdown.add_widget(dd_bt)

        # Update gridlayout height
        self.root.ids.group_grid.height += 40

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
                if self.secondary_options["%s_file" % op]:
                    add_node("%s (Save as new file)" % op,
                             self.main_nodes["proc_sec"])
                else:
                    add_node("%s" % op, self.main_nodes["proc_sec"])
            if self.main_nodes["proc_sec"].is_open is False:
                self.operation_tv.toggle_node(self.main_nodes["proc_sec"])
        else:
            self.main_nodes["proc_sec"].opacity = .2

        # Output file
        clear_nodes(self.main_nodes["main_file"])
        if self.output_file == "" and self.main_operations["conversion"]:
            add_node("[Dependent on input data]", self.main_nodes["main_file"])
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

        # Clear filters and haplotype name
        self.filter_settings = [25, 50]
        self.hap_prefix = "Hap"

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

        # Fade out of additional options
        if self.process_options.opacity == 1:
            self.toggle_process_options()

        # Changes in buttons with dynamic text
        self.process_grid_wgt.ids.conv_formatbt.text = "Fasta"
        self.process_grid_wgt.ids.conv.text = "Select..."

        # Turn switches off
        for switch in self.secondary_operations:
            self.process_options.ids[switch].active = False

        for switch in self.secondary_options:
            self.process_options.ids[switch].active = False

    ########################### PROCESS SCREEN #################################

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

    def show_popup(self, title, content, size_hint=(.9, .9), size=None):
        """
        General purpose method to create a popup widget
        :param title: string. Title of the popup
        :param content: widget object. The contents of the popup widget
        :param size_hint: tuple. Size hint for the widget
        :param size: tuple. The absolute size for the popup. If this argument is
        used, the size_hint will be ignored
        """

        # Ignore size_hint is absolute size is provided
        if size:
            self._popup = Popup(title=title, content=content, size=size,
                                size_hint=(None, None), auto_dismiss=False)
        else:
            self._popup = Popup(title=title, content=content,
                                size_hint=size_hint, auto_dismiss=False)
        self._popup.open()

    def save_file(self, path, filename, idx):
        """
        Adds functionality to the save button in the output file chooser. It
        gathers information on the specified path through filechooser, file
        name through textinput and the widget text when called.

        For now, only one main output file can be provided, so the its path
        is stored in a string attribute.

        :param path: string. complete path
        :param filename: string. file name only
        :param idx: sting. widget id
        """

        # Ensures that empty file names cannot be specified
        while filename != "":

            # Adds output file to storage
            self.output_file = join(path, filename)
            # Renames the output file button text
            self.process_grid_wgt.ids.conv.text = filename
            # Close popup
            self.dismiss_popup()
            # Breaks loop
            break

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
            if wgt.state == "down" and idx not in self.output_formats:
                self.output_formats.append(idx)
            elif wgt.state == "normal" and idx in self.output_formats:
                self.output_formats.remove(idx)

        self.dismiss_popup()

        # Updates the text in the select format button. In case only one format
        # is specified, the text will be that format; if multiple formats are
        # specified, the text will inform the number of selected formats; if no
        # format is specified, a no selected format text will appear
        if len(self.output_formats) == 1:
            self.process_grid_wgt.ids.conv_formatbt.font_size = 15
            self.process_grid_wgt.ids.conv_formatbt.text = \
                self.output_formats[0].title()

        elif len(self.output_formats) == 0:
            self.process_grid_wgt.ids.conv_formatbt.font_size = 14
            self.process_grid_wgt.ids.conv_formatbt.text = "No formats selected"

        else:
            self.process_grid_wgt.ids.conv_formatbt.font_size = 16
            self.process_grid_wgt.ids.conv_formatbt.text = "%s selected" % (
                len(self.output_formats))

        # Updates the Gcoder option depending on whether the nexus output format
        # is the only one selected
        if self.output_formats == ["nexus"]:
            self.process_options.ids.gcoder.disabled = False
        else:
            self.process_options.ids.gcoder.active = False
            self.process_options.ids.gcoder.disabled = True

    def save_filter(self, gap_val, mis_val):
        """
        Stores the information of the FilterDialog
        """

        self.filter_settings = [gap_val,
                                mis_val]

        self.dismiss_popup()

    def save_phylip_extra(self):

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
            self.process_options.ids.zorro.text = "Off"

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

    def dialog_filechooser(self, value):
        """
        Generates a file chooser popup for the user to select an output file
        """

        # Inherits the layout defined in the .kv file under <SaveDialog>
        content = SaveDialog(cancel=self.dismiss_popup)
        # Sets the home path as starting point
        content.ids.sd_filechooser.path = self.home_path

        # Save output file for conversion/concatenation purposes
        # Providing this operation will allow the filechooser widget to
        # know which output file is this
        content.ids.sd_filechooser.text = value

        self.show_popup(title="Choose output file", content=content)

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
        It handles common misktakes, such as using "," instead of "." for
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

    def dialog_text(self, title=""):
        """
        Generates a simple text dialog to capture text input
        """

        content = TextDialog(cancel=self.dismiss_popup)
        content.ids.txt_dlg.text = self.hap_prefix

        self.show_popup(title=title, content=content,
                        size=(200, 150))

    def dialog_warning(self, msg1, msg2):

        content = WarningDialog(cancel=self.dismiss_popup)
        content.ids.warning_l.text = "[b][size=18]%s[/size][/b]\n\n%s" % (msg1,
                                                                          msg2)

        self.show_popup(title="Warning!", content=content, size=(400, 300))

    def dialog_zorro(self):

        content = ZorroDialog(cancel=self.dismiss_popup)

        content.ids.zorro_switch.active = self.secondary_options["zorro"]
        content.ids.zorro_txt.text = self.zorro_suffix

        self.show_popup(title="ZORRO support", content=content,
                        size=(350, 200))

    def save_hap_prefix(self, text_wgt):
        """
        Saves the specified string suffix for haplotypes when collapsing
        :param text_wgt. The widget of the text input to retrieve its text
        property
        """

        self.hap_prefix = text_wgt.text

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
        for k in self.main_operations:
            if op == k:
                self.main_operations[op] = True
            else:
                self.main_operations[k] = False

        # Disables output file button and other conversion/concatenation
        # specific buttons
        if op == "conversion":
            self.process_grid_wgt.ids.conv.disabled = True
            self.process_options.ids.zorro.disabled = True
        else:
            self.process_grid_wgt.ids.conv.disabled = False
            self.process_options.ids.zorro.disabled = False

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
            self.process_grid_wgt.height = self.process_height + (55 * len(
                self.process_options.ids.main_grid.children))

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

    def load_files(self, files):
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

        aln_list = AlignmentList(files)

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

        return aln_list.bad_alignments

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

        if self.alignment_list:
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
                if aln.model:
                    file_inf[file_name]["is_aln"] += " (Concatenated)"

                # Get length of largest sequence if not aligned, or alignment
                # length
                file_inf[file_name]["aln_len"] = aln.locus_length

        return file_inf

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

        # Setting the alignment to use. A deepcopy of the active alignment list
        # is used because it may be possible to do changes in the taxa data set
        # of the AlignmentList object, which should not change the original
        # self.active_alignment_list. This is because when taxa are removed from
        # the alignment list, there is no way to return those taxa to the
        # object
        aln_object = deepcopy(self.active_alignment_list)

        # Update active taxa set of the alignment object
        aln_object = self.update_active_taxaset(aln_object)

        # Concatenation
        if self.main_operations["concatenation"]:
            aln_object = aln_object.concatenate()

            # Concatenation of ZORRO files
            if self.secondary_options["zorro"]:
                zorro_data = data.Zorro(aln_object, self.zorro_suffix)
                zorro_data.write_to_file(self.output_file)

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
            for name, obj in write_aln.items():
                obj.write_to_file(self.output_formats, name,
                                interleave=self.secondary_options["interleave"],
                                partition_file=self.create_partfile)
        else:
            for name, obj in write_aln.items():
                name = name.replace(self.output_file, "")
                obj.write_to_file(self.output_formats, output_suffix=name,
                                interleave=self.secondary_options["interleave"],
                                partition_file=self.create_partfile)


if __name__ == '__main__':
    TriFusionApp().run()