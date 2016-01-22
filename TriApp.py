#!/usr/bin/python2

__author__ = 'diogo'

if __name__ == "__main__":

    from kivy.config import Config

    Config.set("kivy", "log_level", "warning")
    Config.set("kivy", "desktop", 1)
    Config.set("kivy", "exit_on_escape", 0)
    Config.set("graphics", "resizable", 1)
    Config.set("graphics", "fullscreen", 0)
    Config.set("graphics", "height", 700)
    Config.set("graphics", "width", 1000)
    Config.set("input", "mouse", "mouse, disable_multitouch")

    from kivy.base import EventLoop
    EventLoop.ensure_window()

    # Kivy imports
    from kivy.app import App
    from kivy.animation import Animation
    from kivy.uix.image import Image
    from kivy.uix.widget import Widget
    from kivy.uix.checkbox import CheckBox
    from kivy.lang import Builder
    from kivy.properties import ListProperty, DictProperty
    from kivy.clock import Clock
    from kivy.uix.treeview import TreeView, TreeViewLabel

    # Main program imports
    from ortho import protein2dna
    from process.base import Base
    from data.resources.info_data import informative_storage
    from data.resources.background_tasks import *
    from data.resources.custom_widgets import *
    from base.plotter import *

    # Other imports
    import os
    from os.path import dirname, join, exists, expanduser, basename
    from os import sep
    from collections import OrderedDict
    from copy import deepcopy
    from functools import partial
    import matplotlib.patches as patches
    import psutil
    import pickle
    import multiprocessing
    import time
    import re
    import sys
    import logging
    import shutil
    import urllib
    import matplotlib


    from TriFusion import TriFusionApp

    TriFusionApp().run()
