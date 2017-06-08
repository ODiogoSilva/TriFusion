"""
Welcome to the TriFusion API reference guide. This reference guide details
the sub-packages and modules used for each component of TriFusion.

This guide is intended for power users and developers that wish to modify
the source and/or contribute with new code. A brief description of each
main component of TriFusion follows below with the relevant modules and
classes specified. These descriptions can be used to send you on your
way to the sections of interest, where a more detailed documentation
is available.

What is TriFusion
=================

TriFusion is a GUI and command line application designed to streamline the
gathering, processing and visualization of phylogenomic data. It is broadly
divided in these three modules.

Orthology
---------

Provides a pipeline for running OrthoMCL, with the code ported to python
and SQLite instead of the original perl and MySQL. OrthoMCL is the most
popular ortholog detection pipeline, and TriFusion offers an  easy and
intuitive way of running the pipeline while providing the complete
range of options available.

In addition to the search pipeline, TriFusion allows the filtering and visual
exploration of the resulting ortholog groups. In the end, orthologs can
be exported as protein or DNA sequences.

Process
-------

At its core, the Process module is a conversion and concatenation tool that
handles very large sequence alignment data sets. It reads and exports
alignments into several popular formats used in phylogenetics and population
genetics. In addition to these main operations, TriFusion offers a wide
array of manipulations that can be performed on alignment data, such as
filtering, collapsing, creating consensus, etc.

Statistics
----------

Generates a wide array of graphical visualizations and statistical
analyses of alignment data sets.

How can TriFusion be used
=========================

TriFusion can be used as a:

    - Desktop application with graphical interface (TriFusion).
    - Command line application (orthomcl_pipeline, TriSeq and TriStats).
    - Library of high performance classes to parse, modify, export and
      plot alignment data.

Components of TriFusion
=======================

The TriFusion package is the result of multiple modular components that
work in combination to produce the main application. This modularity
means that the graphical interface is separated from the multiple
backends that power its features, and each backend works independently
of each other. This greatly facilitates changing the existing objects
or creating new ones for specific modules without having to worry about
other aspects of the package.

Here is a brief overview of the GUI and backend components.

TriFusion GUI
-------------

The core graphical interface of TriFusion is controlled by two main files:

  - :mod:`trifusion.app`: Contains the :class:`~trifusion.app.TriFusionApp`
    class, with most of the methods and attributes responsible for
    the interactivity and event handling of the TriFusion application.

  - :mod:`trifusion.trifusion.kv`: Contains graphical instructions for the
    main window and the definition of the majority of the Widgets in the
    kivy language.

.. warning:: API documentation of :mod:`trifusion.app` is still under progress.

Main screens
~~~~~~~~~~~~

The graphical instructions in kivy language for the 8 main screens of
TriFusion are defined in the `trifusion/data/screens` directory, each screen
with its own `.kv` file. These files contain only the definition of
the root widgets of each screen. Other Widgets that are later added to
the screen via some method
should be defined in the main `trifusion.trifusion.kv` file.
The initial setup of these screens in performed in the
:func:`~trifusion.app.TriFusionApp.build` method.

Custom widgets
~~~~~~~~~~~~~~

Custom widgets can be created using the `kivy` toolkit. `Kivy` provides
a convenient way of defining graphical instructions to build the Widget's
layout using the kivy language instead of directly with python. Therefore,
the layout of new widgets can be defined in the `trifusion.trifusion.kv` file.
If they need to be used by the python code, they also have to be defined
as a new class in the :mod:`trifusion.data.resources.custom_widgets` module
and imported in the module where they will be used.
This class can be empty (e.g.
:func:`~trifusion.data.resources.custom_widgets.TableCell`) or it can harbor
custom attributes and methods that are useful for that widget
(e.g. :func:`~trifusion.data.resources.custom_widgets.FileChooserL`).

.. warning:: API documentation of
             :mod:`trifusion.data.resources.custom_widgets` is still under
             progress.

Icons and backgrounds
~~~~~~~~~~~~~~~~~~~~~

Icons and background images are stored in the `trifusion/data/backgrounds`
directory. In either python or kivy files, these backgrounds can be
reference with the path relative to the package root. For example, in any
python file, they can be referenced like::

    bn = "data/backgrounds/check_ok.png"

Or in any kivy file::

    background_normal: "data/backgrounds/bt_process.png"

Running tasks in the background
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Time consuming operations are executed in worker threads separated from the
main GUI thread. These background tasks are defined in
:mod:`trifusion.data.resources.background_tasks`. The module documentation
provides more information on how to setup background tasks in TriFusion.

In-App help
~~~~~~~~~~~

Help buttons are spread throughout TriFusion. The help information is defined
in multiple dictionary objects in :mod:`trifusion.data.resources.info_data`.

Orthology backend
-----------------

The orthology search pipeline is defined in
:mod:`trifusion.orthomcl_pipeline`. This module can be used as a CLI program
to execute the pipeline and is also used by TriFusion as a library.

The handling and exploration of ortholog group files, the output of the
orthology search operation, is made in the
:mod:`trifusion.ortho.OrthomclToolbox` module, in particular using the
:class:`~trifusion.ortho.OrthomclToolbox.MultiGroupsLight` and
:class:`~trifusion.ortho.OrthomclToolbox.GroupLight` classes.

.. warning:: API documentation for :mod:`trifusion.orthomcl_pipeline` and
             :mod:`trifusion.ortho.OrthomclToolbox` is still in progress.

Process backend
---------------

The main functionality of the Process module (and the TriSeq CLI) program
are provided by the modules in the :mod:`trifusion.process` sub package.
Classes that handle alignment data are defined in
:mod:`trifusion.process.sequence`, while data set partitions are handled
in the :mod:`trifusion.process.data` module.

Statistics backend
------------------

The generation of plot data in the Statistics screen or by the TriStats
CLI program is a joint effort between the
:class:`~trifusion.process.sequence.AlignmentList` class and the plotting
functions defined in the :mod:`trifusion.base.plotter` module. Briefly,
the methods of the :class:`~trifusion.process.sequence.AlignmentList` class
are responsible for generating the data and plotting instructions, while
the functions in the :mod:`trifusion.base.plotter` module receive that
information and generate the plot.
"""

__version__ = "0.5.9"
__build__ = "080617"
__author__ = "Diogo N. Silva"
__copyright__ = "Diogo N. Silva"
__credits__ = ["Diogo N. Silva", "Tiago F. Jesus", "Fernando Alves"]
__license__ = "GPL3"
__maintainer__ = "Diogo N. Silva"
__email__ = "o.diogosilva@gmail.com"
__status__ = "4 - Beta"
