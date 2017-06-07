"""
The `resources` submodule contains:

`background_tasks`
------------------

The `background_tasks` module is where the functions executed in the worker
thread of TriFusion are defined. All time consuming functions that need to
be executed in the background to prevent the app from freezing should
be defined in this module and then called from the main TriFusion's `app`
file.

'custom_widgets`
----------------

Custom kivy widgets that need to be used in TriFusion's main `app` file
are defined here.

`info_data`
-----------

Contains the information storage for several parts of TriFusion. All
information from the in-app help is defined here in the `informative_storage`
dictionary.

`stats`
-------

An auxiliary file to TriFusion's Statistics screen. It defines a dictionary
object with the information of the rapid switch float layout that is
added to most Statistics plots.

"""