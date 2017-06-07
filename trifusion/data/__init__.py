"""
The `data` sub-package is where the background data, kivy screens definition,
bundled binary and .dll files, background tasks and other general data
is contained.

`backgrounds`
-------------

Stores all the backgrounds and icons of TriFusion.

`resources`
-----------

Stores several resources used by TriFusion, including custom widgets,
definition of background (worker) tasks, information data for in-app help,
and auxiliary statistics data structures:

    - :mod:`~trifusion.data.resources.background_tasks`: Includes all
      TriFusion tasks that are executed on the background.
    - :mod:`~trifusion.data.resources.info_data`: Contains dictionary
      objects with the in-app help text.
    - :mod:`~trifusion.data.resources.custom_widgets`: Contains class
      definition of custom widgets used in TriFusion.
    - :mod:`~trifusion.data.resources.stats`: Dictionary with information
      for the Statistics FloatLayout with the plot type switch.

It also contains the bundled MCL binary for all supported OS and the usearch
.dll file, which are required by TriFusion's Orthology.

`screens`
---------

The definition os the several TriFusion screens in kivy language.
"""