Data set groups
===============

.. note::

    **Data availability for this tutorial**: the medium sized data
    set of 614 genes and 48 taxa that will be used can be downloaded
    `here <https://github.com/ODiogoSilva/TriFusion-tutorials/raw/master/tutorials/Datasets/Process/medium_protein_dataset/medium_protein_dataset.zip>`_.

What are active data sets
-------------------------

Most operations in TriFusion can be applied to either the **total data set**
(all files and taxa currently loaded) or to custom made data sets,
named **active data sets**. When a custom data set is specified, operations
will be applied only on the active files and/or taxa and ignore all others.
These active data sets can defined in TriFusion in several ways and serve
to quickly apply different operations on different sets of files/taxa.

.. figure:: https://github.com/ODiogoSilva/TriFusion-tutorials/raw/master/tutorials/images/dataset_active_example.png
    :alt: pic
    :align: center

    Example of custom active file (left) and taxa (right) data sets.

How to define active data sets
------------------------------

Active data sets can be created/modified in two main ways:

    - `Toggle file/taxa buttons in side panel`_
        - `Mouse click toggling`_
        - `Import selection from file`_
    - `Create data set groups`_
        - `Manual creation in TriFusion`_
        - `Group creation from file`_


Toggle file/taxa buttons in side panel
--------------------------------------

Mouse click toggling
^^^^^^^^^^^^^^^^^^^^

By default, when data is loaded into TriFusion **all** files/taxa are active.
Therefore, the total and active data sets are the same. The quickest way
to modify the active data set is by navigating to ``Menu -> Open/View Data``
and toggle the corresponding file/taxa buttons. ``Shift + click`` is also
supported to select multiple contiguous files/taxa.

.. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/gifs/process_tutorial2_mouse_toggle.gif
    :alt: pic

Active files/taxa will appear with a blue background, while inactive
buttons will have no background. A label below the button list displays
how many files/taxa are currently active.

Import selection from file
^^^^^^^^^^^^^^^^^^^^^^^^^^

When dealing with a larger number of files/taxa it may be more convenient to
provide the active data set through a text file. This should be a simple
text file containing the names of the desired files/taxa in each line. You
can create it yourself, or download an example `from
here <https://github.com/ODiogoSilva/TriFusion-tutorials/raw/master/tutorials/Datasets/Process/medium_protein_dataset/taxa_list.txt>`_.

::

    # Example of a text file for taxa selection in TriFusion
    Agaricus_bisporus
    Botrytis_cinerea
    Coniophora_puteana

::

    # Example of a text file for file selection in TriFusion (note the extension)
    BasidioOnly2585_linsi_missingFilter_concPrep.fasta
    BasidioOnly2685_linsi_missingFilter_concPrep.fasta
    BasidioOnly2686_linsi_missingFilter_concPrep.fasta

Open the ``Menu -> Open/View Data`` side panel and click on the ``+`` button
at the bottom of either the *Files* or *Taxa* tabs. This will open a sub-menu
with several options, one of which is ``Select file/taxa names from .txt``.
Clicking this button will open a file browser where you can provide the
file containing the file/taxa names. Once you select the text file, the
the active file/taxa names will update.

.. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/images/dataset_toggle_fromfile.png
    :alt: pic

.. warning::

    After loading the file, **ONLY** the specified items will become active,
    regardless of the previous active data set. Names that do not match any of
    the files/taxa present in TriFusion will be ignored.

.. note::

    You can also save any active files/taxa on the side panel to a text file
    by clicking the ``Export selected file/taxa names to .txt``.

Create data set groups
----------------------

When the workflow requires the execution of operations to multiple
taxa/files data sets, it is more convenient to define all data set groups and
then use the dropdown menus (see `How to apply data set groups`_ below) to
select the desired active data set. Data set groups can be defined in
TriFusion by navigating to ``Menu > Dataset Groups``.

.. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/images/dataset_creationpanel.png
    :alt: pic

File and taxa groups are sorted into two tabs, like in the ``Open/View Data``
panel, and clicking the ``Set new file/taxa group`` button will start the
creation of the group.

.. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/images/dataset_triage.png
    :alt: pic

Here you can choose to create the data set group either manually in TriFusion,
or by providing the names of the files/taxa in a text file.

Manual creation in TriFusion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::
    This option is discouraged for larger data sets (>500 items). In these
    cases, it is recommended to use the `Group creation from file`_ method.

The creation of groups is the same for both files and taxa. In this tutorial,
we will create a taxa group by clicking in the *Taxa* tab and then the
``Set new taxa group`` button at the bottom of the side panel.
Here, groups can be created by selecting the desired taxa from the
*All taxa* column and using the arrow buttons to move them to the
*Selected taxa* column. Once the group is complete, give it a unique name
and the group is ready to be defined. If you wish to create multiple groups
in one sitting, click the ``Apply`` button to create the group but remain
in the dialog.

.. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/gifs/process_tutorial2_manual_selection.gif
    :alt: pic

Any previously created group will be listed under the *Created groups*
column. These can be selected to move their corresponding taxa to the
*Selected taxa* column and continue a new group definition from there.

Group creation from file
^^^^^^^^^^^^^^^^^^^^^^^^

Here, we only have to provide a text file with the names of the files/taxa
we wish to select for the group. The text file is the same as the one
described in the `Import selection from file`_ example.

::

    # Example of a text file for taxa selection in TriFusion
    Agaricus_bisporus
    Botrytis_cinerea
    Coniophora_puteana

After providing the file with the group names, specify a unique name of the
new data set group, and that's it!

.. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/gifs/process_tutorial2_file_selection.gif
    :alt: pic

How to apply data set groups
----------------------------

Now that we know how to create active data set groups, the final step is how
can they be specified.

Orthology
^^^^^^^^^

When using the **Orthology** module, only the active proteome files are used
for the Orthology search operation.

Process and Statistics
^^^^^^^^^^^^^^^^^^^^^^

For both **Process** and **Statistics** modules, the *active* data set
is selected by default (that is, the file/taxa buttons active in the side
panel). You can change to the *total* data set or to any user made data
set by clicking the group's name in the corresponding dropdown menu.

    Dropdown menu in the **Process** screen:

    .. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/images/process_data_set_selection.png
        :alt: pic

    Dropdown menu in the **Statistics** screen:

    .. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/images/statistics_data_set_selection.png
        :alt: pic