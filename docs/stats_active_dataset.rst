Update the active data set
==========================

.. note::

    **Data availability for this tutorial**: the medium sized data
    set of 614 genes and 48 taxa that will be used can be downloaded
    `here <https://github.com/ODiogoSilva/TriFusion-tutorials/raw/master/tutorials/Datasets/Process/medium_protein_dataset/medium_protein_dataset.zip>`_.

Data exploration analyses
-------------------------

The analyses in the **Statistics** module are not limited to the
**total data** set loaded into TriFusion. You can
**modify the active file/taxa data sets** or create data set groups
in TriFusion
(see tutorial :doc:`dataset_groups`), and then select them in the bottom
of the Statistics side panel.

.. image:: https://github.com/ODiogoSilva/TriFusion-tutorials/raw/master/tutorials/gifs/stats_tutorial1_change_active.gif
    :alt: pic

Following the guidelines in the
:doc:`dataset_groups` tutorial, we created a taxa group of 12 elements that
contains taxa whose name starts with an "A", "B" or "C", named **A_to_C**.
To change the taxa data set to the newly define group, click in the drop down
menu for the taxa data set and select the **A_to_C** option.

.. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/images/stats_selecting_taxa_group.png
    :alt: pic

Now, all selected analyses will use this set of 12 taxa instead of the full
48 taxa data set. If you want to update the currently displayed analyses,
click the **refresh** button next to the data set selection drop down menus.

.. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/images/stats_subset.png
    :alt: pic

Summary statistics
------------------

It is also possible to change the active data set when visualizing the
summary statistics of your data set and it can be particularly useful. For
example, if you suspect that a group of taxa or alignment may be responsible
for a particular large share of variability of missing data, you could
create

