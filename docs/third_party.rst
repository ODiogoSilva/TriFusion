Setup of USEARCH
================

The `USEARCH <http://www.drive5.com/usearch/>`_ software is required to
perform the **Orthology** search operation and to export ortholog
groups into nucleotide sequences. However, due to licensing issues, USEARCH
cannot be bundled with Triusion, so it requires some user intervention to
setup. But don't fret! Everything can be up and running with just a few
simple steps. Moreover, after the initial setup, TriFusion will store
the USEARCH executable internally and use it for all subsequent sessions.

    - Step 1: Download the USEARCH executable for your corresponding operating
      system `here <http://www.drive5.com/usearch/download.html>`_.

    - Step 2: If USEARCH is not reachable by TriFusion, you will see a warning
      like this when you navigate to ``Orthology -> Show additional options
      -> USEARCH``:

      .. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/images/usearch_warning_box.png
          :alt: pic
          :align: center

      Click the ``Fix it`` button, and then the ``Search USEARCH executable``
      button.

      .. image:: https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/images/usearch_file_search.png
          :alt: pic
          :align: center

    - Step 3: Search for the executable you have downloaded in **Step 1** and
      click the ``Save`` button.

And that's it. When a valid USEARCH executable is provided, the previous
warning should be replaced with a green box saying "*USEARCH is installed and
reachable*". You are good to go!
