TriFusion
---------

Making life easier for phylogenomic data gathering, processing and visualization
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Website: http://odiogosilva.github.io/TriFusion/

| |Build Status|
| |Code Health|
| |codecov|
| |PyPI|
| |PyPI|
| |AUR|

What is TriFusion?
------------------

TriFusion is a modern GUI and command line application designed to make
the life of anyone with **proteome** and/or **alignment sequence data**
easier and more pleasurable. Regardless of your experience in
bioinformatics, TriFusion is easy to use and offers a wide array of
powerfull features to help you deal with your data. At the same time, it
was developed to handle the enormous amount of data that is generated
nowadays.

TriFusion is an open source, cross-platform application written in
`Python 2.7 <https://www.python.org/>`__ and using the
`Kivy <https://github.com/kivy/kivy>`__ framework to build the graphical
interface.

What can TriFusion do for you?
------------------------------

Here is an overview of what it can do for you across its three main
modules.

Orthology - Search and explore orthologs across proteomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  **Searches for ortholog sequences** across multiple species.
-  Filters ortholog sequences according to the **gene copy number**
   and/or **number of taxa** present.
-  **Graphical visualization** of ortholog data.
-  Exports your orthologs as **protein or nucleotide sequences**.

`Find out more <https://odiogosilva.github.io/TriFusion/#featurette>`__

Process - Blazing fast processing of alignment files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  **Conversion** or **concatenation** of alignment files into several
   popular formats (`check supported
   formats <https://github.com/ODiogoSilva/TriFusion/wiki/Supported-Process-formats>`__).
-  **Collapse** identical sequences into the same haplotype.
-  Create **consensus** sequences for each alignment with several
   options on how to handle sequence variation.
-  **Filter** either alignments (according to whether they contain or
   exclude certain taxa, to a minimum proportion of taxa, and/or
   variable sites) or alignment columns (according to codon position,
   missing data and gaps).
-  **Code indel patterns** of your alignments into a binary matrix that
   is appended to the alignment.
-  **Revert concatenated alignments** or export sub-regions into
   individual files
-  Set **gene and codon partitions** as well as **substitution models**
   (Nexus format)
-  Create **file/taxa groups** to quickly perform operations on
   different sets of data.
-  It's **fast** and **memory efficient**. Converting 3,093 files with
   376 taxa can be performed in just 30 seconds and using less than 90Mb
   of RAM memory (`check the benchmarks
   table <https://github.com/ODiogoSilva/TriFusion/wiki/Benchmarks>`__).

`Find out more <https://odiogosilva.github.io/TriFusion/#featurette>`__

Statistics - Effortless visual exploration of your data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Provides instant information on overall and per gene **summary
   statistcs**.
-  TriFusion offers **dozens of graphical and statistical options** to
   explore your data:

   -  General information plots.
   -  Polymorphism and sequence variation plots.
   -  Missing data plots.
   -  Outlier plots.

-  Publication ready figures

`Find out more <https://odiogosilva.github.io/TriFusion/#featurette>`__

Installation
------------

Executables binaries
~~~~~~~~~~~~~~~~~~~~

The latest stable release of TriFusion can be installed as a standalone
application using one of the following installers. This **only includes
the GUI component** of TriFusion. If you also want the command line
version, see `Installation from source <#installation-from-source>`__.

Linux
^^^^^

-  Debian package based (`See
   list <https://en.wikipedia.org/wiki/Category:Debian-based_distributions>`__):
   `TriFusion-0.5.0.deb <https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0.deb>`__

-  RPM package based (`See
   list <https://en.wikipedia.org/wiki/Category:RPM-based_Linux_distributions>`__):
   `TriFusion-0.5.0.rpm <https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0.rpm>`__

-  ArchLinux/Manjaro (`See
   list <https://wiki.archlinux.org/index.php/Arch_based_distributions>`__):
   `TriFusion-0.5.0.tar.xz <https://aur.archlinux.org/packages/trifusion-bin/>`__
   is available on AUR.

MacOS
^^^^^

-  `TriFusion-0.5.0.app.zip <https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0-MacOS.app.zip>`__

Windows
^^^^^^^

-  `TriFusion-0.5.0 64bit
   installer <https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0-Win64.msi>`__
-  `TriFusion-0.5.0 32bit
   installer <https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0-Win32.msi>`__

Note for Windows 8.x and 10 users:
''''''''''''''''''''''''''''''''''

Executing the TriFusion installer may generate a warning from
SmartScreen. To continue with the installation, click the "More info"
label and then "Run anyway".

Installation from source
~~~~~~~~~~~~~~~~~~~~~~~~

TriFusion is on `PyPi <https://pypi.python.org/pypi/trifusion/>`__ and
can be easily installed with ``pip``. If you want to use the only the
command line version run this command on any operating system with
``pip``.

::

    pip install trifusion --user

For additional installation information or if you want to install the
complete TriFusion package, `follow these instructions according to your
operating
system <https://github.com/ODiogoSilva/TriFusion/wiki/Install-from-source>`__.

--------------

If you are unconvinced that a terminal version would be useful/pratical,
check out how easy and fast it is to use TriFusion to process 614 Fasta
alignments into phylip and nexus output formats :-):

How to use
----------

Tutorials on how to use TriFusion for its many tasks can be perused
`here <http://odiogosilva.github.io/TriFusion/#tutorials>`__.

Documentation
-------------

You can download TriFusion User Guide
`here <https://github.com/ODiogoSilva/TriFusion/raw/master/docs/manual.pdf>`__.

Citation
--------

When using OrthoMCL to find ortholog clusters, please cite the original
software:

Fischer, S., Brunk, B. P., Chen, F., Gao, X., Harb, O. S., Iodice, J.
B., Shanmugam, D., Roos, D. S. and Stoeckert, C. J. Using OrthoMCL to
Assign Proteins to OrthoMCL-DB Groups or to Cluster Proteomes Into New
Ortholog Groups Current Protocols in Bioinformatics. 2011
35:6.12.1–6.12.19.

We're working on a manuscript for TriFusion now.

.. |Build Status| image:: https://travis-ci.org/ODiogoSilva/TriFusion.svg?branch=master
   :target: https://travis-ci.org/ODiogoSilva/TriFusion
.. |Code Health| image:: https://landscape.io/github/ODiogoSilva/TriFusion/master/landscape.svg?style=flat
   :target: https://landscape.io/github/ODiogoSilva/TriFusion/master
.. |codecov| image:: https://codecov.io/gh/ODiogoSilva/TriFusion/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/ODiogoSilva/TriFusion
.. |PyPI| image:: https://img.shields.io/pypi/pyversions/trifusion.svg
   :target: https://pypi.python.org/pypi/trifusion
.. |PyPI| image:: https://img.shields.io/pypi/v/trifusion.svg
   :target: https://pypi.python.org/pypi/trifusion
.. |AUR| image:: https://img.shields.io/aur/version/trifusion.svg
   :target: https://aur.archlinux.org/packages/trifusion/
