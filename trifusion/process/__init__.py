"""
Introduction to TriFusion's process module
==========================================

The `process` sub-package is the main backend of TriFusion's Process and
Statistics modules and of TriSeq and TriStats CLI programs. The most
important classes are defined in the `sequence` module: `Alignment` and
`AlignmentList`.

What it does
------------

The `process` module contains the classes and functions responsible for
parsing, modifying, writing and plotting alignment data.

Sub modules
-----------

base
~~~~
Contains several methods and function that are inherited or used by
`Alignment` and `AlignmentList` objects, as well as by the TriSeq and
TriStats CLI programs.

data
~~~~
Contains the `Partition` class, used by `Alignment` and `AlignmentList`
classes to handle partitions in the alignments.

error_handling
~~~~~~~~~~~~~~
Contains custom made Exception sub-classes.

sequence
~~~~~~~~
Contains the `Alignment` and `AlignmentList` classes, responsible for the
majority of the heavy lifting when dealing with alignment files. See the
module's documentation for further details.
"""