"""
Introduction to TriFusion's process module
==========================================

The `process` subpackage is the main backend of TriFusion's Process and
Statistics modules and of TriSeq and TriStats CLI programs. The most
important classes are defined in the :mod:`~trifusion.process.sequence` module:
:class:`~trifusion.process.sequence.Alignment` and
:class:`~trifusion.process.sequence.AlignmentList`.

What it does
------------

The `process` module contains the classes and functions responsible for
parsing, modifying, writing and plotting alignment data.

Submodules description
----------------------

:mod:`~trifusion.process.base`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Contains several methods and function that are inherited or used by
:class:`~trifusion.process.sequence.Alignment` and
:class:`~trifusion.process.sequence.AlignmentList` objects, as well as by
the TriSeq and TriStats CLI programs.

:mod:`~trifusion.process.data`
~~~~
Contains the :class:`~trifusion.process.data.Partitions`  class, used by
:class:`~trifusion.process.sequence.Alignment` and
:class:`~trifusion.process.sequence.AlignmentList`
classes to handle partitions in the alignments.

:mod:`~trifusion.process.error_handling`
~~~~~~~~~~~~~~
Contains custom made Exception sub-classes.

:mod:`~trifusion.process.sequence`
~~~~~~~~
Contains the :class:`~trifusion.process.sequence.Alignment`  and
:class:`~trifusion.process.sequence.AlignmentList` classes, responsible
for the majority of the heavy lifting when dealing with alignment files. See
the module's documentation for further details.
"""