"""The `sequence` module of TriFusion contains the main classes handling
alignment sequence data. These are :class:`.Alignment` and
:class:`.AlignmentList`. Here follows a brief explanation of how these
classes work and how to deal with the sqlite database.

`Alignment` class
-----------------

The :class:`.Alignment` class is the main interface with single alignment
files. It can be viewed as the building block of an :class:`.AlignmentList`
object, which can have one or more :class:`.Alignment` objects. It contains
all methods and attributes that pertain to a given alignment and are used to
retrieve information or modify it. **However, it is NOT meant to be used
independently**, but rather within the context of an :class:`.AlignmentList`
object. The data from each alignment is stored in a single sqlite database
during the execution of TriFusion or the TriSeq/TriStats CLI programs. The
connection to this database is automatically handled by
:class:`.AlignmentList` for all :class:`.Alignment` objects included in it.
In this way, we can use the :class:`.AlignmentList` class to handle the
setup of the sqlite3 database, and focus on single alignment data handling
in general in this class.

The main types of methods defined in this class are:

Parsers
~~~~~~~

Parsing methods are defined for each format with `_read_<format>`:

    - :meth:`~.Alignment._read_phylip`: Parses phylip format.
    - :meth:`~.Alignment._read_fasta`: Parses fasta format.
    - :meth:`~.Alignment._read_nexus`: Parses nexus format.
    - :meth:`~.Alignment._read_loci`: Parses pyRAD/ipyrad loci format.
    - :meth:`~.Alignment._read_stockholm`: Parses stockholm format.

They are always called from the :meth:`~.Alignment.read_alignment` method,
and not directly. When an :class:`.Alignment` object is instantiated with a
path to an alignment file, it automatically detects the format of the
alignment and keeps that information on the :attr:`~.Alignment.input_format`
attribute. The :meth:`~.Alignment.read_alignment` method then calls the
parsing method corresponding to that format. That information is stored in a
dictionary::

    parsing_methods = {
        "phylip": self._read_phylip,
        "fasta": self._read_fasta,
        "loci": self._read_loci,
        "nexus": self._read_nexus,
        "stockholm": self._read_stockholm
    }

    # Call the appropriate method
    parsing_methods[self.input_format]()

Each format has its own parsing method (which can be modified directly). To
add a new format, it is necessary to add it to the automatic format
recognition in :meth:`~trifusion.process.base.Base.autofinder`. Then,
create the new parsing method, using the same `_read_<format>` notation and
add it to the `parsing_methods` dictionary in
:meth:`~.Alignment.read_alignment`.

New parsers must insert alignment data into a table in the sqlite database.
This table is automatically created when the :class:`.Alignment` object
is instantiated, and its name is stored in the :attr:`~.Alignment.table_name`
attribute (see :meth:`~.Alignment._create_table`)::

    cur.execute("CREATE TABLE [{}]("
                "txId INT,"
                "taxon TEXT,"
                "seq TEXT)".format(table_name))

To insert data into the database, a taxon id (`txId`), taxon name (`taxon`)
and sequence (`seq`) must be provided. For example::

    cur.execute("INSERT INTO [{}] VALUES (?, ?, ?)".format(self.table_name),
                (0, "spa", "AAA"))

At the end of the parsing, these attributes should be set:

    - :attr:`~.Alignment.taxa_list`.
    - :attr:`~.Alignment.taxa_idx`.
    - :attr:`~.Alignment.locus_length`.
    - :attr:`~.Alignment.partitions`.

Data fetching
~~~~~~~~~~~~~

To facilitate fetching alignment data from the database, several generators
and data retrieval methods are defined:

    - :meth:`~.Alignment.iter_columns`: Iterates over the columns of the
      alignment.
    - :meth:`~.Alignment.iter_columns_uniq`: Iterates over the columns of
      the alignment but yields only unique characters.
    - :meth:`~.Alignment.iter_sequences`: Iterates over each sequence in the
      alignment.
    - :meth:`~.Alignment.iter_alignment`: Iterates over both taxon name
      and sequence in the alignment.
    - :meth:`~.Alignment.get_sequence`: Return the sequence from a particular
      taxon.

These should always use the
:func:`.setup_intable` decorator and defined with the `table_suffix`, `table_name`
and `table` arguments. For instance, the :meth:`.Alignment.iter_sequences` method is a
generator that allows the iteration over the sequences in the alignment
and is defined as::

    @setup_intable
    def iter_sequences(self, table_suffix="", table_name=None, table=None):

When calling these methods, only the `table_suffix` and `table_name` have to
be provided. In fact, the value provided to the `table` argument at calling
time is ignored. The decorator will check the values of both `table_suffix`
and `table_name` and evaluate the database table that will be used. This
final table name will then be provided as the `table` argument value within
the decorator. In this way, these methods can be called like::

    for seq in self.iter_sequences(table_suffix="_collapse"):
        # Do stuff

In this case, the :func:`.setup_intable` decorator will append the
`table_suffix` to the name of the original alignment table. If
:attr:`.Alignment.table_name`="main"`, then the final table name in this
case will be "main_collapse".

Alternatively, we can use `table_name`::

    for seq in self.iter_sequences(table_name="main_2"):
        # Do stuff

In this case, the final table name will be "main_2".

If the final table name does not exist, the method falls back to the original
table name defined by the `table_name` attribute.

Alignment modifiers
~~~~~~~~~~~~~~~~~~~

Methods that perform modifications to the alignment are also defined here.
These include:

    - :meth:`~.Alignment.collapse`: Collapses alignment into unique sequences.
    - :meth:`~.Alignment.consensus`: Merges alignment into a single sequence.
    - :meth:`~.Alignment.filter_codon_positions`: Filters alignment columns
      according to codon position.
    - :meth:`~.Alignment.filter_missing_data`: Filters columns according to
      missing data.
    - :meth:`~.Alignment.code_gaps`: Code alignment's indel patterns as a
      binary matrix.

For example, the :meth:`.Alignment.collapse` method transforms the original
alignment into a new one that contains only unique sequences. An important
factor to take into account with alignment modifying methods, is that it may
be important to preserve the original alignment data for future operations.
In TriFusion, the original alignment must be available at all times since
users may perform any number of process executions in a single session.
Therefore, all methods that can potentially modify the original alignment
need to be decorated with the :func:`.setup_database` function, and must be
defined with at least the `table_in` and `table_out` arguments. The
decorator and these arguments will work together to determine the database's
table from where the data will be fetched, and to where the modified
alignment will be written. For instance, the :meth:`.Alignment.collapse`
method is defined as::

    @setup_database
    def collapse(..., table_in=None, table_out="collapsed",
                 use_main_table=False):

If we want to perform a collapse of the original alignment and store the
modified alignment in a new table, we could call :meth:`.Alignment.collapse`
like::

    collapse(table_out="new_table")

The :func:`.setup_database` decorator interprets `table_in=None` as an
instruction to use the table with the original alignment, and stores the
modified alignment in `"new_table"`.

However, we may want to perform a collapse operation after a previous
modification from other method. In that case, we can specify a
`table_in`::

    collapse(table_in="other_table", table_out="collapse")

One issue with this approach is that we do not know *a priori* which
operations will be requested by the user nor the order. If one execution
performs, say, a `consensus` and a `collapse`, the new table should be
created in `consensus` and then used as input in `collapse`. However,
if only `collapse` is called, then the new table should only be created
there. To solve this issue, the `setup_database` decorator is smart about
its arguments. We can create a sequence of operations with the same
`table_in` and `table_out` arguments::

    new_table = "master_table"
    if "consensus" in operations:
        consensus(table_in=new_table, table_out=new_table)
    if "collapse" in operations:
        collapse(table_in=new_table, table_out=new_table)

In this simple pipeline, the user may perform either a `consensus`,
a `collapse`, or both. When the first method is called,
the :func:`.setup_database` decorator will check if the table provided in
`table_in` exists. In the first called method it will not exist, so instead
of returning an error, it falls back to the original alignment table and
then writes the modified alignment to "master_table". In the second method,
`table_in` already exists, so it fetches alignment data from the
"master_table". This will work whether these methods are called individually
or in combination.

When there is no need to keep the original alignment data (in single
execution of TriSeq, for instance), the special `use_main_table` argument
can be provided to tell the method to use the original table as the input
and output table. If this argument is True, it supersedes any information
provided by `table_in` or `table_out`::

    collapse(use_main_table=True)

Writers
~~~~~~~

Like parsers, writer methods are defined with `_write_<format>`:

    - :meth:`~.Alignment._write_fasta`: Writes to fasta format.
    - :meth:`~.Alignment._write_phylip`: Writes to phylip format.
    - :meth:`~.Alignment._write_nexus`: Writes to nexus format.
    - :meth:`~.Alignment._write_stockholm`: Writes to stockholm format.
    - :meth:`~.Alignment._write_gphocs`: Writes to gphocs format.
    - :meth:`~.Alignment._write_ima2`: Writes to IMa2 format.
    - :meth:`~.Alignment._write_mcmctree`: Writes to MCMCTree format.

They are always called from the :meth:`.Alignment.write_to_file` method,
not directly. When the :meth:`.Alignment.write_to_file` method is called,
a list with the requested output formats is also provided as an argument.
For each format specified in the argument, the corresponding writer method
is called. That method is responsible for fetching the data from the
database and write it to an output file in the appropriate format. The map
between the formats and the methods is stored in a dictionary::

    write_methods = {
        "fasta": self._write_fasta,
        "phylip": self._write_phylip,
        "nexus": self._write_nexus,
        "stockholm": self._write_stockholm,
        "gphocs": self._write_gphocs,
        "ima2": self._write_ima2,
        "mcmctree": self._write_mcmctree
    }

    # Call apropriate method for each output format
    for fmt in output_format:
        write_methods[fmt](output_file, **kwargs)

The `output_file` and a `kwargs` dictionary are provided as arguments to
each of these methods. The `kwargs` dictionary contains all keyword
arguments used when calling `write_to_file` and each writer method fetches
the ones relevant to the format. For instance, in the beginning of the
:meth:`~.Alignment._write_fasta` method, the relevant keyword arguments are
retrieved::

    # Get relevant keyword arguments
    ld_hat = kwargs.get("ld_hat", False)
    interleave = kwargs.get("interleave", False)
    table_suffix = kwargs.get("table_suffix", None)
    table_name = kwargs.get("table_name", None)
    ns_pipe = kwargs.get("ns_pipe", None)
    pbar = kwargs.get("pbar", None)

In addition to the `write_methods` dictionary, a dictionary mapping the
formats and their corresponding file extensions is also defined::

    format_ext = {"ima2": ".txt",
                  "mcmctree": "_mcmctree.phy",
                  "phylip": ".phy",
                  "nexus": ".nex",
                  "fasta": ".fas",
                  "stockholm": ".stockholm",
                  "gphocs": ".txt"}

To add a new output format writer, simply create the method using the
`_write_<format>` notation and include it in the `write_methods` dictionary,
along with the extension in the `format_ext` variable.

`AlignmentList` class
---------------------

The :class:`.AlignmentList` is the main interface between the user and the
alignment data. It may contain one or more :class:`.Alignment` objects,
which are considered the building blocks of the data set. These
:class:`.Alignment` objects are bound by the same sqlite database connection.

How to create an instance
~~~~~~~~~~~~~~~~~~~~~~~~~

An :class:`.AlignmentList` instance can be created with a single argument,
which is a list of paths to alignment files::

    aln_obj = AlignmentList(["file1.fas", "file2.fas"])

Even when no information is provided, an sqlite database connection is
established automatically (generating a "trifusion.sqlite3" file in the
current working directory by default). However, it is possible and advisable
to specify the path to the sqlite database::

    aln_obj = AlignmentList(["file1.fas", "file2.fas"], sql_path=".sql.db")

A connection to the database can also be provided when instantiating the
class, so it's perhaps more useful to see how the
:meth:`~.AlignmentList.__init__` works::

    def __init__(self, alignment_list, sql_db=None, db_cur=None, db_con=None,
                 pbar=None):

        if db_cur and db_con:
            self.con = db_con
            self.cur = db_cur
        elif sql_db:
            self.sql_path = sql_db
        else:
            self.sql_path = "trifusion.sqlite3"

As we can see, if a database connection is provided via the `db_cur` and
`db_con` arguments, the `sql_path` is ignored. If no sqlite information
is provided, the `sql_path` attribute defaults to "trifusion.sqlite3".

Adding alignment data
~~~~~~~~~~~~~~~~~~~~~

Alignment data can be loaded when initializing the class as above and/or
added later via the :meth:`~AlignmentList.add_alignments` and
:meth:`~AlignmentList.add_alignment_files` methods (in fact,
the :meth:`~AlignmentList.add_alignment_files` method is the one called at
:meth:`~.AlignmentList.__init__`)::

    aln_obj.add_alignment_files(["file3.fas", "file4.fas"])

The most important aspects when adding alignment data is how each alignment
is processed and how errors and exceptions are handled. Briefly, the flow is:

  1. Check for file duplications within the loaded file list. Duplications
     are stored in the `duplicate_alignments` attribute.
  2. Check for file duplications between the loaded files and any
     alignment already present. The duplications go to the same attribute
     as above.
  3. For each provided file, create an :class:`.Alignment` object. Errors
     that occur when creating an `Alignment` object are stored in its
     :attr:`~.Alignment.e`
     attribute. This attribute is then checked before adding the alignment to
     the :class:`.AlignmentList` object.
      1. Check for `InputError` exception. These are malformated files
         and are stored in the `bad_alignments` attribute.
      2. Check for `AlignmentUnequalLength` exception. These are sequence
         sets of
         unequal length (not alignments) and are stored in the
         `non_alignments` attribute.
      3. Check for `EmptyAlignment` exception. These are empty alignments
         and are
         stored in the `bad_alignments` attribute.
      4. Check the sequence code of the `Alignment` object. For now,
         the `AlignmentList` object only accepts alignments of the same
         sequence type (e.g. DNA, protein).
      5. If all of the previous checks pass, add the `Alignment` object
         to the `alignments` attribute.
  4. Update the overall taxa names of the full data set after the
     inclusion of the new alignments.

Wrapping `Alignment` methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The majority of the methods defined in the :class:`.Alignment` object can
also be accessed in the :class:`.AlignmentList` object. These are defined
roughly with the same arguments in both classes so that their behavior is
the same. These can be simple wrappers that call the respective
:class:`.Alignment` method for each alignment in the
:attr:`~.AlignmentList.alignments` attribute. For instance,
the :meth:`~.AlignmentList.change_taxon_name` method is simply::

    def change_taxon_name(self, old_name, new_name):

        for alignment_obj in list(self.alignments.values()):
            alignment_obj.change_taxon_name(old_name, new_name)

        self.taxa_names = [new_name if x == old_name else x
                           for x in self.taxa_names]

To avoid duplicating the argument list, the wrapping method can use `args`
and `kwargs` to transfer arguments. This ensures that if the argument list
is modified in the :class:`.Alignment` method, it doesn't need any
modification in the wrapper method. For instance,
the :meth:`~.Alignment.write_to_file` method of :class:`.Alignment` accepts
a large number of positional and keyword arguments, which would be an hassle
to define an maintain in the wrapper method of :class:`.AlignmentList`. So,
the :meth:`~:AlignmentList.write_to_file` method of :class:`.AlignmentList`
is simply defined as::

    def write_to_file(self, output_format, conversion_suffix="",
                      output_suffix="", *args, **kwargs):

The `output_format`, `conversion_suffix` and `output_suffix` are the only
positional arguments required when calling this wrapper. All the remaining
arguments are packed in the `args` and `kwargs` objects and used
normally by the wrapped method.

Exclusive `AlignmentList` methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some methods are exclusive of :class:`.AlignmentList` because they only make
sense to be applied to lists of alignments (e.g.
:meth:`~.AlignmentList.concatenate`). These have more freedom in how they
are defined and called.

Active and inactive datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Taxa and/or alignments can become 'inactive', that is, they are temporarily
removed from their respective attributes, `taxa_names` and `alignments`.
This means that these 'inactive' elements are ignored when performing most
operations. To change the 'active' status of alignments,
the :meth:`~AlignmentList.update_active_alignments` and
:meth:`~AlignmentList.update_active_alignment` methods are available. For
taxa, the :meth:`~AlignmentList.update_taxa_names` method can be used::

    # Set only two active alignments
    self.update_active_alignments(["file1.fas", "file2.fas"])

    # Set only two active taxa
    self.update_taxa_names(["taxonA", "taxonB"])

Note that all these modifications are reversible. 'Inactive' elements are
stored in the :attr:`~.AlignmentList.shelve_alignments` attribute for
alignments, and :attr:`~.AlignmentList.shelved_taxa` for taxa.

Updating `Alignment` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some tasks perform changes to core attributes of :class:`.AlignmentList`,
but they may also be necessary on each :class:`.Alignment` object. For
instance, the :meth:`~.AlignmentList.remove_taxa` method is used to remove a
list of taxa from the :class:`.AlignmentList` object. It is easy to change
only the relevant :class:`.AlignmentList` attributes, but this change also
requires those particular taxa to be removed in all :class:`.Alignment`
objects. For this reason, such methods should be defined in the same way in
both classes. Using the :meth:`~.AlignmentList.remove_taxa` example::

    def remove_taxa(self, taxa_list, mode="remove"):

        # <changes to AlignmentList>

        for alignment_obj in list(self.alignments.values()):
            alignment_obj.remove_taxa(taxa_list, mode=mode)

As you can see, the usage is the same for both methods.

Plot data methods
~~~~~~~~~~~~~~~~~

The :class:`.AlignmentList` object contains all methods that generate data
for plotting in TriFusion's Statistics screen and TriStats CLI program.
However, it's important to establish a separation between the generation of
plot data, and the generation of the plot itself. The
:class:`.AlignmentList` methods only generate the data and relevant
instructions necessary to to draw the plot. This information is then passed
on to the appropriate plot generation functions, which are defined in
:mod:`trifusion.base.plotter`. The reason for this separation of tasks is
that many different alignment analyses are represented by the same plot.

The complete process of how new plots can be added to TriFusion is
described here_. In this section, we provide only a few guidelines on what
to expect from these methods.

All plot data methods must be decorated with the :func:`.check_data` decorator
and take at least a `Namespace` argument. In most cases, no more arguments
are required::

    @check_data
    def gene_occupancy(self, ns=None):

The :func:`.check_data` decorator is responsible for performing checks
before and after executing the method. The `Namespace` argument, `ns`,
is used to allow communication between the main and worker threads of
TriFusion.

Additional keyword arguments may be defined, but in that case they must be
provided in TriFusion when the
:meth:`trifusion.app.TriFusionApp.stats_show_plot` method is called,
using the `additional_args` argument. This object will be passed to the
:func:`~trifusion.data.resources.background_tasks.get_stats_data` function
in :mod:`trifusion.data.resources.background_tasks` and used when calling
the plot data methods::

    if additional_args:
        plot_data = methods[stats_idx](ns=ns, **additional_args)
    else:
        plot_data = methods[stats_idx](ns)

The other requirement of plot data methods, is that they must always return
a single dictionary object. This dictionary must contain at least one
`key:value` with the "data" key and a numpy.array with the plot data as the
value. The other entries in the dictionary are optional and refer to
instructions for plot generation. For example, the
`missing_data_distribution` method returns::

    return {"data": data,
            "title": "Distribution of missing data",
            "legend": legend,
            "ax_names": ["Proportion", "Number of genes"],
            "table_header": ["Bin"] + legend}

The first entry is the mandatory "data" key with the numpy.array `data`.
The other instructions are "title", which sets the title of the plot, "legend",
which provides the labels for the plot legend, "ax_names", which provides
the name of the `x` and `y` axis, and the "table_header", which specified
the header for the table of that plot.

The allowed plot instructions depend on the plot function that will be used
and not all of them need to be specified.

.. _here: https://github.com/ODiogoSilva/TriFusion/wiki/
          Add-Statistics-plot-analysis

Logging progress
~~~~~~~~~~~~~~~~

The majority of :class:`.AlignmentList` methods support the setup and update of
progress indicators that can be used in TriFusion (GUI) and the CLI programs.
In the case of TriFusion, the progress indicator is provided via a
`multiprocessing.Namespace` object that transfers information between the
main thread (where the GUI changes take place) and the working thread.
In the case of the CLI programs, the indicator is provided via a `ProgressBar`
object. In either case, the setup, update and resetting of the progress
indicators is performed by the same methods.

At the beginning of an operation, the :meth:`~.AlignmentList_set_pipes`
method is called::

    self._set_pipes(ns, pbar, total=len(self.alignments))

The `Namespace` object is defined as `ns` and the `ProgressBar` is defined
as `pbar`. Usually, only one of them is provided, depending on whether it
was called from TriFusion or from a CLI program. We also set the total of
the progress indicator. In this case it's the number of alignments. In case
the operation is called from TriFusion using the `Namespace` object,
this method also checks the number of active alignments. If there is only
one active alignment, it sets a Namespace attribute that will silence the
progress logging of the :class:`.AlignmentList` object and receive the
information from the :class:`.Alignment` object instead::

    if ns:
        if len(self.alignments) == 1:
            ns.sa = True
        else:
            ns.sa = False

Then, inside the task, we can update the progress within a loop using
the :meth:`~.AlignmentList._update_pipes` method::

    for p, alignment_obj in enumerate(list(self.alignments.values())):
        self._update_pipes(ns, pbar, value=p + 1,
                           msg="Some message")

Here, we provide the `Namespace` and `ProgressBar` objects as before.
In addition we provide the `value` associated with each iteration.
Optionally, the `msg` argument can be specified, which is used exclusively
by TriFusion to show a custom message.

At the end of a task its good practice to reset the progress indicators
by using the :meth:`~.AlignmentList._reset_pipes` method::

    self._reset_pipes(ns)

Here, only the `Namespace` object is necessary, since the `ProgressBar`
indicator is automatically reset on :meth:`~.AlignmentList._set_pipes`.

Incorporating a kill switch
~~~~~~~~~~~~~~~~~~~~~~~~~~~

All time consuming methods of :class:`.AlignmentList` accept a `Namespace`
object when called from TriFusion, allowing communication between the main
thread and the worker thread. Since python `threads` cannot be forcibly
terminated like `processes`, most methods should listen to a kill switch
flag that is actually an attribute of the `Namespace` object. This kill
switch is already incorporated into the :meth:`~.AlignmentList._set_pipes`,
:meth:`~.AlignmentList._update_pipes` and
:meth:`~.AlignmentList._reset_pipes` methods. ::

    if ns:
        if ns.stop:
            raise KillByUser("")

What it does is to listen to a kill signal from TriFusion's window (it can
be clicking on the "Cancel" button, for example). When this kill signal is
received in the main thread, it sets the `Namespace.stop` attribute to True
in both main and worker threads. In the worker thread, when the `stop`
attribute evaluates to True, a custom `KillByUser` exception is raised,
immediately stopping the task. The exception is then handled in the
:func:`trifusion.data.resources.background_tasks.process_execution` function
and transmitted to the main thread.

Using SQLite
------------

A great deal of the high performance and memory efficiency of the `sequence`
module comes from the use of sqlite to store and manipulate alignment data
on disk rather than RAM. This means that parsing, modifying and writing
alignment data must be done with great care to ensure that only the
essential data is loaded into memory, while minimizing the number of (
expensive) database queries. This has some implications for the methods of
both :class:`.Alignment` and :class:`.AlignmentList` objects with respect to
how parsing, alignment modification and output writing is performed.

Implications for parsing
~~~~~~~~~~~~~~~~~~~~~~~~

When writing or modifying parsing methods it is important to take into
account that alignment files can be very large (> 1Gb) and loading the
entire data into memory should be avoided. Whenever possible, only the data
of a single taxon should be kept in memory before inserting it into the
database and then releasing that memory. For most formats, particularly
leave formats, it's fairly straightforward to do this. However, interleave
formats can fragment the data from each taxon across the entire file. Since
database insertions and updates are expensive, loading the data in each line
can greatly decrease the performance in these formats. Therefore,
it's preferable to read the alignment file once per taxon, join the entire
sequence of that taxon, and then insert it into the database. This ensures
that only sequence data from one taxon is kept in memory at any given time
and only a minimal number of database insertions are performed. It will also
result in the same file being parsed N times, where N is the number of taxa.
However, the decrease in speed is marginal, since most lines are actually
skipped, whereas the efficiency increase is substantial.

Implications for fetching data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Retrieving data from an sqlite database is not as simple as accessing python
native data structure. Therefore, a set of methods and generators have
been defined in the :class:`.Alignment` object to facilitate the interface with
the data in the database (See `Data fetching`_). When some kind of data
is required from the database, it is preferable to modify or create a
dedicated method, instead of interacting with the database directly. This
creates a layer of abstraction between the database and :class:`.Alignment`/
:class:`.AlignmentList` methods that greatly facilitates future updates.

Implications for alignment modification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When performing modification to the alignment data it is important to take
into account that the original alignment data may need to be preserved for
future operations. These methods must be defined and called using
appropriate decorators and arguments that establish the name of the
database table from where information will be retrieved, and the name of the
table where the information will be written (See `Alignment modifiers`_).

Implications for writing files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When writing alignment data into new output files, the same caution of the
alignment parsing is advised. It's quite easy to let the entire alignment
be loaded into RAM, particularly when writing in interleave format.

"""

import numpy as np
import pandas as pd
from collections import Counter, defaultdict, OrderedDict
import itertools
import re
import os
import io
import json
from os.path import join, basename, splitext, exists
from itertools import compress
from threading import Lock
import functools
import sqlite3

# TriFusion imports

try:
    import process
    from process.base import dna_chars, aminoacid_table, iupac, \
        iupac_rev, iupac_conv, Base
    from process.data import Partitions
    from process.data import PartitionException
    from process.error_handling import DuplicateTaxa, KillByUser, \
        InvalidSequenceType, InputError, EmptyAlignment, \
        MultipleSequenceTypes, SingleAlignment
except ImportError:
    import trifusion.process as process
    from trifusion.process.base import dna_chars, aminoacid_table, iupac, \
        iupac_rev, iupac_conv, Base
    from trifusion.process.data import Partitions
    from trifusion.process.data import PartitionException
    from trifusion.process.error_handling import DuplicateTaxa, KillByUser, \
        InvalidSequenceType, InputError, EmptyAlignment, \
        MultipleSequenceTypes, SingleAlignment

# import pickle
# TODO: Create a SequenceSet class for sets of sequences that do not conform
# to an alignment, i.e. unequal length.This would eliminate the problems of
# applying methods designed for alignments to sets of sequences with unequal
# length would allows these sets of sequences to have methods of their own.
# However, the __init__ of the current Alignment object should apply to both
# SequenceSet and Alignment classes. So, I'll have to re-structure the code
# somehow.
# TODO After creating the SequenceSet class, an additional class should be
# used to make the triage of files to either the Alignment or SequenceSet
# classes

# Lock mechanism to prevent concurrent access to sqlite database
lock = Lock()


class LookupDatabase(object):
    """Decorator handling hash lookup table with pre-calculated values.

    This decorator class is used to decorate class methods that calculate
    pairwise sequence similarity. To ensure proper functionality, the
    decorated method should be called first with a single "connect" argument
    and after finishing all calculations, with a single "disconnect" argument.
    These special method callings are necessary to setup and close the
    database storing the calculated values, respectively.

    Parameters
    ----------
    func : function
        Decorated function
    
    Attributes
    ----------
    func : function
        Decorated function
    con : None or sqlite.Connection
        Connection object of the sqlite database
    c : None or sqlite.Cursor
        Cursor object of the sqlite database
    
    Notes
    -----
    The `con` and `c` attributes are initialized when the decorated function
    is called with the single "connect" argument
    (e.g. decorated_func("connect")). When the decorated function is called
    with the single "disconnect" argument (e.g. decorated_func("disconnect")),
    the database changes are committed, the sqlite Connection is closed and
    the Cursor is reset.
    """

    def __init__(self, func):

        self.func = func
        self.con = None
        self.c = None

    def __call__(self, *args):
        """Wraps the call of `func`.
        
        When the decorated method is called, this code is wrapped around
        its execution. It accepts an arbitrary number of positional arguments
        but is is currently designed to decorate the `_get_similarity`
        method of the `AlignmentList` class. There are three expected calling
        modes:
        
            1. decorated_func("connect") : Initializes the Connection and
            Cursor objects of the database
            2. decorated_func(seq1, seq2, locus_lenght) : The typical main
            execution of the method, providing the sequence1 and sequence2
            strings along with the integer with their lenght
            3. decorated_func("disconnect") : Commits changes to database
            and closes Connection and Cursor objects.
        
        The path of this sqlite database is automatically obtained from
        the `AlignmentList` attribute `sql_path`. We use the path to the
        same directory, but with a different name for the database, "pw.db".
        
        Parameters
        ----------
        args : list
            List of positional arguments for decorated method
        
        Notes
        -----
        When the decorated method is called normally (i.e., not with a
        "connect" or "disconnect" argument"), it first creates an hash
        of the argument list (excluding `self`). This hash is used to query
        the database to check if a value for that combination of sequence
        strings and sequence length has already been calculated. If yes,
        the result is immediately  returned without executing the decorated
        method. If no, the method is executed to perform calculations. The
        result of this calculation is then stored in the database and the
        result is returned.

        """

        c_args = args[1:3]

        # Initialize database connection.
        if c_args[0] == "connect":

            # Get path to database based on the AlignmentList db
            tmp_dir = os.path.dirname(args[0].sql_path)

            self.con = sqlite3.connect(join(tmp_dir, "pw.db"))
            self.c = self.con.cursor()
            self.c.execute("PRAGMA synchronous = OFF")

            # Create table, if it does not yet exist
            if not self.c.execute("SELECT name FROM sqlite_master WHERE type="
                                  "'table' AND name='pw_table'").fetchall():
                self.c.execute("CREATE TABLE pw_table (hash integer"
                               ", seq_len integer"
                               ", val float, "
                               "effective_len float,"
                               "PRIMARY KEY (hash, seq_len))")

            return

        elif c_args[0] == "disconnect":

            self.con.commit()
            self.con.close()
            self.c = None

            return

        if self.c:

            h = hash(c_args)
            self.c.execute("SELECT * FROM pw_table WHERE hash = ?"
                           " AND seq_len = ?", (h, args[-1]))
            vals = self.c.fetchone()
            if vals:
                return vals[2], vals[3]
            else:
                value = self.func(*args)
                h = hash(c_args)
                self.c.execute("INSERT INTO pw_table VALUES (?, ?, ?, ?)",
                               (h, args[-1], value[0], value[1]))
                return value

        else:
            return self.func(*args)

    def __repr__(self):
        """ Return the function's docstring. """
        return self.func.__doc__

    def __get__(self, obj, objtype):
        """ Support instance methods """
        return functools.partial(self.__call__, obj)


def check_data(func):
    """Decorator handling the result from AlignmentList plotting methods.
    
    This should decorate all plotting methods from the `AlignmentList` class.
    It can be used to control warnings and handle exceptions when calling
    the plotting methods. Currently, it ensures that a list of methods are
    not executed when the `AlignmentList` instance contains a single
    alignment, and handles cases where the plotting methods return an
    empty array of data. Further control can be added here to prevent
    methods from being executed in certain conditions and handling
    certain outputs of the plotting methods.
    
    The only requirement is that, even when the plotting methods are not
    executed, this should always return a dictionary. In such cases,
    this dictionary should contain a single key:value with
    {"exception":<exception_type_string>}. The <exception_type_string> should
    be a simple string with an informative name that will be handled in
    the `stats_write_plot` of the `TriFusionApp` class.
    
    Parameters
    ----------
    func : function
        Decorated function

    Returns
    -------
    res : dict
        The value returned by `func` is always a dictionary. If the
        decorator prevents the execution of `func` for some reason, ensure
        that a dictionary is always return with a single key:value
        {"exception": <exception_type_string>}
    """

    # functools.wraps(func) allows the correct generation of docstrings
    # by sphinx of the decorated objects.
    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        # Plot methods that should not be allowed to continue with only one
        #  active alignment
        no_single_plot = ["outlier_missing_data", "outlier_missing_data_sp",
                          "outlier_segregating", "outlier_segregating_sp",
                          "outlier_sequence_size", "outlier_sequence_size_sp",
                          "average_seqsize_per_species", "average_seqsize",
                          "sequence_similarity", "sequence_segregation",
                          "length_polymorphism_correlation",
                          "taxa_distribution", "cumulative_missing_genes",
                          "gene_occupancy", "missing_data_distribution",
                          "missing_genes_average"]

        # Calling outlier method with a single alignments should immediately
        # raise an exception
        if len(args[0].alignments) == 1:
            if func.__name__ in no_single_plot:
                return {"exception": "single_alignment"}

        res = func(*args, **kwargs)

        if "data" in res:
            if np.asarray(res["data"]).any():
                return res
            else:
                return {"exception": "empty_data"}
        elif "exception" in res:
            return res

    return wrapper


def setup_database(func):
    """Decorator handling the active database tables.
    
    Decorates methods from the `Alignment` object that use and
    perform modifications to the original alignment data. All methods
    decorated with this must have the keyword arguments `table_in` and
    `table_out` (and `use_main_table`, optionally). The values associated with
    these arguments will determine which tables will be used and created
    before the execution of the decorated method (See Notes for the
    rationale). The strings provided as arguments for `table_in` and
    `table_out` will serve as a suffix to the `Alignment` instance attribute
    `db_idx`.
    
    If `use_main_table` is provided, and is set to True, both `table_in`
    and `table_out` will default to `Alignment.db_idx`. Values provided
    for `table_in` and `table_out` at calling time will be ignored.
    
    The following cases assume `use_main_table` is not provided or set to
    False.
    
    If both `table_in` and `table_out` are provided, `table_out` is modified
    so that `table_out = Alignment.db_idx + table_out`. If `table_out`
    is not provided, it defaults to `Alignment.db_idx`.
    
    If the final `table_out` does not exist in the database, create it.
    In this case, if `table_in` will default to `Alignment.db_idx`. If
    `table_out` already exists in the database, then `table_in=table_out`.
    
    Parameters
    ----------
    func : function
        Decorated function
    
    Notes
    -----
    In order to fully understand the mechanism behind the setup of the
    database tables, one must first know that when an `Alignment` instance
    is created, it generates a table in the database containing the
    original data from the alignment.
    In TriFusion (GUI), this original table MUST NOT be modified, since users
    may want to execute several methods on the same `Alignment` object.
    Therefore, when a particular method needs to modify the original alignment,
    a new temporary table is created to store the modified version
    until the end of the execution. If a second modification is requested,
    it will also be necessary to set the input table as the output table
    of the previous alignment modification.
    In the case of TriSeq (CLI), there is no such requirement,
    so it's much simpler to use the same table  for all modifications.
    
    This decorator greatly simplifies this process in the same way for
    all methods of the `Alignment` object that modify the original alignment
    data. To accomplish this, all decorated method must have the keyword
    arguments: `table_in` and `table_out` (and `use_main_table`, optionally).
    
    For methods called within the execution of TriFusion, the idea is simple.
    Since we don't know which methods will be used by the user, all chained
    methods that will create the same output file can be called with
    `table_in=new_table` and `table_out=new_table`.
    Note that bot arguments have the same value. Whatever is the first method
    being called, it will face the fact that "new_table" does not yet exist.
    In this case, the decorator will create a "new_table" in the database and
    reset `table_in` to the original `Alignment.db_idx`. This ensures that
    the first method will still be able to fetch the alignment data. In the
    following methods, `table_in` and `table_out` will be used based on the
    original values, ensuring that the alignment data is being fetched
    from the last modification made and that the modification chain is
    maintained. In this way, the execution of `Alignment` methods can
    be the same, regardless of the order or number of operations requested
    by the user.
    
    In the case of TriSeq, the methods should set `use_main_table=True`.
    This will always set `table_in=table_out=Alignment.db_idx`, which
    means that the main database will be used for all methods.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        # If use_main_table is set, the table_in will be overwritten by
        # table_out
        if "use_main_table" in kwargs:
            if kwargs["use_main_table"]:
                kwargs["table_out"] = args[0].table_name
                kwargs["table_in"] = args[0].table_name
                return func(*args, **kwargs)

        # Get sqlite database cursor and main table name
        sql_cur = args[0].cur

        # Define sqlite table name based on the main Alignment table name
        # and the provided db_idx argument
        if "table_out" in kwargs:
            table_out = args[0].table_name + kwargs["table_out"]
            kwargs["table_out"] = table_out
        else:
            table_out = args[0].table_name
            kwargs["table_out"] = args[0].table_name

        # Add table name to active table names
        if table_out not in args[0].tables:
            args[0].tables.append(table_out)

        # Check if the table already exists, if not create it.
        if not sql_cur.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND"
                " name='{}'".format(table_out)).fetchall():

            sql_cur.execute("CREATE TABLE [{}]("
                             "txId INT,"
                             "taxon TEXT,"
                             "seq TEXT)".format(table_out))
            # Set the table_in argument (table where the sequences will be
            # fecthed) to the main Alignment table name
            kwargs["table_in"] = args[0].table_name
        else:
            # Since the table_out already exists, sequences will be fetched
            # from this table instead of the main table name
            kwargs["table_in"] = table_out

        return func(*args, **kwargs)

    return wrapper


def setup_intable(func):
    """Decorator handling active database table when fetching data.
    
    This class is mean to decorate methods of the `Alignment` class that
    retrieves alignment data from the database. The requirement is that
    these methods have the positional arguments `table_suffix`,
    `db_idx` and `table`. These are all optional, and only `table_suffix`
    and `db_idx` should be used when calling these methods. The values
    of these two will be used to define the value of `table`, so any
    value provided to this argument when calling the decorated method will
    be ignored. In the end, the `table` variable will contain the final
    table name and only this variable will be used to retrieve the alignment
    data.
    
    `table_suffix` is a suffix that is appended to `Alignment.db_idx`.
    `db_idx` defines a complete table name.
    
    If neither `table_suffix` nor `db_idx` are provided, then `table`
    will default to `Alignment.db_idx`.
    
    If `db_idx` is provided, `table=db_idx`. If `table_suffix`
    is provided, `table=Alignment.db_idx + table_suffix`.
    
    If both are provided, `db_idx` takes precedence over `table_suffix`
    and `table=db_idx`.
    
    In any case, we test the existence of the final `table` value in the
    database. If it does not exist, `table` will revert to
    `Alignment.db_idx` to prevent errors.
    
    Parameters
    ----------
    func : function
        Decorated function
    
    Attributes
    ----------
    func : function
        Decorated function
    """

    # functools.wraps(func) allows the correct generation of docstrings
    # by sphinx of the decorated objects.
    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        # Get sqlite database cursor and main table name
        sql_cur = args[0].cur

        if "db_idx" not in kwargs and "table_suffix" not in kwargs:
            kwargs["table"] = args[0].table_name
            return func(*args, **kwargs)

        try:
            tb_name = kwargs["db_idx"]
        except KeyError:
            tb_name = None

        try:
            tb_suffix = kwargs["table_suffix"]
        except KeyError:
            tb_suffix = ""

        if tb_name:
            table = tb_name
        elif tb_suffix:
            table = args[0].table_name + tb_suffix
        else:
            table = args[0].table_name

        # Check if the specified table exists. If not, revert to
        # self.db_idx. When executing Process operations, the first valid
        # operation may be one that only has a "table_in" argument but several
        # other operations may be defined before with a different table suffix.
        # In such case, the table passed to "table_in" was not yet created and
        # populated in those prior methods. Therefore, here we account for
        # those cases by reverting the table to self.db_idx if table_in
        # was not yet created.
        if not sql_cur.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND"
                " name='{}'".format(table)).fetchall():
            table = args[0].table_name

        kwargs["table"] = table

        return func(*args, **kwargs)

    return wrapper


class AlignmentException(Exception):
    """ Generic Alignment object exception. """
    pass


class AlignmentUnequalLength(Exception):
    """ Raised when sequences in alignment have unequal length. """
    pass


class Alignment(Base):
    """Main interface for single alignment files.
    
    The `Alignment` class is the main interface for single alignment files,
    providing methods that parse, query, retrieve and modify alignment data.
    It is usually created within the context of `AlignmentList`,
    which handles the connection with the sqlite database. Nevertheless,
    `Alignment` instances can be created by providing the `input_alignment`
    argument, which can be one of the following:

    - path to alignment file. In this case, the file is parsed and
      information stored in a new sqlite table.

    - sqlite table name. In this case, the table name must be present
      in the database.

    In either case, the sqlite `Connection` and `Cursor` objects should
    be provided.
    
    When an `Alignment` object is instantiated, it first generates the
    `db_idx` based on the `input_alignment` string, filtering all
    characters that are not alpha numeric. Then, it queries the database
    to check is a table already exists with that name. If yes, it is assumed
    that the alignment data is stored in the provided table name. If there
    is no table with that name, it is assumed that `input_alignment` is a
    path to the alignment file and the regular parsing ensues. An empty
    table is created, the sequence type, format and missing data symbol are
    automatically detected and the alignment is parsed according to the
    detected format.
    
    Parameters
    ----------
    input_alignment : str
        Can be either the path to an alignment file or the name of a table
        in the sqlite database.
    sql_cursor : sqlite3.Cursor
        Cursor object of the sqlite database.
    input_format : str, optional
        File format of `input_alignment`. If `input_alignment` is a
        file path, the format will be automatically detect from the file.
        The value provided with this argument overrides the automatic
        detection's result.
    partitions : `trifusion.process.data.Partitions`, optional
        If provided, it will set the `_partitions` attribute. This should
        be used only when `input_alignment` is a database table name.
    locus_length : int, optional
        Sets the length of the current alignment, stored in the `locus_length`
        attribute. This option should only be used when `input_alignment`
        is a database table name. Otherwise, it is automatically set during
        alignment parsing.
    sequence_code : tuple, optional
        Sets the `sequence_code` attribute with the information on
        (<sequence_type>, <missing data symbol>). This option should only be
         used when `input_alignment` is a database table name. Otherwise,
         it is automatically set during alignment parsing.
    taxa_list : list, optional
        Sets the list attribute `taxa_list` with the names of the taxa
        present in the alignment. This option should only be
        used when `input_alignment` is a database table name. Otherwise,
        it is automatically set during alignment parsing.
    taxa_idx : dict, optional
        Sets the dictionary attribute `_taxa_idx` that maps the taxon names
        to their index in the sqlite database table. This option should only
        be used when `input_alignment` is a database table name. Otherwise,
        it is automatically set during alignment parsing.
    
    Attributes
    ----------
    cur : sqlite3.Cursor
        Cursor object of the sqlite database.
    db_idx : str
        Name of the sqlite database's table storing the sequence
        data.
    partitions : `trifusion.process.data.Partitions`
        Stores the _partitions and substitution model definition for the
        alignment.
    locus_length : int
        Length of the alignment in base pairs or residues.
    restriction_range : str
        Only used when gaps are coded in a binary matrix. Stores a string
        with the range of the restriction-type data that will encode
        gaps and will only be used when nexus is in the output format.
    e : None or Exception
        Stores any exceptions that occur during the parsing of the
        alignment file. It remains None unless something wrong happens.
    taxa_list : list
        List with the active taxon names.
    taxa_idx : dict
        Maps the taxon names to their corresponding index in the sqlite
        database. The index is not retrieved from the position of the taxon
        in `taxa_list` to prevent messing up when taxa are removed from the
        `Alignment` object.
    shelved_taxa : list
        List of ignored taxon names.
    path : str
        Full path to alignment file.
    sname : str
        Basename of the alignment file without the extension.
    name : str
        Basename of the alignment file with the extension.
    sequence_code : tuple
        Contains information on (<sequence type>, <missing data symbol>),
        e.g. ("Protein", "x").
    interleave_data : bool
        Attribute that is set to True when interleave data has been
        created for the alignment data.
    input_format : str
        Format of the input alignment file.
    
    Notes
    -----
    The `Alignment` class is not meant to be used directly (although it is
    quite possible to do so if you handle the database connections before).
    Instead, use the `AlignmentList` class, even if there is only one
    alignment. All `Alignment` methods are available from the `AlignmentList`
    and it is also possible to retrieve specific `Alignment` objects.

    The `Alignment` class was designed to be a lightweight, fast and
    powerful interface between alignment data and a set of manipulation
    and transformation methods.
    For performance and efficiency purposes, all alignment data is stored
    in a sqlite database that prevents the entire alignment from being
    loaded into memory. To facilitate the retrieval and iteration over the
    alignment data, several methods (`iter_columns`, 'iter_sequences`, etc)
    are available to handle the interface with the database and retrieve
    only the necessary information.

    See Also
    --------
    AlignmentList
    AlignmentList.add_alignment_files
    AlignmentList.retrieve_alignment
    """

    def __init__(self, input_alignment, input_format=None, partitions=None,
                 locus_length=None, sequence_code=None,
                 taxa_idx=None, sql_cursor=None,
                 db_idx=None, ignore_db_check=False, temp_dir=""):

        self.cur = sql_cursor

        if isinstance(partitions, Partitions):
            self._partitions = partitions
        else:
            self._partitions = Partitions()
            """
            Initializing a Partitions instance for the current alignment. By
            default, the Partitions instance will assume that the whole Alignment
            object is one single partition. However, if the current alignment is the
            result of a concatenation, or if a _partitions file is provided, the
            _partitions argument can be used. Partitions can be later changed using
            the set_partitions method. Substitution models objects are associated
            with each partition and they default to None
            """

        if not locus_length:
            self.locus_length = 0
        else:
            self.locus_length = locus_length
            """
            The length of the alignment object. Even if the current alignment
            object is partitioned, this will return the length of the entire
            alignment
            """

        self.restriction_range = None
        """
        This option is only relevant when gaps are coded. This will store a
        string with the range of the restriction-type data that will encode
        gaps and will only be used when nexus is in the output format
        """

        self.e = None
        """
        The `e` attribute will store any exceptions that occur during the
        parsing of the alignment object. It remains None unless something
        wrong happens.
        """

        self.shelved_taxa = []
        """
        Attribute that will store shelved taxa. When retrieving taxa from
        the database, this list will be checked and present taxa will
        be ignored
        """

        self._taxa_idx = None
        """
        Attribute that will store the index of the taxa in the sql database.
        This index is stored in a dictionary instead of being retrieved
        from the index position of the taxa_list list because taxa may
        be removed (which would mess up the link between the two indexes)
        """

        self.path = input_alignment
        """
        Attribute with full path to alignment file
        """
        self.sname = basename(splitext(input_alignment)[0])
        """
        Attribute with basename of alignment file without extension
        """
        self.name = basename(input_alignment)
        """
        Attribute with basename of alignment file
        """

        self.sequence_code = sequence_code

        if taxa_idx:
            self._taxa_idx = taxa_idx
        else:
            self._taxa_idx = OrderedDict()

        self.db_idx = db_idx
        """
        Creates a database table for the current alignment. This will be the
        link attribute between the sqlite database and the remaining methods
        that require information on the sequence data. This also removes
        any non alpha numeric characters the table name might have and ensures
        that it starts with a aphabetic character to avoid a database error
        """

        self.master_table = "alignment_data"
        """
        String attribute with the default name of the master table setup by
        the :class:`.AlignmentList` object.
        """

        self.interleave_data = False
        """
        Attribute that is set to True when interleave data has been
        created for the alignment data. Buiding the interleave matrix is a bit
        costly, so when it is requested by the user it is built once and
        stored in the database, and then further usages will use that table.
        """

        self.temp_dir = temp_dir if temp_dir else "."

        if not ignore_db_check:

            # Get alignment format and code. Sequence code is a tuple of
            # (DNA, N) or (Protein, X)
            finder_content = self.autofinder(input_alignment)
            # Handles the case where the input format is invalid and
            # finder_content is an Exception
            if isinstance(finder_content, Exception) is False:
                self.input_format, self.sequence_code = self.autofinder(
                    input_alignment)

                # In case the input format is specified, overwrite the
                # attribute
                if input_format:
                    self.input_format = input_format
                    """
                    Format of the alignment file.
                    """

                # parsing the alignment
                self.read_alignment()
            else:
                # Setting the sequence code attribute for seq type checking
                # in AlignmentList
                self.sequence_code = None
                self.e = finder_content

        # In case there is a table for the provided input_alignment
        else:
            self.input_format = input_format
            if taxa_idx:
                self.store_taxa_idx(self.temp_dir, taxa_idx)
            if partitions:
                self.store_partitions(self.temp_dir, partitions)

    def __iter__(self):
        """Iterator behavior for `Alignment`.
        
        Iterates over taxa and sequence data from the alignment. Taxon names
        in the `shelved_taxa` attribute are ignored. This will always
        retrieve the alignment data from the master database table,
        `db_idx`.
        
        Yields
        ------
        tx : str
            taxon name.
        seq : str
            sequence string.
        """

        for tx, seq in self.cur.execute(
                "SELECT taxon,seq from alignment_data WHERE aln_idx=?",
                (self.db_idx,)):
            if tx not in self.shelved_taxa:
                yield tx, seq

    def _create_table(self, table_name, index=None, cur=None):
        """Creates a new table in the database.
        
        Convenience method that creates a new table in the sqlite database.
        It accepts a custom Cursor object, which overrides the `cur`
        attribute.
        
        Parameters
        ----------
        table_name : str
            Name of the table to be created.
        index : list, optional
            Provide a list with [<name of index>, <columns to be indexed>].
            (e.g., ["concindex", "txId"]).
        cur : sqlite3.Cursor, optional
            Custom Cursor object used to query the database.
        """

        if not cur:
            cur = self.cur

        cur.execute("CREATE TABLE [{}]("
                    "txId INT,"
                    "taxon TEXT,"
                    "seq TEXT,"
                    "aln_idx INT)".format(table_name))

        if index:
            cur.execute("CREATE INDEX {} ON [{}]({})".format(
                index[0], table_name, index[1]))

    def _set_format(self, input_format):
        """Manually set the input format of the `Alignment` object.
        
        Manually sets the input format associated with the `Alignment` object

        Parameters
        ----------
        input_format : str
            The input format.
        """

        self.input_format = input_format

    def get_taxa_idx(self):

        return json.load(open(join(self.temp_dir, str(self.db_idx))),
                         object_pairs_hook=OrderedDict)

    def store_taxa_idx(self, temp_path, tx_idx=None):

        if not self.temp_dir:
            self.temp_dir = temp_path

        if not self.temp_dir:
            self.temp_dir = "."

        tx_idx = tx_idx if tx_idx else self._taxa_idx

        with open(join(self.temp_dir, str(self.db_idx)), "w") as fh:
            json.dump(tx_idx, fh)

        self._taxa_idx = None

    def get_partitions(self):

        with open(join(self.temp_dir, str(self.db_idx) + "_p")) as fh:
            part = Partitions()
            part.__dict__ = json.load(fh, object_pairs_hook=OrderedDict)

        return part

    def store_partitions(self, temp_path, partitions_obj=None):

        if not self.temp_dir:
            self.temp_dir = temp_path if temp_path else "."

        part = partitions_obj if partitions_obj else self._partitions

        with open(join(self.temp_dir, str(self.db_idx) + "_p"), "w") as fh:
            json.dump(part.__dict__, fh)

        self._partitions = None

    @staticmethod
    def _set_pipes(ns=None, pbar=None, total=None, msg=None):
        """Setup of progress indicators for both GUI and CLI task executions.
         
         This handles the setup of the objects responsible for updating
         the progress of task's execution of both TriFusion (GUI) and
         TriSeq (CLI). At the beginning of any given task, these objects
         can be initialized by providing either the Namespace object (`ns`)
         in the case of TriFusion, or the ProgressBar object (`pbar`), in
         the case of TriSeq. Along with one of these objects, the expected
         `total` of the progress should also be provided. The `ns` and
         `pbar` objects are updated at each iteration of a given task,
         and the `total` is used to get a measure of the progress.

         Optionally, a message can be also provided for the Namespace object
         that will be used by TriFusion.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.
        total : int
            Expected total of the task's progress.
        msg : str, optional
            A secondary message that appears in TriFusion's progress dialogs.

        Notes
        -----
        The progress of any given task can be provided by either an
        `Alignment` or `AlignmentList` instance. Generally, the tasks follow
        the progress of the `AlignmentList` instance, unless that instance
        contains only one `Alignment` object. In that case, the progress
        information is piped from the `Alignment` instance. For that reason,
        and for the Namespace (`ns`) object only, this method first checks
        if the `ns.sa` attribute exists. If it does, it means that the parent
        `AlignmentList` instance contains only one `Alignment` object and,
        therefore, the progress information shoud come from here.

        Examples
        --------
        Start a progress counter  for a task that will make 100 iterations::

          self._set_pipes(ns=ns_obj, pbar=pbar_obj, total=100,
          msg="Some message")

        See Also
        --------
        _reset_pipes
        _update_pipes
        """

        # Check if the `sa` attribute exists in the Namespace object.
        # If yes, then the logging will be performed here. If not,
        # The progress attributes will not be set for `ns`.
        try:
            sa = ns.sa
        except AttributeError:
            sa = True

        if ns:
            # Listen for kill switch
            if ns.stop:
                raise KillByUser("")
            # Set attributes only if parent AlignmentList contains a single
            # Alignment object
            if sa:
                ns.total = total
                ns.counter = 0
                ns.msg = msg

        if pbar:
            pbar.max_value = total
            # Reset ProgressBar object counter
            pbar.update(0)

    @staticmethod
    def _update_pipes(ns=None, pbar=None, value=None, msg=None):
        """Update progress indicators for both GUI and CLI task executions.

        This method provides a single interface for updating the progress
        objects `ns` or `pbar`, which should have been initialized at the
        beginning of the task with the `_set_pipes` method.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.
        value : int
            Value of the current progress index
        msg : str, optional
            A secondary message that appears in TriFusion's progress dialogs.

        Notes
        -----
        The progress of any given task can be provided by either an
        `Alignment` or `AlignmentList` instance. Generally, the tasks follow
        the progress of the `AlignmentList` instance, unless that instance
        contains only one `Alignment` object. In that case, the progress
        information is piped from the `Alignment` instance. For that reason,
        and for the Namespace (`ns`) object only, this method first checks
        if the `ns.sa` attribute exists. If it does, it means that the parent
        `AlignmentList` instance contains only one `Alignment` object and,
        therefore, the progress information shoud come from here.

        Examples
        --------
        Update the counter in an iteration of 100::

            for i in range(100):
                self._update_pipes(ns=ns_obj, pbar=pbar_obj, value=i,
                                   msg="Some string")

        See Also
        --------
        _set_pipes
        _reset_pipes
        """

        try:
            sa = ns.sa
        except AttributeError:
            sa = True

        if ns:
            if ns.stop:
                raise KillByUser("")
            if sa:
                ns.counter = value
                ns.msg = msg
        if pbar:
            pbar.update(value)

    @staticmethod
    def _reset_pipes(ns):
        """Reset progress indicators for both GUI and CLI task executions.

        This should be done at the end of any task that initialized the
        progress objects, but it only affects the Namespace object. It
        resets all Namespace attributes to None.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        See Also
        --------
        _set_pipes
        _update_pipes
        """

        try:
            sa = ns.sa
        except AttributeError:
            sa = True

        if ns:
            if ns.stop:
                raise KillByUser("")
            if sa:
                ns.total = ns.counter = ns.msg = ns.sa = None

    def iter_sequences(self, table_name=None):
        """ Generator over sequence strings in the alignment.

        Generator for sequence data of the `Alignment` object.
        Sequence data is retrieved from a database table specified either by
        the `db_idx` or `table_suffix` arguments. `db_idx` will always
        take precedence over `table_suffix` if both are provided. If none
        are provided, the default `Alignment.db_idx` is used. If the
        table name provided by either `db_idx` or `table_suffix` is
        invalid, the default `Alignment.db_idx` is also used.

        Parameters
        ----------
        table_suffix : string
            Suffix of the table from where the sequence data is fetched.
            The suffix is appended to `Alignment.db_idx`.
        table_name : string
            Name of the table from where the sequence data is fetched.
        table : string
            This argument will be automatically setup by the `SetupInTable`
            decorator. Do not use directly.

        Yields
        ------
        seq : string
            Sequence string for a given taxon

        See Also
        --------
        SetupInTable

        Notes
        -----
        All generators ignore data associated with taxa present in the
        `shelved_taxa` attribute.
        """

        table_name = table_name if table_name else self.master_table

        try:
            # Locking mechanism necessary to avoid concurrency issues when
            # accessing the database. This ensures that only one Cursor
            # object is querying the database at any given time
            lock.acquire(True)
            for tx, seq in self.cur.execute(
                    "SELECT taxon,seq "
                    "FROM [{}] "
                    "WHERE aln_idx=?".format(table_name), (self.db_idx, )):
                if tx not in self.shelved_taxa:
                    yield seq
        finally:
            lock.release()

    def iter_alignment(self, table_name):
        """Generator for (taxon, sequence) tuples.

        Generator that yields (taxon, sequence) tuples from the database.
        Sequence data is retrieved from a database table specified either by
        the `db_idx` or `table_suffix` arguments. `db_idx` will always
        take precedence over `table_suffix` if both are provided. If none
        are provided, the default `Alignment.db_idx` is used. If the
        table name provided by either `db_idx` or `table_suffix` is
        invalid, the default `Alignment.db_idx` is also used.

        Parameters
        ----------
        table_suffix : string
            Suffix of the table from where the sequence data is fetched.
            The suffix is append to `Alignment.db_idx`.
        table_name : string
            Name of the table from where the sequence data is fetched.
        table : string
            This argument will be automatically setup by the `SetupInTable`
            decorator. Do not use directly.

        Yields
        ------
        tx : string
            Taxon name.
        seq : string
            Sequence string.

        See Also
        --------
        SetupInTable

        Notes
        -----
        All generators ignore data associated with taxa present in the
        `shelved_taxa` attribute.
        """

        table_name = table_name if table_name else self.master_table

        try:
            # Locking mechanism necessary to avoid concurrency issues when
            # accessing the database. This ensures that only one Cursor
            # object is querying the database at any given time
            lock.acquire(True)
            for tx, seq in self.cur.execute(
                    "SELECT taxon, seq "
                    "FROM [{}] "
                    "WHERE aln_idx=?".format(table_name), (self.db_idx,)):
                if tx not in self.shelved_taxa:
                    yield tx, seq
        finally:
            lock.release()

    def get_sequence(self, taxon, table_name=None):
        """Returns the sequence string for a given taxon.

        Returns the sequence string of the corresponding `taxon`. If the
        taxon does not exist in the table, raises a KeyError.
        The sequence is retrieved from a database table specified either by
        the `db_idx` or `table_suffix` arguments. `db_idx` will always
        take precedence over `table_suffix` if both are provided. If none
        are provided, the default `Alignment.db_idx` is used. If the
        table name provided by either `db_idx` or `table_suffix` is
        invalid, the default `Alignment.db_idx` is also used.

        Parameters
        ----------
        taxon : str
            Name of the taxon.
        table_suffix : string
            Suffix of the table from where the sequence data is fetched.
            The suffix is append to `Alignment.db_idx`.
        table_name : string
            Name of the table from where the sequence data is fetched.
        table : string
            This argument will be automatically setup by the `SetupInTable`
            decorator. Do not use directly.

        Returns
        -------
        seq : str
            Sequence string.

        Raises
        ------
        KeyError
            If the taxon does not exist in the specified table.

        Notes
        -----
        Ignores data associated with taxa present in the `shelved_taxa`
        attribute.
        """

        table_name = table_name if table_name else self.master_table

        tx_idx = self.get_taxa_idx()

        try:
            # Locking mechanism necessary to avoid concurrency issues when
            # accessing the database. This ensures that only one Cursor
            # object is querying the database at any given time
            lock.acquire(True)
            if taxon in tx_idx and \
                    taxon not in self.shelved_taxa:
                seq = self.cur.execute(
                    "SELECT seq "
                    "FROM [{}] "
                    "WHERE txId=? "
                    "AND aln_idx=?".format(table_name),
                    (tx_idx[taxon], self.db_idx)).fetchone()[0]
                return seq

            else:
                raise KeyError
        finally:
            lock.release()

    def shelve_taxa(self, lst):
        """Shelves taxa from `Alignment` methods.

        The taxa provided in the `lst` argument will be ignored in
        subsequent operations. In practice, taxa from the `taxa_list`
        attribute that are present in `lst` will be moved to the
        `shelved_taxa` attribute and ignored in all methods of the class.
        To revert all taxa to the `taxa_list` attribute, simply call this
        method with an empty list.

        Parameters
        ----------
        lst : list
            List with taxa names that should be ignored.
        """

        self.shelved_taxa = [x for x in lst if x in self.get_taxa_idx()]

    def _insert_data(self, txId, taxon, seq):
        """

        Parameters
        ----------
        txId
        taxon
        seq

        Returns
        -------

        """
        try:

            lock.acquire(True)

            self.cur.execute("INSERT INTO alignment_data VALUES (?, ?, ?, ?)",
                             (txId, unicode(taxon), seq, self.db_idx))

        finally:
            lock.release()


    def _read_interleave_phylip(self, ntaxa):
        """ Alignment parser for interleave phylip format.

        Parses an interleave phylip alignment file and stores taxa and
        sequence data in the database. This method is only called from the
        `_read_phylip` method, when the regular phylip parser detects that the
        file is in interleave format.

        Parameters
        ----------
        ntaxa : int
            Number of taxa contained in the alignment file.

        Returns
        -------
        size_list : list
            List of sequence size (int) for each taxon in the alignment file.
            Used to check size consistency at the end of the alignment parsing
            in `read_alignment`.

        See Also
        --------
        _read_phylip

        Notes
        -----
        The phylip interleave format splits the alignment into blocks of
        a certain lenght (usually 90 characters) separated by blank lines.
        This means that to gather the complete sequence of a given taxon,
        the parser has to read the entire file. To prevent the entire
        alignment to be loaded into memory, we actually iterate over a range
        determined by the number of taxa in the alignment. In each iteration,
        we open a file handle and we retrieve only the sequence of
        a particular taxon. This ensures that only sequence data for
        a single taxon is store in memory at any given time. This also means
        that the alignment file has to be read N times, where N = number of
        taxa. However, since the vast majority of lines are actually skipped
        in each iteration, the decrease in speed is marginal, while the
        gains in memory efficient are much larger.
        """

        size_list = []

        for i in xrange(ntaxa):

            # Variable that will store the sequence for the current taxon.
            sequence = []
            # Index identifier of the current taxon
            idx = 0
            # When True, it means that the alignment line has the taxon
            # and sequence. When False, it means that the line contains
            # only sequence.
            taxa_gather = True

            fh = open(self.path)

            # Skip header an any potential blank lines
            header = ""
            while not header:
                header = next(fh)

            for line in fh:

                # At blank lines, reset the idx and set taxa_gather to False.
                # All lines from this point on will only contain sequence.
                if not line.strip():
                    idx = 0
                    taxa_gather = False
                else:
                    # If the idx of the alignment block matches the taxon
                    # from the current iteration, add to sequence variable
                    if idx == i:
                        # Remove the taxon from the gathering
                        if taxa_gather:
                            sequence.append(
                                "".join(line.strip().lower().split()[1:]))
                        else:
                            sequence.append(
                                "".join(line.strip().lower().split()))
                    idx += 1

            seq = "".join(sequence)
            taxa = self._taxa_idx.keys()[i]

            size_list.append(len(seq))

            # Insert data into database
            self._insert_data(i, taxa, seq)

            fh.close()

        return size_list

    def _read_interleave_nexus(self, ntaxa):
        """ Alignment parser for interleave nexus format.

        Parses an interleave nexus alignment file and stores taxa and
        sequence data in the database.This method is only called from the
        `_read_nexus` method, when the regular nexus parser detects that the
        file is in interleave format.

        Parameters
        ----------
        ntaxa : int
            Number of taxa contained in the alignment file.

        Returns
        -------
        size_list : list
            List of sequence size (int) for each taxon in the alignment file.
            Used to check size consistency at the end of the alignment parsing
            in `read_alignment`.

        See Also
        --------
        _read_nexus

        Notes
        -----
        The nexus interleave format splits the alignment into blocks of
        a certain lenght (usually 90 characters) separated by blank lines.
        This means that to gather the complete sequence of a given taxon,
        the parser has to read the entire file. To prevent the entire
        alignment to be loaded into memory, we actually iterate over a range
        determined by the number of taxa in the alignment. In each iteration,
        we open a file handle and we retrieve only the sequence of a
        particular taxon. This ensures that only sequence data for
        a single taxon is store in memory at any given time. This also means
        that the alignment file has to be read N times, where N = number of
        taxa. However, since the vast majority of lines are actually skipped
        in each iteration, the decrease in speed is marginal, while the
        gains in memory efficient are much larger.
        """

        size_list = []

        for i in xrange(ntaxa):

            # Variable that will store the sequence for the current taxon.
            sequence = []
            # counter used to skip the Nexus header and footer
            counter = 0
            # Index identifier of the current taxon
            idx = 0
            taxa = None

            fh = open(self.path)

            for line in fh:

                if line.strip().lower() == "matrix" and counter == 0:
                    counter = 1

                elif line.strip() == ";" and counter == 1:
                    counter = 2

                elif counter == 1:

                    if not line.strip():
                        idx = 0
                    else:
                        if idx == i:
                            if not taxa:
                                taxa = line.strip().split()[0].replace(" ",
                                                                       "")
                                taxa = self.rm_illegal(taxa)

                            sequence.append(
                                "".join(line.strip().lower().split()[1:]))

                        idx += 1

            self._taxa_idx[taxa] = i

            seq = "".join(sequence)
            if not self.locus_length:
                self.locus_length = len(seq)

            size_list.append(len(seq))

            self._insert_data(i, taxa, seq)

            fh.close()

    def _eval_missing_symbol(self, sequence):
        """Evaluates missing data symbol from sequence string.

        This method is performed when sequence data is being parsed from the
        alignment and only executes when the missing data symbol stored in
        `Alignment.sequence_code` is not defined. It attempts to count
        the regular characters used to denote missing data. First, tries to
        find "?" characters, which are set as the symbol if they exist.
        Then, it finds "n" characters and if sequence type is set to DNA,
        this character is set as the symbol.
        Finally, it finds "x" characters and if sequence type is set to
        Protein, this character is set as the symbol.

        Parameters
        ----------
        sequence : str
            Sequence string

        """

        if not self.sequence_code[1]:

            if sequence.count("?"):
                self.sequence_code[1] = "?"

            elif sequence.count("n") and self.sequence_code[0] == "DNA":
                self.sequence_code[1] = "n"

            elif sequence.count("x") and self.sequence_code[0] == "Protein":
                self.sequence_code[1] = "x"

    def _read_phylip(self):
        """Alignment parser for phylip format.

        Parses a phylip alignment file and stored taxa and sequence data in
        the database.

        See Also
        --------
        read_alignment
        """

        fh = open(self.path)

        # Variable storing the lenght of each sequence
        size_list = []

        # Get the number of taxa and sequence length from the file header
        header = fh.readline().split()
        self.locus_length = int(header[1])
        self._partitions.set_length(self.locus_length)
        taxa_num = int(header[0])

        # These three following attributes allow the support for
        # interleave phylip files
        # Flags whether the current line should be parsed for taxon name
        taxa_gather = True
        # Counter that makes the correspondence between the current line
        # and the appropriate taxon
        c = 0

        for line in fh:

            # Ignore empty lines
            if line.strip() == "":
                continue

            # The counter is reset when surpassing the number of
            # expected taxa.
            if c + 1 > taxa_num:
                c = 0
                taxa_gather = False

            try:

                # Here, assume that all lines with taxon names have
                # already been processed, since he taxa_pos variable
                # already has the same number as the expected taxa in the
                #  phylip header
                if not taxa_gather:

                    # Oh boy, this seems like an interleave phylip file.
                    # Redirect parsing to appropriate method
                    # self.cur.execute("DELETE FROM [{}]".format(
                    #     self.db_idx))
                    size_list = self._read_interleave_phylip(taxa_num)
                    break

                # To support interleave phylip, while the taxa_pos
                # variable has not reached the expected number of taxa
                # provided in header[0], treat the line as the first
                # lines of the phylip where the first field is the taxon
                # name
                elif taxa_gather:

                    taxa = line.split()[0].replace(" ", "")
                    taxa = self.rm_illegal(taxa)

                    try:
                        # Joint multiple batches of sequence
                        # Remove any possible whitespace by splitting
                        # according to whitespace. This also ensures that
                        # sequence is always a list when added to
                        # sequence_data
                        sequence = line.strip().lower().split()[1:]
                    except IndexError:
                        sequence = [""]

                    self._taxa_idx[taxa] = c

                    seq = "".join(sequence)

                    # Evaluate missing data symbol if undefined
                    self._eval_missing_symbol(seq)

                    self._insert_data(c, taxa, seq)

                    size_list.append(len(seq))

                    # Add counter for interleave processing
                    c += 1

            except IndexError:
                pass

        # Updating _partitions object
        self._partitions.add_partition(self.name, self.locus_length,
                                       file_name=self.path)
        
        fh.close()

        # Checks the size consistency of the alignment
        if len(set(size_list)) > 1:
            self.e = AlignmentUnequalLength()

    def _read_fasta(self):
        """Alignment parser for fasta format.

        Parses a fasta alignment file and stores taxa and sequence data in
        the database.

        See Also
        --------
        read_alignment
        """

        fh = open(self.path)

        # Variable storing the lenght of each sequence
        size_list = []

        sequence = []
        taxa = None
        idx = 0
        for line in fh:
            if line.strip().startswith(">"):

                if sequence:
                    seq = "".join(sequence)

                    # Evaluate missing data symbol if undefined
                    self._eval_missing_symbol(seq)

                    self._insert_data(idx, taxa, seq)

                    self._taxa_idx[taxa] = idx

                    if not self.locus_length:
                        self.locus_length = len(seq)

                    size_list.append(len(seq))

                    # sequence_data.append((taxa, "".join(sequence)))
                    sequence = []
                    idx += 1

                taxa = line[1:].strip()
                taxa = self.rm_illegal(taxa)

            elif line.strip() != "" and taxa:
                sequence.append(line.strip().lower().
                                replace(" ", "").replace("*", ""))

        if sequence:
            seq = "".join(sequence)

            # Evaluate missing data symbol if undefined
            self._eval_missing_symbol(seq)

            self._insert_data(idx, taxa, seq)

            self._taxa_idx[taxa] = idx

            if not self.locus_length:
                self.locus_length = len(seq)

        size_list.append(len(seq))

        self._partitions.set_length(self.locus_length)

        # Updating _partitions object
        self._partitions.add_partition(self.name, self.locus_length,
                                       file_name=self.path)

        fh.close()

        # Checks the size consistency of the alignment
        if len(set(size_list)) > 1:
            self.e = AlignmentUnequalLength()

    def _read_loci(self):
        """Alignment parser for pyRAD and ipyrad loci format.

        See Also
        --------
        read_alignment
        """

        # Create temporary table
        temp_table = ".locidata"
        self._create_table(temp_table, index=("lociindex", "txId"))

        fh = open(self.path)

        # Variable storing the length of each sequence
        size_list = []

        taxa_list = self.get_loci_taxa(self.path)
        self._taxa_idx = dict((x, p) for p, x in enumerate(taxa_list))

        # Create empty dict
        sequence_data = []

        # Add a counter to name each locus
        locus_c = 1
        # This variable is used in comparison with taxa_list to check which
        # taxa are missing from the current partition. Missing taxa
        # are filled with missing data.
        present_taxa = []

        # Set default missing data symbol as "n"
        if not self.sequence_code[1]:
            self.sequence_code[1] = "n"

        for line in fh:
            # Parse a line with sequence data
            if not line.strip().startswith("//") and line.strip() != "":
                fields = line.strip().split()
                taxon = fields[0].lstrip(">")
                present_taxa.append(taxon)
                cur_seq = fields[1].lower()
                sequence_data.append(
                    (self._taxa_idx[taxon], taxon, cur_seq, self.db_idx))
                size_list.append(len(cur_seq))

            # End of a partition in the loci file
            elif line.strip().startswith("//"):

                # Checks the size consistency of the previous _partitions
                if len(set(size_list)) > 1:
                    self.e = AlignmentUnequalLength()

                # Reset list of sequence size for next partition
                size_list = []

                # Get length of previous partition based on the last
                # sequence
                locus_len = len(fields[1])
                self.locus_length += locus_len

                # Add partition
                self._partitions.add_partition("locus_{}".format(locus_c),
                                               locus_len,
                                               file_name=self.path)
                
                # Adding missing data
                for tx in self._taxa_idx:
                    if tx not in present_taxa:
                        sequence_data.append(
                            (self._taxa_idx[tx],
                             tx,
                             self.sequence_code[1] * locus_len,
                             self.db_idx))

                locus_c += 1

                present_taxa = []

                # Insert previous partition in the database
                self.cur.executemany(
                    "INSERT INTO [{}] VALUES (?, ?, ?, ?)".format(
                        temp_table), sequence_data)

                # Reset sequence data for next partition
                sequence_data = []

        self._partitions.set_length(self.locus_length)

        fh.close()

        # Add temp table to master table
        self.cur.execute(
            "INSERT INTO alignment_data (txId, taxon, seq, aln_idx) "
            "SELECT txId, taxon, GROUP_CONCAT(seq, ''), {} "
            "FROM [{}] "
            "GROUP BY txId".format(self.db_idx, temp_table))

        self.cur.execute("DROP TABLE [{}]".format(temp_table))

        return size_list

    def _read_nexus(self):
        """Alignment parser for nexus format.

        Parses a nexus alignment file and stores taxa and sequence data in
        the database.

        See Also
        --------
        read_alignment
        """

        fh = open(self.path)

        # Variable storing the lenght of each sequence
        size_list = []

        counter = 0
        idx = 0
        ntaxa = None
        interleave = None
        for line in fh:

            # Fetch the number of taxa from nexus header. This will be
            # necessary for efficient interleave parsing
            if not ntaxa:
                try:
                    ntaxa = int(re.search(r".*ntax=(.+?) ",
                                          line).group(1))
                except ValueError:
                    self.e = InputError("Could not recognize number"
                                        " of taxa from nexus header")
                    return
                except AttributeError:
                    pass

            # Fetch the number of sites from nexus header. This will be
            # necessary for efficient interleave parsing
            if not self.locus_length:
                try:
                    self.locus_length = int(
                        re.search(r".*nchar=(.+?)[ ,;]",
                                  line).group(1))
                except ValueError:
                    self.e = InputError("Could not recognize number"
                                        " of sites from nexus header")
                    return
                except AttributeError:
                    pass

            # Fetch the interleave value from the  nexus header. This
            # will be necessary for efficient interleave parsing
            if interleave == None:
                try:
                    interleave = re.search(r"interleave=(.+?) ",
                                           line).group(1).lower()
                    interleave = True if interleave.strip() == "yes" \
                        else False
                except AttributeError:
                    pass

            # Set the counter to beggin data parsing
            if line.strip().lower() == "matrix" and counter == 0:
                counter = 1

            # Stop sequence parser here
            elif line.strip() == ";" and counter == 1:
                counter = 2

                # self.locus_length = len(sequence_data[0][2])
                self._partitions.set_length(self.locus_length)

            # Start parsing here
            elif counter == 1:

                # If neither ntaxa or interleave are defined, return
                # with an error in the 'e' attribute
                if not ntaxa or interleave == None:
                    self.e = InputError("Could not recognize ntax or"
                                        " interleave parameters from"
                                        " nexus header")
                    return

                if interleave == False:

                    # In interleave, this marks the change in matrix blocks
                    if not line.strip():
                        continue

                    # Get taxon name
                    taxa = line.strip().split()[0].replace(" ", "")
                    taxa = self.rm_illegal(taxa)

                    # Update taxa_list and _taxa_idx attributes
                    self._taxa_idx[taxa] = idx

                    # Get sequence string
                    seq = "".join(line.strip().split()[1].split())

                    # Evaluate missing data symbol if undefined
                    self._eval_missing_symbol(seq)

                    size_list.append(len(seq))

                    # Add sequence to sqlite database
                    self._insert_data(idx, taxa, seq)

                    idx += 1

                else:
                    self._read_interleave_nexus(ntaxa)
                    counter = 2

            # If _partitions are specified using the charset command, this
            # section will parse the _partitions
            elif line.strip().startswith("charset"):
                self._partitions.read_from_nexus_string(line,
                                                        file_name=self.path)

            # If substitution models are specified using the lset or prset
            # commands, this will parse the model parameters
            if ("lset" in line.lower() or "prset" in line.lower()) and \
                            counter == 2:
                self._partitions.parse_nexus_model(line)

        # If no _partitions have been added during the parsing of the nexus
        # file, set a single partition
        if self._partitions.partitions == OrderedDict():
            self._partitions.add_partition(self.name, self.locus_length,
                                           file_name=self.path)

        fh.close()

        # Checks the size consistency of the alignment
        if len(set(size_list)) > 1:
            self.e = AlignmentUnequalLength()

    def _read_stockholm(self):
        """Alignment parser for stockholm format.

        Parses a stockholm alignment file and stores taxa and sequence data in
        the database.

        See Also
        --------
        read_alignment
        """

        fh = open(self.path)

        # Variable storing the lenght of each sequence
        size_list = []
        idx = 0

        # Skip first header line
        next(fh)

        for line in fh:
            # Skip header and comments
            if line.startswith("#"):
                pass
            # End of file
            elif line.strip() == "//":
                break
            # Parse sequence data
            elif line.strip() != "":

                taxa = line.split()[0].split("/")[0]
                taxa = self.rm_illegal(taxa)

                try:
                    sequence = line.split()[1].strip().lower()
                    # Evaluate missing data symbol if undefined
                    self._eval_missing_symbol(sequence)
                except IndexError:
                    sequence = ""

                size_list.append(len(sequence))

                self._insert_data(idx, taxa, sequence)

                self._taxa_idx[taxa] = idx

        self.locus_length = len(sequence)
        self._partitions.set_length(self.locus_length)

        # Updating _partitions object
        self._partitions.add_partition(self.name, self.locus_length,
                                       file_name=self.path)

        fh.close()

        # Checks the size consistency of the alignment
        if len(set(size_list)) > 1:
            self.e = AlignmentUnequalLength()

    def read_alignment(self):
        """Main alignment parser method.

        This is the main alignment parsing method that is called when the
        `Alignment` object is instantiated with a file path as the argument.
        Given the alignment format automatically detected in `__init__`,
        it calls the specific method that parses that alignment format.
        After the execution of this method, all attributes of the class will
        be set and the full range of methods can be applied.
        """

        parsing_methods = {
            "phylip": self._read_phylip,
            "fasta": self._read_fasta,
            "loci": self._read_loci,
            "nexus": self._read_nexus,
            "stockholm": self._read_stockholm
        }

        parsing_methods[self.input_format]()

        # If the missing data symbol could not be evaluated during alignment
        # parsing, set the defaults
        default_missing = {"DNA": "n", "Protein": "x"}
        if not self.sequence_code[1]:
            self.sequence_code[1] = default_missing[self.sequence_code[0]]

        if len(self._taxa_idx) != len(set(self._taxa_idx.keys())):
            duplicate_taxa = self.duplicate_taxa(self._taxa_idx.keys())
            self.e = DuplicateTaxa("The following taxa were duplicated in"
                                   " the alignment: {}".format(
                "; ".join(duplicate_taxa)))



    def remove_taxa(self, taxa_list_file, mode="remove"):
        """ Removes taxa from the `Alignment` object.

        Removes taxa based on the  `taxa_list_file` argument from the
        `Alignment` object.

        This method supports a list or path to CSV file as the
        `taxa_list_file` argument. In case of CSV path, it parses the CSV
        file into a list. The CSV file should be a simple text file containing
        taxon names in each line.

        This method also supports two removal modes.
            - remove: removes the specified taxa
            - inverse: removes all but the specified taxa

        Parameters
        ----------
        taxa_list_file : list or string
            List containing taxa names or string with path to a CSV file
            containig the taxa names.
        mode : {"remove", "inverse"}, optional
            The removal mode (default is remove).
        """

        tx_idx = self.get_taxa_idx()

        def remove(list_taxa):
            for tx in list_taxa:
                self.cur.execute(
                    "DELETE FROM alignment_data "
                    "WHERE txId=? AND aln_idx=?",
                    (tx_idx[tx], self.db_idx))
                del tx_idx[tx]

            self.store_taxa_idx(self.temp_dir, tx_idx)

        def inverse(list_taxa):
            for tx in list(tx_idx.keys()):
                if tx not in list_taxa:
                    self.cur.execute(
                        "DELETE FROM alignment_data "
                        "WHERE txId=? AND aln_idx=?",
                        (tx_idx[tx], self.db_idx))
                    del tx_idx[tx]

            self.store_taxa_idx(self.temp_dir, tx_idx)

        # Checking if taxa_list is an input csv file:
        try:
            file_handle = open(taxa_list_file[0])
            taxa_list = self.read_basic_csv(file_handle)

        # If not, then the method's argument is already the final list
        except (IOError, IndexError):
            taxa_list = taxa_list_file

        # Filter taxa_list for the taxa that are actually present in this
        # Alignment object
        taxa_list = [x for x in taxa_list if x in tx_idx]

        if mode == "remove":
            remove(taxa_list)
        if mode == "inverse":
            inverse(taxa_list)

    def change_taxon_name(self, old_name, new_name):
        """Changes the name of a particular taxon.

        Parameters
        ----------
        old_name : str
            Original taxon name.
        new_name : str
            New taxon name.
        """

        tx_idx = self.get_taxa_idx()

        # Change in taxa_list
        if old_name in tx_idx:
            # Change in the database
            self.cur.execute("UPDATE [{}] SET taxon=? WHERE txId=?".format(
                self.db_idx), (new_name, tx_idx[old_name]))
            # Change in taxa_index
            tx_idx[new_name] = tx_idx[old_name]
            del tx_idx[old_name]

        self.store_taxa_idx(self.temp_dir, tx_idx)

    def _check_partitions(self, partition_obj):
        """Consistency check of the `partition_obj` for `Alignment` object.

        Performs checks to ensure that a `partition_obj` is consistent
        with the alignment data of the `Alignment` object.

        Parameters
        ----------
        partition_obj : trifusion.process.data.Partitions
            Partitions object.

        See Also
        --------
        trifusion.process.data.Partitions
        """

        # Checks if total lenght of _partitions matches the lenght of the
        # current alignment

        if partition_obj.counter != self.locus_length:
            return process.data.InvalidPartitionFile("Partitions in partition"
                   " file are inconsistent with current alignment. Last "
                   "partition range: {}; Alignmet length: {}".format(
                partition_obj.counter, self.locus_length
            ))

    def set_partitions(self, partitions):
        """Updates the `partition` attribute.

        Parameters
        ----------
        partitions : trifusion.process.data.Partitions
            Partitions object.
        """

        # Checks partition's consistency
        er = self._check_partitions(partitions)

        if isinstance(er, process.data.InvalidPartitionFile):
            return er
        else:
            self.store_partitions(self.temp_dir, partitions)

    taxa_idx = property(get_taxa_idx)
    partitions = property(get_partitions)


class AlignmentList(Base):
    """Main interface for groups of `Alignment` objects.

    The `AlignmentList` class can be seen as a group of `Alignment` objects
    that is treated as a single cohesive data set. It is the main class
    used to perform alignment operations on TriFusion and TriSeq programs,
    even when there is only a single alignment.

    An instance can be created by simply providing a list of paths to
    alignment files, which can be empty (and alignments added later). It
    is recommended that the path to the sqlite database be specified
    via the `sql_db` argument, but this is optional.

    Individual alignment files can be provided when instantiating the class
    or added later, and each is processed and stored as an `Alignment` object.
    The `AlignmentList` object provides attributes that are relative to the
    global data set. For instance, each `Alignment` object may have their
    individual and unique list of taxa, but the corresponding attribute in the
    `AlignmentList` object contains the list of all taxa that are found
    in the `Alignment` object list.
    
    Several methods of the `Alignment` class are also present in this class
    with the same usage (e.g. `collapse`). These are basically wrappers that
    apply the method to all `Alignment` objects and may have some
    modifications adapted to multiple alignments.
    
    Other methods are exclusive and desgined to deal with multiple alignments,
    such as concatenation (`concatenate`), filter by taxa (`filter_by_taxa`),
    etc.

    Plotting methods are exclusive to `AlignmentList`, even the
    ones that are focused on single alignments. In this case, the
    `Alignment` object can be specified and processed individually directly
    from this class.

    Parameters
    ----------
    alignment_list : list
        List with paths to alignment files. Can be empty.
    sql_db : str, optional
        Path where the sqlite3 database will be created.
    db_cur : sqlite3.Cursor, optional
        Provide a database Cursor object in conjunction with the Connection
        object (`db_con`) to connect to an existing database.
    db_con : sqlite3.Connection, optional
        Provide a database Connection object in conjunction with the Cursor
        object (`db_cur`) to connect to an existing database.
    pbar : ProgressBar, optional
        A ProgressBar object used to log the progress of TriSeq execution.

    Attributes
    ----------
    con : sqlite3.Connection
        Connection object of the sqlite database.
    cur : sqlite3.Cursor
        Cursor object of the sqlite database.
    sql_path : str
        Path to sqlite3 database file.
    alignments : collections.OrderedDict
        Stores the 'active' `Alignment objects` as values, with the
        corresponding key being `Alignment.path`.
    shelve_alignments : collections.OrderedDict
        Stores the 'inactive' or 'shelved' `Alignment` objects.
        `AlignmentList` methods will operate only on the `alignments`
        attribute, unless explicitly stated otherwise. The key:value is the
        same as in the `alignments` attribute.
    bad_alignments : list
        Stores the `Alignment.name` attribute of badly formatted alignments.
    duplicate_alignments : list
        Stores the `Alignment.name` attribute of duplicated alignment paths.
    non_alignments : list
        Stores the `Alignment.name` attribute of sequence sets of unequal
        length.
    taxa_names : list
        Stores the name of all taxa across the `Alignment` object list.
    shelved_taxa : list
        Stores the 'inactive' or 'shelved' taxa. `AlignmentList` methods
        will operate only on the `taxa_names` attribute, unless explicitly
        stated otherwise.
    path_list : list
        List of `Alignment.paths` for all `Alignment` objects.
    filtered_alignments : collections.OrderedDict
        Ordered dictionary of four key:values used to count the number
        of `Alignment` objects that were filtered by four operations. The
        keys are {"By minimum taxa", "By taxa", "By variable sites",
        "By informative sites"}.
    sequence_code : tuple
        Contains information on (<sequence type>, <missing data symbol>),
        e.g. ("Protein", "x").
    gap_symbol : str
        Symbol used to denote gaps in the alignment (default is '-').
    summary_stats : dict
        Dictionary that stores several summary statistics that are calculated
        once the `get_summary_stats` method is executed.
    summary_gene_table : pandas.DataFrame
        DataFrame containing the summary statistics for each `Alignment`
        object. Also populated when the `get_summary_stats` method is
        executed.
    temporary_tables : list
        List of active database tables associated with the `AlignmentList`
        instance.
    partitions : trifusion.process.data.Partitions
        Partitions object that refers to the total `AlignmentList`.
    """

    def __init__(self, alignment_list, sql_db=None, db_cur=None, db_con=None,
                 pbar=None):

        # Create connection and cursor for sqlite database
        # If `db_cur` and `db_con` are both provided, setup the database
        # connection and ignore `sql_db`
        if db_cur and db_con:
            self.con = db_con
            """
            Database Connection object
            """
            self.cur = db_cur
            """
            Database Cursor object
            """
            self.sql_path = sql_db
        # If only `sql_db` is provided, set it as the attribute
        elif sql_db:
            self.sql_path = sql_db
        # If no database information is provided, set a temporary database
        # file in the current working directory
        else:
            self.sql_path = "trifusion.sqlite3"
            """
            Path to sqlite3 database file
            """

        self.master_table = "alignment_data"
        """Name of the table with the original alignment data"""

        if not db_cur and not db_con:
            self.con = sqlite3.connect(self.sql_path, check_same_thread=False)
            self.cur = self.con.cursor()
            self.cur.execute("PRAGMA synchronous = OFF")

        if not self._table_exists(self.master_table):
            self._create_table(self.master_table,
                               index=("main_idx", "aln_idx"))

        self.alignments = OrderedDict()
        """
        Stores the "active" `Alignment` objects for the current
        `AlignmentList`. Keys will be the `Alignment.path` for quick lookup of
        `Alignment` object values
        """

        self.all_alignments = OrderedDict()
        """
        Stores all the alignments of the AlignmentList instance, whether they
        are active or inactive. Do not change unless adding new alignments
        or clearing AlignmentList instance.
        """

        self.alignment_idx = OrderedDict()
        """
        OrderedDict that that maps the "aln_file" column value of the sqlite
        database, to the path/filename  of the corresponding alignment.
        """

        self.shelved_idx = []
        """
        List with shelved alignment idx.
        """

        self._idx = 1

        self.bad_alignments = []
        """
        Attribute list that stores the Alignment.name attribute of badly
        formatted alignments.
        """

        self.duplicate_alignments = []
        """
        Attribute list that stores duplicate Alignment.name.
        """

        self.non_alignments = []
        """
        Attribute list that stores sequence sets of unequal length.
        """

        self.taxa_names = []
        """
        List with the name of the taxa included in the AlignmentList object
        """

        self.shelved_taxa = []
        """
        List with non active taxa
        """

        self.path_list = []
        """
        List of `Alignment.path`
        """

        self.filtered_alignments = OrderedDict([("By minimum taxa", None),
                                               ("By taxa", None),
                                               ("By variable sites", None),
                                               ("By informative sites", None)])
        """
        Records the number of filtered alignments by each individual filter type
        """

        self.sequence_code = None
        """
        Tuple with the AlignmentList sequence code. Either ("DNA", "n") or
        ("Protein", "x")
        """
        self.gap_symbol = "-"

        self.summary_stats = {"genes": 0, "taxa": 0, "seq_len": 0, "gaps": 0,
                              "avg_gaps": [], "missing": 0, "avg_missing": [],
                              "variable": 0, "avg_var": [], "informative": 0,
                              "avg_inf": []}
        """
        Dictionary with summary statistics for the active alignments
        """

        columns = ["genes", "nsites", "taxa", "var", "inf", "gap", "missing"]
        self.summary_gene_table = pd.DataFrame(columns=columns)
        """
        Dictionary with summary statistics for each alignment.
        """

        self.format_ext = {"ima2": ".txt",
                            "mcmctree": "_mcmctree.phy",
                            "phylip": ".phy",
                            "nexus": ".nex",
                            "fasta": ".fas",
                            "stockholm": ".stockholm",
                            "gphocs": ".txt"}
        """Dictionary that stores the suffix for each output format"""

        self.temporary_tables = []
        """
        Lists the currently active tables. This is mainly used for the
        conversion of the consensus alignments into a single Alignment object
        """

        self.interleave_data = False

        self.partition_data = False

        self.partition_count = {}

        self.size = 0
        """
        Integer with the size of the complete data set
        """

        # Set _partitions object
        self.partitions = Partitions()

        # if type(alignment_list[0]) is str:
        if alignment_list:

            self.add_alignment_files(alignment_list, pbar=pbar)

    def __iter__(self):
        """Iterator behavior for `AlignmentList`

        Returns an iterator over the 'active' `Alignment` objects.

        Returns

        _ : iter
            Iterable over `alignments.values()`
        """
        return iter(self.alignments.values())

    def iter_alignments(self, table_name=None, include_txid=False):

        table_name = table_name if table_name else self.master_table

        # Check if table exists and is not empty. In any of these conditions,
        # fallback to the master table
        try:
            if not self.cur.execute(
                "SELECT * FROM {}".format(table_name)).fetchone():
                table_name = self.master_table
        except sqlite3.OperationalError:
            table_name = self.master_table

        try:

            lock.acquire(True)

            for txId, taxon, seq, aln_idx in self.cur.execute(
                    "SELECT txId, taxon, seq, aln_idx "
                    "FROM [{}] "
                    "WHERE aln_idx NOT IN ({})".format(
                        table_name,
                        ", ".join([str(x) for x in self.shelved_idx]))):
                if taxon not in self.shelved_taxa:
                    if include_txid:
                        yield txId, taxon, seq, aln_idx
                    else:
                        yield taxon, seq, aln_idx

        finally:
            lock.release()

    def iter_columns(self, table_name=None, aln_idx=None):

        table_name = table_name if table_name else self.master_table

        # Check if table exists and is not empty. In any of these conditions,
        # fallback to the master table
        try:
            if not self.cur.execute(
                    "SELECT * FROM {}".format(table_name)).fetchone():
                table_name = self.master_table
        except sqlite3.OperationalError:
            table_name = self.master_table

        try:

            lock.acquire(True)

            query = "SELECT " \
                    "GROUP_CONCAT(substr(seq, {pos}, 100000)), " \
                    "aln_idx " \
                    "FROM [{tb}] " \
                    "WHERE {cond} " \
                    "GROUP BY aln_idx"

            if aln_idx:
                cond = "aln_idx={} ".format(aln_idx)
            else:
                shelved = ", ".join([str(x) for x in self.shelved_idx])
                cond = "aln_idx NOT in ({})".format(shelved)

            for p in xrange(0, self.size, 100000):
                for res in ((x.split(","), y) for x, y in self.cur.execute(
                        query.format(pos=p, tb=table_name, cond=cond))):
                    for col in itertools.izip(*res[0]):
                        yield col, res[1]

        finally:
            lock.release()

    def _create_table(self, table_name, index=None, cur=None):
        """Creates a new table in the database.

        Convenience method that creates a new table in the sqlite database.
        It accepts a custom Cursor object, which overrides the `cur`
        attribute. Optionally, an index can also be specified.

        Parameters
        ----------
        table_name : str
            Name of the table to be created.
        index : list, optional
            Provide a list with [<name of index>, <columns to be indexed>].
            (e.g., ["concindex", "txId"]).
        cur : sqlite3.Cursor, optional
            Custom Cursor object used to query the database.
        """

        if not cur:
            cur = self.cur

        cur.execute("CREATE TABLE [{}]("
                    "txId INT,"
                    "taxon TEXT,"
                    "seq TEXT,"
                    "aln_idx INT)".format(table_name))

        if index:
            cur.execute("CREATE INDEX {} ON [{}]({})".format(
                index[0], table_name, index[1]))

    def _table_exists(self, table_name, cur=None):
        """ Checks if a table exists in the database.

        Convenience method that checks if a table exists in the database.
        It accepts a custom Cursor object, which overrides the `cur`
        attribute.

        Parameters
        ----------
        table_name : str
            Name of the table.
        cur: sqlite3.Cursor, optional
            Custom Cursor object used to query the database.

        Returns
        -------
        res : list
            List with the results of a query for 'table' type with
            `db_idx` name. Is empty when the table does not exist.

        Notes
        -----
        This returns a list that will contain one item if the table exists
        in the database. If it doesn't exist, returns an empty list.
        Therefore, this can be used simply like::

            if self._table_exists("my_table"):
                 # Do stuff
        """

        if not cur:
            cur = self.cur

        return cur.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND"
            " name='{}'".format(table_name)).fetchall()

    def _get_idx(self, aln_path):
        """Returns the `align_idx` for a given `aln_path`.

        Returns the alignment identifier number in the sqlite database
        (`aln_idx`) for a given alignment path (`aln_path`).

        Parameters
        ----------
        aln_path

        Returns
        -------

        """

        try:
            idx = [idx for idx, path in self.alignment_idx.items() if
                   path == aln_path][0]
        except IndexError:
            idx = None

        return idx

    @staticmethod
    def _check_killswitch(ns):

        if ns:
            if ns.stop:
                raise KillByUser("")

    def _set_pipes(self, ns=None, pbar=None, total=None, msg=None,
                   ignore_sa=False):
        """Setup of progress indicators for both GUI and CLI task executions.

        This handles the setup of the objects responsible for updating
        the progress of task's execution of both TriFusion (GUI) and
        TriSeq (CLI). At the beginning of any given task, these objects
        can be initialized by providing either the Namespace object (`ns`)
        in the case of TriFusion, or the ProgressBar object (`pbar`), in
        the case of TriSeq. Along with one of these objects, the expected
        `total` of the progress should also be provided. The `ns` and
        `pbar` objects are updated at each iteration of a given task,
        and the `total` is used to get a measure of the progress.

        Optionally, a message can be also provided for the Namespace object
        that will be used by TriFusion.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.
        total : int
            Expected total of the task's progress.
        msg : str, optional
            A secondary message that appears in TriFusion's progress dialogs.

        Notes
        -----
        The progress of any given task can be provided by either an
        `Alignment` or `AlignmentList` instance. Generally, the tasks follow
        the progress of the `AlignmentList` instance, unless that instance
        contains only one `Alignment` object. In that case, the progress
        information is piped from the `Alignment` instance. For that reason,
        and for the Namespace (`ns`) object only, this method checks if
        there is only one 'active' `Alignment` object. If yes, set the
        `Namepsace.sa` flag to True, so that the progress indicators
        are piped from the `Alignment` object instead.

        Examples
        --------
        Start a progress counter  for a task that will make 100 iterations::

            self._set_pipes(ns=ns_obj, pbar=pbar_obj, total=100,
                            msg="Some message")
        """

        if ns:
            if ns.stop:
                raise KillByUser("")

            if len(self.alignments) == 1 and not ignore_sa:
                ns.sa = True
            else:
                ns.sa = False
                ns.total = total
                ns.counter = 0
        if pbar:
            pbar.max_value = total

    @staticmethod
    def _update_pipes(ns=None, pbar=None, value=None, msg=None,
                      ignore_sa=False):
        """Update progress indicators for both GUI and CLI task executions.

        This method provides a single interface for updating the progress
        objects `ns` or `pbar`, which should have been initialized at the
        beginning of the task with the `_set_pipes` method.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.
        value : int
            Value of the current progress index
        msg : str, optional
            A secondary message that appears in TriFusion's progress dialogs.

        Notes
        -----
        The progress of any given task can be provided by either an
        `Alignment` or `AlignmentList` instance. Generally, the tasks follow
        the progress of the `AlignmentList` instance, unless that instance
        contains only one `Alignment` object. In that case, the progress
        information is piped from the `Alignment` instance. For that reason,
        and for the Namespace (`ns`) object only, if the `Namespace.sa`
        attribute is set to True, then to not update the progress. Progress
        is being piped directly from the `Alignment` object.

        Examples
        --------
        Update the counter in an iteration of 100::

            for i in range(100):
                self._update_pipes(ns=ns_obj, pbar=pbar_obj, value=i,
                                   msg="Some string")
        """
        if ns:
            if ns.stop:
                raise KillByUser("")
            if not ns.sa or ignore_sa:
                ns.msg = msg
                ns.counter = value

        if pbar:
            pbar.update(value)

    @staticmethod
    def _reset_pipes(ns):
        """Reset progress indicators for both GUI and CLI task executions.

        This should be done at the end of any task that initialized the
        progress objects, but it only affects the Namespace object. It
        resets all Namespace attributes to None.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        """

        if ns:
            if ns.stop:
                raise KillByUser("")
            ns.counter = ns.total = ns.msg = ns.sa = None

    def aln_names(self):
        """Returns list of basenames of `alignments` file paths."""

        return sorted([basename(x) for x in self.alignments])

    def resume_database(self):
        """Reconnects to the sqlite database.

        The connection to the database can be closes or even severed when
        passing objects between threads. This methods allows the
        reconnection of the database for the `AlignmentList` and
        *all* (even the shelved ones) `Alignment` objects.
        """

        self.con = sqlite3.connect(self.sql_path, check_same_thread=False)
        self.cur = self.con.cursor()

        for aln in self.all_alignments.values():
            aln.cur = self.cur
            aln.con = self.con

    def set_database_connections(self, cur, con):
        """Provides Connection and Cursor to `Alignment` objects.

        Sets the database connections manually for *all* (even the
        shelved ones) `Alignment` objects.

        Parameters
        ----------
        cur : sqlite3.Cursor
            Provide a database Cursor object.
        con : sqlite3.Connection
            Provide a database Connection object.
        """

        self.cur = cur
        self.con = con

        for aln in self.all_alignments.values():
            aln.cur = cur
            aln.con = con

    def get_tables(self):
        """Return list with `db_idx` of *all* `Alignment` objects.

        Returns
        -------
        _ : list
            List with `Alignment.db_idx`.

        """

        return [x.table_name for x in
                self.all_alignments.values()]

    def remove_tables(self, preserve_tables=None, trash_tables=None):
        """Drops tables from the database.

        This method can be used to drop tables from the database.

        If no arguments are provided it will remove *all*  database tables
        associated with the `AlignmentList` object.  If the
        `preserve_tables` argument is provided, the table names present in
        this list will not be dropped. If `trash_tables` is provided,
        only the tables names specified in that list are dropped. If
        both are provided, `trash_tables` will take precedence.

        Parameters
        ----------
        preserve_tables : list, optional
            If provided, the table names in this list will NOT be dropped.
        trash_tables : list, optional
            If provided, only the table names in this list will be dropped.
            Takes precedence over `preserve_tables` if both are provided.
        """

        tables = self.cur.execute(
            "SELECT name FROM sqlite_master WHERE type='table';").fetchall()

        for tb in [x[0] for x in tables if x[0] != self.master_table]:
            self.cur.execute("DROP TABLE [{}]".format(tb))

    def clear_alignments(self):
        """Clears all attributes and data from the `AlignmentList` object."""

        # Drop active databases of the AlignmentList instance, namely consensus
        # and concatenation, if they exist
        for table in self.temporary_tables:
            self.cur.execute("DROP TABLE [{}]".format(table))

        self.cur.execute("DELETE FROM [{}]".format(self.master_table))

        self.alignments = OrderedDict()
        self.all_alignments = OrderedDict()
        self.alignment_idx = OrderedDict()
        self.bad_alignments = []
        self.duplicate_alignments = []
        self.non_alignments = []
        self.taxa_names = []
        self.shelved_taxa = []
        self.path_list = []
        self._idx = 1
        self.size = 0
        self.filtered_alignments = OrderedDict([("By minimum taxa", None),
                                                ("By taxa", None),
                                                ("By variable sites", None),
                                                ("By informative sites", None)])
        self.sequence_code = None
        self.summary_stats = {"genes": 0, "taxa": 0, "seq_len": 0, "gaps": 0,
                              "avg_gaps": [], "missing": 0, "avg_missing": [],
                              "variable": 0, "avg_var": [], "informative": 0,
                              "avg_inf": []}
        columns = ["genes", "nsites", "taxa", "var", "inf", "gap", "missing"]
        self.summary_gene_table = pd.DataFrame(columns=columns)
        self.temporary_tables = []
        self.partitions = Partitions()

    def _reset_summary_stats(self):
        """Resets the `summary_stats` attribute."""

        self.summary_stats = {"genes": 0, "taxa": 0, "seq_len": 0, "gaps": 0,
                              "avg_gaps": [], "missing": 0, "avg_missing": [],
                              "variable": 0, "avg_var": [], "informative": 0,
                              "avg_inf": []}

    def update_active_alignments(self, aln_list=None, all_files=False,
                                 ns=None, pbar=None):
        """Sets the active alignments.

        Sets the 'active' alignments stored in the `alignments` attribute
        and the 'inactive' alignments stored in the `shelve_alignments`
        attribute, based on the `aln_list` argument. Regardless of the
        current 'active' and 'inactive' alignments, calling this method
        with the `aln_list` argument will set only those as the 'active'
        alignments.

        Alternatively, the `all_files` bool argument can be provided to
        set all alignments as 'active'.

        After updating the `alignments`and `shelve_alignments` attributes,
        the `taxa_names` attribute is also updated.

        Parameters
        ----------
        aln_list : list
            List of `Alignment.path` that will be set as 'active'.
        all_files : bool
            If True, all alignments become active. Ignores `aln_list`.
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.
        """
        #
        # self._set_pipes(ns, pbar, total=len(self.alignments.keys() +
        #                                     self.shelve_alignments.keys()))

        if all_files:

            self.alignments = OrderedDict(
                (p, aln) for p, aln in self.all_alignments)
            self.shelved_idx = []

        elif aln_list is not None:

            self.alignments = OrderedDict(
                (p, aln) for p, aln in self.all_alignments.items()
                if p in aln_list)
            self.shelved_idx = [idx for idx, aln in self.alignment_idx.items()
                                if aln.path not in self.alignments]

        # Update taxa names
        self.taxa_names = self._get_taxa_list()

        #Update size
        self.size = sum((x.locus_length for x in self.alignments.values()))

    def update_active_alignment(self, aln_name, direction):
        """Updates the 'active' status of a single alignment.

        Changes the 'active' status of a single alignment in a specified
        direction.

        After updating the `alignments`and `shelve_alignments` attributes,
        the `taxa_names` attribute is also updated.

        Parameters
        ----------
        aln_name : str
            `Alignment.name` to update.
        direction : str
            Direction where `aln_name` will be moved. Can be
            {"shelve", "active"}.
        """

        if direction == "shelve":
            self.shelved_idx.append(self._get_idx(
                self.alignments[aln_name]))
            del self.alignments[aln_name]

        else:
            self.alignments[aln_name] = self.all_alignments[aln_name]
            self.shelved_idx.remove(self._get_idx(self.alignments[aln_name]))

        # Update taxa names
        self.taxa_names = self._get_taxa_list()

        # Update size
        self.size = sum((x.locus_length for x in self.alignments.values()))

    def update_taxa_names(self, taxa_list=None, all_taxa=False):
        """Sets the active taxa.

        Sets the 'active' taxa stored in the `taxa_names` attribute and the
        'inactive' taxa stored in the `shelved_taxa` attribute, based on
        the `taxa_list` argument.

        Alternatively, the `all_taxa` bool argument can be provided to
        set all taxa as 'active'.

        Parameters
        ----------
        taxa_list : list
            List of taxon names that will be set as 'active'.
        all_taxa : bool
            If True, all alignments become active. Ignores `taxa_list`.
        """

        # Activate all taxa
        if all_taxa:
            for tx in self.shelved_taxa:
                self.taxa_names.append(tx)

            self.shelved_taxa = []

        # Activate only taxa specified by taxa_list
        elif taxa_list:

            for tx in self.taxa_names + self.shelved_taxa:
                if tx not in taxa_list and tx in self.taxa_names:
                    try:
                        self.taxa_names.remove(tx)
                        self.shelved_taxa.append(tx)
                    except ValueError:
                        pass

                elif tx in taxa_list and tx in self.shelved_taxa:
                    try:
                        self.taxa_names.append(tx)
                        self.shelved_taxa.remove(tx)
                    except ValueError:
                        pass

        # Update individual Alignment objects
        for aln_obj in self.alignments.values():
            aln_obj.shelve_taxa(self.shelved_taxa)

    def format_list(self):
        """Returns list of unique sequence types from `Alignment` objects.

        Returns
        -------
        _ : list
            List with unique sequence types for `Alignment` objects as
            strings.
        """

        return list(set([x.sequence_code[0] for x in
                         self.alignments.values() if x]))

    def _get_taxa_list(self, only_active=False):
        """Gets the global taxa list from all alignments.

        If called with no arguments, it returns the taxa list from *all*
        alignments, including the 'inactive' ones. The `only_active`
        argument can be set to True to return only the taxa list from the
        'active' files.

        Parameters
        ----------
        only_active : bool
            If True, only returns the taxa list from the 'active' alignments.

        Returns
        -------
        full_taxa : list
            List with taxon names as strings.
        """

        if only_active:
            full_taxa = list(set().union(*[x.taxa_idx.keys() for x in
                                           self.alignments.values()]))
        else:
            full_taxa = list(set().union(*[x.taxa_idx.keys() for x in
                                           self.all_alignments.values()]))

        return full_taxa

    def _get_filename_list(self):
        """Returns list with the `Alignment.name` of alignments.

        Returns
        -------
        _ : list
            List of `Alignment.name` of alignments.
        """
        return [alignment.name for alignment in self.alignments.values()]

    def set_partition_from_alignment(self, alignment_obj, reset=False):
        """Updates `_partitions` object with the provided `Alignment` object.

        Uses the `_partitions` of an `Alignment` object to set the `_partitions`
        attribute of `AlignmentList`

        Parameters
        ----------
        alignment_obj : trifusion.process.sequence.Alignment
            `Alignment` object.
        reset : bool
            If True, clears the current _partitions
        """
        
        if reset:
            self.partitions = Partitions()

        aln_parts = alignment_obj.get_partitions()

        # Update _partitions object
        # NOTE: The use_counter argument is set to False here so that when
        # the locus_range is provided, the counter does not interfere with the
        # ranges.
        if not aln_parts.is_single():
            for k, v in aln_parts:
                self.partitions.add_partition(
                    k, locus_range=v[0], codon=v[1],
                    use_counter=True, file_name=alignment_obj.path,
                    model_cls=aln_parts.models[k])
        else:
            self.partitions.add_partition(
                alignment_obj.name, use_counter=True,
                file_name=alignment_obj.path,
                length=alignment_obj.locus_length,
                model_cls=aln_parts.models[alignment_obj.name])

        self.size = self.partitions.counter

    def add_alignments(self, alignment_obj_list, ignore_paths=False):
        """Add a list of `Alignment` objects to the current `AlignmentList`.

        Adds a list of alignments, already as `Alignment` objects, to the
        current `AlignmentList`. This method performs several checks
        to assure that the `Alignment` objects being added are correct. If
        all checks out, updates all attributes of the class with the new
        alignments.

        The `ignore_paths` argument can be used to prevent `Alignment`
        objects with short paths (only basename, for instance) from being
        added if they have the same basename and another in the `alignments`
        dict.

        Parameters
        ----------
        alignment_obj_list : list
            List with `Alignment` objects.
        ignore_paths : bool
            If True, the `Alignment.path` of the new alignments are not
            checked with the loaded alignments.
        """

        for alignment_obj in alignment_obj_list:

            # Check for badly formatted alignments
            if isinstance(alignment_obj.e, InputError):
                self.bad_alignments.append(alignment_obj.path)
            elif isinstance(alignment_obj.e,
                            AlignmentUnequalLength):
                self.non_alignments.append(alignment_obj.path)

            # if not ignore_paths:
            if not ignore_paths:
                if alignment_obj.path in [x.path for x in
                                          self.alignments.values()]:
                    self.duplicate_alignments.append(alignment_obj.name)
                else:
                    # Get seq code
                    if not self.sequence_code:
                        self.sequence_code = alignment_obj.sequence_code

                    sql_path = self.sql_path if self.sql_path else "."

                    alignment_obj.store_taxa_idx(os.path.dirname(sql_path))
                    alignment_obj.store_partitions(os.path.dirname(sql_path))

                    self.alignments[alignment_obj.path] = alignment_obj
                    self.set_partition_from_alignment(alignment_obj)
            else:
                # Get seq code
                if not self.sequence_code:
                    self.sequence_code = alignment_obj.sequence_code

                sql_path = self.sql_path if self.sql_path else "."

                alignment_obj.store_taxa_idx(os.path.dirname(sql_path))
                alignment_obj.store_partitions(os.path.dirname(sql_path))

                self.alignments[alignment_obj.name] = alignment_obj
                self.set_partition_from_alignment(alignment_obj)
                # self.filename_list.append(alignment_obj.name)

        self.taxa_names = self._get_taxa_list()

    def add_alignment_files(self, file_name_list, pbar=None,
                            ns=None):
        """Adds a list of alignment files to the current `AlignmentList`.

        Adds a list of alignment paths to the current `AlignmentList`. Each
        path in the list is used to instantiate an `Alignment` object and
        several checks are performed to ensure that the alignments are
        correct and compliant with the other members of the `AlignmentList`
        object.

        Parameters
        ----------
        file_name_list : list
            List with paths to sequence alignment files.
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.
        """

        # Check for duplicates among current file list
        for f in [x for x, y in Counter(file_name_list).items() if y > 1]:
            self.duplicate_alignments.append(f)
            file_name_list.remove(f)

        # Check for duplicates between current file list and previous file list
        for i in set(self.path_list).intersection(set(file_name_list)):
            self.duplicate_alignments.append(i)
            file_name_list.remove(i)

        if ns:
            ns.progress = 0

        if pbar:
            pbar.max_value = len(file_name_list)

        for p, aln_path in enumerate(file_name_list):

            # Progress bar update for command line version
            if pbar:
                pbar.update(p + 1)

            if ns:
                ns.progress += 1
                ns.m = "Processing file {}".format(
                    basename(aln_path))

                if ns.stop:
                    raise KillByUser("Child thread killed by user")

            aln_obj = Alignment(aln_path, sql_cursor=self.cur,
                                db_idx=self._idx,
                                temp_dir=os.path.dirname(self.sql_path))

            if isinstance(aln_obj.e, InputError):
                self.bad_alignments.append(aln_obj.path)
            elif isinstance(aln_obj.e, AlignmentUnequalLength):
                self.non_alignments.append(aln_obj.path)
            elif isinstance(aln_obj.e, EmptyAlignment):
                self.bad_alignments.append(aln_obj.path)
            else:

                sql_path = self.sql_path if self.sql_path else "."

                aln_obj.store_taxa_idx(os.path.dirname(sql_path))
                aln_obj.store_partitions(os.path.dirname(sql_path))

                # Get seq code
                if not self.sequence_code:
                    self.sequence_code = aln_obj.sequence_code
                # Check for multiple sequence types. If True,
                # raise Exception
                elif self.sequence_code[0] != aln_obj.sequence_code[0]:
                    raise MultipleSequenceTypes("Multiple sequence "
                        "types detected: {} and {}".format(
                            self.sequence_code[0],
                            aln_obj.sequence_code[0]))

                self.all_alignments[aln_obj.path] = aln_obj
                self.alignments[aln_obj.path] = aln_obj
                self.set_partition_from_alignment(aln_obj)
                self.path_list.append(aln_obj.path)
                self.alignment_idx[self._idx] = aln_obj
                self._idx += 1

        if self.alignments:
            self.taxa_names = self._get_taxa_list()

    def retrieve_alignment(self, name):
        """Return `Alignment` object with a given `name`.

        Retrieves the `Alignment` object with a matching `name` attribute
        from `alignments` or `shelve_alignments`. Returns None if no match
        is found.

        Parameters
        ----------
        name : str
            `name` attribute of an `Alignment` object.

        Returns
        -------
        _ : trifusion.process.sequence.Alignment
            `Alignment` object.
        """

        if name in self.all_alignments:
            return self.all_alignments[name]
        else:
            return None

    def iter_alignment_files(self):
        """Returns an iterable on `Alignment.path` of each alignment.

        Returns
        -------
        _ : iter
            Iterable over the `Alignment.path` attribute of each alignment.
        """

        return iter(alignment.path for alignment in self.alignments.values())

    def write_taxa_to_file(self, file_name="Taxa_list.csv"):
        """Writes the taxa in `taxa_list` to a file.

        Parameters
        ----------
        file_name : string
            Path to the generated file (default is 'Taxa_list.csv').
        """

        output_handle = open(file_name, "w")

        for taxon in self.taxa_names:
            output_handle.write(taxon + "\n")

        output_handle.close()

    def concatenate(self, table_in="", table_out="", ns=None, pbar=None):
        """Concatenates alignments into a single `Alignment` object.

        Concatenates multiple sequence alignments creating a single `Alignment`
        object and the auxiliary Partitions object defining the _partitions
        of the concatenated alignment

        Parameters
        ----------
        table_in : string
            Name of database table containing the alignment data that is
            used for this operation.
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.

        Returns
        -------
        concatenated_alignment : trifusion.process.sequence.Alignment
            The `Alignment` object with the concatenated data and _partitions.

        Notes
        -----
        Since sqlite queries are quite expensive, the previous approach
        of querying the sequence for each taxa AND alignment table was
        dropped. In this approach, sqlite will do all the heavy lifting.
        First, a temporary table with the same definition as all alignment
        tables is populated with data from all taxa and alignments. It is
        important that and index is created on the new txId column, which
        will redifine the txId of individual alignments according to the
        global `taxa_names` attribute. In the first procedure, there
        will be only one query per alignment. When the temporary table is
        complete, some sql operations are used to group sequences from each
        taxon for all alignments and then concatenated them, all in a single
        query. This query returns the concatenated sequences, corrected txId
        and taxon information, ready to populate the final concatenation
        table. This approach is an order o magnitude faster than the
        previous one where thousands of queries could be performed.
        """

        def fill_missing_taxa(aln, cur):
            """Inserts missing sequences from current alignment in database

            Given an :class:`.Alignment` object, this will check for the
            taxa that are present in the overall
            :attr:`.AlignmentList.taxa_names` attribute, but not in the
            individual alignment. For each of those taxa, insert missing
            data sequences in the database.

            Parameters
            ----------
            aln : trifusion.process.sequence.Alignment
                :class:`.Alignment` object.
            cur : sqlite3.Cursor
                Cursor object of the sqlite database.
            """

            # Get list of missing taxa in the provided alignment
            missing_tx = set(self.taxa_names) - set(aln_obj.taxa_idx.keys())

            for tx in missing_tx:
                # Insert missing data in database
                cur.execute(
                    "INSERT INTO [{}] VALUES (?, ?, ?, ?)".format(
                        temp_table),
                    (taxa_idx[tx],
                     tx,
                     aln_obj.sequence_code[1] * aln_obj.locus_length,
                     1))

        # Set progress pipes
        self._set_pipes(ns, pbar, total=len(self.alignments))
        # Progress counter
        c = 1

        # Create table temporary table and create an index
        temp_table = ".concatenation"
        self._create_table(temp_table, index=["concidex", "txId"])

        # Variables that will store the taxa_list and _taxa_idx that will be
        # provided when instantiating the Alignment object
        taxa_idx = dict((tx, idx) for idx, tx in enumerate(self.taxa_names))

        # Create new cursor to insert data into database while the main
        # cursor returns the query
        conc_cur = self.con.cursor()

        # Defining empty alignment object
        aln_obj = None

        prev_idx = ""
        for taxon, seq, aln_idx in self.iter_alignments(table_in):

            # This happens when the alignment changes during the iteration.
            if aln_idx != prev_idx:

                # Update progress
                self._update_pipes(ns, pbar, value=c)
                c += 1

                # If prev_idx is already defined, it means this is the second
                # alignment during the iteration. In that case, fill the table
                # with missing data from the previous alignment.
                if prev_idx:
                    fill_missing_taxa(aln_obj, conc_cur)

                prev_idx = aln_idx
                aln_obj = self.alignment_idx[prev_idx]

            # Add data to temporary table
            conc_cur.execute("INSERT INTO [{}] VALUES (?, ?, ?, ?)".format(
                    temp_table), (taxa_idx[taxon], taxon, seq, 1))

        # Perform the fill missing taxa routine for the last alignment
        # Or when there is only one active alignment.
        if aln_idx:
            aln_obj = self.alignment_idx[aln_idx]
            fill_missing_taxa(aln_obj, conc_cur)

        # Setup final table that will have the concatenation
        table_out = table_out if table_out else self.master_table
        if self._table_exists(table_out):
            self.cur.execute("DROP TABLE [{}]".format(table_out))
        self._create_table(table_out, index=["conc_idx", "aln_idx"])

        # Reset progress information for next loop
        self._reset_pipes(ns)
        self._set_pipes(ns, pbar, total=len(self.taxa_names))

        # Here is where sqlite will do the heavy lifting of grouping sequences
        # from the same taxon across all alignments and then concatenating
        # them. It is important that the column that will be grouped has
        # an index (for perfomance and memory reasons). While seqs are
        # grouped by txId, the GROUP_CONCAT() method is used to return
        # concatenated strings.
        for p, (idx, tx, seq, aln_idx) in enumerate(conc_cur.execute(
                "SELECT txId, taxon, GROUP_CONCAT(seq, ''), aln_idx "
                "FROM [{}] "
                "GROUP BY txId".format(temp_table))):

            self._update_pipes(ns, pbar, value=p + 1,
                               msg="Concatenating taxon {}".format(tx))

            self.cur.execute("INSERT INTO [{}] VALUES (?, ?, ?, ?)".format(
                table_out), (idx, tx, seq, aln_idx))

        # Variable that will store the length of the concatenated alignment
        # and provided it when initializing the Alignment object
        locus_length = len(seq)

        # DROP the temporary table
        self.cur.execute("DROP TABLE [{}]".format(temp_table))

        self._idx += 1
        aln = Alignment(table_out, sql_cursor=self.cur,
                        taxa_idx=taxa_idx,
                        ignore_db_check=True, partitions=self.partitions,
                        sequence_code=self.sequence_code,
                        locus_length=locus_length,
                        db_idx=self._idx,
                        temp_dir=os.path.dirname(self.sql_path))

        # Reset alignment_idx attribute to reflect the single concatenated
        # alignment
        self.alignment_idx[1] = aln
        # Reset the shelved_idx. After the concatenation, the file
        # no longer takes effect
        self.shelved_idx = []

    def filter_min_taxa(self, min_taxa, ns=None, pbar=None):
        """Filters `alignments` by minimum taxa proportion.

        Filters `Alignment` objects based on a minimum taxa representation
        threshold. Alignments with less that the specified minimum taxa
        percentage will be moved to the `filtered_alignments` attribute.

        NOTE: Since this filtering is meant to be performed when executing
        the process operations, it will permanently change the AlignmentList
        object, which means both self.alignments and self._partitions. Not doing
        so and removing/adding the _partitions would create a great deal of
        conflicts that can be easily avoided by simply copying the
        AlignmentList object and modifying this object for the process
        execution

        Parameters
        ----------
        min_taxa : int
            Percentage of minimum allowed taxa representation (e.g. 50).
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.
        """

        self._set_pipes(ns, pbar, total=len(self.alignments))

        self.filtered_alignments["By minimum taxa"] = 0

        filtered_alns = []
        active_alns = []

        for p, (k, aln_obj) in enumerate(list(self.alignments.items())):

            self._update_pipes(
                ns, pbar, value=p + 1, msg="Evaluating file {}".format(
                    basename(aln_obj.name)))

            if len(aln_obj.taxa_idx) < \
                    (float(min_taxa) / 100.) * len(self.taxa_names):
                filtered_alns.append(aln_obj.path)
                self.filtered_alignments["By minimum taxa"] += 1
            else:
                active_alns.append(aln_obj.path)

        # Remove all _partitions at once. Much performance, such speed, wow.
        self.partitions.remove_partition(file_list=filtered_alns)
        # Update active alignments
        self.update_active_alignments(active_alns, ns=ns)

        self._reset_pipes(ns)

    def filter_by_taxa(self, taxa_list, filter_mode, ns=None, pbar=None):
        """Filters `alignments` if they contain or exclude certain taxa.

        Filters the `alignments` attribute by a given taxa list.
        The filtering may be performed to filter alignments
        depending on whether they include or exclude certain taxa.

        Parameters
        ----------
        taxa_list : list
            List of taxa names.
        filter_mode : str
            Filter mode. Can be "contain" or "exclude".
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.
        """

        self._set_pipes(ns, pbar, total=len(self.alignments))

        self.filtered_alignments["By taxa"] = 0

        # Support automatic file detection if taxa_list is a string
        if isinstance(taxa_list, str):
            try:
                file_handle = open(taxa_list)
                taxa_list = self.read_basic_csv(file_handle)
            except IOError:
                pass
        elif len(taxa_list) == 1:
            try:
                file_handle = open(taxa_list[0])
                taxa_list = self.read_basic_csv(file_handle)
            except IOError:
                pass

        # Stores Alignment.path to be filtered
        filtered_alns = []
        # Stpres Ailignment.name of active alignments
        active_alns = []

        for p, (k, aln_obj) in enumerate(list(self.alignments.items())):

            self._update_pipes(
                ns, pbar, value=p + 1, msg="Filtering file {}".format(
                    basename(aln_obj.name)))

            # Filter alignments that do not contain at least all taxa in
            # taxa_list
            if filter_mode.lower() == "contain":
                if set(taxa_list) - set(aln_obj.taxa_idx.keys()) != set():
                    filtered_alns.append(aln_obj.path)
                    self.filtered_alignments["By taxa"] += 1
                else:
                    active_alns.append(aln_obj.path)

            # Filter alignments that contain the taxa in taxa list
            if filter_mode.lower() == "exclude":
                if any((x for x in taxa_list
                        if x in aln_obj.taxa_idx.keys())):
                    filtered_alns.append(aln_obj.path)
                    self.filtered_alignments["By taxa"] += 1
                else:
                    active_alns.append(aln_obj.path)

        # Update _partitions
        self.partitions.remove_partition(file_list=filtered_alns)
        # Update active files
        self.update_active_alignments(active_alns, ns=ns)

        self.taxa_names = self._get_taxa_list(only_active=True)

        # If the resulting alignment is empty, raise an Exception
        if self.alignments == {}:
            if ns:
                ns.exception = "EmptyAlignment"
            raise EmptyAlignment("Alignment is empty after taxa filter")

        self._reset_pipes(ns)

    def filter_codon_positions(self, position_list, table_in=None,
                               table_out=None, ns=None, pbar=None):
        """Filters codon positions in each `Alignment` object.

        This wraps the execution of `filter_codon_positions` method for
        each `Alignment` object in the `alignments` attribute.

        Parameters
        ----------
        position_list : list
            List of three bool elements that correspond to each codon position.
            Ex. [True, True, True] will save all positions while
            [True, True, False] will exclude the third codon position
        table_in : string
            Name of database table containing the alignment data that is
            used for this operation.
        table_out : string
            Name of database table where the final alignment will be inserted
            (default is "filter").
        use_main_table : bool
            If True, both `table_in` and `table_out` are ignore and the main
            table `Alignment.db_idx` is used as the input and output
            table (default is False).
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.

        See Also
        --------
        Alignment.filter_codon_positions
        """

        def index(length, pos):
            """
            index generator
            """
            for _ in range(0, length, 3):
                for j in pos:
                    if j:
                        yield 1
                    else:
                        yield 0

        # Create temporary table
        temp_table = ".codonfilter"
        self._create_table(temp_table)

        # Set progress pipes
        self._set_pipes(ns, pbar, total=len(self.alignments))
        c = 1

        # Reset _partitions
        self.partitions = Partitions()

        # Create out table if does not exist
        if not self._table_exists(table_out):
            self._create_table(table_out)

        # Set temporary cursor to perform database changes while querying
        temp_cur = self.con.cursor()

        # Defining empty alignment object
        aln_obj = None

        prev_idx = ""
        for txId, taxon, seq, aln_idx in self.iter_alignments(
                table_in, include_txid=True):

            # This happens when the alignment changes during the iteration.
            if aln_idx != prev_idx:

                # Update the partition size of the last alignment
                if prev_idx:
                    aln_obj.locus_length = len(final_seq)
                    self.set_partition_from_alignment(aln_obj)

                aln_obj = self.alignment_idx[aln_idx]

                # Update progress
                self._update_pipes(ns, pbar, value=c,
                                   msg="Filtering file {}".format(
                                       aln_obj.name))
                c += 1

                prev_idx = aln_idx

            # Use itertools.compress to remove positions of the string
            # based on position list
            final_seq = "".join(list(itertools.compress(
                seq, index(aln_obj.locus_length, position_list))))

            # Add data to
            temp_cur.execute(
                "INSERT INTO [{}] VALUES (?, ?, ?, ?)".format(temp_table),
                (txId, taxon, final_seq, aln_idx))

        # Perform partition correction on the last loop element
        if aln_obj:
            aln_obj.locus_length = len(final_seq)
            self.set_partition_from_alignment(aln_obj)
        
        # If a previous table_out exist, replace with this new one
        if self._table_exists(table_out):
            self.cur.execute("DROP TABLE [{}];".format(table_out))

        # Replace table_out with collapsed table
        self.cur.execute("ALTER TABLE [{}] RENAME TO [{}]".format(
            temp_table, table_out))

        self._reset_pipes(ns)

    def _filter_terminals(self, table_in, table_out, ns=None):
        """

        Parameters
        ----------
        table_in
        table_out
        ns

        Returns
        -------

        """

        # Set progress pipes
        self._set_pipes(ns, None, total=len(self.alignments))

        # Create temporary table
        temp_table = ".filterterminals"
        self._create_table(temp_table)

        # Create temporary cursor to edit database while querying
        temp_cur = self.con.cursor()

        prev_idx = ""
        c = 1
        for txId, taxon, seq, aln_idx in self.iter_alignments(
                table_in, include_txid=True):

            # This happens when the alignment changes during the iteration
            if aln_idx != prev_idx:

                # Update progress
                self._update_pipes(ns, None, value=c,
                                   msg="Filtering terminals")
                c += 1

                prev_idx = aln_idx

            # Condition where the sequence only has gaps
            if not seq.strip("-"):
                seq = self.sequence_code[1] * len(seq)
                temp_cur.execute(
                    "INSERT INTO [{}] VALUES (?, ?, ?, ?)".format(
                        temp_table), (txId, taxon, seq, aln_idx))

                continue

            seq = list(seq)
            counter, reverse_counter = 0, -1

            while seq[counter] == "-":
                seq[counter] = self.sequence_code[1]
                counter += 1

            while seq[reverse_counter] == "-":
                seq[reverse_counter] = self.sequence_code[1]
                reverse_counter -= 1

            temp_cur.execute(
                "INSERT INTO [{}] VALUES (?, ?, ?, ?)".format(temp_table),
                (txId, taxon, "".join(seq), aln_idx))

        # Check if input and output tables are the same. If they are,
        # drop the old table and replace with this new one
        if self._table_exists(table_out):
            self.cur.execute("DROP TABLE [{}]".format(table_out))

        self.cur.execute("ALTER TABLE [{}] "
                         "RENAME TO [{}]".format(temp_table, table_out))

    def _filter_columns(self, gap_threshold, missing_threshold, table_in,
                        table_out, ns=None, pbar=None):

        # Create pipes
        self._set_pipes(ns, pbar, total=self.size)

        # Create temporary table
        temp_table = ".filtercolumns"
        self._create_table(temp_table)

        # Stores the binary lists that will be used to compress the final
        # sequences
        filtered_res = {}

        # Create temporary cursor to edit database while querying
        temp_cur = self.con.cursor()

        aln_obj = None

        prev_idx = ""
        for p, (column, aln_idx) in enumerate(self.iter_columns(table_in)):

            # This happes when the alignment changes during the iteration
            if aln_idx != prev_idx:

                # Add the binary list from the previous alignment
                if prev_idx:
                    filtered_res[prev_idx] = filtered_cols
                
                aln_obj = self.alignment_idx[aln_idx]
                taxa_number = len(aln_obj.taxa_idx)

                # Reset binary list for next alignment
                filtered_cols = []

                prev_idx = aln_idx

                # Update progress
                self._update_pipes(ns, pbar, value=p + 1, 
                                   msg="Filtering columns")

            cadd = column.count
            
            # Calculating metrics
            gap_proportion = (float(cadd("-")) /
                              taxa_number) * float(100)
            missing_proportion = (float(cadd(self.sequence_code[1])) /
                                  taxa_number) * float(100)
            total_missing_proportion = gap_proportion + missing_proportion

            if gap_proportion <= gap_threshold and \
                    total_missing_proportion <= missing_threshold:

                filtered_cols.append(1)

            else:
                filtered_cols.append(0)

        # Add binary list from last alignment
        if aln_obj:
            filtered_res[aln_idx] = filtered_cols

        self._set_pipes(ns, pbar, total=len(self.alignments))
        c = 1

        aln_obj = None

        # Now we compress the sequences with the binary lists from the
        # previous iteration
        prev_idx = ""
        for p, (txId, taxon, seq, aln_idx) in enumerate(
                self.iter_alignments(table_in, include_txid=True)):

            if prev_idx != aln_idx:

                if prev_idx:
                    # Update partition size
                    aln_obj.locus_length = len(final_seq)
                    self.set_partition_from_alignment(aln_obj)

                self._update_pipes(ns, pbar, value=c,
                                   msg="Compressing sequences")
                c += 1

                aln_obj = self.alignment_idx[aln_idx]

                prev_idx = aln_idx

            final_seq = "".join(compress(seq, filtered_res[aln_idx]))

            temp_cur.execute(
                "INSERT INTO [{}] VALUES (?, ?, ?, ?)".format(temp_table),
                (txId, taxon, final_seq, aln_idx))

        if aln_obj:
            aln_obj.locus_length = len(final_seq)
            self.set_partition_from_alignment(aln_obj)
            
        #Update size
        self.size = sum((x.locus_length for x in self.alignments.values()))

        # Check if input and output tables are the same. If they are,
        # it means that the output table already exists and is being
        # updated
        if self._table_exists(table_out):
            self.cur.execute("DROP TABLE [{}];".format(table_out))

        self.cur.execute("ALTER TABLE [{}] RENAME TO [{}]".format(
            temp_table, table_out))

        self._reset_pipes(ns)

    def filter_missing_data(self, gap_threshold, missing_threshold,
                            table_in=None, table_out=None, ns=None,
                            pbar=None, use_main_table=False):
        """Filters missing data in each `Alignment` object.

        This wraps the execution of `filter_missing_data` method for
        each `Alignment` object in the `alignments` attribute.

        Parameters
        ----------
        gap_threshold : int
            Integer between 0 and 100 defining the percentage above which
            a column with that gap percentage is removed.
        missing_threshold : int
            Integer between 0 and 100 defining the percentage above which
            a column with that gap+missing percentage is removed.
        table_in : string
            Name of database table containing the alignment data that is
            used for this operation.
        table_out : string
            Name of database table where the final alignment will be inserted
            (default is "filter").
        use_main_table : bool
            If True, both `table_in` and `table_out` are ignore and the main
            table `Alignment.db_idx` is used as the input and output
            table (default is False).
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.

        See Also
        --------
        Alignment.filter_missing_data
        """

        if use_main_table:
            table_in = table_out = self.master_table

        # Reset _partitions
        self.partitions = Partitions()

        self._filter_terminals(table_in, table_out, ns)

        self._filter_columns(gap_threshold, missing_threshold,
                             table_in, table_out, ns)

        self._reset_pipes(ns)

    def filter_segregating_sites(self, min_val, max_val, table_in=None,
                                 ns=None, pbar=None):
        """Filters `Alignment` objects according to segregating sites number.

        Filters `alignments` according to whether or not their number
        of segregating sites is within a specified range. `Alignment` objects
        with a number of segregating sites outside the range are move to
        the `shelve_alignments` and the counter of `filtered_alignments` is
        updated.

        Parameters
        ----------
        min_val : int
            Minimum number of segregating sites for the alignment to pass.
            Can be None, in which case there is no lower bound.
        max_val : int
            Maximum number of segregating sites for the alignment to pass.
            Can be None, in which case there is no upper bound.
        table_in : string
            Name of database table containing the alignment data that is
            used for this operation.
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.

        See Also
        --------
        Alignment.filter_segregating_sites
        """

        # Set progress pipes
        self._set_pipes(ns, pbar, total=len(self.alignments))

        # Stores the active Alignment.path after the fileter
        active_alns = []

        missing_symbols = [self.sequence_code[1], "-"]

        prev_idx = ""
        c = 1
        res = None
        collect = False
        for column, aln_idx in self.iter_columns(table_in):

            # This happens when the alignment changes during the iteration.
            if prev_idx != aln_idx:

                if prev_idx:
                    # Get the range result from the previous alignment
                    # Collect the result if collect has not been flagged
                    # (which means the alignment processing reached the
                    # end of the sequence) OR when the collect AND result
                    # are not None, which means that the early collect
                    # was set in the last nucleotide
                    if not collect or (collect and res):
                        res = self._test_range(s, min_val, max_val)
                        if res == "save":
                            active_alns.append(
                                self.alignment_idx[prev_idx].path)

                        res = False

                    collect = False

                # Reset variable sites counter
                s = 0

                # Update progress
                self._update_pipes(ns, pbar, value=c,
                                   msg="Filtering file {}".format(
                                       self.alignment_idx[aln_idx].name))

                prev_idx = aln_idx

                c += 1

            # The collect flag is set to True when the result of the
            # current alignment is known before finishing it.
            if collect:
                # Save the result for the current alignment
                if res == "save":
                    active_alns.append(self.alignment_idx[aln_idx].path)

                res = False
                continue

            v = len([x for x in set(column) if x not in missing_symbols])

            # Skip column without variation
            if v == 1:
                continue

            if v > 1:
                s += 1

            # Add these tests so that the method may exit earlier if the
            # conditions are met, precluding the analysis of the entire
            # alignment
            if min_val and v >= min_val and max_val is None:
                res = "save"
                collect = True

            if max_val and v > max_val:
                res = "shelve"
                collect = True

        if not collect or (collect and res):
            res = self._test_range(s, min_val, max_val)
            if res == "save":
                active_alns.append(self.alignment_idx[aln_idx].path)

        self.filtered_alignments["By variable sites"] = \
            len(self.alignments) - len(active_alns)

        # Update _partitions
        filtered_alns = [x.path for x in self.alignments.values()
                         if x.path not in active_alns]
        self.partitions.remove_partition(file_list=filtered_alns)
        # Update active files
        self.update_active_alignments(active_alns)

        self._reset_pipes(ns)

    @staticmethod
    def _test_range(s, min_val, max_val):
        """

        Parameters
        ----------
        s
        min_val
        max_val

        Returns
        -------

        """

        # If both values were specified, check if s is within range
        if max_val is not None and min_val is not None \
                and s >= min_val and s <= max_val:
            return "save"
        # If only min_val was specified, check if s is higher
        elif max_val is None and s >= min_val:
            return "save"
        # If only max_val was specified, check if s is lower
        elif min_val is None and s <= max_val:
            return "save"
        else:
            return "shelve"

    def filter_informative_sites(self, min_val, max_val, table_in=None,
                                 ns=None, pbar=None):
        """Filters `Alignment` objects according to informative sites number.

        Filters `alignments` according to whether or not their number
        of informative sites is within a specified range. `Alignment` objects
        with a number of informative sites outside the range are move to
        the `shelve_alignments` and the counter of `filtered_alignments` is
        updated.

        Parameters
        ----------
        min_val : int
            Minimum number of informative sites for the alignment to pass.
            Can be None, in which case there is no lower bound.
        max_val : int
            Maximum number of informative sites for the alignment to pass.
            Can be None, in which case there is no upper bound.
        table_in : string
            Name of database table containing the alignment data that is
            used for this operation.
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.

        See Also
        --------
        Alignment.filter_informative_sites
        """

        # Set progress pipes
        self._set_pipes(ns, pbar, total=len(self.alignments))

        # Stores the active Alignment.path after the fileter
        active_alns = []

        missing_symbols = [self.sequence_code[1], "-"]
        
        prev_idx = ""
        c = 1
        res = None
        collect = False
        for column, aln_idx in self.iter_columns(table_in):

            # This happens when the alignment changes during the iteration.
            if prev_idx != aln_idx:

                if prev_idx:
                    # Get the range result from the previous alignment
                    if not collect or (collect and res):
                        res = self._test_range(s, min_val, max_val)
                        if res == "save":
                            active_alns.append(
                                self.alignment_idx[prev_idx].path)

                        res = False
                    
                    collect = False

                # Reset informative sites counter
                s = 0
                
                self._update_pipes(ns, pbar, value=c,
                                   msg="Filtering file {}".format(
                                       self.alignment_idx[aln_idx].name))
                c += 1

                prev_idx = aln_idx

            # The collect flag is set to True when the result of the
            # current alignment is known before finishing it.
            if collect:
                # Save the result for the current alignment
                if res == "save":
                    active_alns.append(self.alignment_idx[aln_idx].path)

                res = False
                continue

            column = Counter(column)

            # Skip columns without variation
            if len(column) == 1:
                continue

            # Remove missing data
            for sym in missing_symbols:
                del column[sym]

            if len([x for x in column.values() if x >= 2]) >= 2:
                s += 1

            # Add these tests so that the method may exit earlier if the
            # conditions are met, precluding the analysis of the entire
            # alignment
            if min_val and s >= min_val and max_val is None:
                res = "save"
                collect = True

            if max_val and s > max_val:
                res = "shelve"
                collect = True

        if not collect or (collect and res):
            res = self._test_range(s, min_val, max_val)
            if res == "save":
                active_alns.append(self.alignment_idx[aln_idx].path)

        self.filtered_alignments["By informative sites"] = \
            len(self.alignments) - len(active_alns)

        # Update _partitions
        filtered_alns = [x.path for x in self.alignments.values()
                         if x.path not in active_alns]
        self.partitions.remove_partition(file_list=filtered_alns)
        # Update active files
        self.update_active_alignments(active_alns)

        self._reset_pipes(ns)

    def remove_taxa(self, taxa_list, mode="remove"):
        """Removes the specified taxa.

        Removes the specified taxa from the global `AlignmentList` attributes
        and from each `Alignment` object.

        This method supports a list or path to CSV file as the
        `taxa_list_file` argument. In case of CSV path, it parses the CSV
        file into a list. The CSV file should be a simple text file containing
        taxon names in each line.

        This method also supports two removal modes.
            - remove: removes the specified taxa
            - inverse: removes all but the specified taxa

        Parameters
        ----------
        taxa_list_file : list or string
            List containing taxa names or string with path to a CSV file
            containig the taxa names.
        mode : {"remove", "inverse"}
            The removal mode (default is remove).

        See Also
        --------
        Alignment.remove_taxa
        """

        # Checking if taxa_list is an input csv file:
        try:
            file_handle = open(taxa_list[0])
            taxa_list = self.read_basic_csv(file_handle)

        # If not, then the method's argument is already the final list
        except (IOError, IndexError):
            pass

        for alignment_obj in list(self.alignments.values()):
            alignment_obj.remove_taxa(taxa_list, mode=mode)

        # Updates taxa names
        if mode == "remove":
            for tx in taxa_list:
                try:
                    self.taxa_names.remove(tx)
                except ValueError:
                    # TODO: log a warning
                    pass
        elif mode == "inverse":
            self.taxa_names = [tx for tx in taxa_list if tx in self.taxa_names]

    def change_taxon_name(self, old_name, new_name):
        """Changes the name of a taxon. """

        for alignment_obj in list(self.alignments.values()):
            alignment_obj.change_taxon_name(old_name, new_name)

        # update taxa names
        self.taxa_names = [new_name if x == old_name else x
                           for x in self.taxa_names]

    def remove_file(self, filename_list):
        """Removes alignments.

        Removes `Alignment` objects based on their `path` attribute.

        Parameters
        ----------
        filename_list : list
            List of `Alignment.path` to be removed.
        """

        # If filename_list corresponds to all files in the current alignment
        # list, dispatch the clear_alignments methods
        if set(self.path_list) - set(filename_list) == set([]):
            self.clear_alignments()

        self.all_alignments = OrderedDict(
            (p, aln) for p, aln in self.all_alignments
            if p not in filename_list)
        self.alignments = OrderedDict(
            (p, aln) for p, aln in self.alignments
            if p not in filename_list)
        self.alignment_idx = OrderedDict(
            (idx, aln) for idx, aln in self.alignment_idx
            if aln.path not in filename_list
        )

        # Updates taxa names
        self.taxa_names = self._get_taxa_list()

        # Update size
        self.size = sum((x.locus_length for x in self.alignments.values()))

    def select_by_taxa(self, taxa_list, mode="strict"):
        """Selects a list of `Alignment` objects by taxa.

        This method is used to select `Alignments` based on whether they
        contain or exclude certain taxa. The selected alignments are
        returned in a list.

        Parameters
        ----------
        taxa_list : list
            List of taxa names.
        mode : str
            Selection mode. Can be:

                - "strict": Selects alignments containing all and only the
                provide taxa.
                - "inclusive": Selects alignments containing the provided
                taxa.
                - "relaxed": Selects alignments if they contain at least
                one of the provided taxa.

        Returns
        -------
        selected_alignments : list
            List of `Alignment` objects.
        """

        selected_alignments = []

        # taxa_list may be a file name (string) or a list containing the name
        # of the taxa. If taxa_list is a file name this code will parse the
        # csv file and return a list of the taxa. Otherwise, the taxa_list
        # variable remains the same.
        try:
            file_handle = open("".join(taxa_list))
            taxa_list = self.read_basic_csv(file_handle)
        except EnvironmentError:
            pass

        for alignment_obj in self.alignments.values():

            alignment_taxa = alignment_obj.taxa_idx.keys()

            # Selected only the alignments with the exact same taxa
            if mode == "strict":
                if set(taxa_list) == set(alignment_taxa):
                    selected_alignments.append(alignment_obj)

            # Selected alignments that include the specified taxa
            if mode == "inclusive":
                if set(taxa_list) - set(alignment_taxa) == set():
                    selected_alignments.append(alignment_obj)

            if mode == "relaxed":
                for taxon in taxa_list:
                    if taxon in alignment_taxa:
                        selected_alignments.append(alignment_obj)
                        continue

        return selected_alignments

    def code_gaps(self, table_out="gaps", table_in=None, use_main_table=False,
                  pbar=None, ns=None):
        """Code gaps in each `Alignment` object.

        This wraps the execution of `code_gaps` method for
        each `Alignment` object in the `alignments` attribute.

        Parameters
        ----------
        table_in : string
            Name of database table containing the alignment data that is
            used for this operation.
        table_out : string
            Name of database table where the final alignment will be inserted
            (default is "gaps").
        use_main_table : bool
            If True, both `table_in` and `table_out` are ignore and the main
            table `Alignment.db_idx` is used as the input and output
            table (default is False).
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.

        See Also
        --------
        Alignment.code_gaps
        """

        def gap_listing(sequence):
            """ Function that parses a sequence string and returns the
            position of indel events. The returned list is composed of
            tuples with the span of each indel """

            gap = "-+"
            gap_list = []

            for gobj in re.finditer(gap, sequence):
                
                gap_list.append((gobj.start(), gobj.end()))
            
            return gap_list

        def gap_binary_generator(sequence, gap_list):
            """ This function contains the algorithm to construct the binary
             state block for the indel events """

            binary_states = []

            for span in gap_list:

                substr = sequence[span[0]:span[1]]

                # Check if the entire span is only gaps
                if set(substr) == {"-"}:

                    if span[0] < 0:
                        prev = "n"
                    else:
                        prev = sequence[span[0] - 1]

                    try:
                        aft = sequence[span[1]]
                    except IndexError:
                        aft = "n"

                    # Check if the boundary positions of the span are gaps
                    # or not. The score 1 will only be given to spans that
                    # contain only gaps and are surrounded by non gaps.
                    # If there are gaps surrounding, then the score is "-"
                    if prev != "-" and aft != "-":
                        binary_states.append("1")
                    else:
                        binary_states.append("-")

                else:
                    binary_states.append("0")

            sequence = "".join([sequence, "".join(binary_states)])

            return sequence

        if use_main_table:
            table_in = table_out = self.master_table

        # Set progress pipes
        self._set_pipes(ns, pbar, total=len(self.alignments))

        # Storage of gaps positions for each alignment
        master_gaps = {}

        # Get the complete gap positions for each alignment and store it in
        # the master_gaps dict
        prev_idx = ""
        for taxon, seq, aln_idx in self.iter_alignments(table_in):

            if prev_idx != aln_idx:

                if prev_idx:
                    master_gaps[prev_idx] = complete_gap_list

                complete_gap_list = []

                prev_idx = aln_idx

            current_list = gap_listing(seq)
            complete_gap_list += [gap for gap in current_list if gap not in
                                  complete_gap_list]

        if aln_idx:
            master_gaps[aln_idx] = complete_gap_list

        # Create temporary table
        temp_table = ".codegaps"
        self._create_table(temp_table)
        c = 1

        # Create temporary cursor for database
        temp_cur = self.con.cursor()

        prev_idx = ""
        for txId, taxon, seq, aln_idx in self.iter_alignments(
                table_in, include_txid=True):

            if prev_idx != aln_idx:

                if prev_idx:
                    if len(master_gaps[prev_idx]):
                        aln_obj.restriction_range = "{}-{}".format(
                            int(aln_obj.locus_length),
                            len(master_gaps[prev_idx]) +
                            int(aln_obj.locus_length))

                aln_obj = self.alignment_idx[aln_idx]

                self._update_pipes(ns, pbar, value=c)
                c += 1

                prev_idx = aln_idx

            final_seq = gap_binary_generator(seq, master_gaps[aln_idx])

            temp_cur.execute(
                "INSERT INTO [.codegaps] VALUES (?, ?, ?, ?)",
                (txId, taxon, final_seq, aln_idx))

        if prev_idx:
            if len(master_gaps[aln_idx]):
                aln_obj.restriction_range = "{}-{}".format(
                    int(aln_obj.locus_length),
                    len(master_gaps[aln_idx]) +
                    int(aln_obj.locus_length))

        if self._table_exists(table_out):
            self.cur.execute("DROP TABLE [{}];".format(table_out))

        # Replace table_out with collapsed table
        self.cur.execute("ALTER TABLE [.codegaps] RENAME TO [{}]".format(
            table_out))

        self._reset_pipes(ns)

        #TODO: must check binary gap generator

    @staticmethod
    def write_loci_correspondence(hap_dict, output_file, dest="./"):
        """Writes the file mapping taxa to unique haplotypes for `collapse`.

        Parameters
        ----------
        hap_dict : dict
            Dictionary mapping the haplotype (key) to a list of taxa (value)
            sharing the same haplotype.
        output_file : str
            Name of the haplotype correspondence file
        dest : str, optional
            Path to directory where the `output_file` will be written.
        """

        output_handle = open(join(dest, output_file + ".haplotypes"), "w")

        for haplotype, taxa_list in sorted(hap_dict.items()):
            output_handle.write("%s: %s\n" % (haplotype, "; ".join(taxa_list)))

        output_handle.close()

    def collapse(self,  write_haplotypes=True, haplotypes_file=None,
                 dest=".", conversion_suffix="", haplotype_name="Hap",
                 table_in=None, table_out="collapsed", use_main_table=False,
                 ns=None, pbar=None):
        """Collapses equal sequences for each `Alignment` object.

        This wraps the execution of `collapse` method for
        each `Alignment` object in the `alignments` attribute.

        Parameters
        ----------
        write_haplotypes : bool
            If True, a file mapping the taxa names name to their respective
            haplotype name will be generated (default is True).
        haplotypes_file : string
            Name of the file mapping taxa names to their respective haplotype.
            Only used when `write_haplotypes` is True. If it not provided and
            `write_haplotypes` is True, the file name will be determined
            from `Alignment.name`.
        dest : string
            Path of directory where `haplotypes_file` will be generated
            (default is ".").
        conversion_suffix : string
            Suffix appended to the haplotypes file. Only used when
            `write_haplotypes` is True and `haplotypes_file` is None.
        haplotype_name : string
            Prefix of the haplotype string. The final haplotype name will be
            `haplotype_name` + <int> (default is "Hap").
        table_in : string
            Name of database table containing the alignment data that is
            used for this operation.
        table_out : string
            Name of database table where the final alignment will be inserted
            (default is "collapsed").
        use_main_table : bool
            If True, both `table_in` and `table_out` are ignore and the main
            table `Alignment.db_idx` is used as the input and output
            table (default is False).
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.

        See Also
        --------
        Alignment.collapse
        """

        if use_main_table:
            table_out = table_in = self.master_table

        # Set progress pipes
        self._set_pipes(ns, pbar, total=len(self.alignments))
        c = 1

        # Create temporary table
        temp_table = ".collapsed"
        self._create_table(temp_table)

        # Set temporary cursor to perform database changes while querying
        temp_cur = self.con.cursor()

        prev_idx = ""
        for taxon, seq, aln_idx in self.iter_alignments(table_in):

            # This happens when the alignment changes during the iteration.
            if aln_idx != prev_idx:

                self._update_pipes(ns, pbar, value=c)
                c += 1

                if prev_idx:
                    # Write the haplotypes file for the previous alignment
                    if write_haplotypes:
                        if not haplotypes_file:
                            hf = aln_obj.sname + conversion_suffix
                        else:
                            hf = haplotypes_file

                        self.write_loci_correspondence(hap_dic, hf, dest)

                aln_obj = self.alignment_idx[aln_idx]

                # Dict that will store the hashes of each sequence. It will
                # be used to check if the current sequence has already been
                # processed. The value will be the haplotype
                hash_list = {}

                # hap_dic will store the haplotype name as key and a list of
                # the taxa with the same sequence as a list value
                hap_dic = {}
                hap_counter = 0

                prev_idx = aln_idx

            cur_hash = hash(seq)
            
            if cur_hash not in hash_list:

                # Create name for new haplotype and update other haplotype
                # related variables
                haplotype = "{}_{}".format(haplotype_name, hap_counter + 1)
                hash_list[cur_hash] = haplotype
                hap_dic[haplotype] = [taxon]

                # Add data to table
                temp_cur.execute(
                    "INSERT INTO [.collapsed] VALUES (?, ?, ?, ?)",
                    (hap_counter, unicode(haplotype), seq, aln_idx))

                hap_counter += 1

            else:

                hap_dic[hash_list[cur_hash]].append(taxon)

        # Write haplotypes for the last alignment
        if aln_obj:
            if write_haplotypes:
                if not haplotypes_file:
                    hf = aln_obj.sname + conversion_suffix
                else:
                    hf = haplotypes_file
                self.write_loci_correspondence(hap_dic, hf, dest)

        # The collapse operation is special in the sense that the former taxon
        # names are no longer valid. Therefore, we drop the previous table and
        # populate a new one with the collapsed data
        if self._table_exists(table_out):
            self.cur.execute("DROP TABLE [{}];".format(table_out))

        # Replace table_out with collapsed table
        self.cur.execute("ALTER TABLE [.collapsed] RENAME TO [{}]".format(
            table_out))

        self._reset_pipes(ns)

    def consensus(self, consensus_type, single_file=False, table_in=None,
                  table_out=None, use_main_table=False, ns=None,
                  pbar=None):
        """Creates consensus sequences for each alignment.

        This wraps the execution of `consensus` method for
        each `Alignment` object in the `alignments` attribute. It has an
        additional argument, `single_file`, which can be used to merge
        the consensus sequence from each alignment into a new `Alignment`
        object.

        Parameters
        ----------
        consensus_type : {"IUPAC", "Soft mask", "Remove", "First sequence"}
            Type of variation handling. See summary above.
        single_file : bool
            If True, the consensus sequence of each `Alignment` object will
            be merged into a single `Alignment` object.
        table_in : string
            Name of database table containing the alignment data that is
            used for this operation.
        table_out : string
            Name of database table where the final alignment will be inserted
            (default is "consensus").
        use_main_table : bool
            If True, both `table_in` and `table_out` are ignore and the main
            table `Alignment.db_idx` is used as the input and output
            table (default is False).
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.
            
        Returns
        consensus_aln : trifusion.process.sequence.Alignment
            Alignment object with all the consensus sequences merged. It's
            only returned when the `single_file` parameter is set to True.
        
        See Also
        --------
        Alignment.consensus
        """

        def add_to_database(txid, tx, s, fidx, aln_name=None):

            # Get first sequence and set the skip flag to True
            temp_cur.execute(
                "INSERT INTO [{}] VALUES (?, ?, ?, ?)".format(
                    temp_table), (txid, tx, s, fidx))

            seq_len = len(s)

            size[0] += seq_len

            if single_file:
                # Update _partitions
                self.partitions.add_partition("consensus", length=seq_len)
                taxa_idx[txid] = aln_name
            else:
                part = Partitions()
                part.add_partition("consensus", length=seq_len)
                self._idx += 1
                aln = Alignment(
                    aln_name,
                    sql_cursor=self.cur,
                    taxa_idx=taxa_idx,
                    ignore_db_check=True, partitions=part,
                    sequence_code=self.sequence_code,
                    locus_length=seq_len,
                    db_idx=self._idx,
                    temp_dir=os.path.dirname(self.sql_path))
                idx_storage[fidx] = aln
                aln_storage[aln_name] = aln

        if use_main_table:
            table_in = table_out = self.master_table
        else:
            table_out = table_out if table_out else self.master_table

        idx_storage = OrderedDict()
        aln_storage = OrderedDict()

        # Create new cursor to insert data into database while the main
        # cursor returns the query
        temp_cur = self.con.cursor()

        # Variables that will store the taxa_list and _taxa_idx to provide
        # to the Alignment object
        if single_file:
            taxa_idx = {}
        else:
            taxa_idx = {"consensus": 0}

        # Reset AlignmentList size
        size = [0]

        # Reset the _partitions object
        self.partitions = Partitions()

        # Create temporary table
        temp_table = ".consensus"
        self._create_table(temp_table, index=["conindex", "aln_idx"])

        # The First sequence mode is unique in the sense that it does not
        # require an iteration over the columns of each alignment.
        if consensus_type == "First sequence":
            skip = False
            prev_idx = ""
            # Set progress pipes
            self._set_pipes(ns, pbar, len(self.alignments))
            c = 1
            for txId, taxon, seq, aln_idx in self.iter_alignments(
                    table_in, include_txid=True):

                # This happens when the alignment changes during the iteration.
                if aln_idx != prev_idx:
                    
                    aln_name = self.alignment_idx[aln_idx].name

                    # Update progress
                    self._update_pipes(
                        ns, pbar, value=c,
                        msg="Processing file {}".format(aln_name))
                    c += 1

                    # Set final aln_idx. When single_file is set to True, all
                    # final aln_idx are 1. Else, the original aln_idx is used.
                    final_idx = 1 if single_file else aln_idx

                    skip = False

                if skip:
                    continue

                add_to_database(0, "consensus", seq, final_idx, aln_name)

                # Set skip flag to True so that all remaining sequences of the
                # current alignment are ignored
                skip = True

        else:
            self._set_pipes(ns, pbar, total=len(self.alignments))
            c = 1

            prev_idx = ""
            for column, aln_idx in self.iter_columns(table_in):

                if aln_idx != prev_idx:

                    if prev_idx:

                        final_seq = "".join(consensus_seq)

                        if single_file:
                            add_to_database(c - 1, aln_name, final_seq,
                                            final_idx, aln_name)
                        else:
                            add_to_database(0, "consensus", final_seq,
                                            final_idx, aln_name)

                    # Empty consensus sequence
                    consensus_seq = []

                    # Set final aln_idx. When single_file is set to True, all
                    # final aln_idx are 1. Else, the original aln_idx is used.
                    final_idx = 1 if single_file else aln_idx

                    aln_name = self.alignment_idx[aln_idx].name

                    self._update_pipes(
                        ns, pbar, value=c,
                        msg="Processing file {}".format(aln_name))
                    c += 1

                    prev_idx = aln_idx

                column = list(set(column))

                # If invariable, include in consensus and skip. Increases speed
                # when there are no missing or gap characters
                if len(column) == 1:
                    consensus_seq.append(column[0])
                    continue

                # Remove missing data and gaps
                column = sorted([x for x in column if
                                 x != self.sequence_code[1] and x != "-"])

                # In case there is only missing/gap characters
                if not column:
                    consensus_seq.append(self.sequence_code[1])
                    continue

                # Check for invariable without missing data
                if len(column) == 1:
                    consensus_seq.append(column[0])
                    continue

                if consensus_type == "IUPAC":
                    # Convert variation into IUPAC code
                    try:
                        consensus_seq.append(iupac["".join(column)])
                    except KeyError:
                        iupac_code = sorted(set(list(
                            itertools.chain.from_iterable(
                                [x if x in dna_chars else iupac_rev[x]
                                 for x in column]))))
                        consensus_seq.append(iupac["".join(iupac_code)])
                    continue

                elif consensus_type == "Soft mask":
                    consensus_seq.append(self.sequence_code[1])
                    continue

                elif consensus_type == "Remove":
                    continue

        if prev_idx:
            aln_name = self.alignment_idx[aln_idx].name
            final_seq = "".join(consensus_seq)
            if single_file:
                add_to_database(c - 1, aln_name, final_seq,
                                final_idx, aln_name)
            else:
                add_to_database(0, "consensus", final_seq,
                                final_idx, aln_name)

        if self._table_exists(table_out):
            self.cur.execute("DROP TABLE [{}];".format(table_out))

        # Replace table_out with collapsed table
        self.cur.execute(
            "ALTER TABLE [{}] RENAME TO [{}]".format(
                temp_table, table_out))

        if single_file:
            self.size = size[0]
            self._idx += 1
            aln = Alignment("consensus", sql_cursor=self.cur,
                            taxa_idx=taxa_idx,
                            ignore_db_check=True, partitions=self.partitions,
                            sequence_code=self.sequence_code,
                            locus_length=self.size,
                            db_idx=self._idx,
                            temp_dir=os.path.dirname(self.sql_path))
            self.alignment_idx = OrderedDict()
            self.alignment_idx[1] = aln
            self.taxa_names = taxa_idx.keys()
        else:
            self.taxa_names = ["consensus"]
            self.alignment_idx = idx_storage
            self.alignments = aln_storage
            self.shelved_idx = []

        self._reset_pipes(ns)

    def reverse_concatenate(self, aln_name=None, table_in=None,
                            table_out=None, pbar=None, ns=None):
        """Reverse a concatenated file according to the _partitions.

        This is basically a wrapper of the `reverse_concatenate` method
        of `Alignment` and is rarely useful. It is preferable to retrieve
        an `Alignment` object, provide some _partitions and then
        perform the reverse concatenation. Since the reverse concatenation
        operation can only be applied to `Alignment` objects, this method
        actually creates a concatenated `Alignment` .

        Parameters
        ----------
        aln_name : str
            :attr:`Alignment.path` attribute with path to alignment file.
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        reverted_alns : AlignmentList
            AlignmentList with the _partitions as `Alignment` objects.
        """
        
        def add_alignment(part_name, part_len, taxa_idx, part):

            # Create Alignment object
            self._idx += 1
            current_aln = Alignment(part_name,
                                    sql_cursor=self.cur,
                                    locus_length=part_len,
                                    sequence_code=self.sequence_code,
                                    taxa_idx=taxa_idx,
                                    partitions=part,
                                    ignore_db_check=True,
                                    db_idx=self._idx,
                                    temp_dir=os.path.dirname(self.sql_path))

            return current_aln

        table_out = table_out if table_out else self.master_table

        # If aln_name is provided, reverse concatenation will be applied
        # on a single alignment. Change active file set to contain only
        # that alignment.
        if aln_name:
            self.update_active_alignments([aln_name])

        # Set progress pipes
        self._set_pipes(ns, pbar, total=len(self.taxa_names), ignore_sa=True)

        taxa_idx_master = defaultdict(dict)

        rev_aln_idx = OrderedDict()
        alns = OrderedDict()

        # Create temporary table
        temp_table = ".reversedata"
        self._create_table(temp_table, index=("revindex", "aln_idx"))

        temp_cur = self.con.cursor()

        # Get corrected partition names for sqlite database
        part_map = self._get_part_names()

        for p, (taxon, seq, aln_idx) in enumerate(
                self.iter_alignments(table_in)):

            self._update_pipes(ns, pbar, value=p + 1, ignore_sa=True,
                               msg="Processing taxon {}".format(taxon))

            # Index for _partitions
            part_idx = 1

            for name, part_range in self.partitions:

                name = part_map[name]

                if part_range[1]:

                    for i in range(3):

                        part_seq = seq[part_range[0][0]:
                                       part_range[0][1] + 1][i::3]

                        if part_seq.replace(self.sequence_code[1], "") == "":
                            part_idx += 1
                            continue

                        # Adapt name of partition for codon variant
                        cname = name + str(i)

                        taxa_idx_master[cname][taxon] = p

                        temp_cur.execute(
                            "INSERT INTO [{}] VALUES (?, ?, ?, ?)".format(
                                temp_table), (p, taxon, part_seq, part_idx))

                        part_idx += 1

                else:

                    part_seq = seq[part_range[0][0]:part_range[0][1] + 1]

                    if part_seq.replace(self.sequence_code[1], "") == "":
                        part_idx += 1
                        continue

                    taxa_idx_master[name][taxon] = p

                    temp_cur.execute(
                        "INSERT INTO [{}] VALUES (?, ?, ?, ?)".format(
                            temp_table), (p, taxon, part_seq, part_idx))

                    part_idx += 1

        if self._table_exists(table_out):
            self.cur.execute("DROP TABLE [{}];".format(table_out))

        self._create_table(table_out, index=("finalrevindex", "aln_idx"))
        self.cur.execute(
            "INSERT INTO [{}] (txId, taxon, seq, aln_idx) "
            "SELECT txId, taxon, seq, aln_idx "
            "FROM [{}] "
            "ORDER BY aln_idx".format(table_out, temp_table))

        part_idx = 1
        for name, part_range in self.partitions:

            name = part_map[name]
            
            if part_range[1]:

                for i in range(3):

                    cname = "{}_codon_{}".format(name, str(i))

                    part_len = (part_range[0][1] - part_range[0][0] + 1) / 3

                    taxa_idx = taxa_idx_master[cname]

                    p = Partitions()
                    p.add_partition(cname, length=part_len)

                    saln = add_alignment(cname, part_len,
                                         taxa_idx=taxa_idx, part=p)

                    alns[cname] = saln
                    rev_aln_idx[part_idx] = saln
                    
                    part_idx += 1

            else:

                part_len = part_range[0][1] - part_range[0][0] + 1

                taxa_idx = taxa_idx_master[name]

                p = Partitions()
                p.add_partition(name, length=part_len)
                saln = add_alignment(name, part_len,
                                     taxa_idx=taxa_idx, part=p)

                alns[name] = saln
                rev_aln_idx[part_idx] = saln
                
                part_idx += 1

        self.alignments = alns
        self.all_alignments = alns
        self.alignment_idx = rev_aln_idx

    def _get_part_names(self):
        """Returns partition name strings compliant with sqlite database

        Returns
        -------

        """

        # Get corrected partition names for sqlite database
        part_map = {}
        for name, _ in self.partitions:
            part_map[name] = "".join([x for x in name if x.isalnum()])

        return part_map

    def _get_interleave_data(self, table_name=None, ns_pipe=None,
                             pbar=None):

        if self.cur.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND"
                " name='.interleavedata'").fetchall():

            self.cur.execute("DELETE FROM [.interleavedata]")

        else:
            self.cur.execute("CREATE TABLE [.interleavedata]("
                             "taxon TEXT,"
                             "seq TEXT,"
                             "slice INT,"
                             "aln_idx INT)")
            self.cur.execute("CREATE INDEX interindex ON "
                            "[.interleavedata](slice, aln_idx)")

        prev_file = ""

        temp_cur = self.con.cursor()

        for taxon, seq, aln_idx in self.iter_alignments(
                table_name=table_name):

            if aln_idx != prev_file:
                prev_file = aln_idx
                aln_obj = self.alignment_idx[aln_idx]
                
            counter = 0

            for i in xrange(90, aln_obj.locus_length, 90):

                temp_cur.execute("INSERT INTO [.interleavedata] VALUES "
                                 " (?, ?, ?, ?)", (taxon, seq[counter:i], i,
                                                   aln_idx))
                counter = i

            try:
                if aln_obj.locus_length % 90:
                    i += 1
                    temp_cur.execute("INSERT INTO [.interleavedata] VALUES "
                                     " (?, ?, ?, ?)",
                                     (taxon, seq[counter:], i, aln_idx))
            except UnboundLocalError:
                temp_cur.execute("INSERT INTO [.interleavedata] VALUES "
                                 " (?, ?, ?, ?)", (taxon, seq, 0, aln_idx))

        return True

    def _get_partition_data(self, table_name, ns=None, pbar=None):
        """

        Parameters
        ----------
        table_name
        ns
        pbar

        Returns
        -------

        """

        partition_table = ".partitiondata"
        if self.cur.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND"
                " name='{}'".format(partition_table)).fetchall():

            self.cur.execute("DELETE FROM [{}]".format(partition_table))

        else:
            self.cur.execute("CREATE TABLE [{}]("
                             "taxon TEXT,"
                             "seq TEXT,"
                             "part_name TEXT,"
                             "part INT,"
                             "aln_idx INT)".format(partition_table))
            self.cur.execute("CREATE INDEX partindex ON "
                            "[{}](part, aln_idx)".format(partition_table))

        temp_cur = self.con.cursor()

        # Get corrected partition names for sqlite database
        part_map = self._get_part_names()

        for taxon, seq, aln_idx in self.iter_alignments(table_name):

            part_idx = 0

            for name, part_range in self.partitions:

                name = part_map[name]

                if name not in self.partition_count:
                    self.partition_count[name] = [[], 0]

                part_seq = seq[part_range[0][0]:part_range[0][1] + 1]

                if part_seq.replace(self.sequence_code[1], "") == "":
                    part_idx += 1
                    continue

                temp_cur.execute(
                    "INSERT INTO [{}] "
                    "VALUES (?, ?, ?, ?, ?)".format(partition_table),
                    (taxon, part_seq, name, part_idx, aln_idx))
                part_idx += 1
                self.partition_count[name][0].append(taxon)
                self.partition_count[name][1] = len(part_seq)

        return True

    def _setup_newfile(self, fh, aln_file, output_dir, suffix, output_file,
                       ns):

        # Close previous file object, if exists
        if fh:
            fh.close()

        # Get path/file name of alignment/partition
        if suffix and not output_file:
            output_file = self.alignment_idx[aln_file].sname + suffix
            if output_dir:
                output_file = join(output_dir, output_file)
                if not exists(output_dir):
                    os.makedirs(output_dir)
        else:
            output_file = output_file

        if exists(output_file):

            # File exists, issue a warning through the appropriate pipe
            if ns:
                if not ns.apply_all:
                    ns.file_dialog = output_file
                    while not ns.status:
                        pass

            # When the dialog has been close, check if the file is to be
            # skipped or overwritten
            if ns:
                if ns.status == "skip":
                    return None, None

        # Reset pipes, if any
        if ns:
            ns.status = None

        # Return new file object
        fh = open(output_file, "w")

        return fh, output_file

    def _write_fasta(self, suffix, output_file, **kwargs):

        ld_hat = kwargs.get("ld_hat", False)
        interleave = kwargs.get("interleave", False)
        output_dir = kwargs.pop("output_dir", None)
        table_name = kwargs.get("table_name", self.master_table)
        ns = kwargs.get("ns_pipe", None)
        pbar = kwargs.get("pbar", None)

        # File object that will be used to write sequence data
        fh = None

        # Stores the string of the last file
        prev_file = ""

        self._set_pipes(ns, pbar, total=len(self.alignments))
        c = 1
        
        for taxon, seq, aln_idx in self.iter_alignments(table_name):

            if aln_idx != prev_file:
                prev_file = aln_idx
                fh, _ = self._setup_newfile(
                    fh, aln_idx, output_dir, suffix, output_file, ns)

                self._update_pipes(ns, pbar, value=c,
                                   msg="Writing files")
                c += 1

                # If fh is set to None, the current file is set to skip
                if not fh:
                    continue

                # Do stuff at the beginning of the file, before inserting
                # Sequence data

                # Get Alignment object containing info of the current alignment
                aln_obj = self.alignment_idx[aln_idx]

                # If LD HAT sub format has been specificed, write the first
                # line containing the number of sequences, sites and
                # genotype phase
                if ld_hat:
                    fh.write("{} {} {}\n".format(
                        len(aln_obj.taxa_idx), aln_obj.locus_length, "2"))

            # Skip row if fh is not set (skipping files)
            if not fh:
                continue

            if ld_hat:
                # Truncate sequence name to 30 characters
                fh.write(">%s\n" % (taxon[:30]))
                # Limit each sequence line to 2000 characters
                if len(seq) > 2000:
                    for i in range(0, len(seq), 2000):
                        fh.write("%s\n" % (seq[i:i + 2000]))
            elif interleave:
                fh.write(">{}\n".format(taxon))
                counter = 0
                for i in range(90, aln_obj.locus_length, 90):
                    fh.write("{}\n".format(seq[counter:i]))
                    counter = i

                fh.write("{}\n".format(seq[counter:]))
            else:
                fh.write(">{}\n{}\n".format(taxon, seq))

        # When all files are skipeed, fh remains None
        if fh:
            fh.close()

    def _write_phylip_partitions(self, aln_obj, partition_file,
                                 output_file, model_phylip):

        # In case there is a concatenated alignment being written
        if not aln_obj.partitions.is_single() and partition_file:
            partition_file = open(output_file.split(".")[0]
                                  + "_part.File", "w")
            for name, lrange in self.partitions:
                # Get model from app if it exists and there are no codon
                # positions
                model = model_phylip if self.partitions.models[name] == \
                                        [None] or len(
                    self.partitions.models[name][1]) > 1 \
                    else self.partitions.models[name][1][0]
                partition_file.write("%s, %s = %s\n" % (
                    model if model else
                    self.sequence_code[0], name,
                    "-".join([str(x + 1) for x in
                              lrange[0]])))
            partition_file.close()

    def _write_phylip(self, suffix, output_file, **kwargs):

        # Get relevant keyword arguments
        interleave = kwargs.get("interleave", False)
        tx_space_phy = kwargs.get("tx_space_phy", 40)
        cut_space_phy = kwargs.get("cut_space_phy", 39)
        phy_truncate_names = kwargs.get("phy_truncate_names", False)
        partition_file = kwargs.get("partition_file", None)
        model_phylip = kwargs.get("model_phylip", None)
        table_name = kwargs.get("table_name", None)
        ns = kwargs.get("ns_pipe", None)
        pbar = kwargs.get("pbar", None)
        output_dir = kwargs.get("output_dir", None)

        # File object that will be used to write sequence data
        fh = None

        # Stores the string of the last file
        prev_file = ""

        # Change taxa space if phy_truncate_names option is set to True
        if phy_truncate_names:
            cut_space_phy = 10

        if interleave:

            if not self.interleave_data:
                self.interleave_data = self._get_interleave_data(
                        table_name, ns, pbar)

            for taxon, seq, p, aln_idx in self.cur.execute(
                "SELECT taxon, seq, slice, aln_idx from [.interleavedata] "
                "ORDER BY aln_idx, slice"):

                if prev_file != aln_idx:
                    fh, of = self._setup_newfile(
                        fh, aln_idx, output_dir, suffix, output_file, ns)
                    prev_file = aln_idx

                    if not fh:
                        continue

                    aln_obj = self.alignment_idx[aln_idx]
                    
                    self._write_phylip_partitions(aln_obj, partition_file,
                                                  of, model_phylip)

                    fh.write("{} {}\n".format(
                        len(aln_obj.taxa_idx) - len(aln_obj.shelved_taxa),
                        aln_obj.locus_length))

                    write_tx = True
                    prev = 90

                if not fh:
                    continue

                if p != prev:
                    fh.write("\n")
                    prev = p
                    write_tx = False

                if write_tx:
                    fh.write("{} {}\n".format(
                        taxon[:cut_space_phy].ljust(tx_space_phy),
                        seq))
                else:
                    fh.write("{}\n".format(seq))

        else:

            for taxon, seq, aln_idx in self.iter_alignments(table_name):

                if aln_idx != prev_file:
                    prev_file = aln_idx
                    fh, of = self._setup_newfile(
                        fh, aln_idx, output_dir, suffix, output_file, ns)

                    if not fh:
                        continue

                    aln_obj = self.alignment_idx[aln_idx]

                    self._write_phylip_partitions(aln_obj, partition_file,
                                                  of, model_phylip)

                    # Get Alignment object containing info of the current
                    # alignment
                    aln_obj = self.alignment_idx[aln_idx]

                    fh.write("{} {}\n".format(
                        len(aln_obj.taxa_idx) - len(aln_obj.shelved_taxa),
                        aln_obj.locus_length))

                if not fh:
                    continue

                fh.write("{} {}\n".format(
                    taxon[:cut_space_phy].ljust(tx_space_phy),
                    seq))

        # When all files are skipeed, fh remains None
        if fh:
            fh.close()

    @staticmethod
    def _write_nexus_partitions(aln_obj, use_charset, fh,
                                use_nexus_models, outgroup_list):

        if use_charset:
            # Writing _partitions, if any
            if not aln_obj.partitions.is_single():
                fh.write("\nbegin mrbayes;\n")
                p = 0
                # Full gene _partitions
                for name, lrange in aln_obj.partitions:
                    # If there are codon _partitions, write those
                    if lrange[1]:
                        for i in lrange[1]:
                            fh.write("\tcharset %s_%s = %s-%s\\3;\n" %
                                           (name, i + 1, i + 1,
                                            lrange[0][1] + 1))
                            p += 1
                    else:
                        fh.write("\tcharset %s = %s-%s;\n" %
                                       (name, lrange[0][0] + 1,
                                        lrange[0][1] + 1))
                        p += 1

                fh.write("\tpartition part = %s: %s;\n\tset "
                               "partition=part;\nend;\n" %
                               (p, ", ".join([name for name in
                                aln_obj.partitions.get_partition_names()])))

        # Write models, if any
        if use_nexus_models:
            if any(aln_obj.partitions.models.values()[0][0]):

                fh.write("\nbegin mrbayes;\n")
                i = 1
                for name, models in aln_obj.partitions.models.items():
                    if models[1]:
                        for m in models[1]:
                            if m:
                                m_str = \
                                    aln_obj.partitions._models["mrbayes"][m]
                                if not aln_obj.partitions.is_single():
                                    fh.write("applyto=({}) "
                                             "lset ".format(i) +
                                             " ".join(m_str) + "\n")
                                else:
                                    fh.write("lset " + " ".join(m_str) + "\n")
                            i += 1
                fh.write("end;\n")

        # In case outgroup taxa are specified
        if outgroup_list is not None:

            # This assures that only the outgroups present in the current
            #  file are written
            compliant_outgroups = [taxon for taxon in outgroup_list
                                   if taxon in aln_obj.taxa_idx]
            if compliant_outgroups is not []:
                fh.write("\nbegin mrbayes;\n\toutgroup %s\nend;\n" %
                               (" ".join(compliant_outgroups)))

    @staticmethod
    def _write_nexus_header(aln_obj, fh, gap, interleave):

        if aln_obj.restriction_range is not None:
            fh.write("#NEXUS\n\nBegin data;\n\tdimensions "
                     "ntax={} nchar={} ;\n\tformat datatype=mixed"
                     "({}:1-{}, restriction:{}) interleave={} "
                     "gap={} missing={} ;\n\tmatrix\n".format(
                        len(aln_obj.taxa_idx) - len(aln_obj.shelved_taxa),
                        aln_obj.locus_length,
                        aln_obj.sequence_code[0],
                        aln_obj.locus_length - 1,
                        aln_obj.restriction_range,
                        "yes" if interleave else "no",
                        gap,
                        aln_obj.sequence_code[1].upper()))
        else:
            fh.write(
                "#NEXUS\n\nBegin data;\n\tdimensions ntax={}"
                " nchar={} ;\n\tformat datatype={} "
                "interleave={} gap={} missing={} ;\n\t"
                "matrix\n".format(
                    len(aln_obj.taxa_idx) - len(aln_obj.shelved_taxa),
                    aln_obj.locus_length,
                    aln_obj.sequence_code[0],
                    "yes" if interleave else "no",
                    gap,
                    aln_obj.sequence_code[1].upper()))

    def _write_nexus(self, suffix, output_file, **kwargs):

        # Get relevant keyword arguments
        interleave = kwargs.get("interleave", False)
        tx_space_nex = kwargs.get("tx_space_nex", 40)
        cut_space_nex = kwargs.get("cut_space_nex", 39)
        gap = kwargs.get("gap", "-")
        use_charset = kwargs.get("use_charset", True)
        use_nexus_models = kwargs.get("use_nexus_models", True)
        outgroup_list = kwargs.get("outgroup_list", None)
        table_name = kwargs.get("table_name", None)
        output_dir = kwargs.get("output_dir", None)
        ns = kwargs.get("ns_pipe", None)
        pbar = kwargs.get("pbar", None)

        # File object that will be used to write sequence data
        fh = None

        # Stores the string of the last file
        prev_file = ""

        if interleave:

            if not self.interleave_data:
                self.interleave_data = self._get_interleave_data(
                        table_name, ns, pbar)

            for taxon, seq, p, aln_idx in self.cur.execute(
                    "SELECT taxon, seq, slice, aln_idx "
                    "FROM [.interleavedata] "
                    "ORDER BY aln_idx, slice"):

                if aln_idx != prev_file:

                    if fh:
                        fh.write(";\n\tend;")

                        self._write_nexus_partitions(aln_obj, use_charset, fh,
                                                     use_nexus_models,
                                                     outgroup_list)

                    prev_file = aln_idx
                    fh, of = self._setup_newfile(
                        fh, aln_idx, output_dir, suffix, output_file, ns)

                    if not fh:
                        continue

                    aln_obj = self.alignment_idx[aln_idx]

                    self._write_nexus_header(aln_obj, fh, gap, interleave)

                    prev = 90

                if not fh:
                    continue

                if p != prev:
                    fh.write("\n")
                    prev = p

                fh.write("{} {}\n".format(
                    taxon[:cut_space_nex].ljust(tx_space_nex),
                    seq.upper()))

        else:

            for taxon, seq, aln_idx in self.iter_alignments(table_name):

                if aln_idx != prev_file:

                    if fh:
                        fh.write(";\n\tend;")

                        self._write_nexus_partitions(aln_obj, use_charset, fh,
                                                     use_nexus_models,
                                                     outgroup_list)

                    prev_file = aln_idx
                    fh, of = self._setup_newfile(
                        fh, aln_idx, output_dir, suffix, output_file, ns)

                    if not fh:
                        continue

                    aln_obj = self.alignment_idx[aln_idx]

                    self._write_nexus_header(aln_obj, fh, gap, interleave)

                if not fh:
                    continue

                fh.write("{} {}\n".format(taxon[:cut_space_nex].ljust(
                    tx_space_nex), seq))

            if fh:
                fh.write(";\n\tend;")
                self._write_nexus_partitions(aln_obj, use_charset, fh,
                                             use_nexus_models,
                                             outgroup_list)

        # When all files are skipeed, fh remains None
        if fh:
            fh.close()

    def _write_stockholm(self, suffix, output_file, **kwargs):

        table_name = kwargs.get("table_name", None)
        ns = kwargs.get("ns_pipe", None)
        output_dir = kwargs.get("output_dir", None)
        pbar = kwargs.get("pbar", None)

        # File object that will be used to write sequence data
        fh = None

        # Stores the string of the last file
        prev_file = ""

        for taxon, seq, aln_idx in self.iter_alignments(table_name):

            if aln_idx != prev_file:

                if fh:
                    fh.write("//\n")

                prev_file = aln_idx
                fh, _ = self._setup_newfile(
                    fh, aln_idx, output_dir, suffix, output_file, ns)

                fh.write("# STOCKHOLM V1.0\n")

            fh.write("{}\t{}\n".format(taxon, seq))

        # When all files are skipeed, fh remains None
        if fh:
            fh.close()

    def _write_gphocs(self, suffix, output_file, **kwargs):

        table_name = kwargs.get("table_name", None)
        ns = kwargs.get("ns_pipe", None)
        pbar = kwargs.get("pbar", None)
        output_dir = kwargs.get("output_dir", None)

        # File object that will be used to write sequence data
        fh = None

        # Check if partition data has already been built in the database.
        # If not, create it
        if not self.partition_data:
            self.partition_data = self._get_partition_data(table_name, ns,
                                                           pbar)

        # Stores the string of the last file
        prev_idx = ""
        prev_part = ""
        for taxon, seq, pname, pidx, aln_idx in self.cur.execute(
                "SELECT taxon, seq, part_name, part, aln_idx "
                "FROM [.partitiondata] "
                "ORDER BY aln_idx, part"):

            if prev_idx != aln_idx:
                fh, of = self._setup_newfile(fh, aln_idx, output_dir,
                                             suffix, output_file, ns)

                # Write file header
                fh.write("{}\n".format(len(self.partitions.partitions)))
                prev_idx = aln_idx

                if not fh:
                    continue

            # Write partition header at the beginning
            if prev_part != pidx:
                fh.write("\n{} {} {}\n".format(
                    pname,
                    len(self.partition_count[pname][0]),
                    self.partition_count[pname][1]))

                prev_part = pidx
            
            # Write sequence data
            fh.write("{}\t{}\n".format(taxon, seq))

        if fh:
            fh.close()

    def _write_ima2(self, suffix, output_file, **kwargs):
        
        tx_space_ima2 = kwargs.get("tx_space_ima2", 10)
        cut_space_ima2 = kwargs.get("cut_space_ima2", 8)
        ima2_params = kwargs.get("ima2_params", None)
        table_name = kwargs.get("table_name", None)
        ns = kwargs.get("ns_pipe", None)
        pbar = kwargs.get("pbar", None)
        output_dir = kwargs.get("output_dir", None)
        
        population_file = ima2_params[0]
        population_tree = ima2_params[1]
        mutational_model = ima2_params[2]
        inheritance_scalar = ima2_params[3]

        # Get information on which species belong to each population from
        # the populations file
        population_handle = open(population_file)
        population_storage = OrderedDict()
        for line in population_handle:
            taxon, population = line.strip().split()
            try:
                population_storage[population.strip()].append(taxon)
            except KeyError:
                population_storage[population.strip()] = [taxon]

        # File object that will be used to write sequence data
        fh = None

        # Check if partition data has already been built in the database.
        # If not, create it
        if not self.partition_data:
            self.partition_data = self._get_partition_data(table_name, ns,
                                                           pbar)

        # Stores the string of the last file
        prev_idx = ""
        prev_part = ""
        for taxon, seq, pname, pidx, aln_idx in self.cur.execute(
                "SELECT taxon, seq, part_name, part, aln_idx "
                "FROM [.partitiondata] "
                "ORDER BY aln_idx, part"):

            if prev_idx != aln_idx:
                fh, of = self._setup_newfile(fh, aln_idx, output_dir,
                                             suffix, output_file, ns)

                fh.write("Input file for IMa2 using %s alignments\n"
                         "{}\n"  # Line with number of loci
                         "{}\n"  # Line with name of populations
                         "{}\n"  # Line with population string
                         "".format(
                            len(self.partitions.partitions),
                            len(population_storage),
                            " ".join(population_storage.keys()),
                            population_tree))

                prev_idx = aln_idx

            if prev_part != pidx:

                # Get number of members for each population for the current
                # partition
                cur_tx = set(self.partition_count[pname][0])
                members = " ".join([str(len(cur_tx & set(x))) for x in
                                    population_storage.values()])

                fh.write("{} {} {} {} {}\n".format(
                    pname,
                    members,
                    self.partition_count[pname][1],
                    mutational_model,
                    inheritance_scalar))

                prev_part = pidx

            fh.write("{}{}\n".format(
                taxon[:cut_space_ima2].ljust(tx_space_ima2), seq))

        if fh:
            fh.close()

    def _write_mcmctree(self, suffix, output_file, **kwargs):

        tx_space_phy = kwargs.get("tx_space_phy", 40)
        cut_space_phy = kwargs.get("cut_space_phy", 39)
        table_name = kwargs.get("table_name", None)
        ns = kwargs.get("ns_pipe", None)
        pbar = kwargs.get("pbar", None)
        output_dir = kwargs.get("output_dir", None)

        # File object that will be used to write sequence data
        fh = None

        # Check if partition data has already been built in the database.
        # If not, create it
        if not self.partition_data:
            self.partition_data = self._get_partition_data(table_name, ns,
                                                           pbar)

        # Stores the string of the last file
        prev_idx = ""
        prev_part = ""
        for taxon, seq, pname, pidx, aln_idx in self.cur.execute(
                "SELECT taxon, seq, part_name, part, aln_idx "
                "FROM [.partitiondata] "
                "ORDER BY aln_idx, part"):

            if prev_idx != aln_idx:
                fh, of = self._setup_newfile(fh, aln_idx, output_dir,
                                             suffix, output_file, ns)
                prev_idx = aln_idx

                if not fh:
                    continue

            if prev_part != pidx:

                fh.write("{} {}\n".format(
                    len(self.partition_count[pname][0]),
                    self.partition_count[pname][1]))

                prev_part = pidx

            fh.write("{}  {}\n".format(
                taxon[:cut_space_phy].ljust(tx_space_phy),
                seq))

    def write_to_file(self, output_format, conversion_suffix="",
                      output_suffix="", *args, **kwargs):
        """Writes `Alignment` objects into files.

        Wrapper of the `write_to_file` method of the `Alignment` object.

        Parameters
        ----------
        output_format : list
            List with the output formats to generate. Options are:
            {"fasta", "phylip", "nexus", "stockholm", "gphocs", "ima2"}.
        output_file : str
            Name of the output file. If using `ns_pipe` in TriFusion,
            it will prompt the user if the output file already exists.
        tx_space_nex : int
            Space (in characters) provided for the taxon name in nexus
            format (default is 40).
        tx_space_phy : int
            Space (in characters) provided for the taxon name in phylip
            format (default is 40).
        tx_space_ima2 : int
            Space (in characters) provided for the taxon name in ima2
            format (default is 10).
        cut_space_nex : int
            Set the maximum allowed character length for taxon names in
            nexus format. Longer  names are truncated (default is 39).
        cut_space_phy : int
            Set the maximum allowed character length for taxon names in
            phylip format. Longer  names are truncated (default is 39).
        cut_space_ima2 : int
            Set the maximum allowed character length for taxon names in
            ima2 format. Longer  names are truncated (default is 8).
        interleave : bool
            Determines whether the output alignment
            will be in leave (False) or interleave (True) format. Not all
            output formats support this option.
        gap : str
            Symbol for alignment gaps (default is '-').
        model_phylip : str
            Substitution model for the auxiliary partition file of phylip
            format, compliant with RAxML.
        outgroup_list : list
            Specify list of taxon names that will be defined as the outgroup
            in the nexus output format. This may be useful for analyses with
            MrBayes or other software that may require outgroups.
        ima2_params : list
            A list with the additional information required for the ima2
            output format. The list should contains the following information:

               1. (str) path to file containing the species and populations.
               2. (str) Population tree in newick format, e.g. (0,1):2.
               3. (str) Mutational model for all alignments.
               4. (str) inheritance scalar.

        use_charset : bool
            If True, _partitions from the `Partitions` object will be written
            in the nexus output format (default is True).
        partition_file : bool
            If True, the auxiliary _partitions file will be written (default
            is True).
        output_dir : str
            If provided, the output file will be written on the specified
            directory.
        phy_truncate_names : bool
            If True, taxon names in phylip format are truncated to a maximum
            of 10 characters (default is False).
        ld_hat : bool
            If True, the fasta output format will include a first line
            compliant with the format of LD Hat and will truncate sequence
            names and sequence length per line accordingly.
        use_nexus_models : bool
            If True, writes the _partitions charset block in nexus format.
        ns_pipe : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.
        pbar : ProgressBar
            A ProgressBar object used to log the progress of TriSeq execution.
        table_name : string
            Name of the table from where the sequence data is fetched.
        """

        ns_pipe = kwargs.get("ns_pipe", None)
        pbar = kwargs.pop("pbar", None)
        output_file = kwargs.pop("output_file", None)

        write_methods = {
            "fasta": self._write_fasta,
            "phylip": self._write_phylip,
            "nexus": self._write_nexus,
            "stockholm": self._write_stockholm,
            "gphocs": self._write_gphocs,
            "ima2": self._write_ima2,
            "mcmctree": self._write_mcmctree
        }

        for fmt in output_format:

            filename = None

            if output_file:
                suffix = None
                filename = output_file + self.format_ext[fmt]
            else:
                suffix = conversion_suffix + output_suffix + \
                    self.format_ext[fmt]

            write_methods[fmt](suffix, filename, **kwargs)

    def get_gene_table_stats(self, active_alignments=None, sortby=None,
                             ascending=True):
        """Returns summary statistics for each `Alignment`.

        Gets summary statistics for each individual gene for the gene table
        view of TriFusion. Summary statistics have already been calculated in
        get_summary_stats, so here we only get the active alignments for
        showing. Therefore, this function should only be called after
        get_summary_stats.

        Parameters
        ----------
        active_alignments : list
            List of `Alignment.name` objects to show
        sortby : str
            Table header used for sorting. Can be {"nsites", "taxa", "var",
            "inf", "gap", "missing"}.
        ascending : bool
            If True, the `sortby` column will be sorted in ascending order
            (default is True).

        Returns
        -------
        summary_gene_table : pandas.DataFrame
            DataFrame with summary statistics for each alignment.
        table : list
            List with the data for creating .csv tables.
        """

        # Set table header
        table = [["Gene name", "Number of sites", "Number of taxa",
                  "Variable sites", "Informative sites", "Gaps",
                  "Missing data"]]
        # Table line keys matching summary_gene_table dict values
        tl = ["nsites", "taxa", "var", "inf", "gap", "missing"]

        # Update the active alignments
        if active_alignments and active_alignments != list(
                self.alignments.keys()):
            self.update_active_alignments(active_alignments)

        # Filter table with active alignments
        summary_gene_table = self.summary_gene_table[
            self.summary_gene_table["genes"].isin(self.aln_names())]

        if not sortby:
            sortby = ["genes"]

        summary_gene_table.sort_values(sortby, ascending=ascending,
                                       inplace=True)

        # Populate table information
        for k in self.summary_gene_table.itertuples():
            # Add table line
            table.append(list(k[1:]))

        return summary_gene_table, table

    # Stats methods
    def get_summary_stats(self, active_alignments=None, ns=None):
        """Calculates summary statistics for the 'active' alignments.

        Creates/Updates summary statistics for the active alignments.

        Parameters
        ----------
        active_alignments : list
            List of `Alignment` objects for which the summary statistics
            will be calculated
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        summary_stats : dict
            Dictionary with the overall summary statistics
        table : list
            List with overall summary statistics for creating .csv tables.
        """

        def add_data():
            # Get values for current alignment for average calculations
            self.summary_stats["avg_gaps"].append(cur_gap)
            self.summary_stats["avg_missing"].append(cur_missing)
            self.summary_stats["avg_var"].append(cur_var)
            self.summary_stats["avg_inf"].append(cur_inf)

            if not any(self.summary_gene_table["genes"] == aln.name):
                # Get row information for current gene in dict format
                row_info = (aln.name,
                            aln.locus_length,
                            len(aln.taxa_idx),
                            cur_var,
                            cur_inf,
                            cur_gap,
                            cur_missing)
                # Create DataFrame from dict
                row_dt = pd.DataFrame(
                    [row_info], index=[c],
                    columns=self.summary_gene_table.columns)
                # Add row DF to main DF
                self.summary_gene_table = pd.concat(
                    [self.summary_gene_table,
                     row_dt])


        # Update active alignments if they changed since last update
        if active_alignments and \
                active_alignments != list(self.alignments.keys()):
            self.update_active_alignments(active_alignments)

        self._set_pipes(ns, None, total=len(self.alignments))
        c = -1

        # Set table header for summary_stats
        table = [["Genes", "Taxa", "Alignment length", "Gaps",
                  "Gaps per gene", "Missing data", "Missing data per gene",
                  "Variable sites", "Variable sites per gene",
                  "Informative sites", "Informative sites per gene"]]
        # Table line keys matching summary_stats for table completion
        tl = ["genes", "taxa", "seq_len", "gaps", "avg_gaps", "missing",
              "avg_missing", "variable", "avg_var", "informative", "avg_inf"]

        # Reset summary_stats
        self._reset_summary_stats()

        # Get number of alignments
        self.summary_stats["genes"] = len(self.alignments)

        # Get number of taxa
        self.summary_stats["taxa"] = len(self.taxa_names)

        # Get statistics that require iteration over alignments
        prev_idx = ""
        for col, aln_idx in self.iter_columns():

            if prev_idx != aln_idx:

                if prev_idx:
                    add_data()

                self._update_pipes(ns, None, value=c)
                c += 1

                # Get current alignment
                aln = self.alignment_idx[aln_idx]
                self.summary_stats["seq_len"] += aln.locus_length

                cur_gap, cur_missing = 0, 0
                cur_var, cur_inf = 0, 0

                prev_idx = aln_idx

                if ns:
                    if ns.stop:
                        raise KillByUser("")

            col = Counter(col)

            # Get missing data and gaps
            if self.sequence_code[1] in col:
                self.summary_stats["missing"] += 1
                cur_missing += 1
            if self.gap_symbol in col:
                self.summary_stats["gaps"] += 1
                cur_gap += 1

            # Get variability information
            # Filter missing data
            for i in [self.sequence_code[1], self.gap_symbol]:
                del col[i]

            # If it's only missing data, ignore
            if col:
                # Get variable sites
                if len(col) > 1:
                    self.summary_stats["variable"] += 1
                    cur_var += 1

                # If any of the remaining sites is present in more than two
                # taxa score the site as informative
                if len([x for x in col.values() if x >= 2]) >= 2:
                    self.summary_stats["informative"] += 1
                    cur_inf += 1
        if aln_idx:
            add_data()

        # Get average values
        for k in ["avg_gaps", "avg_missing", "avg_var", "avg_inf"]:
            self.summary_stats[k] = round(np.mean(self.summary_stats[k]))

        # Get percentage values for missing data
        total = self.summary_stats["seq_len"] * self.summary_stats["taxa"]
        for k in ["gaps", "missing"]:
            n = self.summary_stats[k]
            percentage = round((float(n) / float(total)) * 100, 2)
            self.summary_stats[k] = "{0:,} ({1}%)".format(n, percentage)

        # Get percentage values for variation
        for k in ["variable", "informative"]:
            n = self.summary_stats[k]
            percentage = round((float(n) /
                                float(self.summary_stats["seq_len"])) * 100, 2)
            self.summary_stats[k] = "{0:,} ({1}%)".format(n, percentage)

        # Complete table information
        table.append([self.summary_stats[x] for x in tl])

        return dict(self.summary_stats), table

    @check_data
    def gene_occupancy(self, ns=None):
        """Create data for gene occupancy plot.

        Creates data for an interpolation plot to visualize the amount of
        missing genes in a alignment dataset.

        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title
        """

        data = []

        for alignment in self.alignments.values():
            if ns:
                if ns.stop:
                    raise KillByUser()

            data.append([1 if x in alignment.taxa_idx else 0
                         for x in self.taxa_names])

        data = np.transpose(data)

        return {"data": data,
                "ax_names": ["Genes", "Taxa"],
                "title": "Gene occupancy"}

    @check_data
    def missing_data_distribution(self, ns=None):
        """Create data for overall distribution of missing data plot.

        Creates data for overall distribution of missing data. This will
        calculate the amount of gaps, missing and actual data.

        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "legend": mpl Legend object,

            "table_header": list with headers of table.
        """

        def add_data():
            data_storage["Gaps"].append(np.mean(gaps_g))
            data_storage["Missing data"].append(np.mean(missing_g))
            data_storage["Data"].append(np.mean(data_g))

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        legend = ["Gaps", "Missing data", "Data"]

        data_storage = OrderedDict((x, []) for x in legend)

        prev_idx = ""
        for taxon, seq, aln_idx in self.iter_alignments():

            if prev_idx != aln_idx:
                
                if prev_idx:
                    # Add data to taxa missing in the current alignment
                    missing = (x for x in self.taxa_names
                               if x not in aln.taxa_idx)
                    for tx in missing:
                        gaps_g.append(0)
                        missing_g.append(1)
                        data_g.append(0)

                    add_data()

                self._update_pipes(ns, None, value=c)
                c += 1

                # Get alignment object
                aln = self.alignment_idx[aln_idx]

                gaps_g, missing_g, data_g = [], [], []

                prev_idx = aln_idx

            # Get gaps
            gaps = float(seq.count("-"))
            gaps_g.append(gaps / float(aln.locus_length))
            # Get missing data
            missing = float(seq.count(aln.sequence_code[1]))
            missing_g.append(missing / float(aln.locus_length))
            # Get actual data
            actual_data = (float(aln.locus_length) - gaps - missing) / \
                float(aln.locus_length)
            data_g.append(actual_data)

        if aln_idx:
            add_data()

        data = [x for x in data_storage.values()]

        return {"data": data,
                "title": "Distribution of missing data",
                "legend": legend,
                "ax_names": ["Proportion", "Number of genes"],
                "table_header": ["Bin"] + legend}

    @check_data
    def missing_data_per_species(self, ns=None):
        """Creates data for missing data proportion per species plot.

        For each taxon in `taxon_list` calculate the proportion of gap,
        missing and actual data.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "legend": mpl Legend object,

            "table_header": list with headers of table,

            "normalize": bool, whether data values should be normalized,

            "normalize_factor": int, factor to use in normalization
        """

        # Data for a stacked bar plot. First element for gaps, second for
        # missing, third for actual data
        data_storage = OrderedDict((taxon, [0, 0, 0]) for taxon in
                                   self.taxa_names)
        total_len = 0

        legend = ["Gaps", "Missing", "Data"]

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        prev_idx = ""
        for taxon, seq, aln_idx in self.iter_alignments():

            if prev_idx != aln_idx:

                if prev_idx:
                    # Add data to taxa missing in the current alignment
                    missing = (x for x in self.taxa_names
                               if x not in aln.taxa_idx)
                    for tx in missing:
                        data_storage[tx][1] += aln.locus_length

                self._update_pipes(ns, None, value=c)
                c += 1

                aln = self.alignment_idx[aln_idx]
                total_len += aln.locus_length

                prev_idx = aln_idx

            # Get gaps
            gaps = seq.count("-")
            data_storage[taxon][0] += gaps
            # Get missing
            missing = seq.count(aln.sequence_code[1])
            data_storage[taxon][1] += missing
            # Get actual data
            actual_data = aln.locus_length - gaps - missing
            data_storage[taxon][2] += actual_data

        data_storage = OrderedDict(sorted(data_storage.items(),
                                          key=lambda x: x[1][1] + x[1][0],
                                          reverse=True))

        data = np.array([[float(x[0]) for x in
                          data_storage.values()],
                         [float(x[1]) for x in
                          data_storage.values()],
                         [float(x[2]) for x in
                          data_storage.values()]])

        return {"data": data,
                "title": "Distribution of missing data per species",
                "labels": list(data_storage.keys()),
                "legend": legend,
                "ax_names": [None, "Frequency"],
                "table_header": ["Taxon", "Gaps", "%", "Missing", "%", "Data",
                                 "%"],
                "normalize": True,
                "normalize_factor": total_len}

    @check_data
    def missing_genes_per_species(self, ns=None):
        """Creates data with distribution of missing genes per species plot.

        Calculates, for each taxon, the number of alignments where it is
        present.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "labels": list, label for xticks,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "table_header": list with headers of table
        """

        data_storage = OrderedDict((taxon, 0) for taxon in self.taxa_names)

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        for aln in self.alignments.values():

            self._update_pipes(ns, None, value=c)
            c += 1

            for key in data_storage:
                if key not in aln.taxa_idx:
                    data_storage[key] += 1

        # Sort data in descending order of missing genes
        data_storage = OrderedDict(sorted(data_storage.items(), reverse=True,
                                          key=lambda t: t[1]))

        return {"data": [list(data_storage.values())],
                "labels": list(data_storage.keys()),
                "title": "Distribution of missing genes per species",
                "ax_names": [None, "Frequency"],
                "table_header": ["Taxon", "Missing genes"]
                }

    @check_data
    def missing_genes_average(self, ns=None):
        """Creates histogram data with average missing taxa.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "table_header": list with headers of table
        """

        data = []

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        for aln in self.alignments.values():

            self._update_pipes(ns, None, value=c)
            c += 1

            data.append(len(set(self.taxa_names) - set(aln.taxa_idx.keys())))

        return {"data": data,
                "title": "Distribution of missing taxa",
                "ax_names": ["Number of missing taxa", "Frequency"],
                "table_header": ["Number of missing taxa", "Frequency"]}

    @check_data
    def average_seqsize_per_species(self, ns=None):
        """Creates data for average sequence size per taxon.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "labels": list, label for xticks,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title
        """

        data_storage = OrderedDict((taxon, []) for taxon in self.taxa_names)

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        prev_idx = ""
        for taxon, seq, aln_idx in self.iter_alignments():

            if prev_idx != aln_idx:

                self._update_pipes(ns, None, value=c)
                c += 1

                prev_idx = aln_idx

            data_storage[taxon].append(
                len(seq.replace("-", "").replace(self.sequence_code[1], "")))

        # Adapt y-axis label according to sequence code
        seq_code = self.sequence_code[0]
        ax_ylabel = "Size (bp)" if seq_code == "DNA" else "Size (residues)"

        data_storage = OrderedDict(sorted(data_storage.items(), reverse=True,
                                   key=lambda t: np.mean(t[1])))

        return {"data": list(data_storage.values()),
                "labels": list(data_storage.keys()),
                "title": "Sequence size distribution per species",
                "ax_names": [None, ax_ylabel]}

    @check_data
    def average_seqsize(self, ns=None):
        """Creates data for the average sequence size for the entire data set.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "table_header": list with headers of table
        """

        data_storage = []

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        for aln in self.alignments.values():

            self._update_pipes(ns, None, value=c)
            c += 1

            data_storage.append(aln.locus_length)

        # Adapt y-axis label according to sequence code
        seq_code = aln.sequence_code[0]
        ax_xlabel = "Size (bp)" if seq_code == "DNA" else "Size (residues)"

        return {"data": data_storage,
                "title": "Average sequence size distribution",
                "ax_names": [ax_xlabel, "Frequency"],
                "table_header": [ax_xlabel, "Frequency"]}

    @check_data
    def characters_proportion(self, ns=None):
        """Creates data for the proportion of nt/aa for the data set.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "labels": list, label for xticks,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "table_header": list with headers of table,
        """

        data_storage = Counter()

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        prev_idx = ""
        for _, seq, aln_idx in self.iter_alignments():

            if prev_idx != aln_idx:

                self._update_pipes(ns, None, value=c)
                c += 1

                prev_idx = aln_idx

            data_storage += Counter(
                seq.replace("-", "").replace(self.sequence_code[1], ""))

        # Determine total number of characters
        chars = float(sum(data_storage.values()))

        # Valid characters list
        valid_chars = dna_chars if self.sequence_code[0] == "DNA" else \
            list(aminoacid_table.keys())

        data, xlabels = zip(*[(float(x) / chars, y.upper()) for y, x in
                              data_storage.items() if y in valid_chars])

        title = "Nucleotide proportions" if self.sequence_code[0] == "DNA" \
            else "Amino acid proportions"
        ax_xlabel = "Nucleotide" if self.sequence_code[0] == "DNA" \
            else "Amino acid"

        return {"data": [data],
                "labels": xlabels,
                "title": title,
                "ax_names": [ax_xlabel, "Proportion"],
                "table_header": [ax_xlabel, "Proportion"]}

    @check_data
    def characters_proportion_per_species(self, ns=None):
        """Creates data for the proportion of nt/aa per species.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "labels": list, label for xticks,

            "ax_names": 2 element list with axis labels [x, y],

            "legend": mpl Legend object,

            "title": str with title,

            "table_header": list with headers of table}
        """

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        data_storage = OrderedDict((x, Counter()) for x in self.taxa_names)

        prev_idx = ""
        for taxon, seq, aln_idx in self.iter_alignments():

            if prev_idx != aln_idx:

                self._update_pipes(ns, None, value=c)
                c += 1

                prev_idx = aln_idx

            data_storage[taxon] += Counter(
                seq.replace("-", ""). replace(self.sequence_code[1], ""))

        legend = dna_chars if self.sequence_code[0] == "DNA" else \
            list(aminoacid_table.keys())

        data = [[] for _ in legend]

        for p, char in enumerate(legend):
            for c in data_storage.values():
                chars = float(sum([x for y, x in c.items() if y in legend]))
                data[p].append(float(c[char]) / chars)

        data = np.array(data)

        title = "Nucleotide proportions" if self.sequence_code[0] == "DNA" \
            else "Amino acid proportions"

        ax_ylabel = "Nucleotide" if self.sequence_code[0] == "DNA" \
            else "Amino acid"

        return {"data": data,
                "title": title,
                "labels": list(data_storage.keys()),
                "legend": legend,
                "ax_names": ["Taxa", ax_ylabel],
                "table_header": ["Taxon"] + legend}

    @LookupDatabase
    def _get_similarity(self, seq1=None, seq2=None, aln_len=None):
        """Gets the similarity between two strings and the effective length.

        Compares two sequences and calculates the number of similarities.
        This comparison ignores columns with missing data, so the final
        sequence length that is actually compared has to be adjusted.

        Parameters
        ----------
        seq1, seq2 : str
            Sequence strings.
        aln_len : int
            Original lenght of `seq1` and `seq2`.

        Returns
        -------
        sim : float
            Number of pairwise similarities.
        ef_len : float
            Effective sequence length, without missing data.

        Notes
        -----
        To avoid performing the calculations between identical string pairs,
        this method is decorated with a class that stores the calculation
        result of new string pairs in a database, and fetches that result
        when the same string pair appears later on.
        """

        seq1 = np.array(list(seq1))
        seq2 = np.array(list(seq2))

        ef_len_ = (seq1 != self.sequence_code[1]) & \
            (seq1 != self.gap_symbol) & \
            (seq2 != self.sequence_code[1]) & \
            (seq2 != self.gap_symbol)

        sim_ = (seq1 == seq2) & ef_len_

        sim = np.count_nonzero(sim_)
        ef_len = np.count_nonzero(ef_len_)

        return float(sim), float(ef_len)

    @check_data
    def sequence_similarity(self, ns=None):
        """Creates data for average sequence similarity plot.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y]
        """

        self._get_similarity("connect")

        data = []

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        for aln in self.alignments.values():

            self._update_pipes(ns, None, value=c)
            c += 1

            aln_similarities = []

            for seq1, seq2 in itertools.combinations(aln.iter_sequences(), 2):

                self._check_killswitch(ns)

                sim, total_len = self._get_similarity(seq1, seq2,
                                                      aln.locus_length)

                if total_len:
                    aln_similarities.append(sim / total_len)

            if aln_similarities:
                data.append(np.mean(aln_similarities) * 100)

        self._get_similarity("disconnect")

        return {"data": data,
                "ax_names": ["Similarity (%)", "Frequency"]}

    @check_data
    def sequence_similarity_per_species(self, ns=None):
        """Creates data for sequence similarity per species pair plot.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "labels": list, label for xticks,

            "color_label": str, label for colorbar
        """

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        self._get_similarity("connect")

        # Create matrix for parwise comparisons
        data = [np.empty((len(self.taxa_names), 0)).tolist() for _ in
                range(len(self.taxa_names))]

        taxa_pos = OrderedDict((x, y) for y, x in enumerate(self.taxa_names))

        for aln in self.alignments.values():

            self._update_pipes(ns, None, value=c)
            c += 1

            for tx1, tx2 in itertools.combinations(taxa_pos.keys(), 2):

                self._check_killswitch(ns)

                try:
                    seq1, seq2 = aln.get_sequence(tx1), aln.get_sequence(tx2)
                except KeyError:
                    continue

                sim, l = self._get_similarity(seq1, seq2, aln.locus_length)
                if l:
                    data[taxa_pos[tx1]][taxa_pos[tx2]].append(sim / l)

        data = np.array([[np.mean(y) if y else 0. for y in x] for x in data])
        mask = np.tri(data.shape[0], k=0)
        data = np.ma.array(data, mask=mask)

        self._get_similarity("disconnect")

        return {"data": data,
                "color_label": "Pairwise sequence similarity",
                "labels": list(taxa_pos)}

    @check_data
    def sequence_similarity_gene(self, gene_name, window_size, ns=None):
        """Creates data for sliding window sequence similarity for alignment.

        Retrieves an alignment using `gene_name` and calculates the
        average pair-wise sequence similarity along the alignment length
        using a sliding window approach.

        Parameters
        ----------
        gene_name : str
            `Alignment.name` from an alignment
        window_size : int
            Size of sliding window.
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "window_size": int, size of sliding window,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "table_header": list with headers of table
        """

        aln_obj = self.retrieve_alignment(gene_name)

        data = []

        for i in range(0, aln_obj.locus_length, window_size):

            window_similarities = []

            seqs = np.array([[y for y in x[i:i + window_size]] for x in
                             aln_obj.iter_sequences()])

            for seq1, seq2 in itertools.combinations(seqs, 2):

                s, t = self._get_similarity("".join(seq1), "".join(seq2),
                                            window_size)
                if t:
                    sim = s / t
                else:
                    sim = 0

                window_similarities.append(sim * 100)

            if window_similarities:
                data.append(np.mean(window_similarities))

        return {"data": data,
                "title": "Sequence similarity sliding window for gene\n %s"
                         % basename(gene_name),
                "window_size": window_size,
                "ax_names": ["Sequence (bp)", "Similarity (%)"],
                "table_header": ["Sequence (bp)", "Similarity (%)"]}

    @check_data
    def sequence_segregation(self, ns=None, proportions=False):
        """Creates data for overall distribution of segregating sites.

        Parameters
        ----------
        proportions : bool
            If True, use proportions instead of absolute values
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "table_header": list with headers of table,

            "real_bin_num":
        """
        
        def add_data():
            if proportions:
                data.append(
                    float(segregating_sites / float(aln.locus_length)))
            else:
                data.append(segregating_sites)

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        data = []

        prev_idx = ""
        for column, aln_idx in self.iter_columns():

            if prev_idx != aln_idx:
                
                if prev_idx:
                    add_data()

                self._update_pipes(ns, None, value=c)
                c += 1

                aln = self.alignment_idx[aln_idx]

                segregating_sites = 0

                prev_idx = aln_idx

            # Remove gaps and missing characters
            column = set([x for x in column if x != aln.sequence_code[1]
                          and x != self.gap_symbol])

            if len(column) > 1:
                segregating_sites += 1

        if proportions:
            ax_names = ["Segregating sites", "Percentage"]
            real_bin = False
        else:
            ax_names = ["Segregating sites", "Frequency"]
            real_bin = True

        return {"data": data,
                "ax_names": ax_names,
                "title": "Distribution of segregating sites",
                "table_header": ax_names,
                "real_bin_num": real_bin}

    @check_data
    def sequence_segregation_per_species(self, ns=None):
        """Creates data for a triangular matrix of sequence segregation for
        pairs of taxa.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "labels": list, label for xticks,

            "color_label": str, label for colorbar
        """

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        self._get_similarity("connect")

        # Create matrix for parwise comparisons
        data = [np.empty((len(self.taxa_names), 0)).tolist() for _ in
                range(len(self.taxa_names))]

        taxa_pos = OrderedDict((x, y) for y, x in enumerate(self.taxa_names))

        for aln in self.alignments.values():

            self._update_pipes(ns, None, value=c)
            c += 1

            for tx1, tx2 in itertools.combinations(taxa_pos.keys(), 2):

                self._check_killswitch(ns)

                try:
                    seq1, seq2 = aln.get_sequence(tx1), aln.get_sequence(tx2)
                except KeyError:
                    continue

                s, t = self._get_similarity(seq1, seq2, aln.locus_length)
                if t:
                    aln_diff = t - s
                else:
                    aln_diff = 0

                data[taxa_pos[tx1]][taxa_pos[tx2]].append(aln_diff)

        data = np.array([[np.mean(y) if y else 0. for y in x] for x in data])
        mask = np.tri(data.shape[0], k=0)
        data = np.ma.array(data, mask=mask)

        self._get_similarity("disconnect")

        return {"data": data,
                "labels": list(taxa_pos),
                "color_label": "Segregating sites"}

    @check_data
    def sequence_segregation_gene(self, gene_name, window_size, ns=None):
        """Create data for a sliding window analysis of segregating sites.

        Retrieves an alignment using `gene_name` and calculates the
        number of segregating sites along the alignment length
        using a sliding window approach.

        Parameters
        ----------
        gene_name : str
            `Alignment.name` from an alignment
        window_size : int
            Size of sliding window.
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "window_size": int, size of sliding window,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "table_header": list with headers of table
        """

        aln_obj = self.retrieve_alignment(gene_name)

        data = []

        for i in range(0, aln_obj.locus_length, window_size):

            segregating_sites = 0

            seqs = np.array([[y for y in x[i:i + window_size]] for x in
                             aln_obj.iter_sequences()])

            for column in zip(*seqs):
                column = set([x for x in column if x != aln_obj.sequence_code[1]
                              and x != "-"])

                if len(column) > 1:
                    segregating_sites += 1

            data.append(segregating_sites)

        return {"data": data,
                "title": "Number of segregating sites sliding window for "
                         "gene\n %s" % basename(gene_name),
                "window_size": window_size,
                "ax_names": ["Sequence (bp)", "Segregating sites"],
                "table_header": ["Sequence (bp)", "Segregating sites"]}

    @check_data
    def length_polymorphism_correlation(self, ns=None):
        """Creates data for correlation between alignment length and
        informative sites.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "correlation": bool, whether to perform spearman correlation,

            "title": str with title,

            "table_header": list with headers of table
        """

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 0

        data_length = []
        data_inf = []

        missing_symbols = [self.sequence_code[1], self.gap_symbol]

        prev_idx = ""
        for column, aln_idx in self.iter_columns():

            if prev_idx != aln_idx:

                if prev_idx:
                    data_length.append(aln.locus_length)
                    data_inf.append(informative_sites)

                self._update_pipes(ns, None, value=c)
                c += 1

                aln = self.alignment_idx[aln_idx]
                informative_sites = 0

                prev_idx = aln_idx

            column = Counter(column)
            for sym in missing_symbols:
                del column[sym]

            if len(column) > 1:
                # Remove most common and check the length of the remaining
                del column[column.most_common()[0][0]]
                if sum(column.values()) > 1:
                    informative_sites += 1

        return {"data": [data_length, data_inf],
                "title": "Correlation between alignment length and number of "
                         "variable sites",
                "ax_names": ["Alignment length", "Informative sites"],
                "table_header": ["Alignment length", "Informative sites"],
                "correlation": True}

    @check_data
    def allele_frequency_spectrum(self, ns=None, proportions=False):
        """Creates data for the allele frequency spectrum.

        Calculates the allele frequency spectrum of the entire
        alignment data set. Here, multiple alignments are effectively treated
        as a single one. This method is exclusive of DNA sequence type and
        supports IUPAC ambiguity codes.

        Parameters
        ----------
        proportions : bool
            If True, use proportions instead of absolute values
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "table_header": list with headers of table
        """

        # Make check for sequence type consistency here for TriStats.py. In
        # TriFusion the check is made before calling this method
        if self.sequence_code[0] != "DNA":
            return {"exception": InvalidSequenceType}

        data = []
        missing_symbols = [self.sequence_code[1], self.gap_symbol]

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 0

        prev_idx = ""
        for column, aln_idx in self.iter_columns():

            if prev_idx != aln_idx:

                self._update_pipes(ns, None, value=c)
                c += 1

                prev_idx = aln_idx

            column = Counter(column)

            for sym in missing_symbols:
                del column[sym]
            # Get length
            col_len = sum(column.values())

            # Consider only bi-allelic SNPs
            if len(column) == 2:
                # Remove most common and check the length of the remaining
                del column[column.most_common()[0][0]]
                # Append number of derived alleles
                if proportions:
                    data.append(float(sum(column.values())) /
                                float(col_len))
                else:
                    data.append(sum(column.values()))

        return {"data": data,
                "title": "Allele frequency spectrum",
                "ax_names": ["Derived allele frequency", "Frequency"],
                "table_header": ["Derived allele frequency", "Frequency"],
                "real_bin_num": True}

    @check_data
    def allele_frequency_spectrum_gene(self, gene_name, ns):
        """Creates data for the allele frequency spectrum of an `Alignment`.

        Parameters
        ----------
        gene_name : str
            `Alignment.name` of an alignment.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "table_header": list with headers of table,

            "real_bin_num":
        """

        # Make check for sequence type consistency here for TriStats.py. In
        # TriFusion the check is made before calling this method
        if self.sequence_code[0] != "DNA":
            return {"exception": InvalidSequenceType}

        aln = self.retrieve_alignment(gene_name)
        data = []

        for column, _ in self.iter_columns(aln_idx=aln.db_idx):

            # Remove gaps and missing characters
            column = Counter([iupac_conv[x] for x in column if
                              x != aln.sequence_code[1] and
                              x != self.gap_symbol])

            # Consider only bi-allelic SNPs
            if len(column) == 2:
                # Remove most common and check the length of the remaining
                del column[column.most_common()[0][0]]
                # Append number of derived alleles
                data.append(sum(column.values()))

        return {"data": data,
                "title": "Allele frequency spectrum",
                "ax_names": ["Derived allele frequency", "Frequency"],
                "table_header": ["Derived allele frequency", "Frequency"],
                "real_bin_num": True}

    @check_data
    def taxa_distribution(self, ns=None):
        """Creates data for a distribution of taxa frequency across alignments.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "labels": list, label for xticks,

            "title": str with title,

            "color_label": str, label for colorbar,

            "real_bin_num":
        """

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        data = []

        for aln in self.alignments.values():

            self._update_pipes(ns, None, value=c)
            c += 1

            # Get number of taxa
            data.append(len(aln.taxa_idx))

        return {"data": data,
                "title": "Distribution of taxa frequency",
                "ax_names": ["Number of taxa", "Frequency"],
                "table_header": ["Number of taxa", "Frequency"],
                "real_bin_num": True}

    @check_data
    def cumulative_missing_genes(self, ns=None):
        """Creates data for a distribution of the maximum number of genes
        available for consecutive thresholds of missing data.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "labels": list, label for xticks,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "table_header": list with headers of table}
        """

        size_storage = []
        data = []

        self._set_pipes(ns, None, total=self.alignments)
        c = 1

        # total number of taxa in data set
        taxa = float(len(self.taxa_names))

        for aln in self.alignments.values():

            self._update_pipes(ns, None, value=c)
            c += 1

            # Get number of taxa
            size_storage.append((float(len(aln.taxa_idx)) / taxa) * 100)

        labels = []
        for i in xrange(0, 105, 5):

            # Get percentage
            data.append(len([x for x in size_storage if x > i]))
            labels.append(str(i))

        return {"data": [data],
                "title": "Gene frequency for decreasing values of missing "
                         "data",
                "labels": labels,
                "ax_names": ["Minimum taxa representation", "Gene frequency"],
                "table_header": ["Percentage", "Frequency"]}

    @staticmethod
    def _mad_based_outlier(p, threshold=3.5):
        """Calculate MAD based outliers.

        An outlier detection method based on median absolute deviation. This
        code was adapted from http://stackoverflow.com/a/22357811/1990165.
        The usage of the median is much less biased that the mean and is
        robust to smaller data sets.

        Parameters
        ----------
        p : numpy.array
            Array with data observations.
        threshold : float
            Modified Z-score to use as a threshold.
        """

        if len(p.shape) == 1:
            p = p[:, None]
        # Get sample median.
        median = np.median(p, axis=0)
        diff = np.sum((p - median) ** 2, axis=-1)
        diff = np.sqrt(diff)
        med_abs_deviation = np.median(diff)

        # The 0.6745is a constant that makes values roughly equivalent in
        # units to standard deviations.
        z_score = 0.6745 * diff / med_abs_deviation

        return z_score > threshold

    @check_data
    def outlier_missing_data(self, ns=None):
        """Create data for missing data alignment outliers.

        Get data for outlier detection of genes based on the distribution of
        average missing data per gene. Data points will be based on the
        proportion of missing data symbols out of the possible total. For
        example, in an alignment with three taxa, each with 100 sites,
        the total possible missing data is 300 (100 * 3). Here, missing data
        will be gathered from all taxa and a proportion will be calculated
        based n the total possible

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "outliers": list of outlier data points,

            "outliers_labels": list of outlier labels
        """

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        data_labels = []
        data_points = []

        prev_idx = ""
        for taxon, seq, aln_idx in self.iter_alignments():

            if prev_idx != aln_idx:

                if prev_idx:
                    m_data = float(gn_data) / float(total_len)
                    data_points.append(m_data)
                    data_labels.append(aln.sname)

                self._update_pipes(ns, None, value=c)
                c += 1

                aln = self.alignment_idx[aln_idx]

                total_len = aln.locus_length * len(aln.taxa_idx)
                gn_data = 0

                prev_idx = aln_idx

            gn_data += seq.count(self.sequence_code[1]) + \
                seq.count(self.gap_symbol)

        if prev_idx:
            m_data = float(gn_data) / float(total_len)
            data_points.append(m_data)
            data_labels.append(aln.sname)

        data_points = np.asarray(data_points)
        data_labels = np.asarray(data_labels)

        # Get outliers
        outliers_points = data_points[self._mad_based_outlier(data_points)]
        # Get outlier labels
        outliers_labels = list(data_labels[self._mad_based_outlier(
            data_points)])

        return {"data": data_points,
                "title": "Missing data outlier gene detection",
                "outliers": outliers_points,
                "outliers_labels": outliers_labels,
                "ax_names": ["Proportion of missing data", "Frequency"]}

    @check_data
    def outlier_missing_data_sp(self, ns=None):
        """Create data for missing data taxa outliers.

        Gets data for outlier detection of species based on missing data. For
        this analysis, genes for which a taxa is completely absent will be
        ignored for the calculations of that taxa. The reason for this,
        is that including genes where the taxon is absent would bias the
        outlier detection towards taxa that have low prevalence in the data
        set, even if they have low missing data in the alignments where they
        are present.

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "outliers": list of outlier data points,

            "outliers_labels": list of outlier labels}
        """

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        data = dict((tx, []) for tx in self.taxa_names)

        prev_idx = ""
        for taxon, seq, aln_idx in self.iter_alignments():

            if prev_idx != aln_idx:

                self._update_pipes(ns, None, value=c)
                c += 1

                aln = self.alignment_idx[aln_idx]
                total_len = aln.locus_length
                
                prev_idx = aln_idx

            m_data = float(seq.count(self.sequence_code[1]) +
                           seq.count(self.gap_symbol)) / float(total_len)
            data[taxon].append(m_data)

        # Get average for each taxon
        for tx, vals in data.items():
            data[tx] = np.mean(vals)

        # Prepara data for plotting
        data_points = []
        data_labels = []
        for tx, vals in data.items():
            data_points.append(vals)
            data_labels.append(tx)

        data_points = np.asarray(data_points)
        data_labels = np.asarray(data_labels)

        # Get outliers
        outliers_points = data_points[self._mad_based_outlier(data_points)]
        # Get outlier taxa
        outlier_labels = list(data_labels[self._mad_based_outlier(data_points)])

        return {"data": data_points,
                "title": "Missing data outlier taxa detection",
                "outliers": outliers_points,
                "outliers_labels": outlier_labels,
                "ax_names": ["Proportion of missing data", "Frequency"]}

    @check_data
    def outlier_segregating(self, ns=None):
        """Creates data for segregating site alignment outliers.

        Generates data for the outlier detection of genes based on
        segregating sites. The data will be based on the number of alignments
        columns with a variable number of sites, excluding gaps and missing
        data

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "outliers": list of outlier data points,

            "outliers_labels": list of outlier labels
        """

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        data_points = []
        data_labels = []
        missing_symbols = [self.sequence_code[1], self.gap_symbol]

        prev_idx = ""
        for column, aln_idx in self.iter_columns():

            if prev_idx != aln_idx:

                if prev_idx:
                    # Get proportion of segregating sites for current alignment
                    data_points.append(float(segregating_sites) /
                                       float(aln.locus_length))
                    data_labels.append(aln.name)

                self._update_pipes(ns, None, value=c)
                c += 1

                aln = self.alignment_idx[aln_idx]

                segregating_sites = 0
                
                prev_idx = aln_idx

            column = list(x for x in set(column) if x not in missing_symbols)

            if len(column) > 1:
                segregating_sites += 1

        if prev_idx:
            # Get proportion of segregating sites for current alignment
            data_points.append(float(segregating_sites) /
                               float(aln.locus_length))
            data_labels.append(aln.name)

        data_points = np.asarray(data_points)
        data_labels = np.asarray(data_labels)

        # Get outliers
        outliers_points = data_points[self._mad_based_outlier(data_points)]
        # Get outlier taxa
        outlier_labels = list(data_labels[self._mad_based_outlier(data_points)])

        return {"data": data_points,
                "title": "Sequence variation outlier gene detection",
                "outliers": outliers_points,
                "outliers_labels": outlier_labels,
                "ax_names": ["Proportion of segregating sites", "Frequency"]}

    @check_data
    def outlier_segregating_sp(self, ns=None):
        """Create data for segregating sites taxa outliers.

        Generates data for the outlier detection of species based on their
        average pair-wise proportion of segregating sites. Comparisons
        with gaps or missing data are ignored

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "outliers": list of outlier data points,

            "outliers_labels": list of outlier labels
        """

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        self._get_similarity("connect")

        data = OrderedDict((tx, []) for tx in self.taxa_names)

        for aln in self.alignments.values():

            self._update_pipes(ns, None, value=c)
            c += 1

            for tx1, tx2 in itertools.combinations(data.keys(), 2):

                self._check_killswitch(ns)

                try:
                    seq1, seq2 = aln.get_sequence(tx1), aln.get_sequence(tx2)
                except KeyError:
                    continue

                s, t_len = self._get_similarity(seq1, seq2, aln.locus_length)

                if t_len:
                    s_data = (t_len - s) / t_len
                else:
                    s_data = 0

                data[tx1].append(s_data)
                data[tx2].append(s_data)

        for tx, vals in data.items():
            data[tx] = np.mean(vals)

        # Prepara data for plotting
        data_points = []
        data_labels = []
        for tx, vals in data.items():
            data_points.append(vals)
            data_labels.append(tx)

        data_points = np.asarray(data_points)
        data_labels = np.asarray(data_labels)

        # Get outliers
        outliers_points = data_points[self._mad_based_outlier(data_points)]
        # Get outlier taxa
        outlier_labels = list(data_labels[self._mad_based_outlier(data_points)])

        self._get_similarity("disconnect")

        return {"data": data_points,
                "title": "Sequence variation outlier taxa detection",
                "outliers": outliers_points,
                "outliers_labels": outlier_labels,
                "ax_names": ["Proportion of segregating sites", "Frequency"]}

    @check_data
    def outlier_sequence_size(self, ns=None):
        """Creates data for sequence size alignment outliers.

        Generates data for the outlier detection of genes based on their
        sequence length (excluding missing data).

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "outliers": list of outlier data points,

            "outliers_labels": list of outlier labels
        """

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        data_labels = []
        data_points = []

        prev_idx = ""
        for taxon, seq, aln_idx in self.iter_alignments():

            if prev_idx != aln_idx:
                
                if prev_idx:
                    gn_avg = np.mean(gn_l)
                    data_points.append(gn_avg)
                    data_labels.append(aln.name)

                self._update_pipes(ns, None, value=c)
                c += 1
                
                gn_l = []
                aln = self.alignment_idx[aln_idx]

                prev_idx = aln_idx

            gn_l.append(len(seq.replace(aln.sequence_code[1], "").
                            replace(self.gap_symbol, "")))

        if prev_idx:
            gn_avg = np.mean(gn_l)
            data_points.append(gn_avg)
            data_labels.append(aln.name)

        data_points = np.asarray(data_points)
        data_labels = np.asarray(data_labels)

        # Get outliers
        outliers_points = data_points[self._mad_based_outlier(data_points)]
        # Get outliers labels
        outliers_labels = list(data_labels[self._mad_based_outlier(
            data_points)])

        return {"data": data_points,
                "title": "Sequence size outlier gene detection",
                "outliers": outliers_points,
                "outliers_labels": outliers_labels,
                "ax_names": ["Sequence size", "Frequency"]}

    @check_data
    def outlier_sequence_size_sp(self, ns=None):
        """Create data for sequence size taxa outliers.

        Generates data for the outlier detection of species based on their
        sequence length (excluding missing data)

        Parameters
        ----------
        ns : multiprocesssing.Manager.Namespace
            A Namespace object used to communicate with the main thread
            in TriFusion.

        Returns
        -------
        _ : dict
            "data": `numpy.array` with data for plotting,

            "ax_names": 2 element list with axis labels [x, y],

            "title": str with title,

            "outliers": list of outlier data points,

            "outliers_labels": list of outlier labels
        """

        self._set_pipes(ns, None, total=len(self.alignments))
        c = 1

        data = dict((tx, []) for tx in self.taxa_names)

        prev_idx = ""
        for taxon, seq, aln_idx in self.iter_alignments():

            if prev_idx != aln_idx:

                self._update_pipes(ns, None, value=c)
                c += 1

                aln = self.alignment_idx[aln_idx]

                prev_idx = aln_idx

            s_data = len(seq.replace(self.sequence_code[1], "").
                         replace(self.gap_symbol, ""))
            data[taxon].append(s_data)

        # Get average for each taxon
        for tx, vals in data.items():
            data[tx] = np.mean(vals)

        # Preparing data for plotting
        data_points = []
        data_labels = []
        for tx, vals in data.items():
            data_points.append(vals)
            data_labels.append(tx)

        data_points = np.asarray(data_points)
        data_labels = np.asarray(data_labels)

        # Get outliers
        outliers_points = data_points[self._mad_based_outlier(data_points)]
        # Get outlier taxa
        outlier_labels = list(data_labels[self._mad_based_outlier(data_points)])

        return {"data": data_points,
                "title": "Sequence size outlier taxa detection",
                "outliers": outliers_points,
                "outliers_labels": outlier_labels,
                "ax_names": ["Sequence size", "Frequency"]}

__author__ = "Diogo N. Silva"
__credits__ = ["Diogo N. Silva", "Tiago F. Jesus"]