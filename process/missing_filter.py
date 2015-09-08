#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  
#  Copyright 2012 Unknown <diogo@arch>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  Author: Diogo N. Silva
#  Version: 0.1
#  Last update: 11/02/14

## TODO: Replace the progression prints with calls of the progression class

## TODO: __init__ is uncalled for. Replaced with another method and use it
## in the children objects. This is used for inheritance

from collections import OrderedDict
from itertools import compress

class MissingFilter ():
    """ Contains several methods used to trim and filter missing data from
    alignments. It's mainly used for inheritance """

    def __init__(self, alignment_dict, gap_threshold=50, missing_threshold=75,
                 gap_symbol="-", missing_symbol="n"):
        """ the gap_threshold variable is a cut-off to total_missing_proportion
         and missing_threshold in a cut-off to missing_proportion """

        self.alignment = alignment_dict
        self.gap = gap_symbol
        self.missing = missing_symbol

        # Defining special variables
        self.old_locus_length = None
        self.locus_length = None

        # Defining thresholds
        self.gap_threshold = gap_threshold
        self.missing_threshold = missing_threshold

        # Basic filter
        self.filter_terminals()
        self.filter_columns()

    def filter_terminals(self):
        """ Given an alignment, this will replace the gaps in the extremities
         of the alignment with missing data """

        for taxa, seq in self.alignment.items():

            if not seq.strip(self.gap):
                seq = self.missing * len(seq)
                self.alignment[taxa] = seq
                continue

            trim_seq = list(seq)
            counter, reverse_counter = 0, -1

            while trim_seq[counter] == self.gap:
                trim_seq[counter] = self.missing
                counter += 1

            while trim_seq[reverse_counter] == self.gap:
                trim_seq[reverse_counter] = self.missing
                reverse_counter -= 1

            seq = "".join(trim_seq)

            self.alignment[taxa] = seq

    def filter_columns(self):
        """ Here several missing data metrics are calculated, and based on
         some user defined thresholds, columns with inappropriate missing
         data are removed """

        taxa_number = float(len(self.alignment))
        self.old_locus_length = len(list(self.alignment.values())[0])

        filtered_cols = []

        # Creating the column list variable
        for column in zip(*self.alignment.values()):

            cadd = column.count

            # Calculating metrics
            gap_proportion = (float(cadd(self.gap)) /
                              taxa_number) * float(100)
            missing_proportion = (float(cadd(self.missing)) /
                                  taxa_number) * float(100)
            total_missing_proportion = gap_proportion + missing_proportion

            if total_missing_proportion <= self.gap_threshold or \
                    missing_proportion <= self.missing_threshold:

                filtered_cols.append(1)

            else:
                filtered_cols.append(0)

        self.alignment = OrderedDict((taxa,
                            "".join(list(compress(seq, filtered_cols))))
                            for taxa, seq in self.alignment.items())
        self.locus_length = len(list(self.alignment.values())[0])

__author__ = "Diogo N. Silva"
__copyright__ = "Diogo N. Silva"
__credits__ = ["Diogo N. Silva"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Diogo N. Silva"
__email__ = "o.diogosilva@gmail.com"
__status__ = "Prototype"