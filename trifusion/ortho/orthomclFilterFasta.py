#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os
import re

try:
    from process.error_handling import KillByUser
except ImportError:
    from trifusion.process.error_handling import KillByUser


def orthomcl_filter_fasta(input_dir, min_length, max_stop_percent, db, dest,
                             nm=None):

    def handle_seq(seq, length, stop_cnt):
        is_bad = 0
        stop_percent = ((length - stop_cnt) / length) * 100

        if length < min_length or stop_percent > max_stop_percent:
            bad.write(seq + "\n")
            is_bad = 1
        else:
            good.write(seq + "\n")

        return is_bad

    good = open(os.path.join(dest, "backstage_files", db), "w")
    bad = open(os.path.join(dest, "backstage_files", "poorProteins.txt"), "w")

    filenames = [os.path.join(input_dir, x) for x in os.listdir(input_dir)]

    reject_rates = []

    # Setup progression information
    if nm:
        if nm.stop:
            raise KillByUser("")
        nm.total = len(filenames)
        nm.counter = 0

    for filename in filenames:

        if nm:
            if nm.stop:
                raise KillByUser("")
            nm.counter += 1
            nm.msg = "Filtering file {}".format(os.path.basename(filename))

        if filename.startswith('.'):
            continue

        input_file = open(filename, 'r')
        seq_count = 0
        reject_seq_count = 0
        current_seq = ""
        current_len = 0
        current_stop_cnt = 0

        # process lines of one file
        for line in input_file:

            if nm:
                if nm.stop:
                    raise KillByUser("")

            if line.startswith('>'):
                if current_seq:
                    seq_count += 1
                    reject_seq_count += handle_seq(current_seq,
                                                   current_len,
                                                   current_stop_cnt)
                    current_seq = ""
                    current_len = 0
                    current_stop_cnt = 0
            else:
                line_len = len(line)
                current_len += line_len
                line = re.sub('[^A-Za-z]', '', line)
                current_stop_cnt += line_len - len(line)

            current_seq += line

        reject_seq_count += handle_seq(current_seq,
                                       current_len,
                                       current_stop_cnt)
        seq_count += 1

        # add file stats to reject count if it qualifies
        if reject_seq_count:
            pct = reject_seq_count / seq_count * 100
            if pct > 10:
                reject_rates.append([input_file, pct])

        input_file.close()

    good.close()
    bad.close()


__author__ = "Fernando Alves"

