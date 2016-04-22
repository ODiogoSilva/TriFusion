#!/usr/bin/python2


def mcl_to_groups(prefix, start_id, infile, outfile):

    try:
        start_id = int(start_id)
    except ValueError:
        raise ValueError("StartId is not a number")

    input_file = open(infile, "r")
    out = open(outfile, "w")

    for line in input_file:
        out.write(prefix + str(start_id) + ": " + line)
        start_id += 1
