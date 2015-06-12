__author__ = 'diogo'

from process.sequence import AlignmentList
import logging

def load_proc(aln_list, file_list, nm):
    try:
        if aln_list:
            aln_list.add_alignment_files(file_list, shared_namespace=nm)
            aln_obj = aln_list
        else:
            aln_obj = AlignmentList(file_list, shared_namespace=nm)
        nm.alns = aln_obj
    except:
        logging.exception("Unexpected error when loading input data")
        nm.exception = True
