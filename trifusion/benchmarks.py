#!/usr/bin/python

import os
import shutil
import cProfile
from memory_profiler import profile

from process.sequence import AlignmentList
from tests.data_files import dna_data_fas
# Benchmarks for comparing file I/O vs. sqlite database for process operations

benchmark_dir = "C:\Users\Diogo\Documents\Tests\output"
sql_path = "C:\Users\Diogo\Documents\Tests\output\sequence.db"

class Wrapper(object):

    def __init__(self, func):

        self.func = func

    def __call__(self, *args):

        if not os.path.exists(benchmark_dir):
            os.mkdir(benchmark_dir)

        self.func(*args)

        shutil.rmtree(benchmark_dir)
        print("Test finished")

################################# DATA ########################################

medium = [os.path.join("C:\Users\Diogo\Documents\TriFusion\Test_data\Process"
                                 "\medium_protein_dataset", x) for x in os.listdir("C:\Users\Diogo\Documents\TriFusion\Test_data\Process"
                                 "\medium_protein_dataset")]

basidio9sp = [os.path.join("C:\Users\Diogo\Documents\TriFusion\Test_data\Process\Basidio3k", x) for x in os.listdir("C:\Users\Diogo\Documents\TriFusion\Test_data\Process\Basidio3k")]

large = [os.path.join("C:\Users\Diogo\Documents\TriFusion\Test_data\Process\\7kloci", x) for x in os.listdir("C:\Users\Diogo\Documents\TriFusion\Test_data\Process\\7kloci")]

single_large = "C:\Users\Diogo\Documents\TriFusion\Test_data\Process\single_big_file\c97d5m4p2_MM80_interleave.phy"

single_large_interleave = "C:\Users\Diogo\Documents\TriFusion\Test_data\Process\single_big_file\c97d5m4p2_MM80_interleave.fas"

loci_file = "C:\Users\Diogo\Documents\GitHub\TriFusion\\trifusion\\tests\data\c97d5m4p2.loci"

nexus_file = "C:\Users\Diogo\Documents\GitHub\TriFusion\\trifusion\\tests\data\BaseConcatenation.nex"

stock_file = "C:\Users\Diogo\Documents\GitHub\TriFusion\\trifusion\\tests\data\BaseConc7.stockholm"
################################ Input ########################################

# Test 1: Single file
@Wrapper
@profile
def run_single():
    print("Running Input test: run_single")
    single_file = os.path.join("tests", "data", "BaseConc1.fas")
    AlignmentList([single_file], dest=benchmark_dir)

# Test 2: 7 files
@Wrapper
@profile
def seven_files():
    print("Running Input test: seven_files")
    AlignmentList(dna_data_fas, dest=benchmark_dir)

# Test 3: 614 files
@Wrapper
#@profile
def medium_dataset():
    print("Running Input test: medium_dataset")
    AlignmentList(medium, dest=benchmark_dir, sql_db=sql_path)

# Test 4: 3k files
@Wrapper
#@profile
def basidio9sp_dataset():
    print("Running Input test: basidio9sp_dataset")
    AlignmentList(basidio9sp, dest=benchmark_dir, sql_db=sql_path)

# Test 5: 7k files
@Wrapper
@profile
def large_dataset():
    print("Running Input test: large_dataset")
    AlignmentList(large, dest=benchmark_dir)

# Test 6: One single file
@Wrapper
# @profile
def large_file():
    print("Running Input test: large_file")
    AlignmentList([single_large], dest=benchmark_dir, sql_db=sql_path)

@Wrapper
# @profile
def large_interleave_file():
    print("Running Input test: large_interleave_file")
    AlignmentList([single_large_interleave], dest=benchmark_dir)

# Test 8: Loci file
@Wrapper
#@profile
def loci_dataset():
    print("Running Input test: loci_dataset")
    AlignmentList([loci_file], sql_db=sql_path)

# Test 9: Nexus file
@Wrapper
def nexus_dataset():
    print("Running Input test: nexus_dataset")
    x = AlignmentList([nexus_file], sql_db=sql_path)
    # x.remove_taxa(['999_RAD_original', '3305_RAD_original',
    #                '3536_RAD_original', 'spa', 'spb', 'spc', 'spd'])
    x.change_taxon_name("spa", "BADONG")
    print(x.taxa_names)

# Test10: Stockholm file
@Wrapper
def stock_dataset():
    print("Running Input test: stock_dataset")
    AlignmentList([stock_file], sql_db=sql_path)

################################ Output #######################################

# if not os.path.exists(benchmark_dir):
#     os.mkdir(benchmark_dir)
#
# aln_obj = AlignmentList(medium, dest=benchmark_dir)

#@profile
def concat(aln_obj):
    print("Running Output test: concat")
    x = aln_obj.concatenate(alignment_name="test", dest=benchmark_dir,
                            remove_temp=True)

#@profile
def filter_alns(aln_obj):
    print("Running Output test: filter_alns")
    x = aln_obj.filter_missing_data(50, 25)

@profile
def consensus_alns(aln_obj):
    print("Running Output test: consensus_alns")
    x = aln_obj.consensus(consensus_type="Soft mask", single_file=True)

@profile
def collapse_alns(aln_obj):
    print("Running Output test: collapse_alns")
    aln_obj.collapse(dest=benchmark_dir)

################################ Run ##########################################

#cProfile.run("run_single()")
#run_single()

#cProfile.run("seven_files()")
#seven_files()

# cProfile.run("medium_dataset()")
# medium_dataset()

# cProfile.run("basidio9sp_dataset()")
#basidio9sp_dataset()

#cProfile.run("large_dataset()")
#large_dataset()

# cProfile.run("large_file()")
# large_file()

#cProfile.run("large_interleave_file()")
#large_interleave_file()

# cProfile.run("loci_dataset()")
# loci_dataset()

cProfile.run("nexus_dataset()")
# nexus_dataset()

# cProfile.run("stock_dataset()")
# stock_dataset()

#cProfile.run("concat(aln_obj)")
#concat(aln_obj)

#cProfile.run("filter_alns(aln_obj)")
#filter_alns(aln_obj)

#cProfile.run("consensus_alns(aln_obj)")
#consensus_alns(aln_obj)

#cProfile.run("collapse_alns(aln_obj)")
# collapse_alns(aln_obj)

# shutil.rmtree(benchmark_dir)