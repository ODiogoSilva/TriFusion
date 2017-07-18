Limitations for input files
---------------------------

The **Process** module deals with several input formats and sequence types,
which begs the question of whether there are limitations on the type
of files that can be loaded simultaneously into TriFusion.

**The answer is almost none.**

TriFusion was designed to capture all the details
about your files automatically and to handle any combination you can
throw at it. In the example below, 8 alignment files of nucleotide and
protein sequences in Fasta, Nexus, Phylip and Stockholm formats are loaded
simultaneously. Then, these files are easily concatenated into a single
file with just a few clicks.

.. image:: https://github.com/ODiogoSilva/TriFusion-tutorials/raw/master/tutorials/gifs/process_mix_input.gif
    :alt: pic

Moreover, defining partitions when there are multiple files and sequence
types can be extremely time consuming and error prone to perform manually.
That is why TriFusion handles all of that automatically. Even though we
did not dealt with partitions in the above example, when you open a Nexus
alignment file, you can see that the header and partitions block are
correctly defined without any user intervention::

    #NEXUS
    Begin data;
        dimensions ntax=101 nchar=7030 ;
        format datatype=mixed(dna:1-3934,protein:3935-7030) interleave=no gap=-;

    (... DATA ...)

    begin mrbayes;
        charset DNAfas = 1-668;
        charset DNAnex = 669-1140;
        charset DNAphy = 1141-1808;
        charset DNAstockholm = 1809-2476;
        charset PROTEINphy = 2477-3934;
        charset PROTEINfasta = 3935-4966;
        charset PROTEINnex = 4967-5998;
        charset PROTEINstockholm = 5999-7030;
        partition part = 8: DNAfas, DNAnex, DNAphy, DNAstockholm, PROTEINphy, PROTEINfasta, PROTEINnex, PROTEINstockholm;
        set partition=part;
    end;

