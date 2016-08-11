
"""
Stores long attributes for statistics plot data
"""

"""
stats_compliant contains the information required to setup the
Gene/Species/Average plot switcher in Statistics plots. Detail on this
object are in https://github.com/ODiogoSilva/TriFusion/wiki/Add-Statistics-plot-analysis#configure-triggers-between-plot-types

Briefly, args1, args2 and single_gene: Contain a dictionary with a single
element. The key is "plot_idx" and the value is the plot index string for
that plot type:
    .:args1: Species plot type
   .:args2: Average plot type
   .:single_gene: Single gene plot type
If a certain plot type does not exist, the value should be None
"""
stats_compliant = {
    "Distribution of sequence size":
    {"args1": None,
     "args2": {"plt_idx": "Distribution of sequence size all"},
     "active_bt": "sp",
     "single_gene": None},

    "Distribution of sequence size all":
    {"args1": {"plt_idx": "Distribution of sequence size"},
     "args2": None,
     "active_bt": "avg",
     "single_gene": None},

    "Proportion of nucleotides or residues":
    {"args1": {"plt_idx": "Proportion of nucleotides or residues"
                          " sp"},
     "args2": None,
     "active_bt": "avg",
     "single_gene": None},

    "Proportion of nucleotides or residues sp":
    {"args1": None,
     "args2": {"plt_idx": "Proportion of nucleotides or residues"},
     "active_bt": "sp",
     "single_gene": None},

    "Pairwise sequence similarity":
    {"args1": {"plt_idx": "Pairwise sequence similarity sp"},
    "args2": None,
     "active_bt": "avg",
     "single_gene": {"plt_idx": "Pairwise sequence similarity gn"}},

    "Pairwise sequence similarity sp":
    {"args1": None,
     "args2": {"plt_idx": "Pairwise sequence similarity"},
     "active_bt": "sp",
     "single_gene": {"plt_idx": "Pairwise sequence similarity gn"}},

    "Pairwise sequence similarity gn":
    {"args1": {"plt_idx": "Pairwise sequence similarity sp"},
     "args2": {"plt_idx": "Pairwise sequence similarity"},
     "active_bt": "gene",
     "single_gene": {"plt_idx": "Pairwise sequence similarity gn"}},

    "Allele Frequency Spectrum":
    {"args1": None,
     "args2": None,
     "active_bt": "avg",
     "single_gene": {"plt_idx": "Allele Frequency Spectrum gn"}},

    "Allele Frequency Spectrum prop":
        {"args1": None,
         "args2": None,
         "active_bt": "avg",
         "single_gene": {"plt_idx": "Allele Frequency Spectrum gn"}},

    "Allele Frequency Spectrum gn":
    {"args1": None,
     "args2": {"plt_idx": "Allele Frequency Spectrum"},
     "active_bt": "gene",
     "single_gene": {"plt_idx": "Allele Frequency Spectrum gn"}},

    "Distribution of missing orthologs":
    {"args1": None,
     "args2": {"plt_idx": "Distribution of missing orthologs avg"},
     "active_bt": "sp",
     "single_gene": None},

    "Distribution of missing orthologs avg":
    {"args1": {"plt_idx": "Distribution of missing orthologs"},
     "args2": None,
     "active_bt": "avg",
     "single_gene": None},

    "Distribution of missing data":
    {"args1": {"plt_idx": "Distribution of missing data sp"},
     "args2": None,
     "active_bt": "avg",
     "single_gene": None},

    "Distribution of missing data sp":
    {"args1": None,
     "args2": {"plt_idx": "Distribution of missing data"},
     "active_bt": "sp",
     "single_gene": None},

    "Segregating sites":
    {"args1": {"plt_idx": "Segregating sites sp"},
     "args2": None,
     "active_bt": "avg",
     "single_gene": {"plt_idx": "Segregating sites gn"}},

    "Segregating sites prop":
    {"args1": {"plt_idx": "Segregating sites sp"},
     "args2": None,
     "active_bt": "avg",
     "single_gene": {"plt_idx": "Segregating sites gn"}},

    "Segregating sites sp":
    {"args1": None,
     "args2": {"plt_idx": "Segregating sites"},
     "active_bt": "sp",
     "single_gene": {"plt_idx": "Segregating sites gn"}},

    "Segregating sites gn":
    {"args1": {"plt_idx": "Segregating sites sp"},
     "args2": {"plt_idx": "Segregating sites"},
     "active_bt": "gene",
     "single_gene": {"plt_idx": "Segregating sites gn"}},

    "Missing data outliers":
    {"args1": {"plt_idx": "Missing data outliers sp"},
     "args2": None,
     "active_bt": "avg",
     "single_gene": None},

    "Missing data outliers sp":
    {"args1": None,
     "args2": {"plt_idx": "Missing data outliers"},
     "active_bt": "sp",
     "single_gene": None},

    "Segregating sites outliers":
    {"args1": {"plt_idx": "Segregating sites outliers sp"},
     "args2": None,
     "active_bt": "avg",
     "single_gene": None},

    "Segregating sites outliers sp":
    {"args1": None,
     "args2": {"plt_idx": "Segregating sites outliers"},
     "active_bt": "sp",
     "single_gene": None},

    "Sequence size outliers sp":
    {"args1": None,
     "args2": {"plt_idx": "Sequence size outliers"},
     "active_bt": "sp",
     "single_gene": None},

    "Sequence size outliers":
    {"args1": {"plt_idx": "Sequence size outliers sp"},
     "args2": None,
     "active_bt": "avg",
     "single_gene": None}
}



__author__ = "Diogo N. Silva"
