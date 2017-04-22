## TriFusion

##### Making life easier for phylogenomic data gathering, processing and visualization

Website: http://odiogosilva.github.io/TriFusion/

<img src="https://raw.githubusercontent.com/ODiogoSilva/TriFusion-tutorials/master/tutorials/images/trifusion_home_screen.png"/>

[![Build Status](https://travis-ci.org/ODiogoSilva/TriFusion.svg?branch=master)](https://travis-ci.org/ODiogoSilva/TriFusion)
[![Code Health](https://landscape.io/github/ODiogoSilva/TriFusion/master/landscape.svg?style=flat)](https://landscape.io/github/ODiogoSilva/TriFusion/master)
[![codecov](https://codecov.io/gh/ODiogoSilva/TriFusion/branch/master/graph/badge.svg)](https://codecov.io/gh/ODiogoSilva/TriFusion)
[![PyPI](https://img.shields.io/pypi/pyversions/trifusion.svg)](https://pypi.python.org/pypi/trifusion)
[![PyPI](https://img.shields.io/pypi/v/trifusion.svg)](https://pypi.python.org/pypi/trifusion)
[![AUR](https://img.shields.io/aur/version/trifusion.svg)](https://aur.archlinux.org/packages/trifusion/)

[comment]: <> (<img align="right" height="128" src="https://github.com/ODiogoSilva/TriFusion/blob/43a41005ee8b1f69d7ae04684b0a0e595c527b4f/trifusion/data/backgrounds/trifusion-icon-256.png?raw=true"/>)

## What is TriFusion?

TriFusion is a modern GUI and command line application designed to make the life of anyone with **proteome** and/or **alignment sequence data** easier and more pleasurable. Regardless of your experience in bioinformatics, TriFusion is easy to use and offers a wide array of powerfull features to help you deal with your data. At the same time, it was developed to handle the enormous amount of data that is generated nowadays.

TriFusion is an open source, cross-platform application written in [Python 2.7](https://www.python.org/) and using the [Kivy](https://github.com/kivy/kivy) framework to construct graphical interface.

## What can TriFusion do for you?

Here is an overview of what it can do for you across its three main modules.

### Orthology - Search and explore orthologs across proteomes

 - **Searches for ortholog sequences** across multiple species.
 - Filters ortholog sequences according to the **gene copy number** and/or **number of taxa** present.
 - **Graphical visualization** of ortholog data.
 - Exports your orthologs as **protein or nucleotide sequences**.

 [Find out more]()

### Process - Blazing fast processing of alignment files

 - **Conversion** or **concatenation** of alignment files into several popular formats ([check supported formats]()).
 - **Collapse** identical sequences into the same haplotype.
 - Create **consensus** sequences for each alignment with several options on how to handle sequence variation.
 - **Filter** either alignments (according to whether they contain or exclude certain taxa, to a minimum proportion of taxa, and/or variable sites) or alignment columns (according to codon position, missing data and gaps).
 - **Code indel patterns** of your alignments into a binary matrix that is appended to the alignment.
 - **Revert concatenated alignments** or export sub-regions into individual files
 - Set **gene and codon partitions** as well as **substitution models** (Nexus format)
 - It's **fast** and **memory efficient**. Converting 3,093 files with 376 taxa can be performed in just 30 seconds and using less than 90Mb of RAM memory ([check the benchmarks table](https://github.com/ODiogoSilva/TriFusion/wiki/Benchmarks)).

[Find out more]()

### Statistics - Effortless visual exploration of your data

- Provides instant information on overall and per gene **summary statistcs**.
- TriFusion offers **dozens of graphical and statistical options** to explore your data:
     - General information plots.
     - Polymorphism and sequence variation plots.
     - Missing data plots.
     - Outlier plots.
- Publication ready figures

[Find out more]()

___

### Installation

#### Executables binaries

The latest stable release of TriFusion can be installed as a standalone application using one of the following installers. This only includes the GUI component of TriFusion.

##### Linux

- Debian package based ([See list](https://en.wikipedia.org/wiki/Category:Debian-based_distributions)): [TriFusion-0.5.0.deb](https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0.deb)

- RPM package based ([See list](https://en.wikipedia.org/wiki/Category:RPM-based_Linux_distributions)): [TriFusion-0.5.0.rpm](https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0.rpm)

- ArchLinux/Manjaro ([See list](https://wiki.archlinux.org/index.php/Arch_based_distributions)): [TriFusion-0.5.0.tar.xz](https://aur.archlinux.org/packages/trifusion-bin/) is available on AUR.

##### MacOS

- [TriFusion-0.5.0.app.zip](https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0-MacOS.app.zip)

##### Windows

- [TriFusion-0.5.0 64bit installer](https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0-Win64.msi)
- [TriFusion-0.5.0 32bit installer](https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0-Win32.msi)

###### Note for Windows 8.x and 10 users:

Executing the TriFusion installer may generate a warning from SmartScreen. To continue with the installation, click the "More info" label and then "Run anyway".

#### Installation from source

TriFusion can be installed directly from source (whether from the latest stable release or from the git version), which has its advantages:

- Stay on the bleeding edge of the application's development (Git version)

- Minimize the application size install

- Modify the code any way you want. Contributions are welcome.


[Here are instructions on how to install TriFusion from source](https://github.com/ODiogoSilva/TriFusion/wiki/Install-from-source). (Go ahead, it won't bite.)

### How to use

#### GUI version

##### Windows

TriFusion should be available in the Start Menu, ready to launch.

##### Linux

If TriFusion was installed using package managers, a TriFusion shortcut should have been created and available using your distribution's application launcher (Under the Science category). Alternatively, TriFusion can be executed from the terminal with the command `TriFusion`.

#### Command line versions

Command line versions for each of the three modules in TriFusion are only available when installing the program from source. Once installed, each module can be executed in the terminal using the following commands:

- Orthology search module:

    `orthomcl_pipeline`

- Process module:

    `TriSeq`

- Statistics module:

    `TriStats`

For more information on the command line versions, see the manual.

### Documentation

You can download TriFusion User Guide [here](https://github.com/ODiogoSilva/TriFusion/raw/master/docs/manual.pdf).

### Citation

When using OrthoMCL to find ortholog clusters, please cite the original software:

Fischer, S., Brunk, B. P., Chen, F., Gao, X., Harb, O. S., Iodice, J. B., Shanmugam, D., Roos, D. S. and Stoeckert, C. J. Using OrthoMCL to Assign Proteins to OrthoMCL-DB Groups or to Cluster Proteomes Into New Ortholog Groups Current Protocols in Bioinformatics. 2011 35:6.12.1–6.12.19.

Coming soon!
