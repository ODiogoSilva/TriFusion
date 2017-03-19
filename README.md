## Warning

TriFusion is still in active development, so the releases are quite a bit out of date 
and the install scripts and steps may not be properly working. If you have any issues,
please contact me directly.

## TriFusion

Website: http://odiogosilva.github.io/TriFusion/

[![Build Status](https://travis-ci.org/ODiogoSilva/TriFusion.svg?branch=master)](https://travis-ci.org/ODiogoSilva/TriFusion)
[![Code Health](https://landscape.io/github/ODiogoSilva/TriFusion/master/landscape.svg?style=flat)](https://landscape.io/github/ODiogoSilva/TriFusion/master)
[![codecov](https://codecov.io/gh/ODiogoSilva/TriFusion/branch/master/graph/badge.svg)](https://codecov.io/gh/ODiogoSilva/TriFusion)


###### Streamlining phylogenomic data gathering, processing and visualization

<img align="right" height="256" src="https://github.com/ODiogoSilva/TriFusion/blob/43a41005ee8b1f69d7ae04684b0a0e595c527b4f/trifusion/data/backgrounds/trifusion-icon-256.png?raw=true"/>


TriFusion is a GUI and command line application designed to streamline the workflow of phylogenomic projects. With the dramatic increase in size of data sets for phylogenetics and population genetics, programming has become a crucial tool to gather, process and analyze the data. However, this may still represent a hurdle that precludes the execution of such projects by a broader range of researchers. TriFusion aims to mitigate this issue by providing a user-friendly visual interface that empowers users without any programming knowledge with a set of tools and operations that can handle large sets of data.

TriFusion is an open source, cross-platform application written in [Python 2.7](https://www.python.org/) and using the [Kivy](https://github.com/kivy/kivy) framework to construct graphical interface.

### Installation

#### Executables binaries

The latest stable release of TriFusion can be installed as a standalone application using one of the following installers.

##### Linux

- Debian package based ([See list](https://en.wikipedia.org/wiki/Category:Debian-based_distributions)): [trifusion-bin-0.4.11-linux.deb](https://github.com/ODiogoSilva/TriFusion/releases/download/0.4.11/trifusion-bin-0.4.11-linux.deb)

- RPM package based ([See list](https://en.wikipedia.org/wiki/Category:RPM-based_Linux_distributions)): [trifusion-bin-0.4.11-linux.rpm](https://github.com/ODiogoSilva/TriFusion/releases/download/0.4.11/trifusion-bin-0.4.11-linux.rpm)

- ArchLinux/Manjaro ([See list](https://wiki.archlinux.org/index.php/Arch_based_distributions)): [trifusion-bin-0.4.11](https://aur.archlinux.org/packages/trifusion-bin/) is available on AUR.

##### Windows

- [TriFusion-0.4.12 64bit installer](https://github.com/ODiogoSilva/TriFusion/releases/download/0.4.11/TriFusion-v0.4.12-windows64.msi)
- [TriFusion-0.4.12 32bit installer](https://github.com/ODiogoSilva/TriFusion/releases/download/0.4.11/TriFusion-v0.4.12-windows32.msi)

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
