## TriFusion
###### Streamlining phylogenomic data gathering, processing and visualization

TriFusion is a GUI and command line application designed to streamline the workflow of phylogenomic projects. With the dramatic increase in size of data sets for phylogenetics and population genetics, programming has become a crucial tool to gather, process and analyze the data. However, this may still represent a hurdle that precludes the execution of such projects by a broader range of researchers. TriFusion aims to mitigate this issue by providing a user-friendly visual interface that empowers users without any programming knowledge with a set of tools and operations that can handle large sets of data.

TriFusion is an open source, cross-platform application written in [Python 2.7](https://www.python.org/) and using the [Kivy](https://github.com/kivy/kivy) framework to construct graphical interface.

### What can TriFusion do for you

The application is roughly divided into three main modules: 

- **Orthology**: Detect and visualize of orthologs, co-orthologs and recent paralogs from complete genomes with just a few clicks
- **Process**: Convert or concatenate alignments of Protein/DNA/RNA into several output formats (fasta, phylip, nexus, mcmctree, IMa2) in addition to a number of secondary operations (collapsing, gap coding, filtering, etc.)
- **Statistics**: Visualize your data set using several plotting operations.

For a complete overview of the features for each module, please refer to the manual (Coming soon).

### Executables

[Linux executable (V.0.1.1)](https://github.com/ODiogoSilva/TriFusion/releases/download/v0.1.1/TriFusion-linux_current.zip)

[Windows installer (V.0.1.1)](https://github.com/ODiogoSilva/TriFusion/releases/download/v0.1.1/TriFusion-windows_current.msi)

##### Note for Windows 8.x and 10 users:

Executing the TriFusion installer may generate a warning from SmartScreen. To circumvent this issue, click the "More info" label and then "Run anyway".

### Installation from source

TriFusion is regularly updated with new features and bug fixes. These will eventually be added into the executable versions but if you wish to stay on the bleeding edge of the application's development (at the cost of potential creeping bugs), use the git version. The dependencies required to run TriFusion follow below.

#### Dependencies

Python 2.7 is required to run TriFusion. It should already be present in most Unix operating systems but not on Windows. For Windows users, I recommend installing a python distribution package (such as [Anaconda](https://www.continuum.io/downloads#_windows)) that already comes with several required python libraries. Alternatively, the vanilla Python installer can be downloaded [here](https://www.python.org/downloads/).

Assuming python is already installed, the following libraries are required:

- `kivy` 
- `matplotlib` (included in Anaconda)
- `numpy` (included in Anaconda)
- `psutil` (included in Anaconda)
- `scipy` (included in Anaconda)
- `configparser`
- `seaborn`

These can be easily installed using `pip`:

```
pip install kivy matplotlib numpy psutil scipy seaborn
```

##### Note for Windows users

Installation of Kivy on Windows machines, requires a few more commands that are explained in [kivy's webpage](https://kivy.org/docs/installation/installation-windows.html#installation). Shortly, Kivy's dependencies must be installed beforehand:

```
python -m pip install docutils pygments pypiwin32 kivy.deps.sdl2 kivy.deps.glew
```

### Know issues

When using the executable version of TriFusion make sure its path does not contain non-ASCII characters, otherwise the app will not open. This issue is currently being fixed.

### Documentation

User guide coming soon!

### Citation

When using OrthoMCL to find ortholog clusters, please cite the original software:

Fischer, S., Brunk, B. P., Chen, F., Gao, X., Harb, O. S., Iodice, J. B., Shanmugam, D., Roos, D. S. and Stoeckert, C. J. Using OrthoMCL to Assign Proteins to OrthoMCL-DB Groups or to Cluster Proteomes Into New Ortholog Groups Current Protocols in Bioinformatics. 2011 35:6.12.1–6.12.19.

Coming soon!
