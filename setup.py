try:
    from setuptools import setup
except ImportError:
    import ez_setup
    ez_setup.use_setuptools()
    from setuptools import setup

import trifusion

VERSION = trifusion.__version__

with open('README.rst') as f:
    readme = f.read()


def mcl_data_files():
    """
    Automatically detects the platform and sets the appropriate mcl binary
    """

    data_file = ["trifusion/data/resources/mcl/linux/mcl",
                 "trifusion/data/resources/mcl/MacOS/mcl",
                 "trifusion/data/resources/mcl/windows/64bit/mcl64.exe",
                 "trifusion/data/resources/mcl/windows/64bit/cygwin1.dll",
                 "trifusion/data/resources/mcl/windows/32bit/mcl32.exe",
                 "trifusion/data/resources/mcl/windows/32bit/cygwin1.dll"]

    return data_file


mcl_file = mcl_data_files()

setup(
    name="trifusion",
    version=VERSION,
    packages=["trifusion",
              "trifusion.base",
              "trifusion.data",
              "trifusion.data.backgrounds",
              "trifusion.data.resources",
              "trifusion.data.resources.theme",
              "trifusion.data.screens",
              "trifusion.ortho",
              "trifusion.process",
              "trifusion.progressbar"],
    package_data={"trifusion": ["*.kv"],
                  "trifusion.data.screens": ["*.kv"],
                  "trifusion.data.backgrounds": ["*.png", "*.ico"],
                  "trifusion.data.resources": [".desktop"]},
    install_requires=[
        "seaborn",
        "configparser",
        "pandas",
        "matplotlib",
        "numpy",
        "psutil",
        "scipy",
    ],
    description=("Streamlining phylogenomic data gathering, processing and "
                 "visualization"),
    long_description=readme,
    url="https://github.com/ODiogoSilva/TriFusion",
    author="Diogo N Silva",
    author_email="odiogosilva@gmail.com",
    license="GPL3",
    classifiers=["Development Status :: 4 - Beta",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: GNU General Public License v3 ("
                 "GPLv3)",
                 "Natural Language :: English",
                 "Operating System :: POSIX :: Linux",
                 "Operating System :: MacOS :: MacOS X",
                 "Operating System :: Microsoft :: Windows",
                 "Programming Language :: Python",
                 "Programming Language :: Python :: 2.7",
                 "Topic :: Scientific/Engineering :: Bio-Informatics"],
    scripts=mcl_file,
    entry_points={
        "gui_scripts": [
            "TriFusion = trifusion.TriFusion:gui_exec"
        ],
        "console_scripts": [
            "orthomcl_pipeline = trifusion.orthomcl_pipeline:main",
            "TriSeq = trifusion.TriSeq:main",
            "TriStats = trifusion.TriStats:main"
        ]
    },
)