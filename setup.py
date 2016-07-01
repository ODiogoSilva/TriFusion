import ez_setup
ez_setup.use_setuptools()

from setuptools import setup

with open('README.md') as f:
    readme = f.read()

setup(
    name="TriFusion",
    version="0.3.17",
    packages=["trifusion",
              "trifusion.base",
              "trifusion.data",
              "trifusion.data.backgrounds",
              "trifusion.data.resources",
              "trifusion.data.screens",
              "trifusion.ortho",
              "trifusion.process",
              "trifusion.stats"],
    package_data={"trifusion": ["*.kv"],
                  "trifusion.data.screens": ["*.kv"],
                  "trifusion.data.backgrounds": ["*.png", "*.ico"]},
    description=("Streamlining phylogenomic data gathering, processing and "
                 "visualization"),
    long_description=readme,
    url="https://github.com/ODiogoSilva/TriFusion",
    author="Diogo N Silva",
    author_email="odiogosilva@gmail.com",
    license="GPL3",
    classifiers=["Development Status :: 3 - Alpha",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: GNU General Public License v3 ("
                 "GPLv3)",
                 "Natural Language :: English",
                 "Operating System:: POSIX:: Linux",
                 "Operating System :: Microsoft :: Windows",
                 "Programming Language:: Python:: 2:: Only",
                 "Topic :: Scientific/Engineering :: Bio-Informatics"],
    entry_points={
        "gui_scripts": [
            "TriFusion = trifusion.TriFusion:main"
        ]
    },
)