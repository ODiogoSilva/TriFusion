Installation
============

Binaries and Installers
-----------------------

.. note::

    The installers **only provide the GUI version** of TriFusion. If
    you want to install both the GUI and command line versions, check the
    `Installation from source`_.

The easiest way to install TriFusion is through binaries and installers
provided for Windows, MacOS and Linux.

Windows
^^^^^^^

    - `TriFusion 1.0.0 64bit installer <https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0-Win64.msi>`_
    - `TriFusion 1.0.0 32bit installer <https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0-Win32.msi>`_

MacOS
^^^^^

    - `TriFusion-1.0.0.app <https://github.com/ODiogoSilva/TriFusion/releases/download/0.5.0/TriFusion-v0.5.0-MacOS.app.zip>`_

Linux
^^^^^

    - Ubuntu and Debian based: `TriFusion-1.0.0.deb <https://github.com/ODiogoSilva/TriFusion/releases/download/1.0.0/TriFusion-v1.0.0.deb>`_
    - RPM based: `TriFusion-1.0.0.rpm <https://github.com/ODiogoSilva/TriFusion/releases/download/1.0.0/TriFusion-v1.0.0.rpm>`_
    - General linux: `TriFusion-1.0.0.tar.xz <https://aur.archlinux.org/packages/trifusion-bin/>`_


Installation from source
------------------------

TriFusion is available on `PyPi <https://pypi.python.org/pypi/trifusion>`_
nd can be easily installed with ``pip`` . This
will only install the command line versions of TriFusion
(``TriSeq``, ``TriStats`` and ``orthomcl_pipeline``). Therefore, if you are
interested only in the command line version of TriFusion, and assuming you
have ``python2.7`` and ``pip`` on your system, installing TriFusion is simply:

.. code-block:: bash

    pip install trifusion

The dependencies for the graphical user interface require only a
few extra commands that are provided below for each operating system.

Windows
^^^^^^^

Windows does not come with  a python installation by default. We recommend
using a package manager, such as
`Anaconda <https://www.continuum.io/downloads>`_, which automatically
installs most of the dependencies (**Note that TriFusion requires
python2.7**). After installing python, you will need
to install ``kivy`` by executing the following commands on a command
line prompt:

.. code-block:: bash

    python -m pip install --upgrade pip wheel setuptools
    python -m pip install docutils pygments pypiwin32 kivy.deps.sdl2 kivy.deps.glew
    python -m pip install kivy.deps.gstreamer --extra-index-url https://kivy.org/downloads/packages/simple/
    python -m pip install kivy

Then, install trifusion by typing:

.. code-block:: bash

    pip install trifusion

MacOS (using homebrew)
^^^^^^^^^^^^^^^^^^^^^^

f you do not have `homebrew <https://brew.sh/>`_ yet, you'll need to install
it:

.. code-block:: bash

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Then, to install TriFusion and it's dependencies:

.. code-block:: bash

    brew install sdl2 sdl2_image sdl2_ttf sdl2_mixer
    pip install -I Cython==0.23
    USE_OSX_FRAMEWORKS=0 pip install kivy
    pip install trifusion

Ubuntu (and relatives)
^^^^^^^^^^^^^^^^^^^^^^

On Ubuntu, there are PPAs available for the installation of TriFusion via
``apt-get`` in addition to the ``pip`` installation method.

Via PPA
~~~~~~~

    - Add one of the following PPAs:
        ::

            # Stable release:
            sudo add-apt-repository ppa:o-diogosilva/trifusion
            # Daily release:
            sudo add-apt-repository ppa:o-diogosilva/trifusion-daily

    - Upgrade your package list and install TriFusion:
        ::

            sudo apt-get update && sudo apt-get install trifusion

Via ``pip``
~~~~~~~~~~~

.. code-block:: bash

    sudo apt-get install python-pip build-essential python-dev libsdl2-dev
    pip install cython==0.23
    pip install kivy
    pip install trifusion


Debian
^^^^^^

As with Ubuntu, you may install TriFusion via the available PPAs or with
``pip``.


Via PPA
~~~~~~~

    - Add one of the following PPAs manually to the ``sources.list`` file:
        ::

            # Stable release:
            http://ppa.launchpad.net/o-diogosilva/trifusion/ubuntu trusty main
            # Daily release:
            http://ppa.launchpad.net/o-diogosilva/trifusion-daily/ubuntu trusty main

    - Add the GPG key to your apt keyring:
        ::

             sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys D4F1E8E6


    - Upgrade your package list and install TriFusion:
        ::

            sudo apt-get update && sudo apt-get install trifusion

Via ``pip``
~~~~~~~~~~~

.. code-block:: bash

    sudo apt-get install python-pip build-essential python-dev libsdl2-dev
    pip install cython==0.23
    pip install kivy
    pip install trifusion

RPM based
^^^^^^^^^

.. code-block:: bash

    dnf install python-pip python-devel redhat-rpm-config freeglut-devel SDL* libsdl2-dev
    pip install cython==0.23
    pip install kivy
    pip install trifusion

ArchLinux
^^^^^^^^^

There are three `AUR <https://aur.archlinux.org/>`_ packages for TriFusion:

    - `trifusion <https://aur.archlinux.org/packages/trifusion/>`_:
      The latest release of TriFusion, based on source code.
    - `trifusion-bin <https://aur.archlinux.org/packages/trifusion-bin/>`_:
      The latest release of TriFusion. in binary format. Does
      not require dependencies to be installed, as all the necessary libs
      are bundled with the distributed binary
    - `trifusion-git <https://aur.archlinux.org/packages/trifusion-git/>`_: The bleeding edge version directly from git. Requires
      dependencies to be installed, as it is also source code based.

Just use any `AUR helper <https://wiki.archlinux.org/index.php/AUR_helpers>`_
to handle the packages for you, or download the *PKGBUILD* you require and
use ``makepkg``.

