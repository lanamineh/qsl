Getting the source code
#######################

Clone the git repository by running

.. code-block:: console

   $ git clone git@githome.gq:lanamineh/qsl.git

Documentation
#############

To build this documentation, install python3 and pip. Then run:

.. code-block:: console
		
   $ sudo apt install doxygen python3 python3-pip
   $ python3 -m pip install breathe sphinx_rtd_theme sphinx_copybutton

After installing the prerequisites, run

.. code-block:: console
		
   $ make docs

from the top level directory in the repository. Use a browser to navigate to docs/html/index.html (relative to the root of the repository).

Prerequisites for Building
##########################

There are a few prerequisites for building the core simulator.

CMake (at least version 3.12.4)
*******************************

The project uses the CMake buildsystem. A relatively recent version of CMake is required (at least 3.12.4), which might not be available in the ubuntu repository. The CMake website contains instructions for installing the latest release on the `downloads <https://cmake.org/download/>`_ page. For example, to install version 3.20.0, run the commands below. They will install cmake into a custom folder (``$HOME/opt``), so that you can easily get rid of it later if you want. You can change the installation location to anything using the ``--prefix`` option (you might need sudo rights depending on what you choose).

.. code-block:: console
		
   $ wget https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0-linux-x86_64.sh
   $ chmod u+x cmake-3.20.0-linux-x86_64.sh
   $ ./cmake-3.20.0-linux-x86_64.sh --help # If you want to see your options
   $ mkdir $HOME/opt # Make a custon installation directory
   $ ./cmake-3.20.0-linux-x86_64.sh --prefix=$HOME/opt/ --exclude-subdir

If you choose to put the installation in ``$HOME/opt/``, remember to add ``$HOME/opt/bin`` to your path (add ``export PATH=$PATH:$HOME/opt/bin/`` to your ``.bashrc``).

To check everything worked, run

.. code-block:: console

   $ cmake --version # Should output 3.20.0

If you get a different version, or any errors, you might have installed cmake previously using the package manager. You can uninstall it using

.. code-block:: console

   $ sudo apt remove cmake

This should remove the system installation of cmake, so that your new cmake is used instead.

GCC (at least version 10)
*************************

A recent version of gcc is required (at least gcc 10). At the moment this in not in the ubuntu package manager. You can install it by adding a ppa as follows:

.. code-block:: console

   $ sudo add-apt-repository ppa:ubuntu-toolchain-r/test
   $ sudo apt update
   $ sudo apt install gcc-10
   $ sudo apt install g++-10
   
If you have multiple versions of gcc on your system, and you want to be able to swap between them, a convenient way to do it is using ``update-alternatives``. For example, assuming you also have gcc-9 installed to ``/usr/bin/gcc-9``, you can run:

.. code-block:: console

   $ sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 900 --slave /usr/bin/g++ g++ /usr/bin/g++-9 
   $ sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 1000 --slave /usr/bin/g++ g++ /usr/bin/g++-10

The order of arguments to ``update-alternatives --install`` is:

#. Path to the link exeutable (e.g. /usr/bin/gcc).
#. Name of link group.
#. Path to the true executable being linked (use ``whereis g++-9`` to find it).
#. A priority, which determines the default choice (higher number means higher priority).

The ``--slave`` arguments (which are the same as for ``--install``) make the current g++ version depend on the current gcc version. 

Once you have set up the link group, you can swap between different alternatives using the ``--config`` option:

.. code-block:: console

   $ sudo update-alternatives --config gcc
   There are 2 choices for the alternative gcc (providing /usr/bin/gcc).

   Selection    Path             Priority   Status
   ------------------------------------------------------------
   * 0            /usr/bin/gcc-10   1000      auto mode
     1            /usr/bin/gcc-10   1000      manual mode
     2            /usr/bin/gcc-9    900       manual mode

   Press <enter> to keep the current choice[*], or type selection number:

You can enter a number and press enter to select the compiler version. After running the command, check the compiler version for gcc and g++ using:

.. code-block:: console

   $ gcc --version
   $ g++ --version

The command ``update-alternatives`` can be used any time you want to choose between different versions of a program. See ``man update-alternatives`` for more information.

Building the Core Simulator
###########################

To build the simulator, clone the repository and install the prerequisites as above, and then run the following commands from the top level directory:

.. code-block:: console

   $ mkdir build # Make a folder for the cmake build
   $ cd build/
   $ cmake ..
   $ cmake --build .

Executable files will now be in ``build/bin/`` and the library ``libqsl.so`` will be in ``build/lib/``.
