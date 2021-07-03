# An optimised quantum computer simulator

[![CMake](https://github.com/lanamineh/qsl/actions/workflows/cmake.yml/badge.svg?branch=master)](https://github.com/lanamineh/qsl/actions/workflows/cmake.yml)

This repository contains a quantum computer simulator which is optimised for a low to medium number of qubits (10-20, no more than about 28). The program will be able to simulate:

- a range of single- and two-qubit gates
- the outcome of measurements of the state vector

## What to do first

If you don't have access to the documentation from somewhere else, the first step is to build the documentation. For that, install doxygen and sphinx as follows:

```bash
 sudo apt install doxygen python3 python3-pip
 python3 -m pip install breathe sphinx_rtd_theme sphinx_copybutton
```

Then build the documentation by changing to the top level directory (containing this README.md) and running:

```bash
make docs
```

The documentation will then be located in `docs/html/index.html` (open it in a web browser).

## Build

The project uses cmake. To build the simulator, run the following commands in the top level directory (where this file is):

```bash
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/ ..
cmake --build .
```

Make sure to include the dots at the ends of the cmake commands. Replace ``/path/to/install/`` with any chosen installation destination. For example, if ``~/opt/`` is chosen, the shared library will be installed to ``~/opt/lib/`` and the header files will be installed to ``~/opt/include/``.

To build with clang instead of GCC, add the following options to the clang configuration: ``-DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++`` (making sure the paths to the clang compiler are correct). To build with the intel compiler, use the following configuration: ``-DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc``

After the build has finished, test executables will be in ``build/bin`` and the shared library ``libqsl.so`` will be in ``build/lib``

## Install

To install, after successfully building, change into the ``build/`` directory, and run

```bash
cmake --build . --target install
```

You might need to use sudo if the installation destination requires extra privileges. After the install has finished, check everything worked by making the example:

```bash
cd ../example
make
./example.bin
```

## Source code documentation

To read the (doxygen generated) documentation, clone the repository and run `make docs` inside the top level directory. The documentation will be created in a folder called html (navigate to `html/index.html`.

test change
