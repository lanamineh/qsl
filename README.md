# A C++ Quantum Computer Simulation Library

[![build](https://github.com/lanamineh/qsl/actions/workflows/build.yml/badge.svg)](https://github.com/lanamineh/qsl/actions/workflows/build.yml) [![tests](https://github.com/lanamineh/qsl/actions/workflows/tests.yml/badge.svg)](https://github.com/lanamineh/qsl/actions/workflows/tests.yml) [![codecov](https://codecov.io/gh/lanamineh/qsl/branch/master/graph/badge.svg?token=VYUJ0OZIEZ)](https://codecov.io/gh/lanamineh/qsl) [![Documentation Status](https://readthedocs.org/projects/qsl/badge/?version=latest)](https://qsl.readthedocs.io/en/latest/?badge=latest) 

:warning: **The project is not stable yet. The documentation and library interface is likely to change in the future.** 

The repository contains a library for simulating quantum computers with a low to medium number of qubits. It is designed to be simple to use and fast. The library contains a few different types of simulator, depending on what you want to do. At the moment, the available simulators are as follows:

| Simulator | Description |
|-----------|-------------|
| Default   | an unrestricted quantum simulator without parallelisation.
| OpenMP    | an unrestricted quantum simulator using OpenMP parallelisation.
| Number Preserving | a simulator optimised for a restricted set of number-preserving gates.
| Resizeable | a simulator that can be resized (qubits can be added or removed). 

The following code snippet shows a simple example of the library usage:

```c++
#include <qsl/qubits.hpp>
int main()
{
	qsl::Qubits<qsl::Type::Default, double> q{ 5 }; // Make a 5 qubit simulator
	q.hadamard(2); // Apply a Hadamard gate to the second qubit
	q.controlNot(0,1); // Apply a CNOT gate between qubits 0 and 1
	q.print(); // Print the resulting state
	unsigned outcome = q.measure(3); // Measure qubit 3 and return the result
}
```

Full documentation for how to install and use the library is contained in the [documentation](https://qsl.readthedocs.io). If you want to install the library to give it a go, you can follow the instructions below. 

# Try it out: quick installation

If you want to get the library up and running quickly to try it out, follow these instructions. You need a linux distribution based on **Ubuntu 20.04 LTS** or more recent. If you are using a distribution based on Ubuntu 18.04, you can still install the library, but you will need to install some recent packages not on the package manager. Have a look at the instructions [below](#installation-using-ubuntu-1804-lts) for what to do.

To install the prerequisite `g++-10` compiler and `cmake` build system, run

```bash
sudo apt install gcc-10 g++-10 cmake
```

Now, to build and install the library, clone this repository and change into the top level directory (where this README is).

```bash
mkdir build # Make a build directory
cd build # Change to the build directory
CC=gcc-10 CXX=g++-10 cmake .. # Configure the build to use g++-10
cmake --build . # Build the library
sudo cmake --install . # Install the library to /usr/include/ and /usr/lib/
```

Now you should be able to compile the example above using the following line (assuming its in a file called main.cpp):

```bash
g++-10 -std=c++20 main.cpp -lqsl
```

If you encounter any errors, have a look at the documentation for more detailed installation instructions. 


# How to remove the library

If you want to remove the library that is installed using the steps above, run the following commands (double check you type them correctly!)

```bash
sudo rm -r /usr/include/qsl # Remove the include files
sudo rm /usr/lib/libqsl.so # Remove the shared object file
```

You can also delete the repository you cloned.

# Installation using Ubuntu 18.04 LTS

To build the library, you will need to install gcc-10, g++-10. Run the following commands to install the compilers:

```bash
# This will work on Ubuntu 18.04, but not on Ubuntu 16.04.
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt install gcc-10 g++-10
```

You will need to install cmake version 3.12.4 or more recent. See the [instructions for doing this](https://qsl.readthedocs.io/en/latest/install.html#cmake-at-least-version-3-12-4) in the main documentation.

Now you should be able to build the library using the following commands:

```bash
mkdir build # Make a build directory
cd build # Change to the build directory
CC=gcc-10 CXX=g++-10 cmake .. # Configure the build to use g++-10
cmake --build . # Build the library
sudo cmake --build . --target install # Install the library to /usr/include/ and /usr/lib/
```
