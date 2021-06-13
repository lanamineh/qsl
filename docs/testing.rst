This section is all about verifying that the simulator works correctly, and analysing its performance.

Verification
############

A single class, ``Verify``, checks that the simulators behave correctly. There are two types of verification that can be performed:

* Comparison of the results of a simulator with a known true result
* Comparison of the results from two different simulators
    
The first test is a stronger test for correctness, but the second check is often sufficient (iff all the simulators agree, then it is likely they are all correct).

Some of the operations that the simulators perform are deterministic, while others aren't. For the deterministic operations, such as gate applications, it is possible to check that the result is equal to another given result. For probabilistic operations (such as measurement, which involves a random collapse of the state vector), it is necessary to perform a large number of identical tests and check that the average result is close to another given result.


Benchmarking
############
.. doxygenfile:: bench.hpp

Google qsim
###########

Google wrote a simulator called qsim to help validating the superconducting qubit supremacy experiment in 2019. Their `documentation <https://quantumai.google/qsim>`_ explains the use of the program as a C++ library, which involves copying the contents of the ``lib/`` folder to a ``qsim/`` folder that is in the search path of your C++ compiler (on Linux, ``/usr/local/include/`` will do). Then library files can be included using ``#include <qsim/gates_qsim.h>`` in a C++ file.

		 
Quest
#####

We have used the `quest <https://quest.qtechtheory.org>`_ simulator to verify the results of |project|. Quest is predominatly written in C, so we have wrapped it in a C++ class to make it easier to use.

.. doxygenclass:: Quest


Qiskit
######

Clone the Qiskit Aer repository ``git clone
https://github.com/qiskit/qiskit-aer`` into a source folder. Rename
``CMakeLists.txt`` in the qiskit repository to ``CMakeLists.txt.old`` and
replace with ``CMakeLists.txt.qiskit`` in the qsl repository. Then,

1) Install conan as it is required to make qiskit with ``python3 -m pip install
   conan``.
2) Download the latest version of nlohmann-json and put in your include
   path. (Explain this more!)
3) Make a ``build`` directory inside the qiskit folder.
4) Change into the build directory and run ``cmake ..``. This step does not seem
   to work, will need to fix it.
