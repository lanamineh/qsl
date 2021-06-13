This section contains miscellaneous functions and classes that are used in the implementation of the simulator

Quantum related functions
#########################
.. doxygenfile:: quantum.hpp

Random numbers
##############

Random numbers are used in several contexts in the simulator. Firstly, quantum measurement simulation requires generation of random numbers to determine whether to pick one outcome or another. Secondly, benchmarking and verification functions generate random samples on which to test the simulator.

The Random class is used in the Qubits object to generate the random numbers necessary for simulating measurements.

.. doxygenclass:: qsl::Random
   :members:	      

Other functions are available for generating random quantum states, generating lists of random phases, etc. as follows.
      
.. doxygenfile:: random.hpp
   :sections: func
	      
Complex numbers
###############		 
.. doxygenfile:: complex.hpp

Miscellaneous
#############
.. doxygenfile:: misc.hpp
