The Qubits object
#################

The most important class in the simulator is called **Qubits**, which holds all the information for the current simulation. The object depends on an enum template parameter, Type, which determines the simulator implementation.

.. doxygenenum:: Type

Default Simulator
#################

The most basic simulator is called Qubits<Type::Default>. It does not use any special optimisations (like inline assembly), and does not support threading.

.. doxygenclass:: Qubits< Type::Default, Fp >

Constructing the object
=======================

There are several ways to construct a qubits object, which are summarised in the example below:

.. code-block:: c++

    // Make a simulator object with 2 qubits, using
    // double precision The object is initialised
    // to the all zero state.
    Qubits<Type::Default, double> q1(2);

    // To initialise an object in an arbitrary starting
    // state, first make a std::vector containing the
    // state and then initialise from the std::vector.
    // The constructor checks that the state vector is
    // a power of two. If not, an exception is thrown.
    std::vector<complex<double>> state{ {0,0}, {1,0}, {0,0}, {0,0} };
    Qubits<Type::Default, double> q2(state);

    // You can construct a new Qubits object by
    // copying the state of an already existing
    // object like this
    Qubits<Type::Default, double> q3(q2);

    // You can also assign one qubits object to
    // another, provided they contain the same
    // number of qubits. An exception is thrown
    // if the number of qubits is different.
    q3 = q1;

    // All the above comments apply to Qubits
    // object with different precision, but
    // you cannot mix different precisions
    Qubits<Type::Default, float> q4(3);
    Qubits<Type::Default, float> q5(q4); // Copy q4 to q5

The full documentation for all the constructors is given below.
  
.. doxygengroup:: qubits_constructors
   :content-only:

Quantum Gates
=============

The Qubits class implements a selection of one- and two- qubit gates, which are described below. The gates can be called on an instance of the Qubits object as in the following example:

.. code-block:: c++

   // Make a Qubits object with 3 qubits using double precision
   Qubits<Type::Default, double> q(3);
   
   q.pauliX(0); // Perform pauliX on the zeroth qubit 
   q.rotateX(1, 0.5); // Rotate qubit one about the X axis by 0.5 rad

.. doxygengroup:: qubits_gates
   :content-only:

Measurement
===========

After performing a series of gates, you can measure or otherwise analyse the state of the qubits using the functions below.

.. doxygengroup:: qubits_meas
   :content-only:

Utilities
=========

These functions can be used to get information about the Qubits object.

.. doxygengroup:: qubits_utils
   :content-only:
      

OpenMP-based simulator
######################

This simulator makes use of parallelisation using the OpenMP library. 

.. doxygenclass:: Qubits< Type::Omp, Fp >

Constructing the object
=======================

.. doxygengroup:: qubits_omp_constructors
   :content-only:

Quantum gates and Measurement
=============================

The member functions for applying gates and performing measurement use the same syntax as those of the default simulator object.
      
