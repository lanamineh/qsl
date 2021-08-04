Style Guide
###########

This section contains the style guide for the code. 

General Guidelines
==================

* Do not use ``using namespace std;``, or any other variant, unless it is localised inside a function body and significantly increases readability.

* Try to pass large objects (such as ``std::vector<double>`` or other classes) by constant reference. Pass builtin types (e.g. ``int``, ``double``, etc.) by value. 

Indenting and whitespace
========================

Indents are 4 spaces.

The built-in constructions ``if``, ``for``, ``switch``, etc. should be followed by one space.

Statements which call a function, or define a callable object, should not have a space. For example, if ``foo`` is a function of ``int``, use ``foo(3)``, not ``foo (3)``. Function prototypes and definitions do not have a space after the function name.

Placement of braces
===================

The opening brace for ``if`` statements, ``for`` loops, ``switch`` statements, etc. should be on the same line as the opening statement:

.. code-block:: c++

   if (done == true) {
	return result;
   }

Braces should always be used in these constructions, even if there is only one line in the body of the construction as above. Otherwise, adding new lines to the body might lead to bugs if you forget to add the newly required braces.

In contrast, the opening brace for classes, structs, namespaces, and functions should begin directly below the keyword on the next line, as follows:

.. code-block:: c++

   void printInt(int a)
   {
	std::cout << a << std::endl;
   }

This automatically makes the function or class easier to read without needing to add an unnecessary empty line after the prototype line.
   
Class, function and variable names
==================================

.. todo:: Write me plz


Comments
========

Every file should contain a doxygen comment explaining what is inside the file, in the following format:

.. code-block:: c++

   /**
    * \file foobar.cpp
    * \brief Contains function implementations for foo and bar
    *
    * Put an optional longer comment about the contents of the
    * file here.
    *
    */

Everything in a ``*.hpp`` file (function prototypes, class/struct definitions, namespace definitions, ``using`` statements, etc.) must be fully doxygen commented (using an extended doxygen comment block). Each comment must contain a ``\brief`` and a detailed description if necessary. The comment should be sufficient to understand what the object is used for. Implementation details should not be included, unless they are relevant to the use of the object.

If functions take arguments, each argument must be documented using ``\param``. If the function does not return ``void``, the returned value must be documented using ``\return``.

Apart from the file doxygen comment at the top, ``*.cpp`` and ``*.tpp`` should not contain any comments outside function bodies. Function bodies can include single- or multi-line comments (constructed from repeated lines beginning with ``//``). Do not use doxygen for general comments. Do not to use multiline comment blocks (``/* */``), because this can make commenting out large blocks of code difficult. Only write a comment if it adds information which is not (too) obvious from the C++ code.

The only exception to the above rules is ``\todo`` statements, which can be included anywhere in any source file. A ``\todo`` statement should be written in the following way:

.. code-block:: c++

   ///\todo What needs to be done
		
Laying out Classes and Structs
==============================

Class/struct definitions should follow this convention regarding the use of whitespace:

.. code-block:: c++

   /**
    * \brief An example class
    *
    * Include a description about how to 
    * use the class.
    *
    */
   class Foo
   {
	int a; ///< Brief inline variable description
	double b; ///< Another brief description

	/**
	 * \brief A more complicated object
	 *
	 * Variables which required more description
	 * should have their own doxygen comment like
	 * this.
	 *
	 */
	Bar b;

	// Put private member functions here
	
   public:

	const int c; ///< Public member variables should be const
   
       /**
        * \brief Constructors should be listed before
	* other public member function 
        *
        * Member functions should always have an
	* extended doxygen comment, because it is
	* nearly always necessary to explain what
	* it does in more than one line!
	*
        */
	Foo();

	/**
	* \brief Get the Bar member
	* \return The internal Bar object
	*/
	Bar getBar() const;
   };

Observe the following guidelines about writing a class

* Do not use the keyword ``private``. Private members should be listed first anyway. Do not intermingle ``public`` and ``private`` members.
* Do not include an empty newline before the first class member or after the final class member.
* Public data members should be declared ``const``. If a data member is not const, but must be accessible from outside the class, then write a get method for it. If it is accessible from outside the class, it might be accidentally changed.
* All members of a class should be doxygen commented.
* If you can make a member function ``const``, you should. This allows the compiler to check that the member function does not modify anything in the class. 
* Do not include a member function definition in the class definition. Put them inside a corresponding ``.cpp`` file instead.

Structs can be used for objects where all the members are public,  whose only purpose is to group data togther. Structs should be kept as small as possible, and have no member functions other than constructors. Structs should be laid out and commented as described above, and should not contain the keyword ``public``.

Other names (namespaces, types, etc.)
=====================================

.. todo:: Write me plz

File naming and file structure
==============================

File names should comprise a single lowercase word with no underscores, followed by ``.hpp`` for header files, ``.cpp`` for implementation files or ``.tpp`` for template implementations. No file should exceed 1000 lines.

.. todo:: Is 'single lowercase word' too restrictive?

Header files
************
	  
Header files (``.hpp``) should contain class definitions and function prototypes, along with detailed documentation. In addition to exposing functions to ``.cpp`` files that use the functions, it should be possible for a user of the functions to read the documentation and understand what all the functions and classes do. Unrelated functions should be placed in different header files, and the filename should give some indication about what is in the file. The contents of the header file ``file.hpp`` should be wrapped in a header guard as follows:

.. code-block::

   #ifndef FILE_HPP
   #define FILE_HPP

   // File contents here
   
   #endif

Implementation files
********************
   
Implentation files (``.cpp``) should contain definitions of the class member functions and standalone functions defined in the header files. Implementation files should be split up into multiple files if they become too large. Implementation files do not need extensive comments.

Template files
**************

Template files (``.tpp``) should contain implementations of template functions, analogously to ``.cpp`` files. However, they should be wrapped in header guards and included at the end of their corresponding header file. Template files do not need extensive comments.

Notes
#####

Any general notes from the development should go here. They can be categorised properly later.

* Apparently GCC 10 produces much faster executables than GCC 7 (try compiling test in commit 42d491babc with both and see).


Todos
#####

Important
=========
	  
.. todo:: Write sorting method for sampling, figure out exactly when it is faster than the binary search method, and write a function to swap between the two.

.. todo:: Split up Qubits .cpp file into multiple chunks (e.g. for gates, for measurement, etc.)

.. todo:: Fix cmake build.

.. todo:: Find out why it takes longer to apply gates to qubit 0. Firstly, check that it does actually take longer (and is not an artifact due to 0 being the first qubit that is manipulated). In theory, it should be faster to apply the gates to low qubits.

.. todo:: Redo the variable names in the two qubit gates so that they are consistent with the rest of the code

.. todo:: Decide what to do when function prototypes line wrap (due to long arguments, long return types, templates, etc.)

.. todo:: Move the comments from benchmark1 and benchmark2 to Compare class, and delete benchmark1 and benchmark2 prototypes. Do the same with the other old benchmark functions

.. todo:: Rename ``Time`` class to something else (it conflicts with ``Timer``)

.. todo:: Check whether SingleSim can be implemented as a special case of MultiSim, with the length = 1 (i.e. check that it is not slower to use the vector). If it is fine, then the code for SingleSim could be realised as a special case of MultiSim.

.. todo:: Figure out what graphs we want to plot in python

.. todo:: Check that the mean and variance in the Results class give the right answers.

.. todo:: Simplify the implementation of the functions in the Results class
	  
Less important
==============

.. todo:: Find the proper way to include todos

.. todo:: Sort out the ``Fp_type`` type inside the Qubits class. Is there a way to eliminate it completely? 

.. todo:: Add more gates to the simulators.

.. todo:: Write functions to measure, postselect and sample groups of qubits.

.. todo:: Write a simulator that contains error checking (like bounds checking) for debugging, etc.

.. todo:: Somehow reduce the code repetition in Compare and Time classes

.. todo:: Implement the other functions in the Results class (but before that, decide whether the structure of the Results class is right. This will be obvious after trying to use them e.g. to plot graphs in python).
	  
.. doxygenpage:: todo

License
#######

The code and documentation is covered by the Apache 2.0 License below.

.. literalinclude:: ../LICENSE
