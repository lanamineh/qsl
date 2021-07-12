Gate Conventions
################

The object of the simulator is to include gates which are useful in the simulation of near-term quantum devices.

One-Qubit Gates
***************

The following single qubit quantum gates are implemented by the simulator

* **X-rotation, Y-rotation and Z-rotation**

  These gates, :math:`R_x(\theta)`, :math:`R_y(\theta)` and :math:`R_z(\theta)`, implement rotations of the input state about an axis of the Bloch sphere. The member functions are called ``rotateX``, ``rotateY`` and ``rotateZ``. They are defined as follows:

  .. math::
     \begin{align*}
     R_x(\theta) &= \begin{bmatrix}\cos(\theta/2)&-i\sin(\theta/2)\\-i\sin(\theta/2)&\cos(\theta/2)\end{bmatrix},\\
     R_y(\theta) &= \begin{bmatrix}\cos(\theta/2)&-\sin(\theta/2)\\\sin(\theta/2)&\cos(\theta/2)\end{bmatrix},\\
     R_z(\theta) &= \begin{bmatrix}e^{-i\theta/2}&0\\0&e^{i\theta/2}\end{bmatrix}.
     \end{align*}
   
  Sometimes an equivalent version of the Z rotation, the **phase shift**, is more useful. This is defined as followed:

  .. math::
     \text{phase}(\theta) = e^{i\theta/2}R_z(\theta) = \begin{bmatrix}1&0\\0&e^{i\theta}\end{bmatrix}.


  It is the same as a Z-rotation, up to the global phase :math:`e^{i\theta/2}`, which means interchanging them in a circuit will have no effect on measured outcomes. (However, different amplitudes in the state vector will have different phases depending which is used.) It is accessible using the ``phase`` member function.
  
* **Pauli X, Pauli Y and Pauli Z**

  These gates, :math:`\sigma_x`, :math:`\sigma_y` and :math:`\sigma_z`, are commonly used special cases of the general X-, Y- and Z-rotations, so they have their own member functions: ``pauliX``, ``pauliY`` and ``pauliZ``. They are defined as follows:

  .. math::
     \begin{align}
     \sigma_x &= iR_x(\pi) =\begin{bmatrix}0&1\\1&0\end{bmatrix},\\
     \sigma_y &= iR_y(\pi)= \begin{bmatrix}0&-i\\i&0\end{bmatrix},\\
     \sigma_z &= iR_z(\pi)= \begin{bmatrix}1&0\\0&-1\end{bmatrix}.
     \end{align}

  As before, the global phases of :math:`i` do not matter, so, for example, every time :math:`R_x(\pi)` is used it a circuit, it is valid to replace it with :math:`\sigma_x`.

* **Hadamard**

  The Hadamard gate is defined as follows:

  .. math::
     H = \frac{1}{\sqrt{2}}\begin{bmatrix}1&1\\1&-1\end{bmatrix},\\

  and is accessible using the ``hadamard`` member function. It can be written in terms of X- and Y-rotations as follows: :math:`H = \sigma_xR_y(\pi/2)`. THe Hadamard gate is mainly used for creating an equal superposition of states.

* **Arbitary unitary**

  The user can specify an arbitrary one-qubit unitary gate with the function
  ``unitary``. 
  
Two-Qubit Gates
***************

Controlled operations
---------------------

Most two-qubit operations of interest for near-term devices are controlled operations, which means that the state of one of the qubits (the control) determines whether a single qubit unitary is applied to another qubit (the target). If that unitary that is applied to the target qubit is :math:`U`, then the general form of a controlled operation is as follows

.. math::

   \text{controlled-}U = \begin{bmatrix}I&0\\0&U\end{bmatrix}

(This matrix assumes that the control is qubit 1 and the target is qubit 0. For more information about qubit indexing, see [link].) The simulator implements the following controlled operations:
   
* **Controlled X-rotation, Y-rotation and Z-rotation**

  These use :math:`U=R_x(\theta),R_y(\theta),R_z(\theta)` as defined above. The corresponding function names are ``controlRotateX``, ``controlRotateY``, ``controlRotateZ``.
  
* **Controlled-NOT**

  This is given by :math:`U=\sigma_x`. The corresponding function name is ``controlNot``.

* **Controlled-H**

  This is given by :math:`U=H`. The corresponding function name is ``controlHadamard``.

  
* **Controlled-Phase**

  This is given by :math:`U=\text{phase}(\theta)`. The corresponding function name is ``controlPhase``.
  
* **Controlled-Z**

  This is given by :math:`U=\sigma_z`. The corresponding function name is
  ``controlZ``.

* **Controlled-unitary**

  The user can specify an arbitrary controlled gate with a one-qubit unitary
  with ``controlUnitary``.


Number preserving
-----------------

The full Hilbert space :math:`\mathcal{H}` of :math:`n` qubits can be written as
a direct sum of the Hilbert spaces :math:`\mathcal{H}_i` with a fixed Hamming
weight :math:`i` for :math:`i=0,...,n`. Number preserving gates preserve
these Hilbert spaces. For example, the gate :math:`\sigma_Z` is number
preserving but :math:`\sigma_X` is not.

An arbitrary one-qubit number preserving gate is the ``phase`` gate. An arbitary
two-qubit number preserving gate preserves the :math:`|00\rangle`,
:math:`\{|01\rangle, |10\rangle \}` and :math:`|11\rangle` subspaces and takes
the form

 .. math::
    \begin{bmatrix}
    1 & 0 & 0 & 0 \\
    0 & u_{00} & u_{01} & 0 \\
    0 & u_{10} & u_{11} & 0 \\
    0 & 0 & 0 & e^{i\theta}
    \end{bmatrix}
  
where the 2x2 sub-matrix of :math:`u_{ij}` is a unitary matrix :math:`U`. Number
preserving gates are useful when simulating quantum chemistry on a quantum
computer and for this reason we also provide a number preserved simulator [link,
will explain why useful with NP simulator]. 

* **Number preserved X- and Y-rotations**

  These gates use :math:`U = R_x(\theta)` and :math:`U = R_y(\theta)`. The
  corresponding function names are ``npRotateX`` and ``npRotateY``. In the first
  case, the gate is equivalent to :math:`e^{-i\theta/2(XX+YY)}` and in the
  second :math:`e^{-i\theta/2(YX-XY)}`.
  
* **Swap and fermionic swap**

  Taking :math:`U=\sigma_x` swaps the state of two qubits, the function name for
  this is ``swap``. A fermionic swap gate acts as a swap gate but for fermions,
  this is useful when using fermion-to-qubit mappings when simulating quantum
  chemistry. It can be accessed with ``fswap`` and is equivalent to a swap gate
  followed by a controlled-Z gate. 

* **Number preserved Hadamard**

  This is given by :math:`U=H` and can be accessed with ``npHadamard``.

* **Arbitary number preserved gate**

  The user can specify an arbitrary :math:`U` with ``npUnitary``.
