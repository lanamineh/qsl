/**
 * \file test.cpp
 * \brief Static test of usability of QSL interface
 *
 * This file contains things that are possible using the interface defined in qsl.hpp, for
 * the purpose of checking usability. The code can be compiled and statically analysed,
 * but not linked (because qsl.hpp lacks an implementation).
 */

#include "qsl.hpp"
#include <vector>

int main()
{
    qsl::basic<double> q1{3}; // A simple qsl simulator containing three qubits
    qsl::resize<double, true> q2{5}; // A resizeable simulator of 5 qubits, with debugging
    qsl::number<double, false, qsl::seq> q3{2}; // Two qubit sequential fixed-number simulator

    //qsl::basic<float> q4{{1,0},{1,0}}; // You cannot implicitly initialise a std::vector
    std::vector<std::complex<float>> s1{{1,0},{1,0}}; // |+) state
    qsl::basic<float> q4{s1}; // You can initialise explicitly from a std::vector    
    
    for (unsigned n{0}; n < q1.size(); ++n) {
	q1.h(n); // Apply a Hadamard gate to all qubits in q1
    }
    q2.cnot(1,2); // Apply a CNOT gate, control = 1, target = 2
    q3.z(1); // Apply a Pauli Z gate to qubit 1
    q1.rx(0, 0.3); // Apply an Rx rotation of 0.3 radians to qubit 0

    [[maybe_unused]] std::size_t d{q1.dim()}; // Get the dimension of the state vector

    std::vector<std::complex<double>> s2{{1,0},{0,0},{0,0},{1,0}};
    q2.set_state(s2); // You can set the state to an arbitrary number of qubits
    // q4.set_state(s2); // You cannot implicitly change the floating point type
    
    
    
}
