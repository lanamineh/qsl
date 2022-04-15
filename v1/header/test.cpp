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
#include <iostream>
#include <random>

int main()
{
    // Constructors
    qsl::basic<double> q1{3}; // A simple qsl simulator containing three qubits
    qsl::resize<double, true> q2{5}; // A resizeable simulator of 5 qubits, with debugging
    qsl::number<double, false, qsl::seq> q3{2}; // Two qubit sequential fixed-number simulator

    //qsl::basic<float> q4{{1,0},{1,0}}; // You cannot implicitly initialise a std::vector
    std::vector<std::complex<float>> s1{{1,0},{1,0}}; // |+) state
    qsl::basic<float> q4{s1}; // You can initialise explicitly from a std::vector    

    // Convert one simulator to another
    qsl::basic<double> q5{q3}; // Any simulator type to basic
    qsl::resize<double> q6{q1}; // Any simulator type to resize
    qsl::number<double,true,qsl::omp> q7{q3}; // Any number sim to number sim
    //qsl::basic<float> q8{q3}; // Cannot implicitly change floating point type    

    // Copy-assignment operator
    q5 = q1;

    // Copy constructor (actually using the coversion constructor)
    qsl::basic<double> q8{q1};

    // Move semantics
    qsl::basic<double> q9{qsl::basic<double>{3}};

    // Conversion between floating point types
    //qsl::basic<float> q10{q1}; // This one is not allowed
    [[maybe_unused]] qsl::basic<float> q10{q1}; // This one seems to use the conversion operator, even though it is explicit
    [[maybe_unused]] qsl::basic<float> q12(q1); // This one also works
    //qsl::basic<float> q14 = q1; // This one is not OK for some reason
    // qsl::basic<float> q13{q3}; // Cannot do this -- constructor is explicit
    [[maybe_unused]] qsl::basic<float> q11{4};
    //q11 = q1; // Cannot do this because conversion is explicit
    q11 = static_cast<qsl::basic<float>>(q1);

    [[maybe_unused]] qsl::resize<long double> q15 = static_cast<qsl::resize<long double>>(q6);
    //qsl::resize<long double> q16 = q6; // This one does not work

    [[maybe_unused]] qsl::number<long double> q17 = static_cast<qsl::number<long double>>(q3);
    //qsl::resize<long double> q18 = q3; // This one does not work

    
    // Gates
    for (unsigned n{0}; n < q1.size(); ++n) {
	q1.h(n); // Apply a Hadamard gate to all qubits in q1
    }
    q2.cnot(1,2); // Apply a CNOT gate, control = 1, target = 2
    q3.z(1); // Apply a Pauli Z gate to qubit 1
    q1.rx(0, 0.3); // Apply an Rx rotation of 0.3 radians to qubit 0
    // q3.cnot(0,1); // This should not work
    
    // Basic functions
    [[maybe_unused]] std::size_t d{q1.dim()}; // Get the dimension of the state vector
    [[maybe_unused]] unsigned n1{q2.size()}; // Get the number of qubits
    [[maybe_unused]] unsigned n2{q3.get_ones()}; // Get the number of ones in a fixed-number simulator
    
    std::vector<std::complex<double>> s2{{1,0},{0,0},{0,0},{1,0}};
    q2.set_state(s2); // You can set the state to an arbitrary number of qubits
    // q4.set_state(s2); // You cannot implicitly change the floating point type

    // TODO
    //q3.make_random(); // Set the simulator to a random state

    std::mt19937_64 g1;
    q1.make_random(g1); // Make a random state using your own generator
    q2.make_random(g1);
    q3.make_random(g1);
    q4.make_random(g1);

    struct VeryRandom
    {
	static constexpr unsigned min() { return 0; };
	static constexpr unsigned max() { return 5; };
	unsigned operator() () { return 3; } // Very random
    } g2;
    q1.make_random(g2); // Make a random state using a custom generator
    
    std::cout << q1 << std::endl; // Send any simulator to an ostream
    std::cout << q2 << std::endl;
    std::cerr << q3 << std::endl;
    std::cerr << q4 << std::endl;
    q2.print(); // Print directly to console output (with a trailing newline)
    q4.print(std::cerr); // Pass a custom output stream

    std::cout << qsl::distance(q1, q2) << std::endl;
    std::cout << qsl::distance(q2, q3) << std::endl;
    std::cout << qsl::distance(q1, q3) << std::endl;
    // std::cout << qsl::distance(q2, q4) << std::endl; // This one (correctly)  doesn't work -- double and float conflict

    std::cout << qsl::distance(q4, s1) << std::endl;    
    std::cout << qsl::distance(q1, s2) << std::endl;    
    std::cout << qsl::distance(s1, q4) << std::endl;    
    std::cout << qsl::distance(s2, q1) << std::endl;    
    // std::cout << qsl::distance(q1, s1) << std::endl; // This one won't work  

    std::cout << qsl::fidelity(q1, q2) << std::endl;
    std::cout << qsl::fidelity(q2, q3) << std::endl;
    std::cout << qsl::fidelity(q1, q3) << std::endl;
    // std::cout << qsl::fidelity(q2, q4) << std::endl; // This one (correctly)  doesn't work -- double and float conflict

    std::cout << qsl::fidelity(q4, s1) << std::endl;    
    std::cout << qsl::fidelity(q1, s2) << std::endl;    
    std::cout << qsl::fidelity(s1, q4) << std::endl;    
    std::cout << qsl::fidelity(s2, q1) << std::endl;    
    // std::cout << qsl::fidelity(q1, s1) << std::endl; // This one won't work  

    std::cout << qsl::inner_prod(q1, q2) << std::endl;
    std::cout << qsl::inner_prod(q2, q3) << std::endl;
    std::cout << qsl::inner_prod(q1, q3) << std::endl;
    // std::cout << qsl::inner_prod(q2, q4) << std::endl; // This one (correctly)  doesn't work -- double and float conflict

    std::cout << qsl::inner_prod(q4, s1) << std::endl;    
    std::cout << qsl::inner_prod(q1, s2) << std::endl;    
    std::cout << qsl::inner_prod(s1, q4) << std::endl;    
    std::cout << qsl::inner_prod(s2, q1) << std::endl;    
    // std::cout << qsl::inner_prod(q1, s1) << std::endl; // This one won't work  

    // Compare different floating point types
    std::cout << qsl::inner_prod(q1,static_cast<qsl::basic<double,false>>(q10)) << std::endl;
    
}
