#include <iostream>
#include <map>
#include <qsl/qubits.hpp>

int main ()
{
    // Create a Qubits object with 3 qubits
    std::cout << "Creating Qubits object..." << std::endl;
    unsigned nqubits = 3;
    Qubits<Type::Default, double> q{nqubits};

    // Print the state of the Qubits object
    std::cout << std::endl
	      << "The Qubits object starts in the all-zero state:"
	      << std::endl;
    q.print();
    
    // Do a Pauli-X gate on the 3rd qubit
    // (remember the first qubit is 0)
    std::cout << "Applying Pauli-X to qubit 2:" << std::endl;
    q.pauliX(2);

    // Print the new state
    q.print();

    // You can repeatedly sample from the state of the
    // Qubits object, mimicking what would happen if
    // you re-prepared the state and measured it
    // many times.
    std::cout << "Sampling from that state 10 times..." << std::endl;
    const std::size_t num_samples = 10;
    std::map<std::size_t, std::size_t> samples = q.sampleAll(num_samples);

    // Since the Qubits object is in a computation basis
    // state, you will get the same outcome every time
    for (auto & [outcome, frequency] : samples) {
	std::cout << "Outcome " << outcome << " occured "
		  << frequency << " time(s)" << std::endl;
    }

    // Now apply a rotation to create a superposition state
    std::cout << std::endl
	      << "Applying an X rotation..." << std::endl;
    double angle = M_PI/2; // Radians
    q.rotateX(1, angle);
    q.print();

    // This time, when you sample 10 times, you will get
    // different outcomes
    std::cout << "Sampling from that state a further 10 times..." << std::endl;
    samples = q.sampleAll(num_samples);

    // Print the samples
    for (auto & [outcome, frequency] : samples) {
	std::cout << "Outcome " << outcome << " occured "
		  << frequency << " time(s)" << std::endl;
    }
    
    
}
