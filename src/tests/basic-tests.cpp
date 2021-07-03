// This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <qsl/utils/random.hpp>
#include <qsl/utils/quantum.hpp>
#include <qsl/utils/misc.hpp>
#include <qsl/qubits.hpp>

TEST_CASE( "Qubits object copy constructor", "[constructors]" ) {

    const unsigned num_qubits{ 5 };
    qsl::Qubits<qsl::Type::Default, double> q{ num_qubits };

    // Set q to a random state
    q.setState(qsl::makeRandomState<double>(num_qubits));

    // Construct another object using the copy constructor
    qsl::Qubits<qsl::Type::Default, double> q_copy{ q };    

    // Check that the other properties are also the same
    REQUIRE(q.getNumQubits() == q_copy.getNumQubits());
    
    // Check that both objects have the same internal state
    double distance = qsl::distance(q.getState(), q_copy.getState());
    REQUIRE(std::abs(distance) < 1e-10);

}

TEST_CASE( "Test Fubini Study distance", "[quantum-utils]" ) {

    std::vector<qsl::complex<double>> a{ {1,0}, {0,1}, {1,1} }; 
    REQUIRE( fubiniStudy(a,a) < 1e-13);
}
