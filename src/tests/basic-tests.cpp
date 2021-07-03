// This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <qsl/utils/random.hpp>
#include <qsl/utils/quantum.hpp>
#include <qsl/utils/misc.hpp>
#include <qsl/qubits.hpp>

TEST_CASE ("Test the innerProduct function", "[quantum-utils]")
{
    const std::vector<qsl::complex<double>> a{ {1,1.2}, {1.2,1}, {1,1.5} };
    const std::vector<qsl::complex<double>> b{ {1,0}, {0,-1}, {-1,0} };
    const std::vector<qsl::complex<double>> c{
	{1.2,0}, {1.2,-2.1}, {1.3,2}, {0,0} };

    // Check that inner product throws exception for different sized vectors
    CHECK_THROWS(innerProduct(a,c));	 
		 
}


TEST_CASE ("Test the checkStateSize function", "[quantum-utils]")
{
    // Not a valid state size, should throw exception
    const std::vector<qsl::complex<double>> pretend_state{ {1,0}, {0,1}, {1,1} };
    CHECK_THROWS(checkStateSize(pretend_state));

    // Check that the function returns the correct state size
    const std::vector<qsl::complex<double>> state{
	{1,0}, {0,1}, {1,0}, {0,1} };

    REQUIRE( checkStateSize(state) == 2 );
}

TEST_CASE( "Test Fubini Study distance", "[quantum-utils]" )
{
    std::vector<qsl::complex<double>> a{ {1,0}, {0,1}, {1,1} }; 
    std::vector<qsl::complex<double>> double_a{ {2,0}, {0,2}, {2,2} }; 

    std::vector<qsl::complex<double>> b{ {1.3,2}, {1,1.3}, {1.5,0.2} }; 
    std::vector<qsl::complex<double>> double_b{ {2.6,4}, {2,2.6}, {3,0.4} }; 

    // Test that distance between equal and scaled vectors is zero
    REQUIRE( fubiniStudy(a,a) < 1e-13);
    REQUIRE( fubiniStudy(double_a,a) < 1e-13);

    double distance_a_b = fubiniStudy(a,b);
    double scaled_distance_a_b = fubiniStudy(a,double_b);

    REQUIRE( std::abs(distance_a_b - scaled_distance_a_b) < 1e-13);
    
}
TEST_CASE( "Qubits object constructors", "[constructors]" ) {

    const unsigned num_qubits{ 5 };
    qsl::Qubits<qsl::Type::Default, double> q{ num_qubits };

    // Set q to a random state
    const std::vector<qsl::complex<double>> state
	= qsl::makeRandomState<double>(num_qubits);
    q.setState(state);

    // Construct another object using the copy constructor
    qsl::Qubits<qsl::Type::Default, double> q_copy{ q };    

    // Check that the other properties are also the same
    REQUIRE(q.getNumQubits() == q_copy.getNumQubits());
    
    // Check that both objects have the same internal state
    double distance = qsl::distance(q.getState(), q_copy.getState());
    REQUIRE(std::abs(distance) < 1e-10);

    // Construct another object using initialisation from state vector
    qsl::Qubits<qsl::Type::Default, double> q_from_state{ state };    

    // Check that both objects have the same internal state
    distance = qsl::distance(q.getState(), q_from_state.getState());
    REQUIRE(std::abs(distance) < 1e-10);
      
}
