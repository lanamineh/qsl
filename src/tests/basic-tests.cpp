// This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <qsl/utils/random.hpp>
#include <qsl/utils/quantum.hpp>
#include <qsl/utils/misc.hpp>
#include <qsl/qubits.hpp>

TEST_CASE( "Test Fubini Study distance", "[quantum-utils]" ) {

    std::vector<qsl::complex<double>> a{ {1,0}, {0,1}, {1,1} }; 
    REQUIRE( fubiniStudy(a,a) < 1e-13);
}
