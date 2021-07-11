#include <gtest/gtest.h>
#include <qsl/qubits.hpp>
#include <armadillo>
#include <complex>

TEST(Gates, OneQubitNoArg)
{   
    const unsigned num_qubits = 8;
    const unsigned targ = 4;
    using Fp = double;
    using Sim = qsl::Qubits<qsl::Type::Default, Fp>; 

    // Create armadillo matrices
    arma::SpMat<std::complex<Fp>> pauliX(2, 2);
    pauliX(0, 1) = 1;
    pauliX(1, 0) = 1;
    
    arma::SpMat<std::complex<Fp>> hadamard(2, 2);
    Fp sqrt2 = 1/std::sqrt(2);
    hadamard(0, 0) = sqrt2;
    hadamard(0, 1) = sqrt2;
    hadamard(1, 0) = sqrt2;
    hadamard(1, 1) = -sqrt2;

    // Create list of gates mapped to matrices
    std::vector<std::pair<qsl::Gate<Sim, unsigned>, 
			  arma::SpMat<std::complex<Fp>>>> gates;
    gates.push_back({&Sim::pauliX, pauliX});
    gates.push_back({&Sim::hadamard, hadamard});

    for (const auto & [fn, mat] : gates) {    
	// Make a random state
	Sim q{num_qubits};
	const std::vector<qsl::complex<Fp>> state
	    = qsl::makeRandomState<Fp>(num_qubits);
	q.setState(state);

	// Set an armadillo vector to the same state
	std::size_t dim = 1 << num_qubits;
	arma::Col<std::complex<Fp>> v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    v(i) = std::complex<Fp>{state[i].real, state[i].imag};
	}

	// Apply gate to qubits
	std::invoke(fn, q, targ);

	// Create gate in armadillo
	// Sizes of idenity matrix padding
	std::size_t pre = 1 << targ;
	std::size_t post = 1 << (num_qubits - targ - 1);
	// Tensor to make the gate
	arma::SpMat<std::complex<Fp>> gate =
	    arma::speye<arma::SpMat<std::complex<Fp>>>(pre, pre);
	gate = arma::kron(mat, gate);
	gate = arma::kron(arma::speye<arma::SpMat<std::complex<Fp>>>(post, post),
			  gate);

	// Apply gate in armadillo
	v = gate * v;

	// Read qubit state into armadillo
	std::vector<qsl::complex<Fp>> res = q.getState();
	arma::Col<std::complex<Fp>> qubit_v(dim);
	for (std::size_t i = 0; i < dim; i++) {
	    qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
	}
    
    	EXPECT_TRUE(arma::approx_equal(v, qubit_v, "both", 1e-8, 1e-10));
    }
}
