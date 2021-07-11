#include <gtest/gtest.h>
#include <qsl/qubits.hpp>
#include <armadillo>
#include <complex>

TEST(OneQubitGates, pauliX)
{
    const unsigned num_qubits = 10;
    const unsigned targ = 7;
    using Fp = float;

    // Make a random state
    qsl::Qubits<qsl::Type::Default, Fp> q{num_qubits};
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
    q.pauliX(targ);

    // Create gate in armadillo
    arma::SpMat<Fp> X(2, 2);
    X(0, 1) = 1;
    X(1, 0) = 1;
    // Sizes of idenity matrix padding
    std::size_t pre = 1 << targ;
    std::size_t post = 1 << (num_qubits - targ - 1);
    // Tensor to make the gate
    arma::SpMat<Fp> gate = arma::speye<arma::SpMat<Fp>>(pre, pre);
    gate = arma::kron(X, gate);
    gate = arma::kron(arma::speye<arma::SpMat<Fp>>(post, post), gate);

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
