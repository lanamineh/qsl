/*
 *  Copyright 2021 Lana Mineh and John Scott
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 */
 
/**
 * \file test-utils.hpp
 * \brief Utilities for the test function
 *
 */

#include <armadillo>
#include <qsl/concepts.hpp>

/// Gets the bit in the nth position of val
unsigned getBit(unsigned val, unsigned n);

/// Set the nth bit of val to b
unsigned setBit(unsigned val, unsigned n, unsigned b);

/**
 * \brief Construct an armadillo matrix to apply an arbitrary gate
 */
template<std::floating_point Fp>
arma::SpMat<std::complex<Fp>>
makeMatrix(const arma::Mat<std::complex<Fp>> & gate, unsigned nqubits,
	   const std::vector<unsigned> & indices)
{
    // Check that the gate size matches the length of the indices
    if ((1ULL << indices.size()) != gate.n_rows) {
	throw std::logic_error("Length of indices vector does not "
			       "match the size of the gate in call "
			       "to makeMatrix()");
    }

    const std::size_t dim{ 1ULL << nqubits };
    arma::SpMat<std::complex<Fp>> mat(dim,dim);    
	
    // Write columns of mat. The columns of the matrix
    // correspond to the images of the basis states under
    // the action of the unitary matrix. 
    for (std::size_t col = 0; col < dim; col++) {

	// Get column bitstring and find the value of the ctrl and targ bits.
	// The values of the bits of col in the ctrl and targ positions
	// fix which column of the small matrix is used in populating
	// values of the big matrix at column col.
	//
	// Generalised for standard vector of qubit positions. 
	std::vector<unsigned> vals;
	for (std::size_t k = 0; k < indices.size(); k++) {
	    vals.push_back(getBit(col, indices[k]));
	}
	
	// Write rows of mat. For each row index, 
	for (std::size_t n = 0; n < gate.n_rows; n++) {

	    // Make the row index. The rows that are non-zero in a
	    // particular column are the ones whose bits agree with
	    // the bits in the col, apart from at the ctrl and targ
	    // positions. There, they take every possible ctrl and
	    // targ value.
	    std::size_t row = col;
	    for (std::size_t k = 0; k < indices.size(); k++) {
		row = setBit(row, indices[k], getBit(n,k));
	    }

	    // Make the column index for the small matrix
	    std::size_t m = 0;
	    for (std::size_t k = 0; k < vals.size(); k++) {
		m = setBit(m, k, vals[k]);
	    }

	    mat(row,col) = gate(n,m);
	}
    }
    return mat;
}

/**
 * \brief Convert a standard vector state to an armadillo vector
 */
template<std::floating_point Fp>
arma::Col<std::complex<Fp>> toArmaState(const std::vector<qsl::complex<Fp>> & res)
{
    arma::Col<std::complex<Fp>> qubit_v(res.size());
    for (std::size_t i = 0; i < res.size(); i++) {
	qubit_v(i) = std::complex<Fp>{res[i].real, res[i].imag};
    }
    return qubit_v;
}

/**
 * \brief Convert a Simulator state to an armadillo vector
 */
template<qsl::Simulator Sim>
arma::Col<std::complex<typename Sim::Fp_type>> toArmaState(const Sim & sim)
{
    using Fp = Sim::Fp_type;
    // Read qubit state into armadillo to check the state hasn't changed
    std::vector<qsl::complex<Fp>> res = sim.getState();
    return toArmaState<Fp>(res);
}

/**
 * \brief Make the projector onto a particular outcome of the target qubit
 *
 */
template<std::floating_point Fp>
arma::SpMat<std::complex<Fp>>
projector(unsigned num_qubits, unsigned targ, unsigned outcome)
{
    // Calculate the projector for the outcome, which is
    // the observable that has eigenvalue 0 for |~outcome) and eigenvalue
    // 1 for |outcome). The expectation value of this observable in a
    // state is the probability of getting 1 on measurement.
    arma::Mat<std::complex<Fp>> projector(2,2,arma::fill::zeros);
    if (outcome == 0) {
	projector(0,0) = 1;
    } else if (outcome == 1) {
	projector(1,1) = 1;
    } else {
	throw std::out_of_range("outcome must be 0 or 1 in projector() function");
    }
    arma::SpMat<std::complex<Fp>> M = makeMatrix(projector, num_qubits, {targ});
    return M;
}

/**
 * \brief Calculate the probability of measuring a projector outcome
 *
 * Uses the formula prob = (v|P|v), where v is the state and P is 
 * the projector
 */
template<std::floating_point Fp>
Fp probability(const arma::SpMat<std::complex<Fp>> & P,
	       const arma::Col<std::complex<Fp>> & v)
{
    // Calculate the probability of the projector outcome
    Fp prob = arma::cdot(v, P*v).real();
    return prob;
}

/**
 * \brief Collapse a state v using a projector P
 *
 * The output is the state Pv/|Pv| (i.e. the normalised projected state)
 *
 */
template<std::floating_point Fp>
arma::Col<std::complex<Fp>>
applyProjector(const arma::SpMat<std::complex<Fp>> & P,
	       const arma::Col<std::complex<Fp>> & v)
{
    arma::Col<std::complex<Fp>> state = (P * v)/arma::norm(P * v);
    return state;
}
