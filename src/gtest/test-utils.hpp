/**
 * \file test-utils.hpp
 * \brief Utilities for the test function
 *
 */

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
