/* 
 * Copyright (C) 2020 Lana Mineh and John Scott.
 *
 * This file is part of QSL, the quantum computer simulator.
 *
 * QSL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QSL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QSL.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * \file quantum.hpp
 * \brief General quantum-related functions
 *
 */

#ifndef QUANTUM_HPP
#define QUANTUM_HPP

#include <vector>
#include "complex.hpp"

namespace qsl {

    /**
     * \brief Compute the inner product between two state vectors
     *
     * The standard inner product between two state vectors \f$u\f$ 
     * and \f$v\f$ is given by
     *
     * \f[
     * \langle u | v \rangle = \sum_{n=0}^N u_i^* v_i
     * \f]
     *
     * This is the natural inner product on \f$ \mathbb{C}^N\f$,
     * analogous to the Euclidean scalar product on 
     * \f$ \mathbb{R}^N \f$. The inner product gives rise to a 
     * norm, which measures the magnitude of the vector in 
     * \f$ \mathbb{C}^N\f$.
     *
     * The inner product is also used in the definition of the
     * Fubini-Study metric, which is a distance between rays in
     * the projective space of physical states. 
     *
     */
    template<std::floating_point Fp = double>
    complex<Fp> innerProduct(const std::vector<complex<Fp>> & v,
			     const std::vector<complex<Fp>> & w);

    /**
     * \brief Compute the norm of a complex vector
     *
     * The norm of a complex vector v is given by the following expression 
     *
     * \f[
     * \lVert v\rVert = \sqrt{\langle v|v \rangle} 
     *  = \sqrt{\sum_{n=0}^N \left|v_i\right|^2}.
     * \f]
     *
     * which is derived from the inner product on the complex vector space.
     *
     */
    template<std::floating_point Fp = double>
    Fp norm(const std::vector<complex<Fp>> & v);

    /**
     * \brief Compute the norm of a real vector 
     *
     * The norm of a vector v is given by the following expression 
     *
     * \f[
     * \lVert v\rVert = \sqrt{\langle v|v \rangle} 
     *  = \sqrt{\sum_{n=0}^N \left|v_i\right|^2}.
     * \f]
     */
    template<std::floating_point Fp>
    Fp norm(const std::vector<Fp> & v)
    {
	double norm = 0;
	for(std::size_t n=0; n<v.size(); n++) {
	    norm += v[n] * v[n];
	}

	return std::sqrt(norm);
    }



    /**
     * \brief Compute the Fubini-Study metric between two state vectors
     *
     * The Fubini-Study metric is a distance between rays in complex projective
     * space which can be interpreted as the angle between rays.
     *
     * The definition of the Fubini-Study is
     *
     * \f[
     * \gamma (u ,v ) = \arccos {\sqrt{\frac{\langle u \vert v \rangle \;
     * \langle v \vert u \rangle }{\langle u \vert u \rangle \;\langle 
     * v \vert v \rangle }}}
     * = \arccos\left[{\frac{|\langle u|v \rangle|}{\lVert u\rVert\lVert 
     * v\rVert}}\right]
     * \f]
     *
     * This is a more appropriate distance to use for physical states
     * because it compares them up to a global scale factor (e.g. a global
     * phase).
     *
     * To compare two vectors for equality, it is more appropriate to use
     * the complex Euclidean norm which does not ignore global scale factors.
     *
     * For more information, see 
     * <a href="https://en.wikipedia.org/wiki/Fubini%E2%80%93Study_metric">
     * this wikipedia page
     * </a>
     *
     */
    template<std::floating_point Fp = double>
    Fp fubiniStudy(const std::vector<complex<Fp>> & v,
		   const std::vector<complex<Fp>> & w);

    /**
     * \brief Normalise the state vector.
     */
    template<std::floating_point Fp = double>
    Fp normalise(std::vector<complex<Fp>> &state);

    /**
     * \brief Check whether the state vector is a valid size
     *
     * If it is, return the number of qubits associated to the state
     * vector. If not, throw an exception.
     *
     */
    template<std::floating_point Fp = double>
    unsigned checkStateSize(const std::vector<complex<Fp>> & state);

    /**
     * \brief Check if the state is number preserving and find the number
     * of ones if it is.
     *
     * This function is probably inefficient - needs improving.
     * Also might make it a member function of the class.
     *
     * Throws std::logic_error if state is not number preserving.
     */
    template<std::floating_point Fp = double>
    unsigned checkStateNP(const std::vector<qsl::complex<Fp>> & state);
    
}
    
#endif
