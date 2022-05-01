/**
 * \file qsl.hpp
 * \brief First draft of qsl interface.
 */

#include <concepts>
#include <vector>
#include <complex>
#include <iostream>
#include <filesystem>
#include <random>
#include <map>

namespace qsl
{
    template<typename T>
    struct get_precision {};
    
    template<template<std::floating_point,bool,typename> typename S, std::floating_point F, bool D, typename P>
    struct get_precision<S<F, D, P>>
    {
	using type = F;
    };

    // \todo This might not be necessary if we can find a way around the get_precision_t compile error
    template<std::floating_point F>
    struct get_precision<std::vector<std::complex<F>>>
    {
	using type = F;
    };
    
    template<typename T>
    using get_precision_t = get_precision<T>::type;

    
    template<typename S, typename F>
    concept has_state_vector_of_type = requires (S s) {
    	{s.get_state()} -> std::same_as<std::vector<std::complex<F>>>;
    };

    template<typename S>
    concept has_state_vector = requires (S s) {
    	{s.get_state()};// -> Check that return type is std::vector of std::complex of std::floating_point 
    };
    
    template<typename S1, typename S2>
    concept have_same_precision = std::is_same_v<get_precision_t<S1>, get_precision_t<S2>>;
    
    /// \todo Make these all part of a concept
    /// Run functions without parallelisation
    struct seq {};
    /// Run functions with omp
    struct omp {};
    /// Automatically turn omp on and off depending on number of qubits
    struct opt {};

    /**
     * \brief Random generator
     */
    class gen_t
    {
    public:    
	using result_type = std::mt19937_64::result_type;
	/// Return minimum possible value
	static result_type min();
	/// Return maximum possible value
	static result_type max();
	/// Obtain the next3 sample from the generator
	result_type operator() ();

	/// Instantiate generator with a random seed
	gen_t(); 
	/// Instantiate generator with a specific seed
	explicit gen_t(result_type seed);
	/// Return the seed
	result_type seed() const;
	/// Set a new seed
	void seed(result_type seed);
    };
    /// \todo Figure out what keyword needs to go in front
    extern gen_t gen;

    /**
     * \brief General purpose quantum simulator.
     */
    template<std::floating_point F, bool D = false, typename P = opt>
    class basic
    {
    public:
	//using precision_t = F;
	
	/// Instantiate a simulator with num_qubits qubits in the all zero state
	explicit basic(unsigned num_qubits);
	/// Instantiate a simulator based on the state that is passed in 
	explicit basic(const std::vector<std::complex<F>> & state);

	/// Convert from any simulator of the same type
	/// This is not a  copy constructor (because it is templated), so it will not
	/// cause the move constructors to be implicitly deleted.
	/// TODO concept for simulator
	template<has_state_vector_of_type<F> S>
	explicit basic(const S & s);

	/// Allow explicit cast between different simulator types
	template<std::floating_point F1, bool D1, typename P1>
	explicit operator basic<F1,D1,P1>();
	
	/// Get the number of qubits
	unsigned size() const;
	/// Get the dimension of the Hilbert space
	unsigned dim() const;
	/// Return std::vector of state
	std::vector<std::complex<F>> get_state() const;
	/// Set a state
	void set_state(const std::vector<std::complex<F>> & state);
	/// Reset to the all zero computational basis state
	void reset();
	/// Access state vector elements (read-only). 
	const std::complex<F> & operator[](std::size_t index) const;
	/// Generate random state, need to pass in a generator which will default to our global one.
	void make_random(std::uniform_random_bit_generator auto & g = gen)
	    {
		// In order for this to work, the body of this function must
		// be in the header file (if the user is allowed to pass custom
		// gen types)
		g(); // Just to surpress unused-variable warning for now
	    }
	
	/// Print state vector
	void print(std::ostream & os = std::cout) const;
	/// Save state to json
	void save_json(const std::filesystem::path & path) const;
	/// Load state from json
	void load_json(const std::filesystem::path & path);
	/// \todo Select a compressed format
	/// Save state to compressed format
	void save_hdf5(const std::filesystem::path & path) const;
	/// Load state from compressed format
	void load_hdf5(const std::filesystem::path & path);

	// One-qubit gates
	void rx(unsigned targ, F angle);
	void ry(unsigned targ, F angle);
	void rz(unsigned targ, F angle);  
	void phase(unsigned targ, F angle); 
	void h(unsigned targ);
	void x(unsigned targ);
	void y(unsigned targ);
	void z(unsigned targ);
	void u1(unsigned targ, const std::vector<F> & matrix);
	void u1(unsigned targ, const std::vector<std::complex<F>> & matrix);

	// Controlled gates
	void crx(unsigned ctrl, unsigned targ, F angle);
	void cry(unsigned ctrl, unsigned targ, F angle);
	void crz(unsigned ctrl, unsigned targ, F angle);  
	void cphase(unsigned ctrl, unsigned targ, F angle); 
	void ch(unsigned ctrl, unsigned targ);
	void cnot(unsigned ctrl, unsigned targ);
	void cy(unsigned ctrl, unsigned targ);
	void cz(unsigned ctrl, unsigned targ);
	void cu1(unsigned ctrl, unsigned targ, const std::vector<F> & matrix);
	void cu1(unsigned ctrl, unsigned targ, const std::vector<std::complex<F>> & matrix);

	// Rest of the number gates
	void nrx(unsigned targ1, unsigned targ2, F angle);
	void nry(unsigned targ1, unsigned targ2, F angle);
	void nrz(unsigned targ1, unsigned targ2, F angle);
	void swap(unsigned targ1, unsigned targ2);
	void fswap(unsigned targ1, unsigned targ2);
	void iswap(unsigned targ1, unsigned targ2);
	void nh(unsigned targ1, unsigned targ2);
	void nu1(unsigned targ1, unsigned targ2, const std::vector<F> & matrix);
	void nu1(unsigned targ1, unsigned targ2, const std::vector<std::complex<F>> & matrix);

	/// Arbitrary two-qubit unitary
	void u2(unsigned targ1, unsigned targ2, const std::vector<F> & matrix);
	void u2(unsigned targ1, unsigned targ2, const std::vector<std::complex<F>> & matrix);
	
	/// Calculate probability of specific outcome of a qubit
	F prob(unsigned targ, unsigned outcome) const;
	/// Measure one qubit
	unsigned measure(unsigned targ, std::uniform_random_bit_generator auto g = gen);
	/// Measure all of the qubits
	std::size_t measure_all(std::uniform_random_bit_generator auto g = gen);
	/// Postselect one qubit on a specific outcome
	unsigned postselect(unsigned targ, unsigned outcome) const;
	/// Perform measurement sampling on one qubit
	std::vector<std::size_t> sample(unsigned targ, std::size_t samples,
					std::uniform_random_bit_generator auto g = gen) const;
	/// Perform measurement sampling on all the qubits
	std::map<std::size_t, std::size_t> sample_all(std::size_t samples,
						      std::uniform_random_bit_generator auto g = gen) const;

    };


    /**
     * \brief Resizeable quantum simulator
     */
    template<std::floating_point F, bool D = false, typename P = opt>
    class resize
    {
    public:
	/// Instantiate a simulator with num_qubits qubits in the all zero state
	explicit resize(unsigned num_qubits);
	/// Instantiate a simulator based on the state that is passed in 
	explicit resize(const std::vector<std::complex<F>> & state);

	/// Construct (a copy) from any simulator of the same type
	/// TODO concept for simulator
	template<has_state_vector_of_type<F> S>
	explicit resize(const S & s);

	/// Allow explicit cast between different simulator types
	template<std::floating_point F1, bool D1, typename P1>
	explicit operator resize<F1,D1,P1>();
	
	/// Get the number of qubits
	unsigned size() const;
	/// Get the dimension of the Hilbert space
	unsigned dim() const;
	/// Return std::vector of state
	std::vector<std::complex<F>> get_state() const;
	/// Set a state
	void set_state(const std::vector<std::complex<F>> & state);
	/// Reset to the all zero computational basis state
	void reset();
	/// Access state vector elements (read-only). 
	const std::complex<F> & operator[](std::size_t index) const;
	/// Generate random state, need to pass in a generator which will default to our global one.
	void make_random(std::uniform_random_bit_generator auto & g = gen)
	    {
		// In order for this to work, the body of this function must
		// be in the header file (if the user is allowed to pass custom
		// gen types)
		g(); // Just to surpress unused-variable warning for now
	    }


	/// Add a qubit to the end of the state vector
	void add_qubit();
	/// Add a qubit at index targ
	void add_qubit(unsigned targ);
	/// Trims the state vector to the current qubit size
	void trim();
	
	/// Print state vector
	void print(std::ostream & os = std::cout) const;
	/// Save state to json
	void save_json(const std::filesystem::path & path) const;
	/// Load state from json
	void load_json(const std::filesystem::path & path);
	/// \todo Select a compressed format
	/// Save state to compressed format
	void save_hdf5(const std::filesystem::path & path) const;
	/// Load state from compressed format
	void load_hdf5(const std::filesystem::path & path);

	// One-qubit gates
	void rx(unsigned targ, F angle);
	void ry(unsigned targ, F angle);
	void rz(unsigned targ, F angle);  
	void phase(unsigned targ, F angle); 
	void h(unsigned targ);
	void x(unsigned targ);
	void y(unsigned targ);
	void z(unsigned targ);
	void u1(unsigned targ, const std::vector<F> & matrix);
	void u1(unsigned targ, const std::vector<std::complex<F>> & matrix);

	// Controlled gates
	void crx(unsigned ctrl, unsigned targ, F angle);
	void cry(unsigned ctrl, unsigned targ, F angle);
	void crz(unsigned ctrl, unsigned targ, F angle);  
	void cphase(unsigned ctrl, unsigned targ, F angle); 
	void ch(unsigned ctrl, unsigned targ);
	void cnot(unsigned ctrl, unsigned targ);
	void cy(unsigned ctrl, unsigned targ);
	void cz(unsigned ctrl, unsigned targ);
	void cu1(unsigned ctrl, unsigned targ, const std::vector<F> & matrix);
	void cu1(unsigned ctrl, unsigned targ, const std::vector<std::complex<F>> & matrix);

	// Rest of the number gates
	void nrx(unsigned targ1, unsigned targ2, F angle);
	void nry(unsigned targ1, unsigned targ2, F angle);
	void nrz(unsigned targ1, unsigned targ2, F angle);
	void swap(unsigned targ1, unsigned targ2);
	void fswap(unsigned targ1, unsigned targ2);
	void iswap(unsigned targ1, unsigned targ2);
	void nh(unsigned targ1, unsigned targ2);
	void nu1(unsigned targ1, unsigned targ2, const std::vector<F> & matrix);
	void nu1(unsigned targ1, unsigned targ2, const std::vector<std::complex<F>> & matrix);

	/// Arbitrary two-qubit unitary
	void u2(unsigned targ1, unsigned targ2, const std::vector<F> & matrix);
	void u2(unsigned targ1, unsigned targ2, const std::vector<std::complex<F>> & matrix);

	/// Calculate probability of specific outcome of a qubit
	F prob(unsigned targ, unsigned outcome) const;
	/// Measure one qubit
	unsigned measure(unsigned targ, std::uniform_random_bit_generator auto g = gen);
	/// Measure all of the qubits
	std::size_t measure_all(std::uniform_random_bit_generator auto g = gen);
	/// Postselect one qubit on a specific outcome
	unsigned postselect(unsigned targ, unsigned outcome) const;
	/// Perform measurement sampling on one qubit
	std::vector<std::size_t> sample(unsigned targ, std::size_t samples,
					std::uniform_random_bit_generator auto g = gen) const;
	/// Perform measurement sampling on all the qubits
	std::map<std::size_t, std::size_t> sample_all(std::size_t samples,
						      std::uniform_random_bit_generator auto g = gen) const;

	/// Measure one qubit and remove from state vector
	unsigned measure_out(unsigned targ, std::uniform_random_bit_generator auto g = gen);
	/// Postselect one qubit on a specific outcome and remove from state vector
	unsigned postselect_out(unsigned targ, unsigned outcome) const;
    };


    /**
     * \brief Fixed number quantum simulator.
     */
    template<std::floating_point F, bool D = false, typename P = opt>
    class number
    {
    public:
	/// Instantiate a simulator with num_qubits qubits in the all zero state
	explicit number(unsigned num_qubits);
	/// Instantiate a simulator based on the state that is passed in 
	explicit number(const std::vector<std::complex<F>> & state);
	/// Instantiate a simulator with a specific number of ones
	number(unsigned num_qubits, unsigned num_ones);

	/// Convert from any simulator of the same type
	template<has_state_vector_of_type<F> S>
	explicit number(const S & s);

	/// Allow explicit cast between different simulator types
	template<std::floating_point F1, bool D1, typename P1>
	explicit operator number<F1,D1,P1>();

	/// Get the number of qubits
	unsigned size() const;
	/// Get the dimension of the Hilbert space
	unsigned dim() const;
	/// Return std::vector of state
	std::vector<std::complex<F>> get_state() const;
	/// Set a state
	void set_state(const std::vector<std::complex<F>> & state);
	/// Reset to the all zero computational basis state
	void reset();
	/// Access state vector elements (read-only). 
	const std::complex<F> & operator[](std::size_t index) const;
	/// Generate random state, need to pass in a generator which will default to our global one.
	void make_random(std::uniform_random_bit_generator auto & g = gen)
	    {
		// In order for this to work, the body of this function must
		// be in the header file (if the user is allowed to pass custom
		// gen types)
		g(); // Just to surpress unused-variable warning for now
	    }


	/// Get the number of ones
	unsigned get_ones() const;
	/// Set the number of ones - sets the lowest possible computational basis state
	void set_ones(unsigned num_ones);
	
	/// Print state vector
	void print(std::ostream & os = std::cout) const;
	/// Save state to json
	void save_json(const std::filesystem::path & path) const;
	/// Load state from json
	void load_json(const std::filesystem::path & path);
	/// \todo Select a compressed format
	/// Save state to compressed format
	void save_hdf5(const std::filesystem::path & path) const;
	/// Load state from compressed format
	void load_hdf5(const std::filesystem::path & path);

	// One-qubit gates
	void rz(unsigned targ, F angle);  
	void phase(unsigned targ, F angle); 
	void z(unsigned targ);  

	// Controlled gates
	void crz(unsigned ctrl, unsigned targ, F angle);  
	void cphase(unsigned ctrl, unsigned targ, F angle); 
	void cz(unsigned ctrl, unsigned targ);  

	// Rest of the number gates
	void nrx(unsigned targ1, unsigned targ2, F angle);
	void nry(unsigned targ1, unsigned targ2, F angle);
	void nrz(unsigned targ1, unsigned targ2, F angle);
	void swap(unsigned targ1, unsigned targ2);
	void fswap(unsigned targ1, unsigned targ2);
	void iswap(unsigned targ1, unsigned targ2);
	void nh(unsigned targ1, unsigned targ2);
	void nu1(unsigned targ1, unsigned targ2, const std::vector<F> & matrix);
	void nu1(unsigned targ1, unsigned targ2, const std::vector<std::complex<F>> & matrix);

	/// Calculate probability of specific outcome of a qubit
	F prob(unsigned targ, unsigned outcome) const;
	/// Measure one qubit
	unsigned measure(unsigned targ, std::uniform_random_bit_generator auto g = gen);
	/// Measure all of the qubits
	std::size_t measure_all(std::uniform_random_bit_generator auto g = gen);
	/// Postselect one qubit on a specific outcome
	unsigned postselect(unsigned targ, unsigned outcome) const;
	/// Perform measurement sampling on one qubit
	std::vector<std::size_t> sample(unsigned targ, std::size_t samples,
					std::uniform_random_bit_generator auto g = gen) const;
	/// Perform measurement sampling on all the qubits
	std::map<std::size_t, std::size_t> sample_all(std::size_t samples,
						      std::uniform_random_bit_generator auto g = gen) const;

    };

    /// TODO concept for sim
    /// TODO double check where this should go -- in or out the namespace
    template<typename S>
    std::ostream & operator << (std::ostream & os, const S & s);
    
    /// Calculate the Fubini-Study metric between two simulators/vectors
    template<has_state_vector S1, has_state_vector S2>
    requires have_same_precision<S1, S2> 
    get_precision_t<S1> distance(const S1 & s1, const S2 & s2);

    template<std::floating_point F, has_state_vector_of_type<F> S>
    F distance(const S & s1, const std::vector<std::complex<F>> & s2);

    template<std::floating_point F, has_state_vector_of_type<F> S>
    F distance(const std::vector<std::complex<F>> & s1, const S & s2);

    /// Calculate the fidelity between two simulators/vectors
    template<has_state_vector S1, has_state_vector S2>
    requires have_same_precision<S1, S2> 
    get_precision_t<S1> fidelity(const S1 & s1, const S2 & s2);

    template<std::floating_point F, has_state_vector_of_type<F> S>
    F fidelity(const S & s1, const std::vector<std::complex<F>> & s2);

    template<std::floating_point F, has_state_vector_of_type<F> S>
    F fidelity(const std::vector<std::complex<F>> & s1, const S & s2);

    /// Calculate the inner product between two simulators/vectors - do not normalise vector?
    template<has_state_vector S1, has_state_vector S2>
    requires have_same_precision<S1, S2> 
    get_precision_t<S1> inner_prod(const S1 & s1, const S2 & s2);

    template<std::floating_point F, has_state_vector_of_type<F> S>
    F inner_prod(const S & s1, const std::vector<std::complex<F>> & s2);

    template<std::floating_point F, has_state_vector_of_type<F> S>
    F inner_prod(const std::vector<std::complex<F>> & s1, const S & s2);

    

}

