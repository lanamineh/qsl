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
    /// Run functions without parallelisation
    struct seq;
    /// Run functions with omp
    struct omp;
    /// Automatically turn omp on and off depending on number of qubits
    struct opt;

    /**
     * \brief Restricted to the types that specify a parallelisation level
     *
     * This concept restricts a template parameter to the types that are
     * valid parallelisation levels.
     *
     */
    template<typename T>
    concept parallelisation = std::is_same_v<T,seq> || std::is_same_v<T,omp>
	|| std::is_same_v<T,opt>;

    /// A general type is not real or complex
    template<typename T>
    struct real_or_complex_t : std::false_type {};

    /// Specialisation for real floating point
    template<std::floating_point F>
    struct real_or_complex_t<F&> : std::true_type {};

    /// Specialisation for complex floating point
    template<std::floating_point F>
    struct real_or_complex_t<std::complex<F>&> : std::true_type {};

    /**
     * \brief Concept to check whether a type is a real or complex number
     *
     * This concept is true for a type which is either a reference to a 
     * built-in floating point type (a real number), or a reference to a 
     * std::complex of a floating point type (a complex number). 
     * Otherwise it is false. 
     */ 
    template<typename T>
    concept real_or_complex = real_or_complex_t<T>::value;

    template<typename T>
    concept state_vector = requires (T t) {

	// Operator[] must be valid and return a real or complex number
	// TODO Fix this.
	//{t[0]} -> real_or_complex;

	
	// Must return its size like std::vector
	{t.size()} -> std::same_as<std::size_t>; 
    };
    
    /**
     * \brief Struct for obtaining the precision of a simulator or state vector
     */
    template<typename T>
    struct get_precision;

    /**
     * \brief Obtain the precision of a simulator object
     *
     * The resulting precision is stored in the type member.
     */
    template<template<typename,bool,typename> typename S,
	     std::floating_point F,
	     bool D,
	     parallelisation P>
    struct get_precision<S<F, D, P>>
    {
	using type = F;
    };

    /**
     * \brief Obtain the precision for a complex state vector
     * 
     * Precision is stored in the type member.
     */
    template<std::floating_point F>
    struct get_precision<std::vector<std::complex<F>>>
    {
    	using type = F;
    };

    /**
     * \brief Obtain the precision for a real state vector
     * 
     * Precision is stored in the type member.
     */
    template<std::floating_point F>
    struct get_precision<std::vector<F>>
    {
    	using type = F;
    };

    
    /**
     * \brief Allow get_precision to work with built-in floating-point types too
     */
    template<std::floating_point F>
    struct get_precision<F>
    {
	using type = F;
    };
    
    /**
     * \brief Helper to get the precision type of a simulator or state directly
     *
     * \tparam T The simulator or state vector whose precision is wanted.
     */
    template<typename T>
    using get_precision_t = typename get_precision<T>::type;

    /**
     * \brief Check if two types have the same precision
     *
     * This concept compares the precision of two types, which may be
     * simulators, state vectors (real or complex), or built-in floating
     * point types. The concept returns true if the types have the same
     * floating point type. For a simulator, the floating point type is
     * the floating point type of its internal complex state vector. For
     * a real or complex state vector, it is the floating point type
     * of the underlying vector.
     *
     * \tparam T The first type to check
     * \tparam U Another simulator/state vector/built-in type to compare with T
     */
    template<typename T, typename U>
    concept same_precision = std::is_same_v<get_precision_t<T>,
					    get_precision_t<U>>;

    
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
     *
     * \tparam F The floating point precision.
     * \tparam D Set debugging to be on (true) or off (false). Defaults to off
     *           to prioritise speed. 
     * \tparam P Set parallelisation level. Available options are off (qsl::seq), 
     *           on (qsl::omp), or an automatic optimal selection (qsl::opt).
     */
    template<std::floating_point F, bool D = false, parallelisation P = opt>
    class basic
    {
    public:
	/**
	 * \brief Initialise the class with a specified number of qubits.
	 *
	 * The simulator is initialised in the all-zero state.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 */ 
	explicit basic(unsigned num_qubits);

	/**
	 * \brief Instantiate a simulator based on a real state that is passed in. 
         *        The imaginary part of the state vector will be assumed to be all-zero.
	 *
	 * The input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. 
	 *
	 * \param state A (real) vector containing the initial state for the object.
	 */ 	
	explicit basic(const std::vector<F> & state);

	/**
	 * \brief Instantiate a simulator based on the state that is passed in. 
	 *
	 * The input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step.
	 *
	 * \param state A complex vector containing the initial state for the object.
	 */ 	
	explicit basic(const std::vector<std::complex<F>> & state);

	/**
	 * \brief Convert from any simulator using the same floating point type.
	 *
	 * This is not a copy constructor (because it is templated), so it will not
	 * cause the move constructors to be implicitly deleted.
	 *
	 * \param s A valid qsl simulator object that uses the same floating point type.
	 */
	template<state_vector S>
	explicit basic(const S & s);

	/**
	 * \brief Convert between different floating point types for qsl::basic simulators.
	 */
	template<std::floating_point F1, bool D1, parallelisation P1>
	explicit operator basic<F1,D1,P1>();
	
	/**
	 * \brief Get the number of qubits.
	 *
	 * \return Number of qubits.
	 */
	unsigned qubits() const;

	/**
	 * \brief Get the dimension of the Hilbert space = 2 ^ num_qubits.
	 *
	 * \return Dimension of Hilbert space.
	 */
	std::size_t size() const;

	/**
	 * \brief Return the current state of the qubits.
	 *
	 * \return The state of the qubits as a std::vector.
	 */
	std::vector<std::complex<F>> get_state() const;

	/**
	 * \brief Change the simulator state to the real state that is passed in. 
         *        The imaginary part of the state vector will be assumed to be all-zero.
	 *
	 * The input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. Note that this function
	 * allows for a change in the number of qubits.
	 *
	 * \param state A (real) vector containing the new state for the object.
	 */ 		
	void set_state(const std::vector<F> & state);

	basic & operator= (const std::vector<F> & state);
	
	/**
	 * \brief Change the simulator state to the state that is passed in. 
	 *
	 * The input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. Note that this function
	 * allows for a change in the number of qubits.
	 *
	 * \param state A vector containing the new state for the object.
	 */ 		
	void set_state(const std::vector<std::complex<F>> & state);

	/**
	 * \brief Reset to the all-zero computational basis state.
	 */
	void reset();

	/**
	 * \brief Change number of qubits and reset to the all-zero computational basis state.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 */
	void reset(unsigned num_qubits);
	
	/**
	 * \brief Access state vector elements (read-only).
	 *
	 * This is read-only to avoid accidental tampering with the state vector
	 * that might render the state invalid (for example, setting every element to zero).  
	 *
	 * \param index The state vector index to access.
	 * \return The complex amplitude at index.
	 */ 
	const std::complex<F> & operator[](std::size_t index) const;
	
	/**
	 * \brief Set the state vector to a random state.
	 *
	 * This function allows the user to pass in a specific random number generator,
	 * this is to allow for seeding and reproducibility.
	 *
	 * \param g Optional, uniform random number generator.
	 */ 
	void make_random(std::uniform_random_bit_generator auto & g = gen)
	    {
		// In order for this to work, the body of this function must
		// be in the header file (if the user is allowed to pass custom
		// gen types)
		g(); // Just to surpress unused-variable warning for now
	    }
	
	/**
	 * \brief Pretty print the state vector to an output stream.
	 *
	 * Prints the entire state vector up to 6 qubits (64 elements). Beyond
	 * this, it prints the first and last 20 elements.
	 *
	 * \param os Optional, output stream, defaults to printing to the console with std::cout. 
	 */
	void print(std::ostream & os = std::cout) const;

	/**
	 * \brief Save the state vector and associated metadata to a json file.
	 *
	 * Capped at 20 qubits (16MB at double precision to store the state vector).
	 * \todo Specify which metadata are stored.
	 *
	 * \param path Filename to save to.
	 */
	void save_json(const std::filesystem::path & path) const;

	/**
	 * \brief Load in a state vector from a json file.
	 *
	 * The input file must contain a vector that can be mapped to a valid state vector
	 * (length a power of two, non-zero, does not need to be normalised) in the 'state'
	 * field. Capped at 20 qubits.
	 *
	 * \todo Figure out how we will serialise std::vector<std::complex> so the input
	 *       can be specified here.
	 *
	 * \param path The filename to read from. 
	 */
	void load_json(const std::filesystem::path & path);

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
     * \brief Resizeable quantum simulator. 
     *
     * Unlike qsl::basic, the number of qubits in the simulator is allowed to change throughout
     * the computation. This is useful for computations involving ancillae.
     *
     * \tparam F The floating point precision.
     * \tparam D Set debugging to be on (true) or off (false). Defaults to off
     *           to prioritise speed. 
     * \tparam P Set paralellisation level. Available options are off (qsl::seq), 
     *           on (qsl::omp), or an automatic optimal selection (qsl::opt).
     */
    template<std::floating_point F, bool D = false, parallelisation P = opt>
    class resize
    {
    public:
	/**
	 * \brief Initialise the class with a specified number of qubits.
	 *
	 * The simulator is initialised in the all-zero state.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 */ 
	explicit resize(unsigned num_qubits);

	/**
	 * \brief Instantiate a simulator based on a real state that is passed in. 
         *        The imaginary part of the state vector will be assumed to be all-zero.
	 *
	 * The input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. 
	 *
	 * \param state A (real) vector containing the initial state for the object.
	 */ 	
	explicit resize(const std::vector<F> & state);

	/**
	 * \brief Instantiate a simulator based on the state that is passed in. 
	 *
	 * The input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step.
	 *
	 * \param state A complex vector containing the initial state for the object.
	 */ 	
	explicit resize(const std::vector<std::complex<F>> & state);

	/**
	 * \brief Convert from any simulator using the same floating point type.
	 *
	 * This is not a copy constructor (because it is templated), so it will not
	 * cause the move constructors to be implicitly deleted.
	 *
	 * \param s A valid qsl simulator object that uses the same floating point type.
	 */
	template<state_vector S>
	explicit resize(const S & s);

	/**
	 * \brief Convert between different floating point types for qsl::resize simulators.
	 */
	template<std::floating_point F1, bool D1, parallelisation P1>
	explicit operator resize<F1,D1,P1>();
	
	/**
	 * \brief Get the current number of qubits.
	 *
	 * \return Number of qubits.
	 */
	unsigned qubits() const;
	
	/**
	 * \brief Get the current dimension of the Hilbert space = 2 ^ num_qubits.
	 *
	 * \return Dimension of Hilbert space.
	 */
	std::size_t size() const;
	
	/**
	 * \brief Return the current state of the qubits.
	 *
	 * \return The state of the qubits as a std::vector.
	 */
	std::vector<std::complex<F>> get_state() const;

	/**
	 * \brief Change the simulator state to the real state that is passed in. 
         *        The imaginary part of the state vector will be assumed to be all-zero.
	 *
	 * The input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. 
	 *
	 * \param state A (real) vector containing the new state for the object.
	 */ 		
	void set_state(const std::vector<F> & state);

	/**
	 * \brief Change the simulator state to the state that is passed in. 
	 *
	 * The input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. 
	 *
	 * \param state A vector containing the new state for the object.
	 */ 		
	void set_state(const std::vector<std::complex<F>> & state);

	/**
	 * \brief Reset to the all-zero computational basis state.
	 */
	void reset();

	/**
	 * \brief Change number of qubits and reset to the all-zero computational basis state.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 */
	void reset(unsigned num_qubits);

	/**
	 * \brief Access state vector elements (read-only).
	 *
	 * This is read-only to avoid accidental tampering with the state vector
	 * that might render the state invalid (for example, setting every element to zero).  
	 *
	 * \param index The state vector index to access.
	 * \return The complex amplitude at index.
	 */ 
	const std::complex<F> & operator[](std::size_t index) const;

	/**
	 * \brief Set the state vector to a random state.
	 *
	 * This function allows the user to pass in a specific random number generator,
	 * this is to allow for seeding and reproducibility.
	 *
	 * \param g Optional, uniform random number generator.
	 */ 	
	void make_random(std::uniform_random_bit_generator auto & g = gen)
	    {
		// In order for this to work, the body of this function must
		// be in the header file (if the user is allowed to pass custom
		// gen types)
		g(); // Just to surpress unused-variable warning for now
	    }

	/**
	 * \brief Add a qubit in the zero state to the end of the state vector.
	 */
	void add_qubit();

	/**
	 * \brief Insert a qubit in the zero state at position targ in the state vector.
	 *
	 * Qubits that were in positions above targ will be shifted along i.e. the previous 
	 * qubit that was at position targ will now be indexed as targ+1.
	 *
	 * \param targ The value the inserted qubit will be indexed as.
	 */
	void add_qubit(unsigned targ);
	/// Trims the state vector to the current qubit size
	/**
	 * \brief Trims the internal state vector to the current qubit size.
	 *
	 * This function is useful if memory is an issue. For efficiency, the internal
	 * state vector remains as big as the largest number of qubits it has been
	 * required to simulate previously (but these excess elements are not used in
	 * any computation). This function trims the internal state vector to the same
	 * size as the current dimension.  
	 */
	void trim();
	
	/**
	 * \brief Pretty print the state vector to an output stream.
	 *
	 * Prints the entire state vector up to 6 qubits (64 elements). Beyond
	 * this, it prints the first and last 20 elements.
	 *
	 * \param os Optional, output stream, defaults to printing to the console with std::cout. 
	 */
	void print(std::ostream & os = std::cout) const;

	/**
	 * \brief Save the state vector and associated metadata to a json file.
	 *
	 * Capped at 20 qubits (16MB at double precision to store the state vector).
	 * \todo Specify which metadata are stored.
	 *
	 * \param path Filename to save to.
	 */
	void save_json(const std::filesystem::path & path) const;

	/**
	 * \brief Load in a state vector from a json file.
	 *
	 * The input file must contain a vector that can be mapped to a valid state vector
	 * (length a power of two, non-zero, does not need to be normalised) in the 'state'
	 * field. Capped at 20 qubits.
	 *
	 * \todo Figure out how we will serialise std::vector<std::complex> so the input
	 *       can be specified here.
	 *
	 * \param path The filename to read from. 
	 */
	void load_json(const std::filesystem::path & path);
	
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
     *
     * This simulator stores state vectors where the only non-zero amplitudes are
     * associated to bitstrings with a specified number of ones. This can speed up
     * computations involving quantum chemistry.
     *
     * \tparam F The floating point precision.
     * \tparam D Set debugging to be on (true) or off (false). Defaults to off
     *           to prioritise speed. 
     * \tparam P Set paralellisation level. Available options are off (qsl::seq), 
     *           on (qsl::omp), or an automatic optimal selection (qsl::opt).
     */
    template<std::floating_point F, bool D = false, parallelisation P = opt>
    class number
    {
    public:
	/**
	 * \brief Initialise the class with a specified number of qubits.
	 *
	 * The simulator is initialised in the all-zero state, and the number of ones
	 * is set to zero.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 */ 
	explicit number(unsigned num_qubits);

	/**
	 * \brief Instantiate a simulator based on a real state that is passed in. 
         *        The imaginary part of the state vector will be assumed to be all-zero.
	 *
	 * The input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. The input vector must
	 * be of fixed number.
	 *
	 * \param state A (real) vector containing the initial state for the object.
	 */ 	
	explicit number(const std::vector<F> & state);

	/**
	 * \brief Instantiate a simulator based on the state that is passed in. 
	 *
	 * The input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. The input vector must
	 * be of fixed number.
	 *
	 * \param state A complex vector containing the initial state for the object.
	 */ 	
	explicit number(const std::vector<std::complex<F>> & state);

	/**
	 * \brief Initialise the class with a specified number of qubits and number of ones.
	 *
	 * The simulator is initialised in the lowest indexed computational basis state
	 * with num_ones ones.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 * \param num_ones The number of ones in the fixed number simulator.
	 */ 
	number(unsigned num_qubits, unsigned num_ones);

	/**
	 * \brief Convert from any simulator using the same floating point type.
	 *
	 * This is not a copy constructor (because it is templated), so it will not
	 * cause the move constructors to be implicitly deleted.
	 *
	 * \param s A valid qsl simulator object that uses the same floating point type.
	 */
	template<state_vector S>
	explicit number(const S & s);

	/**
	 * \brief Convert between different floating point types for qsl::number simulators.
	 */
	template<std::floating_point F1, bool D1, parallelisation P1>
	explicit operator number<F1,D1,P1>();

	/**
	 * \brief Get the number of qubits.
	 *
	 * \return Number of qubits.
	 */
	unsigned qubits() const;

	/**
	 * \brief Get the dimension of the Hilbert space = 2 ^ num_qubits.
	 *
	 * \return Dimension of Hilbert space.
	 */	
	std::size_t size() const;

	/**
	 * \brief Return the current state of the qubits.
	 *
	 * \return The state of the qubits as a std::vector.
	 */
	std::vector<std::complex<F>> get_state() const;

	/**
	 * \brief Change the simulator state to the real state that is passed in. 
         *        The imaginary part of the state vector will be assumed to be all-zero.
	 *
	 * The input state vector must be non-zero, fixed number and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. Note that this function
	 * allows for a change in the number of qubits.
	 *
	 * \param state A (real) vector containing the new state for the object.
	 */ 		
	void set_state(const std::vector<F> & state);

	/**
	 * \brief Change the simulator state to the state that is passed in. 
	 *
	 * The input state vector must be non-zero, fixed number and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. Note that this function
	 * allows for a change in the number of qubits.
	 *
	 * \param state A vector containing the new state for the object.
	 */ 		
	void set_state(const std::vector<std::complex<F>> & state);

	/**
	 * \brief Reset to the lowest indexed computational basis state for the given num_ones.
	 */
	void reset();

	/**
	 * \brief Change number of qubits and reset to the all-zero computational basis state.
	 *
	 * The number of ones is also reset to zero.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 */
	void reset(unsigned num_qubits);

	/**
	 * \brief Change number of qubits number of ones.
	 *
	 * The state vector is set to the lowest computational basis
	 * state with the specified number of ones.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 * \param num_qubits The number of ones to set.
	 */
	void reset(unsigned num_qubits, unsigned num_ones);
	    
	/**
	 * \brief Access state vector elements (read-only).
	 *
	 * This is read-only to avoid accidental tampering with the state vector
	 * that might render the state invalid (for example, setting every element to zero).  
	 *
	 * \param index The state vector index to access.
	 * \return The complex amplitude at index.
	 */ 
	const std::complex<F> & operator[](std::size_t index) const;

	/**
	 * \brief Set the state vector to a random state.
	 *
	 * This function allows the user to pass in a specific random number generator,
	 * this is to allow for seeding and reproducibility.
	 *
	 * \param g Optional, uniform random number generator.
	 */ 	
	void make_random(std::uniform_random_bit_generator auto & g = gen)
	    {
		// In order for this to work, the body of this function must
		// be in the header file (if the user is allowed to pass custom
		// gen types)
		g(); // Just to surpress unused-variable warning for now
	    }

	/**
	 * \brief Get the current number of ones.
	 *
	 * \return Number of ones.
	 */
	unsigned ones() const;

	/**
	 * \brief Pretty print the state vector to an output stream.
	 *
	 * Prints the entire state vector up to 6 qubits (64 elements). Beyond
	 * this, it prints the first and last 20 elements.
	 *
	 * \param os Optional, output stream, defaults to printing to the console with std::cout. 
	 */
	void print(std::ostream & os = std::cout) const;
	
	/**
	 * \brief Save the state vector and associated metadata to a json file.
	 *
	 * Capped at 20 qubits (16MB at double precision to store the state vector).
	 * \todo Specify which metadata are stored.
	 *
	 * \param path Filename to save to.
	 */	
	void save_json(const std::filesystem::path & path) const;

	/**
	 * \brief Load in a state vector from a json file.
	 *
	 * The input file must contain a vector that can be mapped to a valid state vector
	 * (length a power of two, non-zero, does not need to be normalised) in the 'state'
	 * field. Capped at 20 qubits.
	 *
	 * \todo Figure out how we will serialise std::vector<std::complex> so the input
	 *       can be specified here.
	 *
	 * \param path The filename to read from. 
	 */
	void load_json(const std::filesystem::path & path);
	
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
    template<state_vector S1, state_vector S2>
    requires same_precision<S1, S2> 
    get_precision_t<S1> distance(const S1 & s1, const S2 & s2);

    /// Calculate the fidelity between two simulators/vectors
    template<state_vector S1, state_vector S2>
    requires same_precision<S1, S2> 
    get_precision_t<S1> fidelity(const S1 & s1, const S2 & s2);

    /// Calculate the inner product between two simulators/vectors
    ///- do not normalise vector?
    template<state_vector S1, state_vector S2>
    requires same_precision<S1, S2> 
    get_precision_t<S1> inner_prod(const S1 & s1, const S2 & s2);    

}

