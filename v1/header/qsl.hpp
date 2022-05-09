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

/// The namespace for all items in QSL
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

    /**
     * \brief Struct for obtaining the precision of a simulator or state vector
     */
    template<typename T>
    struct get_precision;

    /// Get precision of real
    template<std::floating_point F>
    struct get_precision<F>
    {
	using type = F;
    };

    /// Get precision of complex
    template<std::floating_point F>
    struct get_precision<std::complex<F>>
    {
	using type = F;
    };

    /// By default, nothing is a complex type
    template<typename T>
    struct is_complex : std::false_type {};

    /// Only a std::complex<T> for floating-point T is complex 
    template<std::floating_point F>
    struct is_complex<std::complex<F>> : std::true_type {};

    /**
     * \brief Concept to check if a type is real or complex
     *
     * This concept is true for floating-point real or complex types. A real type is
     * a built-in floating-point type. A complex type is a std::complex<F> where F is
     * a real type. The concept is also true for any reference or const reference to
     * a real or complex type.
     */
    template<typename T>
    concept real_or_complex = std::floating_point<std::remove_cvref_t<T>>
	|| is_complex<std::remove_cvref_t<T>>::value;
    
    template<typename T>
    concept state_vector = requires (T t) {

	// The type should have a value_type public typedef like std::vector and
	// other containers.
	typename T::value_type;
	    
	///\todo Check whether the input type to operator[] should be restricted
	// Operator[] must be valid and return a real or complex number
	{t[0]} -> real_or_complex;

	///\todo Check whether using convertible_to is appropriate
	// Must return its size like std::vector
	{t.size()} -> std::unsigned_integral; 
    };
    
    template<typename S>
    concept debug_state_vector = state_vector<S> && requires (S s) {
	{s.debug()} -> std::same_as<bool>;
    };
    
    template<state_vector S>
    struct get_precision<S>
    {
	using type = get_precision<typename S::value_type>::type;
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
     * \brief Random number generator complying with std::uniform_random_bit_generator.
     *
     * This is a wrapper around std::mt19937_64.
     *
     * Testing:
     * - Set a seed, read the seed, check they match (both setting with constructor and seed function).
     * - Check operator() by generating some random values and checking they are between min() and max()
     * - Check operator() is uniform, for example, generate lots of random numbers, check 
     *   in python whether they are uniform, then check if the same numbers are regenerated with the seed.
     * - Check seeding works
     */
    class gen_t
    {
    public:    
	using result_type = std::mt19937_64::result_type;
	/// Return minimum possible value
	static result_type min();
	/// Return maximum possible value
	static result_type max();
	/// Obtain the next sample from the generator
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
     * \brief Turn logging of debugging information on/off 
     *
     * When a simulator is created that has debug on (D = true), this function
     * sets whether a logger which prints debugging information is activated or not.
     * Note that this function switches logging on/off globally. 
     *
     * QSL simulators will only record logging information if debugging is enabled
     * (template parameter D = true in the simulator). If D = false, the state of
     * this global logger is ignored.
     *
     * \param activate Switch logging on (true) or off (false).
     * \param os Stream to output the logging data to, defaults to std::cout.
     *
     * Testing:
     * - Enable logging using a dummy stream (such as a string stream), and then
     *   execute one of each kind of operation that should print to the log. This
     *   will test all the simulator logging functions. Compare the string stream
     *   with correct messages.
     * - Disable logging (in the same context as above, with a stream stream), and
     *   check that operations do not write to the stream.
     * 
     */
    void log(bool activate, std::ostream & os = std::cout);
    
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
	using value_type = std::complex<F>;
	constexpr bool debug() { return D; }

	/**
	 * \brief Construct an empty simulator
	 *
	 * Constructs a simulator with zero qubits. The number of qubits
	 * can be changed using the qsl::basic::reset member function. This 
	 * function is not especially useful for one or two simulators, but
	 * may help in contexts that require construction-before-copy (for 
	 * example, adding a simulator into a map without using emplace).
	 *
	 */
	basic();
	
	/**
	 * \brief Initialise the class with a specified number of qubits.
	 *
	 * The simulator is initialised in the all-zero state.
	 *
	 * In debug mode, if the number of qubits is too large to simulate, a
	 * std::runtime_error is thrown.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 *
	 * Testing:
	 * - In debug mode, check exceptions thrown. E.g. try inputting a large number,
	 *   negative numbers etc.
	 * - Check state vector is in the all-zero state and the correct size.
	 */ 
	explicit basic(unsigned num_qubits);

	/**
	 * \brief Instantiate a simulator based on a std::vector or another simulator
         *        that has the same floating point precision.
	 *
	 * If inputting a std::vector, it must be non-zero and have a length which is a 
	 * power of two. It does not need to be normalised as this function will carry
	 * out a normalisation step.
	 *
	 * In debug mode, if the input state is not a power of two or is all zeros,
	 * a std::invalid_argument is thrown.
	 *
	 * This is not a copy constructor (because it is templated), so it will not
	 * cause the move constructors to be implicitly deleted.
	 *
	 * \param s A valid qsl simulator object with the same floating point precision,
         *          or std::vector<std::complex<F>> or std::vector<F>.
	 *
	 * Testing:
	 * - In debug mode, check exceptions thrown. E.g. input wrong length state vector,
	 *   or one of correct length and all zeros.
	 * - Input un-normalised state and check it becomes normalised.
	 * - Input normalised state and check it stays the same.
	 * - Try all possible input types (std::vectors and simulators) and check state is correct.
	 * - Check number of qubits is correctly calculated.
	 */
	template<state_vector S>
	requires same_precision<F, S>
	explicit basic(const S & s);

	/**
	 * \brief Convert between different floating point types for qsl::basic simulators.
	 *
	 * Testing:
	 * - Check state vector is correct and normalised (this function could introduce small
	 *   floating point errors). 
	 */
	template<std::floating_point F1, bool D1, parallelisation P1>
	explicit operator basic<F1,D1,P1>();
	
	/**
	 * \brief Get the number of qubits.
	 *
	 * \return Number of qubits.
	 *
	 * Testing:
	 * - Instantiate a simulator with a number of qubits and check this function 
	 *   returns the correct number.
	 */
	unsigned qubits() const;

	/**
	 * \brief Get the dimension of the Hilbert space = 2 ^ num_qubits.
	 *
	 * \return Dimension of Hilbert space.
	 *
	 * Testing:
	 * - Instantiate a simulator with a number of qubits and check this function 
	 *   returns the correct state vector size.
	 */
	std::size_t size() const;

	/**
	 * \brief Change the simulator state to the state that is passed in. 
	 *
	 * If a std::vector, the input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. Note that this function
	 * allows for a change in the number of qubits.
	 *
	 * In debug mode, if the input state is not a power of two or is all zeros,
	 * a std::invalid_argument is thrown.
	 *
	 * \param state A vector or qsl simulator containing the new state for the object.
	 *
	 * Testing:
	 * - In debug mode, check exceptions thrown. E.g. input wrong length state vector,
	 *   or one of correct length and all zeros.
	 * - Input un-normalised state and check it becomes normalised.
	 * - Input normalised state and check it stays the same.
	 * - Try all possible input types (std::vectors and simulators) and check state is correct.
	 * - Check number of qubits is correctly calculated.
	 */ 		
	template<state_vector S>
	requires same_precision<F, S>
	basic & operator= (const S & state);

	/**
	 * \brief Reset to the all-zero computational basis state.
	 *
	 * Testing:
	 * - Check state after this function is called is the all-zero state.
	 */
	void reset();

	/**
	 * \brief Change number of qubits and reset to the all-zero computational basis state.
	 *
	 * In debug mode, if the number of qubits is too large to simulate, a
	 * std::runtime_error is thrown.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 *
	 * Testing:
	 * - In debug mode, check exceptions thrown. E.g. try inputting a large number,
	 *   negative numbers etc.
	 * - Check state vector is in the all-zero state and the correct size.
	 */
	void reset(unsigned num_qubits);

	/**
	 * \brief Return the current state of the qubits.
	 *
	 * \return The state of the qubits as a std::vector.
	 *
	 * Testing:
	 * - Instantiate a simulator with a specific state and check this function 
	 *   returns the same state (normalised).
	 */
	std::vector<std::complex<F>> state() const;
	
	/**
	 * \brief Access state vector elements (read-only).
	 *
	 * This is read-only to avoid accidental tampering with the state vector
	 * that might render the state invalid (for example, setting every element to zero).  
	 *
	 * A std::out_of_range error is thrown if index is bigger than size()-1.
	 *
	 * \param index The state vector index to access.
	 * \return The complex amplitude at index.
	 *
	 * Testing:
	 * - In debug mode, check exceptions thrown. E.g. try accessing out of bound elements
	 * - Check output of operator[] in a loop against the .state() method.
	 * - Check the state is normalised -- might be superfluous if we check against
	 *   state() (and check that is normalised).
	 */ 
	const std::complex<F> & operator[](std::size_t index) const;
	
	/**
	 * \brief Set the state vector to a random state.
	 *
	 * This function places the simulator in a random state, chosen uniformly from 
	 * the set of normalised states on the surface of the unit sphere in 
	 * \f$\mathcal{C}^n,\f$ where \f$n\f$ is the number of qubits in the simulator.
	 * (You can call qsl::basic::reset first if you need to change the number of
	 * qubits.)
	 *
	 * This function allows the user to pass in a specific random number generator,
	 * this is to allow for seeding and reproducibility. 
	 *
	 * \param g Optional, uniform random number generator.
	 *
	 * Testing:
	 * - Check that the result of several calls to make random is a normalised
	 *   state vector (e.g. using state()) each time. Pass several different generator 
	 *   types and check that the function produces normalised states in all cases.
	 * - Use a fixed (deterministic, seeded) generator, get the resulting random
	 *   vectors that are generated, and hard code them in the tests. Make this set
	 *   quite large, and verify (externally, using python), that the random states
	 *   are correctly distributed.
	 * - The above tests also check that the states are reproducible.
	 * 
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
	 *
	 * Testing:
	 * - Check that the function puts the correctly formatted string to a
	 *   stringstream. 
	 * - Check boundary cases like zero qubits.
	 */
	void print(std::ostream & os = std::cout) const;

	/**
	 * \brief Save the state vector and associated metadata to a json file.
	 *
	 * Capped at 20 qubits (16MB at double precision to store the state vector).
	 * \todo Specify which metadata are stored.
	 *
	 * A std::runtime_error is thrown if too many qubits are to be stored
	 * or if the file cannot be created. These exceptions are still thrown if
	 * the simulator is not in debug mode.
	 *
	 * \param file Filename to save to.
	 *
	 * Testing:
	 * - Check that exceptions are thrown for more than 20 qubits; if the
	 *   file cannot be created for any reason (maybe try to write to 
	 *   "/qsl.json" for the test, although this would work if they ran it
	 *   as root. Might need some thinking about. Maybe could create a directory
	 *   with no permissions as part of the test.
	 * - Check the right file is written, by creating it and then reading it
	 *   back and checking it has the right contents.
	 */
	void save_json(const std::filesystem::path & file) const;

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
	 * A std::runtime_error is thrown if the file doesn't exist or cannot be
	 * read. A std::invalid_argument is thrown if it does not contain a valid json object,
	 * if it reads a json object that does not contain a 'state' field, or if the state that
	 * is read is invalid. These exceptions occur even when debugging is disabled.
	 *
	 * \param file The filename to read from. 
	 *
	 * Testing:
	 * - Test by loading from fixed (hardcoded) JSON files in the test, and check
	 *   that the data is correct. Include a test for normalised and non-normalised
	 *   states.
	 * - Check exceptions (error if wrong fields, invalid json)
	 * - Check round trip save/load json gives the same state.
	 */
	void load_json(const std::filesystem::path & file);

	/**
	 * \brief Rotate around the x-axis of the Bloch sphere 
	 * \f$ e^{-i\theta X/2} \f$:
	 * 
	 * \f[ 
	 * R_x = \begin{pmatrix}
	 *       \cos(\theta/2) & -i\sin(\theta/2) \\
	 *       -i\sin(\theta/2) & \cos(\theta/2) \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 *
	 * Testing:
	 * - In debug mode, check exceptions thrown e.g. input invalid targ.
	 * - Apply gate to random state (smallish no. of qubits) and check output 
	 *   matches with armadillo. Do this with targ = 0 up to num_qubits-1 
	 *   to make sure edge cases have been checked. Pick random values for angle.
	 */
	void rx(unsigned targ, F angle);

	/**
	 * \brief Rotate around the y-axis of the Bloch sphere 
	 * \f$ e^{-i\theta Y/2} \f$:
	 * 
	 * \f[ 
	 * R_y = \begin{pmatrix}
	 *       \cos(\theta/2) & -\sin(\theta/2) \\
	 *       \sin(\theta/2) & \cos(\theta/2) \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 *
	 * Testing: same as rx gate.
	 */
	void ry(unsigned targ, F angle);

	/**
	 * \brief Rotate around the z-axis of the Bloch sphere 
	 * \f$ e^{-i\theta Z/2} \f$:
	 * 
	 * \f[ 
	 * R_z = \begin{pmatrix}
	 *       e^{-i\theta/2} & 0 \\
	 *       0 & e^{i\theta/2} \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 *
	 * Testing: same as rx gate.
	 */
	void rz(unsigned targ, F angle);  

	/**
	 * \brief Apply a phase shift to qubit targ:
	 *
	 * \f[
	 * \text{phase}(\theta) = \begin{pmatrix}
	 *       1 & 0 \\
	 *       0 & e^{i\theta} \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * Also equivalent to \f$ e^{i\theta/2} R_z(\theta) \f$.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 *
	 * Testing: same as rx gate.
	 */	
	void phase(unsigned targ, F angle); 

	/**
	 * \brief Apply the Hadamard gate to qubit targ:
	 *
	 * \f[ 
	 * H = \frac{1}{\sqrt{2}}\begin{pmatrix}
	 *     1 & 1 \\
	 *     1 & -1 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 *
	 * Testing: 
	 * - In debug mode, check exceptions thrown e.g. input invalid targ.
	 * - Apply gate to random state (smallish no. of qubits) and check output 
	 *   matches with armadillo. Do this with targ = 0 up to num_qubits-1 
	 *   to make sure edge cases have been checked. 
	 */
	void h(unsigned targ);

	/**
	 * \brief Apply the Pauli-X gate to qubit targ:
	 *
	 * \f[ 
	 * X = \begin{pmatrix}
	 *     0 & 1 \\
	 *     1 & 0 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * Also equivalent to \f$ iR_x(\pi) \f$.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 *
	 * Testing: same as h gate.
	 */
	void x(unsigned targ);

	/**
	 * \brief Apply the Pauli-Y gate to qubit targ:
	 *
	 * \f[ 
	 * Y = \begin{pmatrix}
	 *     0 & -i \\
	 *     i & 0 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * Also equivalent to \f$ iR_y(\pi) \f$.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 *
	 * Testing: same as h gate.
	 */
	void y(unsigned targ);

	/**
	 * \brief Apply the Pauli-Z gate to qubit targ:
	 *
	 * \f[ 
	 * Z = \begin{pmatrix}
	 *     1 & 0 \\
	 *     0 & -1 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * Also equivalent to \f$ iR_z(\pi) \f$ or \f$ \text{phase}(\pi) \f$
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 *
	 * Testing: same as h gate.
	 */
	void z(unsigned targ);

	/**
	 * \brief Apply an arbitrary one-qubit (real) unitary to targ.
	 *
	 * The matrix must have orthonormal columns, the columns will be
	 * normalised in this function.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 * A std::invalid_argument is thrown if matrix does not have orthonormal columns
	 * (i.e. cannot be normalised to a unitary matrix).
	 *
	 * \param targ The target qubit.
	 * \param matrix The unitary matrix to apply, in row-major form
	 *
	 * Testing: 
	 * - In debug mode, check exceptions thrown e.g. input invalid targ, put in 
	 *   an invalid matrix such as matrix with det = 0, non-orthogonal columns.
	 * - Apply gate to random state (smallish no. of qubits) and check output 
	 *   matches with armadillo. Do this with targ = 0 up to num_qubits-1 
	 *   to make sure edge cases have been checked. Input a few random unitary matrices.
	 */
	void u1(unsigned targ, const std::vector<F> & matrix);

	/**
	 * \brief Apply an arbitrary one-qubit unitary to targ.
	 *
	 * The matrix must have orthonormal columns, the columns will be
	 * normalised in this function. If all the coefficients are real, use the
	 * real version of u1 for improved performance.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 * A std::invalid_argument is thrown if matrix does not have orthonormal columns
	 * (i.e. cannot be normalised to a unitary matrix).
	 *
	 * \param targ The target qubit.
	 * \param matrix The unitary matrix to apply in row-major form.
	 *
	 * Testing: same as real version of u1.
	 */
	void u1(unsigned targ, const std::vector<std::complex<F>> & matrix);

	/**
	 * \brief Perform a controlled X-rotation on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, \f$ R_x \f$ is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 *
	 * Testing: 
	 * - In debug mode, check exceptions thrown e.g. input invalid targ and ctrl, set
	 *   targ = ctrl.
	 * - Apply gate to random state (smallish no. of qubits) and check output 
	 *   matches with armadillo. Do this with targ and ctrl  = 0 up to num_qubits-1 
	 *   to make sure edge cases have been checked. Pick random values for angle.
	 */
	void crx(unsigned ctrl, unsigned targ, F angle);

	/**
	 * \brief Perform a controlled Y-rotation on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, \f$ R_y \f$ is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 *
	 * Testing: same as crx
	 */
	void cry(unsigned ctrl, unsigned targ, F angle);

	/**
	 * \brief Perform a controlled Z-rotation on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, \f$ R_z \f$ is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 *
	 * Testing: same as crx
	 */
	void crz(unsigned ctrl, unsigned targ, F angle);  

	/**
	 * \brief Perform a controlled phase gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, a phase shift is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 *
	 * Testing: same as crx
	 */
	void cphase(unsigned ctrl, unsigned targ, F angle); 

	/**
	 * \brief Perform a controlled Hadamard gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, H is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 *
	 * Testing: 
	 * - In debug mode, check exceptions thrown e.g. input invalid targ and ctrl, set
	 *   targ = ctrl.
	 * - Apply gate to random state (smallish no. of qubits) and check output 
	 *   matches with armadillo. Do this with targ and ctrl  = 0 up to num_qubits-1 
	 *   to make sure edge cases have been checked. 
	 */
	void ch(unsigned ctrl, unsigned targ);

	/**
	 * \brief Perform a controlled Pauli-X (often referred to as a controlled-Not) 
	 *        gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, X is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 *
	 * Testing: same as ch
	 */
	void cnot(unsigned ctrl, unsigned targ);

	/**
	 * \brief Perform a controlled Pauli-Y gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, Y is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 *
	 * Testing: same as ch
	 */
	void cy(unsigned ctrl, unsigned targ);

	/**
	 * \brief Perform a controlled Pauli-Z gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, Z is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 *
	 * Testing: same as ch
	 */
	void cz(unsigned ctrl, unsigned targ);

	/**
	 * \brief Apply a controlled arbitrary one-qubit (real) unitary.
	 *
	 * The matrix must have orthonormal columns, the columns will be
	 * normalised in this function. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 * A std::invalid_argument is thrown if matrix does not have orthonormal columns
	 * (i.e. cannot be normalised to a unitary matrix).
	 *
	 * \param ctrl The control qubit, U is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param matrix The unitary matrix to apply, in row-major form.
	 *
	 * Testing: 
	 * - In debug mode, check exceptions thrown e.g. input invalid targ, put in 
	 *   an invalid matrix such as matrix with det = 0, non-orthogonal columns.
	 * - Apply gate to random state (smallish no. of qubits) and check output 
	 *   matches with armadillo. Do this with targ = 0 up to num_qubits-1 
	 *   to make sure edge cases have been checked. Input a few random unitary matrices.
	 */
	void cu1(unsigned ctrl, unsigned targ, const std::vector<F> & matrix);

	/**
	 * \brief Apply a controlled arbitrary one-qubit unitary.
	 *
	 * The matrix must have orthonormal columns, the columns will be
	 * normalised in this function. If all the coefficients are real, use the
	 * real version of cu1 for improved performance.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 * A std::invalid_argument is thrown if matrix does not have orthonormal columns
	 * (i.e. cannot be normalised to a unitary matrix).
	 *
	 * \param ctrl The control qubit, U is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param matrix The unitary matrix to apply, in row-major form.
	 *
	 * Testing: same as real version of cu1
	 */
	void cu1(unsigned ctrl, unsigned targ, const std::vector<std::complex<F>> & matrix);

	/**
	 * \brief Perform an X rotation on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace. 
	 *
	 * \f[ 
	 *  nR_x(\theta) = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & \cos(\theta/2) & -i\sin(\theta/2) & 0 \\
	 *             0 & -i\sin(\theta/2) & \cos(\theta/2) & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This is equivalent to applying \f$ e^{-i\theta (XX+YY)/2} \f$. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param angle The angle to rotate by (in radians)
	 *
	 * Testing: same as crx
	 */
	void nrx(unsigned targ1, unsigned targ2, F angle);

	/**
	 * \brief Perform a Y rotation on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace. This is equivalent to applying 
	 * \f$ e^{-i\theta (YX-XY)/2} \f$. 
	 *
	 * \ingroup qubits_gates

	 * \f[ 
	 *  nR_y(\theta) = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & \cos(\theta/2) & -\sin(\theta/2) & 0 \\
	 *             0 & \sin(\theta/2) & \cos(\theta/2) & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param angle Angle to rotate by (in radians).
	 *
	 * Testing: same as crx
	 */
	void nry(unsigned targ1, unsigned targ2, F angle);

	/**
	 * \brief Perform a Z rotation on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace. 
	 *
	 * \f[ 
	 *  nR_y(\theta) = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & e^{-i\theta/2} & 0 & 0 \\
	 *             0 & 0 & e^{i\theta/2} & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This is equivalent to applying TODO fix this -> \f$ e^{-i\theta (YX-XY)/2} \f$.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 * 
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param angle Angle to rotate by (in radians).
	 *
	 * Testing: same as crx
	 */	
	void nrz(unsigned targ1, unsigned targ2, F angle);

	/**
	 * \brief Perform a swap gate on two qubits. 
	 *
	 * \f[ 
	 * SWAP = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 0 & 1 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit to swap.
	 * \param targ2 The second qubit to swap.
	 *
	 * Testing: same as ch
	 */
	void swap(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform a fermionic-swap gate on two qubits. 
	 *
	 * \f[ 
	 * FSWAP = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 0 & 1 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 0 & -1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This corresponds to a swap gate followed by a CZ gate.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit to fswap.
	 * \param targ2 The second qubit to fswap.
	 *
	 * Testing: same as ch
	 */
	void fswap(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform an imaginary-swap gate on two qubits. 
	 *
	 * \f[ 
	 * ISWAP = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 0 & i & 0 \\
	 *             0 & i & 0 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This corresponds to the fixed-number \f[R_x(-\pi/2)\f].
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit to iswap.
	 * \param targ2 The second qubit to iswap.
	 *
	 * Testing: same as ch
	 */
	void iswap(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform a Hadamard gate on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace.  
	 *
	 * \f[ 
	 *  nH = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 1/\sqrt{2} & 1/\sqrt{2} & 0 \\
	 *             0 & 1/\sqrt{2} & -1/\sqrt{2} & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 *
	 * Testing: same as ch
	 */
	void nh(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform an arbitrary fixed-number gate, real matrix coefficients   
	 *
	 * \f[ 
	 *  nU = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & a_0 & a_1 & 0 \\
	 *             0 & a_2 & a_3 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * Use this function when all the coefficients of the unitary matrix 
	 * \f$a_i\f$ above are real.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 * A std::invalid_argument is thrown if a does not have orthonormal columns
	 * (i.e. cannot be normalised to a unitary matrix).
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param a The four coefficients of the matrix, stored in 
	 *          row-major order in a vector (in the order they are indexed above). 
	 *
	 * Testing: same as cu1
	 */
	void nu1(unsigned targ1, unsigned targ2, const std::vector<F> & a);

	/**
	 * \brief Perform an arbitrary fixed-number gate, complex matrix coefficients   
	 *
	 * \f[ 
	 *  nU = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & u_0 & u_1 & 0 \\
	 *             0 & u_2 & u_3 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * Use this function when all the coefficients of the unitary matrix 
	 * \f$u_i\f$ above are complex. If all the coefficients are real, use the
	 * real version of nu1 for improved performance.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 * A std::invalid_argument is thrown if u does not have orthonormal columns
	 * (i.e. cannot be normalised to a unitary matrix).
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param u The four complex coefficients of the matrix, stored in 
	 *          row-major order in a vector (in the order they are indexed above).
	 *
	 * Testing: same as cu1
	 */	
	void nu1(unsigned targ1, unsigned targ2, const std::vector<std::complex<F>> & u);

	/**
	 * \brief Perform an arbitrary two-qubit gate, real coefficients
	 *
	 * \f[ 
	 *  U = \begin{pmatrix}
	 *             a_0 & a_1 & a_2 & a_3 \\
	 *             a_4 & a_5 & a_6 & a_7 \\
	 *             a_8 & a_9 & a_{10} & a_{11} \\
	 *             a_{12} & a_{13} & a_{14} & a_{15}
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * Use this function when all the coefficients of the unitary matrix 
	 * \f$a_i\f$ above are real.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 * A std::invalid_argument is thrown if a does not have orthonormal columns
	 * (i.e. cannot be normalised to a unitary matrix).
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param a The 16 real coefficients of the matrix, stored in row-major 
	 *          order in a vector (in the order they are indexed above). 
	 *
	 * Testing: same as cu1
	 */		
	void u2(unsigned targ1, unsigned targ2, const std::vector<F> & a);

	/**
	 * \brief Perform an arbitrary two-qubit gate
	 *
	 * \f[ 
	 *  U = \begin{pmatrix}
	 *             u_0 & u_1 & u_2 & u_3 \\
	 *             u_4 & u_5 & u_6 & u_7 \\
	 *             u_8 & u_9 & u_{10} & u_{11} \\
	 *             u_{12} & u_{13} & u_{14} & u_{15}
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * Use this function when all the coefficients of the unitary matrix 
	 * \f$u_i\f$ above are complex. If all the coefficients are real, use the
	 * real version of u2 for improved performance.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 * A std::invalid_argument is thrown if u does not have orthonormal columns
	 * (i.e. cannot be normalised to a unitary matrix).
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param u The 16 complex coefficients of the matrix, stored in row-major 
	 *          order in a vector (in the order they are indexed above).
	 *
	 * Testing: same as cu1
	 */		
	void u2(unsigned targ1, unsigned targ2, const std::vector<std::complex<F>> & u);

	/** 
	 * \brief Calculate probability of specific outcome of a qubit
	 * 
	 * Calculate the probability that a specific outcome will result from the
	 * measurement of qubit, without measuring the qubit. The state vector
	 * is not affected by this function.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1. A std::invalid_argument is thrown if outcome is not 0 or 1.
	 *
	 * \param targ The qubit under consideration
	 * \param outcome The outcome whose probability is desired.
	 * \return The probability of the specifed outcome.
	 *
	 * Testing:
	 * - In debug mode, check exceptions thrown e.g. input out of range targ, outcome
	 *   not set to 0 or 1.
	 * - Find prob of random state (smallish no. of qubits) and check output 
	 *   matches with armadillo. Do this with targ = 0 up to num_qubits-1, for 0 and 
	 *   1 outcomes to make sure edge cases have been checked. Can probably be checked
	 *   along with the postselect function.
	 */
	F prob(unsigned targ, unsigned outcome) const;

	/** 
	 * \brief Measure one qubit
	 *
	 * Measure a particular qubit and return the outcome. The function chooses
	 * the outcome at random, depending on the probability defined by the amplitudes
	 * in the state vector.
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1.
	 * 
	 * \param targ The qubit to measure
	 * \param g A generator used as the source of randomness
	 * \return The measured outcome (either 0 or 1)
	 *
	 * Testing:
	 * - In debug mode, check exceptions thrown e.g. input out of range targ.
	 * - Measure random state (smallish no. of qubits) and check output 
	 *   matches with armadillo. Do this with targ = 0 up to num_qubits-1 to make sure 
	 *   edge cases have been checked.
	 * - Check randomness seeding by preparing the same random state and checking we
	 *   get the same outcome everytime.
	 * - Check that this function samples from the actual qubit probablity distribution. 
	 *   Can be done by generating a bunch of outcomes, checking distribution in python
	 *   and then making sure the seeding reproduces these outcomes.
	 */
	unsigned measure(unsigned targ, std::uniform_random_bit_generator auto g = gen);

	/** 
	 * \brief Measure all the qubits 
	 *
	 * Measure all the qubits and return the resulting computational basis state as
	 * an integer, whose bits represent measured outcomes. The least significant bits
	 * in the returned integer represent the lowest-index measurement outcomes. The 
	 * outcome for the qubit measurements are chosen randomly, with probabilities 
	 * depending on the amplitudes in the state vector.
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 * 
	 * \param g A generator used as the source of randomness
	 * \return The measurement outcomes (packed little-endian into an integer)
	 *
	 * Testing:
	 * - Measure random state (smallish no. of qubits) and check output 
	 *   matches with armadillo.
	 * - Check randomness seeding by preparing the same random state and checking we
	 *   get the same outcome everytime.
	 * - Check that this function samples from the actual qubit probablity distribution. 
	 *   Can be done by generating a bunch of outcomes, checking distribution in python
	 *   and then making sure the seeding reproduces these outcomes.
	 */
	std::size_t measure_all(std::uniform_random_bit_generator auto g = gen);

	/** 
	 * \brief Collapse one qubit to a specific outcome
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1. A std::invalid_argument is thrown if outcome is not 0 or 1.
	 * 
	 * \param targ The qubit to measure
	 * \param outcome The desired outcome of targ
	 * \return The measured outcome (either 0 or 1), this is the same as outcome
	 *         but is returned for consistency with the measure function.
	 *
	 * Testing:
	 * - In debug mode, check exceptions thrown e.g. input out of range targ, outcome
	 *   not set to 0 or 1.
	 * - Postselect a random state (smallish no. of qubits) and check output 
	 *   matches with armadillo. Do this with targ = 0 up to num_qubits-1, for 0 and 
	 *   1 outcomes to make sure edge cases have been checked. Can probably be checked
	 *   along with the prob function.
	 */
	unsigned postselect(unsigned targ, unsigned outcome);
	
	/** 
	 * \brief Sample the measurement outcome of one qubit repeatedly
	 *
	 * If only a measurement outcome is required from a measurement operation (e.g.
	 * the state will not be used afterwards), it is more efficient to not collapse 
	 * the state vector. This function samples from the probability distribution of the
	 * qubit many times and is equivalent to preparing the same circuit and measuring
	 * qubit targ samples times. The state vector is not modified.
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1.
	 *
	 * \todo Think about whether measure and sample should agree on outcomes (using same seed).
	 * 
	 * \param targ The qubit to sample
	 * \param samples Number of samples to take
	 * \param g A generator used as the source of randomness
	 * \return A length-two vector where element 0 is the number of times 0 was
	 *         sampled, and element 1 the number of times 1 was sampled.
	 *
	 * Testing:
	 * - In debug mode, check exceptions thrown e.g. input out of range targ.
	 * - Measure random state (smallish no. of qubits) and check output 
	 *   matches with armadillo. Do this with targ = 0 up to num_qubits-1 to make sure 
	 *   edge cases have been checked.
	 * - Check randomness seeding by preparing the same random state and checking we
	 *   get the same samples everytime.
	 * - Check that this function samples from the actual qubit probablity distribution. 
	 *   Can be done by generating a bunch of outcomes, checking distribution in python
	 *   and then making sure the seeding reproduces these outcomes.
	 */
	std::vector<std::size_t> sample(unsigned targ, std::size_t samples,
					std::uniform_random_bit_generator auto g = gen) const;

	/** 
	 * \brief Sample the measurement outcome of all the qubits repeatedly
	 *
	 * This function is useful if at the end of a circuit all the qubits need to be
	 * measured (resulting in a computational basis state as an outcome), 
	 * and the process repeated many times. This function generates measurement samples 
	 * from the probability distribution of all the qubits, and does not modify the
	 * state vector.  
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 * 
	 * \param samples Number of samples to take
	 * \param g A generator used as the source of randomness
	 * \return A map of the measurement outcomes where the keys are the computational
	 *         basis state outcomes (encoded as a bitstring) and the associated values
	 *         are the number of times that outcome was measured.
	 *
	 * Testing:
	 * - Measure random state (smallish no. of qubits) and check output 
	 *   matches with armadillo.
	 * - Check randomness seeding by preparing the same random state and checking we
	 *   get the same outcome everytime.
	 * - Check that this function samples from the actual qubit probablity distribution. 
	 *   Can be done by generating a bunch of outcomes, checking distribution in python
	 *   and then making sure the seeding reproduces these outcomes.
	 */
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
	using value_type = std::complex<F>;
	constexpr bool debug() { return D; }

	/**
	 * \brief Construct an empty simulator
	 *
	 * Constructs a simulator with zero qubits. The number of qubits
	 * can be changed using the qsl::resize::reset member function, or by 
	 * calling the qsl::resize::add_qubit method. This function
	 * is not especially useful for one or two simulators, but may help
	 * in contexts that require construction-before-copy (for example,
	 * adding a simulator into a map without using emplace).
	 *
	 */
	resize();
	
	/**
	 * \brief Initialise the class with a specified number of qubits.
	 *
	 * The simulator is initialised in the all-zero state.
	 *
	 * In debug mode, if the number of qubits is too large to simulate, a
	 * std::runtime_error is thrown.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 */ 
	explicit resize(unsigned num_qubits);

	/**
	 * \brief Instantiate a simulator based on a std::vector or another simulator
         *        that has the same floating point precision.
	 *
	 * If inputting a std::vector, it must be non-zero and have a length which is a 
	 * power of two. It does not need to be normalised as this function will carry
	 * out a normalisation step.
	 *
	 * In debug mode, if the input state is not a power of two or is all zeros,
	 * a std::invalid_argument is thrown.
	 *
	 * This is not a copy constructor (because it is templated), so it will not
	 * cause the move constructors to be implicitly deleted.
	 *
	 * \param s A valid qsl simulator object with the same floating point precision,
         *          or std::vector<std::complex<F>> or std::vector<F>.
	 */
	template<state_vector S>
	requires same_precision<F, S>
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
	std::vector<std::complex<F>> state() const;

	/**
	 * \brief Change the simulator state to the state that is passed in. 
	 *
	 * If a stad::vector, the input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. Note that this function
	 * allows for a change in the number of qubits.
	 *
	 * In debug mode, if the input state is not a power of two or is all zeros,
	 * a std::invalid_argument is thrown.
	 *
	 * \param state A vector or qsl simulator containing the new state for the object.
	 */ 		
	template<state_vector S>
	requires same_precision<F, S>
	resize & operator= (const S & state);

	/**
	 * \brief Reset to the all-zero computational basis state.
	 */
	void reset();

	/**
	 * \brief Change number of qubits and reset to the all-zero computational basis state.
	 *
	 * In debug mode, if the number of qubits is too large to simulate, a
	 * std::runtime_error is thrown.
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
	 * A std::out_of_range error is thrown if index is bigger than size()-1.
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
	 * A std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The value the inserted qubit will be indexed as.
	 */
	void add_qubit(unsigned targ);

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
	 * A std::runtime_error is thrown if too many qubits are to be stored
	 * or if the file cannot be created. These exceptions are still thrown if
	 * the simulator is not in debug mode.
	 *
	 * \param file Filename to save to.
	 */
	void save_json(const std::filesystem::path & file) const;

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
	 * A std::runtime_error is thrown if the file doesn't exist or cannot be
	 * read. A std::invalid_argument is thrown if it does not contain a valid json object,
	 * if it reads a json object that does not contain a 'state' field, or if the state that
	 * is read is invalid. These exceptions occur even when debugging is disabled.
	 *
	 * \param file The filename to read from. 
	 */
	void load_json(const std::filesystem::path & file);

	/**
	 * \brief Rotate around the x-axis of the Bloch sphere 
	 * \f$ e^{-i\theta X/2} \f$:
	 * 
	 * \f[ 
	 * R_x = \begin{pmatrix}
	 *       \cos(\theta/2) & -i\sin(\theta/2) \\
	 *       -i\sin(\theta/2) & \cos(\theta/2) \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */
	void rx(unsigned targ, F angle);

	/**
	 * \brief Rotate around the y-axis of the Bloch sphere 
	 * \f$ e^{-i\theta Y/2} \f$:
	 * 
	 * \f[ 
	 * R_y = \begin{pmatrix}
	 *       \cos(\theta/2) & -\sin(\theta/2) \\
	 *       \sin(\theta/2) & \cos(\theta/2) \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */
	void ry(unsigned targ, F angle);

	/**
	 * \brief Rotate around the z-axis of the Bloch sphere 
	 * \f$ e^{-i\theta Z/2} \f$:
	 * 
	 * \f[ 
	 * R_z = \begin{pmatrix}
	 *       e^{-i\theta/2} & 0 \\
	 *       0 & e^{i\theta/2} \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */
	void rz(unsigned targ, F angle);  

	/**
	 * \brief Apply a phase shift to qubit targ:
	 *
	 * \f[
	 * \text{phase}(\theta) = \begin{pmatrix}
	 *       1 & 0 \\
	 *       0 & e^{i\theta} \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * Also equivalent to \f$ e^{i\theta/2} R_z(\theta) \f$.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */	
	void phase(unsigned targ, F angle); 

	/**
	 * \brief Apply the Hadamard gate to qubit targ:
	 *
	 * \f[ 
	 * H = \frac{1}{\sqrt{2}}\begin{pmatrix}
	 *     1 & 1 \\
	 *     1 & -1 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 */
	void h(unsigned targ);

	/**
	 * \brief Apply the Pauli-X gate to qubit targ:
	 *
	 * \f[ 
	 * X = \begin{pmatrix}
	 *     0 & 1 \\
	 *     1 & 0 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * Also equivalent to \f$ iR_x(\pi) \f$.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 */
	void x(unsigned targ);

	/**
	 * \brief Apply the Pauli-Y gate to qubit targ:
	 *
	 * \f[ 
	 * Y = \begin{pmatrix}
	 *     0 & -i \\
	 *     i & 0 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * Also equivalent to \f$ iR_y(\pi) \f$.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 */
	void y(unsigned targ);

	/**
	 * \brief Apply the Pauli-Z gate to qubit targ:
	 *
	 * \f[ 
	 * Z = \begin{pmatrix}
	 *     1 & 0 \\
	 *     0 & -1 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * Also equivalent to \f$ iR_z(\pi) \f$ or \f$ \text{phase}(\pi) \f$
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 */
	void z(unsigned targ);

	/**
	 * \brief Apply an arbitrary one-qubit (real) unitary to targ.
	 *
	 * The matrix must have orthonormal columns, the columns will be
	 * normalised in this function.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param matrix The unitary matrix to apply, in row-major form.
	 */
	void u1(unsigned targ, const std::vector<F> & matrix);

	/**
	 * \brief Apply an arbitrary one-qubit unitary to targ.
	 *
	 * The matrix must have orthonormal columns, the columns will be
	 * normalised in this function.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param matrix The unitary matrix to apply in row-major form.
	 */
	void u1(unsigned targ, const std::vector<std::complex<F>> & matrix);

	/**
	 * \brief Perform a controlled X-rotation on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, \f$ R_x \f$ is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */
	void crx(unsigned ctrl, unsigned targ, F angle);

	/**
	 * \brief Perform a controlled Y-rotation on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, \f$ R_y \f$ is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */
	void cry(unsigned ctrl, unsigned targ, F angle);

	/**
	 * \brief Perform a controlled Z-rotation on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, \f$ R_z \f$ is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */
	void crz(unsigned ctrl, unsigned targ, F angle);  

	/**
	 * \brief Perform a controlled phase gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, a phase shift is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */
	void cphase(unsigned ctrl, unsigned targ, F angle); 

	/**
	 * \brief Perform a controlled Hadamard gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, H is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 */
	void ch(unsigned ctrl, unsigned targ);

	/**
	 * \brief Perform a controlled Pauli-X (often referred to as a controlled-Not) 
	 *        gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, X is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 */
	void cnot(unsigned ctrl, unsigned targ);

	/**
	 * \brief Perform a controlled Pauli-Y gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, Y is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 */
	void cy(unsigned ctrl, unsigned targ);

	/**
	 * \brief Perform a controlled Pauli-Z gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, Z is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 */
	void cz(unsigned ctrl, unsigned targ);

	/**
	 * \brief Apply a controlled arbitrary one-qubit (real) unitary.
	 *
	 * The matrix must have orthonormal columns, the columns will be
	 * normalised in this function.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, U is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param matrix The unitary matrix to apply, in row-major form.
	 */
	void cu1(unsigned ctrl, unsigned targ, const std::vector<F> & matrix);

	/**
	 * \brief Apply a controlled arbitrary one-qubit unitary.
	 *
	 * The matrix must have orthonormal columns, the columns will be
	 * normalised in this function.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, U is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param matrix The unitary matrix to apply, in row-major form.
	 */
	void cu1(unsigned ctrl, unsigned targ, const std::vector<std::complex<F>> & matrix);

	/**
	 * \brief Perform an X rotation on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace. 
	 *
	 * \f[ 
	 *  nR_x(\theta) = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & \cos(\theta/2) & -i\sin(\theta/2) & 0 \\
	 *             0 & -i\sin(\theta/2) & \cos(\theta/2) & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This is equivalent to applying \f$ e^{-i\theta (XX+YY)/2} \f$. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param angle The angle to rotate by (in radians)
	 */
	void nrx(unsigned targ1, unsigned targ2, F angle);

	/**
	 * \brief Perform a Y rotation on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace. This is equivalent to applying 
	 * \f$ e^{-i\theta (YX-XY)/2} \f$. 
	 *
	 * \ingroup qubits_gates

	 * \f[ 
	 *  nR_y(\theta) = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & \cos(\theta/2) & -\sin(\theta/2) & 0 \\
	 *             0 & \sin(\theta/2) & \cos(\theta/2) & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param angle Angle to rotate by (in radians).
	 */
	void nry(unsigned targ1, unsigned targ2, F angle);

	/**
	 * \brief Perform a Z rotation on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace. 
	 *
	 * \f[ 
	 *  nR_y(\theta) = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & e^{-i\theta/2} & 0 & 0 \\
	 *             0 & 0 & e^{i\theta/2} & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This is equivalent to applying TODO fix this -> \f$ e^{-i\theta (YX-XY)/2} \f$.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 * 
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param angle Angle to rotate by (in radians).
	 */	
	void nrz(unsigned targ1, unsigned targ2, F angle);

	/**
	 * \brief Perform a swap gate on two qubits. 
	 *
	 * \f[ 
	 * SWAP = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 0 & 1 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit to swap.
	 * \param targ2 The second qubit to swap.
	 */
	void swap(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform a fermionic-swap gate on two qubits. 
	 *
	 * \f[ 
	 * FSWAP = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 0 & 1 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 0 & -1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This corresponds to a swap gate followed by a CZ gate.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit to fswap.
	 * \param targ2 The second qubit to fswap.
	 */
	void fswap(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform an imaginary-swap gate on two qubits. 
	 *
	 * \f[ 
	 * ISWAP = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 0 & i & 0 \\
	 *             0 & i & 0 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This corresponds to the fixed-number \f[R_x(-\pi/2)\f].
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit to iswap.
	 * \param targ2 The second qubit to iswap.
	 */
	void iswap(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform a Hadamard gate on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace.  
	 *
	 * \f[ 
	 *  nH = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 1/\sqrt{2} & 1/\sqrt{2} & 0 \\
	 *             0 & 1/\sqrt{2} & -1/\sqrt{2} & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 */
	void nh(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform an arbitrary fixed-number gate, real matrix coefficients   
	 *
	 * \f[ 
	 *  nU = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & a_0 & a_1 & 0 \\
	 *             0 & a_2 & a_3 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * Use this function when all the coefficients of the unitary matrix 
	 * \f$a_i\f$ above are real.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param a The four coefficients of the matrix, stored in 
	 *          row-major order in a vector (in the order they are indexed above). 
	 */
	void nu1(unsigned targ1, unsigned targ2, const std::vector<F> & a);

	/**
	 * \brief Perform an arbitrary fixed-number gate, complex matrix coefficients   
	 *
	 * \f[ 
	 *  nU = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & u_0 & u_1 & 0 \\
	 *             0 & u_2 & u_3 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * Use this function when all the coefficients of the unitary matrix 
	 * \f$u_i\f$ above are complex. If all the coefficients are real, use the
	 * real version of nu1 for improved performance.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param u The four complex coefficients of the matrix, stored in 
	 *          row-major order in a vector (in the order they are indexed above).
	 */	
	void nu1(unsigned targ1, unsigned targ2, const std::vector<std::complex<F>> & u);

	/**
	 * \brief Perform an arbitrary two-qubit gate, real coefficients
	 *
	 * \f[ 
	 *  U = \begin{pmatrix}
	 *             a_0 & a_1 & a_2 & a_3 \\
	 *             a_4 & a_5 & a_6 & a_7 \\
	 *             a_8 & a_9 & a_{10} & a_{11} \\
	 *             a_{12} & a_{13} & a_{14} & a_{15}
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * Use this function when all the coefficients of the unitary matrix 
	 * \f$a_i\f$ above are real.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param a The 16 real coefficients of the matrix, stored in row-major 
	 *          order in a vector (in the order they are indexed above). 
	 */		
	void u2(unsigned targ1, unsigned targ2, const std::vector<F> & a);

	/**
	 * \brief Perform an arbitrary two-qubit gate
	 *
	 * \f[ 
	 *  U = \begin{pmatrix}
	 *             u_0 & u_1 & u_2 & u_3 \\
	 *             u_4 & u_5 & u_6 & u_7 \\
	 *             u_8 & u_9 & u_{10} & u_{11} \\
	 *             u_{12} & u_{13} & u_{14} & u_{15}
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * Use this function when all the coefficients of the unitary matrix 
	 * \f$u_i\f$ above are complex. If all the coefficients are real, use the
	 * real version of u2 for improved performance.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param u The 16 complex coefficients of the matrix, stored in row-major 
	 *          order in a vector (in the order they are indexed above).
	 */		
	void u2(unsigned targ1, unsigned targ2, const std::vector<std::complex<F>> & u);

	/** 
	 * \brief Calculate probability of specific outcome of a qubit
	 * 
	 * Calculate the probability that a specific outcome will result from the
	 * measurement of qubit, without measuring the qubit. The state vector
	 * is not affected by this function.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1. A std::invalid_argument is thrown if outcome is not 0 or 1.
	 *
	 * \param targ The qubit under consideration
	 * \param outcome The outcome whose probability is desired.
	 * \return The probability of the specifed outcome.
	 */
	F prob(unsigned targ, unsigned outcome) const;

	/** 
	 * \brief Measure one qubit
	 *
	 * Measure a particular qubit and return the outcome. The function chooses
	 * the outcome at random, depending on the probability defined by the amplitudes
	 * in the state vector.
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1.
	 * 
	 * \param targ The qubit to measure
	 * \param g A generator used as the source of randomness
	 * \return The measured outcome (either 0 or 1)
	 */
	unsigned measure(unsigned targ, std::uniform_random_bit_generator auto g = gen);

	/** 
	 * \brief Measure all the qubits 
	 *
	 * Measure all the qubits and return the resulting computational basis state as
	 * an integer, whose bits represent measured outcomes. The least significant bits
	 * in the returned integer represent the lowest-index measurement outcomes. The 
	 * outcome for the qubit measurements are chosen randomly, with probabilities 
	 * depending on the amplitudes in the state vector.
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 * 
	 * \param g A generator used as the source of randomness
	 * \return The measurement outcomes (packed little-endian into an integer)
	 */
	std::size_t measure_all(std::uniform_random_bit_generator auto g = gen);

	/** 
	 * \brief Collapse one qubit to a specific outcome
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1. A std::invalid_argument is thrown if outcome is not 0 or 1.
	 * 
	 * \param targ The qubit to measure
	 * \param outcome The desired outcome of targ
	 * \return The measured outcome (either 0 or 1), this is the same as outcome
	 *         but is returned for consistency with the measure function.
	 */
	unsigned postselect(unsigned targ, unsigned outcome);
	
	/** 
	 * \brief Sample the measurement outcome of one qubit repeatedly
	 *
	 * If only a measurement outcome is required from a measurement operation (e.g.
	 * the state will not be used afterwards), it is more efficient to not collapse 
	 * the state vector. This function samples from the probability distribution of the
	 * qubit many times and is equivalent to preparing the same circuit and measuring
	 * qubit targ samples times. The state vector is not modified.
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1.
	 * 
	 * \param targ The qubit to sample
	 * \param samples Number of samples to take
	 * \param g A generator used as the source of randomness
	 * \return A length-two vector where element 0 is the number of times 0 was
	 *         sampled, and element 1 the number of times 1 was sampled.
	 */
	std::vector<std::size_t> sample(unsigned targ, std::size_t samples,
					std::uniform_random_bit_generator auto g = gen) const;

	/** 
	 * \brief Sample the measurement outcome of all the qubits repeatedly
	 *
	 * This function is useful if at the end of a circuit all the qubits need to be
	 * measured (resulting in a computational basis state as an outcome), 
	 * and the process repeated many times. This function generates measurement samples 
	 * from the probability distribution of all the qubits, and does not modify the
	 * state vector.  
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 * 
	 * \param samples Number of samples to take
	 * \param g A generator used as the source of randomness
	 * \return A map of the measurement outcomes where the keys are the computational
	 *         basis state outcomes (encoded as a bitstring) and the associated values
	 *         are the number of times that outcome was measured.
	 */
	std::map<std::size_t, std::size_t> sample_all(std::size_t samples,
						      std::uniform_random_bit_generator auto g = gen) const;

	/** 
	 * \brief Measure one qubit and remove from the state vector.
	 *
	 * Measure a particular qubit and return the outcome. Remove the measured 
	 * qubit from the state vector, reducing the number of qubits by one. The function chooses
	 * the outcome at random, depending on the probability defined by the amplitudes
	 * in the state vector.
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1.
	 * 
	 * \param targ The qubit to measure
	 * \param g A generator used as the source of randomness
	 * \return The measured outcome (either 0 or 1)
	 */
	unsigned measure_out(unsigned targ, std::uniform_random_bit_generator auto g = gen);

	/** 
	 * \brief Collapse one qubit to a specific outcome, and remove from the state vector.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1. A std::invalid_argument is thrown if outcome is not 0 or 1.
	 * 
	 * \param targ The qubit to measure
	 * \param outcome The desired outcome of targ
	 * \return The measured outcome (either 0 or 1), this is the same as outcome
	 *         but is returned for consistency with the measure function.
	 */
	unsigned postselect_out(unsigned targ, unsigned outcome);
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
	using value_type = std::complex<F>;
	constexpr bool debug() { return D; }

	/**
	 * \brief Construct an empty simulator
	 *
	 * Constructs a simulator with zero qubits. The number of qubits
	 * can be changed using the qsl::number::reset member function. This 
	 * function is not especially useful for one or two simulators, but 
	 * may help in contexts that require construction-before-copy (for 
	 * example, adding a simulator into a map without using emplace).
	 *
	 */
	number();
	
	/**
	 * \brief Initialise the class with a specified number of qubits.
	 *
	 * The simulator is initialised in the all-zero state, and the number of ones
	 * is set to zero.
	 *
	 * In debug mode, if the number of qubits is too large to simulate, a
	 * std::runtime_error is thrown.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 */ 
	explicit number(unsigned num_qubits);

	/**
	 * \brief Instantiate a simulator based on a std::vector or another simulator
         *        that has the same floating point precision.
	 *
	 * If inputting a std::vector, it must be non-zero and have a length which is a 
	 * power of two. It does not need to be normalised as this function will carry
	 * out a normalisation step.
	 *
	 * In debug mode, if the input state is not a power of two or is all zeros,
	 * a std::invalid_argument is thrown.
	 *
	 * This is not a copy constructor (because it is templated), so it will not
	 * cause the move constructors to be implicitly deleted.
	 *
	 * \param s A valid qsl simulator object with the same floating point precision,
         *          or std::vector<std::complex<F>> or std::vector<F>.
	 */
	template<state_vector S>
	requires same_precision<F, S>
	explicit number(const S & s);
	
	/**
	 * \brief Initialise the class with a specified number of qubits and number of ones.
	 *
	 * The simulator is initialised in the lowest indexed computational basis state
	 * with num_ones ones.
	 *
	 * In debug mode, if the number of qubits is too large to simulate, a
	 * std::runtime_error is thrown. If num_ones is bigger than num_qubits, 
	 * std::invalid_argument is thrown.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 * \param num_ones The number of ones in the fixed number simulator.
	 */ 
	number(unsigned num_qubits, unsigned num_ones);

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
	std::vector<std::complex<F>> state() const;

	/**
	 * \brief Change the simulator state to the state that is passed in. 
	 *
	 * If a stad::vector, the input state vector must be non-zero and have a length which
	 * is a power of two. It does not need to be normalised as
	 * this function will carry out a normalisation step. Note that this function
	 * allows for a change in the number of qubits.
	 *
	 * In debug mode, if the input state is not a power of two or is all zeros,
	 * a std::invalid_argument is thrown.
	 *
	 * \param state A vector or qsl simulator containing the new state for the object.
	 */ 		
	template<state_vector S>
	requires same_precision<F, S>
	number & operator= (const S & state);

	/**
	 * \brief Reset to the lowest indexed computational basis state for the given num_ones.
	 */
	void reset();

	/**
	 * \brief Change number of qubits and reset to the all-zero computational basis state.
	 *
	 * The number of ones is also reset to zero.
	 *
	 * In debug mode, if the number of qubits is too large to simulate, a
	 * std::runtime_error is thrown.
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
	 * In debug mode, if the number of qubits is too large to simulate, a
	 * std::runtime_error is thrown. If num_ones is bigger than num_qubits, 
	 * std::invalid_argument is thrown.
	 *
	 * \param num_qubits The number of qubits to simulate.
	 * \param num_ones The number of ones to set.
	 */
	void reset(unsigned num_qubits, unsigned num_ones);
	    
	/**
	 * \brief Access state vector elements (read-only).
	 *
	 * This is read-only to avoid accidental tampering with the state vector
	 * that might render the state invalid (for example, setting every element to zero).  
	 *
	 * A std::out_of_range error is thrown if index is bigger than size()-1.
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
	 * A std::runtime_error is thrown if too many qubits are to be stored
	 * or if the file cannot be created. These exceptions are still thrown if
	 * the simulator is not in debug mode.
	 *
	 * \param file Filename to save to.
	 */	
	void save_json(const std::filesystem::path & file) const;

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
	 * A std::runtime_error is thrown if the file doesn't exist or cannot be
	 * read. A std::invalid_argument is thrown if it does not contain a valid json object,
	 * if it reads a json object that does not contain a 'state' field, or if the state that
	 * is read is invalid. These exceptions occur even when debugging is disabled.
	 *
	 * \param file The filename to read from. 
	 */
	void load_json(const std::filesystem::path & file);
	
	/**
	 * \brief Rotate around the z-axis of the Bloch sphere 
	 * \f$ e^{-i\theta Z/2} \f$:
	 * 
	 * \f[ 
	 * R_z = \begin{pmatrix}
	 *       e^{-i\theta/2} & 0 \\
	 *       0 & e^{i\theta/2} \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */
	void rz(unsigned targ, F angle);

	/**
	 * \brief Apply a phase shift to qubit targ:
	 *
	 * \f[
	 * \text{phase}(\theta) = \begin{pmatrix}
	 *       1 & 0 \\
	 *       0 & e^{i\theta} \\
	 *       \end{pmatrix} 
	 * \f]
	 *
	 * Also equivalent to \f$ e^{i\theta/2} R_z(\theta) \f$.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */	
	void phase(unsigned targ, F angle);

	/**
	 * \brief Apply the Pauli-Z gate to qubit targ:
	 *
	 * \f[ 
	 * Z = \begin{pmatrix}
	 *     1 & 0 \\
	 *     0 & -1 \\
	 *     \end{pmatrix} 
	 * \f]
	 *
	 * Also equivalent to \f$ iR_z(\pi) \f$ or \f$ \text{phase}(\pi) \f$
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than num_qubits-1.
	 *
	 * \param targ The target qubit.
	 */
	void z(unsigned targ);  

	/**
	 * \brief Perform a controlled Z-rotation on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, \f$ R_z \f$ is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */	
	void crz(unsigned ctrl, unsigned targ, F angle);  

	/**
	 * \brief Perform a controlled phase gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, a phase shift is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 * \param angle The angle to rotate the qubit by (in radians).
	 */
	void cphase(unsigned ctrl, unsigned targ, F angle); 

	/**
	 * \brief Perform a controlled Pauli-Z gate on two qubits. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ or ctrl is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if ctrl = targ.
	 *
	 * \param ctrl The control qubit, Z is applied on the target qubit
	 *             if this qubit is \f$ |1\rangle \f$.
	 * \param targ The target qubit.
	 */
	void cz(unsigned ctrl, unsigned targ);  

	/**
	 * \brief Perform an X rotation on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace. 
	 *
	 * \f[ 
	 *  nR_x(\theta) = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & \cos(\theta/2) & -i\sin(\theta/2) & 0 \\
	 *             0 & -i\sin(\theta/2) & \cos(\theta/2) & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This is equivalent to applying \f$ e^{-i\theta (XX+YY)/2} \f$. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param angle The angle to rotate by (in radians)
	 */
	void nrx(unsigned targ1, unsigned targ2, F angle);

	/**
	 * \brief Perform a Y rotation on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace. This is equivalent to applying 
	 * \f$ e^{-i\theta (YX-XY)/2} \f$. 
	 *
	 * \ingroup qubits_gates

	 * \f[ 
	 *  nR_y(\theta) = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & \cos(\theta/2) & -\sin(\theta/2) & 0 \\
	 *             0 & \sin(\theta/2) & \cos(\theta/2) & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param angle Angle to rotate by (in radians).
	 */
	void nry(unsigned targ1, unsigned targ2, F angle);

	/**
	 * \brief Perform a Z rotation on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace. 
	 *
	 * \f[ 
	 *  nR_y(\theta) = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & e^{-i\theta/2} & 0 & 0 \\
	 *             0 & 0 & e^{i\theta/2} & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This is equivalent to applying TODO fix this -> \f$ e^{-i\theta (YX-XY)/2} \f$.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 * 
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param angle Angle to rotate by (in radians).
	 */	
	void nrz(unsigned targ1, unsigned targ2, F angle);

	/**
	 * \brief Perform a swap gate on two qubits. 
	 *
	 * \f[ 
	 * SWAP = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 0 & 1 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit to swap.
	 * \param targ2 The second qubit to swap.
	 */
	void swap(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform a fermionic-swap gate on two qubits. 
	 *
	 * \f[ 
	 * FSWAP = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 0 & 1 & 0 \\
	 *             0 & 1 & 0 & 0 \\
	 *             0 & 0 & 0 & -1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This corresponds to a swap gate followed by a CZ gate.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit to fswap.
	 * \param targ2 The second qubit to fswap.
	 */
	void fswap(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform an imaginary-swap gate on two qubits. 
	 *
	 * \f[ 
	 * ISWAP = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 0 & i & 0 \\
	 *             0 & i & 0 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * This corresponds to the fixed-number \f[R_x(-\pi/2)\f].
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit to iswap.
	 * \param targ2 The second qubit to iswap.
	 */
	void iswap(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform a Hadamard gate on the \f$ \{|01\rangle,|10\rangle\}\f$
	 * subspace.  
	 *
	 * \f[ 
	 *  nH = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & 1/\sqrt{2} & 1/\sqrt{2} & 0 \\
	 *             0 & 1/\sqrt{2} & -1/\sqrt{2} & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 */
	void nh(unsigned targ1, unsigned targ2);

	/**
	 * \brief Perform an arbitrary fixed-number gate, real matrix coefficients   
	 *
	 * \f[ 
	 *  nU = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & a_0 & a_1 & 0 \\
	 *             0 & a_2 & a_3 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * Use this function when all the coefficients of the unitary matrix 
	 * \f$a_i\f$ above are real.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param a The four coefficients of the matrix, stored in 
	 *          row-major order in a vector (in the order they are indexed above). 
	 */
	void nu1(unsigned targ1, unsigned targ2, const std::vector<F> & a);

	/**
	 * \brief Perform an arbitrary fixed-number gate, complex matrix coefficients   
	 *
	 * \f[ 
	 *  nU = \begin{pmatrix}
	 *             1 & 0 & 0 & 0 \\
	 *             0 & u_0 & u_1 & 0 \\
	 *             0 & u_2 & u_3 & 0 \\
	 *             0 & 0 & 0 & 1
	 *             \end{pmatrix} 
	 * \f]
	 *
	 * Use this function when all the coefficients of the unitary matrix 
	 * \f$u_i\f$ above are complex. If all the coefficients are real, use the
	 * real version of nu1 for improved performance.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ1 or targ2 is bigger 
	 * than num_qubits-1. A std::invalid_argument is thrown if targ1 = targ2.
	 *
	 * \param targ1 The first qubit.
	 * \param targ2 The second qubit.
	 * \param u The four complex coefficients of the matrix, stored in 
	 *          row-major order in a vector (in the order they are indexed above).
	 */	
	void nu1(unsigned targ1, unsigned targ2, const std::vector<std::complex<F>> & u);

	/** 
	 * \brief Calculate probability of specific outcome of a qubit
	 * 
	 * Calculate the probability that a specific outcome will result from the
	 * measurement of qubit, without measuring the qubit. The state vector
	 * is not affected by this function.
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1. A std::invalid_argument is thrown if outcome is not 0 or 1.
	 *
	 * \param targ The qubit under consideration
	 * \param outcome The outcome whose probability is desired.
	 * \return The probability of the specifed outcome.
	 */
	F prob(unsigned targ, unsigned outcome) const;

	/** 
	 * \brief Measure one qubit
	 *
	 * Measure a particular qubit and return the outcome. The function chooses
	 * the outcome at random, depending on the probability defined by the amplitudes
	 * in the state vector.
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1.
	 * 
	 * \param targ The qubit to measure
	 * \param g A generator used as the source of randomness
	 * \return The measured outcome (either 0 or 1)
	 */
	unsigned measure(unsigned targ, std::uniform_random_bit_generator auto g = gen);

	/** 
	 * \brief Measure all the qubits 
	 *
	 * Measure all the qubits and return the resulting computational basis state as
	 * an integer, whose bits represent measured outcomes. The least significant bits
	 * in the returned integer represent the lowest-index measurement outcomes. The 
	 * outcome for the qubit measurements are chosen randomly, with probabilities 
	 * depending on the amplitudes in the state vector.
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 * 
	 * \param g A generator used as the source of randomness
	 * \return The measurement outcomes (packed little-endian into an integer)
	 */
	std::size_t measure_all(std::uniform_random_bit_generator auto g = gen);

	/** 
	 * \brief Collapse one qubit to a specific outcome
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1. A std::invalid_argument is thrown if outcome is not 0 or 1.
	 * 
	 * \param targ The qubit to measure
	 * \param outcome The desired outcome of targ
	 * \return The measured outcome (either 0 or 1), this is the same as outcome
	 *         but is returned for consistency with the measure function.
	 */
	unsigned postselect(unsigned targ, unsigned outcome);
	
	/** 
	 * \brief Sample the measurement outcome of one qubit repeatedly
	 *
	 * If only a measurement outcome is required from a measurement operation (e.g.
	 * the state will not be used afterwards), it is more efficient to not collapse 
	 * the state vector. This function samples from the probability distribution of the
	 * qubit many times and is equivalent to preparing the same circuit and measuring
	 * qubit targ samples times. The state vector is not modified.
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 *
	 * In debug mode, a std::out_of_range error is thrown if targ is bigger than
	 * num_qubits-1.
	 * 
	 * \param targ The qubit to sample
	 * \param samples Number of samples to take
	 * \param g A generator used as the source of randomness
	 * \return A length-two vector where element 0 is the number of times 0 was
	 *         sampled, and element 1 the number of times 1 was sampled.
	 */
	std::vector<std::size_t> sample(unsigned targ, std::size_t samples,
					std::uniform_random_bit_generator auto g = gen) const;

	/** 
	 * \brief Sample the measurement outcome of all the qubits repeatedly
	 *
	 * This function is useful if at the end of a circuit all the qubits need to be
	 * measured (resulting in a computational basis state as an outcome), 
	 * and the process repeated many times. This function generates measurement samples 
	 * from the probability distribution of all the qubits, and does not modify the
	 * state vector.  
	 *
	 * The source of randomness is a generator object, which provides a means to 
	 * obtain repeatable results, or otherwise control the source of randomness for the
	 * function. By default, qsl::gen is used as the source of randomness. 
	 * 
	 * \param samples Number of samples to take
	 * \param g A generator used as the source of randomness
	 * \return A map of the measurement outcomes where the keys are the computational
	 *         basis state outcomes (encoded as a bitstring) and the associated values
	 *         are the number of times that outcome was measured.
	 */
	std::map<std::size_t, std::size_t> sample_all(std::size_t samples,
						      std::uniform_random_bit_generator auto g = gen) const;
    };

    /**
     * \brief Print the full state vector of a simulator object
     *
     * Print the state vector to an output stream, with no trailing newline.
     * Each amplitude is printed on a separate line.
     *
     * This function behaves differently to s.print(), in that it prints the
     * whole state vector (does not truncate the middle), and also does not 
     * print other information about the simulator.
     *
     * \param os The output stream to print the vector to
     * \param s The simulator whose state is to be printed.
     *
     * Testing:
     * - Check that the function prints the correct state vector to a string
     *   stream, for a number of example simulators (and vectors? see below).
     *   Need to think of a way to check the format properly.
     * 
     * Todo:
     * - This function will match std::vector too probably -- can the user still
     *   override it with their own implmentation if they want? (In the global
     *   namespace?) What happens if their one is in a namespace? Might be better
     *   to restrict this to only simulators, since QSL is not really supposed
     *   to tamper with std::vector.
     */
    template<state_vector S>
    std::ostream & operator << (std::ostream & os, const S & s);

        /**
     * \brief Calculate the Fubini-Study metric between two simulators/state vectors
     *
     * The Fubini-Study metric is a distance between rays in complex projective
     * space which can be interpreted as the angle between rays. This is an 
     * appropriate distance to use for physical states because it compares them 
     * up to a global scale factor (e.g. ignoring any global phases).
     *
     * The definition of the Fubini-Study distance is
     *
     * \f[
     * \gamma (u ,v) = \arccos {\sqrt{\frac{\langle u \vert v \rangle \;
     * \langle v \vert u \rangle }{\langle u \vert u \rangle \;\langle 
     * v \vert v \rangle }}}
     * = \arccos\left[{\frac{|\langle u|v \rangle|}{\lVert u\rVert\lVert 
     * v\rVert}}\right]
     * \f]
     *
     * For more information, see 
     * <a href="https://en.wikipedia.org/wiki/Fubini%E2%80%93Study_metric">
     * the wikipedia page
     * </a>
     *
     * You can use this function between two simulator objects, or a simulator
     * object and a state_vector (meaning it has .size() and operator[] that 
     * returns a real or complex number). If the simulator argument has debugging
     * enabled, then std::invalid_argument is thrown in the following two 
     * scenarios: firstly, if the two arguments are different length state vectors;
     * secondly, if either argument is the zero vector (i.e. not normalisable).
     * If debugging is disabled, then the behaviour in these cases is unspecified.
     *
     * \param u The first state vector to compare
     * \param v The second state vector to compare
     * \return The (real) Fubini-Study distance between u and v
     *
     * Testing: 
     * - Check that the right exceptions are thrown when debugging is enabled
     * - Check that passing in the same vector results in distance zero
     * - Check that passing orthogonal states results in M_PI/2 (todo check this)
     * - Check that the scaling/global phases of the state vectors do not matter
     * - Check all the above with both simulators and vectors in any order
     */
    template<state_vector S1, state_vector S2>
    requires same_precision<S1, S2> && (debug_state_vector<S1> || debug_state_vector<S2>)
    get_precision_t<S1> distance(const S1 & u, const S2 & v);

    /**
     * \brief Calculate the fidelity between two simulators/state vectors
     *
     * The fidelity is a common way to measure the distance between two quantum
     * states. It is not a distance, but relates more closely to the ability to
     * distinguish the two states by measuring them. A fidelity of 1 indicates that
     * the states are equal (corresponds to distance(u,v) = 0), whereas a fidelity
     * 0 indicates that the states are orthogonal.
     *
     * The fidelity between two (normalised) state vectors is defined as 
     *
     * \f[
     * F(u, v) = \lvert \langle u|  v \rangle \rvert^2
     * \f]
     *
     * For more information, see 
     * <a href="https://en.wikipedia.org/wiki/Fidelity_of_quantum_states"> the
     * wikipedia page </a>.
     *
     * You can use this function between two simulator objects, or a simulator
     * object and a state_vector (meaning it has .size() and operator[] that 
     * returns a real or complex number). If the simulator argument has debugging
     * enabled, then std::invalid_argument if the two arguments are different 
     * length state vectors. If debugging is disabled, then the behaviour in this 
     * case is unspecified.
     *
     * This function also normalises the states you pass in, so that the fidelity 
     * calculation is valid.
     *
     * \param u The first state vector to compare
     * \param v The second state vector to compare
     * \return The (real) fidelity between u and v.
     *
     * Testing: 
     * - Check that the right exceptions are thrown when debugging is enabled
     * - Check that passing in the same vector results in fidelity 1
     * - Check that passing orthogonal states results in 0
     * - Check that the scaling/global phases of the state vectors do not matter.
     * - Check all the above with both simulators and vectors in any order
     *
     */
    template<state_vector S1, state_vector S2>
    requires same_precision<S1, S2> && (debug_state_vector<S1> || debug_state_vector<S2>)
    get_precision_t<S1> fidelity(const S1 & u, const S2 & v);

    /**
     * \brief Compute the inner product between two simulators/state vectors
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
     * \f$ \mathbb{R}^N \f$. The state vectors are not normalised in
     * the calculation of this function, which distinguishes it from
     * qsl::fidelity (which normalises the input state vectors) and
     * qsl::distance (for which normalisation does not matter).
     *
     * You can use this function between two simulator objects, or a simulator
     * object and a state_vector (meaning it has .size() and operator[] that 
     * returns a real or complex number). If the simulator argument has debugging
     * enabled, then std::invalid_argument if the two arguments are different 
     * length state vectors. If debugging is disabled, then the behaviour in this 
     * case is unspecified.
     *
     * \param u The first state vector to compare
     * \param v The second state vector to compare
     * \return The (complex) inner product between u and v. 
     *
     * Testing: 
     * - Check that the right exceptions are thrown when debugging is enabled
     * - Check that passing orthogonal states results in 0 (todo check this)
     * - Check that the inner product between equal unit vectors is 1
     * - Check that the outputs from some simulators are normalised (checks sims too).
     * - Check all the above with both simulators and vectors in any order 
     *
     */
    template<state_vector S1, state_vector S2>
    requires same_precision<S1, S2> && (debug_state_vector<S1> || debug_state_vector<S2>)
    std::complex<get_precision_t<S1>> inner_prod(const S1 & u, const S2 & v);    


}

