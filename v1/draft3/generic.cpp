#include <vector>
#include <stdexcept>
#include "typelist.hpp"


/**
 * \brief Base class containing state, used by everything
 */
template<std::floating_point Fp>
class Base
{
protected:
    unsigned num_qubits{0};
    std::size_t dim{0};
    std::vector<Fp> state{};

public:
    Base() = default;
    Base(const unsigned num_qubits_in)
	: num_qubits{num_qubits_in}, dim{1UL << num_qubits}, state(dim)
	{}
    [[nodiscard]] unsigned numQubits() const { return num_qubits; }
    [[nodiscard]] unsigned dimension() const { return dim; }
    void print() const {
	std::cout << num_qubits << " qubits" << std::endl;
	for (std::size_t n = 0; n < dim; n++) {
	    std::cout << n << ": " << state[n] << std::endl;
	}
    }
};

template<std::floating_point Fp, bool debug>
class SequentialStateSetter : virtual private Base<Fp>
{
    using Base<Fp>::num_qubits;
    using Base<Fp>::dim;
    using Base<Fp>::state;
    
public:
    void setState(const std::vector<Fp> & new_state) {
	if constexpr (debug) {
	    if (new_state.size() != dim) {
		throw std::runtime_error("Wrong size state");
	    }
	}
	state = new_state;
    }
};

template<std::floating_point Fp, bool debug>
class OmpStateSetter : virtual private Base<Fp>
{
    using Base<Fp>::num_qubits;
    using Base<Fp>::dim;
    using Base<Fp>::state;
    
public:
    void setState(const std::vector<Fp> & new_state) {
	if constexpr (debug) {
	    if (new_state.size() != dim) {
		throw std::runtime_error("Wrong size state");
	    }
	}
#pragma omp parallel for
	for (std::size_t n = 0; n < dim; n++) {
	    state[n] = new_state[n];
	}
    }
};



template<std::floating_point Fp, typename Debug, typename SeqPar>
class GenericClass : virtual public Base<Fp>,
		     public SeqPar::StateSetterPolicy<Fp,Debug::value>
{
public:
    GenericClass(unsigned num_qubits) : Base<Fp>{num_qubits} {}
};


// Interface layer

/** 
 * \brief Debug policy class -- choose debugging or not
 */
template<bool state>
struct Debug
{
    using Default = Debug<false>;
    static constexpr bool value{state};
};

/**
 * \brief Sequential policy class
 *
 * Include everything here which should be used if a sequential object
 * is desired.
 */
struct Sequential
{
    using Default = Sequential;
    template<std::floating_point Fp, bool debug>
    using StateSetterPolicy = SequentialStateSetter<Fp,debug>;
};

/**
 * \brief OpenMP policy class
 *
 * Include policies for OpenMP
 */
struct OpenMP
{
    using Default = Sequential;
    template<std::floating_point Fp, bool debug>
    using StateSetterPolicy = OmpStateSetter<Fp,debug>;
};

/// Most general definition
template<std::floating_point Fp, typename TL_in, typename TL_parsed>
struct GenericBuilder;

/**
 * \brief A class whose purpose is to deduce a parameter pack from a typelist
 *
 * A parameter pack can only be deduced, so this class deduces the parameter
 * pack inside a typelist.
 */
template<std::floating_point Fp, typename... Args, typename... ParsedArgs>
struct GenericBuilder<Fp, TypeList<Args...>, TypeList<ParsedArgs...>>
{
    using type = GenericClass<Fp, ParsedArgs...>;
};


// Final level
template<std::floating_point Fp, typename... Args>
using Generic = GenericBuilder<Fp, TypeList<Args...>,
			       typename ParseTypeList<TypeList<Args...>,
						      TypeList<Debug<false>, Sequential>
						      >::next>::type;

int main()
{
    //              { Debug<false>, Sequential}
    Generic<double, Sequential, Debug<true>> g{2};
    g.setState({1,0,0,0,2});
    g.print();
}
