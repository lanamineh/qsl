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
    void print() const {
	std::cout << "Using sequential state setter" << std::endl;
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
    void print() const {
	std::cout << "Using OpenMP state setter" << std::endl;
    }
};

template<std::floating_point Fp, typename Debug, typename SeqPar>
class GenericClass : virtual public Base<Fp>,
		     public SeqPar::StateSetterPolicy<Fp,Debug::value>
{
public:
    GenericClass(unsigned num_qubits) : Base<Fp>{num_qubits} {}
    void print() const {
	using T = typename SeqPar::StateSetterPolicy<Fp,Debug::value>;
	T::print();
	Debug::print();
	Base<Fp>::print();
    }
};

// Interface layer

/** 
 * \brief Debug policy class -- choose debugging or not
 */
template<bool state>
struct Debugging
{
    static constexpr bool value{state};
};

struct NoDebug : Debugging<false>
{
    using Default = NoDebug;
    static void print() {
	std::cout << "Debugging is disabled" << std::endl;
    }
};

struct Debug : Debugging<true>
{
    using Default = NoDebug;
    static void print() {
	std::cout << "Debugging is enabled" << std::endl;
    }
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
template<std::floating_point Fp, typename TL_parsed>
struct GenericBuilder;

/**
 * \brief A class whose purpose is to deduce a parameter pack from a typelist
 *
 * A parameter pack can only be deduced, so this class deduces the parameter
 * pack inside a typelist.
 *
 * \todo You don't need this because there won't ever be an unknown collection
 * of arguments that need inheriting, as is written here. So this whole level
 * of indirection can be omitted.
 */
template<std::floating_point Fp, typename... ParsedArgs>
struct GenericBuilder<Fp, TypeList<ParsedArgs...>>
{
    using type = GenericClass<Fp, ParsedArgs...>;
};

using DefaultList = TypeList<NoDebug, Sequential>;

template<typename... Args>
struct Opts : TypeList<Args...>
{
    using parse = ParseTypeList<TypeList<Args...>, DefaultList>::next; 
};

// Final level
template<std::floating_point Fp, typename... Args>
using Generic = GenericBuilder<Fp, typename Opts<Args...>::parse>::type;

int main()
{
    //Generic<double> g{2};
    //Generic<double, Sequential, NoDebug> g{2};
    //Generic<double, NoDebug, Sequential> g{2};
    //Generic<double, Debug> g{2};
    Generic<double, OpenMP, Debug> g{2};
    
    g.setState({1,0,0,0});
    g.print();
}
