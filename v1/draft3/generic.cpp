#include <vector>
#include <iostream>

/**
 * \brief Base class for quantum simulators
 *
 * Contains only methods which do not depend on any template
 * parameter other than the floating-point precision. The class
 * is augmented with other classes which specify other parts of
 * the behaviour of the simulation.
 *
 */
template<std::floating_point Fp>
class Base
{
protected:
    unsigned num_qubits = 0;
    std::size_t dim = 0;
    std::vector<Fp> state{};
    
public:

    /**
     * \brief Default constructor to make the virtual inheritance happy...
     *
     * ...without going to the effort of writing all the non-default
     * constructors in the derived classes, if that is what is needed.
     * Need to figure out exactly what is best to do here.
     */
    Base() = default;
    Base(unsigned num_qubits_in)
	: num_qubits{num_qubits_in}, dim{ 1UL << num_qubits }, state(dim)
	{}
};

/**
 * \brief StateSetter base class to provide generic functions
 *
 * This class contains the code for checking the input state. It 
 * is used as a base class for the StateSetter implementations so
 * that they do not have to duplicate the checking code.
 *
 */
template<std::floating_point Fp>
struct StateSetter : virtual private Base<Fp>
{
    void checkInputState(const std::vector<Fp> & new_state) const {

	// Check the state size
	if (new_state.size() != Base<Fp>::state.size()) {
	    throw std::logic_error("Cannot set state: wrong dimension");
	}

	// Check the input state is normalised? etc.
    }
};

template<bool state> struct Debugging
{
    using Default = Debugging<false>; 
    static constexpr bool value = state;
};

/**
 * \brief Set the state to a particular value
 *
 * This class sets the state sequentially (without OpenMP). A bool is used
 * to specify whether debugging checks should be performed.
 *
 */
template<std::floating_point Fp, typename Debug>
struct SequentialStateSetter : virtual private Base<Fp>, public StateSetter<Fp>
{
    void setState(const std::vector<Fp> & new_state)
    	{
	    if constexpr (Debug::value == true) {
		StateSetter<Fp>::checkInputState(new_state);
	    }
	    Base<Fp>::state = new_state;
    	}
};

/**
 * \brief Same set state function but uses OpenMP
 *
 */
template<std::floating_point Fp, typename Debug>
class OmpStateSetter : virtual private Base<Fp>, public StateSetter<Fp>
{
public:
    void setState(const std::vector<Fp> & new_state)
    	{
	    if constexpr (Debug::value == true) {
		StateSetter<Fp>::checkInputState(new_state);		
	    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
	    for (std::size_t index = 0; index < Base<Fp>::dim; index++) {
		Base<Fp>::state[index] = new_state[index];
	    }
	}
};

/**
 * \brief Set the state to a particular value
 *
 * This class sets the state sequentially (without OpenMP). A bool is used
 * to specify whether debugging checks should be performed.
 *
 */
template<std::floating_point Fp>
struct SequentialStateResetter : virtual private Base<Fp>
{
    void reset()
    	{
	    for (auto & amp : Base<Fp>::state) {
		amp = 0;
	    }
	    Base<Fp>::state[0] = 1;
    	}
};

/// OMP state resetter
template<std::floating_point Fp>
class OmpStateResetter : virtual private Base<Fp>
{
public:
    void reset()
    	{
#ifdef _OPENMP
#pragma omp parallel for
#endif
	    for (auto & amp : Base<Fp>::state) {
		amp = 0;
	    }
	    Base<Fp>::state[0] = 1;
	}
};

template<std::floating_point Fp, typename Debug,
	 class StateSetterPolicy,
	 class StateResetterPolicy>
class BetterBase : virtual private Base<Fp>,
		   public StateSetterPolicy,
		   public StateResetterPolicy 
{
#ifndef _OPENMP
    static_assert(not std::is_same_v<StateSetterPolicy<Fp,Debug>,
		  OmpStateSetter<Fp,Debug>>,
		  "You have requested an OpenMP implementation, but _OPENMP "
		  "is not defined. Did you forget -fopenmp?");
#endif

public:
    std::size_t dimension() const { return Base<Fp>::dim; }
    void print() const
	{
	    std::cout << "State on " << Base<Fp>::num_qubits
		      << " qubits:" << std::endl;
	    for (std::size_t index = 0; index < Base<Fp>::dim; index++) {
		std::cout << index << ": " << Base<Fp>::state[index] << std::endl;
	    }
	}  

    
};

/**
 * \brief A collection of all the settings required for sequential operation
 *
 */
struct Sequential
{
    template<std::floating_point Fp, typename Debug>
    using StateSetterPolicy = SequentialStateSetter<Fp, Debug>;

    template<std::floating_point Fp>
    using StateResetterPolicy = SequentialStateResetter<Fp>;
};

struct Omp
{
    using Default = Sequential; 

    template<std::floating_point Fp, typename Debug>
    using StateSetterPolicy = OmpStateSetter<Fp, Debug>;

    template<std::floating_point Fp>
    using StateResetterPolicy = OmpStateResetter<Fp>;
};

template<typename... Others>
struct TypeList
{};

/// General case
template<typename First, typename... Others>
struct TypeList<First, Others...>
{
    static void print()
	{
	    std::cout << typeid(First).name() << ",";
	    TypeList<Others...>::print();
	}
};

/// Base case
template<>
struct TypeList<>
{
    static void print()
	{
	    std::cout << std::endl;
	}
};


// Check whether type First is in typenames Others...
template<typename First, typename... Others>
static constexpr bool contains {
    std::disjunction_v<std::is_same<First, Others>...>
};

template<typename TL_in, typename TL_out = TypeList<>>
struct UniqueTypeList
{};

template<typename First, typename... Others, typename... Pruned>
struct UniqueTypeList<TypeList<First, Others...>, TypeList<Pruned...>>
{    
    // Whether to push the new type of not
    static constexpr bool push = contains<First, Pruned...>;

    // The new type including First
    using WithFirst = TypeList<Pruned..., First>; 

    // The new type excluding First
    using WithoutFirst = TypeList<Pruned...>; 
    
    // The new typelist
    using NewTypeList = std::conditional_t<push, WithoutFirst, WithFirst>;
    
    // Recursive traversal of the list
    using next = UniqueTypeList<TypeList<Others...>, NewTypeList>::next;
};

template<typename... Pruned>
struct UniqueTypeList<TypeList<>, TypeList<Pruned...>>
{
    // The new type excluding First
    using next = TypeList<Pruned...>; 
};

template<typename Find, typename Replace,
	 typename TL_in, typename TL_out = TypeList<>>
struct ReplaceType
{};

template<typename Find, typename Replace,
	 typename First, typename... Others, typename... Pruned>
struct ReplaceType<Find, Replace, TypeList<First, Others...>, TypeList<Pruned...>>
{    
    // Check whether the next type First is equal to Find
    static constexpr bool replace = std::is_same_v<Find, First>;

    // The next type after replacing First with Replace
    using WithReplace = TypeList<Pruned..., Replace>; 

    // The type after leaving First alone
    using WithFirst = TypeList<Pruned..., First>; 
    
    // The new typelist
    using NewTypeList = std::conditional_t<replace, WithReplace, WithFirst>;
    
    // Recursive traversal of the list
    using next = ReplaceType<Find, Replace, TypeList<Others...>, NewTypeList>::next;
};

template<typename Find, typename Replace, typename... Pruned>
struct ReplaceType<Find, Replace, TypeList<>, TypeList<Pruned...>>
{
    // The new type excluding First
    using next = TypeList<Pruned...>; 
};

template<typename TL_in, typename TL_out/* = TypeList<A_def, B_def>*/>
struct ParseTypeList
{};

template<typename First, typename... Others, typename... Pruned>
struct ParseTypeList<TypeList<First, Others...>, TypeList<Pruned...>>
{
    // The default type to be replaced
    using Default = First::Default;
    
    // The next type after replacing First with Replace
    using NewTypeList = ReplaceType<Default, First, TypeList<Pruned...>>::next;
    
    // Recursive traversal of the list
    using next = ParseTypeList<TypeList<Others...>, NewTypeList>::next;
};

template<typename... Pruned>
struct ParseTypeList<TypeList<>, TypeList<Pruned...>>
{
    // The new type excluding First
    using next = TypeList<Pruned...>; 
};


template<std::floating_point Fp, typename... Args>
class GenericClass : public BetterBase<Fp, Args...>
{

public:
    GenericClass(unsigned num_qubits_in) : Base<Fp>{num_qubits_in} {}
};


// template<typename... Args>
// class GenericClass : public Args...
// { };

template<typename TL_in, typename TL_parsed>
class GenericBuilder;

template<typename... Args, typename... ParsedArgs>
struct GenericBuilder<TypeList<Args...>, TypeList<ParsedArgs...>>
{
    using type = GenericClass<ParsedArgs...>;
};

template<typename... Args>
using Generic = GenericBuilder<TypeList<Args...>,
			       typename ParseTypeList<TypeList<Args...>,
						      TypeList<A_def, B_def>
						      >::next>::type;

int main()
{
    Generic<double, Omp, Debugging<false>> q{2};
    q.print();
    
    q.setState({1,0,0,1});
    q.print();

    q.reset();
    q.print();

}
