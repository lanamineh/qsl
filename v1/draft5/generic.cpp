#include <vector>
#include <iostream>

template<std::floating_point Fp>
class Base
{
protected:
    unsigned num_qubits_{0};
    std::size_t dim_{0};
    std::vector<Fp> state_{};
public:
    Base(unsigned num_qubits)
	: num_qubits_{num_qubits},
	  dim_{1UL << num_qubits_},
	  state_(dim_)
	{
	    state_[0] = 1;
	}
    [[nodiscard]] unsigned numQubits() const { return num_qubits_; }
    [[nodiscard]] unsigned dimension() const { return dim_; }
};

template<typename Previous, typename Base>
concept Derived = std::derived_from<Previous, Base>;

template<std::floating_point Fp, typename T>
class Thing
{

};

template<std::floating_point Fp, Derived<Base<Fp>> Previous, bool debug>
class DefaultStateSetter : public Previous
{
protected:
    using Previous::state_;
    using Previous::dim_;
public:
    using Previous::Previous;
    void setState(const std::vector<Fp> & state);
};

template<std::floating_point Fp, Derived<Base<Fp>> Previous, bool debug>
void DefaultStateSetter<Fp, Previous, debug>::setState(const std::vector<Fp> & state)
{
    if constexpr (debug == true) {
    	if (state.size() != dim_) {
    	    throw std::logic_error("Incorrect dimension "
    				   "when setting state");
    	}
    }
    state_ = state;
}

template<std::floating_point Fp, Derived<Base<Fp>> Previous>
class DefaultPrinter : public Previous
{
protected:
    using Previous::num_qubits_;
    using Previous::dim_;
    using Previous::state_;
public:
    using Previous::Previous;
    
    void print() const
	{
	    std::cout << num_qubits_ << " qubits" << std::endl;
	    for (const auto & amp : state_) {
		std::cout << amp << std::endl;
	    }
	}
};

template<std::floating_point Fp, Derived<Base<Fp>> Previous>
class IndexedPrinter : public Previous
{
protected:
    using Previous::num_qubits_;
    using Previous::dim_;
    using Previous::state_;
public:
    using Previous::Previous;
    void print() const
	{
	    std::cout << num_qubits_ << " qubits" << std::endl;
	    for (std::size_t n = 0; n < dim_; n++) {
		std::cout << n << ": " << state_[n] << std::endl;
	    }
	}
};

// Unfortunately, adding concepts to the template template parameter's
// template parameters does not seem to work. 
template<std::floating_point Fp,
	 template<typename,typename> typename... Args>
struct GenericWrapper;

template<std::floating_point Fp>
struct GenericWrapper<Fp> : Base<Fp>
{
    using Base<Fp>::Base;
};

template<typename Fp,
	 template<typename,typename> typename First,
	 template<typename,typename> typename... Rest>
struct GenericWrapper<Fp, First, Rest...> : First<Fp, GenericWrapper<Fp, Rest...>>
{
    using  First<Fp, GenericWrapper<Fp, Rest...>>::First;
};

template<template<typename,typename,bool> typename Policy, bool debug>
struct StateSetter
{
    template<std::floating_point Fp, Derived<Base<Fp>> Previous>
    using type = Policy<Fp, Previous, debug>;
};

template<template<typename,typename> typename Policy>
struct Printer
{
    template<std::floating_point Fp, Derived<Base<Fp>> Previous>
    using type = Policy<Fp, Previous>;
};

template<std::floating_point Fp, typename... Args>
using Generic = GenericWrapper<Fp, Args::template type...>;
					  
int main()
{
    //DefaultStateSetter<double, DefaultPrinter<double, Base<double>>> def;
    //def.setState({1,0,0});
    //def.print();

    Generic<double,
	    Printer<DefaultPrinter>,
	    StateSetter<DefaultStateSetter, false>> base{3};
    base.setState({1});
    base.print();
    
}
