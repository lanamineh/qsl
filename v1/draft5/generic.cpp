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

template<std::floating_point Fp, typename T>
class Thing
{

};

template<std::floating_point Fp, typename Previous>
requires std::derived_from<Previous, Base<Fp>>
class DefaultStateSetter : public Previous
{
protected:
    using Previous::state_;
public:
    using Previous::Previous;
    void setState(const std::vector<Fp> & state) { state_ = state; }
};

template<typename Previous, typename Base>
concept Derived = std::derived_from<Previous, Base>;

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
struct Generic;

template<std::floating_point Fp>
struct Generic<Fp> : Base<Fp>
{
    using Base<Fp>::Base;
};

template<typename Fp,
	 template<typename,typename> typename First,
	 template<typename,typename> typename... Rest>
struct Generic<Fp, First, Rest...> : First<Fp, Generic<Fp, Rest...>>
{
    using  First<Fp, Generic<Fp, Rest...>>::First;
};

// template<std::floating_point Fp, template<typename,typename> typename... Args>
// struct GenLinear;

// template<std::floating_point Fp>
// struct GenLinear<Fp>
// {
//     using type = Base<Fp>;
// };

// template<std::floating_point Fp,
// 	 template<typename,typename> typename First,
// 	 template<typename,typename> typename... Rest>
// struct GenLinear<Fp, First, Rest...>
// {
//     struct inner : First<Fp, typename GenLinear<Fp, Rest...>::type> {};
//     using type = inner;
// };

template<typename a, typename b>
struct AAA
{
};


int main()
{
    //DefaultStateSetter<double, DefaultPrinter<double, Base<double>>> def;
    //def.setState({1,0,0});
    //def.print();

    Generic<double, DefaultPrinter, DefaultStateSetter> base{3};
    base.setState({1});
    base.print();
    
}
