#include <concepts>
#include <iostream>

// Template definition and partial specialisation (works)
//template<template<typename> typename...> class A;
//template<template<typename> typename T> class A<T>;

// Template definition and partial specialisation with concept (fails)

// Compiles with gcc-10 but not gcc-11
//template<template<std::floating_point> typename...> class A;
//template<template<std::floating_point> typename T> class A<T>;

// Fails with gcc-10 and gcc-11
// template<template<std::floating_point...> typename...> class A;
// template<template<std::floating_point...> typename T> class A<T>;


// template<template<std::floating_point...> typename...> class A;
// template<template<std::floating_point...> typename T> class A<T>;

// Checks for correct behaviour below
template<template<std::floating_point> typename... T>
struct X
{
    static void print() { std::cout << "General" << std::endl; }
};

template<template<std::floating_point> typename T>
struct X<T>
{
    static void print() { std::cout << "Special" << std::endl; }
};


template<std::floating_point Fp>
struct A
{
    using type = Fp;
};

template<std::floating_point Fp>
struct B
{
    using type = Fp;
};


template<std::integral I>
struct R
{
    using type = I;
};


int main()
{
    X<A,B>::print();
    //X<R,B>::print(); Constraint failure
    X<A>::print(); // Partial specialisation 

}
