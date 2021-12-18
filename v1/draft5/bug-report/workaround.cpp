#include <iostream>

// General template
template<template<typename> typename... T>
struct X
{
    static void print() { std::cout << "General" << std::endl; }
};

// Specialisation
template<template<typename> typename T>
struct X<T>
{
    static void print() { std::cout << "Special" << std::endl; }
};

template<std::floating_point Fp> struct A {};
template<std::floating_point Fp>struct B {};
template<std::integral I> struct R {};

int main()
{
    X<A,B>::print();
    X<A>::print(); // Use of partial specialisation 
    //X<R,B>::print(); I want this to fail due to wrong concept in R
}
