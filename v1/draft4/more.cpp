#include <iostream>

template<typename NumType>
struct B
{
    const NumType val = 3.2;
    B(NumType val_in) : val{ val_in} {}
};

template<template<typename> typename T>
struct A
{
    T<long int> t_int = 2.3;
    T<float> t_double = 2.3;
    
    void print() const {
	std::cout << t_int.val << std::endl;
	std::cout << t_double.val << std::endl;
    }
};


int main ()
{
    A<B> a;
    a.print();
}
