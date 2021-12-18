#include <concepts>

// Template definition and partial specialisation (works)
//template<template<typename> typename...> class A;
//template<template<typename> typename T> class A<T>;

template<template<std::floating_point> typename...> class A;
template<template<std::floating_point> typename T> class A<T>;
  
int main()
{

}
