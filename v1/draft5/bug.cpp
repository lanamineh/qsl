#include <concepts>

/////////////////// FAILS with seg fault in GCC, gives error message in clang++
/*
template<template<std::floating_point,std::floating_point> typename... Args>
struct Generic;

template<template<typename...> typename First>
struct Generic<First>
{};
*/

/////////////////// FAILS with incorrect behaviour (?) in gcc,
// compiles fine  with clang++
template<template<std::floating_point...> typename... Args>
struct Generic;

template<template<std::floating_point...> typename First>
struct Generic<First>
{};

/////////////////// FAILS with incorrect behaviour (?) in gcc,
// compiles fine  with clang++
/*
template<template<typename...> typename... Args>
struct Generic;

template<template<typename...> typename First>
struct Generic<First>
{};
*/

/////////////////// FAILS
/*
template<template<std::floating_point,std::floating_point> typename... Args>
struct Generic;

template<template<std::floating_point...> typename First>
struct Generic<First>
{};
*/

/////////////////// FAILS
/*
template<template<std::floating_point,std::floating_point> typename... Args>
struct Generic;

template<template<std::floating_point...> typename First>
struct Generic<First>
{};
*/

/////////////////// WORKS
/*
template<template<typename, typename> typename... Args>
struct Generic;

template<template<typename...> typename First>
struct Generic<First>
{};
*/

/////////////////// WORKS
/*
template<template<typename,typename> typename... Args>
struct Generic;

template<template<typename...> typename First>
struct Generic<First>
{};
*/

/////////////////// WORKS
/*
template<template<std::floating_point> typename... Args>
struct Generic;

template<template<typename> typename First>
struct Generic<First>
{};
*/

int main()
{}
