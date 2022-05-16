/**
 * \file qsl-testing.cpp
 *
 * Contains gtest code.
 */

#include <gtest/gtest.h>

#include "meta.hpp"
#include "qsl.hpp"
#include <armadillo>

/**
 * \brief Construct a std::tuple of simulator objects from a std::tuple of std::tuple
 *        simulator template arguments.
 *
 * You can pass in the type of simulator you want to create e.g. qsl::basic, the function
 * will generate qsl::basic simulator types for the template arguments provided in the tuple T.
 */
template <template<typename, bool, typename> typename S, typename T>
struct make_sims {};

template <template<typename, bool, typename> typename S, typename... F, typename... D, typename... P>
struct make_sims<S, std::tuple<std::tuple<F, D, P>...>>
{
    using type = std::tuple<S<F, D::value, P>...>;
};

template <typename... T>
using make_sims_t = make_sims<T...>::type;



template <typename T>
class ConstructorTests : public testing::Test {};

