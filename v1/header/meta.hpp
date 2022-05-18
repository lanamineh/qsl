/**
 * \file meta.hpp
 *
 * Conatins useful meta-programming utilities.
 */

#ifndef QSL_META_HPP
#define QSL_META_HPP

namespace qsl
{
    /**
     * \brief Concatenate the types of multiple tuples.
     *
     * You can input as many tuples as you want as the template parameters,
     * this meta-function will combine the types of each tuple (in order) into
     * a bigger tuple.
     */
    template <typename T, typename...>
    struct cat
    {
	using type = T;
    };

    template <typename... T1, typename... T2, typename... T3>
    struct cat<std::tuple<T1...>, std::tuple<T2...>, T3...>
    {
	using type = cat<std::tuple<T1..., T2...>, T3...>::type;
    };

    template <typename... T>
    using cat_t = cat<T...>::type;

    /**
     * \brief Prepend a type T to the beginning of every tuple in a tuple of tuples.
     *
     * For example, if the arguments are <T, std::tuple<std::tuple<A, B>, std::tuple<C>>>
     * the resulting type will be std::tuple<std::tuple<T, A, B>, std::tuple<T, C>>.
     */
    template <typename T, typename L>
    struct prepend{};

    template <typename T, typename... L>
    struct prepend<T, std::tuple<L...>>
    {
	using type = std::tuple<cat_t<std::tuple<T>, L>...>;
    };

    template <typename T, typename L>
    using prepend_t = prepend<T, L>::type;

    /**
     * \brief Construct all combinations of one type extracted from each tuple in a list.
     *
     * The input is an arbitrary number of std::tuples where you wish to take one type from
     * each tuple and construct all possible combinations. For example, if the input is
     * <std::tuple<A, B>, std::tuple<C, D>>, the resulting type will be 
     * std::tuple<std::tuple<A, C>, std::tuple<A, D>, std::tuple<B, C>, std::tuple<B, D>>.
     *
     * For qsl, we could modify this to output a std::tuple of all the simulator objects if we wanted.
     */
    template <typename... T>
    struct comb;

    template <typename... T>
    struct comb<std::tuple<T...>>
    {
	// At the bottom level of the recursion, the input std::tuple is split into
	// a std::tuple of separate std::tuples for each type in the original tuple.
	using type = std::tuple<std::tuple<T>...>;
    };

    template <typename... T1, typename... T2>
    struct comb<std::tuple<T1...>, T2...>
    {
	// Construct the combination of all the remaining std::tuples
	using combined = comb<T2...>::type;

	// Prepend each type in the first tuple T1 to every combination
	using type = cat_t<prepend_t<T1, combined>...>;
    };

    template <typename... T>
    using comb_t = comb<T...>::type;

    /**
     * \brief Construct a std::tuple of simulator objects from a std::tuple of std::tuple
     *        simulator template arguments.
     *
     * You can pass in the type of simulator you want to create e.g. qsl::basic, the function
     * will generate qsl::basic simulator types for the template arguments provided in the tuple T.
     */
    template <template<typename, bool, typename> typename S, typename T>
    struct make_sims;

    template <template<typename, bool, typename> typename S, typename... F, typename... D, typename... P>
    struct make_sims<S, std::tuple<std::tuple<F, D, P>...>>
    {
	using type = std::tuple<S<F, D::value, P>...>;
    };

    template <template<typename, bool, typename> typename S, typename T>
    using make_sims_t = make_sims<S, T>::type;

    /**
     * \brief Extract tuple contents into gtest.
     */
    template <typename T>
    struct gtest_typelist;

    template <typename... T>
    struct gtest_typelist<std::tuple<T...>>
    {
	using type = ::testing::Types<T...>;
    };

    template <typename T>
    using gtest_typelist_t = gtest_typelist<T>::type;

    /**
     * \brief Extract the debug value from the simulator
     */
    template <typename S>
    constexpr bool is_debug;

    template <template<typename, bool, typename> typename S, typename F, bool D, typename P>
    constexpr bool is_debug<S<F, D, P>> = D;
    
}

#endif 
