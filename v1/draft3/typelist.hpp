#include <iostream>

template<typename... Others>
struct TypeList
{};

/// TypeList<int,double,double>

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
