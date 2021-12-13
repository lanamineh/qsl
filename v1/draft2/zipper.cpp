// Forward declarations
template < typename... Args >
struct zipper;

// Two types meaning true and false
struct true_type {};
struct false_type {};

// The only purpose of this struct is to be associated with Types...
template < typename... Types >
struct typelist {};


// ===================================================
// is_in < type, typelist<...> >::type
//     is true_type if type is in typelist
//     is false_type if type is not in typelist

// Assume TElement is not in the list unless proven otherwise
template < typename TElement, typename TList >
struct is_in {
  typedef false_type type;
};

// If it matches the first type, it is definitely in the list
template < typename TElement, typename... TTail >
struct is_in < TElement, typelist < TElement, TTail... > >
{
  typedef true_type type;
};

// If it is not the first element, check the remaining list
template < typename TElement, typename THead, typename... TTail >
struct is_in < TElement, typelist < THead, TTail... > >
{
  typedef typename is_in < TElement, typelist < TTail... > >::type type;
};

// ===================================================
// add_unique < TNew, typelist<...> >::type
//     is typelist < TNew, ... > if TNew is not already in the list
//     is typelist <...> otherwise

// Append a type to a type_list unless it already exists
template < typename TNew, typename TList,
  typename Tis_duplicate = typename is_in < TNew, TList >::type
  >
struct add_unique;

// If TNew is in the list, return the list unmodified
template < typename TNew, typename... TList >
struct add_unique < TNew, typelist < TList... >, true_type >
{
  typedef typelist < TList... > type;
};

// If TNew is not in the list, append it
template < typename TNew, typename... TList >
struct add_unique < TNew, typelist < TList... >, false_type >
{
  typedef typelist < TNew, TList... > type;
};

// ===================================================
// process_zipper_arguments < Args... >::type
//     returns a typelist of types to be inherited from.
//
// It performs the following actions:
// a) Unpack zipper<...> and typelist <...> arguments
// b) Ignore values that are already in the list

template < typename... Args >
struct process_zipper_arguments;

// Unpack a zipper in the first argument
template < typename... ZipperArgs, typename... Args >
struct process_zipper_arguments < zipper < ZipperArgs... >, Args... >
{
  typedef typename process_zipper_arguments < ZipperArgs..., Args... >::type type;
};

// Unpack a typelist in the first argument
template < typename... TypeListArgs, typename... Args >
struct process_zipper_arguments < typelist < TypeListArgs... >, Args... >
{
  typedef typename process_zipper_arguments < TypeListArgs..., Args... >::type type;
};

// End the recursion if the list is empty
template < >
struct process_zipper_arguments < >
{
  typedef typelist < > type;
};

// Construct the list of unique types by appending them one by one
template < typename THead, typename... TTail >
struct process_zipper_arguments < THead, TTail... >
{
  typedef typename
    add_unique < THead,
      typename process_zipper_arguments < TTail... >::type
    >::type type;
};


// ===================================================
// The zipper class that you might want


// If the list of types is not yet known, process it.
// The inheritance is ugly, but there is a workaround
template < typename... Args >
struct zipper : zipper < typename process_zipper_arguments < Args... >::type >
{
  // // Instead of inheriting, you can use zipper as a factory.
  // // So this:
  // typedef zipper < meas2, zipper < meas1, meas > > mymeas;
  // // Turns to:
  // typedef typename zipper < meas2, zipper < meas1, meas > >::type mymeas;
  typedef zipper < typename process_zipper_arguments < Args... >::type > type;
};

// If the list of types is known, inherit from each type
template <typename... Args>
struct zipper<typelist< Args... >> : public virtual Args...
{};

// ===================================================
// Short usage demo, replace with your own code

struct meas {
    int i;
};

struct meas2 {
    int j;
};

struct meas3 {
    int k;
};


using meas_type = zipper<meas, meas, meas3>;
using meas_type_2 = zipper<meas2, meas_type, meas2>;

using  nicer_meas_type2 = typename zipper <meas_type2>::type;


int main ( int, char** )
{
    meas * m = new meas_type2;
    meas_type2 n;
    nicer_meas_type2 o;

    return 0;
}
