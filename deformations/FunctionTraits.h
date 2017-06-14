#pragma once

#include <boost/mpl/fold.hpp>
#include <boost/function_types/parameter_types.hpp>
#include <type_traits>

namespace {

template < typename Sequence, typename T > struct add_to_tuple;
template < typename T, typename... TCurrent >
struct add_to_tuple< std::tuple<TCurrent...>, T >
{
    using type = std::tuple<TCurrent..., typename std::decay<T>::type>;
};

}

namespace details {

template < typename T, bool IsAFunction >
struct FunctionTraitsImpl
{
    using Arguments = typename boost::mpl::fold<
        typename boost::function_types::parameter_types<T>,
        std::tuple<>,
        add_to_tuple<boost::mpl::_1, boost::mpl::_2>
    >::type;
};

template < typename T >
struct FunctionTraitsImpl<T, false>
{
    using Arguments = typename boost::mpl::fold<
        typename boost::mpl::pop_front< boost::function_types::parameter_types< decltype(&T::operator()) > >::type,
        std::tuple<>,
        add_to_tuple<boost::mpl::_1, boost::mpl::_2>
    >::type;
};

template < typename T >
struct FunctionTraitsCleaned
  : FunctionTraitsImpl< T, std::is_function<T>::value >
{};

}

//TODO: not clear that it works in all cases.
template < typename T >
struct FunctionTraits
  : details::FunctionTraitsCleaned< typename std::remove_pointer< typename std::remove_reference<T>::type >::type >
{};
