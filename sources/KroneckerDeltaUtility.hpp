#pragma once
#include "KroneckerDelta.hpp"
#include "Momentum.hpp"
#include "IndexWrapper.hpp"
#include <utility>
#include <type_traits>

template <class T>
class is_linearly_combinable{
    struct _true {};
    struct _false {};
    
    template <class _internal>
    static _true test( decltype(&_internal::operator+=) );
    template <class _internal>
    static _false test(...);
    
    public:
    static constexpr bool value = std::is_same_v<decltype(test<T>(0)), _true>;
};

struct has_plus{
    has_plus& operator+=(const has_plus& rhs){
      return *this;
    };
};