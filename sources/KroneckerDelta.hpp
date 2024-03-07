#pragma once
#include <utility>

namespace SymbolicOperators{
    template<typename T>
    struct KroneckerDelta{
        T first{};
        T second{};

        bool isOne() const {
            return first == second;
        };
    };

    template<typename T>
    inline auto make_delta(const T& first, const T& second){
        return KroneckerDelta<T>(first, second);
    }

    template<typename T>
    inline auto make_delta(T&& first, T&& second){
        return KroneckerDelta<T>(std::move(first), std::move(second));
    }

    template<typename T>
    bool operator==(const KroneckerDelta& lhs, const KroneckerDelta& rhs){
        if(lhs.first == rhs.first && lhs.second == rhs.second) return true;
        if(lhs.first == rhs.second && lhs.second == rhs.first) return true;
        return false;
    };
    template<typename T>
    bool operator!=(const KroneckerDelta& lhs, const KroneckerDelta& rhs){
        return !(lhs == rhs);
    };

    template<typename T>
    inline std::ostream& operator<<(std::ostream& os, KroneckerDelta<T> const& delta){
        os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
        return os;
    };
}