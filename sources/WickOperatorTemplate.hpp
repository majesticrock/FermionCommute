#pragma once

#include "Operator.hpp"
#include "WickOperator.hpp"
#include <optional>

namespace SymbolicOperators{
    struct WickOperatorTemplate {
        OperatorType type;
        bool is_sc_type{};

        Momentum momentum_difference;
        IndexWrapper allowed_left;
        IndexWrapper allowed_right;

        // Returns the corresponding WickOperator if construction is possible
        // Otherwise, it returns an empty optional
        std::optional<WickOperator> createFromOperators(const Operator& left, const Operator& right);
    
        private:
        std::optional<WickOperator> _handle_sc_type(const Operator& left, const Operator& right);
        std::optional<WickOperator> _handle_num_type(const Operator& left, const Operator& right);
        
    };
}