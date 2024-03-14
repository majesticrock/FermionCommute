#include "WickOperatorTemplate.hpp"
#include <algorithm>
#include <iterator>
#include "Momentum.hpp"
#include "KroneckerDelta.hpp"

namespace SymbolicOperators{
    std::optional<WickOperator> WickOperatorTemplate::_handle_sc_type(const Operator& left, const Operator& right){
        // c_{-k} c_{k} or c_{k}^+ c_{-k}^+
        const Momentum& base{ left.isDaggered ? left.momentum : right.momentum };
        Momentum momentum_diff = left.momentum + right.momentum;

        KroneckerDelta<Momentum> momentum_delta{this->momentum_difference, momentum_diff};
        std::vector<KroneckerDelta<Index>> index_delta;
        index_delta.reserve(this->allowed_left.size() + this->allowed_right.size());
        std::transform(this->allowed_left.begin(), this->allowed_left.end(), left.indizes.begin(), std::back_inserter(index_delta),
            [](const Index& left, const Index& right) { 
                return make_delta(left, right);
            });
        std::transform(this->allowed_right.begin(), this->allowed_right.end(), right.indizes.begin(), std::back_inserter(index_delta),
            [](const Index& left, const Index& right) { 
                return make_delta(left, right);
            });
    }
    std::optional<WickOperator> WickOperatorTemplate::_handle_num_type(const Operator& left, const Operator& right){
        // c_{k}^+ c_{k}
        Momentum momentum_diff = left.momentum - right.momentum;

    }

    std::optional<WickOperator> WickOperatorTemplate::createFromOperators(const Operator& left, const Operator& right){
        if(this->is_sc_type){
            if (left.isDaggered != right.isDaggered)
                return std::nullopt;
            return this->_handle_sc_type(left, right);
        }

        if(left.isDaggered == right.isDaggered)
            return std::nullopt;
        return this->_handle_num_type(left, right);
    }
}