#pragma once

#include "../Momentum.hpp"
#include "../Operator.hpp"

namespace SymbolicOperators{
    inline const Momentum base_k(momentum_pairs{ {1, 'k'} }, false);

	inline const Operator c_k(base_k, SpinUp, false);
	inline const consttexpr Operator c_minus_k(-base_k, SpinDown, false);

    inline const Operator c_k_dagger(base_k, SpinUp, true);
	inline const Operator c_minus_k_dagger(-base_k, SpinDown, true);
}