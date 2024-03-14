#pragma once
#include "KroneckerDelta.hpp"
#include "Momentum.hpp"
#include "IndexWrapper.hpp"
#include "../../Utility/sources/defines_arithmetic_operators.hpp"

namespace SymbolicOperators {
	template <class T>
	bool reduce_deltas(std::vector<KroneckerDelta<T>>& deltas) {
		if constexpr (Utility::is_linearly_combinable_v<T>()) {
			for (auto& delta : deltas) {
				delta.first -= delta.second;
				delta.second = T{};
			}
		}
	}
}