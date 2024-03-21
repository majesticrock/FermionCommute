#pragma once
#include "KroneckerDelta.hpp"
#include "Momentum.hpp"
#include "IndexWrapper.hpp"
#include "../../Utility/sources/defines_arithmetic_operators.hpp"

namespace SymbolicOperators {
	template <class T>
	void remove_delta_squared(std::vector<KroneckerDelta<T>>& deltas) {
		if constexpr (Utility::is_linearly_combinable_v<T>()) {
			for (auto& delta : deltas) {
				delta.first -= delta.second;
				delta.second = T{};
			}
		}
		// Remove duplicates as delta^2 = delta
		for (size_t i = 0U; i < deltas.size(); ++i) {
			for (auto jt = deltas.begin() + i + 1; jt != deltas.end();) {
				if (deltas[i] == *jt) {
					jt = deltas.erase(jt);
				}
				else {
					++jt;
				}
			}
		}
	}

	template <class T>
	void remove_delta_is_one(std::vector<KroneckerDelta<T>>& deltas) {
		for(auto it = deltas.begin(); it != deltas.end();){
			if(it->isOne()){
				it = deltas.erase(it);
			}else{
				++it;
			}
		}
	}
}