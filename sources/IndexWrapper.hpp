#pragma once
#include <iostream>
#include "../../Utility/sources/VectorWrapper.hpp"

namespace SymbolicOperators {
	enum Index { SpinUp = 0, SpinDown, Sigma, SigmaPrime, UndefinedIndex };
	// Returns true if the index represents a variable and false otherwise
	// Example: If the index is SpinUp it is fixed, i.e., non-mutable
	constexpr bool is_mutable(const Index idx) {
		return (idx < 2);
	}

	std::ostream& operator<<(std::ostream& os, const Index index);

	struct IndexWrapper : public Utility::VectorWrapper<Index> {
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->_vector;
		}

		IndexWrapper() = default;
		IndexWrapper(Index _spin) : Utility::VectorWrapper<Index>(1U, _spin) {};
		IndexWrapper(const std::vector<Index>& _indizes)
			: Utility::VectorWrapper<Index>(_indizes) {};
		IndexWrapper(std::vector<Index>&& _indizes)
			: Utility::VectorWrapper<Index>(std::move(_indizes)) {};
	};

	std::ostream& operator<<(std::ostream& os, const IndexWrapper& indizes);
}