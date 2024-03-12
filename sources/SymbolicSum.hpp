#pragma once

#include "../../Utility/sources/VectorWrapper.hpp"
#include "IndexWrapper.hpp"

namespace SymbolicOperators {
	template<class SumIndex>
	struct SymbolicSum : public Utility::VectorWrapper<SumIndex>
	{
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->_vector;
		}

		SymbolicSum() = default;
		SymbolicSum(SumIndex sum_index) : Utility::VectorWrapper<SumIndex>(1U, sum_index) {};
		SymbolicSum(const std::vector<SumIndex>& _indizes)
			: Utility::VectorWrapper<SumIndex>(_indizes) {};
		SymbolicSum(std::vector<SumIndex>&& _indizes)
			: Utility::VectorWrapper<SumIndex>(std::move(_indizes)) {};
	};

	template<class SumIndex>
	std::ostream& operator<<(std::ostream& os, SymbolicSum<SumIndex> const& sum) {
		if (sum.empty()) return os;
		os << "\\sum_{ ";
		for (const auto& index : sum) {
			os << index << " ";
		}
		os << "} ";
		return os;
	}

	typedef SymbolicSum<Index> IndexSum;
	typedef SymbolicSum<char> MomentumSum;
}