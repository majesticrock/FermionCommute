#pragma once
#include <iostream>
#include "Momentum.hpp"
#include "IndexWrapper.hpp"

namespace SymbolicOperators {
    enum OperatorType { Number_Type, CDW_Type, SC_Type, Eta_Type, Undefined_Type };

	std::ostream& operator<<(std::ostream& os, const OperatorType op);

	typedef std::pair<Momentum, Momentum> pair_of_momenta;
	struct WickOperator {
		OperatorType type;
		bool isDaggered;
		Momentum momentum;
		IndexWrapper indizes;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& type;
			ar& isDaggered;
			ar& momentum;
			ar& indizes;
		}

		WickOperator(const OperatorType& _type, const bool _isDaggered, const Momentum& _momentum, const IndexWrapper& _indizes = IndexWrapper());
		WickOperator(const OperatorType& _type, const bool _isDaggered, const Momentum& _momentum, const Index _index);
		WickOperator();

		inline bool usesIndex(const Index index) const noexcept {
			for (const auto& idx : this->indizes) {
				if (idx == index) return true;
			}
			return false;
		};
		inline bool dependsOn(char momentum) const noexcept {
			return this->momentum.isUsed(momentum) != -1;
		}
	};

	std::ostream& operator<<(std::ostream& os, const WickOperator& op);
	std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops);
}