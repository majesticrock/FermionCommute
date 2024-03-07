#pragma once
#include <iostream>
#include "Momentum.hpp"

namespace SymbolicOperators {
    enum OperatorType { Number_Type, CDW_Type, SC_Type, Eta_Type, Undefined_Type };

	inline std::ostream& operator<<(std::ostream& os, const OperatorType op) {
		switch (op) {
		case SC_Type:
			os << "f";
			break;
		case Eta_Type:
			os << "\\eta";
			break;
		case CDW_Type:
			os << "g";
			break;
		case Number_Type:
			os << "n";
			break;
		default:
			os << "ERROR_OPERATOR";
		}
		return os;
	};

	typedef std::pair<Momentum, Momentum> pair_of_momenta;
	struct WickOperator {
		OperatorType type;
		bool isDaggered;
		Momentum momentum;
		std::vector<std::string> indizes;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& type;
			ar& isDaggered;
			ar& momentum;
			ar& indizes;
		}

		WickOperator(const OperatorType& _type, const bool _isDaggered, const Momentum& _momentum, const std::vector<std::string>& _indizes = std::vector<std::string>());
		WickOperator(const OperatorType& _type, const bool _isDaggered, const Momentum& _momentum, const std::string& _index);
		WickOperator();

		inline bool usesIndex(const std::string& index) const noexcept {
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