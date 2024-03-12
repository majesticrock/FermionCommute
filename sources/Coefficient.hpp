#pragma once
#include "Operator.hpp"
#include "IndexWrapper.hpp"

namespace SymbolicOperators {
	struct Coefficient {
		std::string name;
		Momentum momentum;
		// Contains all indizes, standard: first index = spin, all others arbitrary, e.g. orbitals, bands etc
		IndexWrapper indizes;
		// if Coeff(k) = Coeff(-k)
		bool translationalInvariance = true;
		// if Coeff(k+Q) = -Coeff(k)
		bool Q_changes_sign{};
		bool isDaggered{};

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& name;
			ar& momentum;
			ar& indizes;
			ar& isDaggered;
			ar& translationalInvariance;
			ar& Q_changes_sign;
		}

		Coefficient();
		explicit Coefficient(std::string _name);
		Coefficient(std::string _name, const Momentum& _momentum, const IndexWrapper& _indizes, bool _Q_changes_sign = false, bool _isDaggered = false);
		Coefficient(std::string _name, char _momentum, bool add_Q, const IndexWrapper& _indizes, bool _Q_changes_sign = false, bool _isDaggered = false);
		Coefficient(std::string _name, const Momentum& _momentum, bool _Q_changes_sign = false, bool _isDaggered = false);
		Coefficient(std::string _name, char _momentum, bool add_Q = false, bool _Q_changes_sign = false, bool _isDaggered = false);

		inline bool usesIndex(const Index index) const noexcept {
			for (const auto& idx : indizes) {
				if (idx == index) return true;
			}
			return false;
		}
		inline bool dependsOnMomentum() const noexcept {
			return this->momentum.momentum_list.size() > 0;
		};
		inline bool dependsOn(char momentum) const noexcept {
			return this->momentum.isUsed(momentum) != -1;
		}
		inline bool dependsOnTwoMomenta() const noexcept {
			return this->momentum.momentum_list.size() == 2U;
		};
	};

	inline bool operator==(const Coefficient& lhs, const Coefficient& rhs) {
		if (lhs.name != rhs.name) return false;
		if (lhs.momentum != rhs.momentum) return false;
		if (lhs.isDaggered != rhs.isDaggered) return false;
		return (lhs.indizes == rhs.indizes);
	}
	inline bool operator!=(const Coefficient& lhs, const Coefficient& rhs) {
		return !(lhs == rhs);
	}
}