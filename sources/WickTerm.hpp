#pragma once
#include "Term.hpp"
#include <algorithm>

namespace SymbolicOperators {
	enum OperatorType { Number_Type, CDW_Type, SC_Type, Eta_Type,Undefined_Type};

	inline std::ostream& operator<<(std::ostream& os, const OperatorType& op) {
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

	class Term;

	struct WickTerm
	{
		int multiplicity;
		std::vector<Coefficient> coefficients;
		std::vector<char> sum_momenta;
		std::vector<std::string> sum_indizes;
		std::vector<WickOperator> operators;

		// symbolises the Kronecker delta
		std::vector<pair_of_momenta> delta_momenta;
		std::vector<std::pair<std::string, std::string>> delta_indizes;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& multiplicity;
			ar& coefficients;
			ar& sum_momenta;
			ar& sum_indizes;
			ar& operators;
			ar& delta_momenta;
			ar& delta_indizes;
		}

		std::vector<Operator> temporary_operators;

		explicit WickTerm(const Term* base);
		explicit WickTerm(const Term& base);
		WickTerm();

		inline bool includesType(const OperatorType operator_type) const {
			return std::any_of(this->operators.begin(), this->operators.end(),
				[operator_type](const WickOperator& op) { return op.type == operator_type; });
		};
		inline bool hasSingleCoefficient() const noexcept {
			return this->coefficients.size() == 1U;
		};
		inline bool usesIndex(const std::string& index) const noexcept {
			for (const auto& op : operators) {
				if (op.usesIndex(index)) return true;
			}
			for (const auto& coeff : coefficients) {
				if (coeff.usesIndex(index)) return true;
			}
			return false;
		}
		inline bool isIdentity() const noexcept {
			return this->operators.empty();
		}
		inline bool isBilinear() const noexcept {
			return this->operators.size() == 1U;
		};
		inline bool isQuartic() const noexcept {
			return this->operators.size() == 2U;
		}
		// Returns this->multiplicity, but always as a double
		inline double getFactor() const noexcept {
			return static_cast<double>(this->multiplicity);
		}
		// Returns the position of the first operator that depends on 'momentum'
		// Returns -1 if no operators depend on 'momentum'.
		inline int whichOperatorDependsOn(char momentum) const noexcept {
			for (int i = 0U; i < operators.size(); ++i)
			{
				if (operators[i].dependsOn(momentum)) return i;
			}
			return -1;
		}
		inline const Coefficient& getFirstCoefficient() const {
			assert(!(this->coefficients.empty()));
			return this->coefficients.front();
		};

		inline bool handled() const noexcept {
			if (this->temporary_operators.empty()) return true;
			return !(this->operators.empty());
		}
		bool swapToWickOperators(std::vector<WickTerm>& reciever);
		// returns false if there is atleast one delta
		// or a combination of deltas, that can never be achieved
		// for example delta_k,k+Q, as k can never be equal to k+Q
		bool setDeltas();
		// May call setDeltas. If setDeltas returns false this functions also returns false
		// In all other cases it returns true
		bool computeSums();
		void discardZeroMomenta();
		void renameSums();
		void sort();
		void applySpinSymmetry();
		void applyTranslationalSymmetry();
		void applyPhaseSymmetry();
	};

	inline bool operator==(const WickOperator& lhs, const WickOperator& rhs) {
		if (lhs.type != rhs.type) return false;
		if (lhs.isDaggered != rhs.isDaggered) return false;
		if (lhs.momentum != rhs.momentum) return false;
		return (lhs.indizes == rhs.indizes);
	};
	inline bool operator!=(const WickOperator& lhs, const WickOperator& rhs) {
		return !(lhs == rhs);
	};
	inline bool operator==(const WickTerm& lhs, const WickTerm& rhs) {
		if (lhs.coefficients != rhs.coefficients) return false;
		if (lhs.sum_indizes != rhs.sum_indizes) return false;
		if (lhs.sum_momenta != rhs.sum_momenta) return false;

		if (lhs.delta_indizes != rhs.delta_indizes) return false;
		if (lhs.delta_momenta != rhs.delta_momenta) return false;

		if (lhs.operators.size() != rhs.operators.size()) return false;
		// The Wick "operators" are actually just numbers
		// therefore I might be interested to implement permutations as well...
		for (size_t i = 0; i < lhs.operators.size(); i++)
		{
			if (lhs.operators[i] != rhs.operators[i]) return false;
		}
		return true;
	};
	inline bool operator!=(const WickTerm& lhs, const WickTerm& rhs) {
		return !(lhs == rhs);
	};

	void wicks_theorem(const Term& term, std::vector<WickTerm>& reciever);
	void clearEtas(std::vector<WickTerm>& terms);
	void cleanWicks(std::vector<WickTerm>& terms);

	std::ostream& operator<<(std::ostream& os, const WickOperator& op);
	std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops);
	std::ostream& operator<<(std::ostream& os, const WickTerm& term);
	std::ostream& operator<<(std::ostream& os, const std::vector<WickTerm>& terms);
}