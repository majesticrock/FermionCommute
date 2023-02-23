#pragma once
#include "Term.hpp"

namespace SymbolicOperators {
	typedef std::pair<Momentum, Momentum> pair_of_momenta;
	struct WickOperator {
		std::string type;
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

		WickOperator(const std::string& _type, const bool _isDaggered, const Momentum& _momentum, const std::vector<std::string>& _indizes = std::vector<std::string>());
		WickOperator(const std::string& _type, const bool _isDaggered, const Momentum& _momentum, const std::string& _index);
		WickOperator();
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
		WickTerm();

		inline bool isIdentity() const {
			return this->operators.empty();
		}
		inline bool handled() const {
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

	void cleanWicks(std::vector<WickTerm>& terms);

	std::ostream& operator<<(std::ostream& os, const WickOperator& op);
	std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops);
	std::ostream& operator<<(std::ostream& os, const WickTerm& term);
	std::ostream& operator<<(std::ostream& os, const std::vector<WickTerm>& terms);
}