#pragma once
#include <algorithm>

#include "Operator.hpp"
#include "Coefficient.hpp"
#include "WickTerm.hpp"

namespace SymbolicOperators {
	typedef std::pair<Momentum, Momentum> pair_of_momenta;
	inline bool pair_equal_allow_permutation(const std::pair<std::string, std::string>& left, const std::pair<std::string, std::string>& right) {
		if (left.first == right.first && left.second == right.second) return true;
		if (left.first == right.second && left.second == right.first) return true;
		return false;
	};
	inline bool pair_equal_allow_permutation(std::pair<Momentum, Momentum> left, std::pair<Momentum, Momentum> right) {
		if (left.first.add_Q) {
			left.first.add_Q = false;
			left.second.add_Q = !(left.second.add_Q);
		}
		if (right.first.add_Q) {
			right.first.add_Q = false;
			right.second.add_Q = !(right.second.add_Q);
		}
		if (left.first == right.first && left.second == right.second) return true;

		if (right.second.add_Q) {
			right.second.add_Q = false;
			right.first.add_Q = !(right.first.add_Q);
		}
		if (left.first == right.second && left.second == right.first) return true;
		return false;
	};

	class Term {
	private:
		std::vector<Coefficient> coefficients;
		std::vector<char> sum_momenta;
		std::vector<std::string> sum_indizes;
		std::vector<Operator> operators;
		// symbolises the Kronecker delta
		std::vector<pair_of_momenta> delta_momenta;
		std::vector<std::pair<std::string, std::string>> delta_indizes;

		friend struct WickTerm;
	public:
		int multiplicity;
		Term(int _multiplicity, Coefficient _coefficient, const std::vector<char>& _sum_momenta,
			const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(int _multiplicity, Coefficient _coefficient, const std::vector<char>& _sum_momenta,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(int _multiplicity, Coefficient _coefficient, const std::vector<std::string>& _sum_indizes,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(int _multiplicity, Coefficient _coefficient,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(int _multiplicity, const std::vector<char>& _sum_momenta, const std::vector<std::string>& _sum_indizes,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(int _multiplicity, const std::vector<char>& _sum_momenta,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(int _multiplicity, const std::vector<std::string>& _sum_indizes,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		explicit Term(int _multiplicity,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term();

		inline bool isIdentity() const {
			return this->operators.empty();
		}
		void print() const;
		inline void flipSign() {
			this->multiplicity *= -1;
		}
		inline const std::vector<Operator>& getOperators() const {
			return this->operators;
		}

		bool setDeltas();
		bool computeSums();
		void discardZeroMomenta();
		void sort();
		// Unifies the sum indizes
		void renameSums();
		//void wick(std::vector<WickTerm>& reciever) const;

		// Checks for equality of everything except of multiplicity
		inline bool isEqual(const Term& other) const {
			if (this->coefficients.size() != other.coefficients.size()) return false;
			for (size_t i = 0; i < coefficients.size(); i++)
			{
				if (this->coefficients[i] != other.coefficients[i]) return false;
			}
			if (this->sum_indizes.size() != other.sum_indizes.size()) return false;
			if (this->sum_momenta.size() != other.sum_momenta.size()) return false;
			for (size_t i = 0; i < sum_indizes.size(); i++)
			{
				if (this->sum_indizes[i] != other.sum_indizes[i]) return false;
			}
			for (size_t i = 0; i < sum_momenta.size(); i++)
			{
				if (this->sum_momenta[i] != other.sum_momenta[i]) return false;
			}

			if (this->delta_indizes.size() != other.delta_indizes.size()) return false;
			if (this->delta_momenta.size() != other.delta_momenta.size()) return false;
			for (size_t i = 0; i < delta_indizes.size(); i++)
			{
				if (this->delta_indizes[i] != other.delta_indizes[i]) return false;
			}
			for (size_t i = 0; i < delta_momenta.size(); i++)
			{
				if (this->delta_momenta[i] != other.delta_momenta[i]) return false;
			}

			if (this->operators.size() != other.operators.size()) return false;
			for (size_t i = 0; i < this->operators.size(); i++)
			{
				if (this->operators[i] != other.operators[i]) return false;
			}
			return true;
		};

		std::string toStringWithoutPrefactor() const;

		friend void normalOrder(std::vector<Term>& terms);
		friend void commutator(std::vector<Term>& reciever, const Term& left, const Term& right);
		friend std::ostream& operator<<(std::ostream& os, const Term& term);

		inline void hermitianConjugate() {
			std::reverse(this->operators.begin(), this->operators.end());
			for (auto& op : this->operators) {
				op.hermitianConjugate();
			}
		};
		inline void renameMomenta(char what, char to) {
			for (auto& op : operators)
			{
				op.momentum.replaceOccurances(what, Momentum(to));
			}
		};
	};
	void commutator(std::vector<Term>& reciever, const std::vector<Term>& left, const std::vector<Term>& right);
	inline void commutator(std::vector<Term>& reciever, const Term& left, const std::vector<Term>& right) {
		const std::vector<Term> buffer = { left };
		commutator(reciever, buffer, right);
	};
	inline void commutator(std::vector<Term>& reciever, const std::vector<Term>& left, const Term& right) {
		const std::vector<Term> buffer = { right };
		commutator(reciever, left, buffer);
	};

	inline bool operator==(const Term& lhs, const Term& rhs) {
		return lhs.isEqual(rhs);
	};
	inline bool operator!=(const Term& lhs, const Term& rhs) {
		return !(lhs.isEqual(rhs));
	};

	std::ostream& operator<<(std::ostream& os, const Coefficient& coeff);
	std::ostream& operator<<(std::ostream& os, const std::vector<Coefficient>& coeffs);
	std::ostream& operator<<(std::ostream& os, const std::vector<Term>& terms);

	void cleanUp(std::vector<Term>& terms);
	inline void hermitianConjugate(std::vector<Term>& terms) {
		for (auto& t : terms) {
			t.hermitianConjugate();
		}
	};
	inline void renameMomenta(std::vector<Term>& terms, char what, char to) {
		for (auto& t : terms) {
			t.renameMomenta(what, to);
		}
	};
	inline std::string toStringWithoutPrefactor(const std::vector<Term>& terms) {
		std::string ret = "";
		for (size_t i = 0; i < terms.size(); i++)
		{
			if (terms[i].multiplicity < 0) {
				ret += "-";
			}
			else if (i > 0) {
				ret += "+";
			}
			ret += terms[i].toStringWithoutPrefactor();
		}
		return ret;
	}
}