#define append_vector(a, b) a.insert(a.end(), b.begin(), b.end())

#include "Term.hpp"
#include <cmath>
#include <sstream>

namespace SymbolicOperators {
	Term::Term(int _multiplicity, Coefficient _coefficient, const std::vector<char>& _sum_momenta, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators)
		: coefficients(1, _coefficient), sum_momenta(_sum_momenta), sum_indizes(_sum_indizes), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(int _multiplicity, Coefficient _coefficient, const std::vector<char>& _sum_momenta, const std::vector<Operator>& _operators)
		: coefficients(1, _coefficient), sum_momenta(_sum_momenta), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(int _multiplicity, Coefficient _coefficient, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators)
		: coefficients(1, _coefficient), sum_indizes(_sum_indizes), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(int _multiplicity, Coefficient _coefficient, const std::vector<Operator>& _operators)
		: coefficients(1, _coefficient), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(int _multiplicity, const std::vector<char>& _sum_momenta, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators)
		: coefficients(), sum_momenta(_sum_momenta), sum_indizes(_sum_indizes), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(int _multiplicity, const std::vector<char>& _sum_momenta, const std::vector<Operator>& _operators)
		: coefficients(), sum_momenta(_sum_momenta), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(int _multiplicity, const std::vector<std::string>& _sum_indizes, const std::vector<Operator>& _operators)
		: coefficients(), sum_indizes(_sum_indizes), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(int _multiplicity, const std::vector<Operator>& _operators)
		: coefficients(), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term()
		: coefficients(), operators(), multiplicity(0) {}

	void Term::print() const {
		std::cout << *this << std::endl;
	}

	void Term::setDeltas()
	{
		// Set all deltas up to the same notation
		for (auto& delta : delta_momenta) {
			if (delta.first.add_Q) {
				delta.first.add_Q = false;
				delta.second.add_Q = (delta.second.add_Q != true);
			}
			if (delta.second.momentum_list.size() == 1 && delta.first.momentum_list.size() > 1) {
				std::swap(delta.first, delta.second);
			}
			if (delta.first.momentum_list.size() == 1 && delta.first.momentum_list[0].first < 0) {
				delta.first.flipMomentum();
				delta.second.flipMomentum();
			}
			else if (delta.first.momentum_list.size() > 1 && delta.second.momentum_list.size() > 1) {
				bool foundCandidate = false;
				int index = 0;
				delta.second -= delta.first;
				delta.first.momentum_list.clear();

				for (auto m : sum_momenta)
				{
					index = delta.second.isUsed(m);
					if (index >= 0) {
						foundCandidate = true;
						if (abs(delta.second.momentum_list[index].first) == 1) {
							break;
						}
					}
				}
				if (!foundCandidate) index = 0;

				if (delta.second.momentum_list[index].first > 0) {
					delta.second.flipMomentum();
				}
				delta.first.momentum_list.push_back(delta.second.momentum_list[index]);
				delta.first.flipMomentum();
				if (abs(delta.first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta.first << std::endl;
				delta.second.momentum_list.erase(delta.second.momentum_list.begin() + index);
			}

			if (abs(delta.first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta.first << std::endl;
			for (auto& op : operators) {
				op.momentum.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
			}
			for (auto& coeff : coefficients) {
				coeff.momentum.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
			}
		}
		for (auto& delta : delta_indizes) {
			for (auto& op : operators) {
				for (auto it = op.indizes.begin(); it != op.indizes.end(); ++it)
				{
					if (delta.first == UP || delta.first == DOWN) {
						if (*it == delta.second) {
							*it = delta.first;
						}
					}
					else {
						if (*it == delta.first) {
							*it = delta.second;
						}
					}
				}
			}
		}

		// Remove delta^2
		for (int i = 0; i < delta_momenta.size(); i++)
		{
			for (int j = i + 1; j < delta_momenta.size(); j++)
			{
				if (pair_equal_allow_permutation(delta_momenta[i], delta_momenta[j])) {
					delta_momenta.erase(delta_momenta.begin() + j);
					--i;
					break;
				}
			}
		}
		for (int i = 0; i < delta_indizes.size(); i++)
		{
			for (int j = i + 1; j < delta_indizes.size(); j++)
			{
				if (pair_equal_allow_permutation(delta_indizes[i], delta_indizes[j])) {
					delta_indizes.erase(delta_indizes.begin() + j);
					--i;
					break;
				}
			}
		}
	}

	void Term::computeSums() {
		auto changeAllIndizes = [&](const std::string& replaceWhat, const std::string& replaceWith) {
			for (auto& op : operators) {
				for (auto it = op.indizes.begin(); it != op.indizes.end(); ++it)
				{
					if (*it == replaceWhat) {
						*it = replaceWith;
					}
				}
			}
			for (auto& coeff : coefficients) {
				for (auto it = coeff.indizes.begin(); it != coeff.indizes.end(); ++it)
				{
					if (*it == replaceWhat) {
						*it = replaceWith;
					}
				}
			}
		};

		for (int i = 0; i < sum_indizes.size(); i++)
		{
			for (int j = 0; j < delta_indizes.size(); j++)
			{
				if (delta_indizes[j].first == sum_indizes[i]) {
					changeAllIndizes(sum_indizes[i], delta_indizes[j].second);
					sum_indizes.erase(sum_indizes.begin() + i);
					delta_indizes.erase(delta_indizes.begin() + j);
					--i;
					break;
				}
				else if (delta_indizes[j].second == sum_indizes[i]) {
					changeAllIndizes(sum_indizes[i], delta_indizes[j].first);
					sum_indizes.erase(sum_indizes.begin() + i);
					delta_indizes.erase(delta_indizes.begin() + j);
					--i;
					break;
				}
			}
		}

		auto changeAllMomenta = [&](const char replaceWhat, const Momentum& replaceWith) {
			for (auto& op : operators) {
				op.momentum.replaceOccurances(replaceWhat, replaceWith);
			}
			for (auto& coeff : coefficients) {
				coeff.momentum.replaceOccurances(replaceWhat, replaceWith);
			}
		};

		for (int i = 0; i < sum_momenta.size(); i++)
		{
			for (int j = 0; j < delta_momenta.size(); j++)
			{
				if (delta_momenta[j].first.momentum_list[0].second == sum_momenta[i]) {
					changeAllMomenta(sum_momenta[i], delta_momenta[j].second);
					if (abs(delta_momenta[j].first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta_momenta[j].first << std::endl;

					sum_momenta.erase(sum_momenta.begin() + i);
					delta_momenta.erase(delta_momenta.begin() + j);
					--i;
					break;
				}
				else {
					int index = delta_momenta[j].second.isUsed(sum_momenta[i]);
					if (index < 0) continue;

					Momentum buffer(delta_momenta[j].second.momentum_list[index].second, delta_momenta[j].second.momentum_list[index].first);
					if (abs(buffer.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << buffer << std::endl;
					delta_momenta[j].second.momentum_list.erase(delta_momenta[j].second.momentum_list.begin() + index);
					delta_momenta[j].second -= delta_momenta[j].first;

					if (buffer.momentum_list[0].first > 0) {
						delta_momenta[j].second.flipMomentum();
					}
					changeAllMomenta(sum_momenta[i], delta_momenta[j].second);

					sum_momenta.erase(sum_momenta.begin() + i);
					delta_momenta.erase(delta_momenta.begin() + j);
					--i;
					break;
				}
			}
		}
	}

	void Term::discardZeroMomenta() {
		for (auto& op : operators) {
			for (auto it = op.momentum.momentum_list.begin(); it != op.momentum.momentum_list.end();) {
				if (it->first == 0) {
					it = op.momentum.momentum_list.erase(it);
				}
				else {
					++it;
				}
			}
		}
		for (auto& coeff : coefficients) {
			for (auto it = coeff.momentum.momentum_list.begin(); it != coeff.momentum.momentum_list.end();) {
				if (it->first == 0) {
					it = coeff.momentum.momentum_list.erase(it);
				}
				else {
					++it;
				}
			}
		}
	}

	void Term::sort()
	{
		for (auto& coeff : coefficients) {
			if (coeff.translationalInvariance && coeff.momentum.momentum_list.size() > 0) {
				if (coeff.momentum.momentum_list[0].first < 0) {
					coeff.momentum.flipMomentum();
				}
			}
		}
		for (int i = 0; i < operators.size(); i++)
		{
			for (int j = i + 1; j < operators.size(); j++)
			{
				if (operators[i].isDaggered == operators[j].isDaggered) {
					if (operators[i].isDaggered) {
						// c^+ c^+
						if (operators[j].indizes[0] == UP && operators[i].indizes[0] != UP) {
							std::swap(operators[i], operators[j]);
							if (abs(i - j) % 2 != 0) flipSign();
						}
						else if (operators[i].indizes[0] == DOWN && operators[j].indizes[0] != DOWN) {
							std::swap(operators[i], operators[j]);
							if (abs(i - j) % 2 != 0) flipSign();
						}
					}
					else {
						// c c
						if (operators[j].indizes[0] == DOWN && operators[i].indizes[0] != DOWN) {
							std::swap(operators[i], operators[j]);
							if (abs(i - j) % 2 != 0) flipSign();
						}
						else if (operators[i].indizes[0] == UP && operators[j].indizes[0] != UP) {
							std::swap(operators[i], operators[j]);
							if (abs(i - j) % 2 != 0) flipSign();
						}
					}
				}
			}
		}
		for (int i = 0; i < operators.size(); i++)
		{
			for (int j = i + 1; j < operators.size(); j++)
			{
				if (operators[i].isDaggered != operators[j].isDaggered) continue;
				if (operators[i].indizes[0] != operators[j].indizes[0]) continue;

				if (operators[i].momentum.momentum_list[0].second > operators[j].momentum.momentum_list[0].second) {
					std::swap(operators[i], operators[j]);
					if (abs(i - j) % 2 != 0) flipSign();
				}
			}
		}
	}

	void Term::renameSums()
	{
		const char name_list[3] = { 'q', 'p', 'r' };
		const char buffer_list[3] = { ':', ';', '|' };
		for (int i = 0; i < sum_momenta.size(); i++)
		{
			if (i >= 3) {
				std::cerr << "More than 3 momenta, time to implement this..." << std::endl;
				break;
			}
			if (sum_momenta[i] == name_list[i]) continue;

			for (auto& op : operators) {
				op.momentum.replaceOccurances(sum_momenta[i], Momentum(buffer_list[i]));
			}
			for (auto& coeff : coefficients) {
				coeff.momentum.replaceOccurances(sum_momenta[i], Momentum(buffer_list[i]));
			}
			sum_momenta[i] = name_list[i];
		}

		for (int i = 0; i < sum_momenta.size(); i++)
		{
			for (auto& op : operators) {
				op.momentum.replaceOccurances(buffer_list[i], Momentum(name_list[i]));
			}
			for (auto& coeff : coefficients) {
				coeff.momentum.replaceOccurances(buffer_list[i], Momentum(name_list[i]));
			}
		}
	}

	std::string Term::toStringWithoutPrefactor() const
	{
		std::ostringstream os;
		if (!this->sum_indizes.empty()) {
			os << "\\sum_{ ";
			for (const auto& index : this->sum_indizes) {
				os << index << " ";
			}
			os << "}";
		}
		if (!this->sum_momenta.empty()) {
			os << "\\sum_{ ";
			for (const auto& momentum : this->sum_momenta) {
				os << momentum << " ";
			}
			os << "}";
		}
		os << this->coefficients << " ";
		for (const auto& delta : this->delta_momenta) {
			os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
		}
		for (const auto& delta : this->delta_indizes) {
			os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
		}

		if (this->isIdentity()) {
			os << " \\mathbb{1} ";
			return os.str();
		}
		for (const auto& op : this->operators) {
			os << op << " ";
		}
		return os.str();
	}

	void normalOrder(std::vector<Term>& terms) {
		for (int t = 0; t < terms.size();) {
		normalOder_outerLoop:
			if (t >= terms.size()) break;
			int n = terms[t].operators.size();
			int new_n;
			while (n > 1) {
				new_n = 0;
				for (int i = 1; i < terms[t].operators.size(); i++)
				{
					if (!(terms[t].operators[i - 1].isDaggered) && (terms[t].operators[i].isDaggered)) {
						bool other_deltas = false;
						new_n = i;
						// Swap cc^+
						terms[t].flipSign();
						std::swap(terms[t].operators[i - 1], terms[t].operators[i]);

						// Add a new term where cc^+ is replaced by the appropriate delta
						Term new_term(terms[t]);
						new_term.flipSign();
						if (new_term.operators[i - 1].indizes.size() != new_term.operators[i].indizes.size()) {
							std::cerr << "Operators do not have the same index count." << std::endl;
							throw;
						}

						if ((new_term.operators[i - 1].indizes[0] == UP || new_term.operators[i - 1].indizes[0] == DOWN)
							&& (new_term.operators[i].indizes[0] == UP || new_term.operators[i].indizes[0] == DOWN)) {
							if (new_term.operators[i - 1].indizes[0] != new_term.operators[i].indizes[0]) {
								// Case: one up the other down, then free anticommutation
								continue;
							}
						}
						else {
							new_term.delta_indizes.push_back(
								std::make_pair(new_term.operators[i - 1].indizes[0], new_term.operators[i].indizes[0]));
						}
						for (int c = 1; c < new_term.operators[i - 1].indizes.size(); c++)
						{
							// if the indizes are not the same we emplace a delta
							// otherwise no action is required
							if (new_term.operators[i - 1].indizes[c] != new_term.operators[i].indizes[c]) {
								other_deltas = true;
								new_term.delta_indizes.push_back(
									std::make_pair(new_term.operators[i - 1].indizes[c], new_term.operators[i].indizes[c]));
							}
						}
						if (new_term.operators[i - 1].momentum != new_term.operators[i].momentum) {
							other_deltas = true;
							new_term.delta_momenta.push_back(
								std::make_pair(new_term.operators[i - 1].momentum, new_term.operators[i].momentum)
							);
						}

						new_term.operators.erase(new_term.operators.begin() + i - 1, new_term.operators.begin() + i + 1);
						if (other_deltas) terms.push_back(new_term);
					}
					else if (terms[t].operators[i - 1] == terms[t].operators[i]) {
						// two identical fermion operators = 0
						terms.erase(terms.begin() + t);
						goto normalOder_outerLoop;
					}
				}
				n = new_n;
			}
			++t;
		}
	}

#define fill_reciever(x) reciever[0].x = left.x; append_vector(reciever[0].x, right.x); reciever[1].x = left.x; append_vector(reciever[1].x, right.x);
	void commutator(std::vector<Term>& reciever, const Term& left, const Term& right)
	{
		reciever.resize(2);
		reciever[0] = left;
		reciever[0].multiplicity *= right.multiplicity;
		append_vector(reciever[0].operators, right.operators);
		reciever[1] = right;
		reciever[1].multiplicity *= left.multiplicity;
		append_vector(reciever[1].operators, left.operators);
		reciever[1].flipSign();

		fill_reciever(coefficients);
		fill_reciever(sum_momenta);
		fill_reciever(sum_indizes);
		fill_reciever(delta_momenta);
		fill_reciever(delta_indizes);

		normalOrder(reciever);
	}
	void commutator(std::vector<Term>& reciever, const std::vector<Term>& left, const std::vector<Term>& right)
	{
		reciever.reserve(2 * left.size() * right.size());
		std::vector<Term> reciever_buffer(2);

		for (int i = 0; i < left.size(); i++)
		{
			for (int j = 0; j < right.size(); j++)
			{
				commutator(reciever_buffer, left[i], right[j]);
				append_vector(reciever, reciever_buffer);
			}
		}
	}

	std::ostream& operator<<(std::ostream& os, const Term& term)
	{
		if (term.multiplicity > 0) {
			os << "+";
		}
		os << term.multiplicity << " \\cdot ";
		if (!term.sum_indizes.empty()) {
			os << "\\sum_{ ";
			for (const auto& index : term.sum_indizes) {
				os << index << " ";
			}
			os << "}";
		}
		if (!term.sum_momenta.empty()) {
			os << "\\sum_{ ";
			for (const auto& momentum : term.sum_momenta) {
				os << momentum << " ";
			}
			os << "}";
		}
		os << term.coefficients << " ";
		for (const auto& delta : term.delta_momenta) {
			os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
		}
		for (const auto& delta : term.delta_indizes) {
			os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
		}

		if (term.isIdentity()) {
			os << " \\mathbb{1} ";
			return os;
		}
		for (const auto& op : term.operators) {
			os << op << " ";
		}
		return os;
	}
	std::ostream& operator<<(std::ostream& os, const std::vector<Term>& terms)
	{
		for (std::vector<Term>::const_iterator it = terms.begin(); it != terms.end(); ++it)
		{
			os << "\t&" << *it;
			if (it != terms.end() - 1) {
				os << " \\\\";
			}
			os << "\n";
		}
		return os;
	}

	void cleanUp(std::vector<Term>& terms)
	{
		for (auto& term : terms) {
			term.setDeltas();
			term.computeSums();
			term.discardZeroMomenta();
			term.setDeltas();
			term.discardZeroMomenta();
			term.renameSums();
			term.sort();
		}

		// remove duplicates
		for (int i = 0; i < terms.size(); i++)
		{
			for (int j = i + 1; j < terms.size(); j++)
			{
				if (terms[i] == terms[j]) {
					terms[i].multiplicity += terms[j].multiplicity;
					terms.erase(terms.begin() + j);
					--i;
					break;
				}
			}
		}

		// removes any terms that have a 0 prefactor
		for (auto it = terms.begin(); it != terms.end();)
		{
			if (it->multiplicity == 0) {
				it = terms.erase(it);
			}
			else {
				++it;
			}
		}
	}

	void Term::wick(std::vector<WickTerm>& reciever) const
	{
		if (this->operators.size() > 4) std::cerr << "Wick for n>4 not yet implemented!" << std::endl;
		if (this->isIdentity()) {
			reciever.push_back(WickTerm(this));
		}
		else {
			for (int i = 1; i < operators.size(); i++)
			{
				WickTerm buffer(this);
				if ((i % 2) == 0) {
					buffer.multiplicity *= -1;
				}
				buffer.temporary_operators.reserve(buffer.temporary_operators.size() + 2);
				buffer.temporary_operators.push_back(operators[0]);
				buffer.temporary_operators.push_back(operators[i]);

				if (this->operators.size() > 2) {
					std::vector<Operator> copy_these_operators = this->operators;
					copy_these_operators.erase(copy_these_operators.begin() + i);
					copy_these_operators.erase(copy_these_operators.begin());

					buffer.temporary_operators.reserve(buffer.temporary_operators.size() + 2);
					buffer.temporary_operators.push_back(copy_these_operators[0]);
					buffer.temporary_operators.push_back(copy_these_operators[1]);
				}

				reciever.push_back(buffer);
			}
		}

		for (int i = 0; i < reciever.size();)
		{
			if (reciever[i].handled()) {
				++i;
				continue;
			}
			if (!(reciever[i].swapToWickOperators(reciever))) {
				reciever.erase(reciever.begin() + i);
			}
			else {
				++i;
			}
		}
	}
}