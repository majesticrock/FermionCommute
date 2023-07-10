#define append_vector(a, b) a.insert(a.end(), b.begin(), b.end())

#include "Term.hpp"
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

	bool Term::setDeltas()
	{
		for (auto& delta : delta_momenta)
		{
			for (auto it = delta.first.momentum_list.begin(); it != delta.first.momentum_list.end(); )
			{
				int index = delta.second.isUsed(it->second);
				if (index < 0) {
					++it;
					continue;
				}

				int remainder = delta.second.momentum_list[index].first - it->first;
				if (remainder == 0) {
					delta.second.momentum_list.erase(delta.second.momentum_list.begin() + index);
					it = delta.first.momentum_list.erase(it);
					continue;
				}

				delta.second.momentum_list[index].first = remainder;
				it = delta.first.momentum_list.erase(it);
			}
			if (delta.first.momentum_list.size() == 0) {
				if (delta.second.momentum_list.size() == 0) continue;
				std::swap(delta.first, delta.second);
			}
			if (delta.first.add_Q) {
				delta.first.add_Q = false;
				delta.second.add_Q = !(delta.second.add_Q);
			}
			if (delta.first.momentum_list.front().first < 0) {
				delta.first.flipMomentum();
				delta.second.flipMomentum();
			}
			if (delta.first.momentum_list.size() > 1 && delta.second.momentum_list.size() == 0) {
				delta.second.momentum_list.push_back(delta.first.momentum_list[1]);
				delta.second.flipMomentum();
				delta.first.momentum_list.erase(delta.first.momentum_list.begin() + 1);
			}
		}

		for (auto it = delta_momenta.begin(); it != delta_momenta.end(); )
		{
			if (it->first.momentum_list.size() == 0 && it->second.momentum_list.size() == 0) {
				// 0 = Q can never be achieved
				if (it->first.add_Q != it->second.add_Q) return false;
				it = delta_momenta.erase(it);
			}
			else {
				++it;
			}
		}

		// Set all deltas up to the same notation
		for (auto& delta : delta_momenta) {
			for (auto& delta2 : delta_momenta) {
				for (auto it = delta2.first.momentum_list.begin(); it != delta2.first.momentum_list.end();) {
					int pos = delta2.second.isUsed(it->second);
					if (pos < 0) { ++it; continue; }
					it->first -= delta2.second.momentum_list[pos].first;
					if (it->first == 0) {
						it = delta2.first.momentum_list.erase(it);
						delta2.second.momentum_list.erase(delta2.second.momentum_list.begin() + pos);
						continue;
					}
					++it;
				}
			}

			if (delta.first.momentum_list.size() == 0) {
				if (delta.second.momentum_list.size() == 0) continue;
				if (delta.second.momentum_list.size() == 1) {
					std::swap(delta.first, delta.second);
				}
				else {
					delta.first.momentum_list.push_back(delta.second.momentum_list.back());
					if (delta.first.momentum_list.front().first > 0) {
						delta.second.flipMomentum();
					}
					else {
						delta.first.flipMomentum();
					}
					delta.second.momentum_list.pop_back();
				}
			}
			if (delta.second.momentum_list.size() == 1 && delta.first.momentum_list.size() > 1) {
				std::swap(delta.first, delta.second);
			}
			if (delta.first.momentum_list.size() > 1 && delta.second.momentum_list.size() > 1) {
				bool foundCandidate = false;
				int index = 0;
				delta.second -= delta.first;
				delta.first.momentum_list.clear();

				for (auto m : sum_momenta)
				{
					index = delta.second.isUsed(m);
					if (index >= 0) {
						foundCandidate = true;
						if (std::abs(delta.second.momentum_list[index].first) == 1) {
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
				if (std::abs(delta.first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta.first << std::endl;
				delta.second.momentum_list.erase(delta.second.momentum_list.begin() + index);
			}
			if (delta.first.add_Q) {
				delta.first.add_Q = false;
				delta.second.add_Q = !(delta.second.add_Q);
			}
			if (delta.first.momentum_list.size() == 1 && delta.first.momentum_list[0].first < 0) {
				delta.first.flipMomentum();
				delta.second.flipMomentum();
			}

			if (std::abs(delta.first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta.first << std::endl;
			for (auto& op : operators) {
				op.momentum.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
			}
			for (auto& coeff : coefficients) {
				coeff.momentum.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
			}
			for (auto& delta2 : delta_momenta) {
				if (delta2 == delta) continue;
				delta2.first.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
				delta2.second.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
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

				auto delta_buffer = delta_momenta[j];
				delta_buffer.first.flipMomentum();
				delta_buffer.second.flipMomentum();
				if (pair_equal_allow_permutation(delta_momenta[i], delta_buffer)) {
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

		// Erase delta_k,k etc
		for (auto it = delta_momenta.begin(); it != delta_momenta.end();)
		{
			// k = k + Q can never be achieved
			if (it->first.differsOnlyInQ(it->second)) return false;

			if (it->first == it->second) {
				it = delta_momenta.erase(it);
			}
			else {
				++it;
			}
		}
		for (auto it = delta_indizes.begin(); it != delta_indizes.end();)
		{
			// UP can never be DOWN and vice versa
			if (it->first == UP && it->second == DOWN) return false;
			if (it->first == DOWN && it->second == UP) return false;

			if (it->first == it->second) {
				it = delta_indizes.erase(it);
			}
			else {
				++it;
			}
		}
		return true;
	}

	bool Term::computeSums() {
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
			for (auto it = delta_momenta.begin(); it != delta_momenta.end();) {
				it->first.replaceOccurances(replaceWhat, replaceWith);
				it->second.replaceOccurances(replaceWhat, replaceWith);
				++it;
			}
		};

		for (int i = 0; i < sum_momenta.size(); i++)
		{
			for (int j = 0; j < delta_momenta.size(); j++)
			{
				if (delta_momenta[j].first.momentum_list[0].second == sum_momenta[i]) {
					changeAllMomenta(sum_momenta[i], delta_momenta[j].second);
					if (!(delta_momenta[j].first.momentum_list.empty())) {
						if (std::abs(delta_momenta[j].first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta_momenta[j].first << std::endl;
					}

					sum_momenta.erase(sum_momenta.begin() + i);
					delta_momenta.erase(delta_momenta.begin() + j);
					--i;
					if (!(setDeltas())) return false;
					break;
				}
				else {
					int index = delta_momenta[j].second.isUsed(sum_momenta[i]);
					if (index < 0) continue;

					Momentum buffer(delta_momenta[j].second.momentum_list[index].second, delta_momenta[j].second.momentum_list[index].first);
					if (std::abs(buffer.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << buffer << std::endl;
					delta_momenta[j].second.momentum_list.erase(delta_momenta[j].second.momentum_list.begin() + index);
					delta_momenta[j].second -= delta_momenta[j].first;

					if (buffer.momentum_list[0].first > 0) {
						delta_momenta[j].second.flipMomentum();
					}
					changeAllMomenta(sum_momenta[i], delta_momenta[j].second);

					sum_momenta.erase(sum_momenta.begin() + i);
					delta_momenta.erase(delta_momenta.begin() + j);
					--i;
					if (!(setDeltas())) return false;
					break;
				}
			}
		}
		return true;
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
			if (coeff.Q_changes_sign && coeff.momentum.add_Q) {
				coeff.momentum.add_Q = false;
				flipSign();
			}
		}

		size_t new_n;
		size_t n = operators.size();
		while (n > 1U) {
			new_n = 0U;
			for (size_t i = 1U; i < n; ++i)
			{
				if (operators[i].isDaggered != operators[i - 1].isDaggered) continue;
				if (operators[i].isDaggered) {
					// c^+ c^+
					if (operators[i].indizes[0] == UP && operators[i - 1].indizes[0] != UP) {
						std::swap(operators[i], operators[i - 1]);
						flipSign();
						new_n = i;
					}
					else if (operators[i - 1].indizes[0] == DOWN && operators[i].indizes[0] != DOWN) {
						std::swap(operators[i], operators[i - 1]);
						flipSign();
						new_n = i;
					}
				}
				else {
					// c c
					if (operators[i].indizes[0] == DOWN && operators[i - 1].indizes[0] != DOWN) {
						std::swap(operators[i], operators[i - 1]);
						flipSign();
						new_n = i;
					}
					else if (operators[i - 1].indizes[0] == UP && operators[i].indizes[0] != UP) {
						std::swap(operators[i], operators[i - 1]);
						flipSign();
						new_n = i;
					}
				}
			}
			n = new_n;
		}

		n = operators.size();
		while (n > 1U) {
			new_n = 0U;
			for (size_t i = 1U; i < n; ++i)
			{
				if (operators[i].isDaggered != operators[i - 1].isDaggered) continue;
				if (operators[i].indizes[0] != operators[i - 1].indizes[0]) continue;

				if (operators[i - 1].momentum.momentum_list[0].second > operators[i].momentum.momentum_list[0].second) {
					std::swap(operators[i], operators[i - 1]);
					flipSign();
					new_n = i;
				}
			}
			n = new_n;
		}

		// check whether we can swap the sign of each momentum in the coefficients
		for (const auto& coeff : coefficients) {
			if (!(coeff.translationalInvariance)) return;
			if (coeff.momentum.momentum_list.size() > 1) return;
		}

		for (const auto& sum_mom : sum_momenta) {
			bool first_occurance = true;
			for (auto& op : operators) {
				int i = op.momentum.isUsed(sum_mom);
				if (i > -1) {
					if (first_occurance) {
						if (op.momentum.momentum_list[i].first < 0) {
							first_occurance = false;
						}
						else {
							break;
						}
					}
					op.momentum.momentum_list[i].first *= -1;
				}
			}
		}
	}

	void Term::renameSums()
	{
		constexpr char name_list[3] = { 'q', 'p', 'r' };
		constexpr char buffer_list[3] = { ':', ';', '|' };
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

		if (sum_indizes.size() == 1 && sum_indizes.front() == "\\sigma'") {
			sum_indizes.front() = "\\sigma";
			for (auto& op : operators) {
				for (auto& index : op.indizes) {
					if (index == "\\sigma'") index = "\\sigma";
				}
			}
			for (auto& coeff : coefficients) {
				for (auto& index : coeff.indizes) {
					if (index == "\\sigma'") index = "\\sigma";
				}
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
		for (const auto& delta : delta_momenta) {
			os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
		}
		for (const auto& delta : delta_indizes) {
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
			size_t n = terms[t].operators.size();
			size_t new_n;
			while (n > 1U) {
				new_n = 0U;
				for (size_t i = 1U; i < terms[t].operators.size(); ++i)
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
							throw std::invalid_argument("Operators do not have the same index count.");
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
						else {
							other_deltas = true;
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
		for (std::vector<Term>::iterator it = terms.begin(); it != terms.end();) {
			//std::cout << count++ << " of " << terms.size() << ":&\t" << *it << "\\\\" << std::endl;
			if (!(it->setDeltas())) {
				it = terms.erase(it);
				continue;
			}
			it->discardZeroMomenta();
			if (!(it->computeSums())) {
				it = terms.erase(it);
				continue;
			}
			if (!(it->setDeltas())) {
				it = terms.erase(it);
				continue;
			}
			it->discardZeroMomenta();
			it->renameSums();
			it->sort();

			++it;
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
}