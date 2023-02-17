#include "WickTerm.hpp"
#include <map>

namespace SymbolicOperators {
	bool WickTerm::setDeltas()
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

	void WickTerm::computeSums()
	{
		auto changeAllIndizes = [&](const std::string replaceWhat, const std::string replaceWith) {
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
			for (auto& delta : delta_indizes) {
				if (delta.first == replaceWhat) {
					delta.first = replaceWith;
				}
				if (delta.second == replaceWhat) {
					delta.second = replaceWith;
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

		auto changeAllMomenta = [&](const char replaceWhat, const Momentum replaceWith) {
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

	void WickTerm::discardZeroMomenta()
	{
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

	void WickTerm::renameSums()
	{
		const char name_list[3] = { 'q', 'p', 'r' };
		const char buffer_list[3] = { ':', ';', '|' };
		for (int i = 0; i < sum_momenta.size(); i++)
		{
			if (i >= 3) {
				throw std::invalid_argument("More than 3 momenta, time to implement this...");
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

		for (const auto& sum : sum_momenta)
		{
			for (auto& op : operators) {
				int index = op.momentum.isUsed(sum);
				if (index < 0) continue;
				if (op.momentum.momentum_list.size() == 1) break;

				Momentum buffer = op.momentum;
				buffer.flipMomentum();
				buffer.momentum_list[index].first *= -1;
				buffer.momentum_list[index].second = buffer_list[0];

				for (auto& op2 : operators) {
					op2.momentum.replaceOccurances(sum, buffer);
					op2.momentum.replaceOccurances(buffer_list[0], Momentum(sum));
				}
				for (auto& coeff : coefficients) {
					coeff.momentum.replaceOccurances(sum, buffer);
					coeff.momentum.replaceOccurances(buffer_list[0], Momentum(sum));
				}
			}
		}
		this->discardZeroMomenta();
	}

	void WickTerm::sort()
	{
		std::map<std::string, int> sort_map;
		sort_map["n"] = 0;
		sort_map["g"] = 1;
		sort_map["f"] = 2;
		sort_map["\\eta"] = 3;

		for (auto& op : operators) {
			if (op.type == "g" && op.momentum.add_Q) {
				op.momentum.add_Q = false;
				op.isDaggered = !(op.isDaggered);
			}
		}

		for (int i = 0; i < operators.size(); i++)
		{
			for (int j = i + 1; j < operators.size(); j++)
			{
				if (sort_map.find(operators[i].type)->second > sort_map.find(operators[j].type)->second) {
					std::swap(operators[i], operators[j]);
				}
				else if (operators[i].type == operators[j].type) {
					if (operators[i].momentum.momentum_list[0].second > operators[j].momentum.momentum_list[0].second) {
						std::swap(operators[i], operators[j]);
					}
				}
			}
		}

		for (auto& delta : delta_momenta) {
			if (delta.first.momentum_list.size() == 1 && delta.second.momentum_list.size() == 1) {
				// This comparison is well defined because we save the momentum as char i.e. byte
				// which is easily comparable
				if (delta.first.momentum_list[0].second > delta.second.momentum_list[0].second) {
					std::swap(delta.first, delta.second);
					if (delta.first.momentum_list[0].first < 0) {
						delta.first.flipMomentum();
						delta.second.flipMomentum();
					}
					if (delta.first.add_Q) {
						delta.first.add_Q = false;
						delta.second.add_Q = !(delta.second.add_Q);
					}
				}
			}
		}
	}

	void WickTerm::applySpinSymmetry()
	{
		// Expectation values for spin up and down are the same
		for (auto& op : operators) {
			if (op.indizes.size() > 0) {
				op.indizes[0] = UP;
			}
		}
	}

	void WickTerm::applyTranslationalSymmetry()
	{
		// Expectation values for k and -k are the same
		for (auto& op : operators) {
			if (op.momentum.momentum_list[0].first < 0) {
				op.momentum.flipMomentum();
			}
		}
	}

	void WickTerm::applyPhaseSymmetry()
	{
		// g^+ = g and f^+ = f
		for (auto& op : operators) {
			if (op.type == "f" || op.type == "g") {
				op.isDaggered = false;
			}
		}
	}

	void cleanWicks(std::vector<WickTerm>& terms)
	{
		for (std::vector<WickTerm>::iterator it = terms.begin(); it != terms.end();) {
			//std::cout << count++ << " of " << terms.size() << ":&\t" << *it << "\\\\" << std::endl;
			if (!(it->setDeltas())) {
				it = terms.erase(it);
				continue;
			}
			it->discardZeroMomenta();
			it->computeSums();
			if (!(it->setDeltas())) {
				it = terms.erase(it);
				continue;
			}
			it->discardZeroMomenta();
			it->renameSums();
			it->sort();

			it->applyPhaseSymmetry();
			it->applySpinSymmetry();
			it->applyTranslationalSymmetry();
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

		for (size_t i = 0; i < terms.size(); i++)
		{
			for (size_t j = i + 1; j < terms.size(); j++)
			{
				if (terms[i].delta_momenta.size() == 0 && terms[j].delta_momenta.size() > 0) {
					std::swap(terms[i], terms[j]);
				}
				if (terms[i].delta_momenta.size() > 0 && terms[j].delta_momenta.size() > 0) {
					if (terms[i].delta_momenta.size() > terms[j].delta_momenta.size()) {
						std::swap(terms[i], terms[j]);
					}
					else if (terms[i].delta_momenta.size() == terms[j].delta_momenta.size()) {
						if (terms[i].delta_momenta[0].second.add_Q && !(terms[j].delta_momenta[0].second.add_Q)) {
							std::swap(terms[i], terms[j]);
						}
					}
				}
			}
		}
	}
}