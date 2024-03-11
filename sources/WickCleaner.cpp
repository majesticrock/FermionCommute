#include "WickTerm.hpp"
#include <map>

namespace SymbolicOperators {
	bool WickTerm::setDeltas()
	{
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

			if (delta.first.momentum_list.size() == 1 && delta.first.momentum_list[0].first < 0) {
				delta.first.flipMomentum();
				delta.second.flipMomentum();
			}
			if (delta.first.add_Q) {
				delta.first.add_Q = false;
				delta.second.add_Q = !(delta.second.add_Q);
			}

			if (abs(delta.first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta.first << std::endl;
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
					if (delta.first == SpinUp || delta.first == SpinDown) {
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
				if (delta_momenta[i] == delta_momenta[j]) {
					delta_momenta.erase(delta_momenta.begin() + j);
					--i;
					break;
				}

				auto delta_buffer = delta_momenta[j];
				delta_buffer.first.flipMomentum();
				delta_buffer.second.flipMomentum();
				if (delta_momenta[i] == delta_buffer) {
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
				if (delta_indizes[i] == delta_indizes[j]) {
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
			// SpinUp can never be SpinDown and vice versa
			if (it->first == SpinUp && it->second == SpinDown) return false;
			if (it->first == SpinDown && it->second == SpinUp) return false;

			if (it->first == it->second) {
				it = delta_indizes.erase(it);
			}
			else {
				++it;
			}
		}
		return true;
	}

	bool WickTerm::computeSums()
	{
		auto changeAllIndizes = [&](const Index replaceWhat, const Index replaceWith) {
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
					if (abs(delta_momenta[j].first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta_momenta[j].first << std::endl;
					changeAllMomenta(sum_momenta[i], delta_momenta[j].second);

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
					if (!(setDeltas())) return false;
					break;
				}
			}
		}
		return true;
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
		constexpr char name_list[3] = { 'q', 'p', 'r' };
		constexpr char buffer_list[3] = { ':', ';', '|' };
		for (int i = 0; i < sum_momenta.size(); i++)
		{
			if (i >= 3) {
				throw std::invalid_argument("More than 3 momenta, time to implement this...");
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
				if (buffer.momentum_list[index].first > 0) buffer.flipMomentum();
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
		discardZeroMomenta();
	}

	void WickTerm::sort()
	{
		for (auto& delta : delta_momenta) {
			if (delta.first.momentum_list.size() == 1 && delta.second.momentum_list.size() == 1) {
				// This comparison is well defined because we save the momentum as char i.e. byte
				// which is easily comparable
				if (delta.first.momentum_list[0].second < delta.second.momentum_list[0].second) {
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
				for (auto& op : operators) {
					op.momentum.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
				}
				for (auto& coeff : coefficients) {
					coeff.momentum.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
				}
			}
		}

		for (auto& op : operators) {
			if (op.type == CDW_Type && op.momentum.add_Q) {
				op.momentum.add_Q = false;
				op.isDaggered = !(op.isDaggered);
			}
		}

		for (int i = 0; i < operators.size(); i++)
		{
			for (int j = i + 1; j < operators.size(); j++)
			{
				if (operators[i].type > operators[j].type) {
					std::swap(operators[i], operators[j]);
				}
				else if (operators[i].type == operators[j].type) {
					if (operators[i].momentum.momentum_list[0].second > operators[j].momentum.momentum_list[0].second) {
						std::swap(operators[i], operators[j]);
					}
					else if (operators[i].momentum.momentum_list[0].second == operators[j].momentum.momentum_list[0].second) {
						if (operators[i].momentum.add_Q && !(operators[j].momentum.add_Q)) {
							std::swap(operators[i], operators[j]);
						}
					}
				}
			}
		}

		for (auto& coeff : coefficients) {
			coeff.momentum.sort();
			if (coeff.translationalInvariance && coeff.momentum.momentum_list.size() > 0) {
				if (coeff.momentum.momentum_list[0].first < 0) {
					coeff.momentum.flipMomentum();
				}
			}
			if (coeff.Q_changes_sign && coeff.momentum.add_Q) {
				coeff.momentum.add_Q = false;
				this->multiplicity *= -1;
			}

			for (const auto& sum : sum_momenta) {
				int idx = coeff.momentum.isUsed(sum);
				if (idx < 0) continue;

				if (coeff.momentum.momentum_list[idx].first < 0) {
					for (auto& op : operators) {
						for (auto& mom : op.momentum.momentum_list)
						{
							if (mom.second == sum) mom.first *= -1;
						}
					}
					for (auto& coeff2 : coefficients) {
						for (auto& mom : coeff2.momentum.momentum_list)
						{
							if (mom.second == sum) mom.first *= -1;
						}
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
				op.indizes[0] = SpinUp;
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
			if (op.type == SC_Type || op.type == CDW_Type) {
				op.isDaggered = false;
			}
		}
	}

	void clearEtas(std::vector<WickTerm>& terms)
	{
		for (auto it = terms.begin(); it != terms.end();) {
			bool isEta = false;
			for (const auto& op : it->operators) {
				if (op.type == Eta_Type) {
					isEta = true;
					break;
				}
			}
			if (isEta) {
				it = terms.erase(it);
			}
			else {
				++it;
			}
		}
	}

	void cleanWicks(std::vector<WickTerm>& terms)
	{
		// Assuming (for now) that all <eta> = 0
		clearEtas(terms);
		for (auto& term : terms) {
			for (std::vector<Coefficient>::iterator it = term.coefficients.begin(); it != term.coefficients.end();) {
				if (it->name == "") {
					it = term.coefficients.erase(it);
				}
				else {
					++it;
				}
			}
		}
		for (std::vector<WickTerm>::iterator it = terms.begin(); it != terms.end();) {
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

			//it->applyPhaseSymmetry();
			//it->applySpinSymmetry();
			it->applyTranslationalSymmetry();

			for (auto jt = it->sum_indizes.begin(); jt != it->sum_indizes.end();)
			{
				if (it->usesIndex(*jt)) {
					++jt;
				}
				else {
					// We are assuming there are only spin indizes here (spin 1/2)
					// If another kind of index arises I have to readress this section.
					it->multiplicity *= 2;
					jt = it->sum_indizes.erase(jt);
				}
			}
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

		std::map<std::string, int> coeff_map;
		coeff_map["\\epsilon_0"] = 0;
		coeff_map["\\frac{U}{N}"] = 1;
		coeff_map["\\tilde{V}"] = 2;

		// Sort terms
		for (size_t i = 0; i < terms.size(); i++)
		{
			for (size_t j = i + 1; j < terms.size(); j++)
			{
				if (terms[i].delta_momenta.size() == 0 && terms[j].delta_momenta.size() > 0) {
					std::swap(terms[i], terms[j]);
				}
				if (terms[i].delta_momenta.size() > 0 && terms[j].delta_momenta.size() > 0) {
					if (terms[i].delta_momenta.size() < terms[j].delta_momenta.size()) {
						std::swap(terms[i], terms[j]);
					}
					else if (terms[i].delta_momenta.size() == terms[j].delta_momenta.size()) {
						if (terms[i].delta_momenta[0].second.add_Q && !(terms[j].delta_momenta[0].second.add_Q)) {
							std::swap(terms[i], terms[j]);
						}
						else if (terms[i].coefficients.size() > 0) {
							if (coeff_map.find(terms[j].coefficients[0].name)->second < coeff_map.find(terms[i].coefficients[0].name)->second) {
								std::swap(terms[i], terms[j]);
							}
						}
					}
				}
				else if (terms[i].delta_momenta.size() == 0 && terms[j].delta_momenta.size() == 0) {
					if (terms[i].coefficients.size() > 0) {
						if (coeff_map.find(terms[j].coefficients[0].name)->second < coeff_map.find(terms[i].coefficients[0].name)->second) {
							std::swap(terms[i], terms[j]);
						}
					}
				}
			}
		}
	}
}