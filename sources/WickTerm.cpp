#include "WickTerm.hpp"
#include <variant>

#define LEFT temporary_operators[i]
#define RIGHT temporary_operators[i + 1]
#define L_SPIN temporary_operators[i].indizes[0]
#define R_SPIN temporary_operators[i + 1].indizes[0]

namespace SymbolicOperators {
	WickTerm::WickTerm(const Term* base)
		: multiplicity(base->multiplicity), coefficients(base->coefficients), sum_momenta(base->sum_momenta),
		sum_indizes(base->sum_indizes), operators(), delta_momenta(base->delta_momenta),
		delta_indizes(base->delta_indizes), temporary_operators()
	{
	}
	WickTerm::WickTerm(const Term& base)
		: multiplicity(base.multiplicity), coefficients(base.coefficients), sum_momenta(base.sum_momenta),
		sum_indizes(base.sum_indizes), operators(), delta_momenta(base.delta_momenta),
		delta_indizes(base.delta_indizes), temporary_operators()
	{
	}
	WickTerm::WickTerm()
		: multiplicity(0), coefficients(), sum_momenta(), sum_indizes(),
		operators(), delta_momenta(), delta_indizes(), temporary_operators()
	{
	}

	void wick_processor(const std::vector<Operator>& remaining, std::vector<WickTerm>& reciever_list, std::variant<WickTerm, Term> buffer)
	{
		if (remaining.empty()) {
			reciever_list.push_back(std::get<WickTerm>(buffer));
			return;
		}
		for (size_t i = 1U; i < remaining.size(); ++i)
		{
			if (std::holds_alternative<Term>(buffer)) {
				WickTerm temp(std::get<Term>(buffer));
				buffer = temp;
			}
			if ((i % 2) == 0) {
				std::get<WickTerm>(buffer).multiplicity *= -1;
			}
			std::get<WickTerm>(buffer).temporary_operators.reserve(std::get<WickTerm>(buffer).temporary_operators.size() + 2);
			std::get<WickTerm>(buffer).temporary_operators.push_back(remaining[0]);
			std::get<WickTerm>(buffer).temporary_operators.push_back(remaining[i]);

			std::vector<Operator> copy_operators = remaining;
			copy_operators.erase(copy_operators.begin() + i);
			copy_operators.erase(copy_operators.begin());
			wick_processor(copy_operators, reciever_list, buffer);

			// delete last two elements, as they are to be updated in the next iteration
			std::get<WickTerm>(buffer).temporary_operators.pop_back();
			std::get<WickTerm>(buffer).temporary_operators.pop_back();
			if ((i % 2) == 0) {
				std::get<WickTerm>(buffer).multiplicity *= -1;
			}
		}
	}

	void wicks_theorem(const Term& term, std::vector<WickTerm>& reciever)
	{
		//if (this->operators.size() > 4) throw std::length_error("Wick for n>4 not yet implemented!");
		if (term.isIdentity()) {
			reciever.push_back(WickTerm(term));
		}
		else {
			std::vector<WickTerm> buffer_list;
			{
				size_t value = 1;
				// Computes the double factorial; total number of products in wicks theorem = (2n - 1)!!
				// Needs to stay signed in case of size() being 0
				for (long n = 2 * term.getOperators().size() - 1; n > 0; n -= 2)
				{
					value *= n;
				}
				buffer_list.reserve(value);
			}
			wick_processor(term.getOperators(), buffer_list, term);
			reciever.insert(reciever.end(), buffer_list.begin(), buffer_list.end());
		}

		for (size_t i = 0U; i < reciever.size();)
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

	bool WickTerm::swapToWickOperators(std::vector<WickTerm>& reciever)
	{
		this->operators.reserve(temporary_operators.size() / 2);
		WickTerm this_copy = *this;
		std::vector<WickTerm> these_copies;
		these_copies.push_back(*this);

		auto setDeltas = [&](const Operator& left, const Operator& right, size_t index) {
			Momentum copy_momentum = left.momentum;
			if (left.isDaggered == right.isDaggered) copy_momentum.flipMomentum();

			for (size_t k = 1U; k < left.indizes.size(); ++k)
			{
				these_copies[index].delta_indizes.push_back(make_delta(left.indizes[k], right.indizes[k]));
				this_copy.delta_indizes.push_back(these_copies[index].delta_indizes.back());
			}

			these_copies[index].delta_momenta.push_back(make_delta(copy_momentum, right.momentum));
			copy_momentum.add_Q = !(copy_momentum.add_Q);
			this_copy.delta_momenta.push_back(make_delta(copy_momentum, right.momentum));
			};

		for (size_t i = 0U; i < temporary_operators.size(); i += 2U)
		{
			size_t copies_size = these_copies.size();
			for (size_t j = 0U; j < copies_size; ++j)
			{
				this_copy = these_copies[j];
				if (LEFT.isDaggered == RIGHT.isDaggered) {
					if (L_SPIN == R_SPIN) return false;
					if (LEFT.isDaggered) { // c^+ c^+
						if (L_SPIN == SpinDown) {
							throw std::invalid_argument("c^+ c^+: Left spin is down while right isn't. Did you forget to sort the terms?");
						}
						if (L_SPIN != SpinUp) {
							these_copies[j].delta_indizes.push_back(make_delta(L_SPIN, SpinUp));
							this_copy.delta_indizes.push_back(these_copies[j].delta_indizes.back());
						}
						if (R_SPIN != SpinDown) {
							these_copies[j].delta_indizes.push_back(make_delta(R_SPIN, SpinDown));
							this_copy.delta_indizes.push_back(these_copies[j].delta_indizes.back());
						}
						// Due to the dagger we need to swap left and right
						setDeltas(RIGHT, LEFT, j);

						these_copies[j].operators.push_back(WickOperator(SC_Type, true, LEFT.momentum));
						if (LEFT.indizes.size() > 1) {
							these_copies[j].operators.back().indizes = std::vector<Index>(LEFT.indizes.begin() + 1, LEFT.indizes.end());
						}
						this_copy.operators.push_back(these_copies[j].operators.back());
						this_copy.operators.back().type = Eta_Type;
					}
					else { // cc
						if (L_SPIN == SpinUp) {
							throw std::invalid_argument("c^+ c^+: Left spin is down while right isn't. Did you forget to sort the terms?");
						}
						if (L_SPIN != SpinDown) {
							these_copies[j].delta_indizes.push_back(make_delta(L_SPIN, SpinDown));
							this_copy.delta_indizes.push_back(these_copies[j].delta_indizes.back());
						}
						if (R_SPIN != SpinUp) {
							these_copies[j].delta_indizes.push_back(make_delta(R_SPIN, SpinUp));
							this_copy.delta_indizes.push_back(these_copies[j].delta_indizes.back());
						}
						setDeltas(LEFT, RIGHT, j);

						these_copies[j].operators.push_back(WickOperator(SC_Type, false, RIGHT.momentum));
						if (RIGHT.indizes.size() > 1) {
							these_copies[j].operators.back().indizes = std::vector<Index>(RIGHT.indizes.begin() + 1, RIGHT.indizes.end());
						}
						this_copy.operators.push_back(these_copies[j].operators.back());
						this_copy.operators.back().type = Eta_Type;
					}
				}
				else {
					// c^+ c
					if (L_SPIN == SpinUp && R_SPIN == SpinDown) return false;
					if (L_SPIN == SpinDown && R_SPIN == SpinUp) return false;

					if (L_SPIN != R_SPIN) {
						these_copies[j].delta_indizes.push_back(make_delta(L_SPIN, R_SPIN));
						this_copy.delta_indizes.push_back(these_copies[j].delta_indizes.back());
					}
					// Left and right are swapped due to the definition of g
					setDeltas(RIGHT, LEFT, j);

					these_copies[j].operators.push_back(WickOperator(Number_Type, false, LEFT.momentum, LEFT.indizes));
					this_copy.operators.push_back(these_copies[j].operators.back());
					this_copy.operators.back().type = CDW_Type;
					if (this_copy.operators.back().momentum.add_Q) {
						this_copy.operators.back().momentum.add_Q = false;
						this_copy.operators.back().isDaggered = true;
					}
				}
				these_copies.push_back(this_copy);
			}
		}
		delta_indizes = these_copies.back().delta_indizes;
		delta_momenta = these_copies.back().delta_momenta;
		this->operators = these_copies.back().operators;
		these_copies.pop_back();

		reciever.insert(reciever.end(), these_copies.begin(), these_copies.end());
		return true;
	}

	std::ostream& operator<<(std::ostream& os, const WickTerm& term)
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
	std::ostream& operator<<(std::ostream& os, const std::vector<WickTerm>& terms)
	{
		for (std::vector<WickTerm>::const_iterator it = terms.begin(); it != terms.end(); ++it)
		{
			os << "\t&" << *it;
			if (it != terms.end() - 1) {
				os << " \\\\";
			}
			os << "\n";
		}
		return os;
	}
}