#include "WickTerm.hpp"

#define LEFT temporary_operators[i]
#define RIGHT temporary_operators[i + 1]
#define L_SPIN temporary_operators[i].indizes[0]
#define R_SPIN temporary_operators[i + 1].indizes[0]

namespace SymbolicOperators {
	WickTerm::WickTerm(const Term* base)
		: multiplicity(base->multiplicity), coefficients(base->coefficients), sum_momenta(base->sum_momenta), sum_indizes(base->sum_indizes),
		operators(), delta_momenta(base->delta_momenta), delta_indizes(base->delta_indizes), temporary_operators()
	{
	}
	WickTerm::WickTerm()
		: multiplicity(0), coefficients(), sum_momenta(), sum_indizes(),
		operators(), delta_momenta(), delta_indizes(), temporary_operators()
	{
	}

	bool WickTerm::swapToWickOperators(std::vector<WickTerm>& reciever)
	{
		this->operators.reserve(temporary_operators.size() / 2);
		WickTerm this_copy = *this;
		std::vector<WickTerm> these_copies;
		these_copies.push_back(*this);

		auto setDeltas = [&](const Operator& left, const Operator& right, bool sc_type, int index) {
			Momentum copy_momentum = left.momentum;
			if (sc_type) copy_momentum.flipMomentum();

			for (int k = 1; k < left.indizes.size(); k++)
			{
				these_copies[index].delta_indizes.push_back(std::make_pair(left.indizes[k], right.indizes[k]));
				this_copy.delta_indizes.push_back(these_copies[index].delta_indizes.back());
			}

			these_copies[index].delta_momenta.push_back(std::make_pair(copy_momentum, right.momentum));
			copy_momentum.add_Q = !(copy_momentum.add_Q);
			this_copy.delta_momenta.push_back(std::make_pair(copy_momentum, right.momentum));
		};

		for (int i = 0; i < temporary_operators.size(); i += 2)
		{
			int copies_size = these_copies.size();
			for (int j = 0; j < copies_size; j++)
			{
				this_copy = these_copies[j];
				if (LEFT.isDaggered == RIGHT.isDaggered) {
					if (L_SPIN == R_SPIN) return false;
					if (LEFT.isDaggered) { // c^+ c^+
						if (L_SPIN == DOWN) {
							throw std::invalid_argument("c^+ c^+: Left spin is down while right isn't. Did you forget to sort the terms?");
						}
						if (L_SPIN != UP) {
							these_copies[j].delta_indizes.push_back(std::make_pair(L_SPIN, UP));
							this_copy.delta_indizes.push_back(these_copies[j].delta_indizes.back());
						}
						if (R_SPIN != DOWN) {
							these_copies[j].delta_indizes.push_back(std::make_pair(R_SPIN, DOWN));
							this_copy.delta_indizes.push_back(these_copies[j].delta_indizes.back());
						}
						// Due to the dagger we need to swap left and right
						setDeltas(RIGHT, LEFT, true, j);

						these_copies[j].operators.push_back(WickOperator("f", true, LEFT.momentum));
						if (LEFT.indizes.size() > 1) {
							these_copies[j].operators.back().indizes = std::vector<std::string>(LEFT.indizes.begin() + 1, LEFT.indizes.end());
						}
						this_copy.operators.push_back(these_copies[j].operators.back());
						this_copy.operators.back().type = "\\eta";
					}
					else { // cc
						if (L_SPIN == UP) {
							throw std::invalid_argument("c^+ c^+: Left spin is down while right isn't. Did you forget to sort the terms?");
						}
						if (L_SPIN != DOWN) {
							these_copies[j].delta_indizes.push_back(std::make_pair(L_SPIN, DOWN));
							this_copy.delta_indizes.push_back(these_copies[j].delta_indizes.back());
						}
						if (R_SPIN != UP) {
							these_copies[j].delta_indizes.push_back(std::make_pair(R_SPIN, UP));
							this_copy.delta_indizes.push_back(these_copies[j].delta_indizes.back());
						}
						setDeltas(LEFT, RIGHT, true, j);

						these_copies[j].operators.push_back(WickOperator("f", false, RIGHT.momentum));
						if (RIGHT.indizes.size() > 1) {
							these_copies[j].operators.back().indizes = std::vector<std::string>(RIGHT.indizes.begin() + 1, RIGHT.indizes.end());
						}
						this_copy.operators.push_back(these_copies[j].operators.back());
						this_copy.operators.back().type = "\\eta";
					}
				}
				else {
					// c^+ c
					if (L_SPIN == UP && R_SPIN == DOWN) return false;
					if (L_SPIN == DOWN && R_SPIN == UP) return false;

					if (L_SPIN != R_SPIN) {
						these_copies[j].delta_indizes.push_back(std::make_pair(L_SPIN, R_SPIN));
						this_copy.delta_indizes.push_back(these_copies[j].delta_indizes.back());
					}
					// Left and right are swapped due to the definition of g
					setDeltas(RIGHT, LEFT, false, j);

					these_copies[j].operators.push_back(WickOperator("n", false, LEFT.momentum, LEFT.indizes));
					this_copy.operators.push_back(these_copies[j].operators.back());
					this_copy.operators.back().type = "g";
					if (this_copy.operators.back().momentum.add_Q) {
						this_copy.operators.back().momentum.add_Q = false;
						this_copy.operators.back().isDaggered = true;
					}
				}
				these_copies.push_back(this_copy);
			}
		}
		this->delta_indizes = these_copies.back().delta_indizes;
		this->delta_momenta = these_copies.back().delta_momenta;
		this->operators = these_copies.back().operators;
		these_copies.pop_back();

		reciever.insert(reciever.end(), these_copies.begin(), these_copies.end());
		return true;
	}

	std::ostream& operator<<(std::ostream& os, const WickOperator& op)
	{
		os << "\\langle " << op.type << "_{ " << op.momentum << ", ";
		for (const auto& index : op.indizes) {
			os << index << " ";
		}
		os << "}";
		if (op.isDaggered) {
			os << "^\\dagger";
		}
		os << " \\rangle";
		return os;
	}
	std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops)
	{
		for (const auto& op : ops) {
			os << op << " ";
		}
		return os;
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

	WickOperator::WickOperator(const std::string& _type, const bool _isDaggered, const Momentum& _momentum, const std::vector<std::string>& _indizes)
		: type(_type), isDaggered(_isDaggered), momentum(_momentum), indizes(_indizes) {}
	WickOperator::WickOperator(const std::string& _type, const bool _isDaggered, const Momentum& _momentum, const std::string& _index)
		: type(_type), isDaggered(_isDaggered), momentum(_momentum), indizes(1, _index) {}
	WickOperator::WickOperator()
		: type(""), isDaggered(false), momentum(), indizes() {}
}