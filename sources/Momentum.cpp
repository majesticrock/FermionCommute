#include "Momentum.hpp"

namespace SymbolicOperators {
	Momentum::Momentum()
		: momentum_list(), add_Q(false) {}
	Momentum::Momentum(const char value, int plus_minus, bool Q)
		: momentum_list(1, std::make_pair(plus_minus, value)), add_Q(Q) {}
	Momentum::Momentum(const momentum_pairs& _momenta, bool Q)
		: momentum_list(_momenta), add_Q(Q) {}

	void Momentum::sort()
	{
		for (size_t i = 0; i < momentum_list.size(); i++)
		{
			for (size_t j = i + 1; j < momentum_list.size(); j++)
			{
				// Comparing two chars is easy
				if (momentum_list[i].second > momentum_list[j].second) {
					std::swap(momentum_list[i], momentum_list[j]);
				}
			}
		}
	}

	Momentum& Momentum::operator+=(const Momentum& rhs)
	{
		this->add_Q = (rhs.add_Q != this->add_Q);
		bool foundOne = false;
		for (size_t i = 0; i < rhs.momentum_list.size(); i++)
		{
			foundOne = false;
			for (size_t j = 0; j < this->momentum_list.size(); j++)
			{
				if (rhs.momentum_list[i].second == this->momentum_list[j].second) {
					foundOne = true;
					this->momentum_list[j].first += rhs.momentum_list[i].first;
					if (this->momentum_list[j].first == 0) {
						this->momentum_list.erase(this->momentum_list.begin() + j);
					}
					break;
				}
			}
			if (!foundOne) {
				this->momentum_list.push_back(rhs.momentum_list[i]);
			}
		}

		return *this;
	}
	Momentum& Momentum::operator-=(const Momentum& rhs)
	{
		this->add_Q = (rhs.add_Q != this->add_Q);
		bool foundOne = false;
		for (size_t i = 0; i < rhs.momentum_list.size(); i++)
		{
			foundOne = false;
			for (size_t j = 0; j < this->momentum_list.size(); j++)
			{
				if (rhs.momentum_list[i].second == this->momentum_list[j].second) {
					foundOne = true;
					this->momentum_list[j].first -= rhs.momentum_list[i].first;
					if (this->momentum_list[j].first == 0) {
						this->momentum_list.erase(this->momentum_list.begin() + j);
					}
					break;
				}
			}
			if (!foundOne) {
				this->momentum_list.push_back(std::make_pair(-rhs.momentum_list[i].first, rhs.momentum_list[i].second));
			}
		}

		return *this;
	}
	void Momentum::addInPlace(const Momentum& rhs)
	{
		(*this) += rhs;
	}

	void Momentum::replaceOccurances(const char replaceWhat, const Momentum& replaceWith)
	{
		for (const auto& x : replaceWith.momentum_list) {
			if (x.second == replaceWhat) {
				throw std::invalid_argument("You are trying to replace a momentum with itself. This has undefined behaviour!");
			}
		}
		for (size_t i = 0; i < momentum_list.size(); ++i) {
			if (momentum_list[i].second == replaceWhat) {
				auto buffer = replaceWith;
				buffer.multiplyMomentum(momentum_list[i].first);
				this->momentum_list.erase(momentum_list.begin() + i);

				(*this) += buffer;
			}
		}
	}

	std::ostream& operator<<(std::ostream& os, const Momentum& momentum)
	{
		if (momentum.momentum_list.empty()) {
			if (momentum.add_Q) {
				os << "Q";
			}
			else {
				os << "0";
			}
			return os;
		}
		for (momentum_pairs::const_iterator it = momentum.momentum_list.begin(); it != momentum.momentum_list.end(); ++it)
		{
			if (it != momentum.momentum_list.begin() && it->first > 0) {
				os << " +";
			}
			else if (it->first == -1) {
				os << " -";
			}
			else if (std::abs(it->first) != 1) {
				os << it->first;
			}
			os << it->second;
		}
		if (momentum.add_Q) {
			os << " + Q";
		}
		return os;
	}
}