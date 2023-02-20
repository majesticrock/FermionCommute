#pragma once
#include <string>
#include <iostream>
#include "Momentum.hpp"

namespace SymbolicOperators {
	const std::string UP = "\\uparrow";
	const std::string DOWN = "\\downarrow";

	struct Operator {
		Momentum momentum;
		// Contains all indizes, standard: first index = spin, all others arbitrary, e.g. orbitals, bands etc
		std::vector<std::string> indizes;
		bool isDaggered;

		Operator(const Momentum& _momentum, const std::vector<std::string>& _indizes, bool _isDaggered);
		Operator(const Momentum& _momentum, const std::string& index, bool _isDaggered);
		Operator(const momentum_pairs& _momentum, const std::string& index, bool _isDaggered);
		Operator(char _momentum, bool add_Q, const std::vector<std::string>& _indizes, bool _isDaggered);
		Operator(char _momentum, int sign, bool add_Q, const std::string& index, bool _isDaggered);

		inline void hermitianConjugate() {
			this->isDaggered = !(this->isDaggered);
		}
	};

	inline bool operator==(const Operator& lhs, const Operator& rhs) {
		if (lhs.isDaggered != rhs.isDaggered) {
			return false;
		}
		if (lhs.indizes.size() != rhs.indizes.size()) {
			// Should never occur, but better to check anyways
			std::cout << "Warning: Number of indizes between two operators does not match!" << std::endl;
			return false;
		}
		for (size_t i = 0; i < lhs.indizes.size(); i++)
		{
			if (lhs.indizes[i] != rhs.indizes[i]) {
				return false;
			}
		}
		return (lhs.momentum == rhs.momentum);
	}
	inline bool operator!=(const Operator& lhs, const Operator& rhs) {
		return !(lhs == rhs);
	}

	std::ostream& operator<<(std::ostream& os, const Operator& op);
}