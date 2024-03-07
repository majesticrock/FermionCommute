#pragma once
#include "Momentum.hpp"
#include <IndexWrapper.hpp>

namespace SymbolicOperators {
	//const std::string UP = "\\uparrow";
	//const std::string DOWN = "\\downarrow";

	struct Operator {
		Momentum momentum;
		// Contains all indizes, standard: first index = spin, all others arbitrary, e.g., orbitals, bands etc
		IndexWrapper indizes;
		bool isDaggered;

		Operator(const Momentum& _momentum, const IndexWrapper _indizes, bool _isDaggered);
		Operator(const momentum_pairs& _momentum, const IndexWrapper _indizes, bool _isDaggered);
		Operator(char _momentum, bool add_Q, const IndexWrapper _indizes, bool _isDaggered);
		Operator(char _momentum, int sign, bool add_Q, const IndexWrapper _indizes, bool _isDaggered);

		inline void hermitianConjugate() {
			this->isDaggered = !(this->isDaggered);
		}
	};

	inline bool operator==(const Operator& lhs, const Operator& rhs) {
		if (lhs.isDaggered != rhs.isDaggered) return false;
		if(lhs.indizes != rhs.indizes) return false;
		return (lhs.momentum == rhs.momentum);
	}
	inline bool operator!=(const Operator& lhs, const Operator& rhs) {
		return !(lhs == rhs);
	}

	std::ostream& operator<<(std::ostream& os, const Operator& op);
}