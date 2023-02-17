#include "Operator.hpp"

namespace SymbolicOperators {
	std::ostream& operator<<(std::ostream& os, const Operator& op)
	{
		os << "c_{ " << op.momentum << ", ";
		for (const auto& index : op.indizes) {
			os << index << " ";
		}
		os << "}";
		if (op.isDaggered) {
			os << "^\\dagger ";
		}
		return os;
	}

	Operator::Operator(const Momentum& _momentum, const std::vector<std::string>& _indizes, bool _isDaggered)
		: momentum(_momentum), indizes(_indizes), isDaggered(_isDaggered) {}

	Operator::Operator(const Momentum& _momentum, const std::string& index, bool _isDaggered)
		: momentum(_momentum), indizes(1, index), isDaggered(_isDaggered) {}

	Operator::Operator(const momentum_pairs& _momentum, const std::string& index, bool _isDaggered)
		: momentum(_momentum), indizes(1, index), isDaggered(_isDaggered) {}

	Operator::Operator(char _momentum, bool add_Q, const std::vector<std::string>& _indizes, bool _isDaggered)
		: momentum(_momentum, add_Q), indizes(_indizes), isDaggered(_isDaggered) {}

	Operator::Operator(char _momentum, int sign, bool add_Q, const std::string& index, bool _isDaggered)
		: momentum(_momentum, sign, add_Q), indizes(1, index), isDaggered(_isDaggered) {}
}