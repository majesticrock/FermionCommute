#include "Coefficient.hpp"

namespace SymbolicOperators {
	std::ostream& operator<<(std::ostream& os, const Coefficient& coeff)
	{
		os << coeff.name;
		if (!coeff.indizes.empty()) {
			os << "_{ ";
			for (const auto& index : coeff.indizes) {
				os << index << " ";
			}
			os << "}";
		}
		if (coeff.isDaggered) {
			os << "^*";
		}
		if (!coeff.momentum.momentum_list.empty()) {
			os << " ( " << coeff.momentum << " )";
		}
		return os;
	}
	std::ostream& operator<<(std::ostream& os, const std::vector<Coefficient>& coeffs) {
		for (std::vector<Coefficient>::const_iterator it = coeffs.begin(); it != coeffs.end(); ++it)
		{
			os << (*it) << " ";
		}
		return os;
	}

	Coefficient::Coefficient()
		: name(""), momentum(), indizes(), isDaggered(false) {}
	Coefficient::Coefficient(std::string _name)
		: name(_name), momentum(), indizes(), isDaggered(false) {}
	Coefficient::Coefficient(std::string _name, const Momentum& _momentum, const std::vector<std::string>& _indizes, bool _isDaggered)
		: name(_name), momentum(_momentum), indizes(_indizes), isDaggered(_isDaggered) {}
	Coefficient::Coefficient(std::string _name, char _momentum, bool add_Q, const std::vector<std::string>& _indizes, bool _isDaggered)
		: name(_name), momentum(_momentum, add_Q), indizes(_indizes), isDaggered(_isDaggered) { }
	Coefficient::Coefficient(std::string _name, const Momentum& _momentum, bool _isDaggered)
		: name(_name), momentum(_momentum), indizes(), isDaggered(_isDaggered) {}
	Coefficient::Coefficient(std::string _name, char _momentum, bool add_Q, bool _isDaggered)
		: name(_name), momentum(_momentum, 1, add_Q), indizes(), isDaggered(_isDaggered) { }
}