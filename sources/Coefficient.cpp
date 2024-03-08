#include "Coefficient.hpp"

namespace SymbolicOperators {
	std::ostream& operator<<(std::ostream& os, const Coefficient& coeff)
	{
		os << coeff.name;
		if (!coeff.indizes.empty()) {
			os << "_{ " << coeff.indizes << "}";
		}
		if (coeff.isDaggered) {
			os << "^*";
		}
		if (!coeff.momentum.momentum_list.empty() || coeff.momentum.add_Q) {
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
		: name(""), momentum(), indizes(), Q_changes_sign(false), isDaggered(false) {}
	Coefficient::Coefficient(std::string _name)
		: name(_name), momentum(), indizes(), Q_changes_sign(false), isDaggered(false) {}
	Coefficient::Coefficient(std::string _name, const Momentum& _momentum, const IndexWrapper& _indizes, bool _Q_changes_sign, bool _isDaggered)
		: name(_name), momentum(_momentum), indizes(_indizes), Q_changes_sign(_Q_changes_sign), isDaggered(_isDaggered) {}
	Coefficient::Coefficient(std::string _name, char _momentum, bool add_Q, const IndexWrapper& _indizes, bool _Q_changes_sign, bool _isDaggered)
		: name(_name), momentum(_momentum, add_Q), indizes(_indizes), Q_changes_sign(_Q_changes_sign), isDaggered(_isDaggered) { }
	Coefficient::Coefficient(std::string _name, const Momentum& _momentum, bool _Q_changes_sign, bool _isDaggered)
		: name(_name), momentum(_momentum), indizes(), Q_changes_sign(_Q_changes_sign), isDaggered(_isDaggered) {}
	Coefficient::Coefficient(std::string _name, char _momentum, bool add_Q, bool _Q_changes_sign, bool _isDaggered)
		: name(_name), momentum(_momentum, 1, add_Q), indizes(), Q_changes_sign(_Q_changes_sign), isDaggered(_isDaggered) { }
}