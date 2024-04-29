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
		os << coeff.momenta << " ";
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
		: name(""), momenta(), indizes(), Q_changes_sign(false), isDaggered(false) {}
	Coefficient::Coefficient(std::string _name)
		: name(_name), momenta(), indizes(), Q_changes_sign(false), isDaggered(false) {}
	Coefficient::Coefficient(std::string _name, const Momentum& _momentum, const IndexWrapper& _indizes, bool _Q_changes_sign, bool _isDaggered)
		: name(_name), momenta(_momentum), indizes(_indizes), Q_changes_sign(_Q_changes_sign), isDaggered(_isDaggered) {}
	Coefficient::Coefficient(std::string _name, const Momentum& _momentum, bool _Q_changes_sign, bool _isDaggered)
		: name(_name), momenta(_momentum), indizes(), Q_changes_sign(_Q_changes_sign), isDaggered(_isDaggered) {}
	Coefficient::Coefficient(std::string _name, const MomentumList& _momenta, const IndexWrapper& _indizes, bool _Q_changes_sign, bool _isDaggered)
		: name(_name), momenta(_momenta), indizes(), Q_changes_sign(_Q_changes_sign), isDaggered(_isDaggered) { }
}