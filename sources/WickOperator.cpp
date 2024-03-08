#include "WickOperator.hpp"

namespace SymbolicOperators {
    WickOperator::WickOperator(const OperatorType& _type, const bool _isDaggered, const Momentum& _momentum, const IndexWrapper& _indizes)
		: type(_type), isDaggered(_isDaggered), momentum(_momentum), indizes(_indizes) {}
	WickOperator::WickOperator(const OperatorType& _type, const bool _isDaggered, const Momentum& _momentum, const Index _index)
		: type(_type), isDaggered(_isDaggered), momentum(_momentum), indizes(_index) {}
	WickOperator::WickOperator()
		: type(Undefined_Type), isDaggered(false), momentum(), indizes() {}

	std::ostream& operator<<(std::ostream& os, const OperatorType op)
	{
		switch (op) {
		case SC_Type:
			os << "f";
			break;
		case Eta_Type:
			os << "\\eta";
			break;
		case CDW_Type:
			os << "g";
			break;
		case Number_Type:
			os << "n";
			break;
		default:
			os << "ERROR_OPERATOR";
		}
		return os;
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
}