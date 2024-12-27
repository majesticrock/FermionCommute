#include "HubbardDispersions.hpp"

namespace mrock::SymbolicOperators {
	std::vector<std::vector<Term>> HubbardDispersions::XP_basis() const
	{
		std::vector<std::vector<Term>> basis = Hubbard::XP_basis();
		for (auto& basis_term : basis) {
			if (basis_term.front().operators.front().is_daggered) {
				basis_term.front().operators.front().momentum += Momentum('x');
			}
			else {
				basis_term.front().operators.front().momentum += Momentum('x', -1);
			}
			if (basis_term.size() == 2U) {
				basis_term[1] = basis_term.front();
			}
			else {
				basis_term.push_back(basis_term.front());
			}
			basis_term[1].hermitianConjugate();
		}
		return basis;
	}
	std::vector<std::vector<Term>> HubbardDispersions::STD_basis() const
	{
		std::vector<std::vector<Term>> ret = Hubbard::STD_basis();
		for (auto& _v : ret) {
			for (auto& v : _v) {
				if (v.operators.front().is_daggered) {
					v.operators.front().momentum += Momentum('x');
				}
				else {
					v.operators.front().momentum += Momentum('x', -1);
				}
			}
		}
		return ret;
	}

	std::string HubbardDispersions::get_subfolder() const
	{
		return "hubbard/dispersions/";
	}
}