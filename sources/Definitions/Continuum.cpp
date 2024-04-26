#include "Continuum.hpp"

namespace SymbolicOperators {
	std::vector<Term> Continuum::hamiltonian() const
	{
		const Term H_T(1, Coefficient("\\epsilon_0", 'q'), SumContainer{ MomentumSum({ 'q' }), Sigma },
			std::vector<Operator>({
				Operator('q', 1, false, Sigma, true), Operator('q', 1, false, Sigma, false)
				}));

		const Term H_U(1, Coefficient("U", 'q'), SumContainer{ MomentumSum({'q', 'p'}) }, std::vector<Operator>({
			c_k_dagger.with_momentum('q'), c_minus_k_dagger.with_momentum('q'),
			c_minus_k.with_momentum('p'), c_k.with_momentum('p')
			}));
		return std::vector<Term>();
	}
	std::vector<WickOperatorTemplate> Continuum::templates() const
	{
		return std::vector<WickOperatorTemplate>();
	}
	std::vector<std::vector<Term>> Continuum::XP_basis() const
	{
		return std::vector<std::vector<Term>>();
	}
	std::vector<std::vector<Term>> Continuum::STD_basis() const
	{
		return std::vector<std::vector<Term>>();
	}
}