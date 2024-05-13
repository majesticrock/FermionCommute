#include "Continuum.hpp"

namespace SymbolicOperators {
	std::vector<Term> Continuum::hamiltonian() const
	{
		const Term H_T(1, Coefficient("\\epsilon_0", Momentum('q')), SumContainer{ MomentumSum({ 'q' }), Sigma },
			std::vector<Operator>({
				Operator('q', 1, false, Sigma, true), Operator('q', 1, false, Sigma, false)
				}));

		const Term H_U(-1, Coefficient("U", MomentumList({ 'r', 'q' })), MomentumSum({ 'r', 'p', 'q' }), std::vector<Operator>({
			Operator('r', 1, false, SpinUp, true), Operator('p', 1, false, SpinDown, true),
			Operator(Momentum("r+p-q"), SpinDown, false),
			Operator(momentum_pairs({ std::make_pair(1, 'q') }), SpinUp, false),
			}));

		return { H_T, H_U };
	}
	std::vector<WickOperatorTemplate> Continuum::templates() const
	{
		return {
			WickOperatorTemplate{ {IndexComparison{false, SpinDown, SpinUp}}, Momentum(), SC_Type, true },
			WickOperatorTemplate{ {IndexComparison{true}}, Momentum(), Number_Type, false }
		};
	}
	std::vector<std::vector<Term>> Continuum::XP_basis() const
	{
		return {
			// 0: f + f^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k, c_k })),
				Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
				}),
			// 1/2: n_up/down
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k }))
				}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k }))
				}),
			// 3: f - f^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k, c_k })),
				Term(-1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
				})
		};
	}
	std::vector<std::vector<Term>> Continuum::STD_basis() const
	{
		return {
			// f, f^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k, c_k }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
			}),
			// n_up/down
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k }))
			})
		};
	}
	std::string Continuum::get_subfolder() const
	{
		return "continuum/";
	}
}