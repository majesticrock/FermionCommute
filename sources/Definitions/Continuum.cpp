#include "Continuum.hpp"

namespace SymbolicOperators {
	std::vector<Term> Continuum::hamiltonian() const
	{
		const Term H_Kin(1, Coefficient("\\epsilon_0", Momentum('q')), SumContainer{ MomentumSum({ 'q' }), Index::Sigma },
			std::vector<Operator>({
				Operator('q', 1, false, Index::Sigma, true), Operator('q', 1, false, Index::Sigma, false)
				}));

		const Term H_Ph(-1, Coefficient("g", MomentumList({ 'q', 'p' })), SumContainer{ MomentumSum({'q', 'p'}) }, std::vector<Operator>({
			c_k_dagger.with_momentum('q'), c_minus_k_dagger.with_momentum('q'),
			c_minus_k.with_momentum('p'), c_k.with_momentum('p')
			}));

		const Term H_C(IntFractional(1, 2), Coefficient("V", Momentum('q')),
			SumContainer{ MomentumSum({ 'r', 'p', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator('r', 1, false, Index::Sigma, true),
				Operator('p', 1, false, Index::SigmaPrime, true),
				Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), Index::SigmaPrime, false),
				Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), Index::Sigma, false),
				}));

		/* SC channel only
		const Term H_C(IntFractional(1, 1), Coefficient("V", Momentum('q')),
			SumContainer{ MomentumSum({ 'p', 'q' }) },
			std::vector<Operator>({
				Operator('p', 1, false, Index::SpinUp, true),
				Operator('p', -1, false, Index::SpinDown, true),
				Operator(momentum_pairs({ std::make_pair(-1, 'p'), std::make_pair(-1, 'q') }), Index::SpinDown, false),
				Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(1, 'q') }), Index::SpinUp, false),
				}));
		*/

		const Term H_BG(-1, Coefficient("\\rho"), SumContainer{ MomentumSum({ 'q' }), Index::Sigma },
			std::vector<Operator>({
				Operator(Momentum("q"), Index::Sigma, true), Operator('q', 1, false, Index::Sigma, false)
				}));

		return { H_Kin, H_Ph, H_C, H_BG };
	}
	std::vector<WickOperatorTemplate> Continuum::templates() const
	{
		return {
			WickOperatorTemplate{ {IndexComparison{false, Index::SpinDown, Index::SpinUp}}, Momentum(), SC_Type, true },
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
				// 1: n_up + down
				std::vector<Term>({
					Term(1, std::vector<Operator>({ c_k_dagger, c_k })),
					Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k }))
					}),
				// 2: f - f^+
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
				// n_up + down
				std::vector<Term>({
					Term(1, std::vector<Operator>({ c_k_dagger, c_k })),
					Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k }))
				})
		};
	}
	std::vector<std::unique_ptr<WickSymmetry>> Continuum::symmetries() const
	{
		std::vector<std::unique_ptr<WickSymmetry>> ret;
		ret.push_back(std::make_unique<SpinSymmetry>());
		ret.push_back(std::make_unique<TranslationalSymmetry>());
		ret.push_back(std::make_unique<PhaseSymmetry<SC_Type>>());
		return ret;
	}
	std::string Continuum::get_subfolder() const
	{
		return "continuum/";
	}
}