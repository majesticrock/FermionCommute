#include "Hubbard.hpp"

namespace mrock::symbolic_operators {
	std::vector<Term> Hubbard::hamiltonian() const
	{
		const Term H_T(1, Coefficient("\\epsilon_0", Momentum('q')), SumContainer{ MomentumSum({ 'q' }), Index::Sigma },
			std::vector<Operator>({
				Operator('q', 1, false, Index::Sigma, true), Operator('q', 1, false, Index::Sigma, false)
				}));

		const Term H_U(1, Coefficient("\\frac{U}{N}"), MomentumSum({ 'r', 'p', 'q' }), std::vector<Operator>({
			Operator('r', 1, false, Index::SpinUp, true), Operator('p', 1, false, Index::SpinDown, true),
			Operator(momentum_symbols({ MomentumSymbol(1, 'p'), MomentumSymbol(-1, 'q') }), Index::SpinDown, false),
			Operator(momentum_symbols({ MomentumSymbol(1, 'r'), MomentumSymbol(1, 'q') }), Index::SpinUp, false),
			}));

		const Term H_V(1, Coefficient("\\tilde{V}", Momentum('q'), true),
			SumContainer{ MomentumSum({ 'r', 'p', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator('r', 1, false, Index::Sigma, true),
				Operator('p', 1, false, Index::SigmaPrime, true),
				Operator(momentum_symbols({ MomentumSymbol(1, 'p'), MomentumSymbol(-1, 'q') }), Index::SigmaPrime, false),
				Operator(momentum_symbols({ MomentumSymbol(1, 'r'), MomentumSymbol(1, 'q') }), Index::Sigma, false),
				}));

		return { H_T, H_U, H_V };//
	}
	std::vector<WickOperatorTemplate> Hubbard::templates() const
	{
		return {
			WickOperatorTemplate{ {IndexComparison{false, Index::SpinDown, Index::SpinUp}}, Momentum(), SC_Type, true },
			WickOperatorTemplate{ {IndexComparison{false, Index::SpinDown, Index::SpinUp}}, Momentum(momentum_symbols(), true), Eta_Type, true },
			WickOperatorTemplate{ {IndexComparison{true}}, Momentum(), Number_Type, false },
			WickOperatorTemplate{ {IndexComparison{true}}, Momentum(momentum_symbols(), true), CDW_Type, false }
		};
	}
	std::vector<std::vector<Term>> Hubbard::XP_basis() const
	{
		return {
			// 0: phi + phi^+
			//std::vector<Term>({
			//	Term(1, std::vector<Operator>({ c_k_down_dagger, c_k })),
			//	Term(1, std::vector<Operator>({ c_k_dagger, c_k_down }))
			//	}),
			// 0: f + f^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k, c_k })),
				Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
				}),
			// 1: eta + eta^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_Q, c_k })),
				Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_Q_dagger }))
				}),
			// 2/3: g_up/down +
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q })),
				Term(1, std::vector<Operator>({ c_k_Q_dagger, c_k }))
				}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k_Q })),
				Term(1, std::vector<Operator>({ c_minus_k_Q_dagger, c_minus_k }))
				}),
			// 4: transversal magnon, hermitian
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q_down })),
				Term(1, std::vector<Operator>({ c_k_Q_down_dagger, c_k }))
				}),
			// 5/6: n_up/down
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k }))
				}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k }))
				}),
			// 7: f - f^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k, c_k })),
				Term(-1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
				}),
			// 8: eta - eta^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_Q, c_k })),
				Term(-1, std::vector<Operator>({ c_k_dagger, c_minus_k_Q_dagger }))
				}),
			// 9/10: g_up/down -
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q })),
				Term(-1, std::vector<Operator>({ c_k_Q_dagger, c_k }))
				}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k_Q })),
				Term(-1, std::vector<Operator>({ c_minus_k_Q_dagger, c_minus_k }))
				}),
			// 11: transversal magnon, antihermitian
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q_down })),
				Term(-1, std::vector<Operator>({ c_k_Q_down_dagger, c_k }))
				})
		};
	}
	std::vector<std::vector<Term>> Hubbard::STD_basis() const
	{
		return {
			// 0/1: f, f^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k, c_k }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
			}),
			// 2/3: n_up/down
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k }))
			}),
			// 4/5: g_up/down
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k_Q }))
			}),
			// 6/7: eta, eta^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_Q, c_k }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_Q_dagger }))
			}),
			// 8/9: transversal magnon
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q_down }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_Q_down_dagger, c_k }))
			})
		};
	}
	std::vector<std::unique_ptr<WickSymmetry>> Hubbard::symmetries() const
	{
		std::vector<std::unique_ptr<WickSymmetry>> ret;
		ret.reserve(2);
		//ret.push_back(std::make_unique<Index::SpinSymmetry>());
		ret.push_back(std::make_unique<InversionSymmetry>());
		ret.push_back(std::make_unique<PhaseSymmetry<SC_Type, CDW_Type>>());
		return ret;
	}
	std::string Hubbard::get_subfolder() const
	{
		return "hubbard/";
	}
}