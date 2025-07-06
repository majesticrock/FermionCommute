#include "LatticeCUT.hpp"

namespace mrock::symbolic_operators {
    std::vector<Term> LatticeCUT::hamiltonian() const
    {
        const Term H_Kin(1, Coefficient("\\epsilon_0", Momentum('q')), SumContainer{ MomentumSum({ 'q' }), Index::Sigma },
			std::vector<Operator>({
				Operator('q', 1, false, Index::Sigma, true), Operator('q', 1, false, Index::Sigma, false)
				})
            );

        const Term H_Ph(-1, Coefficient::RealInversionSymmetric("g", MomentumList({ 'q', 'p' }), std::function<void(Coefficient&)>([](Coefficient& coeff){ coeff.momenta.sort(); })),
			SumContainer{ MomentumSum({ 'p', 'q' }) },
			std::vector<Operator>({
				c_k_dagger.with_momentum('q'), c_minus_k_dagger.with_momentum('q'),
				c_minus_k.with_momentum('p'), c_k.with_momentum('p') })
			);

        const Term H_U(1, Coefficient("U"), MomentumSum({ 'r', 'p', 'q' }), 
            std::vector<Operator>({
			    Operator('r', 1, false, Index::SpinUp, true), Operator('p', 1, false, Index::SpinDown, true),
			    Operator(momentum_symbols({ MomentumSymbol(1, 'p'), MomentumSymbol(-1, 'q') }), Index::SpinDown, false),
			    Operator(momentum_symbols({ MomentumSymbol(1, 'r'), MomentumSymbol(1, 'q') }), Index::SpinUp, false),
			}));

        return { H_Kin, H_Ph, H_U };
    }

    std::vector<WickOperatorTemplate> LatticeCUT::templates() const
    {
        return {
			WickOperatorTemplate{ {IndexComparison{false, Index::SpinDown, Index::SpinUp}}, Momentum(), SC_Type, true },
			WickOperatorTemplate{ {IndexComparison{true}}, Momentum(), Number_Type, false }
		};
    }

    std::vector<std::vector<Term>> LatticeCUT::XP_basis() const
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

    std::vector<std::vector<Term>> LatticeCUT::STD_basis() const
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

    std::vector<std::unique_ptr<WickSymmetry>> LatticeCUT::symmetries() const
    {
        std::vector<std::unique_ptr<WickSymmetry>> ret;
		ret.push_back(std::make_unique<SpinSymmetry>());
		ret.push_back(std::make_unique<InversionSymmetry>());
		ret.push_back(std::make_unique<PhaseSymmetry<SC_Type>>());
		return ret;
    }

    std::string LatticeCUT::get_subfolder() const
    {
        return "lattice_cut/";
    }
}