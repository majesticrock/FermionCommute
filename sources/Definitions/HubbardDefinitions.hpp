#pragma once
#include "StandardOperators.hpp"
#include "../Term.hpp"
#include "../WickTerm.hpp"
#include "../WickOperatorTemplate.hpp"

namespace SymbolicOperators {
    /*struct HubbardDefinitions{
		const Momentum base_k_Q = Momentum(momentum_pairs{ {1, 'k'} }, true);
		const Momentum base_x = Momentum(momentum_pairs{ {1, 'x'} }, false);

		const Operator c_k(base_k, SpinUp, false);
		const Operator c_minus_k(-base_k, SpinDown, false);

		const Operator c_k_dagger(base_k, SpinUp, true);
		const Operator c_minus_k_dagger(-base_k, SpinDown, true);

		const Operator c_k_Q(base_k_Q, SpinUp, false);
		const Operator c_minus_k_Q(-base_k_Q, SpinDown, false);

		const Operator c_k_Q_dagger(base_k_Q, SpinUp, true);
		const Operator c_minus_k_Q_dagger(-base_k_Q, SpinDown, true);

		// transversal magnon
		const Operator c_k_Q_down_dagger(base_k_Q, SpinDown, true);
		const Operator c_k_Q_down(base_k_Q, SpinDown, false);

		const Term H_T(1, Coefficient("\\epsilon_0", 'q'), SumContainer{ MomentumSum({ 'q' }), Sigma },
			op_vec({
				Operator('q', 1, false, Sigma, true), Operator('q', 1, false, Sigma, false)
				}));

		const Term H_U(1, Coefficient("\\frac{U}{N}"), MomentumSum({ 'r', 'p', 'q' }), op_vec({
			Operator('r', 1, false, SpinUp, true), Operator('p', 1, false, SpinDown, true),
			Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), SpinDown, false),
			Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), SpinUp, false),
			}));

		const Term H_V(1, Coefficient("\\tilde{V}", Momentum('q'), true),
			SumContainer{ MomentumSum({ 'r', 'p', 'q' }), IndexSum({ Sigma, SigmaPrime }) },
			op_vec({
				Operator('r', 1, false, Sigma, true),
				Operator('p', 1, false, SigmaPrime, true),
				Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), SigmaPrime, false),
				Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), Sigma, false),
				}));

		const term_vec H = { H_T, H_U, H_V };

		const std::vector<WickOperatorTemplate> templates = {
			WickOperatorTemplate{ {IndexComparison{false, SpinDown, SpinUp}}, Momentum(), SC_Type, true },
			WickOperatorTemplate{ {IndexComparison{false, SpinDown, SpinUp}}, Momentum(momentum_pairs(), true), Eta_Type, true },
			WickOperatorTemplate{ {IndexComparison{true}}, Momentum(), Number_Type, false },
			WickOperatorTemplate{ {IndexComparison{true}}, Momentum(momentum_pairs(), true), CDW_Type, false }
		};
	};*/
}