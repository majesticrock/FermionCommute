#include "DefinitionsBase.hpp"

namespace mrock::symbolic_operators {
	const Momentum base_k = Momentum('k');
	const Momentum base_x = Momentum('x');
	const Momentum base_k_Q{ momentum_pairs{ {1, 'k'} }, true };

	const Operator c_k = Operator{ base_k, Index::SpinUp, false };
	const Operator c_minus_k = Operator{ -base_k, Index::SpinDown, false };
	const Operator c_k_dagger = Operator{ base_k, Index::SpinUp, true };
	const Operator c_minus_k_dagger = Operator{ -base_k, Index::SpinDown, true };

	const Operator c_k_Q{ base_k_Q, Index::SpinUp, false };
	const Operator c_minus_k_Q{ -base_k_Q, Index::SpinDown, false };
	const Operator c_k_Q_dagger{ base_k_Q, Index::SpinUp, true };
	const Operator c_minus_k_Q_dagger{ -base_k_Q, Index::SpinDown, true };
	const Operator c_k_Q_down_dagger{ base_k_Q, Index::SpinDown, true };
	const Operator c_k_Q_down{ base_k_Q, Index::SpinDown, false };

	const Operator c_k_down_dagger{ base_k, Index::SpinDown, true };
	const Operator c_k_down{ base_k, Index::SpinDown, false };
}