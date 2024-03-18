#include "WickOperatorTemplate.hpp"
#include <algorithm>
#include <iterator>
#include "Momentum.hpp"
#include "KroneckerDelta.hpp"

namespace SymbolicOperators {
	TemplateResult WickOperatorTemplate::_handle_sc_type(const Operator& left, const Operator& right) const {
		// c_{-k-q} c_{k} or c_{k}^+ c_{-k-q}^+
		const Operator& base{ left.isDaggered ? left : right };
		const Operator& other{ left.isDaggered ? right : left };
		// q
		const Momentum momentum_diff{ -(base.momentum + other.momentum) };

		TemplateResult result(1U, this->type, base.momentum);
 		result.results.front().op.isDaggered = left.isDaggered;
		result.momentum_delta = make_delta(this->momentum_difference, momentum_diff);

		for (size_t i = 0U; i < indexComparison.size(); ++i)
		{
			if (indexComparison[i].any_identical) {
				result.add_index_delta(make_delta(base.indizes[i], other.indizes[i]));
			}
			else {
				const size_t previous_size{ result.create_branch() };
				result.add_index_delta_range(make_delta(base.indizes[i], indexComparison[i].base), 0U, previous_size);
				result.add_index_delta_range(make_delta(other.indizes[i], indexComparison[i].other), 0U, previous_size);

				// c^+ c^+ can be swapped for the cost of a sign
				result.add_index_delta_range(make_delta(base.indizes[i], indexComparison[i].other), previous_size, previous_size);
				result.add_index_delta_range(make_delta(other.indizes[i], indexComparison[i].base), previous_size, previous_size);
				result.operation_on_range([](TemplateResult::SingleResult& res) { res.factor *= -1; }, previous_size, previous_size);
				result.operation_on_range([&other](TemplateResult::SingleResult& res) { res.op.momentum = other.momentum; }, previous_size, previous_size);
			}
		}

		return result;
	}
	TemplateResult WickOperatorTemplate::_handle_num_type(const Operator& left, const Operator& right) const {
		// c_{k}^+ c_{k+q}
		// q
		const Momentum momentum_diff = right.momentum - left.momentum;
		KroneckerDelta<Momentum> momentum_delta{ this->momentum_difference, momentum_diff };

		std::vector<KroneckerDelta<Index>> index_delta;
		
		TemplateResult result(1U, this->type, left.momentum);
 		result.results.front().op.isDaggered = false;
		result.momentum_delta = make_delta(this->momentum_difference, momentum_diff);

	for (size_t i = 0U; i < indexComparison.size(); ++i)
		{
			if (indexComparison[i].any_identical) {
				result.add_index_delta(make_delta(left.indizes[i], right.indizes[i]));
			}
			else {
				const size_t previous_size{ result.create_branch() };
				result.add_index_delta_range(make_delta(left.indizes[i], indexComparison[i].base), 0U, previous_size);
				result.add_index_delta_range(make_delta(right.indizes[i], indexComparison[i].other), 0U, previous_size);
			}
		}

		return result;
	}

	TemplateResult WickOperatorTemplate::createFromOperators(const Operator& left, const Operator& right) const {
		if (this->is_sc_type) {
			if (left.isDaggered != right.isDaggered)
				return {};
			return this->_handle_sc_type(left, right);
		}

		if (left.isDaggered == right.isDaggered)
			return {};
		// The input needs to be normal ordered.
		// This means that if right is not daggered, left cannot be daggered either
		assert(! left.isDaggered);
		return this->_handle_num_type(left, right);
	}
}