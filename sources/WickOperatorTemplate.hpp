#pragma once

#include "Operator.hpp"
#include "WickOperator.hpp"
#include "KroneckerDelta.hpp"
#include <optional>
#include <algorithm>

namespace SymbolicOperators {
	struct IndexComparison {
		bool any_identical;
		Index base{ UndefinedIndex };
		Index other{ UndefinedIndex };
	};

	struct TemplateResult {
		struct SingleResult {
			int factor{};
			WickOperator op;
			std::vector<KroneckerDelta<Index>> index_deltas;

			// TODO: std::remove_if
			void clear_delta_equals_one(){
				for(auto it = this->index_deltas.begin(); it != this->index_deltas.end();){
					if(it->first == it->second){
						it = this->index_deltas.erase(it);
					} 
					else {
						++it;
					}
				}
			};
			// TODO:: std::any_of
			bool contains_impossible_delta() const {
				for(const auto& delta : this->index_deltas){
					if(!is_mutable(delta.first) && !is_mutable(delta.second) && delta.first != delta.second){
						return true;
					}
				}
				return false;
			};
		};
		std::vector<SingleResult> results;
		KroneckerDelta<Momentum> momentum_delta;

		TemplateResult() = default;
		TemplateResult(size_t initial_size, OperatorType operator_type, const Momentum& base_momentum)
			: results(initial_size) {
			for (auto& result : results)
			{
				result.op.type = operator_type;
				result.op.momentum = base_momentum;
				result.factor = 1;
			}
		}

		template<class Operation>
		void operation_on_range(const Operation& operation, size_t begin, size_t n) {
			for (size_t i = begin; i < begin + n; ++i)
			{
				operation(results[i]);
			}
		}
		template<class Operation>
		void operation_on_each(const Operation& operation) {
			for (auto& res : results)
			{
				operation(res);
			}
		}
		inline void add_index_delta_range(const KroneckerDelta<Index>& index, size_t begin, size_t n) {
			operation_on_range([&index](SingleResult& res) { res.index_deltas.push_back(index); }, begin, n);
		}
		inline void add_index_delta(const KroneckerDelta<Index>& index) {
			operation_on_each([&index](SingleResult& res) { res.index_deltas.push_back(index); });
		}
		inline size_t create_branch() {
			const size_t current_size{ results.size() };
			results.reserve(2 * current_size);
			std::copy_n(results.begin(), current_size, std::back_inserter(results));
			return current_size;
		}

		void clear_impossible(){
			for(auto it = this->results.begin(); it != this->results.end();){
				if(it->contains_impossible_delta()){
					it = this->results.erase(it);
				}
				else{
					++it;
				}
			}
		};
	};

	struct WickOperatorTemplate {
		std::vector<IndexComparison> indexComparison;
		Momentum momentum_difference;
		OperatorType type;
		bool is_sc_type{};

		// Returns the corresponding WickOperator if construction is possible
		// Otherwise, it returns an empty optional
		TemplateResult createFromOperators(const Operator& left, const Operator& right) const;

	private:
		TemplateResult _handle_sc_type(const Operator& left, const Operator& right) const;
		TemplateResult _handle_num_type(const Operator& left, const Operator& right) const;
	};
}