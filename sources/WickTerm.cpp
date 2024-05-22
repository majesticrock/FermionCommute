#include "WickTerm.hpp"
#include "../../Utility/sources/RangeUtility.hpp"
#include "../../Utility/sources/Numerics/MathFunctions.hpp"
#include "KroneckerDeltaUtility.hpp"
#include "../../Utility/sources/StringUtility.hpp"
#include <variant>
#include <numeric>
#include <cctype>
#include <cassert>

#define LEFT temporary_operators[i]
#define RIGHT temporary_operators[i + 1]
#define L_SPIN temporary_operators[i].indizes[0]
#define R_SPIN temporary_operators[i + 1].indizes[0]

namespace SymbolicOperators {
	WickTerm::WickTerm(const Term* base)
		: multiplicity(base->multiplicity), coefficients(base->coefficients), sums(base->sums), operators(),
		delta_momenta(base->delta_momenta), delta_indizes(base->delta_indizes), temporary_operators()
	{
	}
	WickTerm::WickTerm(const Term& base)
		: multiplicity(base.multiplicity), coefficients(base.coefficients), sums(base.sums), operators(),
		delta_momenta(base.delta_momenta), delta_indizes(base.delta_indizes), temporary_operators()
	{
	}

	WickTerm::WickTerm(const WickTerm& base, const TemplateResult::SingleResult& result)
		: multiplicity(result.factor* base.multiplicity), coefficients(base.coefficients), sums(base.sums), operators(base.operators),
		delta_momenta(base.delta_momenta), delta_indizes(base.delta_indizes), temporary_operators()
	{
		this->operators.push_back(result.op);
		this->delta_indizes.insert(this->delta_indizes.end(), result.index_deltas.begin(), result.index_deltas.end());
	}

	void WickTerm::string_parser(std::string&& expression){
		size_t forward_pos{ expression.find(' ') };
		if(forward_pos == std::string::npos) 
			forward_pos = expression.size();
		
		const std::string sub = expression.substr(0U, forward_pos);
		const size_t sub_delimiter = sub.find(':');

		if(sub_delimiter == std::string::npos)
			throw std::invalid_argument("Did not find ':' in " + expression);

		if(sub.substr(0U, sub_delimiter) == "sum"){
			const std::string type = expression.substr(sub_delimiter + 1, expression.find('{', sub_delimiter) - sub_delimiter - 1);
			const std::vector<std::string> argument_list = Utility::extract_elements(expression);

			if(type == "index"){
				this->sums.spins.reserve(argument_list.size());
				for(const auto& arg : argument_list){
					this->sums.spins.push_back(string_to_index.at(arg));
				}
			}
			else if(type == "momentum"){
				this->sums.momenta.reserve(argument_list.size());
				for(const auto& arg : argument_list){
					assert(arg.size() == 1U);
					this->sums.momenta.push_back(arg.front());
				}
			}
			else {
				throw std::invalid_argument("Sum type not recognized " + type + " in expression " + expression);
			}
		}
		else if(sub.substr(0U, sub_delimiter) == "delta"){
			const std::string type = expression.substr(sub_delimiter + 1, expression.find('{', sub_delimiter) - sub_delimiter - 1);
			const std::vector<std::string> argument_list = Utility::extract_elements(expression);
			assert(argument_list.size() == 2U);

			if(type == "index"){
				this->delta_indizes.push_back(make_delta(string_to_index.at(argument_list[0]), string_to_index.at(argument_list[1])));
			}
			else if(type == "momentum"){
				this->delta_momenta.push_back(make_delta(Momentum(argument_list[0]), Momentum(argument_list[1])));
			} 
			else {
				throw std::invalid_argument("Delta type not recognized " + type + " in expression " + expression);
			}
		}
		else if(sub.substr(0U, sub_delimiter) == "c"){
			this->coefficients.push_back(Coefficient::parse_string(sub.substr(sub_delimiter + 1)));
		} 
		else if(sub.substr(0U, sub_delimiter) == "o"){
			this->operators.push_back(WickOperator(sub.substr(sub_delimiter + 1)));
		}
		else{
			throw std::invalid_argument("Did parse expression <" + expression + "> at <" + sub + "> with delimiter " + std::to_string(sub_delimiter));
		}
	}

	WickTerm::WickTerm(const std::string& expression) : multiplicity(1)
	{
		// Syntax
		// [factor] [index_sum] [momentum_sum] [coefficients...] [momentum_deltas...] [index_deltas...] [operators...]
		/* factor needs to be an integer
		*  Sums must be "sum:index{index1,index2,...}" or "sum:momentum{momentum_name1,momentum_name2,...}"
		*  coefficient must be "c:name{Momentum_expression1,...;index1,index2,...}"
		*  deltas must be "delta:momentum{Momentum_expression,Momentum_expression}" or "delta:index{Index,Index}"
		*  operators must be "o:type{Momentum_expression;index1,index2,...}(^+)"
		*/ 

		size_t pos{};
		if(std::isdigit(expression.front()) || expression.front() == '-' || expression.front() == '+') {
			pos = expression.find(' ');
			if(pos != std::string::npos) {
				this->multiplicity = std::stoi(expression.substr(0U, pos));
			}
		}
		size_t new_pos{};
		++pos;
		while(new_pos != std::string::npos) {
			new_pos = expression.find(' ', pos);
			string_parser(expression.substr(pos, new_pos - pos));
			pos = new_pos + 1;
		}
	}

	void wick_processor(const std::vector<Operator>& remaining, WickTermCollector& reciever_list, std::variant<WickTerm, Term> buffer)
	{
		if (remaining.empty()) {
			reciever_list.push_back(std::get<WickTerm>(buffer));
			return;
		}
		for (size_t i = 1U; i < remaining.size(); ++i)
		{
			if (std::holds_alternative<Term>(buffer)) {
				WickTerm temp(std::get<Term>(buffer));
				buffer = temp;
			}
			if ((i % 2) == 0) {
				std::get<WickTerm>(buffer).multiplicity *= -1;
			}
			std::get<WickTerm>(buffer).temporary_operators.reserve(std::get<WickTerm>(buffer).temporary_operators.size() + 2);
			std::get<WickTerm>(buffer).temporary_operators.push_back(remaining[0]);
			std::get<WickTerm>(buffer).temporary_operators.push_back(remaining[i]);

			std::vector<Operator> copy_operators = remaining;
			copy_operators.erase(copy_operators.begin() + i);
			copy_operators.erase(copy_operators.begin());
			wick_processor(copy_operators, reciever_list, buffer);

			// delete last two elements, as they are to be updated in the next iteration
			std::get<WickTerm>(buffer).temporary_operators.pop_back();
			std::get<WickTerm>(buffer).temporary_operators.pop_back();
			if ((i % 2) == 0) {
				std::get<WickTerm>(buffer).multiplicity *= -1;
			}
		}
	}

	WickTermCollector prepare_wick(const std::vector<Term>& terms)
	{
		WickTermCollector prepared_wick;
		const size_t estimated_size = std::accumulate(terms.begin(), terms.end(), size_t{}, [](size_t current, const Term& term) {
			return current + Utility::Numerics::double_factorial(term.getOperators().size());
			});

		prepared_wick.reserve(estimated_size);
		for (const auto& term : terms) {
			if (term.isIdentity()) {
				prepared_wick.push_back(WickTerm(term));
			}
			else {
				wick_processor(term.getOperators(), prepared_wick, term);
			}
		}

		return prepared_wick;
	}

	void wicks_theorem(const std::vector<Term>& terms, const std::vector<WickOperatorTemplate>& operator_templates, WickTermCollector& reciever)
	{
		WickTermCollector prepared_wick = prepare_wick(terms);

		for (auto& w_term : prepared_wick) {
			Utility::append_if(reciever, identifyWickOperators(w_term, operator_templates), [](const WickTerm& wick) {
				return !(is_always_zero(wick.delta_indizes) || is_always_zero(wick.delta_momenta));
				});
		}
	}

	void WickTerm::includeTemplateResult(const TemplateResult::SingleResult& result) {
		this->delta_indizes.insert(this->delta_indizes.begin(), result.index_deltas.begin(), result.index_deltas.end());
		this->operators.push_back(result.op);
		this->multiplicity *= result.factor;
	}

	WickTermCollector identifyWickOperators(const WickTerm& source, const std::vector<WickOperatorTemplate>& operator_templates)
	{
		WickTermCollector ret;
		ret.push_back(source);
		ret.front().temporary_operators.clear();

		for (size_t i = 0U; i < source.temporary_operators.size(); i += 2U)
		{
			std::vector<TemplateResult> template_results;
			for (const auto& operator_template : operator_templates) {
				auto template_result = operator_template.createFromOperators(source.temporary_operators[i], source.temporary_operators[i + 1U]);
				if (template_result)
					template_results.push_back(std::move(template_result));
			}

			const size_t current_size = ret.size();
			const size_t number_additional_elements = std::accumulate(template_results.begin(), template_results.end(), size_t{}, [](size_t current, const TemplateResult& tr) {
				return current + tr.results.size();
				});
			if (number_additional_elements > 1U) {
				Utility::duplicate_n_inplace(ret, number_additional_elements - 1U);
			}

			size_t template_result_it{};
			size_t old_it{};
			for (const auto& tr : template_results) {
				old_it = template_result_it;
				for (const auto& tr_result : tr.results) {
					for (auto it = ret.begin(); it != ret.begin() + current_size; ++it) {
						(it + template_result_it * current_size)->includeTemplateResult(tr_result);
					}
					++template_result_it;
				}
				std::for_each(ret.begin() + old_it * current_size, ret.begin() + current_size * template_result_it, [&tr](WickTerm& ret_element) {
					if (!tr.momentum_delta.isOne())
						ret_element.delta_momenta.push_back(tr.momentum_delta);
					});
			}
		}

		return ret;
	}

	std::ostream& operator<<(std::ostream& os, const WickTerm& term)
	{
		if (term.multiplicity > 0) {
			os << "+";
		}
		os << term.multiplicity << " \\cdot ";
		os << term.sums;
		os << term.coefficients << " ";
		for (const auto& delta : term.delta_momenta) {
			os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
		}
		for (const auto& delta : term.delta_indizes) {
			os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
		}

		if (term.isIdentity()) {
			os << " \\mathbb{1} ";
			return os;
		}
		for (const auto& op : term.operators) {
			os << op << " ";
		}
		return os;
	}
	std::ostream& operator<<(std::ostream& os, const WickTermCollector& terms)
	{
		for (WickTermCollector::const_iterator it = terms.begin(); it != terms.end(); ++it)
		{
			os << "\t&" << *it;
			if (it != terms.end() - 1) {
				os << " \\\\";
			}
			os << "\n";
		}
		return os;
	}
	WickTermCollector& operator+=(WickTermCollector& lhs, const WickTerm& rhs)
	{
		for (auto it = lhs.begin(); it != lhs.end(); ++it) {
			if (*it == rhs) {
				it->multiplicity += rhs.multiplicity;
				if (it->multiplicity == 0)
					lhs.erase(it);
				return lhs;
			}
		}
		lhs.push_back(rhs);
		return lhs;
	}
	WickTermCollector& operator-=(WickTermCollector& lhs, const WickTerm& rhs)
	{
		for (auto it = lhs.begin(); it != lhs.end(); ++it) {
			if (*it == rhs) {
				it->multiplicity -= rhs.multiplicity;
				if (it->multiplicity == 0)
					lhs.erase(it);
				return lhs;
			}
		}
		lhs.push_back(rhs);
		return lhs;
	}

	WickTermCollector& operator+=(WickTermCollector& lhs, const WickTermCollector& rhs)
	{
		for (const auto& term : rhs) {
			lhs += term;
		}
		return lhs;
	}
	WickTermCollector& operator-=(WickTermCollector& lhs, const WickTermCollector& rhs)
	{
		for (const auto& term : rhs) {
			lhs -= term;
		}
		return lhs;
	}
}