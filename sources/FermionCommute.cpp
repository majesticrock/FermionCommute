#include <fstream>
#include <sstream>
#include "../../Utility/sources/LaTeXOutput.hpp"
#include "Definitions/Hubbard.hpp"
#include "WickTerm.hpp"
#include <memory>

using namespace SymbolicOperators;
using term_vec = std::vector<Term>;
using op_vec = std::vector<Operator>;

int main(int argc, char** argv) {
	if (argc < 2) {
		std::cerr << "Which basis?" << std::endl;
		return 1;
	}
	const std::string EXECUTION_TYPE = argv[1];

	const std::unique_ptr<StandardOperators> model = std::make_unique<Hubbard>();

	const term_vec H = model->hamiltonian();
	const std::vector<WickOperatorTemplate> templates = model->templates();

	if (EXECUTION_TYPE == "test") {
		WickTerm wick;
		wick.multiplicity = 1;
		wick.temporary_operators = { Hubbard::c_minus_k_Q, Hubbard::c_k_Q,
			Hubbard::c_k_Q_dagger, StandardOperators::c_k };
		auto wick_results = identifyWickOperators(wick, templates);
		std::cout << "Testing on: $" << wick.temporary_operators << "$\n\n";
		std::cout << "Pre clean:\n\n" << Utility::as_LaTeX(wick_results, "align*") << std::endl;
		cleanWicks(wick_results);
		std::cout << "Post clean:\n\n" << Utility::as_LaTeX(wick_results, "align*") << std::endl;

		return 0;
	}

	const std::string name_prefix = EXECUTION_TYPE == "XP" ? "XP_" : "";
	const bool debug = EXECUTION_TYPE == "debug";
	std::vector<term_vec> basis;
	if (EXECUTION_TYPE == "XP") {
		basis = model->XP_basis();
	}
	else if (EXECUTION_TYPE == "STD") {
		basis = model->STD_basis();
	}
	else if (debug) {
		basis = {
			// 0: f + f^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ Operator(Momentum({{-1, 'k'}, {-1, 'x'}}), SpinDown, false), StandardOperators::c_k}))
			})
		};
	}

	std::vector<term_vec> basis_daggered(basis);
	for (auto& t : basis_daggered) {
		hermitianConjugate(t);
		renameMomenta(t, 'k', 'l');
		if (debug) {
			renameMomenta(t, 'x', 'y');
		}
	}
	for (size_t i = 0U; i < basis.size(); ++i)
	{
		term_vec commute_with_H;
		commutator(commute_with_H, H, basis[i]);
		cleanUp(commute_with_H);
		if (debug)
			std::cout << "\\begin{align*}\n\t[ H, " << toStringWithoutPrefactor(basis[i]) << " ] ="
			<< commute_with_H << "\\end{align*}" << std::endl;

		for (size_t j = 0U; j < basis.size(); ++j)
		{
			term_vec terms;
			commutator(terms, basis_daggered[j], commute_with_H);
			cleanUp(terms);

			if (debug)
				std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
				<< ", [H, " << toStringWithoutPrefactor(basis[i]) << " ]] =" << terms << "\\end{align*}" << std::endl;

			WickTermCollector wicks;
			wicks_theorem(terms, templates, wicks);
			cleanWicks(wicks);

			if (debug) {
				std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
					<< ", [H, " << toStringWithoutPrefactor(basis[i]) << " ]] =" << wicks << "\\end{align*}" << std::endl;
			}

			// serialization
			if (!debug) {
				// create an output file stream and a text archive to serialize the vector
				std::ofstream ofs("../commutators/" + name_prefix + "wick_M_" + std::to_string(j) + "_" + std::to_string(i) + ".txt");
				boost::archive::text_oarchive oa(ofs);
				oa << wicks;
				ofs.close();
			}

			terms.clear();
			wicks.clear();
			commutator(terms, basis_daggered[j], basis[i]);
			cleanUp(terms);
			wicks_theorem(terms, templates, wicks);
			cleanWicks(wicks);

			if (debug)
				std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
				<< ", " << toStringWithoutPrefactor(basis[i]) << " ] =" << wicks << "\\end{align*}" << std::endl;
			// serialization
			if (!debug) {
				// create an output file stream and a text archive to serialize the vector
				std::ofstream ofs("../commutators/" + name_prefix + "wick_N_" + std::to_string(j) + "_" + std::to_string(i) + ".txt");
				boost::archive::text_oarchive oa(ofs);
				oa << wicks;
				ofs.close();
			}
		}
	}

	/* Output code if needed
	*
	Term right(1, Coefficient(), op_vec({
		c_minus_k, c_k
		}));
	Term left(1, Coefficient(), op_vec({
		c_l_dagger, c_minus_l_dagger
		}));

	std::cout << "\\begin{align*}\n\t\\langle [" << left.toStringWithoutPrefactor() << ", " << right.toStringWithoutPrefactor() << "] \\rangle = " << wicks << "\\end{align*}" << std::endl;
	std::cout << "\\begin{align*}\n\t\\langle [" << left.toStringWithoutPrefactor() << ", [ H, " << right.toStringWithoutPrefactor() << "] ] \\rangle = " << wicks << "\\end{align*}" << std::endl;
	*/

	/* Example of how to to read the output
	// create an input file stream and a text archive to deserialize the vector
	std::ifstream ifs("wick_terms.txt");
	boost::archive::text_iarchive ia(ifs);
	WickTermCollector deserialized_terms;
	ia >> deserialized_terms;
	ifs.close();
	*/
	return 0;
}