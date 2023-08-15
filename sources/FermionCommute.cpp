#include "Term.hpp"
#include "WickTerm.hpp"
#include <fstream>
#include <sstream>

using namespace SymbolicOperators;
using term_vec = std::vector<Term>;
using op_vec = std::vector<Operator>;

int main(int argc, char** argv) {
	const Operator c_k('k', 1, false, UP, false);
	const Operator c_minus_k('k', -1, false, DOWN, false);

	const Operator c_k_dagger('k', 1, false, UP, true);
	const Operator c_minus_k_dagger('k', -1, false, DOWN, true);

	const Operator c_k_Q('k', 1, true, UP, false);
	const Operator c_minus_k_Q('k', -1, true, DOWN, false);

	const Operator c_k_Q_dagger('k', 1, true, UP, true);
	const Operator c_minus_k_Q_dagger('k', -1, true, DOWN, true);

	const Term H_T(1, Coefficient("\\epsilon_0", 'q'), std::vector<char>({ 'q' }), std::vector<std::string>({ "\\sigma" }), op_vec({
		Operator('q', 1, false, "\\sigma", true), Operator('q', 1, false, "\\sigma", false)
		}));

	const Term H_U(1, Coefficient("\\frac{U}{N}"), std::vector<char>({ 'r', 'p', 'q' }), op_vec({
		Operator('r', 1, false, UP, true), Operator('p', 1, false, DOWN, true),
		Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), DOWN, false),
		Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), UP, false),
		}));

	const Term H_V(1, Coefficient("\\tilde{V}", Momentum('q'), true), std::vector<char>({ 'r', 'p', 'q' }), std::vector<std::string>({ "\\sigma", "\\sigma'" }),
		op_vec({
			Operator('r', 1, false, "\\sigma", true),
			Operator('p', 1, false, "\\sigma'", true),
			Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), "\\sigma'", false),
			Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), "\\sigma", false),
			}));

	const term_vec H = { H_T, H_U, H_V };//

	if (argc < 2) {
		std::cerr << "Which basis?" << std::endl;
		return 1;
	}

	std::vector<term_vec> basis;
	if (std::strcmp(argv[1], "XP") == 0) {
		basis = {
			// 0: f + f^+
			term_vec({
				Term(1, { c_minus_k, c_k }),
				Term(1, { c_k_dagger, c_minus_k_dagger })
			}),
			// 1: eta + eta^+
			term_vec({
				Term(1, { c_minus_k_Q, c_k }),
				Term(1, { c_k_dagger, c_minus_k_Q_dagger })
			}),
			// 2/3: g_up/down +
			term_vec({
				Term(1, { c_k_dagger, c_k_Q }),
				Term(1, { c_k_Q_dagger, c_k })
			}),
			term_vec({
				Term(1, { c_minus_k_dagger, c_minus_k_Q }),
				Term(1, { c_minus_k_Q_dagger, c_minus_k })
			}),
			// 4/5: n_up/down
			term_vec({
				Term(1, { c_k_dagger, c_k })
			}),
			term_vec({
				Term(1, { c_minus_k_dagger, c_minus_k })
			}),
			// 6: f - f^+
			term_vec({
				Term(1, { c_minus_k, c_k }),
				Term(-1, { c_k_dagger, c_minus_k_dagger })
			}),
			// 7: eta - eta^+
			term_vec({
				Term(1, { c_minus_k_Q, c_k }),
				Term(-1, { c_k_dagger, c_minus_k_Q_dagger })
			}),
			// 8/9: g_up/down -
			term_vec({
				Term(1, { c_k_dagger, c_k_Q }),
				Term(-1, { c_k_Q_dagger, c_k })
			}),
			term_vec({
				Term(1, { c_minus_k_dagger, c_minus_k_Q }),
				Term(-1, { c_minus_k_Q_dagger, c_minus_k })
			})/**/
		};
	}
	else {
		basis = {
			// f, f^+
			term_vec({
				Term(1, { c_minus_k, c_k })
			}),
			term_vec({
				Term(1, { c_k_dagger, c_minus_k_dagger })
			}),
			// n_up/down
			term_vec({
				Term(1, { c_k_dagger, c_k })
			}),
			term_vec({
				Term(1, { c_minus_k_dagger, c_minus_k })
			}),
			// g_up/down
			term_vec({
				Term(1, { c_k_dagger, c_k_Q })
			}),
			term_vec({
				Term(1, { c_minus_k_dagger, c_minus_k_Q })
			}),
			// eta, eta^+
			term_vec({
				Term(1, { c_minus_k_Q, c_k })
			}),
			term_vec({
				Term(1, { c_k_dagger, c_minus_k_Q_dagger })
			})/**/
		};
	}

	std::vector<term_vec> basis_daggered(basis);
	for (auto& t : basis_daggered) {
		hermitianConjugate(t);
		renameMomenta(t, 'k', 'l');
	}
	for (size_t i = 0; i < basis.size(); i++)
	{
		//std::cout << toStringWithoutPrefactor(basis[i]) << " &\\, and \\, " << toStringWithoutPrefactor(basis_daggered[i]) << " \\\\" << std::endl;
		//continue;
		term_vec commute_with_H;
		commutator(commute_with_H, H, basis[i]);
		cleanUp(commute_with_H);

		//std::cout << "\\begin{align*}\n\t[ H, " << toStringWithoutPrefactor(basis[i]) << " ] ="
		//	<< commute_with_H << "\\end{align*}" << std::endl;

		for (size_t j = 0; j < basis.size(); j++)
		{
			std::cout << "\\subsection{" << j << ", " << i << "}" << std::endl;
			term_vec terms;
			commutator(terms, basis_daggered[j], commute_with_H);
			cleanUp(terms);

			//std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
			//	<< ", [H, " << toStringWithoutPrefactor(basis[i]) << " ]] =" << terms << "\\end{align*}" << std::endl;

			std::vector<WickTerm> wicks;
			for (const auto& term : terms) {
				wicks_theorem(term, wicks);
			}
			cleanWicks(wicks);
			std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
				<< ", [H, " << toStringWithoutPrefactor(basis[i]) << " ]] =" << wicks << "\\end{align*}" << std::endl;
			// serialization
			{
				// create an output file stream and a text archive to serialize the vector
				std::ofstream ofs("../commutators/wick_M_" + std::to_string(j) + "_" + std::to_string(i) + ".txt");
				boost::archive::text_oarchive oa(ofs);
				oa << wicks;
				ofs.close();
			}

			terms.clear();
			wicks.clear();
			commutator(terms, basis_daggered[j], basis[i]);
			cleanUp(terms);
			for (const auto& term : terms) {
				wicks_theorem(term, wicks);
			}
			cleanWicks(wicks);
			//std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
			//	<< ", " << toStringWithoutPrefactor(basis[i]) << " ] =" << wicks << "\\end{align*}" << std::endl;
			// serialization
			{
				// create an output file stream and a text archive to serialize the vector
				std::ofstream ofs("../commutators/wick_N_" + std::to_string(j) + "_" + std::to_string(i) + ".txt");
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
	std::vector<WickTerm> deserialized_terms;
	ia >> deserialized_terms;
	ifs.close();
	*/
	return 0;
}