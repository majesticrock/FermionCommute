#include "Term.hpp"
#include "WickTerm.hpp"
#include <fstream>
#include <sstream>

using namespace SymbolicOperators;
using term_vec = std::vector<Term>;
using op_vec = std::vector<Operator>;

int main(int argc, char** argv) {
	const Operator c_k('k', 1, false, SpinUp, false);
	const Operator c_minus_k('k', -1, false, SpinDown, false);

	const Operator c_k_dagger('k', 1, false, SpinUp, true);
	const Operator c_minus_k_dagger('k', -1, false, SpinDown, true);

	const Operator c_k_Q('k', 1, true, SpinUp, false);
	const Operator c_minus_k_Q('k', -1, true, SpinDown, false);

	const Operator c_k_Q_dagger('k', 1, true, SpinUp, true);
	const Operator c_minus_k_Q_dagger('k', -1, true, SpinDown, true);

	// transversal magnon
	const Operator c_k_Q_down_dagger('k', 1, true, SpinDown, true);
	const Operator c_k_Q_down('k', 1, true, SpinDown, false);

	const Term H_T(1, Coefficient("\\epsilon_0", 'q'), MomentumSum({ 'q' }), Sigma, op_vec({
		Operator('q', 1, false, Sigma, true), Operator('q', 1, false, Sigma, false)
		}));

	const Term H_U(1, Coefficient("\\frac{U}{N}"), MomentumSum({ 'r', 'p', 'q' }), op_vec({
		Operator('r', 1, false, SpinUp, true), Operator('p', 1, false, SpinDown, true),
		Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), SpinDown, false),
		Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), SpinUp, false),
		}));

	const Term H_V(1, Coefficient("\\tilde{V}", Momentum('q'), true), MomentumSum({ 'r', 'p', 'q' }), IndexSum({ Sigma, SigmaPrime }),
		op_vec({
			Operator('r', 1, false, Sigma, true),
			Operator('p', 1, false, SigmaPrime, true),
			Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), SigmaPrime, false),
			Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), Sigma, false),
			}));

	const term_vec H = { H_T, H_U, H_V };//

	if (argc < 2) {
		std::cerr << "Which basis?" << std::endl;
		return 1;
	}

	const std::string name_prefix = std::strcmp(argv[1], "XP") == 0 ? "XP_" : "";
	std::vector<term_vec> basis;
	if (std::strcmp(argv[1], "XP") == 0) {
		basis = {
			// 0: f + f^+
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k, c_k })),
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
			}),
			// 1: eta + eta^+
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_Q, c_k })),
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_minus_k_Q_dagger }))
			}),
			// 2/3: g_up/down +
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_k_Q })),
				Term(1, Coefficient(), std::vector<Operator>({ c_k_Q_dagger, c_k }))
			}),
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_dagger, c_minus_k_Q })),
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_Q_dagger, c_minus_k }))
			}),
			// 4: transversal magnon, hermitian
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_k_Q_down })),
				Term(1, Coefficient(), std::vector<Operator>({ c_k_Q_down_dagger, c_k }))
			}),
			// 5/6: n_up/down
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_k }))
			}),
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_dagger, c_minus_k }))
			}),
			// 7: f - f^+
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k, c_k })),
				Term(-1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
			}),
			// 8: eta - eta^+
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_Q, c_k })),
				Term(-1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_minus_k_Q_dagger }))
			}),
			// 9/10: g_up/down -
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_k_Q })),
				Term(-1, Coefficient(), std::vector<Operator>({ c_k_Q_dagger, c_k }))
			}),
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_dagger, c_minus_k_Q })),
				Term(-1, Coefficient(), std::vector<Operator>({ c_minus_k_Q_dagger, c_minus_k }))
			}),
			// 11: transversal magnon, antihermitian
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_k_Q_down })),
				Term(-1, Coefficient(), std::vector<Operator>({ c_k_Q_down_dagger, c_k }))
			})
		};
	}
	else {
		basis = {
			// f, f^+
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k, c_k }))
			}),
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
			}),
			// n_up/down
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_k }))
			}),
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_dagger, c_minus_k }))
			}),
			// g_up/down
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_k_Q }))
			}),
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_dagger, c_minus_k_Q }))
			}),
			// eta, eta^+
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_Q, c_k }))
			}),
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_minus_k_Q_dagger }))
			}),
			// transversal magnon
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_k_Q_down }))
			}),
			std::vector<Term>({
				Term(1, Coefficient(), std::vector<Operator>({ c_k_Q_down_dagger, c_k }))
			}),
		};
	}

	std::vector<term_vec> basis_daggered(basis);
	for (auto& t : basis_daggered) {
		hermitianConjugate(t);
		renameMomenta(t, 'k', 'l');
	}
	for (size_t i = 0U; i < basis.size(); ++i)
	{
		//std::cout << toStringWithoutPrefactor(basis[i]) << " &\\, and \\, " << toStringWithoutPrefactor(basis_daggered[i]) << " \\\\" << std::endl;
		//continue;
		term_vec commute_with_H;
		commutator(commute_with_H, H, basis[i]);
		cleanUp(commute_with_H);

		//std::cout << "\\begin{align*}\n\t[ H, " << toStringWithoutPrefactor(basis[i]) << " ] ="
		//	<< commute_with_H << "\\end{align*}" << std::endl;

		for (size_t j = 0U; j < basis.size(); ++j)
		{
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
			//if (i > 9) {
			//	std::cout << "\\subsection{" << j << ", " << i << "}" << std::endl;
			//	std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
			//		<< ", [H, " << toStringWithoutPrefactor(basis[i]) << " ]] =" << wicks << "\\end{align*}" << std::endl;
			//}
			// serialization
			{
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
			for (const auto& term : terms) {
				wicks_theorem(term, wicks);
			}
			cleanWicks(wicks);
			//if (i > 9)
			//	std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
			//	<< ", " << toStringWithoutPrefactor(basis[i]) << " ] =" << wicks << "\\end{align*}" << std::endl;
			// serialization
			{
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
	std::vector<WickTerm> deserialized_terms;
	ia >> deserialized_terms;
	ifs.close();
	*/
	return 0;
}