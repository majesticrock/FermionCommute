#include "Term.hpp"
#include <fstream>
#include <sstream>
using namespace SymbolicOperators;

int main(int argc, char** argv) {
	Operator c_k('k', 1, false, UP, false);
	Operator c_minus_k('k', -1, false, DOWN, false);

	Operator c_k_dagger('k', 1, false, UP, true);
	Operator c_minus_k_dagger('k', -1, false, DOWN, true);

	Operator c_k_Q('k', 1, true, UP, false);
	Operator c_minus_k_Q('k', -1, true, DOWN, false);

	Operator c_k_Q_dagger('k', 1, true, UP, true);
	Operator c_minus_k_Q_dagger('k', -1, true, DOWN, true);

	Term H_T(1, Coefficient("\\epsilon_0", 'q'), std::vector<char>({ 'q' }), std::vector<std::string>({ "\\sigma" }), std::vector<Operator>({
		Operator('q', 1, false, "\\sigma", true), Operator('q', 1, false, "\\sigma", false)
		}));

	Term H_U(1, Coefficient("\\frac{U}{N}"), std::vector<char>({ 'r', 'p', 'q' }), std::vector<Operator>({
		Operator('r', 1, false, UP, true), Operator('p', 1, false, DOWN, true),
		Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), DOWN, false),
		Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), UP, false),
		}));

	Term H_V(1, Coefficient("\\tilde{V}", Momentum('q'), true), std::vector<char>({'r', 'p', 'q'}), std::vector<std::string>({"\\sigma", "\\sigma'"}),
		std::vector<Operator>({
			Operator('r', 1, false, "\\sigma", true),
			Operator('p', 1, false, "\\sigma'", true),
			Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), "\\sigma'", false),
			Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), "\\sigma", false),
		}));

	std::vector<Term> H = { H_V };

	std::vector<Term> basis = {
		// f, f^+
		Term(1, Coefficient(), std::vector<Operator>({ c_minus_k, c_k }))
		//Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_minus_k_dagger })),
		//// g_up/down
		//Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_k_Q })),
		//Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_dagger, c_minus_k_Q })),
		//// n_up/down
		//Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_k })),
		//Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_dagger, c_minus_k })),
		//// eta, eta^+
		//Term(1, Coefficient(), std::vector<Operator>({ c_minus_k_Q, c_k })),
		//Term(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_minus_k_Q_dagger }))
	};
	std::vector<Term> basis_daggered(basis);
	for (auto& t : basis_daggered) {
		t.hermitianConjugate();
		t.renameMomenta('k', 'l');
	}
	//{
	//	Term test(1, Coefficient(), std::vector<Operator>({ c_k_dagger, c_minus_k_dagger, c_minus_k, c_k }));
	//	std::vector<WickTerm> wicks;
	//	wicks_theorem(test, wicks);
	//	cleanWicks(wicks);
	//	std::cout << "\\begin{align*}\n\t[ " << ", [H, " << test << " ]] =" << wicks << "\\end{align*}" << std::endl;
	//	
	//}
	//return 0;
	for (size_t i = 0; i < basis.size(); i++)
	{
		std::vector<Term> commute_with_H;
		commutator(commute_with_H, H, basis[i]);
		cleanUp(commute_with_H);
		//std::cout << "H: \\begin{align*}\n\t" << commute_with_H << "\\end{align*}" << std::endl;
		for (size_t j = 0; j < basis.size(); j++)
		{
			//if (i != 2) continue;
			//if (j != 0) continue;
			std::cout << "\\subsection{" << j << ", " << i << "}" << std::endl;
			std::vector<Term> terms;
			commutator(terms, basis_daggered[j], commute_with_H);
			cleanUp(terms);
			//std::cout << "full: \\begin{align*}\n\t" << terms << "\\end{align*}" << std::endl;
			std::vector<WickTerm> wicks;
			for (const auto& term : terms) {
				wicks_theorem(term, wicks);
			}
			cleanWicks(wicks);
			std::cout << "\\begin{align*}\n\t[ " << basis_daggered[j].toStringWithoutPrefactor() << ", [H, " << basis[i].toStringWithoutPrefactor() << " ]] =" << wicks << "\\end{align*}" << std::endl;
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
			std::cout << "\\begin{align*}\n\t[ " << basis_daggered[j].toStringWithoutPrefactor() << ", " << basis[i].toStringWithoutPrefactor() << " ] =" << wicks << "\\end{align*}" << std::endl;
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
	Term right(1, Coefficient(), std::vector<Operator>({
		c_minus_k, c_k
		}));
	Term left(1, Coefficient(), std::vector<Operator>({
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