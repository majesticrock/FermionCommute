#include <fstream>
#include <sstream>
#include <mrock/utility/LaTeXOutput.hpp>
#include "Definitions/Hubbard.hpp"
#include "Definitions/HubbardDispersions.hpp"
#include "Definitions/Continuum.hpp"
#include <mrock/symbolic_operators/Wick.hpp>
#include <memory>
#include <filesystem>
#include <boost/archive/binary_oarchive.hpp>

using namespace mrock::symbolic_operators;
using term_vec = std::vector<Term>;
using op_vec = std::vector<Operator>;

void remove_all_x(term_vec& terms) {
	for (auto& term : terms) {
		term.remove_momentum_contribution('x');
	}
}
void remove_all_x(WickTermCollector& terms) {
	for (auto& term : terms) {
		term.remove_momentum_contribution('x');
	}
}

template <typename T>
std::vector<T> joinVectors(const std::vector<T>& vec1, const std::vector<T>& vec2) {
	std::vector<T> result;
	result.reserve(vec1.size() + vec2.size());
	result.insert(result.end(), vec1.begin(), vec1.end());
	result.insert(result.end(), vec2.begin(), vec2.end());
	return result;
}


std::unique_ptr<DefinitionsBase> get_model(std::string const& model_type) {
	if (model_type == "hubbard") {
		return std::make_unique<Hubbard>();
	}
	else if (model_type == "continuum") {
		return std::make_unique<Continuum>();
	}
	else if (model_type == "hubbard_dispersions") {
		return std::make_unique<HubbardDispersions>();
	}
	else {
		throw std::invalid_argument("Model not recognized! " + model_type);
	}
}

int main(int argc, char** argv) {
	const std::string save_folder = "../commutators/";
	/* WickTerm parse_test("1 sum:momentum{p,q} c:V{p;} o:n{k-p-3x;up} o:f{k+l;}");
	std::cout << parse_test << "    " << parse_test.coefficients.size() << std::endl;
	return 0; */
	constexpr bool print = false;
	if (argc < 3) {
		std::cerr << "Syntax: ./build/main <XP/std> <model>" << std::endl;
		return 1;
	}
	const std::string EXECUTION_TYPE = argv[1];
	const std::string MODEL_TYPE = argv[2];

	if (EXECUTION_TYPE == "test") {
		Hubbard hubbard;

		std::vector<std::vector<Term>> base = hubbard.STD_basis();
		std::vector<std::vector<Term>> disp = base;
		for (auto& _v : disp) {
			for (auto& v : _v) {
				if (v.operators.front().is_daggered) {
					v.operators.front().momentum += Momentum('x');
				}
				else {
					v.operators.front().momentum += Momentum('x', -1);
				}
			}
		}

		std::vector<std::vector<Term>> base_daggered(base);
		std::vector<std::vector<Term>> disp_daggered(disp);
		for (auto& vec : base_daggered) {
			hermitianConjugate(vec);
			rename_momenta(vec, 'k', 'l');
		}
		for (auto& vec : disp_daggered) {
			hermitianConjugate(vec);
			rename_momenta(vec, 'k', 'l');
		}


		const Term H_U(1, Coefficient("\\frac{U}{N}"), MomentumSum({ 'r', 'p', 'q' }), std::vector<Operator>({
			Operator('r', 1, false, Index::SpinUp, true), Operator('p', 1, false, Index::SpinDown, true),
			Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), Index::SpinDown, false),
			Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), Index::SpinUp, false),
			}));
		const std::vector<Term> H = { H_U };
		const int inner_idx = 0;
		const int outer_idx = static_cast<int>(!inner_idx);

		term_vec commute_with_H_base;
		commutator(commute_with_H_base, H, base[inner_idx]);
		cleanUp(commute_with_H_base);

		term_vec commute_with_H_disp;
		commutator(commute_with_H_disp, disp[inner_idx], H);
		cleanUp(commute_with_H_disp);
		
		if (true) {
			term_vec joined = joinVectors(commute_with_H_base, commute_with_H_disp);
			remove_all_x(joined);
			cleanUp(joined);

			std::cout << "Single commutator:\n" << joined << std::endl; // Up to here, everything works
		}

		{
			term_vec base_double;
			commutator(base_double, base_daggered[outer_idx], commute_with_H_base);
			cleanUp(base_double);

			term_vec disp_double;
			commutator(disp_double, disp_daggered[outer_idx], commute_with_H_disp);
			cleanUp(disp_double);
			
			term_vec joined = joinVectors(base_double, disp_double);
			cleanUp(joined);
			//std::cout << "joined:\n" << joined << std::endl;

			auto templates = hubbard.templates();
			auto symmetries = hubbard.symmetries();

			WickTermCollector wicks;
			wicks_theorem(joined, templates, wicks);
			clearEtas(wicks);
			cleanWicks(wicks, symmetries);
			remove_all_x(wicks);
			cleanWicks(wicks, symmetries);

			std::cout << "Double commutator:\n" << wicks << std::endl;
		}
		return 0;
	}

	const std::string name_prefix = EXECUTION_TYPE == "XP" ? "XP_" : "";
	const bool debug = EXECUTION_TYPE == "debug";

	const std::unique_ptr<DefinitionsBase> model = get_model(MODEL_TYPE);
	const std::string sub_folder = model->get_subfolder();
	if (!debug)
		std::filesystem::create_directories(save_folder + sub_folder);

	const term_vec H = model->hamiltonian();
	const std::vector<WickOperatorTemplate> templates = model->templates();

	std::vector<term_vec> basis;
	if (EXECUTION_TYPE == "XP") {
		basis = model->XP_basis();
	}
	else if (EXECUTION_TYPE == "std") {
		basis = model->STD_basis();
	}
	else if (debug) {
		basis = model->STD_basis();
	}
	else {
		std::cerr << "Execution type not recognized! Accepted are: 'XP' and 'std'" << std::endl;
		return 1;
	}
	const auto symmetries = model->symmetries();

	std::vector<term_vec> basis_daggered(basis);
	for (auto& t : basis_daggered) {
		hermitianConjugate(t);
		rename_momenta(t, 'k', 'l');
		//if (debug) {
		//	rename_momenta(t, 'x', 'y');
		//}
	}

	if (print) std::cout << "\\begin{align*}\n\t H =" << H << "\\end{align*}" << std::endl;

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
			if (print) {
				std::cout << "\\subsection{" << i << "." << j << "}" << std::endl;
			}
			term_vec terms;
			commutator(terms, basis_daggered[j], commute_with_H);
			cleanUp(terms);

			if (debug)
				std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
				<< ", [H, " << toStringWithoutPrefactor(basis[i]) << " ]] =" << terms << "\\end{align*}" << std::endl;

			WickTermCollector wicks;
			wicks_theorem(terms, templates, wicks);
			clearEtas(wicks);
			cleanWicks(wicks, symmetries);
			
			for (auto& wickterm : wicks) {
				if (wickterm.coefficients.front().name == "\\rho") {
					wickterm.sums.push_back(Index::SigmaPrime);
					wickterm.sums.push_back('r');
					wickterm.coefficients.front() = Coefficient::parse_string("V{0;}");
					wickterm.operators.push_back(WickOperator("n{r;sigma'}"));
				}
			}
			cleanWicks(wicks, symmetries);

			if (debug || print) {
				std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
					<< ", [H, " << toStringWithoutPrefactor(basis[i]) << " ]] =" << wicks << "\\end{align*}" << std::endl;
			}

			// serialization
			if (!debug) {
				// create an output file stream and a text archive to serialize the vector
				std::ofstream ofs(save_folder + sub_folder + name_prefix + "wick_M_" + std::to_string(j) + "_" + std::to_string(i) + ".bin", std::ios::binary);
				boost::archive::binary_oarchive oa(ofs);
				oa << wicks;
				ofs.close();
			}

			terms.clear();
			wicks.clear();
			commutator(terms, basis_daggered[j], basis[i]);
			cleanUp(terms);
			wicks_theorem(terms, templates, wicks);
			clearEtas(wicks);
			cleanWicks(wicks, symmetries);

			if (debug || print)
				std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
				<< ", " << toStringWithoutPrefactor(basis[i]) << " ] =" << wicks << "\\end{align*}" << std::endl;
			// serialization
			if (!debug) {
				// create an output file stream and a text archive to serialize the vector
				std::ofstream ofs(save_folder + sub_folder + name_prefix + "wick_N_" + std::to_string(j) + "_" + std::to_string(i) + ".bin", std::ios::binary);
				boost::archive::binary_oarchive oa(ofs);
				oa << wicks;
				ofs.close();
			}
		}
	}

	return 0;
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