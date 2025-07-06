#include <fstream>
#include <sstream>
#include <mrock/utility/LaTeXOutput.hpp>
#include "Definitions/Hubbard.hpp"
#include "Definitions/HubbardDispersions.hpp"
#include "Definitions/Continuum.hpp"
#include "Definitions/LatticeCUT.hpp"
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

const std::string hubbard_type = "hubbard";
const std::string continuum_type = "continuum";
const std::string hubbard_dispersions_type = "hubbard_dispersions";
const std::string lattice_cut_type = "lattice_cut";

std::unique_ptr<DefinitionsBase> get_model(std::string const& model_type) {
	if (model_type == hubbard_type) {
		return std::make_unique<Hubbard>();
	}
	else if (model_type == continuum_type) {
		return std::make_unique<Continuum>();
	}
	else if (model_type == hubbard_dispersions_type) {
		return std::make_unique<HubbardDispersions>();
	}
	else if (model_type == lattice_cut_type) {
		return std::make_unique<LatticeCUT>();
	}
	else {
		throw std::invalid_argument("Model not recognized! " + model_type);
	}
}

void boson_test(){
	Term H_diag(1, Coefficient("\\lambda"), SumContainer{ MomentumSum({'k'}), Index::Sigma},
        std::vector<Operator>({
            Operator::Boson(Momentum('k'), Index::Sigma, true),
            Operator::Boson(Momentum('k'), Index::Sigma, false)
        }));

	Term H_bogo(1, Coefficient::HoneyComb("\\gamma", Momentum('k'), true), SumContainer{ MomentumSum({'k'}) },
		std::vector<Operator>({
            Operator::Boson(Momentum('k'), Index::TypeA, false),
            Operator::Boson(Momentum('k', -1), Index::TypeB, false)
        }));

	Term H_bogo_conjugate(1, Coefficient::HoneyComb("\\gamma", Momentum('k', -1), false), SumContainer{ MomentumSum({'k'}) },
		std::vector<Operator>({
            Operator::Boson(Momentum('k'), Index::TypeA, true),
            Operator::Boson(Momentum('k', -1), Index::TypeB, true)
        }));

	term_vec H = {H_diag, H_bogo, H_bogo_conjugate};

    Term input(1, std::vector<Operator>({
        Operator::Boson(Momentum('q'), Index::TypeA, false),
        Operator::Boson(Momentum('q', -1), Index::TypeB, false)
    }));

    term_vec commutation_result = commutator(H, input);
	clean_up(commutation_result);

	std::cout << "H = " << H << "\\\\" << std::endl;
    std::cout << "[ H, " << input << "] ="  << commutation_result << std::endl;
}

int main(int argc, char** argv) {
	// Remove comment for boson test
	//boson_test();
	//return 0;

	const std::string save_folder = "../commutators/";
	/* WickTerm parse_test("1 sum:momentum{p,q} c:V{p;} o:n{k-p-3x;up} o:f{k+l;}");
	std::cout << parse_test << "    " << parse_test.coefficients.size() << std::endl;
	return 0; */
	constexpr bool print = true;
	constexpr bool print_terms = false;
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
			hermitian_conjugate(vec);
			rename_momenta(vec, 'k', 'l');
		}
		for (auto& vec : disp_daggered) {
			hermitian_conjugate(vec);
			rename_momenta(vec, 'k', 'l');
		}


		const Term H_U(1, Coefficient("\\frac{U}{N}"), MomentumSum({ 'r', 'p', 'q' }), std::vector<Operator>({
			Operator('r', 1, false, Index::SpinUp, true), Operator('p', 1, false, Index::SpinDown, true),
			Operator(momentum_symbols({ MomentumSymbol(1, 'p'), MomentumSymbol(-1, 'q') }), Index::SpinDown, false),
			Operator(momentum_symbols({ MomentumSymbol(1, 'r'), MomentumSymbol(1, 'q') }), Index::SpinUp, false),
			}));
		const std::vector<Term> H = { H_U };
		const int inner_idx = 0;
		const int outer_idx = static_cast<int>(!inner_idx);

		term_vec commute_with_H_base = commutator(H, base[inner_idx]);
		clean_up(commute_with_H_base);

		term_vec commute_with_H_disp = commutator(disp[inner_idx], H);
		clean_up(commute_with_H_disp);
		
		if (true) {
			term_vec joined = joinVectors(commute_with_H_base, commute_with_H_disp);
			remove_all_x(joined);
			clean_up(joined);

			std::cout << "Single commutator:\n" << joined << std::endl; // Up to here, everything works
		}

		{
			term_vec base_double = commutator(base_daggered[outer_idx], commute_with_H_base);
			clean_up(base_double);

			term_vec disp_double = commutator(disp_daggered[outer_idx], commute_with_H_disp);
			clean_up(disp_double);
			
			term_vec joined = joinVectors(base_double, disp_double);
			clean_up(joined);
			//std::cout << "joined:\n" << joined << std::endl;

			auto templates = hubbard.templates();
			auto symmetries = hubbard.symmetries();

			WickTermCollector wicks;
			wicks_theorem(joined, templates, wicks);
			clear_etas(wicks);
			clean_wicks(wicks, symmetries);
			remove_all_x(wicks);
			clean_wicks(wicks, symmetries);

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
		hermitian_conjugate(t);
		rename_momenta(t, 'k', 'l');
		//if (debug) {
		//	rename_momenta(t, 'x', 'y');
		//}
	}

	if (print) std::cout << "\\begin{align*}\n\t H =" << H << "\\end{align*}" << std::endl;

	for (size_t i = 0U; i < basis.size(); ++i)
	{
		term_vec commute_with_H = commutator(H, basis[i]);
		clean_up(commute_with_H);
		if (debug)
			std::cout << "\\begin{align*}\n\t[ H, " << to_string_without_prefactor(basis[i]) << " ] ="
			<< commute_with_H << "\\end{align*}" << std::endl;

		for (size_t j = 0U; j < basis.size(); ++j)
		{
			if (print) {
				std::cout << "\\subsection{" << i << "." << j << "}" << std::endl;
			}
			term_vec terms = commutator(basis_daggered[j], commute_with_H);
			clean_up(terms);

			if (print_terms || debug)
				std::cout << "\\begin{align*}\n\t[ " << to_string_without_prefactor(basis_daggered[j])
				<< ", [H, " << to_string_without_prefactor(basis[i]) << " ]] =" << terms << "\\end{align*}" << std::endl;

			WickTermCollector wicks;
			wicks_theorem(terms, templates, wicks);
			clear_etas(wicks);
			clean_wicks(wicks, symmetries);
			
			for (auto& wickterm : wicks) {
				if (MODEL_TYPE != continuum_type) break;

				Coefficient& current_coeff = wickterm.coefficients.front();
				if (current_coeff.name == "\\rho") {
					wickterm.sums.push_back(Index::SigmaPrime);
					wickterm.sums.push_back('q');
					current_coeff = Coefficient::parse_string("V{0;}");
					wickterm.operators.push_back(WickOperator("n{q;sigma'}"));
					wickterm.multiplicity *= 2;
				}
				if (current_coeff.name == "\\epsilon_{C.Fock}") {
					wickterm.sums.push_back('q');
					current_coeff = Coefficient::parse_string("V{q;}");
					wickterm.operators.push_back(WickOperator("n{k+q;up}"));
					wickterm.multiplicity *= -2;
				}

				if (current_coeff.name == "\\mu_{Ph}") {
					wickterm.sums.push_back(Index::SigmaPrime);
					wickterm.sums.push_back('q');
					current_coeff = Coefficient::parse_interaction_string("U_\\\\mathrm\\{CUT\\}{k,q,0;}");
					wickterm.operators.push_back(WickOperator("n{q;sigma'}"));
					wickterm.multiplicity *= 2;
				}
				if (current_coeff.name == "\\epsilon_{Phock}") {
					wickterm.sums.push_back('q');
					current_coeff = Coefficient::parse_interaction_string("U_\\\\mathrm\\{CUT\\}{k,k+q,q;}");
					wickterm.operators.push_back(WickOperator("n{k+q;up}"));
					wickterm.multiplicity *= -2;
				}

				if (current_coeff.name == "U_\\mathrm{CUT}") {
					assert(current_coeff.momenta.size() == 3U);
					if (current_coeff.momenta[0] == -current_coeff.momenta[1]) {
						// U_CUT (-k, k, x), i.e., BCS channel
						Momentum _first = current_coeff.momenta[0] + current_coeff.momenta[2];
						Momentum _second = current_coeff.momenta[1] - current_coeff.momenta[2];
						// We expect always something like -k - (-k + l) = l
						assert(_first.size() == 1U && _second.size() == 1U);
						assert(_first == -_second);

						current_coeff.momenta[0] = Momentum(current_coeff.momenta[2].front());
						current_coeff.momenta[1] = Momentum(current_coeff.momenta[2].back());
						current_coeff.momenta.pop_back();
						current_coeff.name = "g";
						current_coeff.inversion_symmetry = true;
						current_coeff.custom_symmetry = std::function<void(Coefficient&)>([](Coefficient& _c) {
								_c.momenta.sort();
							});
						wickterm.multiplicity *= -1;
					}
					else if (current_coeff.momenta.back().empty()) {
						// U_CUT (k, l, 0) = const = -M^2 / omega_D =: -G, i.e., Hartree channel
						current_coeff.momenta = MomentumList();
						current_coeff.name = "G";
						wickterm.multiplicity *= -1;
					}
					else {
						Momentum _first = current_coeff.momenta[0] + current_coeff.momenta[2];
						Momentum _second = current_coeff.momenta[1] - current_coeff.momenta[2];
						if (_first == current_coeff.momenta[1] && _second == current_coeff.momenta[0]) {
							// U_CUT(k, l, l - k) or U_CUT(-k, -l, k - l)
							current_coeff.momenta.pop_back();
							current_coeff.momenta[0] = _second;
							current_coeff.momenta[1] = _first;
							current_coeff.inversion_symmetry = true;
							current_coeff.name = "\\tilde{g}";
						}
					} 
				}
			}
			clean_wicks(wicks, symmetries);

			if (debug || print) {
				std::cout << "\\begin{align*}\n\t\\langle [ " << to_string_without_prefactor(basis_daggered[j])
					<< ", [H, " << to_string_without_prefactor(basis[i]) << " ]] \\rangle =" << wicks << "\\end{align*}" << std::endl;
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
			terms = commutator(basis_daggered[j], basis[i]);
			clean_up(terms);
			wicks_theorem(terms, templates, wicks);
			clear_etas(wicks);
			clean_wicks(wicks, symmetries);

			if (debug || print)
				std::cout << "\\begin{align*}\n\t[ " << to_string_without_prefactor(basis_daggered[j])
				<< ", " << to_string_without_prefactor(basis[i]) << " ] =" << wicks << "\\end{align*}" << std::endl;
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

std::cout << "\\begin{align*}\n\t\\langle [" << left.to_string_without_prefactor() << ", " << right.to_string_without_prefactor() << "] \\rangle = " << wicks << "\\end{align*}" << std::endl;
std::cout << "\\begin{align*}\n\t\\langle [" << left.to_string_without_prefactor() << ", [ H, " << right.to_string_without_prefactor() << "] ] \\rangle = " << wicks << "\\end{align*}" << std::endl;
*/

/* Example of how to to read the output
// create an input file stream and a text archive to deserialize the vector
std::ifstream ifs("wick_terms.txt");
boost::archive::text_iarchive ia(ifs);
WickTermCollector deserialized_terms;
ia >> deserialized_terms;
ifs.close();
*/