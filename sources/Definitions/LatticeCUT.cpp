#include "LatticeCUT.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace mrock::symbolic_operators {
std::vector<Term> LatticeCUT::hamiltonian() const {
    const Term H_Kin(1, Coefficient("\\epsilon_0", Momentum('K')), SumContainer{MomentumSum({'K'}), Index::Sigma},
                     std::vector<Operator>(
                         {Operator('K', 1, false, Index::Sigma, true), Operator('K', 1, false, Index::Sigma, false)}));

    const Term H_Ph(-1,
                    Coefficient::RealInversionSymmetric(
                        "g", MomentumList({'K', 'P'}),
                        std::function<void(Coefficient&)>([](Coefficient& coeff) { coeff.momenta.sort(); })),
                    SumContainer{MomentumSum({'K', 'P'}), IndexSum{}},
                    std::vector<Operator>({c_k_dagger.with_momentum('K'), c_minus_k_dagger.with_momentum('K'),
                                           c_minus_k.with_momentum('P'), c_k.with_momentum('P')}));

    const Term H_U(1, Coefficient("U"), SumContainer{MomentumSum({'P', 'K'}), IndexSum{}},
                   std::vector<Operator>({c_k_dagger.with_momentum('K'), c_minus_k_dagger.with_momentum('K'),
                                          c_minus_k.with_momentum('P'), c_k.with_momentum('P')}));

    return {H_Kin, H_Ph, H_U};
}

std::vector<WickOperatorTemplate> LatticeCUT::templates() const {
    return {WickOperatorTemplate{{SC_Comparison}, Momentum(), OperatorType::SC},
            WickOperatorTemplate{{Num_Comparison}, Momentum(), OperatorType::Number}};
}

std::vector<std::vector<Term>> LatticeCUT::XP_basis() const {
    return {// 0: f + f^+
            std::vector<Term>({Term(1, std::vector<Operator>({c_minus_k, c_k})),
                               Term(1, std::vector<Operator>({c_k_dagger, c_minus_k_dagger}))}),
            // 1: n_up + down
            std::vector<Term>({Term(1, std::vector<Operator>({c_k_dagger, c_k})),
                               Term(1, std::vector<Operator>({c_minus_k_dagger, c_minus_k}))}),
            // 2: f - f^+
            std::vector<Term>({Term(1, std::vector<Operator>({c_minus_k, c_k})),
                               Term(-1, std::vector<Operator>({c_k_dagger, c_minus_k_dagger}))})};
}

std::vector<std::vector<Term>> LatticeCUT::STD_basis() const {
    return {// f, f^+
            std::vector<Term>({Term(1, std::vector<Operator>({c_minus_k, c_k}))}),
            std::vector<Term>({Term(1, std::vector<Operator>({c_k_dagger, c_minus_k_dagger}))}),
            // n_up + down
            std::vector<Term>({Term(1, std::vector<Operator>({c_k_dagger, c_k})),
                               Term(1, std::vector<Operator>({c_minus_k_dagger, c_minus_k}))})};
}

std::vector<std::unique_ptr<WickSymmetry>> LatticeCUT::symmetries() const {
    std::vector<std::unique_ptr<WickSymmetry>> ret;
    ret.push_back(std::make_unique<SpinSymmetry>());
    ret.push_back(std::make_unique<InversionSymmetry>());
    ret.push_back(std::make_unique<PhaseSymmetry<OperatorType::SC>>());
    return ret;
}

std::string LatticeCUT::get_subfolder() const {
    return "lattice_cut/";
}
}  // namespace mrock::symbolic_operators