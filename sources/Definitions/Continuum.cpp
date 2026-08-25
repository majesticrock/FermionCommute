#include "Continuum.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

// #define COULOMB_ONLY_SC_CHANNEL
#define PHONON_ONLY_SC_CHANNEL

#ifdef PHONON_ONLY_SC_CHANNEL
#define PHONON_HAMILTONIAN H_Ph
#else
#define PHONON_HAMILTONIAN H_Ph, H_Phock, H_Phartree
#endif

#ifdef COULOMB_ONLY_SC_CHANNEL
#define COULOMB_HAMILTONIAN H_C
#else
#define COULOMB_HAMILTONIAN H_C, H_BG  //, H_C_Fock
#endif

namespace mrock::symbolic_operators {
std::vector<Term> Continuum::hamiltonian() const {
    const Term H_Kin(1, Coefficient("\\epsilon_0", Momentum('K')), SumContainer{MomentumSum({'K'}), Index::Sigma},
                     std::vector<Operator>(
                         {Operator('K', 1, false, Index::Sigma, true), Operator('K', 1, false, Index::Sigma, false)}));

#ifndef PHONON_ONLY_SC_CHANNEL
    const Term H_Ph(
        IntFractional(1, 2), Coefficient::RealInteraction("U_\\mathrm{CUT}", MomentumList({'K', 'P', 'Q'})),
        SumContainer{MomentumSum({'K', 'P', 'Q'}), IndexSum({Index::Sigma, Index::SigmaPrime})},
        std::vector<Operator>(
            {Operator(Momentum("K+Q"), Index::Sigma, true), Operator(Momentum("P-Q"), Index::SigmaPrime, true),
             Operator(Momentum('P'), Index::SigmaPrime, false), Operator(Momentum('K'), Index::Sigma, false)}));

    const Term H_Phock(-IntFractional(1, 2), Coefficient("\\epsilon_{Phock}", Momentum('Q')),
                       SumContainer{MomentumSum({'Q'}), Index::Sigma},
                       std::vector<Operator>({Operator('Q', 1, false, Index::Sigma, true),
                                              Operator('Q', 1, false, Index::Sigma, false)}));

    const Term H_Phartree(-IntFractional(1, 2), Coefficient("\\mu_{Ph}"),
                          SumContainer{MomentumSum({'Q'}), Index::Sigma},
                          std::vector<Operator>({Operator(Momentum("q"), Index::Sigma, true),
                                                 Operator('Q', 1, false, Index::Sigma, false)}));
#else
    const Term H_Ph(-1,
                    Coefficient::RealInversionSymmetric(
                        "g", MomentumList({'K', 'P'}),
                        std::function<void(Coefficient&)>([](Coefficient& coeff) { coeff.momenta.sort(); })),
                    SumContainer{MomentumSum({'K', 'P'}), IndexSum{}},
                    std::vector<Operator>({c_k_dagger.with_momentum('K'), c_minus_k_dagger.with_momentum('K'),
                                           c_minus_k.with_momentum('P'), c_k.with_momentum('P')}));
#endif

#ifndef COULOMB_ONLY_SC_CHANNEL
    const Term H_C(IntFractional(1, 2), Coefficient("V", Momentum('Q')),
                   SumContainer{MomentumSum({'K', 'P', 'Q'}), IndexSum({Index::Sigma, Index::SigmaPrime})},
                   std::vector<Operator>(
                       {Operator('K', 1, false, Index::Sigma, true), Operator('P', 1, false, Index::SigmaPrime, true),
                        Operator(std::vector<MomentumSymbol>({MomentumSymbol(1, 'P'), MomentumSymbol(-1, 'Q')}),
                                 Index::SigmaPrime, false),
                        Operator(std::vector<MomentumSymbol>({MomentumSymbol(1, 'K'), MomentumSymbol(1, 'Q')}),
                                 Index::Sigma, false)}));

    // const Term H_C_Fock(-IntFractional(1, 2), Coefficient("\\epsilon_{C.Fock}", Momentum('K')), SumContainer{
    // MomentumSum({ 'K' }), Index::Sigma }, 	std::vector<Operator>({ 		Operator('K', 1, false,
    // Index::Sigma, true),
    // Operator('K', 1, false, Index::Sigma, false)
    //		}));

    const Term H_BG(-IntFractional(1, 2), Coefficient("\\rho"), SumContainer{MomentumSum({'K'}), Index::Sigma},
                    std::vector<Operator>(
                        {Operator(Momentum('K'), Index::Sigma, true), Operator('K', 1, false, Index::Sigma, false)}));

#else
    const Term H_C(1, Coefficient("V", Momentum('P')), SumContainer{MomentumSum({'K', 'P'}), IndexSum{}},
                   std::vector<Operator>({
                       Operator('P', 1, false, Index::SpinUp, true),
                       Operator('P', -1, false, Index::SpinDown, true),
                       Operator(std::vector<MomentumSymbol>({MomentumSymbol(-1, 'K'), MomentumSymbol(-1, 'P')}),
                                Index::SpinDown, false),
                       Operator(std::vector<MomentumSymbol>({MomentumSymbol(1, 'K'), MomentumSymbol(1, 'P')}),
                                Index::SpinUp, false),
                   }));
#endif
    return {H_Kin, PHONON_HAMILTONIAN, COULOMB_HAMILTONIAN};
}
std::vector<WickOperatorTemplate> Continuum::templates() const {
    return {WickOperatorTemplate{{SC_Comparison}, Momentum(), OperatorType::SC},
            WickOperatorTemplate{{Num_Comparison}, Momentum(), OperatorType::Number}};
}
std::vector<std::vector<Term>> Continuum::XP_basis() const {
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
std::vector<std::vector<Term>> Continuum::STD_basis() const {
    return {// f, f^+
            std::vector<Term>({Term(1, std::vector<Operator>({c_minus_k, c_k}))}),
            std::vector<Term>({Term(1, std::vector<Operator>({c_k_dagger, c_minus_k_dagger}))}),
            // n_up + down
            std::vector<Term>({Term(1, std::vector<Operator>({c_k_dagger, c_k})),
                               Term(1, std::vector<Operator>({c_minus_k_dagger, c_minus_k}))})};
}
std::vector<std::unique_ptr<WickSymmetry>> Continuum::symmetries() const {
    std::vector<std::unique_ptr<WickSymmetry>> ret;
    ret.push_back(std::make_unique<SpinSymmetry>());
    ret.push_back(std::make_unique<InversionSymmetry>());
    ret.push_back(std::make_unique<PhaseSymmetry<OperatorType::SC>>());
    return ret;
}
std::string Continuum::get_subfolder() const {
    return "continuum/";
}
}  // namespace mrock::symbolic_operators