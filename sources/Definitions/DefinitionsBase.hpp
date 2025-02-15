#pragma once

#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/Operator.hpp>
#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/WickOperatorTemplate.hpp>
#include <mrock/symbolic_operators/WickSymmetry.hpp>

namespace mrock::symbolic_operators {
	extern const Momentum base_k;
	extern const Momentum base_x;
	extern const Momentum base_k_Q;

	extern const Operator c_k; // c_{k up}
	extern const Operator c_minus_k; // c_{-k down}

	extern const Operator c_k_dagger; // c_{k up}^+
	extern const Operator c_minus_k_dagger; // c_{-k down}^+

	extern const Operator c_k_Q;
	extern const Operator c_minus_k_Q;

	extern const Operator c_k_Q_dagger;
	extern const Operator c_minus_k_Q_dagger;

	// transversal magnon
	extern const Operator c_k_Q_down_dagger;
	extern const Operator c_k_Q_down;

	// Ferromagnetic modes
	extern const Operator c_k_down_dagger;
	extern const Operator c_k_down;

	struct DefinitionsBase {
		virtual std::vector<Term> hamiltonian() const = 0;
		virtual std::vector<WickOperatorTemplate> templates() const = 0;
		virtual std::vector<std::vector<Term>> XP_basis() const = 0;
		virtual std::vector<std::vector<Term>> STD_basis() const = 0;
		virtual std::vector<std::unique_ptr<WickSymmetry>> symmetries() const = 0;

		virtual std::string get_subfolder() const = 0;
	};
}