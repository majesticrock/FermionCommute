#pragma once
#include "Hubbard.hpp"

#include <string>
#include <vector>

namespace mrock::symbolic_operators {
	struct HubbardDispersions : public Hubbard {
		virtual std::vector<TermCollector> XP_basis() const override;
		virtual std::vector<TermCollector> STD_basis() const override;

		virtual std::string get_subfolder() const override;
	};
}