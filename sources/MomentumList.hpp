#pragma once
#include "Momentum.hpp"
#include "../../Utility/sources/VectorWrapper.hpp"
#include <algorithm>

namespace SymbolicOperators {
	struct MomentumList : public Utility::VectorWrapper<Momentum>
	{
	private:
		using _parent = Utility::VectorWrapper<Momentum>;

	public:
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->_vector;
		};

		MomentumList() : _parent() {};
		explicit MomentumList(const char value, int plus_minus = 1, bool Q = false)
			: _parent{ Momentum(value, plus_minus, Q) } {};
		explicit MomentumList(const momentum_pairs& _momenta, bool Q = false)
			: _parent{ Momentum(_momenta, Q) } {};
		MomentumList(const Momentum& momentum)
			: _parent{ momentum } {};
		MomentumList(const Momentum& first, const Momentum& second)
			: _parent{ first, second } {};

		inline MomentumList& operator*=(const int rhs) {
			for (auto& mom : _vector) {
				mom *= rhs;
			}
			return *this;
		};
		inline void multiplyMomentum(int factor) {
			(*this) *= factor;
		};
		inline void flipMomentum() {
			(*this) *= -1;
		};

		void replaceOccurances(const char replaceWhat, const Momentum& replaceWith);
		void remove_zeros();
		void flip_single(char momentum);
	};

	std::ostream& operator<<(std::ostream& os, const MomentumList& momenta);
}