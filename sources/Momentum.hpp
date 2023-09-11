#pragma once
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/string.hpp>

namespace SymbolicOperators {
	template <typename T>
	inline bool operator==(const std::vector<T>& lhs, const std::vector<T>& rhs) {
		if (lhs.size() != rhs.size()) return false;
		for (size_t i = 0; i < lhs.size(); i++)
		{
			if (lhs[i] != rhs[i]) return false;
		}
		return true;
	}
	template <typename T>
	inline bool operator!=(const std::vector<T>& lhs, const std::vector<T>& rhs) {
		return !(lhs == rhs);
	}

	typedef std::vector<std::pair<int, char>> momentum_pairs;
	struct Momentum {
		// total momentum is then:    sum_i pair_i.first * pair_i.second
		momentum_pairs momentum_list;
		bool add_Q;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& momentum_list;
			ar& add_Q;
		}

		Momentum();
		explicit Momentum(const char value, int plus_minus = 1, bool Q = false);
		explicit Momentum(const momentum_pairs& _momenta, bool Q = false);

		void sort();

		Momentum& operator+=(const Momentum& rhs);
		Momentum& operator-=(const Momentum& rhs);

		inline void multiplyMomentum(int factor) {
			for (auto& m : momentum_list) {
				m.first *= factor;
			}
		};
		inline void flipMomentum() {
			multiplyMomentum(-1);
		};
		// Returns the position in the momentum_list array of the momentum given in <value>
		// Return -1 if the value is not contained in momentum_list
		inline int isUsed(const char value) const noexcept {
			for (int i = 0; i < momentum_list.size(); ++i) {
				if (momentum_list[i].second == value) return i;
			}
			return -1;
		};
		inline bool differsOnlyInQ(Momentum rhs) const {
			if (rhs.add_Q == this->add_Q) return false;
			rhs.add_Q = this->add_Q;
			if (*this != rhs) return false;
			return true;
		};

		void addInPlace(const Momentum& rhs);
		void replaceOccurances(const char replaceWhat, const Momentum& replaceWith);
		inline bool operator==(const Momentum& rhs) const {
			if (this->add_Q != rhs.add_Q) return false;
			if (this->momentum_list.size() != rhs.momentum_list.size()) return false;
			bool foundOne = true;
			for (size_t i = 0; i < this->momentum_list.size(); i++)
			{
				foundOne = false;
				for (size_t j = 0; j < rhs.momentum_list.size(); j++)
				{
					if (this->momentum_list[i] == rhs.momentum_list[j])
						foundOne = true;
				}
				if (!foundOne) return false;
			}
			return true;
		};
		inline bool operator!=(const Momentum& rhs) const {
			return !(*this == rhs);
		};
	};

	inline Momentum operator+(Momentum lhs, const Momentum& rhs) {
		lhs += rhs;
		return lhs;
	}
	inline Momentum operator-(Momentum lhs, const Momentum& rhs) {
		lhs -= rhs;
		return lhs;
	}

	std::ostream& operator<<(std::ostream& os, const Momentum& momentum);
}