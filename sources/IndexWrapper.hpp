#pragma once
#include <iostream>

namespace SymbolicOperators{
    enum Index{ UP, DOWN, Sigma, SigmaPrime, UndefinedIndex };

	inline std::ostream& operator<<(std::ostream& os, const Index) {
		switch (Spin)
		{
		case UP:
			os << "\\uparrow";
			break;
		case DOWN:
			os << "\\downarrow";
			break;
		case Sigma:
			os << "\\sigma";
			break;
		case SigmaPrime:
			os << "\\sigma'";
			break;
		default:
			os << "ERROR_INDEX";
			break;
		}
		return os;
	};

    struct IndexWrapper{
        std::vector<Index> indizes;

		inline Index& operator[](size_t i) {
			return indizes[i];
		};
		inline Index operator[](size_t i) const {
			return indizes[i];
		};
		inline auto begin() {
			return indizes.begin();
		}
		inline auto begin() const {
			return indizes.begin();
		}
		inline auto end() {
			return indizes.end();
		}
		inline auto end() const {
			return indizes.end();
		}

        IndexWrapper() = default;
        IndexWrapper(Index _spin) : indizes(1U, _spin) {};
        IndexWrapper(const std::vector<Index>& _indizes)
            : indizes(_indizes) {};
    };

	std::ostream& operator<<(std::ostream& os, const IndexWrapper& indizes);
	inline bool operator==(const IndexWrapper& lhs, const IndexWrapper& rhs){
		return lhs.indizes == rhs.indizes;
	};
	inline bool operator!=(const IndexWrapper& lhs, const IndexWrapper& rhs){
		return !(lhs == rhs);
	};
}