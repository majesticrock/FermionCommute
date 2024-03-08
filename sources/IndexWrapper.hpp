#pragma once
#include <iostream>
#include <vector>

namespace SymbolicOperators{
    enum Index{ UP, DOWN, Sigma, SigmaPrime, UndefinedIndex };

	std::ostream& operator<<(std::ostream& os, const Index index);

    struct IndexWrapper{
        std::vector<Index> indizes;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& indizes;
		}

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
		inline bool empty() const noexcept {
			return indizes.empty();
		}
		inline size_t size() const noexcept {
			return indizes.size();
		}
		inline void push_back(const Index index) {
			indizes.push_back(index);
		}
		inline Index front() const {
			return indizes.front();
		}
		inline Index& front() {
			return indizes.front();
		}
		inline Index back() const {
			return indizes.back();
		}
		inline Index& back() {
			return indizes.back();
		}
		template <class iterator>
		inline auto erase(iterator pos) {
			return indizes.erase(pos);
		}
		template <class iterator>
		inline auto erase(iterator first, iterator last) {
			return indizes.erase(first, last);
		}
		template <class iterator>
		inline auto insert(iterator pos, const Index value) {
			return indizes.insert(pos, value);
		}
		template <class iterator, class input_iterator>
		inline auto insert(iterator pos, input_iterator first, input_iterator last) {
			return indizes.insert(pos, first, last);
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