#include "IndexWrapper.hpp"

namespace SymbolicOperators{
    std::ostream& operator<<(std::ostream& os, const IndexWrapper& indizes){
        for(const auto& idx : indizes){
            os << idx << " ";
        }
        return os;
    }
}