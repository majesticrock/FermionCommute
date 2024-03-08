#include "IndexWrapper.hpp"

namespace SymbolicOperators{
    std::ostream& operator<<(std::ostream& os, const Index index)
	{
		switch (index)
		{
		case UP:
			os << "\\uparrow";
			break;
		case DOWN:
			os << "\\downarrow";
			break;
		case Sigma:
			os << Sigma;
			break;
		case SigmaPrime:
			os << SigmaPrime;
			break;
		default:
			os << "ERROR_INDEX";
			break;
		}
		return os;
	}

    std::ostream& operator<<(std::ostream& os, const IndexWrapper& indizes){
        for(const auto& idx : indizes){
            os << idx << " ";
        }
        return os;
    }
}