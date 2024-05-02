#include "TermLoader.hpp"
#include <fstream>
#include <boost/archive/text_iarchive.hpp>

namespace SymbolicOperators {
    void TermLoader::load(std::string const& folder, bool use_XP, int n_terms, int start_at/*=0*/){
        M.resize(n_terms * n_terms);
		N.resize(n_terms * n_terms);
		const int name_offset = (start_at < 0) ? 0 : start_at;

		for (int i = 0; i < n_terms; i++)
		{
			for (int j = 0; j < n_terms; j++)
			{
				const std::string suffix = std::to_string(j + name_offset) + "_" + std::to_string(i + name_offset) + ".txt";
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(folder + (use_XP ? "XP_" : "") + "wick_M_" + suffix);
					boost::archive::text_iarchive ia(ifs);
					this->M[j * n_terms + i].clear();
					ia >> this->M[j * n_terms + i];
					ifs.close();
				}
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(folder + (use_XP ? "XP_" : "") + "wick_N_" + suffix);
					boost::archive::text_iarchive ia(ifs);
					this->N[j * n_terms + i].clear();
					ia >> this->N[j * n_terms + i];
					ifs.close();
				}
			}
		}
    }
}