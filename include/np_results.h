#include <string>
#include <unordered_map>

#include <Eigen/Dense>

#include "membertrix.h"

typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

typedef std::unordered_map<int, data_id_t> ground_truth_t;

class Results {
    private:
        const membertrix & _membertrix;

        ground_truth_t & _ground_truth;

        double _rand_index;

        double _adjusted_rand_index;

        int _verbosity;
    public:
        Results(const membertrix &membertrix, ground_truth_t & ground_truth);

        matrix_t & calculateContingencyMatrix();

        void calculateSimilarity();

        void write(const std::string & workspace, const std::string & path, const std::string & basename);
};
