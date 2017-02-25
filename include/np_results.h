#include <string>
#include <unordered_map>

#include <Eigen/Dense>

#include "common.h"
#include "membertrix.h"
#include "clustering_performance.h"

typedef int cluster_id_t;

typedef std::unordered_map<cluster_id_t, data_id_t> ground_truth_t;

class Results {
    private:
        const membertrix & _membertrix;

        ground_truth_t & _ground_truth;

        int _verbosity;

        clustering_performance _clustering_performance;
    public:
        Results(const membertrix &membertrix, ground_truth_t & ground_truth);

        matrix_t & calculateContingencyMatrix();

        void write(const std::string & workspace, const std::string & path, const std::string & basename);
};
