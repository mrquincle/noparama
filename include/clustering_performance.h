#pragma once

#include <common.h>

class clustering_performance {
    public:
        clustering_performance();

        matrix_t & calculateContingencyMatrix(std::vector<int> & A, std::vector<int> & B);

        void calculateSimilarity(matrix_t & frequencies);

        void write(std::string fname);

    private:
        double _rand_index;

        double _adjusted_rand_index;

        double _purity;
};
