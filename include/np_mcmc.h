#pragma once

#include <vector>
#include <random>

#include <membertrix>

#include <np_data.h>

#include <np_update_clusters.h>
#include <np_update_cluster_population.h>

class MCMC {
	private:

		membertrix _membertrix;

		dataset_t _dataset;

		distribution_t _nonparametrics;

		UpdateClusters _update_clusters;

		UpdateClusterPopulation _update_cluster_population;

		std::default_random_engine _generator;
	public:
		MCMC(UpdateClusters & update_clusters, UpdateClusterPopulation & update_cluster_population);

		void run(dataset_t & dataset);
};
