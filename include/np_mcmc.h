#pragma once

#include <vector>
#include <random>

#include <membertrix>

#include <np_data.h>

#include <np_update_clusters.h>
#include <np_update_cluster_population.h>

class MCMC {
	private:
		//! Membertrix object itself
		membertrix _membertrix;

		//! Reference to UpdateCluster object
		UpdateClusters & _update_clusters;

		//! Reference to UpdateClusterPopulation object
		UpdateClusterPopulation & _update_cluster_population;

		//! Todo: set generator
		std::default_random_engine & _generator;
	
		//! Reference to nonparametrics distribution_t object
		distribution_t & _nonparametrics;
		
		// verbosity
		char _verbosity;

	public:
		MCMC(
				UpdateClusters & update_clusters, 
				UpdateClusterPopulation & update_cluster_population,
				std::default_random_engine & generator,
				distribution_t & nonparametrics
			);

		void run(dataset_t & dataset);
};
