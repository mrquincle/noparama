#pragma once

#include <vector>
#include <random>

#include <membertrix.h>

#include <np_data.h>

#include <np_init_clusters.h>
#include <np_update_clusters.h>
#include <np_update_cluster_population.h>

class MCMC {
	private:
		//! Todo: set generator
		std::default_random_engine & _generator;
	
		//! Membertrix object itself
		membertrix _membertrix;
		
		//! Reference to InitCluster object
		InitClusters & _init_clusters;

		//! Reference to UpdateCluster object
		UpdateClusters & _update_clusters;

		//! Reference to UpdateClusterPopulation object
		UpdateClusterPopulation & _update_cluster_population;

		// verbosity
		char _verbosity;

	public:
		MCMC(
				std::default_random_engine & generator,
				InitClusters & init_clusters, 
				UpdateClusters & update_clusters, 
				UpdateClusterPopulation & update_cluster_population
			);

		/**
		 * Run the MCMC method (can be one or multiple chains) for T steps.
		 */
		void run(dataset_t & dataset, int T);

		/**
		 * The results of the MCMC algorithm is a membership matrix.
		 */
		const membertrix & getMembershipMatrix() const;
};
