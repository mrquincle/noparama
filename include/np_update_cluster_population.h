#pragma once

#include <random>
#include <iostream>

#include <membertrix.h>

/**
 * This class UpdateClusterPopulation deletes, adds, and adjusts clusters. This in contrast with UpdateCluster which 
 * only adjusts the parameters assigned to a cluster and will leave the number of clusters invariant. 
 */
class UpdateClusterPopulation {
	protected:
		std::default_random_engine _generator;
	
		std::uniform_real_distribution<double> _distribution;

		// verbosity
		char _verbosity;

	public:
		UpdateClusterPopulation():
			_distribution(0.0, 1.0) {
		}

		/*!
		 * Update the cluster population. The observation has to be deleted beforehand.
		 * @param[inout] cluster_matrix				Cluster-observation membership matrix
		 * @param[in] data_ids						Observations to be considered for existing and new cluster
		 */
		virtual void update(
				membertrix & cluster_matrix,
				std::vector<data_id_t> data_ids
			) {
			std::cout << "General update" << std::endl;
		}
};
