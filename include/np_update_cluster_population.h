#pragma once

#include <random>

#include <membertrix>

#include <iostream>

/**
 * This class UpdateClusterPopulation deletes, adds, and adjusts clusters. This in contrast with UpdateCluster which 
 * only adjusts the parameters assigned to a cluster and will leave the number of clusters invariant. 
 */
class UpdateClusterPopulation {
	protected:
		std::default_random_engine _generator;
	
		std::uniform_real_distribution<double> _distribution;

	public:
		UpdateClusterPopulation():
			_distribution(0.0, 1.0) {
		}

		/*!
		 * Update the cluster population. The observation has to be deleted beforehand.
		 * @param[inout] cluster_matrix				Cluster-observation membership matrix
		 * @param[in] data_id						Observation to be considered for existing and new cluster
		 */
		virtual void update(
				membertrix & cluster_matrix,
				data_id_t data_id
			) {
			std::cout << "General update" << std::endl;
		}
};