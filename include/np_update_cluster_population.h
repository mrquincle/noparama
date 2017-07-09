#pragma once

#include <random>
#include <iostream>

#include <membertrix.h>
#include <np_data.h>

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
				const data_ids_t & data_ids
			) = 0;
		
		/*!
		 * MCMC statistics. How many times are steps rejected, etc.
		 */
		virtual void printStatistics() = 0; 
};
