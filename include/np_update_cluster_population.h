#pragma once

#include <random>

#include <statistics/distribution.h>

#include <membertrix>

/**
 * This class UpdateClusterPopulation deletes, adds, and adjusts clusters. This in contrast with UpdateCluster which 
 * only adjusts the parameters assigned to a cluster and will leave the number of clusters invariant. 
 */
class UpdateClusterPopulation {
	private:
		std::default_random_engine _generator;
	
		std::uniform_real_distribution<double> _distribution;

		distribution_t _likelihood;
		
		distribution_t _nonparametrics;
		
		distribution_t _prior;

	public:
		/*!
		 * Construct update method for cluster population.
		 * @param[in] likelihood					Likelihood function to be used in update()
		 * @param[in] nonparametrics 				Sufficient statistics of the nonparametric prior (e.g. Dirichlet)
		 * @param[in] prior							Prior to be used in update()
		 */
		UpdateClusterPopulation(
				distribution_t & likelihood,
				distribution_t & nonparametrics, 
				distribution_t & prior
			);

		/*!
		 * Update the cluster population. The observation has to be deleted beforehand.
		 * @param[inout] cluster_matrix				Cluster-observation membership matrix
		 * @param[in] data_id						Observation to be considered for existing and new cluster
		 */
		void update(
				membertrix & cluster_matrix,
				data_id_t data_id
			);
};
