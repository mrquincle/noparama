#pragma once

#include <random>

#include <statistics/distribution.h>

#include <membertrix>

/**
 * This class NealAlgorithm2 updates the cluster population as ...
 *
 */
class NealAlgorithm2: public UpdateClusterPopulation {
	private:
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
		NealAlgorithm2(
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
