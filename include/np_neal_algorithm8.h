#pragma once

#include <random>

#include <statistics/distribution.h>
#include <membertrix.h>
#include <np_update_cluster_population.h>
#include <statistics/dirichlet.h>

struct statistics_t {
	int new_clusters_events;
};

/**
 * This class NealAlgorithm8 updates the cluster population assuming a dirichlet process as nonparametric prior.
 *
 */
class NealAlgorithm8: public UpdateClusterPopulation {
	private:
		//! Random generator to sample values between 0 and 1 for the Metropolis Hasting step.
		std::default_random_engine _generator;

		//! The likelihood 
		distribution_t & _likelihood;
	
		/*! 
		 * The dirichlet process, stored by reference.
		 * If not stored by reference, the constructor would make a copy, loosing all information on the derived class.
		 */
		dirichlet_process & _nonparametrics;

		/*!
		 * The same value as _nonparametrics.getSuffies().alpha. Will be set in the constructor.
		 */
		double _alpha;

		// verbosity
		char _verbosity;

		// statistics
		statistics_t _statistics;
		
	public:
		/*!
		 * Construct update method for cluster population.
		 *
		 * @param[in] generator                  Random generator
		 * @param[in] likelihood                 Likelihood function to be used in update()
		 * @param[in] nonparametrics             Sufficient statistics of the nonparametric prior (e.g. Dirichlet)
		 */
		NealAlgorithm8(
			random_engine_t & generator,
			distribution_t & likelihood,
			dirichlet_process & nonparametrics
		);

		/*!
		 * Update the cluster population. The observation has to be deleted beforehand.
		 *
		 * @param[inout] cluster_matrix          Cluster-observation membership matrix
		 * @param[in] data_id                    Observation to be considered for existing and new cluster
		 */
		void update(
			membertrix & cluster_matrix,
			data_id_t data_id
		);

		/*!
		 * MCMC statistics. How many times are steps rejected, etc.
		 */
		void printStatistics();
};
