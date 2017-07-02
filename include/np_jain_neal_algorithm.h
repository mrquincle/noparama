#pragma once

#include <random>

#include <statistics/distribution.h>
#include <membertrix.h>
#include <np_update_cluster_population.h>
#include <statistics/dirichlet.h>

struct statistics_t {
	// number of new clusters
	int split_attempts;
	int split_cluster_events_accept;
	int split_cluster_events_reject;
	int split_target_likelihood_zero;
	int split_source_likelihood_zero;
	int split_likelihood_both_zero;
	int split_likelihood_both_nonzero;
	int split_likelihood_both_nonzero_accept;
	int split_likelihood_both_nonzero_reject;

	int merge_attempts;
	int merge_cluster_events_accept;
	int merge_cluster_events_reject;
	int merge_target_likelihood_zero;
	int merge_source_likelihood_zero;
	int merge_likelihood_both_zero;
	int merge_likelihood_both_nonzero;
	int merge_likelihood_both_nonzero_accept;
	int merge_likelihood_both_nonzero_reject;
};

/**
 * This class JainNealAlgorithm updates the cluster population assuming a dirichlet process as nonparametric prior.
 *
 */
class JainNealAlgorithm: public UpdateClusterPopulation {
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
		
		bool split(membertrix & cluster_matrix, data_ids_t data_ids, std::vector<cluster_id_t> & prev_clusters);
		bool merge(membertrix & cluster_matrix, data_ids_t data_ids, std::vector<cluster_id_t> & prev_clusters);
	public:
		/*!
		 * Construct update method for cluster population.
		 *
		 * @param[in] generator                  Random generator
		 * @param[in] likelihood                 Likelihood function to be used in update()
		 * @param[in] nonparametrics             Sufficient statistics of the nonparametric prior (e.g. Dirichlet)
		 */
		JainNealAlgorithm(
			random_engine_t & generator,
			distribution_t & likelihood,
			dirichlet_process & nonparametrics
		);

		/*!
		 * Update the cluster population. The observation has to be deleted beforehand.
		 *
		 * @param[inout] cluster_matrix          Cluster-observation membership matrix
		 * @param[in] data_ids                   Observations to be considered for existing and new cluster
		 */
		void update(
			membertrix & cluster_matrix,
			std::vector<data_id_t> data_ids
		);

		/*!
		 * MCMC statistics. How many times are steps rejected, etc.
		 */
		void printStatistics();
};
