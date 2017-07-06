#pragma once

#include <random>

#include <statistics/distribution.h>
#include <membertrix.h>
#include <np_update_cluster_population.h>
#include <statistics/dirichlet.h>

typedef enum { simple_random_split, sams_prior, sams_random_walk, random_mixing } split_method_t;

typedef struct step {
	int attempts;
	int cluster_events_accept;
	int cluster_events_reject;
	int target_likelihood_zero;
	int source_likelihood_zero;
	int likelihood_both_zero;
	int likelihood_both_nonzero;
	int likelihood_both_nonzero_accept;
	int likelihood_both_nonzero_reject;
} step_t;

struct statistics_t {
	step_t split;
	step_t merge;
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

		//! Store pointer to membership matrix locally
		membertrix * _cluster_matrix;

		/*!
		 * The same value as _nonparametrics.getSuffies().alpha. Will be set in the constructor.
		 */
		double _alpha;

		// verbosity
		char _verbosity;

		// statistics
		statistics_t _statistics;

		double ratioStateProb(bool split, int nc0, int nc1);

		double ratioProposal(bool split, int nc);

		void checkLikelihoods(double l1, double l2, step_t statistics_step, double &l21, bool &accept, bool &overwrite);

		void propose_split(data_ids_t & move, data_ids_t & remain, data_id_t data_i, data_id_t data_j, 
				cluster_id_t current_cluster_id, cluster_t & new_cluster, split_method_t split_method);

		// split a single cluster
		bool split(data_ids_t data_ids, cluster_id_t current_cluster);

		// merge the given clusters
		bool merge(data_ids_t data_ids, std::vector<cluster_id_t> & current_clusters);
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
