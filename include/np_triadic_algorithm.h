#pragma once

#include <random>

#include <statistics/distribution.h>
#include <membertrix.h>
#include <np_update_cluster_population.h>
#include <statistics/dirichlet.h>
#include <np_statistics.h>

/**
 * This class TriadicAlgorithm updates the cluster population assuming a dirichlet process as nonparametric prior.
 *
 */
class TriadicAlgorithm: public UpdateClusterPopulation {
	private:
		//! Random generator to sample values between 0 and 1 for the Metropolis Hasting step.
		std::default_random_engine _generator;

		//! The likelihood 
		distribution_t & _likelihood;

		//! Splitting method
		split_method_t _split_method;

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

		double _beta;

		// verbosity
		char _verbosity;

		// statistics
		statistics_t _statistics;

		double ratioStateProb(bool split, const std::vector<int> &more, const std::vector<int> &less);

		double ratioProposal(bool split, int N, int C);

		void propose_merge(std::vector<data_ids_t> &pdata, const data_ids_t &data_ids, 
				cluster_ids_t &cluster_ids, split_method_t split_method);

		void propose_split(std::vector<data_ids_t> &pdata, const data_ids_t &data_ids, 
				cluster_ids_t &cluster_ids, cluster_t *new_cluster, split_method_t split_method);

		// split a single cluster
		bool split(const data_ids_t & data_ids, cluster_ids_t & cluster_ids);

		// merge the given clusters
		bool merge(const data_ids_t & data_ids, cluster_ids_t & clusters_ids);
	public:
		/*!
		 * Construct update method for cluster population.
		 *
		 * @param[in] generator                  Random generator
		 * @param[in] likelihood                 Likelihood function to be used in update()
		 * @param[in] nonparametrics             Sufficient statistics of the nonparametric prior (e.g. Dirichlet)
		 */
		TriadicAlgorithm(
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
			const data_ids_t & data_ids
		);

		/*!
		 * MCMC statistics. How many times are steps rejected, etc.
		 */
		void printStatistics();
};
