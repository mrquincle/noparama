#pragma once

#include <random>

#include <Eigen/Dense>

#include <statistics/distribution.h>

#include <membertrix.h>

#include <np_cluster.h>

/**
 * This class InitClusters initiates clusters using hyperparameters in the form of a non-parametric Bayesian prior.
 */
class InitClusters {
	private:
		std::default_random_engine _generator;

		distribution_t &  _nonparametrics;

		std::uniform_real_distribution<double> _distribution;

		// verbosity
		char _verbosity;
	public:
		/*!
		 * Constructor for the InitClusters class. 
		 *
		 * @param[in] nonparametrics             Sufficient statistics of the nonparametric prior (e.g. Dirichlet)
		 */
		InitClusters(
				random_engine_t & generator,
				distribution_t & nonparametrics
				);

		/*!
		 * Init all cluster parameters. 
		 *
		 * @param[inout] cluster_matrix          Cluster parameters that will be updated
		 */
		void init(
				membertrix & cluster_matrix,
				int K
			);
};

