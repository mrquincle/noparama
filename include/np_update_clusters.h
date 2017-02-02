#pragma once

#include <random>

#include <Eigen/Dense>

#include <statistics/distribution.h>

#include <np_cluster.h>

/**
 * This class UpdateClusters only adjusts clusters. This in contrast with UpdateClusterPopulation which also
 * deletes and removes clusters and does not leave the number of clusters invariant. 
 */
class UpdateClusters {
	private:
		std::default_random_engine _generator;
	
		std::uniform_real_distribution<double> _distribution;

	protected:

		Suffies & propose();

	public:
		/*!
		 * Constructor for the UpdateCluster class. 
		 *
		 * The class needs access to (1) the likelihood of cluster parameters given observations and (2) the 
		 * nonparametric prior (e.g. the parameters of a Dirichlet Process).
		 *
		 * Mathematical, this class needs access to the likelihood:
		 * 	- L(theta|x0,...xN)
		 *
		 * And it should be possible to sample parameters from the nonparametric prior:
		 * 	- theta ~ p(theta|lambda)
		 *
		 * The cluster parameters are summarized through theta, the hyper parameters are summarized through lambda.
		 *
		 * @param[in] likelihood					Likelihood function to be used in update()
		 * @param[in] nonparametrics 				Sufficient statistics of the nonparametric prior (e.g. Dirichlet)
		 */
		UpdateClusters(
				distribution_t & likelihood,
				distribution_t & nonparametrics
			);

		/*!
		 * Update all cluster parameters. The number of clusters will not change.
		 *
		 * The constructor already sets the relevant (parameters of) probability distributions. The input for this
		 * function are only the current cluster parameters and some parameters for the update procedure such as the
		 * number of steps.
		 *
		 * @param[inout] clusters					Cluster parameters that will be updated
		 * @param[in] number_mh_steps				Number of Metropolis Hastings steps
		 */
		void update(
				clusters_t & clusters, 
				int number_mh_steps
			);
};

