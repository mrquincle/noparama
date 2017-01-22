typedef MatrixXd t_cluster;
typedef std:set<t_cluster> t_cluster_s;

/**
 * This class UpdateClusters only adjusts clusters. This in contrast with UpdateClusterPopluation which also
 * deletes and removes clusters and does not leave the number of clusters invariant. 
 */
class UpdateClusters {
	private:
		std::default_random_engine _generator;
	
		std::uniform_real_distribution<double> _distribution;

	public:
		/*!
		 * Construct update method for clusters.
		 * @param[in] likelihood					Likelihood function to be used in update()
		 * @param[in] pred							Posterior predictive to be used in update()
		 * @param[in] prior							Prior to be used in update()
		 */
		UpdateClusters(
				Likelihood & likelihood,
				PosteriorPredictive & pred,
				Prior & prior
			);

		/*!
		 * Update all cluster parameters. The number of clusters will not change.
		 * @param[inout] clusters					Cluster parameters that will be update.
		 * @param[in] nonparametrics 				Sufficient statistics of the nonparametric prior (e.g. Dirichlet)
		 * @param[in] number_mh_steps				Number of Metropolis Hastings steps
		 */
		void update(
				t_cluster_population & clusters, 
				t_nonparametrics & nonparametrics, 
				int number_mh_steps
			);
};

