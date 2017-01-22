typedef MatrixXd t_cluster;
typedef std:set<t_cluster> t_cluster_population;

/**
 * This class UpdateClusterPopulation deletes, adds, and adjusts clusters. This in contrast with UpdateCluster which 
 * only adjusts the parameters assigned to a cluster and will leave the number of clusters invariant. 
 */
class UpdateClusterPopulation {
	private:
		std::default_random_engine _generator;
	
		std::uniform_real_distribution<double> _distribution;

	public:
		/*!
		 * Construct update method for cluster population.
		 * @param[in] likelihood					Likelihood function to be used in update()
		 * @param[in] pred							Posterior predictive to be used in update()
		 * @param[in] prior							Prior to be used in update()
		 */
		UpdateClusterPopulation(
				Likelihood & likelihood,
				PosteriorPredictive & pred,
				Prior & prior
				);

		/*!
		 * Update the cluster population. The observation has to be deleted beforehand.
		 * @param[inout] clusters					Cluster parameters
		 * @param[in] observation					Observation to be considered for existing and new cluster
		 * @param[in] nonparametrics 				Sufficient statistics of the nonparametric prior (e.g. Dirichlet)
		 * @param[in] sample_pdf					Number of MH-steps
		 */
		void update(
				t_cluster_population & clusters, 
				t_data & observation,
				t_nonparametrics & nonparametrics, 
				t_sample_pdf & sample_pdf
				);
};
