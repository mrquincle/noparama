#include <np_update_clusters.h>
#include <np_likelihood.h>
#include <np_posterior_predictive.h>

/**
 * Initiate update cluster class and set likelihood and posterior predictive.
 */
UpdateClusters::UpdateClusters(
		Likelihood & likelihood,
		PosteriorPredictive & pred
		): 
			_likelihood(likelihood),
			_predictive_posterior(pred) {
}

/** 
 * Given a single observation update all cluster parameters 
 * @param clusters					Cluster parameters
 * @param nonparametrics 			Sufficient statistics of the nonparametric prior (e.g. Dirichlet)
 * @param number_mh_steps			Number of MH-steps
 */
UpdateClusters::update(
		t_cluster_population & clusters, 
		t_data & observation
		t_nonparametrics & nonparametrics, 
		t_base_distribution & base_distribution, 
		t_prior & prior,
		int number_mh_steps
		) {

	for (int t = 0; t < number_mh_steps; ++t)
		for (auto cluster: clusters) {

		}
	}

}

