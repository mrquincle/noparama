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
			_predictive_posterior(pred) 
{
}

void UpdateClusters::update(
			t_cluster_population & clusters, 
			t_nonparametrics & nonparametrics, 
			t_prior & prior,
			int number_mh_steps
		) 
{
	double likelihood;
	double proposed_likelihood;

	// This update procedure exists out of a number of MH steps
	for (int t = 0; t < number_mh_steps; ++t) {

		for (auto cluster: clusters) {
			// current likelihood	
			likelihood = _likelihood(prior, cluster.data(), cluster.SufficientStatistics());
			
			// proposed likelihood
			SufficientStatistics & proposed_sufficient_statistics = sample_pdf(nonparametrics.base_distribution);
			proposed_likelihood = _likelihood(prior, cluster.data(), proposed_sufficient_statistics);
			
			if (!likelihood) {
				// likelihood can be zero at start, in that case accept any proposal
				cluster.setSufficientStatistics(proposed_sufficient_statistics);
				continue;
			} 

			double alpha = proposed_likelihood / likelihood;
	
			double reject = _distribution(_generator);

			if (reject < alpha) {
				// accept new cluster parameters
				cluster.setSufficientStatistics(proposed_sufficient_statistics);
			} else {
				// reject, keep cluster as is
			}

		}
	}

}

