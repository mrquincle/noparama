#include <np_update_clusters.h>

#include <np_sample_pdf.h>

/**
 * Initiate update cluster class and set likelihood and posterior predictive.
 */
UpdateClusters::UpdateClusters(
			distribution_t & likelihood,
			distribution_t & nonparametrics
		): 
			_likelihood(likelihood),
{
}

Suffies & UpdateClusters::propose() {
	// proposed parameters (can also be Brownian walk)
	return sample_pdf(_nonparametrics.base_distribution);
}

void UpdateClusters::update(
			clusters_t & clusters, 
			int number_mh_steps
		) 
{
	double likelihood;
	double proposed_likelihood;

	// This update procedure exists out of a number of MH steps
	for (int t = 0; t < number_mh_steps; ++t) {

		for (auto cluster: clusters) {
			// current likelihood for this cluster, compare with new sample
			// reuse _likelihood object
			_likelihood.init(cluster.Suffies());
			likelihood = _likelihood.get(cluster.data());
			
			Suffies & proposed_suffies = propose();
			_likelihood.init(proposed_suffies);
			proposed_likelihood = _likelihood.get(cluster.data());
			
			if (!likelihood) {
				// likelihood can be zero at start, in that case accept any proposal
				cluster.setSuffies(proposed_suffies);
				continue;
			} 

			double alpha = proposed_likelihood / likelihood;
	
			double reject = _distribution(_generator);

			if (reject < alpha) {
				// accept new cluster parameters
				cluster.setSuffies(proposed_suffies);
			} else {
				// reject, keep cluster as is
			}

		}
	}

}

