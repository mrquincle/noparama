#include <np_update_clusters.h>

/**
 * Initiate update cluster class and set likelihood and posterior predictive.
 */
UpdateClusters::UpdateClusters(
			distribution_t & likelihood,
			distribution_t & nonparametrics
		): 
			_likelihood(likelihood),
			_nonparametrics(nonparametrics),
			_distribution(0.0, 1.0)
{
}

Suffies & UpdateClusters::propose() {
	// proposed parameters (can also be Brownian walk)
	return  _nonparametrics.sample(_generator);
//	return sample_pdf(_nonparametrics.base_distribution);
}

void UpdateClusters::update(
			membertrix & cluster_matrix,
			int number_mh_steps
		) 
{
	double likelihood;
	double proposed_likelihood;

	// This update procedure exists out of a number of MH steps
	for (int t = 0; t < number_mh_steps; ++t) {

		auto clusters = cluster_matrix.getClusters();
		for (int k = 0; k < (int)clusters.size(); ++k) {
			// current likelihood for this cluster, compare with new sample
			// reuse _likelihood object
			_likelihood.init(clusters[k]->getSuffies());
			likelihood = _likelihood.probability(cluster_matrix.getData(k));
			
			Suffies & proposed_suffies = propose();
			_likelihood.init(proposed_suffies);
			proposed_likelihood = _likelihood.probability(cluster_matrix.getData(k));
			
			if (!likelihood) {
				// likelihood can be zero at start, in that case accept any proposal
				clusters[k]->setSuffies(proposed_suffies);
				continue;
			} 

			double alpha = proposed_likelihood / likelihood;
	
			double reject = _distribution(_generator);

			if (reject < alpha) {
				// accept new cluster parameters
				clusters[k]->setSuffies(proposed_suffies);
			} else {
				// reject, keep cluster as is
			}

		}
	}

}

