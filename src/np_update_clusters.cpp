#include <np_update_clusters.h>

#include <pretty_print.hpp>

using namespace std;

/**
 * Initiate update cluster class and set likelihood and posterior predictive.
 */
UpdateClusters::UpdateClusters(
			random_engine_t & generator,
			distribution_t & likelihood,
			distribution_t & nonparametrics
		): 
			_generator(generator),
			_likelihood(likelihood),
			_nonparametrics(nonparametrics),
			_distribution(0.0, 1.0)
{
}

Suffies * UpdateClusters::propose() {
	// proposed parameters (can also be Brownian walk)
	return  _nonparametrics(_generator);
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

		fout << "step " << t << endl;

		auto clusters = cluster_matrix.getClusters();
		//for (int k = 0; k < (int)clusters.size(); ++k) {
		for (auto cluster_pair: clusters) {
			auto const &key = cluster_pair.first;
			auto const &cluster = cluster_pair.second;

			if (cluster_matrix.empty(key)) continue;
			// current likelihood for this cluster, compare with new sample
			// reuse _likelihood object
			_likelihood.init(cluster->getSuffies());
			likelihood = _likelihood.probability(cluster_matrix.getData(key));
			fout << "likelihood cluster " << key << " is " << likelihood << endl;
			
			Suffies *proposed_suffies = propose();
			_likelihood.init(*proposed_suffies);
			proposed_likelihood = _likelihood.probability(cluster_matrix.getData(key));
			
			if (!likelihood) {
				// likelihood can be zero at start, in that case accept any proposal
				cluster->setSuffies(*proposed_suffies);
				continue;
			} 

			double alpha = proposed_likelihood / likelihood;
	
			double reject = _distribution(_generator);

			if (reject < alpha) {
				// accept new cluster parameters
				cluster->setSuffies(*proposed_suffies);
			} else {
				// reject, keep cluster as is
			}

		}
	}

}

