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
	_verbosity = 4;
}

Suffies * UpdateClusters::propose() {
	// proposed parameters (can also be Brownian walk)
	return  _nonparametrics(_generator);
}

/**
 * Update cluster parameters. This function should leave the clusters invariant. The membership matrix should not
 * be different 
 */
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

		// here getClusters got corrupted...
		auto & clusters = cluster_matrix.getClusters();
		for (auto cluster_pair: clusters) {
			auto const &key = cluster_pair.first;
			
			fout << "Update cluster " << key << " parameters" << endl;

			auto const &cluster = cluster_pair.second;

			if (cluster_matrix.empty(key)) {
				// shouldn't happen normally
				foutvar(7) << "This cluster (" << key << ") does not exist!" << endl;
				continue;
			}

			// current likelihood for this cluster, compare with new sample and reuse the _likelihood object
			_likelihood.init(cluster->getSuffies());
			dataset_t *dataset = cluster_matrix.getData(key);
			likelihood = _likelihood.probability(*dataset);
			fout << "likelihood cluster " << key << " is " << likelihood << endl;
			
			Suffies *proposed_suffies = propose();
			_likelihood.init(*proposed_suffies);
			proposed_likelihood = _likelihood.probability(*dataset);
			
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

