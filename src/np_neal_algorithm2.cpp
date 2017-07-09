#include <np_neal_algorithm2.h>

NealAlgorithm2::NealAlgorithm2(
			distribution_t & likelihood,
			distribution_t & nonparametrics, 
			distribution_t & prior,
			PosteriorPredictive & pred
		): 
			_likelihood(likelihood),
			_nonparametrics(nonparametrics),
			_prior(prior),
			_predictive_posterior(pred),
			_distribution(0.0, 1.0)
{
}

void NealAlgorithm2::update(
			membertrix & cluster_matrix,
			const data_ids_t & data_ids
		) {

	assert (data_ids.size() == 1);
	data_id_t data_id = data_ids[0];
	
	// remove observation under consideration from cluster
	for (auto index: data_ids) {
		fout << "Retract observation " << index << " from cluster matrix" << endl;
		cluster_matrix.retract(index);
	}

	// Create temporary data structure for K clusters
	size_t K = clusters.length();

	std::vector<double> weighted_likelihood(K);

	data_t & observation = cluster_matrix.getData(data_id);

	// Check all statistics given observation
	int i = 0;
	for (auto cluster: cluster_matrix.getClusters(); ++i) {
		
		weighted_likelihood[i] = _likelihood(_prior, observation, cluster.getSuffies()) *
			cluster.count();
	}

	weighted_predictive_posterior = _predictive_posterior(observation, _nonparametrics.base_distribution) * 
		_nonparametrics.alpha;

	// Calculate parameters for new cluster
	Suffies & suffies = sample_pdf(_nonparametrics.base_distribution);
	
	cluster_t * new_cluster = new cluster_t(suffies);

	double sum_likelihoods = sum(likelihoods);
	
	std::vector<double> cumsum_likelihood;

	double normalization_constant = sum_likelihoods + weighted_predictive_posterior;

	double prob_new_sample = weighted_predictive_posterior / normalization_constant;

	double reject = _distribution(_generator);

	if (reject < prob_new_sample) {
		// note, acceptance of new cluster does not depend on its likelihood
		cluster_matrix.addCluster(new_cluster);
		cluster_matrix.assign(index, data_id);
	} else {
		// pick existing cluster with weights 
		double pick = (reject - prob_new_sample);

		std::partial_sum(weighted_likelihood.begin(), weighted_likelihood.end(), cumsum_likelihood.begin());

		auto lower = std::lower_bound(cumsum_likelihood.begin(), cumsum_likelihood.end(), pick);
		size_t index = std::distance(cumsum_likelihood.begin(), lower);

		cluster_matrix.reassign(index, data_id);
	}
}
