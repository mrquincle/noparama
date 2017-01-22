#include <np_update_cluster_population.h>
#include <np_likelihood.h>
#include <np_posterior_predictive.h>
#include <np_sample_pdf.h>

UpdateClusterPopulation::UpdateClusterPopulation(
		Likelihood & likelihood,
		PosteriorPredictive & pred
		): 
			_likelihood(likelihood),
			_predictive_posterior(pred) {
}

/**
 * Update the cluster population. The observation has to be deleted beforehand.
 */
void UpdateClusterPopulation::update(
		t_cluster_population & clusters, 
		t_data & observation
		t_nonparametrics & nonparametrics, 
		t_base_distribution & base_distribution, 
		t_sample_pdf & sample_pdf,
		t_prior & prior,
		) {

	// Create temporary data structure for K clusters
	size_t K = clusters.length();

	std::vector<double> weighted_likelihood(K);

	// Check all statistics given observation
	int i = 0;
	for (auto cluster: clusters; ++i) {
		
		weighted_likelihood[i] = _likelihood(prior, observation, cluster.getSufficientStatistics()) *
			cluster.count();
	}

	weighted_predictive_posterior = _predictive_posterior(observation, base_distribution) * alpha;

	// Calculate parameters for new cluster
	SufficientStatistics & sufficient_statistics = sample_pdf(base_distribution);
	Cluster new_cluster(sufficient_statistics);

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);

	double sum_likelihoods = sum(likelihoods);
	
	std::vector<double> cumsum_likelihood;

	double normalization_constant = sum_likelihoods + weighted_predictive_posterior;

	double prob_new_sample = weighted_predictive_posterior / normalization_constant;

	double reject = distribution(generator);

	if (reject < prob_new_sample) {
		// note, acceptance of new cluster does not depend on its likelihood
		clusters.push_back(new_cluster);
	} else {
		// pick existing cluster with weights 
		double pick = (reject - prob_new_sample);

		std::partial_sum(weighted_likelihood.begin(), weighted_likelihood.end(), cumsum_likelihood.begin());

		auto lower = std::lower_bound(cumsum_likelihood.begin(), cumsum_likelihood.end(), pick);
		size_t index = std::distance(cumsum_likelihood.begin(), lower);

		clusters[i].add(observation);
	}
}
