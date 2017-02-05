#include <np_neal_algorithm8.h>

#include <statistics/dirichlet.h>

NealAlgorithm8::NealAlgorithm8(
			distribution_t & likelihood,
			distribution_t & nonparametrics
		): 
			_likelihood(likelihood),
			_nonparametrics(nonparametrics)
{
	_alpha = ((dirichlet_distribution&)nonparametrics).getSuffies().alpha;
}

void NealAlgorithm8::update(
			membertrix & cluster_matrix,
			data_id_t data_id
		) {

	// existing clusters
	auto clusters = cluster_matrix.getClusters();

	// Create temporary data structure for K clusters
	size_t K = clusters.size();

	// Calculate parameters for M new clusters
	int M = 3;
	clusters_t new_clusters(M);
	for (int m = 0; m < M; ++m) {
		Suffies & suffies = _nonparametrics.sample(_generator);
		cluster_t *temp = new cluster_t(suffies);
		new_clusters[m] = temp;
	}	

	data_t & observation = cluster_matrix.getDatum(data_id);

	// Get (weighted) likelihoods for each existing and new cluster given the above observation
	std::vector<double> weighted_likelihood(K+M);
	for (int k = 0; k < (int)K; ++k) {
		_likelihood.init(clusters[k]->getSuffies());
		weighted_likelihood[k] = _likelihood.probability(observation) * cluster_matrix.count(k);
	}
	for (int m = 0; m < M; ++m) {
//		weighted_likelihood[i] = _likelihood(_prior, observation, cluster.getSuffies()) *
//			_alpha / M;
	}

	// sample uniformly from the vector "weighted_likelihood" according to the weights in the vector
	// hence [0.5 0.25 0.25] will be sampled in the ratio 2:1:1
	double pick = _distribution(_generator);
	std::vector<double> cumsum_likelihood;
	std::partial_sum(weighted_likelihood.begin(), weighted_likelihood.end(), cumsum_likelihood.begin());
	auto lower = std::lower_bound(cumsum_likelihood.begin(), cumsum_likelihood.end(), pick);
	size_t index = std::distance(cumsum_likelihood.begin(), lower);

	if (index >= K) {
		// pick new cluster
		cluster_t &new_cluster = *new cluster_t(new_clusters[index-K]->getSuffies());
		cluster_id_t cluster_index = cluster_matrix.addCluster(new_cluster);
		cluster_matrix.assign(cluster_index, data_id);
	} else {
		// pick existing cluster with given index
		cluster_matrix.assign(index, data_id);
	}

	// deallocate new clusters that have not been used
	for (int m = 0; m < M; ++m) {
		delete new_clusters[m];
	}
}
