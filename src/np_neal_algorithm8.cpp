#include <np_neal_algorithm8.h>

#include <iostream>
#include <iomanip>
#include <pretty_print.hpp>

using namespace std;

NealAlgorithm8::NealAlgorithm8(
			random_engine_t & generator,
			distribution_t & likelihood,
			dirichlet_process & nonparametrics
		): 
			_generator(generator),
			_likelihood(likelihood),
			_nonparametrics(nonparametrics)
{
	_alpha = _nonparametrics.getSuffies().alpha;

	_verbosity = 5;

	fout << "likelihood: " << _likelihood.getSuffies() << endl;
		
	_statistics.new_clusters_events = 0;
}

void NealAlgorithm8::update(
			membertrix & cluster_matrix,
			data_id_t data_id
		) {
	static int step = 0;
	fout << "Update step " << ++step << " in algorithm 8 for data item " << data_id << endl;

	// existing clusters
	// if data item removed one of the clusters, this is still in there with a particular weight
	auto clusters = cluster_matrix.getClusters();

	// Create temporary data structure for K clusters
	size_t K = clusters.size();
	std::vector<int> cluster_ids(K);

	fout << "Calculate likelihood for K=" << K << " existing clusters" << endl;

	fout << "Calculate parameters for M new clusters" << endl;
	// Calculate parameters for M new clusters
	int M = 3;
	clusters_t new_clusters(M);
	for (int m = 0; m < M; ++m) {
		Suffies *suffies = _nonparametrics.sample_base(_generator);
		cluster_t *temp = new cluster_t(*suffies);
		new_clusters[m] = temp;
		fout << "Suffies generated: " << new_clusters[m]->getSuffies() << endl;
	}	
	
	data_t & observation = *cluster_matrix.getDatum(data_id);
	fout << "Current observation " << data_id << " under consideration: " << observation << endl;

	// Get (weighted) likelihoods for each existing and new cluster given the above observation
	std::vector<double> weighted_likelihood(K+M);

	int k = 0;
	for (auto cluster_pair: clusters) {
		auto const &key = cluster_pair.first;
		auto const &cluster = cluster_pair.second;
	
		cluster_ids[k] = key;

		fout << "Obtained suffies from cluster " << key << ": " << cluster->getSuffies() << endl;

		// here _likelihood is sudden a niw distribution...
		_likelihood.init(cluster->getSuffies());
		fout << "calculate probability " << endl;
		double ll = _likelihood.probability(observation);
		fout << "which is: " << ll << endl;
		weighted_likelihood[k] = _likelihood.probability(observation) * cluster_matrix.count(key);
		/*
		fout << "Cluster ";
		cluster_matrix.print(key, std::cout);
		cout << '\t' << _likelihood.probability(observation) << endl;
		*/
		k++;
		
	}
	for (int m = 0; m < M; ++m) {
		_likelihood.init(new_clusters[m]->getSuffies());
		weighted_likelihood[K+m] = _likelihood.probability(observation) * _alpha / (double)M;
		fout << "New cluster " << m << 
			"[~#" << _alpha/M << "]: " << \
			_likelihood.probability(observation) << \
			new_clusters[m]->getSuffies() << endl;
	}

	// sample uniformly from the vector "weighted_likelihood" according to the weights in the vector
	// hence [0.5 0.25 0.25] will be sampled in the ratio 2:1:1
	fout << "Pick a cluster given their weights" << endl;
	double pick = _distribution(_generator);
	fout << "Pick " << pick << endl;
	std::vector<double> cumsum_likelihood(weighted_likelihood.size());
	fout << "Weighted likelihood: " << weighted_likelihood << endl;
	std::partial_sum(weighted_likelihood.begin(), weighted_likelihood.end(), cumsum_likelihood.begin());
	fout << "Cumulative weights: " << cumsum_likelihood << endl;

	pick = pick * cumsum_likelihood.back();
	fout << "Pick scaled with maximum cumsum: " << pick << endl;

	auto lower = std::lower_bound(cumsum_likelihood.begin(), cumsum_likelihood.end(), pick);
	size_t index = std::distance(cumsum_likelihood.begin(), lower);

	assert (index < cumsum_likelihood.size());
	fout << "Pick value just below pick: " << index << endl;
	
	if (index >= K) {
		// pick new cluster
		fout << "New cluster: " << index-K << endl;
		cluster_t *new_cluster = new cluster_t(new_clusters[index-K]->getSuffies());
		fout << "Add to membership matrix" << endl;
		cluster_id_t cluster_index = cluster_matrix.addCluster(new_cluster);
		fout << "Assign data to cluster with id " << cluster_index << endl;
		cluster_matrix.assign(cluster_index, data_id);
		_statistics.new_clusters_events++;
	} else {
		// pick existing cluster with given cluster_id
		fout << "Existing cluster" << endl;
		
		cluster_id_t cluster_id = cluster_ids[index];

		fout << "Assign data to cluster with id " << cluster_id << endl;
		assert (cluster_matrix.count(cluster_id) != 0);
		
		cluster_matrix.assign(cluster_id, data_id);
	}
	
	// deallocate new clusters that have not been used
	fout << "Deallocate temporary clusters" << endl;
	for (int m = 0; m < M; ++m) {
		delete new_clusters[m];
	}

	// check
	assert(cluster_matrix.assigned(data_id));
}

void NealAlgorithm8::printStatistics() {
	int verbosity = _verbosity;
	_verbosity = 0;
	fout << "Statistics:" << endl;
	fout << " # of new cluster events: " << _statistics.new_clusters_events << endl;
	_verbosity = verbosity;
} 
